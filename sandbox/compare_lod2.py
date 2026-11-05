"""Compare three LOD2 methods and two point clouds with Höganäs reference models.

Run from an environment with dtcc-core and its mesher installed:
    python sandbox/compare_lod2.py

Seven building-only VTUs share absolute EPSG:3008 / RH2000 coordinates. The
reference is not used for fitting. See compare_lod2.md for provenance and metrics.
"""
from pathlib import Path

DATA_DIR = Path("~/scratch/data/höganäs-2026/Export_DTCC").expanduser()
OUTPUT_DIR = Path(__file__).resolve().parent / "output" / "compare_lod2"
# Local EPSG:3008 coordinates. Edit this box to compare another area.
BOUNDS = (91450.0, 6231125.0, 91750.0, 6231425.0)
LOCAL_CRS = "EPSG:3008"
LM_CRS = "EPSG:3006"
TERRAIN_MARGIN = 20.0
MIN_BUILDING_AREA = 15.0
MIN_BUILDING_DETAIL = 0.5
SAMPLE_SPACING = 1.0  # Metres, fixed XY grid shared by all comparisons.
METHODS = ("catalogue", "ransac", "roofer")
ROOFER_EXECUTABLE = "roofer"  # Or an absolute path, e.g. "~/.local/bin/roofer".
ROOFER_TIMEOUT_SECONDS = 90

import csv
import json
import shutil
import subprocess
import tempfile
import time
import warnings
from collections import Counter
from copy import deepcopy

import dtcc_core as dtcc
import fiona
import laspy
import meshio
import numpy as np
import rasterio
import shapely
from pyproj import CRS, Transformer
from rasterio.windows import from_bounds
from scipy.spatial import cKDTree
from shapely.geometry import Polygon, box, shape
from shapely.ops import unary_union

from dtcc_core.reproject.reproject import reproject_array, reproject_surface


def log(message):
    print(f"[compare_lod2] {message}", flush=True)


def check_crs(crs, path):
    """Reject mismatched inputs; only LM coordinates are reprojected."""
    if crs is None:
        raise ValueError(f"Missing CRS: {path}")
    crs = CRS.from_user_input(crs)
    horizontal = crs.sub_crs_list[0] if crs.is_compound else crs
    # Esri's raster WKT rounds the GRS80 flattening; its normalized projection
    # is identical although pyproj's strict CRS equality rejects that rounding.
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', UserWarning)
        equivalent = horizontal.to_dict() == CRS(LOCAL_CRS).to_dict()
    equivalent &= horizontal.datum.name == CRS(LOCAL_CRS).datum.name
    if crs.is_compound and crs.sub_crs_list[1].to_epsg() != 5613:
        raise ValueError(f"Unexpected vertical CRS in {path}: {crs.sub_crs_list[1]}")
    if not equivalent:
        raise ValueError(f"Unexpected CRS in {path}: {crs}")


def polygon_rings(geometry):
    """Fiona exposes Multipatch faces through nested collections and polygons."""
    if geometry.type == "Polygon":
        yield geometry.coordinates
    elif geometry.type == "MultiPolygon":
        yield from geometry.coordinates
    elif geometry.type == "GeometryCollection":
        for child in geometry.geometries:
            yield from polygon_rings(child)
    else:
        raise ValueError(f"Unsupported reference geometry: {geometry.type}")


def triangulate_reference(rings):
    """Triangulate in the dominant plane, retaining the original XYZ vertices.

    No extrusion, smoothing, hole filling, or vertical shifting of the reference.
    Nonplanar polygons retain their vertices but have ambiguous interior surfaces.
    """
    rings = [np.asarray(r, dtype=float) for r in rings]
    rings = [r[:-1] if np.array_equal(r[0], r[-1]) else r for r in rings]
    # Exact A-B-A excursions enclose zero area (one supplied roof has a 4 mm
    # spur). Remove only those retraced edges, never repair crossing surfaces.
    for k, ring in enumerate(rings):
        ring = ring[np.any(ring != np.roll(ring, 1, axis=0), axis=1)]
        while len(ring) >= 3:
            spikes = np.flatnonzero(np.all(np.roll(ring, 1, axis=0) == np.roll(ring, -1, axis=0), axis=1))
            if not len(spikes):
                break
            i = int(spikes[0])
            ring = np.delete(ring, [i, (i+1) % len(ring)], axis=0)
        if len(ring) < 3:
            raise ValueError('Degenerate reference ring')
        rings[k] = ring
    outer = rings[0] - rings[0][0]
    normal = np.cross(outer, np.roll(outer, -1, axis=0)).sum(axis=0)
    length = np.linalg.norm(normal)
    if length < 1e-10:
        raise ValueError("Degenerate reference face")
    axes = [i for i in range(3) if i != np.argmax(abs(normal))]
    xyz = np.vstack(rings)
    uv = xyz[:, axes] - xyz[0, axes]
    cuts = np.cumsum([len(r) for r in rings])
    loops = np.split(uv, cuts[:-1])
    polygon = Polygon(loops[0], loops[1:])
    if not polygon.is_valid or polygon.area < 1e-10:
        raise ValueError("Invalid projected reference polygon")
    triangles = shapely.constrained_delaunay_triangles(polygon)
    coords = np.array([np.asarray(t.exterior.coords)[:3] for t in triangles.geoms])
    if not len(coords):
        raise ValueError("Empty reference triangulation")
    distance, faces = cKDTree(uv).query(coords.reshape(-1, 2))
    if distance.max() > 1e-7:
        raise ValueError("Reference triangulator introduced a vertex")
    if abs(sum(t.area for t in triangles.geoms) - polygon.area) > max(1e-8, polygon.area * 1e-8):
        raise ValueError("Reference triangulation lost coverage")
    deviation = float(np.max(abs((xyz - xyz[0]) @ (normal / length))))
    return xyz, faces.reshape(-1, 3), deviation


def concatenate(parts):
    vertices, faces, offset = [], [], 0
    for v, f in parts:
        vertices.append(v)
        faces.append(f + offset)
        offset += len(v)
    if not parts:
        raise ValueError("No triangles to export")
    return np.vstack(vertices), np.vstack(faces)


def save_vtu(name, vertices, faces, *, allow_empty=False, **fields):
    if not np.isfinite(vertices).all():
        raise ValueError(f"Nonfinite coordinates in {name}")
    xyz = vertices[faces]
    area2 = np.linalg.norm(np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0]), axis=1)
    if (not len(faces) and not allow_empty) or np.any(area2 <= 1e-12):
        raise ValueError(f"Empty or degenerate mesh: {name}")
    path = OUTPUT_DIR / f"{name}.vtu"
    meshio.write_points_cells(path, vertices, [("triangle", faces)],
                              cell_data={key: [np.asarray(value)] for key, value in fields.items()})
    log(f"Saved {path.name}: {len(vertices):,} vertices, {len(faces):,} triangles")


def load_reference(bounds):
    groups = {}
    inactive = 0
    for filename, kind in (("roofs_buildings_export.shp", 0),
                           ("walls_buildings_export.shp", 1),
                           ("new_consolidated_buildings_export.shp", 2)):
        path = DATA_DIR / "3D-Byggnader" / filename
        with fiona.open(path) as source:
            check_crs(source.crs_wkt, path)
            for feature in source:
                properties = dict(feature.properties)
                if properties.get("Aktualitet", 1) == 0:
                    inactive += 1
                    continue
                key = f"consolidated:{feature.id}" if kind == 2 else f"original:{properties['BDG_ID']}"
                group = groups.setdefault(key, dict(key=key, faces=[], records=[]))
                group['records'].append(dict(file=filename, feature_id=feature.id, properties=properties))
                group['faces'].extend((rings, kind) for rings in polygon_rings(feature.geometry))
    # Whole features only, so boundary buildings are not cut or reconstructed from partial roofs.
    region = box(*bounds)
    groups = [g for g in groups.values() if region.covers(box(*_group_bounds(g)))]
    parts, owners, kinds, records, failures, deviations = [], [], [], [], [], []
    for owner, group in enumerate(groups):
        records.append(dict(reference_id=owner, key=group['key'], records=group['records']))
        for face_id, (rings, kind) in enumerate(group['faces']):
            try:
                vertices, faces, deviation = triangulate_reference(rings)
            except ValueError as exc:
                failures.append(dict(key=group['key'], face=face_id, source_kind=kind, reason=str(exc)))
                continue
            parts.append((vertices, faces))
            owners.extend([owner] * len(faces))
            kinds.extend([kind] * len(faces))
            deviations.append(deviation)
    vertices, faces = concatenate(parts)
    save_vtu('reference', vertices, faces, reference_id=owners, source_kind=kinds)
    return vertices, faces, dict(buildings=len(groups), inactive_features_excluded=inactive,
        omitted_faces=failures, nonplanar_faces_over_1cm=int(np.sum(np.asarray(deviations) > .01)),
        maximum_planarity_deviation_m=max(deviations), sources=records)


def _group_bounds(group):
    points = np.vstack([np.asarray(r)[:, :2] for rings, _ in group['faces'] for r in rings])
    return (*points.min(axis=0), *points.max(axis=0))


def local_pointcloud(bounds):
    points, classifications = [], []
    counts = Counter()
    for path in sorted((DATA_DIR / 'LAS').glob('*.las')) + sorted((DATA_DIR / 'LAS').glob('*.laz')):
        with laspy.open(path) as source:
            check_crs(source.header.parse_crs(), path)
            for chunk in source.chunk_iterator(2_000_000):
                keep = ((chunk.x >= bounds[0]) & (chunk.x <= bounds[2]) &
                        (chunk.y >= bounds[1]) & (chunk.y <= bounds[3]))
                cls = np.asarray(chunk.classification)
                labels, numbers = np.unique(cls[keep], return_counts=True)
                counts.update(dict(zip(map(int, labels), map(int, numbers))))
                keep &= ~np.isin(cls, [7, 18]) & ~np.asarray(chunk.withheld, dtype=bool)
                points.append(np.column_stack((chunk.x[keep], chunk.y[keep], chunk.z[keep])))
                classifications.append(cls[keep])
    if not points or not sum(map(len, points)):
        raise ValueError("No local LAS points within the requested bounds")
    pc = dtcc.model.PointCloud(points=np.vstack(points), classification=np.concatenate(classifications))
    log(f"Local cloud: {len(pc.points):,} usable points; original class counts {dict(counts)}")
    return pc, dict(original_class_counts=dict(counts), usable_points=len(pc.points))


def load_terrain(bounds):
    path = DATA_DIR / 'Raster' / 'hojdraster_export.tif'
    with rasterio.open(path) as source:
        check_crs(source.crs, path)
        if not box(*source.bounds).covers(box(*bounds)):
            raise ValueError("Comparison bounds plus terrain margin exceed raster coverage")
        window = from_bounds(*bounds, transform=source.transform).round_offsets().round_lengths()
        data = source.read(1, window=window, masked=True)
        if np.ma.getmaskarray(data).any() or not np.isfinite(data).all():
            raise ValueError("Common terrain contains missing values")
        return dtcc.model.Raster(data=np.asarray(data, dtype=float), georef=source.window_transform(window),
                           crs=LOCAL_CRS, nodata=source.nodata)


def ground_alignment(cloud, terrain):
    """Report a datum sanity check without applying a fitted Z correction."""
    ground = cloud.points[np.asarray(cloud.classification) == 2]
    ground = ground[::max(1, len(ground)//20_000)]
    col, row = ~terrain.georef * (ground[:, 0], ground[:, 1])
    inside = (col >= 0) & (row >= 0) & (col < terrain.data.shape[1]) & (row < terrain.data.shape[0])
    ground, col, row = ground[inside], col[inside].astype(int), row[inside].astype(int)
    if not len(ground):
        return dict(samples=0, median_raster_minus_cloud_m=None)
    difference = terrain.data[row, col]-ground[:, 2]
    return dict(samples=len(ground), median_raster_minus_cloud_m=float(np.median(difference)),
                p90_absolute_difference_m=float(np.quantile(abs(difference), .9)))


def prepare_city(buildings, cloud, terrain):
    buildings = dtcc.builder.extract_roof_points(deepcopy(buildings), cloud)
    buildings = dtcc.builder.compute_building_heights(buildings, terrain, overwrite=True)
    city = dtcc.model.City()
    city.add_terrain(deepcopy(terrain))
    city.add_buildings(buildings, remove_outside_terrain=False)
    return city


def reconstruct(city, method, name):
    log(f"Computing LOD2 buildings: {method}, {name}")
    audit = {}
    start = time.monotonic()
    try:
        meshes = dtcc.builder.build_city_surface_mesh(deepcopy(city), lod='lod2', lod2_method=method,
            isolated=True, merge_meshes=False, stage_audit=audit,
            min_building_detail=MIN_BUILDING_DETAIL, min_building_area=MIN_BUILDING_AREA,
            report_mesh_quality=False)
    finally:
        write_json(f'{name}_audit.json', audit)
    buildings = meshes[1:]  # The public isolated pipeline returns terrain first.
    vertices, faces = concatenate([(m.vertices, m.faces) for m in buildings])
    catalogue = audit['roof_catalogue']
    roof_start = len(catalogue['patch_markers'])
    roof = np.concatenate([np.asarray(m.markers) >= roof_start for m in buildings])
    owners = np.concatenate([next(f.values for f in m.fields if f.name == 'building') for m in buildings])
    save_vtu(name, vertices, faces, building=owners, is_roof=roof.astype(np.uint8))
    summary = dict(seconds=time.monotonic()-start, buildings=len(buildings), triangles=len(faces),
                   outcomes=dict(Counter(row.get('final_outcome', 'unknown') for row in catalogue['sources'])),
                   retained_nonflat_fits=catalogue['retained_nonflat_fits'],
                   block_fallbacks=len(catalogue['block_fallbacks']),
                   shapes=dict(Counter(row.get('shape') or 'unfitted' for row in catalogue['sources'])))
    log(f"{name}: {summary}")
    return vertices, faces[roof], summary


def resolve_roofer():
    """Check the standalone CLI before downloading inputs or changing outputs."""
    setup = ('Roofer is required for this comparison. Install it on macOS/Linux with:\n'
        'curl -fsSL https://raw.githubusercontent.com/3DBAG/roofer/refs/heads/main/distribution/install.sh | ROOFER_VERSION=1.0.0 sh\n'
        'export PATH="$HOME/.local/bin:$PATH"\n'
        'roofer --version\n'
        'Or set ROOFER_EXECUTABLE to its installed path at the top of this script.\n'
        'Windows instructions: https://innovation.3dbag.nl/roofer/getting_started.html')
    executable = shutil.which(str(Path(ROOFER_EXECUTABLE).expanduser()))
    if executable is None:
        raise SystemExit(setup)
    try:
        version = subprocess.run([executable, '--version'], check=True,
                                 capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError) as exc:
        raise SystemExit(f'Cannot run Roofer: {exc}\n{setup}') from exc
    return dict(executable=str(Path(executable).resolve()),
                version=(version.stdout or version.stderr).strip())


def read_roofer(folder, source_id):
    """Read one building's CityJSONSeq LOD2.2 skins in absolute coordinates."""
    parts, roof_flags, attributes = [], [], {}
    found = False
    for path in sorted(Path(folder).glob('*.jsonl')):
        transform = None
        for line in path.read_text().splitlines():
            if not line.strip():
                continue
            data = json.loads(line)
            transform = data.get('transform', transform)
            objects = data.get('CityObjects', {})
            if not objects:
                continue
            roots = {key for key, obj in objects.items() if obj['type'] == 'Building'}
            if roots != {source_id}:
                raise ValueError(f'Unexpected Roofer building ID: {sorted(roots)}; expected {source_id}')
            found = True
            attributes.update(objects[source_id].get('attributes', {}))
            selected = {source_id}
            pending = [source_id]
            while pending:
                key = pending.pop()
                children = list(objects[key].get('children', []))
                children += [k for k, obj in objects.items() if key in obj.get('parents', [])]
                for child in children:
                    if child not in objects:
                        raise ValueError(f'Missing Roofer building part: {child}')
                    if child not in selected:
                        selected.add(child)
                        pending.append(child)
            vertices = np.asarray(data.get('vertices', []), dtype=float).reshape(-1, 3)
            if transform is not None:
                vertices = vertices*np.asarray(transform['scale']) + np.asarray(transform['translate'])
            for key, obj in objects.items():
                if key not in selected:
                    continue
                for geometry in obj.get('geometry', []):
                    if str(geometry.get('lod')) != '2.2':
                        continue
                    semantics = geometry['semantics']
                    if geometry['type'] == 'Solid':
                        shells, labels = geometry['boundaries'], semantics['values']
                    elif geometry['type'] == 'MultiSurface':
                        shells, labels = [geometry['boundaries']], [semantics['values']]
                    else:
                        raise ValueError(f'Unsupported Roofer geometry: {geometry["type"]}')
                    for shell, surface_labels in zip(shells, labels, strict=True):
                        for face, label in zip(shell, surface_labels, strict=True):
                            if label is None:
                                raise ValueError('Roofer surface has no semantic label')
                            kind = semantics['surfaces'][label]['type']
                            if kind == 'GroundSurface':
                                continue
                            if kind not in ('RoofSurface', 'WallSurface'):
                                raise ValueError(f'Unsupported Roofer surface: {kind}')
                            v, f, _ = triangulate_reference([vertices[ring] for ring in face])
                            parts.append((v, f))
                            roof_flags.extend([kind == 'RoofSurface']*len(f))
    if not found:
        raise ValueError(f'No Roofer result for building {source_id}')
    if not parts or not any(roof_flags):
        return None, attributes
    vertices, faces = concatenate(parts)
    return (vertices, faces, np.asarray(roof_flags, dtype=bool)), attributes


def reconstruct_roofer(city, name, tool):
    """Run Roofer on the same observations, preserving failures as missing output."""
    start = time.monotonic()
    raw_root = (OUTPUT_DIR/f'{name}_raw').resolve()
    raw_root.mkdir(parents=True, exist_ok=True)
    audit = dict(tool=tool, footprint_policy='original LM outlines, Roofer simplification disabled',
        point_policy='DTCC extracted observations labelled class 6; density ceiling 10000/m2',
        floor_policy='common terrain ground_height attribute', buildings=[])
    parts, owners, roof_flags = [], [], []
    for index, building in enumerate(city.buildings):
        cloud = building.point_cloud
        points = np.asarray(cloud.points) if cloud is not None else np.empty((0, 3))
        row = dict(building=index, source_id=building.id, lm_id=building.attributes.get('objektidentitet'),
                   points=len(points), outcome='missing')
        audit['buildings'].append(row)
        if not len(points):
            row['reason'] = 'no roof observations'
            continue
        directory = Path(tempfile.mkdtemp(prefix=f'b{index:03d}_', dir=raw_root))
        row['raw_directory'] = str(directory)
        try:
            header = laspy.LasHeader(point_format=3, version='1.2')
            header.scales = np.full(3, .001)
            header.offsets = np.floor(points.min(axis=0))
            header.add_crs(CRS(LOCAL_CRS))
            las = laspy.LasData(header)
            las.x, las.y, las.z = points.T
            las.classification = np.full(len(points), 6, dtype=np.uint8)
            las.write(directory/'points.las')
            with fiona.open(directory/'footprint.gpkg', 'w', driver='GPKG', crs=LOCAL_CRS,
                schema=dict(geometry='Polygon', properties={'id': 'str', 'ground': 'float'})) as dst:
                dst.write(dict(geometry=building.lod0.to_polygon(simplify=0).__geo_interface__,
                    properties=dict(id=building.id, ground=float(building.attributes['ground_height']))))
            command = [tool['executable'], '-j', '1', '--no-simplify', '--no-clip-terrain',
                '--ceil-point-density', '10000', '--id-attribute', 'id',
                '--h-terrain-attribute', 'ground', '--h-terrain-strategy', 'user',
                '--lod11-fallback-time', str(int(ROOFER_TIMEOUT_SECONDS*2000/3)),
                str(directory/'points.las'), str(directory/'footprint.gpkg'), str(directory/'result')]
            row['command'] = command
            with (directory/'roofer.log').open('w') as log_file:
                process = subprocess.run(command, check=True, timeout=ROOFER_TIMEOUT_SECONDS,
                                         stdout=log_file, stderr=subprocess.STDOUT, cwd=directory)
            row['returncode'] = process.returncode
            mesh, row['attributes'] = read_roofer(directory/'result', building.id)
            if mesh is None:
                row['reason'] = 'no returned roof geometry'
                continue
            v, f, roof = mesh
            parts.append((v, f))
            owners.extend([index]*len(f))
            roof_flags.extend(roof)
            mode = row['attributes'].get('rf_extrusion_mode')
            row.update(outcome='reconstructed' if mode == 'standard' else 'fallback', triangles=len(f))
        except subprocess.TimeoutExpired as exc:
            row.update(outcome='timeout', error=str(exc))
        except (OSError, subprocess.CalledProcessError, ValueError, KeyError, IndexError) as exc:
            row.update(outcome='failed', error=str(exc))
            if isinstance(exc, subprocess.CalledProcessError):
                row['returncode'] = exc.returncode
    write_json(f'{name}_audit.json', audit)
    vertices, faces = concatenate(parts) if parts else (np.empty((0, 3)), np.empty((0, 3), dtype=int))
    roof = np.asarray(roof_flags, dtype=bool)
    save_vtu(name, vertices, faces, allow_empty=True,
             building=np.asarray(owners, dtype=int), is_roof=roof.astype(np.uint8))
    summary = dict(seconds=time.monotonic()-start, buildings=len(parts), triangles=len(faces),
                   outcomes=dict(Counter(row['outcome'] for row in audit['buildings'])))
    log(f'{name}: {summary}')
    return vertices, faces[roof], summary


def top_heights(vertices, faces, xy):
    """Upper-envelope vertical intersections, independent of triangle winding.

    Vertical walls have no projected area. Ground caps lie below roofs. This
    intentionally compares roof surfaces, not buried walls or nearest 3D points.
    """
    triangles = vertices[faces]
    polygons = shapely.polygons(triangles[:, :, :2])
    valid = shapely.area(polygons) > 1e-8
    triangles, polygons = triangles[valid], polygons[valid]
    point_ids, triangle_ids = shapely.STRtree(polygons).query(shapely.points(xy), predicate='intersects')
    heights = np.full(len(xy), -np.inf)
    if len(point_ids):
        t = triangles[triangle_ids]
        delta = xy[point_ids]-t[:, 0, :2]
        a, b = t[:, 1, :2]-t[:, 0, :2], t[:, 2, :2]-t[:, 0, :2]
        det = a[:, 0]*b[:, 1]-a[:, 1]*b[:, 0]
        u = (delta[:, 0]*b[:, 1]-delta[:, 1]*b[:, 0])/det
        v = (a[:, 0]*delta[:, 1]-a[:, 1]*delta[:, 0])/det
        z = t[:, 0, 2]+u*(t[:, 1, 2]-t[:, 0, 2])+v*(t[:, 2, 2]-t[:, 0, 2])
        np.maximum.at(heights, point_ids, z)
    heights[~np.isfinite(heights)] = np.nan
    return heights


def error_stats(error):
    if not len(error):
        return dict(samples=0, mae_m=None, rmse_m=None, bias_m=None)
    return dict(samples=len(error), mae_m=float(np.mean(abs(error))),
                rmse_m=float(np.sqrt(np.mean(error**2))), bias_m=float(np.mean(error)))


def compare(buildings, reference, variants, bounds):
    x = np.arange(bounds[0]+SAMPLE_SPACING/2, bounds[2], SAMPLE_SPACING)
    y = np.arange(bounds[1]+SAMPLE_SPACING/2, bounds[3], SAMPLE_SPACING)
    xy = np.stack(np.meshgrid(x, y), axis=-1).reshape(-1, 2)
    polygons = [b.lod0.to_polygon(simplify=0) for b in buildings]
    xy = xy[shapely.contains_xy(unary_union(polygons), xy[:, 0], xy[:, 1])]
    truth = top_heights(*reference, xy)
    predictions = {name: top_heights(*mesh, xy) for name, mesh in variants.items()}
    # All variants use exactly the same horizontal samples.
    common = np.isfinite(truth)
    for z in predictions.values():
        common &= np.isfinite(z)
    result = dict(lm_footprint_samples=len(xy), reference_samples=int(np.isfinite(truth).sum()),
                  common_samples=int(common.sum()), sample_spacing_m=SAMPLE_SPACING, variants={})
    rows = []
    for name, z in predictions.items():
        result['variants'][name] = dict(error_stats(z[common]-truth[common]),
            generated_samples=int(np.isfinite(z).sum()),
            pairwise_common_samples=int((np.isfinite(z) & np.isfinite(truth)).sum()))
        for building, polygon in zip(buildings, polygons):
            inside = shapely.contains_xy(polygon, xy[:, 0], xy[:, 1])
            keep = common & inside
            rows.append(dict(variant=name, model_building_id=building.id,
                lm_building_id=building.attributes.get('objektidentitet', building.id),
                footprint_samples=int(inside.sum()), reference_samples=int((inside & np.isfinite(truth)).sum()),
                **error_stats(z[keep]-truth[keep])))
    with (OUTPUT_DIR/'roof_errors.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return result


def _json_default(value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_json(filename, value):
    with (OUTPUT_DIR/filename).open('w') as f:
        json.dump(value, f, indent=2, default=_json_default)


def main():
    if not METHODS or len(set(METHODS)) != len(METHODS) or set(METHODS)-{'catalogue', 'ransac', 'roofer'}:
        raise ValueError('METHODS must select unique catalogue, ransac, and/or roofer methods')
    if ROOFER_TIMEOUT_SECONDS <= 0:
        raise ValueError('ROOFER_TIMEOUT_SECONDS must be positive')
    roofer = resolve_roofer() if 'roofer' in METHODS else None
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    bounds = BOUNDS
    if len(bounds) != 4 or not np.isfinite(bounds).all() or bounds[0] >= bounds[2] or bounds[1] >= bounds[3]:
        raise ValueError('BOUNDS must be (xmin, ymin, xmax, ymax) in local coordinates')
    if SAMPLE_SPACING <= 0 or TERRAIN_MARGIN <= 0:
        raise ValueError('SAMPLE_SPACING and TERRAIN_MARGIN must be positive')
    padded = tuple(v + (-TERRAIN_MARGIN if i < 2 else TERRAIN_MARGIN) for i, v in enumerate(bounds))
    write_json('comparison.json', dict(status='preparing', bounds=bounds, data_dir=DATA_DIR))
    log(f"Block: {bounds}; outputs {LOCAL_CRS}, absolute RH2000 Z (no height shift)")
    terrain = load_terrain(padded)
    reference_v, reference_f, reference_info = load_reference(bounds)
    log(f"Reference: {reference_info['buildings']} building records, {len(reference_info['omitted_faces'])} omitted faces")
    write_json('reference_sources.json', reference_info)
    download_bounds = Transformer.from_crs(LOCAL_CRS, LM_CRS, always_xy=True).transform_bounds(*padded)
    log('Reading LM footprints and point cloud through the DTCC data service/cache')
    buildings = dtcc.io.data.download_footprints(dtcc.Bounds(*download_bounds))
    for b in buildings:
        b.add_geometry(reproject_surface(b.lod0, LM_CRS, LOCAL_CRS), dtcc.model.GeometryType.LOD0)
    region = box(*bounds)
    buildings = [b for b in buildings if region.covers(b.lod0.to_polygon(simplify=0))
                 and b.lod0.to_polygon(simplify=0).area >= MIN_BUILDING_AREA]
    if not buildings:
        raise ValueError('No complete LM footprints within comparison bounds')
    log(f"Using the same {len(buildings)} complete LM footprints for all reconstructions")
    lm = dtcc.io.data.download_pointcloud(dtcc.Bounds(*download_bounds))
    lm_points = reproject_array(lm.points, LM_CRS, LOCAL_CRS)
    keep = ~np.isin(lm.classification, [7, 18])
    lm = dtcc.model.PointCloud(points=lm_points[keep], classification=lm.classification[keep])
    local, local_info = local_pointcloud(padded)
    manifest = dict(status='running', data_dir=DATA_DIR, bounds=bounds, horizontal_crs=LOCAL_CRS, vertical_crs='EPSG:5613 (RH2000)',
        lm_horizontal_crs=LM_CRS, vertical_operation='None; local RH2000 documented in läsmig.txt; LM datum follows product specification (LAZ header lacks CRS)',
        terrain='Common supplied height raster; reference models not used in fitting',
        lm_buildings=[dict(model_id=b.id, lm_id=b.attributes.get('objektidentitet'),
                           bounds=b.lod0.to_polygon(simplify=0).bounds) for b in buildings],
        lm_points=len(lm.points), local_cloud=local_info,
        reference=reference_info, variants={}, roof_observations={}, methods=list(METHODS), roofer=roofer,
        ground_alignment={name: ground_alignment(pc, terrain) for name, pc in [('lm', lm), ('local', local)]})
    write_json('comparison.json', manifest)
    variants = {}
    for cloud_name, cloud in [('lm', lm), ('local', local)]:
        log(f"Extracting roof observations: {cloud_name}, {len(cloud.points):,} cloud points")
        city = prepare_city(buildings, cloud, terrain)
        manifest['roof_observations'][cloud_name] = {b.id: len(b.point_cloud.points) for b in city.buildings}
        for method in METHODS:
            name = f'{method}_{cloud_name}'
            vertices, faces, summary = (reconstruct_roofer(city, name, roofer) if method == 'roofer'
                                        else reconstruct(city, method, name))
            variants[name] = (vertices, faces)
            manifest['variants'][name] = summary
            write_json('comparison.json', manifest)
    manifest['comparison'] = compare(buildings, (reference_v, reference_f), variants, bounds)
    manifest['status'] = 'complete'
    write_json('comparison.json', manifest)
    log(json.dumps(manifest['comparison'], indent=2))
    log('Complete: open reference.vtu and the selected method_cloud.vtu files together in ParaView.')


if __name__ == '__main__':
    main()
