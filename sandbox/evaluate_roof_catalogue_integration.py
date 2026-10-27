"""Reproduce complete-city and fixed-case catalogue evaluation from saved evidence.

Run with the dtcc-core environment. The input is the trusted local 341-building
roof benchmark cache (pickled models, train/test NPZs and manifests). No Roofer
installation is used. All fits use train only; holdouts are scored afterwards.
"""
import argparse
from collections import Counter
from contextlib import contextmanager
import html
import hashlib
import platform
import json
from pathlib import Path
import pickle
import time
from unittest.mock import patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
from scipy.spatial import cKDTree
import shapely
from shapely.geometry import box

from dtcc_core import builder
from dtcc_core.builder.geometry_builders import meshes
from dtcc_core.builder.meshing.roof_validation import validate_shell
from dtcc_core.builder.roofs import fit_roof
from dtcc_core.model import Bounds, City, GeometryType, PointCloud, Raster, Terrain


def dump(value, path):
    path.write_text(json.dumps(value, indent=2, allow_nan=False))


def terrain_raster(manifest, output):
    """Nearest classified ground observations on a 4 m north-to-south grid."""
    saved = output / f"{manifest['area']}-terrain.npz"
    if saved.exists():
        data = np.load(saved)['data']
    else:
        import laspy
        ground = []
        xmin, ymin, xmax, ymax = manifest['bounds']
        for filename in manifest['lidar_files']:
            if not Path(filename).exists():
                raise FileNotFoundError(f'Local terrain evidence unavailable: {filename}')
            with laspy.open(filename) as reader:
                for chunk in reader.chunk_iterator(1000000):
                    xyz = np.c_[np.asarray(chunk.x), np.asarray(chunk.y), np.asarray(chunk.z)]
                    keep = ((np.asarray(chunk.classification)==2) &
                            (xyz[:, 0]>=xmin-10) & (xyz[:, 0]<=xmax+10) &
                            (xyz[:, 1]>=ymin-10) & (xyz[:, 1]<=ymax+10))
                    ground.append(xyz[keep])
        ground = np.vstack(ground)
        if not len(ground):
            raise ValueError('No classified ground observations available')
        x = np.linspace(xmin, xmax, 126)
        y = np.linspace(ymax, ymin, 126)
        xx, yy = np.meshgrid(x, y)
        nearest = cKDTree(ground[:, :2]).query(np.c_[xx.ravel(), yy.ravel()])[1]
        data = ground[nearest, 2].reshape(xx.shape)
        np.savez_compressed(saved, data=data)
    raster = Raster(data=data)
    raster.set_bounds(Bounds(*manifest['bounds']))
    return raster


def city_of(buildings, raster):
    city = City()
    terrain = Terrain()
    terrain.add_geometry(raster, GeometryType.RASTER)
    city.add_terrain(terrain)
    city.add_buildings(buildings)
    return city


def crop_raster(raster, bounds):
    xmin, ymin, xmax, ymax = bounds
    xs = np.linspace(xmin, xmax, 32)
    ys = np.linspace(ymax, ymin, 32)
    bxmin, bymin, bxmax, bymax = meshes._raster_bounds_tuple(raster)
    data = np.array([[raster.get_value(np.clip(x, bxmin+.001, bxmax-.001),
                                      np.clip(y, bymin+.001, bymax-.001)) for x in xs] for y in ys])
    result = Raster(data=data)
    result.set_bounds(Bounds(*bounds))
    return result


def build(city, *, catalogue, mesher=None):
    audit, captured = {}, []
    original = meshes._build_city_surface_mesh_from_ground_mesh
    def capture(**kwargs):
        result = original(**kwargs)
        captured.append(result)
        return result
    start = time.perf_counter()
    error, result = None, None
    try:
        with patch.object(meshes, '_build_city_surface_mesh_from_ground_mesh', capture):
            result = builder.build_city_surface_mesh(city, lod='lod2' if catalogue else 'lod1',
                stage_audit=audit, report_mesh_quality=False, max_mesh_size=8.,
                building_mesh_triangle_size=4., min_building_area=15., mesher=mesher)
    except Exception as exc:
        error = str(exc)
        if captured:
            result = captured[-1]  # Retain invalid assemblies for diagnosis only.
    return result, dict(error=error, seconds=time.perf_counter()-start, audit=audit)


def mesh_score(mesh, points, roof_start, roof_end):
    if mesh is None or not len(points):
        return None
    markers = np.asarray(mesh.markers)
    xyz = mesh.vertices[mesh.faces[(markers >= roof_start) & (markers < roof_end)]]
    if not len(xyz):
        return None
    xy = xyz[:, :, :2]
    tree = shapely.STRtree(shapely.polygons(np.concatenate([xy, xy[:, :1]], axis=1)))
    pi, ti = tree.query(shapely.points(points[:, :2]), predicate='intersects')
    predictions = np.full(len(points), np.nan)
    if len(pi):
        tri = xyz[ti]
        slopes = np.linalg.solve(tri[:, 1:, :2]-tri[:, :1, :2],
                                 (tri[:, 1:, 2]-tri[:, :1, 2])[..., None])[..., 0]
        z = tri[:, 0, 2] + np.sum(slopes*(points[pi, :2]-tri[:, 0, :2]), axis=1)
        predictions[pi] = z
    finite = np.isfinite(predictions)
    e = abs(predictions[finite]-points[finite, 2])
    return dict(point_coverage=float(finite.mean()), n=int(finite.sum()),
                median_ae=float(np.median(e)) if len(e) else None,
                p90_ae=float(np.quantile(e, .9)) if len(e) else None)


def render(ax, mesh, points, title, origin, count=None):
    ax.set_title(title, fontsize=10)
    if mesh is not None:
        xyz = mesh.vertices[mesh.faces] - np.r_[origin, 0.]
        markers = np.asarray(mesh.markers)
        # Fit-independent wall/roof classification by triangle normals also
        # works with the default pipeline's scalar region count.
        cross = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
        walls = np.linalg.norm(cross[:, :2], axis=1) > abs(cross[:, 2])*100
        terrain = markers < 0
        colors = np.tile([.25, .60, .82, 1.], (len(xyz), 1))
        colors[walls] = [.84, .70, .49, 1.]
        colors[terrain] = [.71, .76, .66, .48]
        if len(points):
            colors[:, 3] = .18
            colors[terrain, 3] = .08
        ax.add_collection3d(Poly3DCollection(xyz, facecolors=colors,
                           edgecolors=(.12,.16,.20,.28), linewidths=.22))
        ax.auto_scale_xyz(xyz[:, :, 0].ravel(), xyz[:, :, 1].ravel(), xyz[:, :, 2].ravel())
    else:
        ax.text2D(.1, .45, 'ASSEMBLY FAILED\nNo shell returned', transform=ax.transAxes, color='darkred')
    if len(points):
        ax.scatter(points[:, 0]-origin[0], points[:, 1]-origin[1], points[:, 2],
                   s=1.5, c=points[:, 2], cmap='viridis', alpha=.65, depthshade=False)
    ax.view_init(elev=28, azim=-55)
    ax.set_box_aspect((1, 1, .6))
    ax.set_xlabel('x offset (m)', fontsize=7)
    ax.set_ylabel('y offset (m)', fontsize=7)
    ax.set_zlabel('elevation (m)', fontsize=7)
    ax.tick_params(labelsize=6)


def write_report(output, summary, gallery, areas):
    """Render saved measurements, retaining failed and interrupted experiments."""
    extra = output/'lund-triangle-audit.json'
    extra_section = ''
    if extra.exists():
        check = json.loads(extra.read_text())
        summary['additional_triangle_check'] = check
        dump(summary, output/'summary.json')
        extra_section = ('<section><h2>Prior incomplete full-Lund Triangle check (before corrections)</h2>'
            '<p>The optional-backend run was stopped after approximately 92 seconds total, '
            'without a returned triangulation. This is an incomplete full-city check, '
            'separate from the passing synthetic acceptance tests for both backends.</p>'
            '<p><a href="lund-triangle-audit.json">Retained stopped-run diagnostics</a></p></section>')
    sections = []
    def error_text(score):
        if not score or score.get('median_ae') is None:
            return 'unavailable'
        return f"{score['median_ae']:.3f} / {score['p90_ae']:.3f} m"
    city_sections = []
    for area, result in areas.items():
        cat = result['catalogue']
        attempts = cat.get('assemblies', [])
        failures = [a.get('error') for a in attempts if not a.get('valid')]
        shell = attempts[-1].get('shell', {}) if attempts else {}
        contract = shell.get('quality_contract', {})
        worst = shell.get('mesh_audit', {}).get('worst_element_face')
        patch_count = len(cat.get('patch_markers', []))
        marker = worst.get('marker') if worst else None
        owner = (cat['patch_markers'][int(marker)%patch_count]
                 if marker is not None and 0<=marker<2*patch_count else None)
        inspection = dict(warning=shell.get('quality_warning'), contract=contract,
                          worst_face=worst, owner=owner)
        if marker is not None and 0<=marker<patch_count:
            vertices = np.asarray(worst['vertices'])
            inspection['wall_panel_width_m'] = float(np.linalg.norm(vertices[:,None,:2]-vertices[None,:,:2],axis=2).max())
            inspection['limitation'] = ('A mandatory roof breakpoint creates a narrow wall panel. '
                'The shared boundary-only wall zipper retains a quality warning; '
                'topology validation does not certify engineering mesh quality.')
        quality = html.escape(json.dumps(inspection, indent=2))
        source_rows = []
        for source in cat.get('sources', []):
            reasons = source['diagnostics'].get('rejection_reasons', []) + [
                f['reason'] for f in source.get('block_fallbacks', [])]
            source_rows.append('<tr>'+''.join(f'<td>{html.escape(str(v))}</td>' for v in (
                source.get('source_key',source.get('source_keys')), source.get('shape'), source['outcome'],
                source.get('final_outcome'), source.get('retained_fit_area',0), '; '.join(reasons)))+'</tr>')
        source_table = ('<details><summary>Every assembled source: fit retention and fallback reasons</summary>'
            '<table><tr><th>Source</th><th>Provisional model</th><th>Fit outcome</th><th>Final outcome</th><th>Retained fit area (m²)</th><th>Fallback reasons</th></tr>'
            + ''.join(source_rows)+'</table></details>')
        city_sections.append(f'<section><h2>Complete {area} surface</h2>'
            f'<p>Runtime: {result["seconds"]:.2f} s. Topology valid: {shell.get("topology_valid",False)}. '
            f'Minimum triangle quality: {shell.get("minimum_quality")}; quality warning: {shell.get("quality_warning",False)}.</p>'
            f'<details><summary>Quality contract and thin-triangle inspection</summary><pre>{quality}</pre></details>'
            f'<p>Retained fits: {cat.get("retained_fits",0)}; retained non-flat fits: {cat.get("retained_nonflat_fits",0)}.</p>'
            f'<p>{len(cat.get("blocks", []))} cleaned blocks; {len(cat.get("block_fallbacks", []))} whole-block fallbacks. '
            f'Final error: {html.escape(str(result["error"]))}.</p>'
            +source_table+f'<p>Rejected assembly attempts retained: {html.escape(str(failures))}</p>'
            f'<p><a href="{area}-audit.json">Full city audit including every fallback reason</a> · '
            f'<a href="{area}-catalogue.vtu">Final assembled VTU</a></p></section>')
    for number, g in enumerate(gallery):
        info = json.dumps({k:v for k,v in g.items() if k!='image'}, indent=2)
        score_by = {r['index']:r for r in g['case_holdout']}
        table_rows = []
        for fit in g['fits']:
            score = score_by[fit['index']]
            rejection = '; '.join(fit['fit'].get('rejection_reasons', []))
            table_rows.append('<tr>'+''.join(f'<td>{html.escape(str(value))}</td>' for value in (
                fit['index'], fit['shape'], 'accepted' if fit['fit']['accepted'] else 'rejected → LOD1',
                error_text(fit['heldout_provisional']), error_text(score['lod1']),
                error_text(score['catalogue']), rejection))+'</tr>')
        table = '<table><tr><th>Source</th><th>Provisional shape</th><th>Fit outcome</th><th>Provisional holdout median / p90</th><th>LOD1 holdout median / p90</th><th>Assembled catalogue holdout median / p90</th><th>Fit rejection</th></tr>'+''.join(table_rows)+'</table>'
        shell = g['shell'] or {}
        mesh_status = (f"Topology valid: {shell.get('valid',False)}; "
            f"minimum quality: {shell.get('minimum_quality')}; "
            f"quality warning: {shell.get('quality_warning')}; "
            f"block fallbacks: {len(g['fallback_blocks'])}.")
        sections.append(f'<section><h2>{html.escape(str(g["case"]))}</h2><p>{html.escape(mesh_status)}</p>{table}<img src="{g["image"]}"><details><summary>Fit quality, mesh validity and fallback diagnostics</summary><pre>{html.escape(info)}</pre></details><a href="case-{number:02d}.json">Full audit</a></section>')
    report = '''<!doctype html><meta charset="utf-8"><title>Integrated roof catalogue evaluation</title><style>body{font:16px system-ui;max-width:1500px;margin:32px auto;background:#fafafa;color:#182c3a}img{width:100%}section{background:white;padding:20px;margin:24px 0}pre{white-space:pre-wrap;font-size:12px}h1,h2{font-weight:600}summary{cursor:pointer}table{border-collapse:collapse;font-size:13px}td,th{padding:8px;border:1px solid #ccd6dc;text-align:left}</style><h1>Integrated roof catalogue — complete surfaces</h1><p>All thirteen fixed cases are retained, including both five-building groups and failures. Blue: roofs; tan: walls; green: terrain. Right: held-out LiDAR. A geometrically valid mesh does not establish a good roof fit. Provisional fitted errors and the final assembled errors are reported separately; rejected fits and whole-block fallbacks remain visible.</p>'''
    report += '<p>Provisional fits: '+str(summary.get('provisional_fits'))+'; accepted fits: '+str(summary['accepted_fits'])+' of '+str(summary['buildings'])+' sources. Backend: '+str(summary.get('mesher'))+'.</p><details><summary>Complete measurements and fallback records</summary><pre>'+html.escape(json.dumps(summary, indent=2))+'</pre></details>'+''.join(city_sections)+extra_section+''.join(sections)
    status_path = output/'run-status.json'
    if status_path.exists():
        run = json.loads(status_path.read_text())
        report = report.replace('<h1>Integrated roof catalogue — complete surfaces</h1>',
            '<h1>Integrated roof catalogue — complete surfaces</h1>'
            + '<p>Supervised run: '+html.escape(run['status'])
            + f"; {run['seconds']:.2f} s; peak RSS {run['peak_rss_bytes']/1024**2:.1f} MiB. "
            + 'Budget: 600 s wall time, 2 GiB RSS. <a href="run-status.json">Resource record</a>.</p>')
    other = '../catalogue-integration-triangle/report.html' if summary.get('mesher')=='dtcc_mesher' else '../catalogue-integration/report.html'
    report = report.replace('<h1>Integrated roof catalogue — complete surfaces</h1>',
        '<h1>Integrated roof catalogue — complete surfaces</h1>'
        + f'<p><a href="{other}">Comparison report for the other backend</a> · <a href="results.json">All 341 per-source measurements and final outcomes</a></p>')
    (output/'report.html').write_text(report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--mesher', choices=['dtcc_mesher','triangle'], default='dtcc_mesher')
    parser.add_argument('--wall-baseline', type=Path, help='Preserved wall-comparison baseline root; run old/new/LOD1 evaluation')
    args = parser.parse_args()
    if args.wall_baseline is not None:
        from roof_wall_evaluation import run
        run(args)
        return
    args.output.mkdir(parents=True, exist_ok=True)
    from dtcc_core.builder import _dtcc_builder
    import scipy
    import dtcc_mesher
    repository = Path(__file__).resolve().parents[1]
    sources = list((repository/'dtcc_core/builder/roofs').rglob('*.py'))
    sources += [repository/'dtcc_core/builder/meshing'/name for name in (
        'catalogue_surface.py', 'roof_regions.py', 'roof_interfaces.py', 'roof_validation.py', 'roof_coverage.py', 'flat_mesh_backends.py', 'dtcc_mesher_backend.py')]
    sources += [repository/'dtcc_core/cpp/include/MeshBuilder.h',
                repository/'dtcc_core/cpp/dtcc_builder.cpp',
                repository/'dtcc_core/builder/geometry_builders/meshes.py']
    dump(dict(input=str(args.input.resolve()), cases=str(args.cases.resolve()),
        python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
        shapely=shapely.__version__, dtcc_mesher=getattr(dtcc_mesher, '__version__', None),
        native_extension=str(_dtcc_builder.__file__),
        native_sha256=hashlib.sha256(Path(_dtcc_builder.__file__).read_bytes()).hexdigest(),
        production_sha256={str(p.relative_to(repository)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        reconstruction_options=dict(lod='lod2', max_mesh_size=8.,
            building_mesh_triangle_size=4., min_building_area=15., min_building_detail=.5, mesher=args.mesher),
        holdout_used_for_fitting=False), args.output/'experiment.json')
    cases = json.loads(args.cases.read_text())
    dump(cases, args.output/'cases.json')
    originals, rasters, area_meshes, rows, areas = {}, {}, {}, [], {}
    for area in ['gothenburg', 'lund']:
        originals[area] = pickle.load(open(args.input/area/'buildings.pkl', 'rb'))
        manifest = json.loads((args.input/area/'manifest.json').read_text())
        rasters[area] = terrain_raster(manifest, args.output)
        for index, b in enumerate(originals[area]):
            cloud = np.load(args.input/area/f'b{index:03d}'/'points.npz')
            b.add_geometry(PointCloud(points=cloud['train']), GeometryType.POINT_CLOUD)
            fitted = fit_roof(b.lod0.to_polygon(), cloud['train'], b.attributes['ground_height'], source_key=b.id)
            error = abs(fitted.model.predict(cloud['test'][:, :2])-cloud['test'][:, 2]) if fitted.model and len(cloud['test']) else np.array([])
            rows.append(dict(area=area, index=index, id=b.id, fit=fitted.diagnostics,
                shape=fitted.model.shape if fitted.model else None,
                heldout_provisional=dict(n=len(error), median_ae=float(np.median(error)) if len(error) else None,
                    p90_ae=float(np.quantile(error, .9)) if len(error) else None)))
        city = city_of(originals[area], rasters[area])
        mesh, result = build(city, catalogue=True, mesher=args.mesher)
        dump(result, args.output/f'{area}-audit.json')
        if mesh:
            mesh.save(args.output/f'{area}-catalogue.vtu')
        area_meshes[area] = mesh
        areas[area] = dict(error=result['error'], seconds=result['seconds'],
                          catalogue=result['audit'].get('roof_catalogue', {}))
        print(area, 'city complete', result['error'], flush=True)
        mapping = result['audit'].get('roof_catalogue', {}).get('patch_markers', [])
        for row in [r for r in rows if r['area']==area]:
            pts = np.load(args.input/area/f"b{row['index']:03d}"/'points.npz')['test']
            row['assembled_sources'] = [s for s in result['audit'].get('roof_catalogue',{}).get('sources',[]) if s.get('source_id')==row['id']]
            row['final_outcome'] = (row['assembled_sources'][0].get('final_outcome')
                if row['assembled_sources'] else 'no_retained_cleaned_source_cell')
            row['final_fallback_reasons'] = [reason for source in row['assembled_sources']
                for reason in source['diagnostics'].get('rejection_reasons', [])]
            row['final_fallback_reasons'] += [fallback['reason'] for source in row['assembled_sources']
                for fallback in source.get('block_fallbacks', [])]
            if not row['assembled_sources']:
                row['final_fallback_reasons'] = ['Source is absent from the retained cleaned source cells; no correspondence was invented']
            row['assembled_holdout'] = mesh_score(mesh, pts, len(mapping), 2*len(mapping)) if mapping else None
        dump(rows, args.output/'results.json')
    gallery = []
    for number, case in enumerate(cases):
        area = case['area']
        indices = case.get('indices', [case.get('index')])
        buildings = [originals[area][i] for i in indices]
        bounds = shapely.union_all([b.lod0.to_polygon() for b in buildings]).bounds
        extent = box(*bounds).buffer(6).intersection(box(*meshes._raster_bounds_tuple(rasters[area]))).bounds
        raster = crop_raster(rasters[area], extent)
        city = city_of(buildings, raster)
        flat, flat_result = build(city, catalogue=False, mesher=args.mesher)
        cat, cat_result = build(city, catalogue=True, mesher=args.mesher)
        dump(dict(case=case, lod1=flat_result, catalogue=cat_result), args.output/f'case-{number:02d}.json')
        for name, mesh in [('lod1', flat), ('catalogue', cat)]:
            if mesh:
                mesh.save(args.output/f'case-{number:02d}-{name}.vtu')
        points = np.vstack([np.load(args.input/area/f'b{i:03d}'/'points.npz')['test'] for i in indices])
        diagnostics = cat_result['audit'].get('roof_catalogue', {})
        patch_count = len(diagnostics.get('patch_markers', []))
        case_scores = []
        flat_count = (int(np.max(flat.markers))+1)//2 if flat is not None else 0
        for i, source in zip(indices, buildings):
            heldout = np.load(args.input/area/f'b{i:03d}'/'points.npz')['test']
            case_scores.append(dict(index=i,
                catalogue=mesh_score(cat, heldout, patch_count, 2*patch_count) if patch_count else None,
                lod1=mesh_score(flat, heldout, flat_count, 2*flat_count) if flat_count else None))
        origin = np.asarray(extent[:2])
        fig = plt.figure(figsize=(17, 6))
        render(fig.add_subplot(131, projection='3d'), flat, np.empty((0,3)), 'Default LOD1: assembled roofs, walls, terrain', origin)
        render(fig.add_subplot(132, projection='3d'), cat, np.empty((0,3)), 'Catalogue: assembled shell' + (' (FAILED)' if cat_result['error'] else ''), origin)
        render(fig.add_subplot(133, projection='3d'), cat, points, 'Catalogue + outer holdout LiDAR', origin)
        fit_labels = ', '.join(f"{r['index']}: {r['shape'] or 'sparse'} "
            + ('accepted' if r['fit']['accepted'] else 'rejected → LOD1')
            for r in rows if r['area']==area and r['index'] in indices)
        shell_info = diagnostics.get('assemblies', [{}])[-1].get('shell', {})
        fig.suptitle(f"{area} {indices}: {case['reason']}\n{fit_labels}; "
                     f"topology valid={shell_info.get('valid', False)}, quality warning={shell_info.get('quality_warning',False)}, block fallbacks={len(diagnostics.get('block_fallbacks', []))}", fontsize=10)
        # Match view scales so apparent differences come from geometry.
        axes = fig.axes
        xlim, ylim = axes[0].get_xlim(), axes[0].get_ylim()
        zlo = min(ax.get_zlim()[0] for ax in axes)
        zhi = max(ax.get_zlim()[1] for ax in axes)
        for ax in axes:
            ax.set_xlim(xlim); ax.set_ylim(ylim); ax.set_zlim(zlo, zhi)
        fig.tight_layout()
        image = f'case-{number:02d}.png'
        fig.savefig(args.output/image, dpi=140)
        plt.close(fig)
        fitting = [r for r in rows if r['area']==area and r['index'] in indices]
        diagnostics = cat_result['audit'].get('roof_catalogue', {})
        final = diagnostics.get('assemblies', [{}])[-1]
        gallery.append(dict(case=case, image=image, error=cat_result['error'],
            fits=fitting, case_holdout=case_scores, fallback_blocks=diagnostics.get('block_fallbacks', []),
            shell=final.get('shell'), lod1_error=flat_result['error']))
        print('case', number, area, indices, cat_result['error'], flush=True)
    summary = dict(mesher=args.mesher, buildings=len(rows), provisional_fits=sum(r['shape'] is not None for r in rows), accepted_fits=sum(r['fit']['accepted'] for r in rows),
        provisional_shapes=dict(Counter(r['shape'] for r in rows)),
        complete_cities={a:dict(error=r['error'], seconds=r['seconds'],
            block_fallbacks=len(r['catalogue'].get('block_fallbacks', [])),
            retained_fits=r['catalogue'].get('retained_fits',0),
            retained_nonflat_fits=r['catalogue'].get('retained_nonflat_fits',0),
            individual_fallbacks=sum(s['outcome']=='lod1_fallback' for s in r['catalogue'].get('sources',[])),
            fallback_reasons=r['catalogue'].get('block_fallbacks',[]),
            final_shell=r['catalogue'].get('assemblies', [{}])[-1].get('shell')) for a, r in areas.items()},
        cases=len(cases), failed_cases=sum(g['error'] is not None for g in gallery),
        evaluation_holdout_used_for_fitting=False,
        terrain='Nearest classified local ground on 4 m grids; case views resample this raster')
    dump(summary, args.output/'summary.json')
    dump(gallery, args.output/'gallery.json')
    write_report(args.output, summary, gallery, areas)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
