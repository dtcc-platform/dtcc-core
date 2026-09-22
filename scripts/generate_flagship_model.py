"""Build the reproducible Delft flagship development model (no network access).

Run from the checkout:
    python scripts/generate_flagship_model.py
Or from scripts/: python generate_flagship_model.py
Default source and output paths are relative to the checkout, not the working directory.
See docs/flagship-model.md for the pinned source, attribution and Python recipes.
"""

import argparse
from collections import Counter
import copy
import hashlib
import gzip
import json
from pathlib import Path
from time import perf_counter
from datetime import datetime, timedelta, timezone

import numpy as np

from dtcc_core import builder, io
from dtcc_core.datasets import load_model_package
from dtcc_core.datasets.schema import DatasetContext
from dtcc_core.model import (
    Bounds, Building, City, Field, Grid, LineString,
    Mesh, MultiLineString, MultiSurface, Object, Point, PointCloud, Raster,
    SemanticRegion, SensorCollection, Solid, Surface, Terrain,
    VolumeGrid, VolumeMesh, exchange,
)

NS = 'https://github.com/dtcc-platform/dtcc-core/schemas/model#'
ANCHOR_ID = 'NL.IMBAG.Pand.0503100000000030'
CRS = 'EPSG:7415'  # Amersfoort / RD New + NAP height, metre axes.
NAP = 'https://www.opengis.net/def/crs/EPSG/0/5709'
SYNTHETIC = 'Synthetic DTCC development data; not a survey, observation or solver result.'
CREDIT = '© 3DBAG by tudelft3d and 3DGI'
VERSION = '3.0.0'
TIMES = (0., 60., 120., 180.)
PROFILES = {'standard': (8., 25.), 'stress': (8., 20.)}  # district surface / volume spacing
DOMAIN = tuple(json.loads(Path(__file__).with_name('flagship-sources.json').read_text())['domain_rd'])
LANDMARKS = {'Nieuwe Kerk': (84493.542, 447578.225),
             'Oude Kerk': (84215., 447622.999),
             'EWI tower': (85397.406, 446091.568)}


def timestamp(seconds):
    return (datetime(2026, 9, 13, 12, tzinfo=timezone.utc) +
            timedelta(seconds=seconds)).isoformat().replace('+00:00', 'Z')


class UrbanScene:
    """One synthetic terrain/field authority, in RD New metres and NAP heights.

    Footprint prisms approximate obstacles; roof details are not CFD boundaries.
    The smooth wake/vortex/plume construction is illustrative, not a PDE solution.
    All consumers sample this same function, including sensors and streamlines.
    """

    def __init__(self, city, bounds):
        from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
        from shapely import Polygon, STRtree, box
        self.bounds = bounds
        self.origin = np.asarray([bounds.xmin, bounds.ymin, 0.])
        self.obstacles = []
        controls = []
        for building in city.buildings:
            footprint = building.get_geometry(lod='0')
            if footprint is None:
                continue
            for surface in footprint.surfaces:
                polygon = Polygon(surface.vertices[:, :2], [h[:, :2] for h in surface.holes])
                if not polygon.is_valid or polygon.is_empty:
                    raise ValueError(f'Invalid obstacle footprint: {building.id}')
                floor = float(np.mean(surface.vertices[:, 2]))
                self.obstacles.append((polygon, floor, float(building.bounds.zmax)))
                controls.extend(np.column_stack([surface.vertices[:, :2],
                                                  np.full(len(surface.vertices), floor)]))
        self.base = float(np.median([o[1] for o in self.obstacles]))
        controls.extend([[xx, yy, self.base] for xx in (bounds.xmin, bounds.xmax)
                         for yy in (bounds.ymin, bounds.ymax)])
        controls = np.asarray(controls)
        _, indices = np.unique(controls[:, :2], axis=0, return_index=True)
        controls = controls[indices]
        self._linear_ground = LinearNDInterpolator(controls[:, :2], controls[:, 2])
        self._nearest_ground = NearestNDInterpolator(controls[:, :2], controls[:, 2])
        self._footprints = STRtree([o[0] for o in self.obstacles])
        self._floors = np.array([o[1] for o in self.obstacles])
        self._roofs = np.array([o[2] for o in self.obstacles])
        self._centres = np.array([o[0].centroid.coords[0] for o in self.obstacles])
        self._radii = np.array([max(3., np.sqrt(o[0].area)/2) for o in self.obstacles])
        # Finite 8-radius support bounds wake work to nearby buildings. The
        # smooth taper makes this explicit approximation continuous at its edge.
        self._wake_tree = STRtree([box(*(c-8*r), *(c+8*r))
                                  for c, r in zip(self._centres, self._radii)])
        water = []
        for obj in objects(city):
            if obj.semantic_type == NS+'WaterBody':
                for surface in obj.get_geometry(id='water_surface').surfaces:
                    water.append(Polygon(surface.vertices[:, :2], [h[:, :2] for h in surface.holes]))
        self._water = STRtree(water)

    def water_mask(self, xy):
        from shapely import points
        mask = np.zeros(len(xy), dtype=bool)
        indices, _ = self._water.query(points(xy[:, :2]), predicate='intersects')
        mask[indices] = True
        return mask

    def ground(self, xy):
        from shapely import points
        xy = np.asarray(xy)
        z = self._linear_ground(xy[:, :2])
        missing = ~np.isfinite(z)
        z[missing] = self._nearest_ground(xy[missing, :2])
        z[self.water_mask(xy)] = -.6  # Authored canal bed, not a measured depth.
        samples, obstacles = self._footprints.query(points(xy[:, :2]), predicate='intersects')
        z[samples] = self._floors[obstacles]
        return z

    def air(self, xyz, ground=None):
        from shapely import points
        if ground is None:
            ground = self.ground(xyz[:, :2])
        valid = xyz[:, 2] > ground + .05
        samples, obstacles = self._footprints.query(points(xyz[:, :2]), predicate='intersects')
        inside = (xyz[samples, 2] >= self._floors[obstacles]) & (xyz[samples, 2] <= self._roofs[obstacles])
        valid[samples[inside]] = False
        return valid

    def sample(self, xyz, seconds=0.):
        """Return scalar/vector arrays; invalid air samples are NaN with a mask.

        Vector components are projected easting, northing, up (not true north).
        Work is chunked to bound temporary arrays on the stress configuration.
        """
        if len(xyz) > 16384:
            chunks = [self.sample(xyz[i:i+16384], seconds) for i in range(0, len(xyz), 16384)]
            return {key: np.concatenate([c[key] for c in chunks]) for key in chunks[0]}
        phase = 2*np.pi*seconds/240.
        angle = .25 + .35*np.sin(phase)
        direction = np.array([np.cos(angle), np.sin(angle)])
        cross = np.array([-direction[1], direction[0]])
        local = xyz - self.origin
        ground = self.ground(xyz[:, :2])
        height = np.maximum(xyz[:, 2] - ground, 0.)
        speed = 2.5 + 1.4*np.log1p(height/3.)
        velocity = np.column_stack([speed[:, None]*direction, np.zeros(len(xyz))])
        shelter = np.zeros(len(xyz))
        from shapely import points
        samples, obstacles = self._wake_tree.query(points(xyz[:, :2]))
        relative = xyz[samples, :2] - self._centres[obstacles]
        along, across = relative @ direction, relative @ cross
        radius = self._radii[obstacles]
        h = np.maximum(self._roofs[obstacles]-self._floors[obstacles], 1.)
        taper = np.maximum(0., 1.-(np.linalg.norm(relative, axis=1)/(8*radius))**2)**2
        wake = taper*np.exp(-((along-1.8*radius)/(2.4*radius))**2 -
                    (across/(radius+.15*np.maximum(along, 0)))**2 - (height[samples]/(1.5*h))**2)
        np.maximum.at(shelter, samples, wake)
        # Use the strongest local shelter, so dense blocks do not accumulate
        # unrealistically large reverse velocities from overlapping prisms.
        velocity[:, :2] *= (1-.8*shelter[:, None])
        turn = .35*wake*np.sin(across/radius+phase)*np.cos(along/(2*radius))
        np.add.at(velocity[:, 0], samples, turn*cross[0])
        np.add.at(velocity[:, 1], samples, turn*cross[1])
        np.add.at(velocity[:, 2], samples, .12*wake*np.sin(along/radius)*np.exp(-height[samples]/(2*h)))
        # A broad, moving recirculation cell adds signed components away from walls.
        cx = .65*self.bounds.width + 35*np.sin(phase)
        cy = .55*self.bounds.height + 25*np.cos(phase)
        dx, dy = local[:, 0]-cx, local[:, 1]-cy
        vortex = 5*np.exp(-(dx*dx+dy*dy)/10000 - height/70)
        velocity[:, 0] -= vortex*dy/45
        velocity[:, 1] += vortex*dx/45
        velocity[:, 2] += .6*np.sin(local[:, 0]/35+phase)*np.cos(local[:, 1]/45)*np.exp(-height/40)
        source = self.origin[:2] + [.25*self.bounds.width, .45*self.bounds.height]
        relative = xyz[:, :2]-source
        along, across = relative @ direction, relative @ cross
        spread = 20 + .15*np.maximum(along, 0)
        pulse = 80 + (seconds % 240)*3.
        concentration = 90*np.exp(-((along-pulse)/180)**2 - (across/spread)**2 -
                                   ((height-5-.045*np.maximum(along, 0))/9)**2)
        temperature = 291 + 4*np.exp(-height/18)*(1-shelter) + 1.5*np.sin(phase+local[:, 0]/90)
        magnitude = np.linalg.norm(velocity, axis=1)
        pressure = 101325 - 12*height + 5*(speed**2-magnitude**2)
        temperature -= 2.5*self.water_mask(xyz[:, :2])*np.exp(-height/20)
        valid = self.air(xyz, ground)
        result = {'velocity': velocity.astype('float32'), 'speed': magnitude.astype('float32'),
                  'air_temperature': temperature.astype('float32'), 'pressure': pressure,
                  'tracer_concentration': concentration.astype('float32')}
        for values in result.values():
            values[~valid] = np.nan
        result['air_mask'] = valid
        result['flow_zone'] = np.where(~valid, 0, np.where(shelter > .35, 2, 1)).astype('uint8')
        return result


FIELD_UNITS = {'velocity': 'm/s', 'speed': 'm/s', 'air_temperature': 'K',
               'pressure': 'Pa', 'tracer_concentration': 'ug/m3', 'air_mask': '1', 'flow_zone': '1'}


def sampled_fields(scene, xyz, seconds=0., association='sample'):
    return [field(name, values, association, FIELD_UNITS[name], 3 if name == 'velocity' else 1)
            for name, values in scene.sample(xyz, seconds).items()]


def objects(root):
    yield root
    for group in root.children.values():
        for child in group:
            yield from objects(child)


def geometries(geometry):
    """Include shape children, which can also own fields and transforms."""
    yield geometry
    for child in getattr(geometry, 'surfaces', []):
        yield from geometries(child)
    for child in getattr(geometry, 'linestrings', []):
        yield from geometries(child)


def attach(owner, geometry, id, *, lod=None, role=None):
    for child in geometries(geometry):
        if isinstance(child, Raster):
            child.crs = CRS
        elif hasattr(child, 'transform'):
            child.transform.srs = CRS
    owner.add_geometry(geometry, id=id, lod=lod, role=role)
    return geometry


def feature(city, id, kind=None, cls=Object, **attributes):
    obj = cls(id=id, attributes={'name': id, 'synthetic': True,
                                'source': SYNTHETIC, **attributes})
    if kind:
        obj.semantic_type = NS + kind
    obj.transform.srs = CRS
    city.add_child(obj)
    return obj


def field(name, values, association, unit, dim=1):
    return Field(name=name, values=np.asarray(values), association=association,
                 unit=unit, dim=dim, description=SYNTHETIC)


def region(kind, indices, *, parent=None, **attributes):
    return SemanticRegion(semantic_type=NS + kind, id=kind.lower(),
                          indices=np.asarray(indices, dtype=np.int64),
                          parent=parent, attributes=attributes)


def rectangle(x, y, z, width, depth):
    return Surface(vertices=np.array([[x, y, z], [x+width, y, z],
                                     [x+width, y+depth, z], [x, y+depth, z]]))


def box(x, y, z, width, depth, height):
    v = np.array([[x, y, z], [x+width, y, z], [x+width, y+depth, z],
                  [x, y+depth, z], [x, y, z+height], [x+width, y, z+height],
                  [x+width, y+depth, z+height], [x, y+depth, z+height]])
    # Outward-oriented bottom, roof, front, right, back, left.
    rings = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4],
             [1, 2, 6, 5], [2, 3, 7, 6], [3, 0, 4, 7]]
    return Solid(surfaces=[Surface(vertices=v[r].copy()) for r in rings],
                 shells=[np.arange(6, dtype=np.int64)])


def horizontal_samples(bounds, spacing):
    """Cell-centre XY samples in row-major (northing, easting) order."""
    nx, ny = int(np.ceil(bounds.width/spacing)), int(np.ceil(bounds.height/spacing))
    xx = bounds.xmin + (np.arange(nx)+.5)*bounds.width/nx
    yy = bounds.ymin + (np.arange(ny)+.5)*bounds.height/ny
    X, Y = np.meshgrid(xx, yy)
    return np.column_stack([X.ravel(), Y.ravel()]), (ny, nx)


def add_streamlines(frame, scene, seconds):
    """Integrate short steady streamlines through the same snapshot velocity."""
    b = scene.bounds
    seeds = np.array([[b.xmin+8, yy, zz] for yy in np.linspace(b.ymin+10, b.ymax-10, 24)
                      for zz in (scene.base+5, scene.base+22)])
    lines = [[] for _ in seeds]
    active = np.ones(len(seeds), dtype=bool)
    points = seeds.copy()
    for _ in range(260):
        ids = np.flatnonzero(active)
        if not len(ids):
            break
        velocity = scene.sample(points[ids], seconds)['velocity']
        for i, v in zip(ids, velocity):
            if not np.all(np.isfinite(v)):
                active[i] = False
                continue
            lines[i].append(points[i].copy())
            points[i] += .6*v
            active[i] = (b.xmin <= points[i, 0] <= b.xmax and b.ymin <= points[i, 1] <= b.ymax
                         and b.zmin <= points[i, 2] <= b.zmax)
    paths = []
    for points in lines:
        if len(points) < 2:
            continue
        points = np.asarray(points)
        line = LineString(vertices=points)
        values = scene.sample(points, seconds)
        line.fields = [field(name, values[name], 'vertex', FIELD_UNITS[name])
                       for name in ('speed', 'tracer_concentration')]
        paths.append(line)
    attach(frame, MultiLineString(linestrings=paths), 'streamlines', role='steady_streamlines')


def add_numerics(city, scene, surface_spacing, volume_spacing):
    b = scene.bounds
    terrain = feature(city, 'synthetic-terrain', cls=Terrain,
                      description='Interpolated building-base elevations and authored canal bed; synthetic terrain')
    xy, shape = horizontal_samples(b, surface_spacing)
    ground = scene.ground(xy)
    vertices = np.column_stack([xy, ground])
    ny, nx = shape
    j, i = np.meshgrid(np.arange(ny-1), np.arange(nx-1), indexing='ij')
    k = (j*nx+i).ravel()
    faces = np.concatenate([np.column_stack([k, k+1, k+nx+1]),
                            np.column_stack([k, k+nx+1, k+nx])]).astype('int32')
    mesh = Mesh(vertices=vertices, faces=faces)
    temperature = scene.sample(vertices+[0, 0, 2])['air_temperature']
    mesh.fields = [field('ground_temperature', temperature, 'vertex', 'K')]
    attach(terrain, mesh, 'tin', lod='1', role='synthetic_terrain')
    # A complete coarse triangulation is preferable to dropping arbitrary fine
    # triangles in the overview. The full-resolution TIN remains authoritative.
    stride = max(1, int(np.ceil(25/surface_spacing)))
    rows = np.unique(np.r_[np.arange(0, ny, stride), ny-1])
    columns = np.unique(np.r_[np.arange(0, nx, stride), nx-1])
    coarse = vertices.reshape(ny, nx, 3)[rows[:, None], columns].reshape(-1, 3)
    cj, ci = np.meshgrid(np.arange(len(rows)-1), np.arange(len(columns)-1), indexing='ij')
    ck, cn = (cj*len(columns)+ci).ravel(), len(columns)
    coarse_faces = np.concatenate([np.column_stack([ck, ck+1, ck+cn+1]),
                                   np.column_stack([ck, ck+cn+1, ck+cn])]).astype('int32')
    attach(terrain, Mesh(vertices=coarse, faces=coarse_faces), 'overview', lod='0', role='preview_terrain')
    cloud = PointCloud(points=vertices.copy(), classification=np.full(len(vertices), 2, dtype='uint8'),
                       intensity=(np.arange(len(vertices)) % 65536).astype('uint16'),
                       return_number=np.ones(len(vertices), dtype='uint8'),
                       num_returns=np.ones(len(vertices), dtype='uint8'))
    cloud.classification.reshape(shape)[-3:, :3] = 9
    cloud.fields = [field('confidence', (.8+.15*np.cos((xy[:, 0]-b.xmin)/60)).astype('float32'), 'sample', '1')]
    attach(terrain, cloud, 'ground_samples', role='synthetic_samples')
    # Samples lie exactly at raster cell centres. The normal DEM builder records
    # its method and preserves the deliberately withheld 3x3 corner as no-data.
    dem = builder.build_terrain_dem(cloud, surface_spacing, unit='m', vertical_reference=NAP,
        bounds=Bounds(xmin=b.xmin, ymin=b.ymin, xmax=b.xmax, ymax=b.ymax),
        window_size=0, radius=.1, hole_fill='none', source=SYNTHETIC)
    attach(terrain, dem.get_geometry(id='dem'), 'dem', role='elevation')
    terrain.attributes['elevation_rasters'] = dem.attributes['elevation_rasters']

    simulation = feature(city, 'synthetic-flow-domain', solver=None,
        description='Building-aware synthetic wind, thermal and tracer scene',
        timestamp=timestamp(0), time_seconds=0., period_seconds=240.,
        vector_basis=['projected_easting', 'projected_northing', 'up'],
        flow_zone_labels={'0': 'solid or ground', '1': 'open air', '2': 'wake'},
        construction='UrbanScene.sample: analytic shear, local footprint wakes, vortex and advecting Gaussian pulse',
        tetrahedralization='Six tetrahedra per Cartesian grid box, shared vertices and consistent body diagonal; no boundary fitting',
        limitations=['Not a fluid solver; no conservation or no-slip guarantee',
                     'Footprint prisms approximate roofs; volume cells are not boundary conforming',
                     'NaN fields and air_mask distinguish solid/ground samples from valid zero values'],
        snapshot_ids=[f'synthetic-flow-t{int(t):03d}' for t in TIMES])
    grid = Grid(width=nx, height=ny)
    grid.bounds = Bounds(xmin=b.xmin, ymin=b.ymin, zmin=scene.base,
                         xmax=b.xmax, ymax=b.ymax, zmax=scene.base)
    gy, gx = np.gradient(ground.reshape(shape), b.height/ny, b.width/nx)
    runoff = np.maximum(.0, 2-3*(ground-ground.min()))
    grid.fields = [field('surface_runoff', runoff.astype('float32'), 'cell', 'mm/h'),
                   field('runoff_direction', np.column_stack([-gx.ravel(), -gy.ravel()]).astype('float32'),
                         'cell', '1', 2)]
    attach(simulation, grid, 'runoff_grid', role='analysis_grid')

    nxv, nyv, nzv = [int(np.ceil(size/volume_spacing)) for size in (b.width, b.height, b.depth)]
    volume = VolumeGrid(width=nxv, height=nyv, depth=nzv)
    volume.bounds = copy.deepcopy(b)
    # Match VolumeGrid.coordinates(): (northing, easting, elevation) array order.
    xv = np.linspace(b.xmin, b.xmax, nxv+1)
    yv = np.linspace(b.ymin, b.ymax, nyv+1)
    zv = np.linspace(b.zmin, b.zmax, nzv+1)
    X, Y, Z = np.meshgrid((xv[:-1]+xv[1:])/2, (yv[:-1]+yv[1:])/2, (zv[:-1]+zv[1:])/2)
    centres = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    volume.fields = sampled_fields(scene, centres, association='cell')
    attach(simulation, volume, 'air_grid', role='simulation_volume')

    # Same domain as the structured grid; six positively oriented tetrahedra/cube.
    # Store only t=0 in 3D. Time snapshots use dense slices, avoiding duplicated volumes.
    X, Y, Z = np.meshgrid(xv, yv, zv, indexing='ij')
    xyz = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    I, J, K = np.meshgrid(np.arange(nxv), np.arange(nyv), np.arange(nzv), indexing='ij')
    a = ((I*(nyv+1)+J)*(nzv+1)+K).ravel()
    dx, dy = (nyv+1)*(nzv+1), nzv+1
    corners = a[:, None]+np.array([0, dx, dy, dx+dy, 1, dx+1, dy+1, dx+dy+1])
    pattern = np.array([[0, 1, 3, 7], [0, 3, 2, 7], [0, 2, 6, 7],
                        [0, 6, 4, 7], [0, 4, 5, 7], [0, 5, 1, 7]])
    cells = corners[:, pattern].reshape(-1, 4).astype('int32')
    values = scene.sample(xyz)
    tetra = VolumeMesh(vertices=xyz, cells=cells)
    tetra.fields = [field(name, values[name], 'vertex', FIELD_UNITS[name], 3 if name == 'velocity' else 1)
                    for name in ('velocity', 'pressure', 'air_mask')]
    # A six-component symmetric velocity outer product exercises tensor-shaped
    # data with explicit ordering. It is not labelled as a physical stress tensor.
    u, v, w = values['velocity'].T
    tetra.fields.append(field('velocity_dyadic', np.column_stack([u*u, v*v, w*w, u*v, u*w, v*w]),
                              'vertex', 'm2/s2', 6))
    simulation.attributes['field_components'] = {'velocity_dyadic': ['xx', 'yy', 'zz', 'xy', 'xz', 'yz'],
                                                'runoff_direction': ['easting', 'northing']}
    cell_values = scene.sample(xyz[cells].mean(axis=1))
    tetra.markers = cell_values['flow_zone'].astype('int32')
    tetra.fields.extend([field('tracer_concentration', cell_values['tracer_concentration'], 'cell', 'ug/m3'),
                         field('cell_air_mask', cell_values['air_mask'], 'cell', '1')])
    attach(simulation, tetra, 'tetrahedra', role='simulation_volume')

    # Ordinary native Objects carry explicit snapshot metadata. No private wire
    # extensions or implicit tensor/time axes are introduced.
    candidates = np.flatnonzero(scene.air(np.column_stack([xy, ground+2])))
    station_indices = candidates[np.linspace(0, len(candidates)-1, 8, dtype=int)]
    for seconds in TIMES:
        frame = feature(simulation, f'synthetic-flow-t{int(seconds):03d}',
                        time_seconds=seconds, timestamp=timestamp(seconds),
                        scenario='changing_wind', period_seconds=240.,
                        slice_heights_agl_m=[2., 15.], sample_shape=list(shape))
        frame.relations = {'samples': [simulation.id]}
        for height in (2., 15.):
            samples = PointCloud(points=np.column_stack([xy, ground+height]))
            samples.fields = sampled_fields(scene, samples.points, seconds)
            attach(frame, samples, f'height_{int(height)}m', role='terrain_following_samples')
        focus_bounds = Bounds(**city.attributes['flagship_focus_bounds'])
        fine_spacing = 2.
        fine_xy, fine_shape = horizontal_samples(focus_bounds, fine_spacing)
        fine_ground = scene.ground(fine_xy)
        fine = PointCloud(points=np.column_stack([fine_xy, fine_ground+2.]))
        fine.fields = sampled_fields(scene, fine.points, seconds)
        attach(frame, fine, 'focus_2m', role='fine_terrain_following_samples')
        frame.attributes.update(focus_sample_shape=list(fine_shape), focus_spacing_m=fine_spacing)
        # An actual vertical plane, with independent spatial resolution.
        z = np.arange(b.zmin+surface_spacing/2, b.zmax, surface_spacing)
        x = xy[:shape[1], 0]
        XX, ZZ = np.meshgrid(x, z)
        section = PointCloud(points=np.column_stack([XX.ravel(),
                            np.full(XX.size, b.ymin+.5*b.height), ZZ.ravel()]))
        section.fields = sampled_fields(scene, section.points, seconds)
        attach(frame, section, 'vertical_section', role='field_slice')
        frame.attributes['vertical_section_shape'] = list(XX.shape)
        add_streamlines(frame, scene, seconds)
        stations = feature(frame, f'synthetic-sensors-t{int(seconds):03d}', cls=SensorCollection,
                           timestamp=timestamp(seconds), time_seconds=seconds)
        # Stable locations, selected from valid ground-level samples at t=0.
        for i, index in enumerate(station_indices):
            sensor = feature(stations, f'synthetic-sensor-{i}-t{int(seconds):03d}',
                             station_id=f'station-{i}', timestamp=timestamp(seconds))
            point = Point(x=float(xy[index, 0]), y=float(xy[index, 1]), z=float(ground[index]+2))
            point.fields = sampled_fields(scene, np.array([[point.x, point.y, point.z]]), seconds)
            attach(sensor, point, 'location')
            sensor.relations = {'samples': [frame.id]}
        frame.relations['observed_by'] = [s.id for s in stations.stations()]
    return simulation


def enrich(city, bounds, *, mesh_real=True, detail='standard'):
    """Build a coherent world-coordinate scene, retaining source representations."""
    surface_spacing, volume_spacing = PROFILES[detail]
    all_bounds = copy.deepcopy(bounds)
    all_bounds.zmin = min(bounds.zmin, -.6)-1
    all_bounds.zmax = max(bounds.zmax+30, 60.)
    scene = UrbanScene(city, all_bounds)
    add_numerics(city, scene, surface_spacing, volume_spacing)
    derived = []
    if mesh_real:
        # Dense facade/roof fields cover the reference neighbourhood.
        # Selected source representations and attributes remain untouched.
        for building in city.buildings:
            if building.id not in city.attributes.get('reference_building_ids', []):
                continue
            for part in building.building_parts:
                solid = part.get_geometry(lod='2.2')
                if solid is None:
                    continue
                mesh = solid.mesh(mesher='dtcc_mesher', triangle_size=2. if detail == 'standard' else 1.)
                triangles = mesh.vertices[mesh.faces]
                normals = np.cross(triangles[:, 1]-triangles[:, 0], triangles[:, 2]-triangles[:, 0])
                area = np.linalg.norm(normals, axis=1)/2
                if not np.all(area > 0):
                    raise ValueError(f'Degenerate derived triangle: {part.id}')
                mesh.normals = normals/(2*area[:, None])
                mesh.markers = np.full(len(mesh.faces), -1, dtype='int32')
                for i, r in enumerate(mesh.regions):
                    mesh.markers[r.indices] = i
                sun = np.array([-.4, -.3, .8660254])
                centres = triangles.mean(axis=1)
                # Footprint-prism shadow rays, independent of terrain or roof certification.
                shaded = np.zeros(len(centres), dtype=bool)
                for distance in (2., 8., 20., 45., 90.):
                    shaded |= ~scene.air(centres + distance*sun)
                irradiance = 120+780*np.maximum(mesh.normals@sun, 0)*~shaded
                mesh.fields = [Field(name='triangle_area', values=area, unit='m2', association='face',
                                     description='Derived from triangulated source geometry'),
                               field('solar_irradiance', irradiance.astype('float32'), 'face', 'W/m2'),
                               field('surface_temperature', (289+.015*irradiance).astype('float32'), 'face', 'K'),
                               field('in_shadow', shaded, 'face', '1')]
                attach(part, mesh, 'flagship_boundary_mesh', lod='2.2', role='derived_boundary_mesh')
                # The ordinary preview outlines polygons with holes. A coarse
                # triangulation displays their roofs/openings without drawing
                # every analysis triangle or simplifying away those openings.
                attach(part, solid.mesh(mesher='dtcc_mesher'), 'flagship_preview_mesh',
                       lod='2.2', role='preview_boundary_mesh')
                derived.append({'object_id': part.id, 'representation_id': 'flagship_boundary_mesh',
                                'preview_representation_id': 'flagship_preview_mesh', 'source_lod': '2.2',
                                'method': 'dtcc_mesher triangulation; approximate prism shadows'})
    city.attributes['derived_representations'] = derived
    city.attributes['flagship_sampling'] = {'detail': detail, 'surface_spacing_max_m': surface_spacing,
        'volume_spacing_max_m': volume_spacing, 'focus_surface_spacing_m': 2., 'boundary_triangle_target_m': 2. if detail == 'standard' else 1., 'times_seconds': list(TIMES),
        'time_encoding': 'Snapshot Objects with time_seconds and timestamp; not a native time-series axis',
        'coordinates': 'World coordinates, identity affines, EPSG:7415; vectors follow projected axes'}
    city.calculate_bounds()
    return city


def load_sources(directory):
    """Admit pinned tiles, preserve the old reference patch, select whole buildings.

    Complete building envelopes must lie in the square: edge buildings are not
    clipped into invalid solids. Sources are checked before any output is written.
    """
    manifest = json.loads(Path(__file__).with_name('flagship-sources.json').read_text())
    sources = {source['file']: source for source in manifest['sources']}
    source_metadata = {}

    def read_source(filename):
        source = sources[filename]
        path = directory/source['file']
        if not path.is_file():
            raise ValueError(f'Missing source {path}; run scripts/download_flagship_sources.py')
        raw = path.read_bytes()
        if hashlib.sha256(raw).hexdigest() != source['sha256']:
            raise ValueError(f'Source checksum mismatch: {path}')
        document = json.loads(gzip.decompress(raw) if path.suffix == '.gz' else raw)
        if document.get('type') == 'CityJSON':
            source_metadata[source['file']] = document.get('metadata', {}).copy()
            # The Core CityJSON profile admits CRS and extent only. Preserve
            # supplier contact/title/version metadata on the resulting City.
            document['metadata'] = {k: v for k, v in document.get('metadata', {}).items()
                                    if k in ('referenceSystem', 'geographicalExtent')}
        return document

    old = io.load_3dbag(read_source('3dbag.city.json'), extent_policy='recompute')
    anchor = next(b for b in old.buildings if b.id == ANCHOR_ID).bounds
    centre = np.array([(anchor.xmin+anchor.xmax)/2, (anchor.ymin+anchor.ymax)/2])
    def distance(b):
        bb = b.bounds
        return (float(np.sum((np.array([(bb.xmin+bb.xmax)/2, (bb.ymin+bb.ymax)/2])-centre)**2)), b.id)
    reference = sorted(old.buildings, key=distance)[:64]
    del old
    city = City(id='dtcc-flagship-delft-v3')
    city.transform.srs = CRS
    chosen = {b.id: b for b in reference}
    origin = {b.id: '3dbag.city.json' for b in reference}
    x0, y0, x1, y1 = DOMAIN
    for filename in sources:
        if not filename.endswith('.gz'):
            continue
        tile = io.load_3dbag(read_source(filename), extent_policy='recompute')
        for building in tile.buildings:
            b = building.bounds
            if x0 <= b.xmin and y0 <= b.ymin and b.xmax <= x1 and b.ymax <= y1:
                if building.id not in chosen:
                    chosen[building.id] = building
                    origin[building.id] = filename
        print(f'Admitted {filename}: {len(chosen)} buildings selected', flush=True)
    if 'NL.IMBAG.Pand.0503100000030264' not in chosen:
        raise ValueError('The EWI high-rise complex must be wholly included')
    city.add_children([chosen[id] for id in sorted(chosen)])
    city.calculate_bounds()
    if not (x0 <= city.bounds.xmin and y0 <= city.bounds.ymin and
            city.bounds.xmax <= x1 and city.bounds.ymax <= y1):
        raise ValueError('Reference buildings extend beyond the selected square')
    bounds = Bounds(xmin=x0, ymin=y0, xmax=x1, ymax=y1,
                    zmin=city.bounds.zmin, zmax=city.bounds.zmax)
    # Use the existing patch envelope, with 30 m context, for fine field samples.
    xy = np.array([[b.bounds.xmin, b.bounds.ymin, b.bounds.xmax, b.bounds.ymax] for b in reference])
    focus = {'xmin': float(xy[:, 0].min()-30), 'ymin': float(xy[:, 1].min()-30),
             'xmax': float(xy[:, 2].max()+30), 'ymax': float(xy[:, 3].max()+30)}
    city.attributes.update(reference_building_ids=sorted(b.id for b in reference),
                           building_source_files=origin, source_metadata=source_metadata, flagship_focus_bounds=focus,
                           landmarks_rd={name: list(xy) for name, xy in LANDMARKS.items()})
    return city, bounds, manifest, read_source('water.geojson')


def add_water(city, document, bounds):
    """Real BGT plan outlines; explicitly synthetic NAP elevation and temperature."""
    from shapely.geometry import shape, box
    domain = box(bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax)
    if any(link['rel'] == 'next' for link in document.get('links', [])):
        raise ValueError('Water source is incomplete (unread next page)')
    for item in document['features']:
        attrs = item['properties']
        if attrs.get('eind_registratie') or attrs.get('termination_date'):
            continue
        polygon = shape(item['geometry'])
        if not polygon.is_valid:
            raise ValueError(f'Invalid BGT water outline {item["id"]}')
        clipped = polygon.intersection(domain)
        polygons = [clipped] if clipped.geom_type == 'Polygon' else list(clipped.geoms)
        polygons = [p for p in polygons if p.geom_type == 'Polygon' and not p.is_empty and p.area > 0]
        if not polygons:
            continue
        obj = Object(id='bgt-water-'+str(item['id']), semantic_type=NS+'WaterBody',
                     attributes={'source': 'BGT / PDOK', 'source_attributes': attrs,
                                 'elevation_note': 'Synthetic water plane at -0.5 m NAP; only XY outlines are surveyed'})
        obj.transform.srs = CRS
        surface = MultiSurface(surfaces=[Surface().from_polygon(p, -.5) for p in polygons])
        surface.regions = [region('WaterSurface', np.arange(len(polygons)))]
        surface.fields = [field('water_temperature', np.full(len(polygons), 289., dtype='float32'),
                                'face', 'K')]
        attach(obj, surface, 'water_surface', lod='1', role='water_outline')
        city.add_child(obj)


def inventory(city):
    native, semantic, geometry_types, lods, fields = Counter(), Counter(), Counter(), Counter(), []
    region_types, elements = Counter(), Counter()
    object_list = list(objects(city))
    for obj in object_list:
        native[type(obj).__name__] += 1
        semantic[obj.semantic_type or type(obj).__name__] += 1
        for rep_id, rep in obj.geometry.items():
            geometry_types[type(rep.geometry).__name__] += 1
            g = rep.geometry
            for name in ('vertices', 'faces', 'cells', 'points'):
                if hasattr(g, name):
                    elements[type(g).__name__+'.'+name] += len(getattr(g, name))
            if isinstance(g, (Grid, VolumeGrid)):
                elements[type(g).__name__+'.cells'] += g.num_cells
            if isinstance(g, Raster):
                elements['Raster.pixels'] += g.data.size
            if rep.lod is not None:
                lods[rep.lod] += 1
            for geom in geometries(rep.geometry):
                region_types.update(r.semantic_type.rsplit('#', 1)[-1] for r in getattr(geom, 'regions', []))
                for f in getattr(geom, 'fields', []):
                    finite = f.values[np.isfinite(f.values)]
                    fields.append({'object_id': obj.id, 'representation_id': rep_id,
                                   'name': f.name, 'association': f.association, 'unit': f.unit,
                                   'shape': list(f.values.shape), 'dtype': str(f.values.dtype), 'dim': f.dim,
                                   'missing_components': int(f.values.size-finite.size),
                                   'component_range': [float(finite.min()), float(finite.max())] if finite.size else None})
    return {'objects': len(object_list), 'native_object_types': dict(sorted(native.items())),
            'semantic_types': dict(sorted(semantic.items())),
            'representation_types': dict(sorted(geometry_types.items())),
            'lods': dict(sorted(lods.items())), 'semantic_regions': dict(sorted(region_types.items())),
            'fields': fields, 'element_counts': dict(elements),
            'coverage': {'associations': sorted({f['association'] for f in fields}),
                         'field_dimensions': sorted({f['dim'] for f in fields}),
                         'field_dtypes': sorted({f['dtype'] for f in fields}),
                         'known_limits': ['Time and component labels use Object metadata; no native time/tensor axis',
                                          'Detached synthetic architecture and affine specimen omitted from this city scene',
                                          'FieldSlice and StreamlineCollection wrappers are outside canonical exchange',
                                          'No DeSO administrative-area example or exhaustive dtype matrix']},
            'relation_targets': sum(len(ids) for o in object_list for ids in o.relations.values())}


def preview(city, path):
    """Use the same saved-data inspection view available to a local user."""
    import matplotlib
    matplotlib.use('Agg')
    from inspect_flagship_model import dashboard
    import matplotlib.pyplot as plt
    figure = dashboard(city, interactive=False)
    figure.savefig(path, dpi=150)
    plt.close(figure)


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-dir', type=Path, default=root/'data/flagship-source')
    parser.add_argument('--output', type=Path, default=root/'data/flagship',
                        help='Output directory (default: checkout/data/flagship)')
    parser.add_argument('--detail', choices=PROFILES, default='standard',
                        help='standard: 8 m district / 25 m volume / 2 m focus; stress: 8 / 20 / 2 m, with 1 m boundary meshes')
    args = parser.parse_args()
    start = perf_counter()
    city, bounds, manifest, water = load_sources(args.source_dir)
    chosen = city.buildings
    # Keep all source LoDs in the reference patch. Elsewhere retain footprints
    # and the finest available source representation, bounding the native file
    # without simplifying roof geometry or losing the source attributes.
    reference_ids = set(city.attributes['reference_building_ids'])
    for building in chosen:
        if building.id in reference_ids:
            continue
        for obj in objects(building):
            records = [(id, rep) for id, rep in obj.geometry.items() if rep.lod != '0']
            finest = max(records, key=lambda item: float(item[1].lod or 0))[0] if records else None
            obj.geometry = {id: rep for id, rep in obj.geometry.items() if rep.lod == '0' or id == finest}
    original = {b.id: hashlib.sha256(exchange.dumps(b)).hexdigest() for b in chosen}
    city.attributes['flagship'] = {'version': VERSION, 'name': 'Delft flagship development model',
        'sources': manifest['sources'], 'source_release': manifest['release'],
        'credit': CREDIT, 'license': 'https://creativecommons.org/licenses/by/4.0/',
        'attribution_url': 'https://docs.3dbag.nl/en/copyright/',
        'domain_rd': list(DOMAIN),
        'selection': 'Whole building envelopes contained in 2 km square; original 64 buildings take precedence by BAG ID',
        'real_building_ids': sorted(b.id for b in chosen),
        'changes': 'Spatial subset; all LoDs in reference patch, finest source LoD elsewhere; real BGT water outlines; reference-patch solar meshes; district and fine-patch fields',
        'synthetic_notice': SYNTHETIC,
        'coordinates': 'EPSG:7415; metre coordinates and NAP elevations',
        'limitations': ['No full-source geometric validity certification or repair',
                       'Synthetic terrain, water elevations and fields, not observations or solver results',
                       'Reference patch uses the older pinned sample, other buildings use v20250903',
                       'Whole buildings crossing domain edges are omitted',
                       'Wind obstacles and shadows use footprint prisms; tetrahedra are not boundary conforming']}
    add_water(city, water, bounds)
    enrich(city, bounds, detail=args.detail)
    print('Fields and reference boundary meshes generated', flush=True)
    generated_seconds = perf_counter()-start
    args.output.mkdir(parents=True, exist_ok=True)
    city.dataset_context = DatasetContext(
        identity={'name': city.id, 'title': 'Delft flagship development model', 'version': VERSION},
        metadata={'crs': [CRS], 'geographic_coverage': 'Delft, Netherlands',
                  'description': 'Real buildings with explicitly synthetic urban fields'},
        presentation={'headline': 'Delft: an urban field laboratory',
                      'summary': 'Explore changing wind, tracer transport, surface heating and model semantics.',
                      'view_hints': {'default_snapshot_id': 'synthetic-flow-t000',
                                     'default_representation': 'height_2m'},
                      'warnings': [SYNTHETIC]},
        request={'dataset_name': city.id, 'parameters': city.attributes['flagship_sampling']},
        provenance={'sources': manifest['sources']+[{'description': SYNTHETIC}],
                    'processing_steps': [{'operation': 'generate_flagship_model.py', 'version': VERSION,
                                          'source_selection': sorted(b.id for b in chosen)}]})
    path = args.output/'flagship.dtcc'
    start_save = perf_counter()
    city.save(path)
    save_seconds = perf_counter()-start_save
    start_load = perf_counter()
    restored = io.load_model(path)
    load_seconds = perf_counter()-start_load
    payload = path.read_bytes()
    assert exchange.dumps(restored) == payload
    # The principal file has no mixed-frame envelopes, including after decoding.
    assert restored.bounds.width == 2000 and restored.bounds.height == 2000
    assert restored.dataset_context is None
    # Remove the derived representations in the independent decode, proving all
    # selected source/enrichment facts retained their exact native values.
    for obj in objects(restored):
        obj.geometry.pop('flagship_boundary_mesh', None)
        obj.geometry.pop('flagship_preview_mesh', None)
    for obj in objects(restored):
        if obj.id in original:
            assert hashlib.sha256(exchange.dumps(obj)).hexdigest() == original[obj.id], obj.id
    restored.buildings[0].attributes['storeys_above_ground'] = True
    try:
        restored.save(path)
    except ValueError as error:
        assert 'storeys_above_ground' in str(error)
    else:
        raise AssertionError('Invalid storey count was accepted')
    assert path.read_bytes() == payload
    del restored
    package = city.export(args.output/'flagship.dtccpkg', canonical=True)
    packaged = load_model_package(package.path)
    assert exchange.dumps(packaged) == payload
    assert packaged.dataset_context == city.dataset_context
    del packaged
    start_preview = perf_counter()
    preview(city, args.output/'preview.png')
    preview_seconds = perf_counter()-start_preview
    report = inventory(city)
    import resource
    import sys
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        peak_rss *= 1024
    report.update(schema_version=city.schema_version, wire_version=exchange.VERSION,
                  native_bytes=len(payload), native_sha256=hashlib.sha256(payload).hexdigest(),
                  real_buildings=len(chosen), original_source_preserved=True,
                  native_and_package_exact=True, rejected_save_preserved_file=True,
                  domain_rd=list(DOMAIN), reference_buildings=len(city.attributes['reference_building_ids']),
                  water_bodies=sum(o.semantic_type == NS+'WaterBody' for o in objects(city)),
                  tallest_buildings=sorted([{'id': b.id, 'roof_elevation_nap_m': b.bounds.zmax}
                                            for b in city.buildings], key=lambda b: -b['roof_elevation_nap_m'])[:10],
                  sampling=city.attributes['flagship_sampling'],
                  timings_seconds={'generation': round(generated_seconds, 3), 'save': round(save_seconds, 3),
                                   'load': round(load_seconds, 3), 'preview': round(preview_seconds, 3)},
                  peak_process_rss_bytes=peak_rss,
                  elapsed_seconds=round(perf_counter()-start, 3))
    (args.output/'inventory.json').write_text(json.dumps(report, indent=2)+'\n')
    (args.output/'README.md').write_text((Path(__file__).resolve().parents[1]/'docs/flagship-model.md').read_text())
    print(json.dumps({key: report[key] for key in ('objects', 'native_bytes', 'timings_seconds',
                                                  'peak_process_rss_bytes', 'elapsed_seconds')}, indent=2))
    print(f'Inspect: python {root/"scripts/inspect_flagship_model.py"} {path}')


if __name__ == '__main__':
    main()
