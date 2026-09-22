"""Inspect the saved flagship using Core plotting and Matplotlib (no web services).

Run without arguments for a city/field dashboard with a snapshot slider. Use
--view city for a native Core preview, fine for local wind, or tetra for a mesh cut. --save writes
a PNG without opening a window. All displayed values come from the loaded model.
"""

import argparse
import copy
from pathlib import Path

import numpy as np

from dtcc_core import io


def objects(root):
    yield root
    for group in root.children.values():
        for child in group:
            yield from objects(child)


def physical_scene(city):
    """A read-only display selection excluding numerical domains and snapshots."""
    def select(obj):
        view = copy.copy(obj)
        representation = ('overview' if obj.id == 'synthetic-terrain' else
                          'flagship_preview_mesh' if 'flagship_preview_mesh' in obj.geometry else None)
        if representation:
            view.geometry = {representation: obj.geometry[representation]}
        view.children = {key: [select(child) for child in group if child.id != 'synthetic-flow-domain']
                         for key, group in obj.children.items()}
        return view
    return select(city)


def values(geometry, name):
    return next(f.values for f in geometry.fields if f.name == name)


def credit(figure):
    color = '#526577' if sum(figure.get_facecolor()[:3]) > 1.5 else '#bacad3'
    figure.text(.99, .012, '© 3DBAG by tudelft3d and 3DGI | CC BY 4.0 | '
                'BGT / PDOK water outlines | Fields and elevations of terrain/water are synthetic',
                ha='right', fontsize=8, color=color)


def footprint_patches(city, origin, *, water=False):
    """Compound paths preserve islands/courtyards in the plan views."""
    from matplotlib.path import Path as MplPath
    from matplotlib.patches import PathPatch
    from shapely import Polygon
    from shapely.geometry.polygon import orient
    patches, heights = [], []
    owners = [o for o in objects(city) if o.semantic_type and o.semantic_type.endswith('#WaterBody')] if water else city.buildings
    for obj in owners:
        geometry = obj.get_geometry(id='water_surface') if water else obj.get_geometry(lod='0')
        if geometry is None:
            continue
        for surface in geometry.surfaces:
            polygon = orient(Polygon(surface.vertices[:, :2], [h[:, :2] for h in surface.holes]))
            rings = [np.asarray(r.coords)-origin for r in [polygon.exterior, *polygon.interiors]]
            paths = [MplPath(r, [MplPath.MOVETO]+[MplPath.LINETO]*(len(r)-2)+[MplPath.CLOSEPOLY]) for r in rings]
            patches.append(PathPatch(MplPath.make_compound_path(*paths)))
            heights.append(obj.bounds.zmax-float(np.mean(surface.vertices[:, 2])))
    return patches, heights


def city_map(city, ax):
    from matplotlib.collections import PatchCollection
    features = {o.id: o for o in objects(city)}
    b = features['synthetic-flow-domain'].get_geometry(id='air_grid').bounds
    origin = np.array([b.xmin, b.ymin])
    patches, _ = footprint_patches(city, origin, water=True)
    ax.add_collection(PatchCollection(patches, facecolor='#73bde0', edgecolor='none'))
    patches, heights = footprint_patches(city, origin)
    buildings = PatchCollection(patches, cmap='cividis', edgecolor='#273b46', linewidth=.15)
    buildings.set_array(np.asarray(heights))
    buildings.set_clim(0, 100)
    ax.add_collection(buildings)
    for name, xy in city.attributes.get('landmarks_rd', {}).items():
        xy = np.asarray(xy)-origin
        ax.annotate(name, xy, xytext=({'Oude Kerk': (-80, -15), 'Nieuwe Kerk': (22, 15), 'EWI tower': (-90, 15)}.get(name, (8, 12))),
                    textcoords='offset points', fontsize=8, weight='bold',
                    arrowprops={'arrowstyle': '-', 'color': '#263c4b'},
                    bbox={'facecolor': 'white', 'alpha': .85, 'edgecolor': 'none', 'pad': 2})
    focus = city.attributes.get('flagship_focus_bounds')
    if focus:
        from matplotlib.patches import Rectangle
        ax.add_patch(Rectangle((focus['xmin']-b.xmin, focus['ymin']-b.ymin),
                     focus['xmax']-focus['xmin'], focus['ymax']-focus['ymin'],
                     fill=False, edgecolor='#ba4357', linewidth=1, linestyle='--'))
    ax.set(xlim=(0, b.width), ylim=(0, b.height), aspect='equal',
           xlabel='Easting from domain origin [m]', ylabel='Northing from domain origin [m]',
           title=f'Delft · {len(city.buildings):,} real buildings · roof height [m]')
    ax.set_facecolor('#eef0e9')
    ax.figure.colorbar(buildings, ax=ax, fraction=.035, pad=.025)


def fine_view(city, frame):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    samples = frame.get_geometry(id='focus_2m')
    shape = tuple(frame.attributes['focus_sample_shape'])
    b = city.attributes['flagship_focus_bounds']
    origin = np.array([b['xmin'], b['ymin']])
    extent = (0, b['xmax']-b['xmin'], 0, b['ymax']-b['ymin'])
    fig, ax = plt.subplots(figsize=(12, 9))
    artist = ax.imshow(values(samples, 'speed').reshape(shape), origin='lower', extent=extent,
                       vmin=0, vmax=10, cmap='viridis')
    patches, _ = footprint_patches(city, origin)
    ax.add_collection(PatchCollection(patches, facecolor='none', edgecolor='#223545', linewidth=.5))
    stride = max(1, shape[1]//35)
    xy = samples.points[:, :2].reshape(*shape, 2)[::stride, ::stride]-origin
    velocity = values(samples, 'velocity').reshape(*shape, 3)[::stride, ::stride]
    ax.quiver(xy[..., 0], xy[..., 1], velocity[..., 0], velocity[..., 1],
              angles='xy', scale_units='xy', scale=.4, width=.002)
    ax.set(xlim=extent[:2], ylim=extent[2:], xlabel='Easting from focus origin [m]',
           ylabel='Northing from focus origin [m]', title=f'Fine reference patch · {frame.attributes["focus_spacing_m"]:g} m samples · {frame.attributes["timestamp"]}')
    fig.colorbar(artist, ax=ax, label='Wind speed [m/s] at 2 m AGL')
    credit(fig)
    return fig


def tetra_view(city):
    """Intersect stored tetrahedral connectivity; no regenerated mesh or fields."""
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    simulation = next(o for o in objects(city) if o.id == 'synthetic-flow-domain')
    mesh = simulation.get_geometry(id='tetrahedra')
    plane = (mesh.bounds.ymin+mesh.bounds.ymax)/2 + .173  # Avoid an exactly coincident grid face.
    y = mesh.vertices[mesh.cells, 1]
    selected = np.flatnonzero((y.min(axis=1) < plane) & (y.max(axis=1) > plane))
    pressure = values(mesh, 'pressure')
    tracer = values(mesh, 'tracer_concentration')
    polygons, pvalues, cvalues = [], [], []
    edges = ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3))
    for index in selected:
        ids = mesh.cells[index]
        vertices = mesh.vertices[ids]
        points, data = [], []
        for a, b in edges:
            if (vertices[a, 1] < plane) == (vertices[b, 1] < plane):
                continue
            t = (plane-vertices[a, 1])/(vertices[b, 1]-vertices[a, 1])
            point = vertices[a]+t*(vertices[b]-vertices[a])
            points.append(point[[0, 2]])
            data.append((1-t)*pressure[ids[a]]+t*pressure[ids[b]])
        if len(points) < 3:
            continue
        points = np.array(points)
        centred = points-points.mean(axis=0)
        order = np.argsort(np.arctan2(centred[:, 1], centred[:, 0]))
        points[:, 0] -= mesh.bounds.xmin
        polygons.append(points[order])
        pvalues.append(float(np.mean(data)))
        cvalues.append(tracer[index])
    fig, axes = plt.subplots(2, 1, figsize=(15, 8), constrained_layout=True)
    for ax, array, label, cmap in zip(axes, [pvalues, cvalues],
            ['Interpolated vertex pressure, averaged over cut polygon [Pa]', 'Cell tracer concentration [µg/m³]'],
            ['viridis', 'magma']):
        colors = plt.get_cmap(cmap).with_extremes(bad='#c2cad1')
        collection = PolyCollection(polygons, array=np.ma.masked_invalid(array), cmap=colors,
                                    edgecolors='#40515c', linewidths=.15)
        ax.add_collection(collection)
        ax.set(xlim=(0, mesh.bounds.width), ylim=(mesh.bounds.zmin, mesh.bounds.zmax),
               xlabel='Easting from domain origin [m]', ylabel='Elevation [m NAP]', title=label)
        fig.colorbar(collection, ax=ax, fraction=.025, pad=.015)
    fig.suptitle(f'{len(mesh.cells):,} stored tetrahedra · cut at N = {plane:.3f} m · vertical scale exaggerated\n'
                 'Cartesian mesh, not building conforming · grey = masked samples · synthetic fields')
    return fig


def dashboard(city, *, interactive=True, frame_index=0):
    """Inspect four stored snapshots with shared, fixed colour scales.

    Maps use local easting/northing offsets for legibility; the Core 3D preview
    retains the native CRS. Wind arrows show actual projected XY components.
    Heights marked AGL follow the synthetic terrain. The section uses NAP Z.
    """
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider

    features = {o.id: o for o in objects(city)}
    simulation = features['synthetic-flow-domain']
    frames = [features[id] for id in simulation.attributes['snapshot_ids']]
    bounds = simulation.get_geometry(id='air_grid').bounds
    origin = np.array([bounds.xmin, bounds.ymin])
    extent = (0, bounds.width, 0, bounds.height)
    max_speed = max(float(np.nanmax(values(f.get_geometry(id='height_2m'), 'speed'))) for f in frames)
    max_speed = np.ceil(max_speed)
    max_tracer = 10*np.ceil(max(float(np.nanmax(values(f.get_geometry(id='height_15m'),
                                              'tracer_concentration'))) for f in frames)/10)
    figure = plt.figure(figsize=(16, 11), facecolor='#f4f7fa')
    layout = figure.add_gridspec(2, 2, left=.045, right=.96, bottom=.13, top=.90,
                               wspace=.24, hspace=.29)
    scene_ax = figure.add_subplot(layout[0, 0])
    city_map(city, scene_ax)
    wind_ax = figure.add_subplot(layout[0, 1])
    tracer_ax = figure.add_subplot(layout[1, 0])
    section_ax = figure.add_subplot(layout[1, 1])

    def footprints(ax):
        from matplotlib.collections import PatchCollection
        patches, _ = footprint_patches(city, origin)
        ax.add_collection(PatchCollection(patches, facecolor='none', edgecolor='#223545',
                                         linewidth=.3, alpha=.8))

    first = frames[frame_index]
    shape = tuple(first.attributes['sample_shape'])
    wind = first.get_geometry(id='height_2m')
    tracer = first.get_geometry(id='height_15m')
    section = first.get_geometry(id='vertical_section')
    section_shape = tuple(first.attributes['vertical_section_shape'])
    section_z = section.points[:, 2].reshape(section_shape)[:, 0]
    dz = section_z[1]-section_z[0] if len(section_z) > 1 else bounds.depth
    section_extent = (0, bounds.width, section_z[0]-dz/2, section_z[-1]+dz/2)
    artists = [wind_ax.imshow(values(wind, 'speed').reshape(shape), origin='lower', extent=extent,
                             cmap='viridis', vmin=0, vmax=max_speed, interpolation='nearest'),
               tracer_ax.imshow(values(tracer, 'tracer_concentration').reshape(shape), origin='lower',
                                extent=extent, cmap='magma', vmin=0, vmax=max_tracer, interpolation='nearest'),
               section_ax.imshow(values(section, 'air_temperature').reshape(section_shape), origin='lower',
                                 extent=section_extent, cmap='inferno', vmin=289, vmax=297,
                                 interpolation='nearest', aspect='auto')]
    for ax, artist, label in zip((wind_ax, tracer_ax, section_ax), artists,
                                ('Wind speed [m/s]', 'Tracer [µg/m³]', 'Air temperature [K]')):
        figure.colorbar(artist, ax=ax, fraction=.035, pad=.025, label=label)
        ax.set_facecolor('#dbe2e8')
        ax.set_xlabel('Easting from domain origin [m]')
        ax.tick_params(labelsize=8)
    for ax in (wind_ax, tracer_ax):
        footprints(ax)
        ax.set_ylabel('Northing from domain origin [m]')
        ax.set_xlim(extent[:2]); ax.set_ylim(extent[2:])
    stride = max(1, shape[1]//28)
    xy = wind.points[:, :2].reshape(*shape, 2)[::stride, ::stride]-origin
    vectors = values(wind, 'velocity').reshape(*shape, 3)[::stride, ::stride]
    arrows = wind_ax.quiver(xy[..., 0], xy[..., 1], vectors[..., 0], vectors[..., 1],
                           angles='xy', scale_units='xy', scale=.08, width=.0025, color='#14252f')
    wind_ax.quiverkey(arrows, .88, 1.035, 5, '5 m/s', labelpos='E', coordinates='axes')
    wind_ax.set_title('Wind at 2 m AGL', fontsize=12, loc='left')
    tracer_ax.set_title('Advecting tracer at 15 m above terrain', fontsize=12, loc='left')
    section_ax.set_title('Vertical temperature section through the neighbourhood', fontsize=12, loc='left')
    section_ax.set_ylabel('Elevation [m NAP]')
    heading = figure.suptitle('', x=.045, ha='left', fontsize=19, weight='bold', color='#142d40')
    figure.text(.045, .922, 'EPSG:7415 · RD New + NAP · Grey cells are solid/ground or missing data · '
                 'Synthetic illustration, not a solver result', fontsize=10, color='#526577')

    def update(index):
        frame = frames[int(index)]
        heading.set_text(f'DTCC / DELFT FLAGSHIP     {frame.attributes["timestamp"]}')
        for artist, rep, field_name, dimensions in zip(artists,
                ('height_2m', 'height_15m', 'vertical_section'),
                ('speed', 'tracer_concentration', 'air_temperature'), (shape, shape, section_shape)):
            artist.set_data(values(frame.get_geometry(id=rep), field_name).reshape(dimensions))
        vector = values(frame.get_geometry(id='height_2m'), 'velocity').reshape(*shape, 3)[::stride, ::stride]
        arrows.set_UVC(vector[..., 0], vector[..., 1])
        figure.canvas.draw_idle()

    update(frame_index)
    if interactive:
        slider = Slider(figure.add_axes([.18, .06, .62, .022]), 'Snapshot', 0, len(frames)-1,
                        valinit=frame_index, valstep=1)
        slider.on_changed(update)
        # Matplotlib widgets need a live reference for callbacks.
        figure._flagship_slider = slider
    else:
        figure.text(.045, .055, 'Open scripts/inspect_flagship_model.py for the time slider and 3D city, fine patch and tetrahedral section views.',
                    fontsize=10, color='#526577')
    credit(figure)
    return figure


def main():
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', nargs='?', type=Path, default=root/'data/flagship/flagship.dtcc')
    parser.add_argument('--view', choices=['dashboard', 'city', 'wind', 'solar', 'dem', 'streamlines', 'map', 'fine', 'tetra'], default='dashboard')
    parser.add_argument('--frame', type=int, choices=range(4), default=0)
    parser.add_argument('--save', type=Path, help='Save a PNG using a noninteractive backend')
    args = parser.parse_args()
    if not args.model.is_file():
        parser.error(f'Model not found: {args.model}; run generate_flagship_model.py first')
    if args.save:
        import matplotlib
        matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    city = io.load_city(args.model)
    if args.view == 'dashboard':
        figure = dashboard(city, interactive=not args.save, frame_index=args.frame)
    else:
        features = {o.id: o for o in objects(city)}
        simulation = features['synthetic-flow-domain']
        frame = features[simulation.attributes['snapshot_ids'][args.frame]]
        if args.view == 'map':
            figure, ax = plt.subplots(figsize=(12, 11))
            city_map(city, ax)
        elif args.view == 'fine':
            figure = fine_view(city, frame)
            ax = figure.axes[0]
        elif args.view == 'tetra':
            figure = tetra_view(city)
            ax = figure.axes[0]
        elif args.view == 'city':
            ax = physical_scene(city).plot(show=False, max_elements=1000000)
        elif args.view == 'wind':
            ax = frame.plot(representation='height_2m', field='velocity', show=False, max_elements=50000)
        elif args.view == 'streamlines':
            ax = frame.plot(representation='streamlines', show=False, max_elements=30000)
        elif args.view == 'solar':
            # Include every building instead of letting an early large facade
            # consume the generic preview budget and omit later buildings.
            count = sum(len(f.values) for obj in features.values() for rep in obj.geometry.values()
                        for f in getattr(rep.geometry, 'fields', []) if f.name == 'solar_irradiance')
            ax = city.plot(field='solar_irradiance', show=False, max_elements=max(1, count))
        else:
            ax = city.plot(representation='dem', show=False, max_elements=50000)
        figure = ax.figure
        credit(figure)
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(args.save, dpi=150)
        print(args.save.resolve())
    else:
        plt.show()


if __name__ == '__main__':
    main()
