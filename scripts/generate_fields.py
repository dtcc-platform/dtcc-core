"""Generate synthetic time-dependent volume fields for visualization testing.

Run from the installed Core checkout:
    python scripts/generate_fields.py

Each time has a VolumeMesh and VolumeGrid .dtcc snapshot and matching .vtu
files. Open the two .pvd collections in ParaView to animate the series.
The domain follows the saved flagship, including its expanded volume height.
See docs/synthetic-fields.md for the field construction and time encoding.
"""

import argparse
import json
from pathlib import Path
import xml.etree.ElementTree as ET

import meshio
import numpy as np

from dtcc_core import io
from dtcc_core.model import City, Field, VolumeGrid, VolumeMesh, exchange


DEFAULT_NX = 128
DEFAULT_NY = 128
DEFAULT_NZ = 36
DEFAULT_NT = 21
DEFAULT_DURATION = 10.0
DEFAULT_OUTPUT = Path(__file__).resolve().parents[1] / "data" / "fields"
DEFAULT_FLAGSHIP = Path(__file__).resolve().parents[1] / "data" / "flagship" / "flagship.dtcc"


def flagship_domain(path):
    """Read the saved flagship's spatial domain and coordinate reference."""
    city = io.load_model(path, expected_type=City)
    bounds = city.bounds.copy()
    lengths = np.array([bounds.width, bounds.height, bounds.depth])
    if not np.isfinite(lengths).all() or np.any(lengths <= 0):
        raise ValueError("Flagship bounds must have finite, positive extents in all three directions")
    if not city.transform.srs or not np.array_equal(city.transform.affine, np.eye(4)):
        raise ValueError("Flagship must declare a coordinate reference and use an identity root transform")
    return bounds, city.transform.srs


def make_geometries(nx, ny, nz, bounds, crs):
    """Build a common lattice, with six conforming, positive tets per cube."""
    grid = VolumeGrid(width=nx, height=ny, depth=nz)
    grid.bounds = bounds.copy()
    grid.transform.srs = crs
    vertices = grid.coordinates()
    # VolumeGrid coordinates flatten (y, x, z), with z varying fastest.
    dx, dy = nz + 1, (nx + 1) * (nz + 1)
    base = (np.arange(ny, dtype=np.int32)[:, None, None] * dy
            + np.arange(nx, dtype=np.int32)[None, :, None] * dx
            + np.arange(nz, dtype=np.int32)[None, None, :]).ravel()
    # VTK hexahedron order: four counterclockwise bottom corners, then top.
    offsets = np.array([0, dx, dx + dy, dy, 1, dx + 1, dx + dy + 1, dy + 1],
                       dtype=np.int32)
    hexahedra = base[:, None] + offsets
    pattern = np.array([[0, 1, 2, 6], [0, 2, 3, 6], [0, 3, 7, 6],
                        [0, 7, 4, 6], [0, 4, 5, 6], [0, 5, 1, 6]])
    tetra = VolumeMesh(vertices=vertices, cells=hexahedra[:, pattern].reshape(-1, 4))
    tetra.transform.srs = crs
    return tetra, grid, hexahedra


def sample_fields(xyz, time_seconds, period_seconds):
    """Sample synthetic flow/pressure in normalized [-pi, pi] coordinates.

    The polynomial/sine XY flow extends the proposed 2D field. Helical waves
    and three orbiting vortex/jet cores add structure in all three directions.
    Pressure combines traveling waves, moving peaks and velocity magnitude.
    These are visualization fixtures, not solutions of fluid equations.
    """
    phase = 2 * np.pi * (time_seconds / period_seconds)
    x, y, z = xyz.T
    xp, yp = x + .55 * np.sin(phase), y + .55 * np.sin(phase)
    velocity = np.column_stack([
        -np.sin(yp + .7*z + phase) + .1*(xp*xp - 2*xp*yp)
        + .75*np.sin(2*z - phase) + .55*np.cos(2*y + phase),
        np.sin(xp - .7*z - phase) + .1*(yp*yp - 2*xp*yp)
        + .75*np.sin(2*x + phase) + .55*np.cos(2*z - phase),
        .75*np.sin(2*y - phase) + .55*np.cos(2*x + phase)
        + .4*np.sin(x + y + z - 2*phase),
    ])
    pressure = (1.1*np.sin(2*x + phase)*np.cos(2*y - phase)*np.cos(z)
                + .6*np.cos(x - 2*z - 2*phase)*np.sin(2*y + phase))
    axes = np.array([[0., 0., 1.], [1., 0., .4], [0., 1., -.4]])
    axes /= np.linalg.norm(axes, axis=1)[:, None]
    for i, axis in enumerate(axes):
        angle = phase + i * 2*np.pi/3
        centre = np.array([1.5*np.cos(angle), 1.5*np.sin(angle),
                           .9*np.sin(phase + i)])
        relative = xyz - centre
        radius = .85 + .15*np.sin(angle)
        core = np.exp(-np.sum(relative*relative, axis=1) / radius**2)
        strength = (-1 if i == 1 else 1) * 3*(1 + .25*np.sin(angle))
        velocity += core[:, None] * (
            strength*np.cross(axis, relative)/radius + .7*np.cos(angle)*axis)
        pressure += (-2.8 if i != 1 else 3.2) * core * (1 + .3*np.sin(angle))
    pressure = 50 * (pressure - .12*np.sum(velocity*velocity, axis=1))
    return [
        Field(name="velocity", unit="m/s", dim=3, association="vertex",
              values=velocity.astype(np.float32),
              description="Synthetic projected (easting, northing, up) velocity; waves and moving vortex/jet cores."),
        Field(name="pressure", unit="Pa", dim=1, association="vertex",
              values=pressure.astype(np.float32),
              description="Synthetic signed gauge pressure; waves, moving peaks and speed contribution."),
    ]


def write_pvd(path, frames, geometry):
    """Write relative VTU references with explicit times for ParaView."""
    root = ET.Element("VTKFile", type="Collection", version="0.1", byte_order="LittleEndian")
    collection = ET.SubElement(root, "Collection")
    for frame in frames:
        ET.SubElement(collection, "DataSet", timestep=repr(frame["time_seconds"]),
                      group="", part="0", file=frame[geometry]["vtu"])
    ET.indent(root)
    ET.ElementTree(root).write(path, encoding="utf-8", xml_declaration=True)


def generate(nx, ny, nz, nt, duration, output_dir, bounds, crs):
    """Write one pair of snapshots at a time, reusing fixed geometry."""
    tetra, grid, hexahedra = make_geometries(nx, ny, nz, bounds, crs)
    origin = np.array([bounds.xmin, bounds.ymin, bounds.zmin])
    lengths = np.array([bounds.width, bounds.height, bounds.depth])
    # Keep the analytic features distributed across the projected city domain.
    sample_xyz = 2*np.pi*((tetra.vertices - origin) / lengths - .5)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"{grid.num_vertices:,} vertices; {len(tetra.cells):,} tetrahedra; "
          f"{grid.num_cells:,} grid cells; {nt} snapshots", flush=True)
    frames = []
    for index in range(nt):
        seconds = duration * (index / (nt - 1)) if nt > 1 else 0.
        fields = sample_fields(sample_xyz, seconds, duration)
        tetra.fields = fields
        grid.fields = fields
        frame = {"index": index, "time_seconds": seconds}
        for name, geometry in (("tet", tetra), ("grid", grid)):
            stem = f"fields_{name}_{index:04d}"
            dtcc_name, vtu_name = f"{stem}.dtcc", f"{stem}.vtu"
            io.save_model(geometry, output_dir / dtcc_name)
            if name == "tet":
                io.save_volume_mesh(tetra, output_dir / vtu_name)
            else:
                meshio.write(output_dir / vtu_name, meshio.Mesh(
                    tetra.vertices, [("hexahedron", hexahedra)],
                    point_data={field.name: field.values for field in fields}))
            frame[name] = {"dtcc": dtcc_name, "vtu": vtu_name}
        frames.append(frame)
        print(f"  Snapshot {index + 1}/{nt}: t={seconds:g} s", flush=True)
    # Publish time indexes after all the referenced snapshots have been written.
    for name in ("tet", "grid"):
        write_pvd(output_dir / f"fields_{name}.pvd", frames, name)
    manifest = {
        "synthetic": True,
        "description": "Analytic visualization fixtures; not a fluid simulation.",
        "time_encoding": "Separate DTCC geometry snapshots; no native Field time axis.",
        "period_seconds": duration,
        "subdivisions": {"nx": nx, "ny": ny, "nz": nz},
        "bounds": {name: getattr(grid.bounds, name)
                   for name in ("xmin", "ymin", "zmin", "xmax", "ymax", "zmax")},
        "crs": crs,
        "field_association": "vertex",
        "velocity_components": ["easting", "northing", "up"],
        "frames": frames,
    }
    (output_dir / "fields-series.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Wrote {output_dir.resolve()}\n"
          "Open fields_tet.pvd and fields_grid.pvd in ParaView.")


def positive_int(value):
    try:
        result = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a positive integer") from None
    if result < 1:
        raise argparse.ArgumentTypeError("must be a positive integer")
    return result


def positive_float(value):
    try:
        result = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be finite and positive") from None
    if not np.isfinite(result) or result <= 0:
        raise argparse.ArgumentTypeError("must be finite and positive")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    for name, default_count in (("nx", DEFAULT_NX), ("ny", DEFAULT_NY), ("nz", DEFAULT_NZ)):
        parser.add_argument(f"-{name}", f"--{name}", type=positive_int, default=default_count,
                            help=f"number of cube subdivisions along {name[-1]} (default: {default_count})")
    parser.add_argument("-nt", "--nt", type=positive_int, default=DEFAULT_NT,
                        help=f"snapshot count including endpoints; 1 writes only t=0 (default: {DEFAULT_NT})")
    parser.add_argument("-T", "--duration", type=positive_float, default=DEFAULT_DURATION,
                        help=f"duration and loop period in seconds (default: {DEFAULT_DURATION:g})")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT,
                        help=f"output directory (default: {DEFAULT_OUTPUT})")
    parser.add_argument("--flagship-model", type=Path, default=DEFAULT_FLAGSHIP,
                        help=f"saved City supplying bounds and CRS (default: {DEFAULT_FLAGSHIP})")
    args = parser.parse_args()
    vertices = (args.nx + 1) * (args.ny + 1) * (args.nz + 1)
    cubes = args.nx * args.ny * args.nz
    # Native mesh arrays: float64 XYZ, float32 velocity/pressure, int32 tets.
    if 40*vertices + 96*cubes >= exchange.MAX_PROTOBUF_BYTES:
        parser.error("mesh arrays alone exceed the .dtcc Protobuf limit (<2 GiB); reduce -nx/-ny/-nz")
    try:
        bounds, crs = flagship_domain(args.flagship_model)
        generate(args.nx, args.ny, args.nz, args.nt, args.duration, args.output_dir, bounds, crs)
    except (OSError, ValueError) as exc:
        parser.exit(1, f"Unable to generate fields: {exc}\n")


if __name__ == "__main__":
    main()
