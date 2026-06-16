# This demo builds LOD2 buildings for three areas in Kungsbacka.

from dtcc_core.builder.meshing import merge_meshes
from dtcc_core.io import download_footprints, download_pointcloud, save_mesh
from dtcc_core.model import City

from kungsbacka_domains import AREAS, OUTPUT_ROOT


def building_mesh(building):
    geometry = building.lod2 if building.lod2 is not None else building.lod1
    if geometry is None:
        return None
    return geometry.mesh(weld=True, snap=0.005)


def city_building_mesh(city):
    meshes = []
    for building in city.buildings:
        mesh = building_mesh(building)
        if mesh is not None and len(mesh.faces) > 0:
            meshes.append(mesh)
    if len(meshes) == 0:
        raise ValueError("No LOD2 or LOD1 fallback building meshes were generated.")
    return merge_meshes(meshes, weld=True, snap=0.005)


output_dir = OUTPUT_ROOT / "meshes_kungsbacka_lod2"
output_dir.mkdir(parents=True, exist_ok=True)


for area in AREAS:
    # Build LOD2 buildings (watertight LOD2 where supported, LOD1 fallback)
    city = City()
    city.bounds = area["bounds"]
    city.add_buildings(download_footprints(area["bounds"]))
    city.add_pointcloud(download_pointcloud(area["bounds"]))
    city.build_lod2_buildings(calculate_heights=True, log_rejections=True)

    lod2_count = sum(1 for building in city.buildings if building.lod2 is not None)
    print(
        f"{area['name']}: buildings={len(city.buildings)} "
        f"lod2={lod2_count} lod1_fallback={len(city.buildings) - lod2_count}"
    )

    mesh = city_building_mesh(city)
    save_mesh(mesh, output_dir / f"lod2_mesh_{area['name']}.xdmf")
    save_mesh(mesh, output_dir / f"lod2_mesh_{area['name']}.stl")
