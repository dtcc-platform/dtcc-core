from dtcc_core.builder.geometry.create import mesh as create_mesh
import dtcc_viewer


sphere = create_mesh.create_sphere(
    (0, 0, 0), 10.0, resolution=20, scale=(1.0, 1.0, 1.0)
)

cylinder = create_mesh.create_cylinder((0, 0, 0), 10.0, 20.0, resolution=8, axis=2)

cylinder.view()
