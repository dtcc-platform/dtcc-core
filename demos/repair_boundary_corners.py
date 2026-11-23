"""Inspect a tetrahedral corner mesh before and after subdivision."""

import dtcc_core as dtcc

box = dtcc.BoxMesh(1, 1, 1)
cells = box.cells[:1].copy()
vertices = box.vertices[cells[0]]
cells[0] = [0, 1, 2, 3]
original = dtcc.model.VolumeMesh(vertices=vertices, cells=cells)
repaired = original.repair_boundary_corners()

print(
    f"Original: {original.num_cells} cell, {original.find_boundary_corners().num_corners} corners"
)
print(
    f"Repaired: {repaired.num_cells} cells, {repaired.find_boundary_corners().num_corners} corners"
)
original.save("output/boundary_corners_original.vtu")
repaired.save("output/boundary_corners_repaired.vtu")
original.view()
repaired.view()
