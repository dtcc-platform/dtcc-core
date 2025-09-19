# Mesh Quality Issues (from building footprint mesh)

When generating triangular meshes directly from building footprints, common issues that break simulations include:

- **Long skinny triangles** (high aspect ratio) - numerically unstable for FEM/CFD.  
- **Non-manifold geometry** - edges shared by too many or too few triangles.  
- **Overlapping or intersecting triangles** - folded surfaces confuse solvers.  
- **Disconnected islands** - isolated building patches not linked to the main domain.  
- **Large variation in triangle size** - unstable timesteps and solver convergence.  

---

# Open Source Quality-Aware Meshing Tools (Mainstream)

| Software  | License   | 2D   | 3D  | Scripting         | Quality constraints             | Adaptivity                       |
|-----------|-----------|------|-----|-------------------|---------------------------------|----------------------------------|
| Triangle  | Freeware* | Yes  | No  | CLI only          | Min angle, area, CDT            | No                               |
| Gmsh      | GPL       | Yes  | Yes | Python, C++       | Size fields, gradation, opt.    | Yes (refine, coarsen, smooth)    |
| JIGSAW    | MPL 2.0   | Yes  | Yes | Python            | Min angle, size, gradation      | Yes (dynamic)                    |
| CGAL      | LGPL/GPL  | Yes  | Yes | C++ (PyGalmesh)   | Min angle, size, gradation      | Yes (refine)                     |
| TetGen    | GPL v2    | Yes* | Yes | CLI, APIs         | Min angle, volume               | Limited                          |
| Netgen    | LGPL      | Yes  | Yes | Python, C++       | Min angle, size                 | Yes (adaptive)                   |
| MMG       | LGPL      | Yes  | Yes | CLI, C, Fortran   | Aspect ratio, size improvement  | Yes (adaptive remeshing)         |
| Salome    | LGPL      | Yes  | Yes | Python, GUI       | Multiple algos (e.g. Netgen)    | Yes                              |

---

# Lightweight Libraries (from ecosystem)

| Software        | License     | 2D Meshing                                                   | 3D Meshing                          | Scripting (Python, C++, CLI)      | Quality Constraints (generation-time)               |
|-----------------|-------------|--------------------------------------------------------------|-------------------------------------|-----------------------------------|-----------------------------------------------------|
| poly2tri        | BSD-3-Clause| Constrained Delaunay of polygons; supports holes & Steiner pts| —                                   | C++ library (CMake); no CLI        | None built-in; indirect via segmentation/Steiner    |
| delaunator-cpp  | MIT         | Unconstrained Delaunay triangulation of 2D point sets        | —                                   | C++ header-only; no CLI            | None; external refinement needed                    |
| CDT (artem-ogre/CDT) | MPL-2.0 | Constrained & conforming Delaunay; preserves input segments  | —                                   | C++ lib; optional Python bindings (PythonCDT) | No native quality switches; refinement loops needed |
| libDistMesh     | GPL-2.0     | DistMesh algorithm for unstructured triangular meshes via size field h(x) | **Yes**: unstructured tetrahedral   | C++ lib (depends on Eigen & Qhull) | Quality via size field h(x); not explicit min-angle |

---

# Notes
- *Triangle*: free for non-commercial use, not OSI-compliant open source.  
- *TetGen*: mainly 3D but can handle 2D polygons.  
- **JIGSAW**: most permissive modern option (MPL 2.0).  
- **libDistMesh**: not a CDT, relies on iterative node relaxation + size field.  
- Lightweight libs (poly2tri, delaunator, CDT, libDistMesh) are **lower-level building blocks** compared to full-featured meshing platforms (Gmsh, CGAL, Netgen).  

