# Boundary Handling in Tiled Mesh Building

## Problem Statement

When tile boundaries cut through buildings, the mesh needs to be split and then stitched back together:

```
Building spanning boundary:

Tile 1 (0-1000m)      Tile 2 (1000-2000m)
┌─────────────┤Building│───────────────┐
│             │(990-1010m)│             │
└─────────────┴───────────┴─────────────┘
              ↑
           Boundary at 1000m
```

## Current Implementation

### Strategy: Overlap-Based Inclusion

1. **Expand tile bounds by overlap distance** (default 50m)
   ```
   Tile 1 bounds: 0-1000m
   Tile 1 expanded: 0-1050m  (includes buffer)

   Tile 2 bounds: 1000-2000m
   Tile 2 expanded: 950-2000m  (includes buffer)
   ```

2. **Include buildings in both tiles if they overlap expanded bounds**
   ```python
   # Building at 990-1010m overlaps both:
   # - Tile 1 expanded (0-1050m) ✓
   # - Tile 2 expanded (950-2000m) ✓
   # Result: Building meshed in BOTH tiles
   ```

3. **Mesh stitching via weld during merge**
   ```python
   merged = merge_meshes(batch_meshes, weld=True)
   # Removes duplicate vertices within default tolerance
   ```

### Advantages
- ✓ Simple implementation
- ✓ Handles boundary buildings automatically
- ✓ No special case code for boundaries
- ✓ Works with existing merge infrastructure

### Disadvantages
- ✗ Inefficient (same building meshed twice)
- ✗ Relies on numerical precision for welding
- ✗ Snap parameter not used
- ✗ Potential for mesh discontinuities
- ✗ No explicit boundary optimization
- ✗ Memory overhead for duplicate data

---

## Issues and Solutions

### Issue 1: Weld Tolerance Not Specified

**Problem:**
```python
merged = merge_meshes(batch_meshes, weld=True)
# What tolerance does weld use?
# Default might be too loose or too tight
```

**Solution:**
```python
# Use snap parameter to specify tolerance
merged = merge_meshes(
    batch_meshes,
    weld=True,
    snap=0.01  # Merge vertices within 0.01 units
)
```

**Recommended values:**
- `snap=0.001` - Very strict, for precise geometry
- `snap=0.01` - Typical, handles normal floating point errors
- `snap=0.1` - Loose, for coarse meshes
- `snap=1.0` - Very loose, not recommended for buildings

### Issue 2: Duplicate Mesh Data

**Problem:**
```
Building at boundary processed twice:

Tile 1 mesh: 500MB (includes boundary building)
Tile 2 mesh: 500MB (includes same boundary building)
             +─────────────────────────
Total:      1000MB (duplicate 50MB for boundary building)

Memory during merge: 1000MB (both in memory)
After merge:        950MB (weld removes duplicates)
```

**Current solution:**
- Accept the duplication
- Weld removes it during merge
- Temporary memory spike, but acceptable

**Better solution (future):**
- Process boundary vertices only once
- Share between tiles via reference
- More complex implementation

### Issue 3: Mesh Quality at Boundaries

**Problem:**
Triangle density might differ on each side of boundary:

```
Tile 1 triangulation:        Tile 2 triangulation:
  └─────────────────┐       ┌──────────────────┘
    (fine triangles)│       │(coarse triangles)
                    └───┬───┘
                   Boundary mismatch!
```

**Current solution:**
- Same meshing parameters used for both tiles
- Same building surfaces triangulated identically
- Weld should match vertices correctly

**Potential issue:**
- Floating point differences in triangulation
- Different mesh algorithms could produce different triangles
- Vertices might not align exactly

### Issue 4: Overlap Size Parameter

**Current default: `tile_overlap=50.0`**

Is 50m enough? It depends on mesh resolution:

```
Mesh triangle size: 5m
Overlap: 50m (10 triangles)
Building near boundary might extend up to 25m into overlap
Result: Building properly included in both tiles ✓

Mesh triangle size: 50m
Overlap: 50m (1 triangle)
Small building might be cut at boundary
Result: Building might be missed ✗
```

**Better approach:**
```python
def calculate_optimal_overlap(max_mesh_size, min_buffer_triangles=2):
    """
    Calculate overlap size based on mesh resolution.
    Ensure at least N triangles of overlap buffer.
    """
    return max_mesh_size * min_buffer_triangles + some_margin
```

---

## Improved Implementation

### 1. Add Snap Parameter Control

**File:** `tiled_mesh_builder.py`

```python
def _merge_meshes_from_disk(
    store: SharedMeshStore,
    batch_size: int = 2,
    progress: bool = True,
    snap: float = 0.01,  # NEW: Snap tolerance
) -> Mesh:
    """
    Merge mesh tiles from disk with boundary snapping.

    Args:
        snap: Vertex snap distance (default 0.01 units)
    """
    while len(tiles) > 1:
        for batch in batches:
            merged = merge_meshes_func(
                batch_meshes,
                weld=True,
                snap=snap  # USE SNAP PARAMETER!
            )
```

And propagate through:

```python
def build_city_mesh_tiled(
    ...
    snap_distance: float = 0.01,  # NEW
    ...
):
    ...
    result_mesh = _merge_meshes_from_disk(
        store,
        progress=progress,
        snap=snap_distance  # Pass it through
    )
```

### 2. Add Boundary Validation

```python
def validate_mesh_boundaries(mesh, tolerance=0.1):
    """
    Check for mesh discontinuities at boundaries.
    """
    # Find boundary edges (edges with only one adjacent face)
    boundary_edges = find_boundary_edges(mesh)

    # For each boundary edge, check for:
    # 1. Hanging vertices (vertices without proper connections)
    # 2. Gaps (distance between vertices > tolerance)
    # 3. Overlaps (vertices in wrong order)

    issues = []
    for edge in boundary_edges:
        v1, v2 = edge
        dist = distance(v1, v2)
        if dist > tolerance:
            issues.append(f"Gap at boundary: {dist}")

    return issues
```

### 3. Explicit Overlap Calculation

```python
def calculate_tile_overlap(max_mesh_size, min_buffer_size=2):
    """
    Calculate optimal overlap based on mesh parameters.

    Args:
        max_mesh_size: Maximum triangle size in mesh
        min_buffer_size: Minimum number of triangle-widths to overlap

    Returns:
        Recommended overlap distance
    """
    return max_mesh_size * min_buffer_size + 10.0  # 10m minimum margin
```

---

## Testing Boundary Handling

### Test 1: Building at Exact Boundary

```python
city = create_test_city()
# Place building exactly at tile boundary: x=1000

# Tile 1: 0-1000m + 50m overlap = 0-1050m ✓ (includes building)
# Tile 2: 1000-2000m - 50m overlap = 950-2000m ✓ (includes building)

mesh = build_city_mesh(city, use_tiling=True, tile_size=1000.0)

# Verify:
# 1. No gaps in mesh
# 2. No duplicate vertices (or properly welded)
# 3. Building mesh continuous
```

### Test 2: Multiple Buildings at Boundary

```python
city = create_test_city()
# Multiple buildings at x=995, x=1000, x=1005

mesh = build_city_mesh(city, use_tiling=True, tile_overlap=100.0)

# Verify:
# 1. All buildings properly meshed
# 2. No missing triangles
# 3. Mesh topology valid
```

### Test 3: Different Snap Tolerances

```python
for snap_tolerance in [0.001, 0.01, 0.1, 1.0]:
    mesh = build_city_mesh(
        city,
        use_tiling=True,
        snap_distance=snap_tolerance
    )
    # Verify mesh validity
    assert is_valid_mesh(mesh)
```

### Test 4: Mesh Quality Metrics

```python
def measure_boundary_quality(mesh, boundary_vertices):
    """
    Measure mesh quality near boundaries.
    """
    # Triangle aspect ratio near boundaries
    # Edge length distribution
    # Vertex degree (how many edges per vertex)
    # Angle quality
```

---

## Recommended Parameters

### For Precise Geometry
```python
mesh = build_city_mesh(
    city,
    use_tiling=True,
    tile_size=1000.0,
    tile_overlap=100.0,  # Large overlap
    snap_distance=0.001,  # Tight tolerance
    max_mesh_size=2.0,    # Fine triangles
)
```

### For Balanced Performance
```python
mesh = build_city_mesh(
    city,
    use_tiling=True,
    tile_size=1000.0,
    tile_overlap=50.0,   # Default
    snap_distance=0.01,  # Default
    max_mesh_size=10.0,  # Normal
)
```

### For Large Cities with Limited Resources
```python
mesh = build_city_mesh(
    city,
    use_tiling=True,
    tile_size=500.0,     # Smaller tiles
    tile_overlap=25.0,   # Reduced overlap
    snap_distance=0.1,   # Loose tolerance
    max_mesh_size=20.0,  # Coarser mesh
)
```

---

## Known Limitations

1. **Duplicate Processing**
   - Buildings at boundaries processed twice
   - Temporary memory overhead
   - Weld removes duplicates but creates temporary spike

2. **Numerical Precision**
   - Floating point errors might prevent exact welding
   - Snap parameter helps but not foolproof
   - Very tight tolerances (< 1mm) might fail

3. **Mesh Quality Variation**
   - Mesh quality might vary near boundaries
   - Different algorithms might produce different results
   - No guarantee of smooth transition

4. **Complex Building Shapes**
   - Very thin buildings might be missed
   - Concave buildings might have issues
   - Non-convex shapes might have edge cases

---

## Future Improvements

### 1. Intelligent Boundary Processing

Process boundary vertices once, share between tiles:

```python
def process_boundary_vertices(tile1, tile2, overlap_region):
    """
    Identify vertices in overlap region.
    Process them once and share between tiles.
    """
    # Find vertices in overlap
    # Create unique list
    # Reference from both tiles
```

### 2. Adaptive Mesh Refinement at Boundaries

Increase triangle density near boundaries:

```python
def refine_boundary_mesh(mesh, boundary_region, refinement_factor=2):
    """Increase mesh density near boundaries"""
```

### 3. Better Validation

```python
def validate_boundary_stitching(mesh, tolerance=0.01):
    """Check for gaps, overlaps, and discontinuities"""
    # Verify:
    # - No hanging vertices
    # - No gaps
    # - No overlapping triangles
    # - Proper face winding
```

### 4. Automatic Overlap Calculation

```python
optimal_overlap = calculate_overlap_from_mesh_quality(
    max_mesh_size,
    building_complexity,
    required_quality
)
```

---

## Summary

**Current approach:**
- ✓ Works for most cases
- ✓ Simple implementation
- ✗ Not optimal
- ✗ Potential boundary issues

**Improvements needed:**
- [ ] Use snap parameter explicitly
- [ ] Add boundary validation
- [ ] Better overlap calculation
- [ ] Mesh quality verification
- [ ] Documentation of limitations

**Action items:**
- Add `snap_distance` parameter to API
- Add boundary validation function
- Add tests for boundary cases
- Document limitations clearly
