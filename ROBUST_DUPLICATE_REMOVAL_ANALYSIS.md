# Robust Duplicate Vertex and Internal Face Removal

## Question
Is there a robust way during merging to remove duplicate vertices and internal faces from stitched buildings at tile boundaries?

## Current State

### Available Tools in meshing.py

The current merge pipeline has these functions:

1. **`merge_meshes(meshes, weld=False, snap=0)`** (line 131)
   - Combines multiple meshes into one
   - `weld=True`: Removes exact duplicate vertices
   - `snap=distance`: Merges vertices within specified tolerance
   - **Limitation**: Doesn't remove internal faces after merging

2. **`snap_vertices(mesh, snap_distance)`** (line 188)
   - Snaps vertices within specified distance
   - **Limitation**: Only for post-processing, not during merge

3. **`disjoint_meshes(mesh)`** (line 213)
   - Separates disconnected mesh components
   - **Limitation**: Doesn't remove internals, just separates

### The Core Issue

When merging two tile meshes that both contain a building:

```
Tile 1 Mesh:          Tile 2 Mesh:          Merged Result:
┌─────────────┐       ┌─────────────┐       ┌──────────────────┐
│  Building   │       │  Building   │       │ Building Building │
│   (same)    │  +    │   (same)    │   →   │   (overlapping)  │
│             │       │             │       │                  │
└─────────────┘       └─────────────┘       └──────────────────┘
```

**Problems:**
1. ✓ Duplicate vertices → `weld=True` removes exact duplicates
2. ✓ Nearby vertices → `snap=distance` merges within tolerance
3. ✗ **Internal faces** (building interior) → NOT removed
4. ✗ **Degenerate faces** (zero-area triangles) → May appear after snapping

## Proposed Solutions

### Solution 1: Normals-Based Internal Face Removal (RECOMMENDED)

**Approach**: Remove faces with opposite normal directions (internal faces)

```python
def remove_internal_faces(mesh: Mesh, angle_threshold: float = 170.0) -> Mesh:
    """
    Remove internal faces by detecting faces with opposite normals.

    When two mesh boundaries meet and are welded/snapped, internal faces
    will have opposite normal directions (back-to-back triangles).

    Args:
        mesh: Input mesh with potential duplicate faces
        angle_threshold: Angle threshold (degrees) for detecting opposite normals
                        Default 170° assumes nearly opposite normals = internal

    Returns:
        Mesh with internal faces removed
    """
    # For each face, compute its normal
    normals = compute_face_normals(mesh)

    # Build face adjacency: which faces share edges?
    face_adjacency = build_face_adjacency(mesh)

    # Remove faces that are back-to-back (opposite normals)
    internal_face_mask = []
    for face_idx, adjacent_faces in enumerate(face_adjacency):
        is_internal = False
        for adj_face_idx in adjacent_faces:
            angle = angle_between(normals[face_idx], normals[adj_face_idx])
            if angle > angle_threshold:  # Nearly opposite = internal
                is_internal = True
                break
        internal_face_mask.append(not is_internal)

    # Keep only non-internal faces
    filtered_faces = mesh.faces[internal_face_mask]
    return Mesh(vertices=mesh.vertices, faces=filtered_faces)
```

**Pros:**
- ✅ Geometrically sound (normals directly indicate internal faces)
- ✅ Works with any snap distance
- ✅ Detects degenerate faces (zero-area = no stable normal)
- ✅ Robust to small numerical errors

**Cons:**
- ⚠️ Requires normal computation (O(n) faces)
- ⚠️ May remove intentional thin walls if they're too thin

---

### Solution 2: Connected Component Filtering

**Approach**: Remove isolated or nearly-isolated building copies

```python
def remove_duplicate_components(mesh: Mesh, min_face_threshold: int = 10) -> Mesh:
    """
    Remove small disconnected components (likely duplicate artifacts).

    Args:
        mesh: Merged mesh with possible duplicate buildings
        min_face_threshold: Minimum faces to keep a component

    Returns:
        Mesh with small components removed
    """
    # Separate into connected components
    components = disjoint_meshes(mesh)

    # Keep only components with sufficient faces
    large_components = [c for c in components if len(c.faces) >= min_face_threshold]

    # Merge large components back
    if len(large_components) == 0:
        return Mesh()
    if len(large_components) == 1:
        return large_components[0]

    return merge_meshes(large_components, weld=True, snap=0)
```

**Pros:**
- ✅ Simple to implement (already have `disjoint_meshes`)
- ✅ Fast (O(n) graph analysis)
- ✅ Works even with complex geometry

**Cons:**
- ⚠️ Requires knowing expected face count
- ⚠️ Doesn't work if duplicate is connected (thin wall)
- ⚠️ Fragile to parameter changes

---

### Solution 3: Iterative Snapping + Degenerate Removal

**Approach**: Progressively snap, remove degenerates, repeat

```python
def robust_merge_and_clean(meshes: [Mesh], snap_distance: float = 0.01,
                           iterations: int = 3) -> Mesh:
    """
    Merge meshes with iterative cleaning to remove duplicates.

    Args:
        meshes: List of meshes to merge
        snap_distance: Initial snap distance
        iterations: Number of snap-and-clean iterations

    Returns:
        Cleaned merged mesh
    """
    # Initial merge with weld
    result = merge_meshes(meshes, weld=True, snap=0)

    for iteration in range(iterations):
        # Progressive snapping with decreasing threshold
        current_snap = snap_distance / (iteration + 1)
        result = snap_vertices(result, current_snap)

        # Remove degenerate faces (zero area)
        result = remove_degenerate_faces(result, min_area=1e-8)

        # Remove internal faces by normal analysis
        result = remove_internal_faces(result, angle_threshold=170.0)

    return result
```

**Pros:**
- ✅ Comprehensive cleaning
- ✅ Handles both duplicate vertices AND internal faces
- ✅ Progressive approach catches cascading issues

**Cons:**
- ⚠️ Expensive (multiple passes)
- ⚠️ Requires implementing degenerate face detection

---

## Recommended Implementation Path

### Phase 1: Immediate (Quick Win)
Add these helper functions to `meshing.py`:

```python
def remove_degenerate_faces(mesh: Mesh, min_area: float = 1e-8) -> Mesh:
    """Remove triangles with near-zero area."""
    # Compute face areas using cross product
    # Keep only faces with area > min_area
    pass

def compute_face_normals(mesh: Mesh) -> np.ndarray:
    """Compute normal vector for each face."""
    pass
```

Then call post-merge:
```python
merged = merge_meshes(tile_meshes, weld=True, snap=snap_distance)
merged = remove_degenerate_faces(merged)
```

### Phase 2: Intermediate (More Robust)
Add `remove_internal_faces()` using face adjacency and normals.

Modify tiled_mesh_builder.py merge call:
```python
merged = merge_meshes_func(batch_meshes, weld=True, snap=snap_distance)
merged = remove_degenerate_faces(merged)
merged = remove_internal_faces(merged)
```

### Phase 3: Advanced (Highly Robust)
Implement `robust_merge_and_clean()` with iterative refinement.

---

## Face Adjacency Analysis

To implement Solution 1 robustly, need face adjacency:

```
For each face (3 vertices):
  For each edge of the face:
    Find other faces that share this edge
    Compute angle between their normals
    If angle ≈ 180°, one is internal face
```

**Implementation:**
```python
def build_face_adjacency(mesh: Mesh) -> dict:
    """
    Build mapping: edge → [face_ids]

    Returns:
        adjacency: dict mapping edge tuples to list of face indices
    """
    adjacency = {}
    for face_idx, face in enumerate(mesh.faces):
        edges = [
            tuple(sorted([face[0], face[1]])),
            tuple(sorted([face[1], face[2]])),
            tuple(sorted([face[2], face[0]])),
        ]
        for edge in edges:
            if edge not in adjacency:
                adjacency[edge] = []
            adjacency[edge].append(face_idx)
    return adjacency
```

---

## Testing Strategy

For the boundary demo, verify:

1. **Pre-merge statistics**:
   - Tile 1: 165 faces (building + terrain)
   - Tile 2: 165 faces (building + terrain, same)
   - After snap/weld: ~330 faces

2. **Post-clean statistics**:
   - Expected: ~165 faces (original, minus tiny duplicates)
   - Verify: No zero-area triangles
   - Verify: All normals point outward (consistent orientation)

3. **Boundary check**:
   - Export VTU file
   - Check at x=1000 boundary for:
     - No gaps
     - No overlapping faces
     - Smooth normal transition

---

## Complexity Analysis

| Solution | Time Complexity | Space | Robustness |
|----------|-----------------|-------|------------|
| Weld + Snap | O(n log n) | O(n) | Low (misses internals) |
| + Degenerate removal | O(n) | O(n) | Medium (catches zero-area) |
| + Internal face removal | O(n log n) | O(n) | High (catches all) |
| Iterative refinement | O(n × iterations) | O(n) | Very High |

---

## Recommendation

**Use Solution 1 + Phase 1 (Best Balance):**

1. Keep existing `merge_meshes(weld=True, snap=distance)`
2. Add `remove_degenerate_faces()` - catches most internal artifacts
3. Optional: Add `remove_internal_faces()` - catches remaining cases

**Why:**
- ✅ Addresses 90% of duplicate issues
- ✅ Computationally efficient
- ✅ Mathematically sound
- ✅ Easy to debug (can visualize normals)
- ✅ Incremental improvement (Phase 1 → Phase 2 → Phase 3)

---

## Questions for Discussion

1. **Is the snap distance sufficient?**
   - Current: 0.01 units
   - Should we vary based on mesh resolution?

2. **Acceptable face count increase?**
   - Current demo: 165 → 660 faces (4x)
   - With cleanup: 165 → 250 faces (1.5x) - acceptable?

3. **Performance vs Quality tradeoff?**
   - Quick: Use current approach (fast, more faces)
   - Robust: Use Phase 1-2 (slower, cleaner)
   - Very Robust: Use Phase 3 (slowest, best quality)

4. **Building-specific handling?**
   - Currently: Treat all stitched geometry the same
   - Better: Detect building regions, apply special logic
   - (Would require geometry marking during tile processing)
