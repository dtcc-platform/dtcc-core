# MPI Scalability Analysis for Tiled Mesh Building

## Question
Could this be done with MPI on thousands of cores? What's the scaling potential?

---

## TL;DR: YES, with Excellent Scaling Potential ✅

The current architecture is **highly suitable for MPI** and can scale to **thousands of cores** with **near-linear speedup** for tile processing, followed by **logarithmic scaling** for the merge phase.

**Expected Performance:**
- **1000 cores**: Process 1000 tiles in ~same time as 1 tile (perfect speedup)
- **10,000 cores**: Process entire countries/continents feasibly
- **Limitation**: Merge phase bottleneck (but manageable with hierarchical approach)

---

## Current Architecture Analysis

### Existing Implementation (Multiprocessing on Single Node)

```python
# Line 232: Uses Python multiprocessing.Pool
with Pool(n_workers) as pool:
    for tile_id, vertices, faces, markers, bounds in pool.imap_unordered(_process_mesh_tile, tile_args):
        store.add_tile(tile_id, bounds, vertices, faces, markers)
```

**Current Limitations:**
- ❌ Shared memory (single machine)
- ❌ Limited by node RAM/cores
- ❌ Typical limit: 64-128 cores per node
- ❌ No inter-node communication

**Current Strengths:**
- ✅ Embarrassingly parallel tile processing
- ✅ Disk-backed storage (HDF5)
- ✅ Independent worker tasks
- ✅ No inter-tile dependencies during processing

---

## MPI Architecture Proposal

### Phase 1: Distributed Tile Processing (Embarrassingly Parallel)

```python
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

# Master process (rank 0) distributes tiles
if rank == 0:
    # Create tiles
    all_tiles = create_tile_grid(city, tile_size, overlap)

    # Distribute tiles to workers
    tiles_per_worker = len(all_tiles) // size
    my_tiles = all_tiles[rank * tiles_per_worker:(rank + 1) * tiles_per_worker]
else:
    # Worker receives tile assignment
    my_tiles = None

my_tiles = comm.scatter(all_tiles, root=0)

# Each worker processes assigned tiles independently
for tile in my_tiles:
    mesh = process_tile(tile)

    # Write to distributed file system (Lustre, GPFS, etc.)
    # Use MPI-IO for parallel writes
    write_tile_to_parallel_hdf5(mesh, tile_id, rank)
```

**Scaling Characteristics:**
- **Perfect linear scaling**: O(N / P) where N=tiles, P=cores
- **No communication overhead** during processing
- **Memory per rank**: O(1) tile (constant)
- **Theoretical limit**: Limited only by I/O bandwidth

---

### Phase 2: Hierarchical Distributed Merge

Current merge is sequential (bottleneck):
```python
# Line 323: Sequential merge iterations
while len(tiles) > 1:
    # Merge pairs: 1000 → 500 → 250 → ... → 1
    pass
```

**MPI Solution: Parallel Tree Reduction**

```
Iteration 0: 1000 tiles across 1000 ranks
    Rank 0 + Rank 1 → merged_0
    Rank 2 + Rank 3 → merged_1
    ...
    → 500 merged tiles

Iteration 1: 500 tiles across 500 ranks
    Rank 0 + Rank 250 → merged_0
    ...
    → 250 merged tiles

Iteration log2(N): Final merge
    Rank 0 has complete mesh
```

**Implementation:**
```python
def mpi_hierarchical_merge(local_mesh, comm):
    """
    Parallel tree reduction for mesh merging.

    Complexity: O(log2(P)) iterations, each O(mesh_size)
    """
    rank = comm.Get_rank()
    size = comm.Get_size()

    active_ranks = size
    step = 1

    while active_ranks > 1:
        # Pair ranks: (0,1), (2,3), (4,5), ...
        if rank % (2 * step) == 0:
            # Receive from partner
            if rank + step < size:
                partner_rank = rank + step
                partner_mesh = comm.recv(source=partner_rank, tag=0)

                # Merge with local mesh
                local_mesh = merge_meshes([local_mesh, partner_mesh],
                                         weld=True, snap=0.01)
        else:
            # Send to partner
            partner_rank = rank - step
            comm.send(local_mesh, dest=partner_rank, tag=0)
            break  # This rank is done

        active_ranks = active_ranks // 2
        step *= 2

    # Rank 0 has final result
    if rank == 0:
        return local_mesh
    else:
        return None
```

**Scaling Characteristics:**
- **Time Complexity**: O(log2(P)) communication rounds
- **Data Transfer**: Each mesh transmitted once
- **Memory**: O(mesh_size) per active rank
- **Bottleneck**: Network bandwidth between nodes

---

## Scaling Analysis

### Best Case: City with 10,000 Tiles

#### Scenario 1: Current (Single Node, 64 cores)
```
Tile processing: 10,000 tiles / 64 cores = 156 tiles per core
  @ 2 seconds/tile = 312 seconds = 5.2 minutes

Merge phase (sequential):
  log2(10,000) ≈ 14 iterations
  @ 5 seconds/iteration = 70 seconds

Total: 382 seconds ≈ 6.4 minutes
```

#### Scenario 2: MPI (1000 cores across 16 nodes)
```
Tile processing: 10,000 tiles / 1000 cores = 10 tiles per core
  @ 2 seconds/tile = 20 seconds

Merge phase (parallel tree):
  log2(1000) ≈ 10 iterations
  @ 2 seconds/iteration (parallel) = 20 seconds

Total: 40 seconds (9.5x speedup!)
```

#### Scenario 3: MPI (10,000 cores across 160 nodes)
```
Tile processing: 10,000 tiles / 10,000 cores = 1 tile per core
  @ 2 seconds/tile = 2 seconds (perfect scaling!)

Merge phase (parallel tree):
  log2(10,000) ≈ 14 iterations
  @ 2 seconds/iteration = 28 seconds

Total: 30 seconds (12.7x speedup!)
```

---

## Performance Bottlenecks and Solutions

### Bottleneck 1: I/O Bandwidth
**Problem**: 10,000 ranks writing to shared filesystem simultaneously

**Solutions:**
1. **MPI-IO with collective writes**
   ```python
   from mpi4py import MPI
   fh = MPI.File.Open(comm, filename, MPI.MODE_CREATE | MPI.MODE_WRONLY)
   fh.Write_at_all(offset, data)  # Collective, optimized I/O
   ```

2. **Hierarchical I/O** (Node-local aggregation)
   ```
   Rank 0-63 (Node 1) → Local SSD → 1 write to parallel FS
   Rank 64-127 (Node 2) → Local SSD → 1 write to parallel FS
   ```
   Reduces I/O contention by 64x

3. **Burst Buffers** (Modern HPC systems)
   - Write to fast local storage first
   - Flush to parallel FS asynchronously

### Bottleneck 2: Merge Communication Overhead
**Problem**: Network bandwidth for large mesh transfers

**Solutions:**
1. **Compressed transfers**
   ```python
   import zlib
   compressed_mesh = zlib.compress(pickle.dumps(mesh))
   comm.send(compressed_mesh, dest=partner_rank)
   ```
   Reduces transfer by 5-10x

2. **Topology-aware merging**
   ```
   Merge spatially adjacent tiles first (same node/switch)
   Reduces cross-node communication
   ```

3. **Progressive merging**
   ```
   Start merge while tiles still processing
   Pipeline: Process → Merge → Process → Merge
   ```

### Bottleneck 3: Load Balancing
**Problem**: Some tiles have more buildings (longer processing)

**Solutions:**
1. **Dynamic task assignment** (MPI master-worker pattern)
   ```python
   # Master distributes tiles on-demand
   while tiles_remaining:
       worker_rank = comm.recv(source=MPI.ANY_SOURCE, tag=REQUEST)
       comm.send(next_tile, dest=worker_rank, tag=TASK)
   ```

2. **Adaptive tiling**
   ```
   Dense urban areas → smaller tiles (faster processing)
   Rural areas → larger tiles (fewer workers needed)
   ```

---

## Scalability Limits

### Strong Scaling (Fixed Problem Size)

For a city with N tiles on P cores:

| Cores (P) | Tiles (N) | Speedup | Efficiency | Limit |
|-----------|-----------|---------|------------|-------|
| 1 | 10,000 | 1.0x | 100% | - |
| 10 | 10,000 | 9.5x | 95% | Comm overhead |
| 100 | 10,000 | 85x | 85% | Load imbalance |
| 1,000 | 10,000 | 650x | 65% | I/O bandwidth |
| 10,000 | 10,000 | 3,500x | 35% | Merge dominates |

**Practical Strong Scaling Limit**: ~1,000 cores for typical city

**Reason**: When P ≈ N, merge phase dominates (can't parallelize beyond log2(N) levels)

### Weak Scaling (Problem Size Grows with Cores)

Keep tiles per core constant (e.g., 10 tiles/core):

| Cores (P) | Tiles (N) | Time | Efficiency | Limit |
|-----------|-----------|------|------------|-------|
| 1 | 10 | 20s | 100% | - |
| 10 | 100 | 21s | 95% | - |
| 100 | 1,000 | 23s | 87% | - |
| 1,000 | 10,000 | 28s | 71% | Merge comm |
| 10,000 | 100,000 | 38s | 53% | I/O limit |

**Practical Weak Scaling Limit**: >10,000 cores feasible

**Reason**: Merge phase is O(log P), grows slowly

---

## Real-World MPI Implementation Roadmap

### Phase 1: MPI Tile Processing (Easiest, Biggest Impact)
**Effort**: Medium (2-3 weeks)
**Benefit**: 10-100x speedup

```python
# Replace multiprocessing.Pool with MPI
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# Distribute tiles across ranks
my_tiles = distribute_tiles(all_tiles, rank, comm.Get_size())

for tile in my_tiles:
    mesh = process_tile(tile)
    write_to_parallel_hdf5(mesh, tile_id, rank)

# Barrier: Wait for all tiles to finish
comm.Barrier()
```

### Phase 2: Parallel I/O with MPI-IO (Reduces Bottleneck)
**Effort**: Medium (2 weeks)
**Benefit**: 5-10x I/O improvement

```python
# Use MPI-IO for efficient parallel writes
fh = MPI.File.Open(comm, "tiles.h5", MPI.MODE_CREATE | MPI.MODE_WRONLY)
fh.Write_at_all(offset, mesh_data)
fh.Close()
```

### Phase 3: Hierarchical Parallel Merge (Full Scaling)
**Effort**: High (4-6 weeks)
**Benefit**: Scales to 1000+ cores

```python
# Tree reduction merge
final_mesh = mpi_hierarchical_merge(my_merged_tiles, comm)
```

### Phase 4: Advanced Optimizations (Optional)
**Effort**: High (ongoing)
**Benefit**: Incremental improvements

- Topology-aware merging
- Asynchronous I/O (overlap compute + I/O)
- GPU acceleration for meshing
- Adaptive load balancing

---

## Comparison: Multiprocessing vs MPI

| Feature | Multiprocessing (Current) | MPI (Proposed) |
|---------|---------------------------|----------------|
| **Scaling limit** | 64-128 cores (single node) | 1,000-10,000 cores (cluster) |
| **Memory** | Shared memory | Distributed memory |
| **Communication** | Shared memory (fast) | Network (slower) |
| **I/O** | Local disk | Parallel filesystem |
| **Fault tolerance** | Process failure = retry | Rank failure = checkpoint/restart |
| **Setup complexity** | Low (works out of box) | High (MPI environment required) |
| **Code changes** | None | Moderate (replace Pool with MPI calls) |

---

## Cost-Benefit Analysis

### For Small Cities (< 1000 tiles)
**Current approach is sufficient.**
- 64-core node completes in minutes
- MPI overhead > benefit

### For Large Cities (1,000-10,000 tiles)
**MPI provides 5-20x speedup.**
- Example: Stockholm metropolitan area
- Current: ~1 hour on 64 cores
- MPI: ~5 minutes on 1000 cores

### For Country/Continental Scale (>10,000 tiles)
**MPI is essential.**
- Example: Entire Sweden or California
- Current: Infeasible (>24 hours)
- MPI: ~1 hour on 10,000 cores

---

## Recommended Next Steps

1. **Prototype MPI tile processing** (Phase 1)
   - Replace `multiprocessing.Pool` with MPI scatter/gather
   - Test on small cluster (16-64 cores)
   - Measure speedup

2. **Benchmark I/O performance**
   - Test HDF5 with MPI-IO
   - Identify I/O bottlenecks
   - Optimize write patterns

3. **Implement hierarchical merge** (Phase 3)
   - Start with simple tree reduction
   - Test communication overhead
   - Optimize for network topology

4. **Scale testing**
   - Test on 100, 500, 1000 cores
   - Measure strong and weak scaling
   - Identify breaking points

---

## Conclusion

✅ **Yes, this architecture is highly MPI-friendly**

✅ **Scaling to 1,000+ cores is feasible** with expected 10-100x speedup

✅ **Tile processing scales perfectly** (embarrassingly parallel)

✅ **Merge phase scales logarithmically** (acceptable bottleneck)

✅ **Practical limit: 1,000-10,000 cores** depending on problem size

**Key Success Factors:**
- Independent tile processing (no inter-tile dependencies)
- Disk-backed storage (already designed for distributed FS)
- Hierarchical merge strategy (reduces communication)
- Mature MPI ecosystem (mpi4py, HDF5 parallel, etc.)

**Next Action:** Implement Phase 1 (MPI tile processing) as proof-of-concept on a small HPC system.
