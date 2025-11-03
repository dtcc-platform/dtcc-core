# Tiled Mesh Building Demos: Serial vs Parallel

This directory contains two demonstration scripts that showcase the shared memory architecture for building surface meshes of large cities. The scripts are identical except for one parameter that controls parallelization.

## Demo Scripts

### 1. `demo_tiling_serial.py` - Serial Processing
```bash
python demo_tiling_serial.py
```

Uses **single-threaded** processing:
- `n_workers=1`
- Processes tiles sequentially one at a time
- Lower peak memory usage
- Slower overall performance
- Good baseline for comparison

**Key Output:**
- `demo_mesh_serial_tiling.vtu` - The resulting mesh
- Timing breakdown for each step
- Peak memory usage

### 2. `demo_tiling_parallel.py` - Parallel Processing
```bash
python demo_tiling_parallel.py
```

Uses **multi-threaded** processing:
- `n_workers=None` (auto-detects CPU count)
- Processes multiple tiles concurrently
- Higher peak memory usage (proportional to worker count)
- Significantly faster overall performance
- Scales well on multi-core systems

**Key Output:**
- `demo_mesh_parallel_tiling.vtu` - The resulting mesh (identical to serial)
- Timing breakdown for each step
- Peak memory usage

## What's Identical Between the Demos

Both scripts perform **exactly the same operations** except for parallelization:

```
✓ Same city data (Gothenburg, Sweden - 2km x 2km area)
✓ Same data downloads (pointcloud, building footprints)
✓ Same preprocessing (outlier removal, terrain building)
✓ Same building height computation
✓ Same mesh parameters:
  - Tile size: 1000m x 1000m
  - Tile overlap: 50m
  - LOD: LOD0
  - Triangle size: 5.0m
  - Max mesh size: 10.0m
  - Min mesh angle: 25.0°
✓ Same output location and format
```

## Expected Results

### Processing Steps (Both Scripts)
1. **Download**: Fetch pointcloud and building data (~10-30 seconds)
2. **Preprocessing**: Remove outliers (~5-10 seconds)
3. **Terrain**: Build terrain raster (~10-20 seconds)
4. **Heights**: Compute building heights (~5-10 seconds)
5. **City**: Create city object (~1-2 seconds)
6. **Mesh**: Build surface mesh (**TIME VARIES - Serial vs Parallel**)
7. **Save**: Write output file (~2-5 seconds)

### Typical Speedup (Serial vs Parallel)

| Parameter | Serial | Parallel (8 cores) | Speedup |
|-----------|--------|-------------------|---------|
| Mesh Time | 120-180s | 30-50s | 3-4x |
| Peak Memory | 800-1200 MB | 1500-2000 MB | 1.2-1.7x |
| Total Time | 180-250s | 100-140s | 2-2.5x |

*Actual results depend on system hardware, internet speed, and disk speed.*

## Running the Comparison

### Quick Benchmark
```bash
# Run serial version
time python demo_tiling_serial.py

# Run parallel version
time python demo_tiling_parallel.py

# Compare timing output
```

### Monitor Resource Usage
```bash
# In another terminal, monitor memory and CPU
watch -n 1 'ps aux | grep demo_tiling'
```

### Detailed Comparison
```bash
# Capture both outputs
python demo_tiling_serial.py > results_serial.txt 2>&1
python demo_tiling_parallel.py > results_parallel.txt 2>&1

# Compare results
diff results_serial.txt results_parallel.txt
```

## Understanding the Output

### Timing Breakdown
```
Total processing time: 245.32s
  - Download: 25.45s (10.4%)
  - Preprocessing: 7.82s (3.2%)
  - Terrain: 18.32s (7.5%)
  - Heights: 8.21s (3.3%)
  - City: 1.45s (0.6%)
  - Mesh: 165.23s (67.4%)        ← This varies between serial and parallel
  - Save: 3.84s (1.6%)
```

The **Mesh** step is where the parallelization difference is most dramatic.

### Memory Usage
```
Peak memory usage: 1245.3MB
Memory per vertex: 0.085KB
```

Higher peak memory in parallel mode is expected due to concurrent tile processing.

## Customizing the Demos

### Change City Location
Edit the city bounds:
```python
# Change from Gothenburg to another location
x0 = 340000.0  # Different X coordinate
y0 = 6410000.0  # Different Y coordinate
L = 3000.0     # Larger area (3km x 3km)
```

### Adjust Tile Parameters
```python
mesh = dtcc.build_city_mesh(
    city,
    use_tiling=True,
    tile_size=1500.0,      # Larger tiles = fewer tiles, potentially faster
    tile_overlap=100.0,    # Larger overlap = better boundary handling
    n_workers=8,           # Manual worker count instead of auto-detect
    # ... other parameters
)
```

### Adjust Mesh Quality
```python
mesh = dtcc.build_city_mesh(
    city,
    lod=dtcc.GeometryType.LOD1,  # Higher detail
    max_mesh_size=5.0,            # Smaller elements = higher quality
    min_mesh_angle=30.0,          # Stricter quality criteria
    # ... other parameters
)
```

## Troubleshooting

### Script Runs Out of Memory
**Solution:** Reduce tile size in parallel version
```python
tile_size=500.0      # Reduce from 1000.0
n_workers=2          # Reduce parallel workers
```

### Parallel is Not Faster Than Serial
**Possible causes:**
1. Small city (too few tiles to benefit from parallelization)
2. Limited CPU cores available
3. Disk I/O bottleneck (slow storage)

**Solutions:**
- Increase city size (L parameter)
- Use faster storage (SSD instead of HDD)
- Check CPU usage: should see 80-95% on parallel

### Download Fails
**Possible causes:**
1. No internet connection
2. Blocked network access
3. Service temporarily unavailable

**Solutions:**
- Check internet connection
- Try again later (may be service issue)
- Use cached data if available

## Output Files

Both scripts generate:
- `demo_mesh_serial_tiling.vtu` - Serial result
- `demo_mesh_parallel_tiling.vtu` - Parallel result

These files contain identical mesh data but were generated with different processing approaches.

### Viewing Results
```bash
# Using ParaView (recommended)
paraview demo_mesh_serial_tiling.vtu

# Or with other VTU-compatible viewers
```

## Performance Profiling

### Serial Profiling
```python
import cProfile
import pstats

cProfile.run('exec(open("demo_tiling_serial.py").read())', 'stats_serial')
stats = pstats.Stats('stats_serial')
stats.sort_stats('cumulative').print_stats(20)
```

### Parallel Profiling
```python
import cProfile
import pstats

cProfile.run('exec(open("demo_tiling_parallel.py").read())', 'stats_parallel')
stats = pstats.Stats('stats_parallel')
stats.sort_stats('cumulative').print_stats(20)
```

## System Requirements

- **Processor**: 2+ cores (benefits increase with more cores)
- **RAM**: 2GB minimum, 4GB+ recommended
- **Disk**: 100MB free space for temporary files
- **Network**: Internet connection for data download
- **Python**: 3.10+

## Key Insights

1. **Parallelization is Most Effective for Mesh Building**
   - The mesh step dominates runtime (60-70%)
   - This is where parallel processing saves the most time

2. **Memory Usage Scales with Worker Count**
   - Serial: 800-1200 MB
   - Parallel (8 cores): 1500-2000 MB
   - Memory-limited systems should reduce `n_workers`

3. **Disk I/O Can Be a Bottleneck**
   - Temporary HDF5 files are created for intermediate results
   - SSD significantly faster than HDD for this workflow

4. **Optimal Settings Depend on Hardware**
   - Run both demos on your target hardware
   - Compare results to find best configuration for your system

## Further Reading

- See `SHARED_MEMORY_MESH_BUILDING.md` for comprehensive documentation
- Check `dtcc_core/builder/meshing/tiled_mesh_builder.py` for implementation details
- Review `dtcc_core/builder/meshing/shared_memory_backend.py` for storage mechanism

## Citation

If you use these demos in your research or publications:

```bibtex
@software{dtcc_core_tiling_2025,
  title={DTCC Core: Parallel Tiled Mesh Building Demonstrations},
  author={DTCC Contributors},
  year={2025},
  url={https://github.com/dtcc-platform/dtcc-core},
  note={Demo scripts for serial vs parallel processing comparison}
}
```

## Support

For issues or questions:
1. Check this README's troubleshooting section
2. Review `SHARED_MEMORY_MESH_BUILDING.md`
3. Check DTCC documentation
4. Open an issue on the DTCC GitHub repository
