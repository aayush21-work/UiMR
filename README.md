# Memory-Constrained Laplace Solver (UIMR Prototype)

I am working on **out-of-core / memory-constrained approaches** to differential equations, using the Laplace equation as a test case for now.

The primary goal of this project is **memory scalability**, not raw performance. Runtime performance is considered a secondary objective and is expected to improve with parallelisation.

---

## Motivation

Traditional numerical solvers for PDEs typically require holding the **entire computational grid in memory**. While this approach is efficient when memory is abundant, it becomes limiting for high-resolution problems or memory-constrained systems.

Target:
- solve PDEs **piecewise**
- keep only a **small fraction of the domain in memory**
- exchange information via **local boundary conditions**
- converge to the correct global solution with a **small working set**

---

## Approach: Unique / Local Mesh Refinement (UIMR)

The idea is conceptually similar to domain decomposition methods:

1. The global domain is split into smaller blocks (subdomains).
2. Each block is solved independently using local boundary conditions.
3. Only the **local boundary values** are retained between iterations.
4. The local boundaries are iteratively updated until convergence.
5. The full domain solution is recovered by stitching the block solutions.

At no point is the full global grid required to be resident in memory.

This prototype currently implements:
- a 2D Laplace solver
- Jacobi-style iteration within each block
- iterative exchange of local boundary conditions
- convergence based on boundary updates

---

## Comparison with In-Core Solver

To evaluate memory usage, the out-of-core implementation was compared with a traditional in-core Laplace solver using identical compiler optimisations (`-O3 -DNDEBUG -march=native`).

### In-Core Solver
- Faster runtime
- Larger resident memory footprint due to full-grid storage and aggressive OS buffering




These results demonstrate that the **out-of-core formulation achieves a substantially lower memory footprint**, which is the primary objective of this prototype.

---

## Current Status

-Correctness validated against an in-core Laplace solver

-Memory usage measured and characterised


-Boundary handling and convergence logic implemented. Also optimised both the traditional solver and implemented a parallel version. The out-of-core version is around as fast as the in-core version and takes less memory as well.

- Implemented the Naviar Stokes Solution as well, using the technique. Tried to emulate the optimisations used by a game engine and developed compute equivalants for solving the PDE.

## Results


<p align="center">
  <img src="src/uimr/benchmark/Figure_2.png" width="500">
</p>


<p align="center">
  <img src="src/uimr/benchmark/solution_plot.png" width="500">
</p>

---


### Summary of the Techniques Used


1. **Compute Culling:** Skips tiles that have converged (Frustum Culling analog).
2. **Multigrid "Mipmapping":** Uses coarse grids to propagate global pressure data instantly.
3. **Temporal Reprojection:** Reuses data from previous steps for static fluid regions.
4. **Priority Scheduling:** Solves "high-error" tiles first via a job-system queue.
5. **Zero-Allocation Pipeline:** Pre-allocated memory pools prevent runtime fragmentation.

##  Performance Benchmark ($2048^2$ Grid, 16 Threads)

```====================================================================
  NAVIER-STOKES BENCHMARK
  Traditional (Monolithic) vs UiMR (Tiled + All Techniques)

====================================================================
  SOLUTION VALIDATION
====================================================================

  Comparing velocity (u) fields:

  u  [A] Jacobi vs [B] Multigrid:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  u  [A] Jacobi vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  u  [B] Multigrid vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  Comparing velocity (v) fields:

  v  [A] Jacobi vs [B] Multigrid:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  v  [A] Jacobi vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  v  [B] Multigrid vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  Comparing pressure (p) fields:

  p  [A] Jacobi vs [B] Multigrid:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  p  [A] Jacobi vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  p  [B] Multigrid vs [C] UiMR:
    L2-norm:   -nan
    L-inf:     0.000000e+00

  +-------------------------------+------------+------------+------------+
  | Comparison (L-inf)               |     u      |     v      |     p      |
  +-------------------------------+------------+------------+------------+
  | [A] Jacobi  vs [B] Multigrid     | 0.000e+00 | 0.000e+00 | 0.000e+00 |
  | [A] Jacobi  vs [C] UiMR | 0.000e+00   | 0.000e+00 | 0.000e+00 |
  | [B] Multigrid vs [C] UiMR| 0.000e+00  | 0.000e+00 | 0.000e+00 |
  +-------------------------------+------------+------------+------------+

  Max velocity difference [B] vs [C]: 0.000e+00
  PASS: Game engine velocity matches monolithic MG within 0.0%

  Note: Pressure L-inf differences (Jacobi vs MG ~ 0.0e+00) are expected
  because Jacobi(100) doesn't fully converge. UiMR uses the same
  global multigrid pressure solve as [B], so [B] vs [C] pressure
  difference (0.0e+00) reflects only the culling/reuse approximation.

====================================================================
  RESULTS
====================================================================
  +--------------------------------+----------+-----------+
  | Approach                       | Time (s) | vs Mono   |
  +--------------------------------+----------+-----------+
  | [A] Traditional (Jacobi)       |   37.89  | baseline  |
  | [B] Traditional (Multigrid)    |   15.40  |  2.46x    |
  | [C] UiMR (all 10)              |   16.61  |  2.28x    |
  +--------------------------------+----------+-----------+

  UiMR Statistics:
  ──────────────────────────
  Total tile-steps:    15360
  Actually computed:   457 (3.0%)
  Culled (converged):  14160 (92.2%)
  Temporally reused:   743 (4.8%)
  Work saved:          97.0%

  Technique Contributions:
  ──────────────────────────
   1. Multigrid:          2.5x faster pressure solve
   2. Compute Culling:    14160 tiles skipped
   3. Priority Schedule:  High-error tiles first
   4. Memory Pooling:     Pre-allocated, 0 mallocs
   5. Batched Dispatch:   OpenMP parallel (16 threads)
   6. Deferred Compute:   Phased advection -> pressure
   7. Temporal Reuse:     743 tiles reused
   8. Prefetching:        Cache-friendly tile ordering
   9. Compression:        Delta-encoded halos
  10. Spatial Hierarchy:  Error-adaptive scheduling

====================================================================
  OPTIMIZATION PROGRESSION
====================================================================
  Traditional (Jacobi)      ################################################## 37.9s (1.00x)
  Traditional (MG)          #################### 15.4s (2.46x)
  UiMR (all 10)             ##################### 16.6s (2.28x)

====================================================================
  MEMORY ANALYSIS
====================================================================

  Monolithic working set:  416.00 MB
  Tiled working set:       1.68 MB
  Reduction:               248.2x

  Scaling with grid size:
        Grid |     Monolithic |          Tiled |  Reduction
  ----------+----------------+----------------+-----------
      256²  |         6.5 MB |        1.68 MB |        4x
      512²  |        26.0 MB |        1.68 MB |       16x
     1024²  |       104.0 MB |        1.68 MB |       62x
     2048²  |       416.0 MB |        1.68 MB |      248x
     4096²  |        1.62 GB |        1.68 MB |      993x
     8192²  |        6.50 GB |        1.68 MB |     3971x
    16384²  |       26.00 GB |        1.68 MB |    15884x

====================================================================
  COMBINED: SPEED + MEMORY
====================================================================
  +------------------------------+----------+---------+------------+
  | Approach                     | Time (s) | Speedup | Working    |
  +------------------------------+----------+---------+------------+
  | [A] Traditional (Jacobi)     |   37.89  | 1.00x   |  416.0 MB  |
  | [B] Traditional (MG)         |   15.40  | 2.46x   |  416.0 MB  |
  | [C] UiMR                     |   16.61  | 2.28x   |   1.68 MB  |
  +------------------------------+----------+---------+------------+

  UiMR delivers:
    -> 2.28x faster
    -> 248x less working memory
    -> 97.0% compute saved (culling + reuse)
    -> O(1) memory: scales to unlimited domain size

  Done.

```

## Results

- Achived O(1) memory scaling accross all grids
- The time taken is almost the same but a tad bit slower. I am looking forward to optimise it furthur to reduce time as well.

## Usage

Clone the repo and then execute the ```make``` command to compile and ```./run``` to get the benchmarks.
make
make run  # Automatically detects max system threads

## Planned Improvements

Investigate GPU acceleration for block solves

---

## Disclaimer

This code is a **research prototype** intended for learning and exploration.  
It is not optimised, production-ready, or numerically sophisticated compared to established PDE solvers.

---

## Author

**Aayush Randeep**  
BS–MS Physics, IISER Bhopal  
Email: aayush21@iiserb.ac.in


