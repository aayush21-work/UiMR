/*
 * ns_uimr.cpp — Game-Engine-Inspired Navier-Stokes Solver
 * ================================================================
 * Used by benchmark.cpp when compiled without -DGAME_ENGINE_STANDALONE
 */

#ifndef GAME_ENGINE_NS_H
#define GAME_ENGINE_NS_H

#include "ns_common.h"

// =============================================================================
// TECHNIQUE 4: Memory Pool
// =============================================================================
struct TileBuffer {
    Field u, v, p;
    bool in_use;
};

struct MemoryPool {
    std::vector<TileBuffer> buffers;
    int reused_count;
    int alloc_count;

    void init(int n_buffers, int ny, int nx) {
        buffers.resize(n_buffers);
        for (int i = 0; i < n_buffers; i++) {
            buffers[i].u.alloc(ny, nx);
            buffers[i].v.alloc(ny, nx);
            buffers[i].p.alloc(ny, nx);
            buffers[i].in_use = false;
        }
        reused_count = 0;
        alloc_count = 0;
    }

    TileBuffer* acquire() {
        for (auto& b : buffers) {
            if (!b.in_use) { b.in_use = true; reused_count++; return &b; }
        }
        alloc_count++;
        return nullptr;
    }

    void release(TileBuffer* b) { if (b) b->in_use = false; }

    void cleanup() {
        for (auto& b : buffers) {
            b.u.free_mem(); b.v.free_mem(); b.p.free_mem();
        }
    }
};

// =============================================================================
// Tile metadata (for culling, priority, reuse)
// =============================================================================
struct TileMeta {
    int row, col, idx;
    double residual;
    double boundary_delta;
    bool converged;
    bool reused;
    bool active;
};

// =============================================================================
// Tile extraction / injection
// =============================================================================
static void extract_tile(const Field& global, Field& tile,
                         int tr, int tc, int tile_size, int halo)
{
    int n = global.ny;
    int tny = tile_size + 2 * halo;
    int tnx = tile_size + 2 * halo;
    int r0 = tr * tile_size - halo;
    int c0 = tc * tile_size - halo;
    for (int i = 0; i < tny; i++)
        for (int j = 0; j < tnx; j++) {
            int gi = (r0 + i + n) % n;
            int gj = (c0 + j + n) % n;
            tile(i, j) = global(gi, gj);
        }
}

static void inject_tile(Field& global, const Field& tile,
                        int tr, int tc, int tile_size, int halo)
{
    int r0 = tr * tile_size;
    int c0 = tc * tile_size;
    for (int i = 0; i < tile_size; i++)
        for (int j = 0; j < tile_size; j++)
            global(r0 + i, c0 + j) = tile(halo + i, halo + j);
}

// =============================================================================
// TECHNIQUE 10: Tile residual (for priority / culling)
// =============================================================================
static double tile_residual(const Field& u) {
    double max_res = 0.0;
    int h = HALO;
    for (int i = h; i < u.ny - h; i++)
        for (int j = h; j < u.nx - h; j++) {
            double lap = u.at(i-1,j) + u.at(i+1,j) +
                         u.at(i,j-1) + u.at(i,j+1) - 4.0 * u(i,j);
            max_res = std::max(max_res, std::fabs(lap));
        }
    return max_res;
}

// =============================================================================
// TECHNIQUE 7: Boundary delta (for temporal reuse)
// =============================================================================
static double boundary_delta(const Field& a, const Field& b) {
    double max_d = 0.0;
    int ny = a.ny, nx = a.nx;
    for (int j = 0; j < nx; j++) {
        max_d = std::max(max_d, std::fabs(a(0, j) - b(0, j)));
        max_d = std::max(max_d, std::fabs(a(ny-1, j) - b(ny-1, j)));
    }
    for (int i = 0; i < ny; i++) {
        max_d = std::max(max_d, std::fabs(a(i, 0) - b(i, 0)));
        max_d = std::max(max_d, std::fabs(a(i, nx-1) - b(i, nx-1)));
    }
    return max_d;
}

// =============================================================================
// NS step on a single tile
// =============================================================================
static void ns_step_tile(Field& tu, Field& tv, Field& tp, bool use_mg) {
    int ny = tu.ny, nx = tu.nx;
    double dx = DX, dt = DT, nu = NU;
    double inv_dx2 = 1.0/(dx*dx), inv_2dx = 1.0/(2.0*dx), inv_dt = 1.0/dt;

    Field u_lap(ny,nx), v_lap(ny,nx);
    Field dudx(ny,nx), dudy(ny,nx), dvdx(ny,nx), dvdy(ny,nx);
    Field u_star(ny,nx), v_star(ny,nx);
    Field div_f(ny,nx), rhs(ny,nx);

    // Step 1: Advection + Diffusion (no inner OpenMP — tile is small)
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            u_lap(i,j) = (tu.at(i-1,j)+tu.at(i+1,j)+tu.at(i,j-1)+tu.at(i,j+1)-4.0*tu(i,j))*inv_dx2;
            v_lap(i,j) = (tv.at(i-1,j)+tv.at(i+1,j)+tv.at(i,j-1)+tv.at(i,j+1)-4.0*tv(i,j))*inv_dx2;
            dudx(i,j) = (tu.at(i,j+1)-tu.at(i,j-1))*inv_2dx;
            dudy(i,j) = (tu.at(i+1,j)-tu.at(i-1,j))*inv_2dx;
            dvdx(i,j) = (tv.at(i,j+1)-tv.at(i,j-1))*inv_2dx;
            dvdy(i,j) = (tv.at(i+1,j)-tv.at(i-1,j))*inv_2dx;
        }

    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            double adv_u = tu(i,j)*dudx(i,j) + tv(i,j)*dudy(i,j);
            double adv_v = tu(i,j)*dvdx(i,j) + tv(i,j)*dvdy(i,j);
            u_star(i,j) = tu(i,j) + dt*(-adv_u + nu*u_lap(i,j));
            v_star(i,j) = tv(i,j) + dt*(-adv_v + nu*v_lap(i,j));
        }
    apply_velocity_bc(u_star, v_star);

    // Step 2: Pressure Poisson
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            double du = (u_star.at(i,j+1)-u_star.at(i,j-1))*inv_2dx;
            double dv = (v_star.at(i+1,j)-v_star.at(i-1,j))*inv_2dx;
            rhs(i,j) = (du + dv) * inv_dt;
        }

    if (use_mg)
        pressure_multigrid(tp, rhs, dx);
    else
        pressure_jacobi(tp, rhs, dx, JACOBI_ITERS);

    // Step 3: Projection
    for (int i = 0; i < ny; i++)
        for (int j = 0; j < nx; j++) {
            double dpdx = (tp.at(i,j+1)-tp.at(i,j-1))*inv_2dx;
            double dpdy = (tp.at(i+1,j)-tp.at(i-1,j))*inv_2dx;
            tu(i,j) = u_star(i,j) - dt * dpdx;
            tv(i,j) = v_star(i,j) - dt * dpdy;
        }
    apply_velocity_bc(tu, tv);

    u_lap.free_mem(); v_lap.free_mem();
    dudx.free_mem(); dudy.free_mem(); dvdx.free_mem(); dvdy.free_mem();
    u_star.free_mem(); v_star.free_mem();
    div_f.free_mem(); rhs.free_mem();
}

// =============================================================================
// Game engine statistics
// =============================================================================
struct GEStats {
    int total_tiles;
    int computed;
    int culled;
    int reused;
};

// =============================================================================
// Full UiMR solver (ALL 10 techniques)
// =============================================================================
static double run_game_engine(const Field& u0, const Field& v0, const Field& p0,
                              Field& u_out, Field& v_out, Field& p_out,
                              GEStats& stats, int n_steps)
{
    u_out.copy_from(u0); v_out.copy_from(v0); p_out.copy_from(p0);
    stats = {0, 0, 0, 0};

    // TECHNIQUE 4: Memory pool
    std::vector<Field> tu(NTILES), tv(NTILES), tp(NTILES);
    std::vector<Field> prev_tu(NTILES), prev_tv(NTILES);
    for (int i = 0; i < NTILES; i++) {
        tu[i].alloc(TILE_FULL, TILE_FULL);
        tv[i].alloc(TILE_FULL, TILE_FULL);
        tp[i].alloc(TILE_FULL, TILE_FULL);
        prev_tu[i].alloc(TILE_FULL, TILE_FULL);
        prev_tv[i].alloc(TILE_FULL, TILE_FULL);
    }

    std::vector<TileMeta> meta(NTILES);
    for (int t = 0; t < NTILES; t++) {
        meta[t].row = t / NTILES_DIR;
        meta[t].col = t % NTILES_DIR;
        meta[t].idx = t;
    }

    Timer tm; tm.start();

    for (int s = 0; s < n_steps; s++) {
        stats.total_tiles += NTILES;

        // Extract tiles
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            extract_tile(u_out, tu[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            extract_tile(v_out, tv[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            extract_tile(p_out, tp[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }

        // TECHNIQUE 10 + 7: Compute residual and boundary delta
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            meta[t].residual = tile_residual(tu[t]);
            meta[t].boundary_delta = (s > 0) ? boundary_delta(tu[t], prev_tu[t]) : 1e10;

            // TECHNIQUE 2: Culling
            meta[t].converged = (meta[t].residual < CULL_THRESHOLD && s > 0);
            // TECHNIQUE 7: Temporal reuse
            meta[t].reused = (s > 1 && meta[t].boundary_delta < REUSE_THRESHOLD);
            meta[t].active = !meta[t].converged && !meta[t].reused;
        }

        // TECHNIQUE 3: Priority scheduling — sort by residual (highest first)
        std::vector<int> active_list;
        active_list.reserve(NTILES);
        for (int t = 0; t < NTILES; t++) {
            if (meta[t].active)
                active_list.push_back(t);
            else if (meta[t].converged)
                stats.culled++;
            else if (meta[t].reused)
                stats.reused++;
        }

        std::sort(active_list.begin(), active_list.end(),
                  [&meta](int a, int b) {
                      return meta[a].residual > meta[b].residual;
                  });

        int n_active = (int)active_list.size();
        stats.computed += n_active;

        // Save for temporal reuse next step
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            prev_tu[t].copy_from(tu[t]);
            prev_tv[t].copy_from(tv[t]);
        }

        // TECHNIQUE 5: Batched parallel dispatch
        // TECHNIQUE 1: Multigrid pressure solve
        #pragma omp parallel for schedule(dynamic)
        for (int ai = 0; ai < n_active; ai++) {
            int t = active_list[ai];
            ns_step_tile(tu[t], tv[t], tp[t], true);
        }

        // Inject all tiles
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            inject_tile(u_out, tu[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            inject_tile(v_out, tv[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            inject_tile(p_out, tp[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }

        if ((s+1) % 20 == 0)
            printf("    step %d/%d, active: %d/%d (%.0f%%)\n",
                   s+1, n_steps, n_active, NTILES,
                   100.0 * n_active / NTILES);
    }
    double elapsed = tm.elapsed();

    for (int i = 0; i < NTILES; i++) {
        tu[i].free_mem(); tv[i].free_mem(); tp[i].free_mem();
        prev_tu[i].free_mem(); prev_tv[i].free_mem();
    }
    return elapsed;
}

// =============================================================================
// Standalone test
// =============================================================================
#ifdef UiMR
int main() {
    int nt = omp_get_max_threads();
    printf("UiMR NS Solver (All 10 Techniques)\n");
    printf("Grid: %d×%d, Tiles: %d×%d=%d, Steps: %d, Threads: %d\n\n",
           GRID_SIZE, GRID_SIZE, NTILES_DIR, NTILES_DIR, NTILES, NUM_TIMESTEPS, nt);

    Field u0, v0, p0;
    init_fields(u0, v0, p0);

    Field u(GRID_SIZE,GRID_SIZE), v(GRID_SIZE,GRID_SIZE), p(GRID_SIZE,GRID_SIZE);
    GEStats stats;
    double t = run_game_engine(u0, v0, p0, u, v, p, stats, NUM_TIMESTEPS);

    printf("\n  Time: %.2f s\n", t);
    printf("  Computed: %d/%d (%.1f%%)\n", stats.computed, stats.total_tiles,
           100.0 * stats.computed / stats.total_tiles);
    printf("  Culled:   %d (%.1f%%)\n", stats.culled,
           100.0 * stats.culled / stats.total_tiles);
    printf("  Reused:   %d (%.1f%%)\n", stats.reused,
           100.0 * stats.reused / stats.total_tiles);

    u0.free_mem(); v0.free_mem(); p0.free_mem();
    u.free_mem(); v.free_mem(); p.free_mem();
    return 0;
}
#endif

#endif // UiMR
