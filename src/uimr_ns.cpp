#ifndef GAME_ENGINE_NS_H
#define GAME_ENGINE_NS_H

#include "ns_common.h"

// =============================================================================
// Tile metadata
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
// Tile extraction / injection (clamped, not periodic)
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
            int gi = std::max(0, std::min(r0 + i, n - 1));
            int gj = std::max(0, std::min(c0 + j, n - 1));
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
// Tile residual for culling
// =============================================================================
static double tile_max_velocity(const Field& tu, const Field& tv) {
    double max_v = 0.0;
    int h = HALO;
    for (int i = h; i < tu.ny - h; i++)
        for (int j = h; j < tu.nx - h; j++)
            max_v = std::max(max_v, tu(i,j)*tu(i,j) + tv(i,j)*tv(i,j));
    return std::sqrt(max_v);
}

// =============================================================================
// Boundary delta for temporal reuse
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
// Advection on a single tile — updates interior only
// =============================================================================
static void advect_tile(const Field& tu, const Field& tv,
                        Field& tu_star, Field& tv_star)
{
    int ny = tu.ny, nx = tu.nx;
    int h = HALO;
    double dx = DX, dt = DT, nu = NU;
    double inv_dx2 = 1.0/(dx*dx), inv_2dx = 1.0/(2.0*dx);

    // Copy halos through
    tu_star.copy_from(tu);
    tv_star.copy_from(tv);

    for (int i = h; i < ny - h; i++)
        for (int j = h; j < nx - h; j++) {
            double u_lap = (tu(i-1,j)+tu(i+1,j)+tu(i,j-1)+tu(i,j+1)-4.0*tu(i,j))*inv_dx2;
            double v_lap = (tv(i-1,j)+tv(i+1,j)+tv(i,j-1)+tv(i,j+1)-4.0*tv(i,j))*inv_dx2;
            double dudx = (tu(i,j+1)-tu(i,j-1))*inv_2dx;
            double dudy = (tu(i+1,j)-tu(i-1,j))*inv_2dx;
            double dvdx = (tv(i,j+1)-tv(i,j-1))*inv_2dx;
            double dvdy = (tv(i+1,j)-tv(i-1,j))*inv_2dx;

            tu_star(i,j) = tu(i,j) + dt*(-tu(i,j)*dudx - tv(i,j)*dudy + nu*u_lap);
            tv_star(i,j) = tv(i,j) + dt*(-tu(i,j)*dvdx - tv(i,j)*dvdy + nu*v_lap);
        }
}

// =============================================================================
// Projection on a single tile — updates interior only
// Computes u
// =============================================================================
static void project_tile(Field& tu, Field& tv,
                         const Field& tu_star, const Field& tv_star,
                         const Field& tp)
{
    int ny = tu.ny, nx = tu.nx;
    int h = HALO;
    double dt = DT, inv_2dx = 1.0/(2.0*DX);

    for (int i = h; i < ny - h; i++)
        for (int j = h; j < nx - h; j++) {
            double dpdx = (tp(i,j+1)-tp(i,j-1))*inv_2dx;
            double dpdy = (tp(i+1,j)-tp(i-1,j))*inv_2dx;
            tu(i,j) = tu_star(i,j) - dt * dpdx;
            tv(i,j) = tv_star(i,j) - dt * dpdy;
        }
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
// Full Game Engine solver
//
// Each timestep:
//   1. Extract tiles from global u, v
//   2. Advect active tiles → u_star, v_star (tiled, parallel, with culling)
//   3. Inject u_star back to global
//   4. Global pressure Poisson solve (multigrid on full grid)
//   5. Extract p tiles
//   6. Project active tiles: u = u_star - dt∇p (tiled, parallel)
//   7. Inject u, v back to global; apply domain BCs
// =============================================================================
static double run_game_engine(const Field& u0, const Field& v0, const Field& p0,
                              Field& u_out, Field& v_out, Field& p_out,
                              GEStats& stats, int n_steps)
{
    int n = GRID_SIZE;
    u_out.copy_from(u0); v_out.copy_from(v0); p_out.copy_from(p0);
    stats = {0, 0, 0, 0};

    // Tile buffers
    std::vector<Field> tu(NTILES), tv(NTILES), tp(NTILES);
    std::vector<Field> tu_star(NTILES), tv_star(NTILES);
    std::vector<Field> prev_tu(NTILES), prev_tv(NTILES);
    for (int i = 0; i < NTILES; i++) {
        tu[i].alloc(TILE_FULL, TILE_FULL);
        tv[i].alloc(TILE_FULL, TILE_FULL);
        tp[i].alloc(TILE_FULL, TILE_FULL);
        tu_star[i].alloc(TILE_FULL, TILE_FULL);
        tv_star[i].alloc(TILE_FULL, TILE_FULL);
        prev_tu[i].alloc(TILE_FULL, TILE_FULL);
        prev_tv[i].alloc(TILE_FULL, TILE_FULL);
    }

    std::vector<TileMeta> meta(NTILES);
    for (int t = 0; t < NTILES; t++) {
        meta[t].row = t / NTILES_DIR;
        meta[t].col = t % NTILES_DIR;
        meta[t].idx = t;
    }

    // Global temporaries for pressure solve
    Field u_star_global(n, n), v_star_global(n, n);
    Field dudx_g(n, n), dvdy_g(n, n), rhs_g(n, n);

    Timer tm; tm.start();

    for (int s = 0; s < n_steps; s++) {
        stats.total_tiles += NTILES;

        // ── Phase 1: Extract u, v tiles ──
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            extract_tile(u_out, tu[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            extract_tile(v_out, tv[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }

        // ── Culling / reuse decision ──
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            meta[t].residual = tile_max_velocity(tu[t], tv[t]);
            meta[t].boundary_delta = (s > 0) ? boundary_delta(tu[t], prev_tu[t]) : 1e10;
            meta[t].converged = (meta[t].residual < CULL_THRESHOLD && s > 0);
            meta[t].reused = (s > 1 && meta[t].boundary_delta < REUSE_THRESHOLD);
            meta[t].active = !meta[t].converged && !meta[t].reused;
        }

        std::vector<int> active_list;
        active_list.reserve(NTILES);
        for (int t = 0; t < NTILES; t++) {
            if (meta[t].active) active_list.push_back(t);
            else if (meta[t].converged) stats.culled++;
            else if (meta[t].reused) stats.reused++;
        }

        std::sort(active_list.begin(), active_list.end(),
                  [&meta](int a, int b) { return meta[a].residual > meta[b].residual; });

        int n_active = (int)active_list.size();
        stats.computed += n_active;

        // Save for next step's reuse check
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            prev_tu[t].copy_from(tu[t]);
            prev_tv[t].copy_from(tv[t]);
        }

        // ── Phase 2: Tiled advection (only active tiles) ──
        #pragma omp parallel for schedule(dynamic)
        for (int ai = 0; ai < n_active; ai++) {
            int t = active_list[ai];
            advect_tile(tu[t], tv[t], tu_star[t], tv_star[t]);
        }
        // For culled/reused tiles, u_star = u (no change)
        for (int t = 0; t < NTILES; t++) {
            if (!meta[t].active) {
                tu_star[t].copy_from(tu[t]);
                tv_star[t].copy_from(tv[t]);
            }
        }

        // ── Phase 3: Inject u_star to global ──
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            inject_tile(u_star_global, tu_star[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            inject_tile(v_star_global, tv_star[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }
        apply_velocity_bc(u_star_global, v_star_global);

        // ── Phase 4: Global pressure solve (multigrid) ──
        compute_ddx(u_star_global, dudx_g, DX);
        compute_ddy(v_star_global, dvdy_g, DX);
        double inv_dt = 1.0 / DT;
        #pragma omp parallel for schedule(static) collapse(2)
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                rhs_g(i,j) = (dudx_g(i,j) + dvdy_g(i,j)) * inv_dt;

        pressure_multigrid(p_out, rhs_g, DX);

        // ── Phase 5: Extract p tiles + u_star tiles for projection ──
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            extract_tile(p_out, tp[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            extract_tile(u_star_global, tu_star[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            extract_tile(v_star_global, tv_star[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }

        // ── Phase 6: Tiled projection (only active tiles) ──
        #pragma omp parallel for schedule(dynamic)
        for (int ai = 0; ai < n_active; ai++) {
            int t = active_list[ai];
            project_tile(tu[t], tv[t], tu_star[t], tv_star[t], tp[t]);
        }
        // For culled/reused tiles, u = u_star (already set)
        for (int t = 0; t < NTILES; t++) {
            if (!meta[t].active) {
                tu[t].copy_from(tu_star[t]);
                tv[t].copy_from(tv_star[t]);
            }
        }

        // ── Phase 7: Inject u, v back to global + BCs ──
        #pragma omp parallel for schedule(static)
        for (int t = 0; t < NTILES; t++) {
            inject_tile(u_out, tu[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
            inject_tile(v_out, tv[t], meta[t].row, meta[t].col, TILE_SIZE, HALO);
        }
        apply_velocity_bc(u_out, v_out);

        if ((s+1) % 20 == 0)
            printf("    step %d/%d, active: %d/%d (%.0f%%)\n",
                   s+1, n_steps, n_active, NTILES,
                   100.0 * n_active / NTILES);
    }
    double elapsed = tm.elapsed();

    // Cleanup
    for (int i = 0; i < NTILES; i++) {
        tu[i].free_mem(); tv[i].free_mem(); tp[i].free_mem();
        tu_star[i].free_mem(); tv_star[i].free_mem();
        prev_tu[i].free_mem(); prev_tv[i].free_mem();
    }
    u_star_global.free_mem(); v_star_global.free_mem();
    dudx_g.free_mem(); dvdy_g.free_mem(); rhs_g.free_mem();
    return elapsed;
}

// =============================================================================
// Standalone test
// =============================================================================
#ifdef GAME_ENGINE_STANDALONE
int main() {
    int nt = omp_get_max_threads();
    printf("Game Engine NS Solver (All 10 Techniques)\n");
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

#endif // GAME_ENGINE_NS_H
