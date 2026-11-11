// Offline experiments. Each writes CSVs to the output directory for scripts/plot.py.
//   pg_experiments <search|dstar|hybrid|rrt|track|all> [--out results]
#include <chrono>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "pg/dstar_lite.hpp"
#include "pg/hybrid_astar.hpp"
#include "pg/rrt_star.hpp"
#include "pg/search.hpp"
#include "pg/tracking.hpp"

using namespace pg;
using Clock = std::chrono::steady_clock;

namespace {

double ms_since(Clock::time_point t0) {
  return std::chrono::duration<double, std::milli>(Clock::now() - t0).count();
}

Cell random_free(const Grid& g, std::mt19937_64& rng) {
  std::uniform_int_distribution<int> X(0, g.width() - 1), Y(0, g.height() - 1);
  for (;;) {
    const Cell c{X(rng), Y(rng)};
    if (g.free(c.x, c.y)) return c;
  }
}

// ---------------------------------------------------------------- 1. optimal grid search
int exp_search(const std::string& out) {
  std::FILE* f = std::fopen((out + "/search_compare.csv").c_str(), "w");
  std::fprintf(f, "map,size,algo,found,cost,expansions,time_ms,cost_ratio\n");
  std::mt19937_64 rng(42);
  int mismatches = 0, bound_violations = 0, n_queries = 0;
  for (int m = 0; m < 60; ++m) {
    const int size = m < 30 ? 256 : 512;
    const Grid g = Grid::random_rects(size, size, 1.0, size / 2, 2, size / 12, 1000 + m);
    const GridGraph8 G(g);
    const Cell s = random_free(g, rng), t = random_free(g, rng);
    auto t0 = Clock::now();
    const auto dj = astar(G, s, t, ZeroHeuristic{});
    const double t_dj = ms_since(t0);
    if (!dj.found) continue;
    ++n_queries;
    auto row = [&](const char* algo, const SearchResult<Cell>& r, double ms) {
      std::fprintf(f, "%d,%d,%s,%d,%.6f,%ld,%.4f,%.6f\n", m, size, algo, r.found, r.cost, r.expansions, ms,
                   r.cost / dj.cost);
    };
    row("dijkstra", dj, t_dj);
    t0 = Clock::now();
    const auto as = astar(G, s, t, OctileHeuristic{});
    row("astar", as, ms_since(t0));
    for (double w : {1.5, 3.0}) {
      t0 = Clock::now();
      const auto wa = astar(G, s, t, OctileHeuristic{}, w);
      row(w == 1.5 ? "wastar_1.5" : "wastar_3", wa, ms_since(t0));
      if (wa.cost > w * dj.cost + 1e-9) ++bound_violations;
    }
    t0 = Clock::now();
    const auto jp = jps(g, s, t);
    row("jps", jp, ms_since(t0));
    bool valid = false;
    const double jc = path_cost(g, jp.path, &valid);
    if (std::abs(as.cost - dj.cost) > 1e-9 || std::abs(jp.cost - dj.cost) > 1e-9 || !valid ||
        std::abs(jc - jp.cost) > 1e-9)
      ++mismatches;
  }
  std::fclose(f);
  std::printf("search: %d queries, %d optimality mismatches, %d weighted-A* bound violations\n", n_queries,
              mismatches, bound_violations);
  return mismatches + bound_violations;
}

// ---------------------------------------------------------------- 2. replanning in unknown maps
// The robot knows only the map border. It senses a disc of radius 6 cells each step, and
// replans either incrementally (D* Lite) or from scratch with A* whenever the belief changed.
int exp_dstar(const std::string& out) {
  std::FILE* f = std::fopen((out + "/dstar.csv").c_str(), "w");
  std::fprintf(f, "map,size,steps,replans,dstar_expansions,astar_expansions,dstar_ms,astar_ms,path_len,mismatch\n");
  std::mt19937_64 rng(7);
  const int R = 6;
  int total_mismatch = 0;
  for (int m = 0; m < 30; ++m) {
    const int S = m < 15 ? 160 : 400;
    const Grid truth = Grid::random_rects(S, S, 1.0, S * S / 365, 3, 18, 500 + m);
    Cell start = random_free(truth, rng), goal = random_free(truth, rng);
    if (!astar(GridGraph8(truth), start, goal, OctileHeuristic{}).found) { --m; continue; }
    Grid belief(S, S, 1.0);
    DStarLite ds(belief, start, goal);
    auto t0 = Clock::now();
    ds.compute_shortest_path();
    double ds_ms = ms_since(t0), as_ms = 0.0;
    long as_exp = 0;
    {
      t0 = Clock::now();
      as_exp += astar(GridGraph8(ds.belief()), start, goal, OctileHeuristic{}).expansions;
      as_ms += ms_since(t0);
    }
    Cell pos = start;
    int steps = 0, replans = 0, mismatch = 0;
    double len = 0.0;
    std::FILE* trace = m == 0 ? std::fopen((out + "/dstar_trace.csv").c_str(), "w") : nullptr;
    if (trace) std::fprintf(trace, "x,y\n%d,%d\n", pos.x, pos.y);
    while (pos != goal && steps < 4 * S * S) {
      std::vector<std::pair<Cell, bool>> changes;
      for (int dy = -R; dy <= R; ++dy)
        for (int dx = -R; dx <= R; ++dx) {
          if (dx * dx + dy * dy > R * R) continue;
          const Cell c{pos.x + dx, pos.y + dy};
          if (!truth.inside(c.x, c.y)) continue;
          if (truth.occupied(c.x, c.y) != ds.belief().occupied(c.x, c.y)) changes.push_back({c, truth.occupied(c.x, c.y)});
        }
      if (!changes.empty()) {
        ++replans;
        t0 = Clock::now();
        ds.update_cells(changes);
        ds.compute_shortest_path();
        ds_ms += ms_since(t0);
        t0 = Clock::now();
        const auto a = astar(GridGraph8(ds.belief()), pos, goal, OctileHeuristic{});
        as_ms += ms_since(t0);
        as_exp += a.expansions;
        if (std::abs(a.cost - ds.cost_to_goal()) > 1e-9) ++mismatch;  // incremental == from scratch
      }
      const Cell nxt = ds.next_step();
      if (nxt == pos) break;
      len += (nxt.x != pos.x && nxt.y != pos.y) ? kSqrt2 : 1.0;
      pos = nxt;
      ds.move_to(pos);
      ++steps;
      if (trace) std::fprintf(trace, "%d,%d\n", pos.x, pos.y);
    }
    if (trace) {
      std::fclose(trace);
      truth.save_pgm(out + "/dstar_truth.pgm");
      ds.belief().save_pgm(out + "/dstar_belief.pgm");
      std::FILE* sg = std::fopen((out + "/dstar_endpoints.csv").c_str(), "w");
      std::fprintf(sg, "sx,sy,gx,gy\n%d,%d,%d,%d\n", start.x, start.y, goal.x, goal.y);
      std::fclose(sg);
    }
    total_mismatch += mismatch + (pos != goal);
    std::fprintf(f, "%d,%d,%d,%d,%ld,%ld,%.3f,%.3f,%.2f,%d\n", m, S, steps, replans, ds.expansions(), as_exp, ds_ms,
                 as_ms, len, mismatch + (pos != goal));
  }
  std::fclose(f);
  std::printf("dstar: 30 maps, %d mismatches/failures\n", total_mismatch);
  return total_mismatch;
}

// ---------------------------------------------------------------- 3. Hybrid A* parking
int exp_hybrid(const std::string& out) {
  const Grid lot = Grid::parking_lot(0.1);
  lot.save_pgm(out + "/parking_lot.pgm");
  const Vehicle veh;
  const Pose start{6.0, 30.0, 0.0};
  const Pose goal{39.25, 4.85, -M_PI / 2};  // nose-in to the empty stall, rear axle 4.85 m up
  struct Variant { const char* name; HybridConfig cfg; };
  std::vector<Variant> vs;
  HybridConfig base;
  vs.push_back({"full", base});
  HybridConfig nodub = base;
  nodub.use_dubins_heuristic = false;
  vs.push_back({"no_dubins_heuristic", nodub});
  HybridConfig noholo = base;
  noholo.use_holonomic_heuristic = false;
  vs.push_back({"no_holonomic_heuristic", noholo});
  HybridConfig noshot = base;
  noshot.analytic_expansion = false;
  vs.push_back({"no_analytic_expansion", noshot});
  HybridConfig fwd = base;
  fwd.allow_reverse = false;
  vs.push_back({"forward_only", fwd});
  std::FILE* f = std::fopen((out + "/hybrid_metrics.csv").c_str(), "w");
  std::fprintf(f, "variant,found,expansions,time_ms,length_m,direction_switches,analytic_finish\n");
  for (const auto& v : vs) {
    const HybridResult r = hybrid_astar(lot, veh, start, goal, v.cfg);
    std::fprintf(f, "%s,%d,%ld,%.1f,%.2f,%d,%d\n", v.name, r.found, r.expansions, r.time_ms, r.length,
                 r.direction_switches, r.analytic_finish);
    std::printf("hybrid %-24s found=%d expansions=%6ld time=%7.1f ms length=%.1f m switches=%d\n", v.name, r.found,
                r.expansions, r.time_ms, r.length, r.direction_switches);
    if (std::string(v.name) == "full" && r.found) {
      std::FILE* p = std::fopen((out + "/hybrid_path.csv").c_str(), "w");
      std::fprintf(p, "x,y,th,dir\n");
      for (std::size_t i = 0; i < r.path.size(); ++i)
        std::fprintf(p, "%.3f,%.3f,%.4f,%d\n", r.path[i].x, r.path[i].y, r.path[i].th, r.dir[i]);
      std::fclose(p);
    }
  }
  std::fclose(f);
  return 0;
}

// ---------------------------------------------------------------- 4. RRT* vs Informed RRT*
int exp_rrt(const std::string& out) {
  const Grid world = Grid::random_rects(400, 400, 0.25, 60, 8, 40, 77);
  const Grid cs = world.inflate(0.5);
  cs.save_pgm(out + "/rrt_cspace.pgm");
  const Vec2 start{3.0, 3.0}, goal{96.0, 96.0};
  const auto ref = astar(GridGraph8(cs), cs.to_cell(start.x, start.y), cs.to_cell(goal.x, goal.y), OctileHeuristic{});
  std::FILE* f = std::fopen((out + "/rrt_history.csv").c_str(), "w");
  std::fprintf(f, "variant,seed,iter,cost\n");
  std::FILE* s = std::fopen((out + "/rrt_summary.csv").c_str(), "w");
  std::fprintf(s, "variant,seed,found,first_iter,final_cost,nodes,time_ms\n");
  for (int informed = 0; informed <= 1; ++informed) {
    for (int seed = 1; seed <= 20; ++seed) {
      RrtConfig c;
      c.iterations = 6000;
      c.informed = informed;
      c.seed = static_cast<std::uint64_t>(seed);
      const RrtResult r = rrt_star(cs, start, goal, c);
      const char* name = informed ? "informed_rrt_star" : "rrt_star";
      for (const auto& [it, cost] : r.history) std::fprintf(f, "%s,%d,%d,%.4f\n", name, seed, it, cost);
      std::fprintf(s, "%s,%d,%d,%d,%.4f,%d,%.1f\n", name, seed, r.found, r.first_solution_iter, r.best_cost, r.nodes,
                   r.time_ms);
      if (seed == 1) {
        std::FILE* p = std::fopen((out + "/rrt_path_" + name + ".csv").c_str(), "w");
        std::fprintf(p, "x,y\n");
        for (const auto& q : r.path) std::fprintf(p, "%.3f,%.3f\n", q.x, q.y);
        std::fclose(p);
      }
    }
  }
  std::fclose(f);
  std::fclose(s);
  std::FILE* rf = std::fopen((out + "/rrt_reference.csv").c_str(), "w");
  std::fprintf(rf, "astar_grid_cost_m\n%.4f\n", ref.cost * cs.resolution());
  std::fclose(rf);
  std::printf("rrt: done (8-connected A* reference %.2f m)\n", ref.cost * cs.resolution());
  return 0;
}

// ---------------------------------------------------------------- 5. tracking with estimation
int exp_track(const std::string& out) {
  const RefPath path = RefPath::road(20.0, 2.0, 1.5);
  std::FILE* pf = std::fopen((out + "/road.csv").c_str(), "w");
  std::fprintf(pf, "s,x,y,th,kappa,v\n");
  for (std::size_t i = 0; i < path.s.size(); ++i)
    std::fprintf(pf, "%.2f,%.3f,%.3f,%.5f,%.5f,%.3f\n", path.s[i], path.x[i], path.y[i], path.th[i], path.kappa[i], path.v[i]);
  std::fclose(pf);

  std::FILE* f = std::fopen((out + "/tracking_metrics.csv").c_str(), "w");
  std::fprintf(f, "tracker,state_source,seed,completed,rms_lat_m,max_lat_m,rms_head_deg,steer_rate_rms,"
                  "rms_pos_err_m,max_pos_err_m,gps_outliers,gps_rejected\n");
  for (const char* tr : {"pure_pursuit", "stanley", "lqr"}) {
    for (int mode = 0; mode < 3; ++mode) {
      const char* src = mode == 0 ? "truth" : mode == 1 ? "ekf" : "ekf_no_gate";
      for (int seed = 1; seed <= 20; ++seed) {
        ScenarioConfig c;
        c.tracker = tr;
        c.use_ekf = mode != 0;
        c.gps_gate = mode == 1;
        c.seed = static_cast<std::uint64_t>(seed);
        if (seed == 1) c.log_path = out + "/track_" + tr + "_" + src + ".csv";
        const ScenarioResult r = run_scenario(path, c);
        std::fprintf(f, "%s,%s,%d,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%ld,%ld\n", tr, src, seed, r.completed, r.rms_lat,
                     r.max_lat, r.rms_head * 180.0 / M_PI, r.steer_rate_rms, r.rms_pos_err, r.max_pos_err,
                     r.gps_outliers, r.gps_rejected);
      }
    }
  }
  std::fclose(f);

  // Filter consistency: 50 Monte Carlo runs, NEES per 10 Hz epoch averaged across runs.
  std::vector<std::vector<double>> nees;
  std::vector<double> nis_all;
  for (int seed = 100; seed < 150; ++seed) {
    ScenarioConfig c;
    c.tracker = "lqr";
    c.seed = static_cast<std::uint64_t>(seed);
    c.gps_outlier_prob = 0.0;  // consistency is about the noise model, not the outlier model
    const auto r = run_scenario(path, c);
    nees.push_back(r.nees);
    nis_all.insert(nis_all.end(), r.nis.begin(), r.nis.end());
  }
  std::size_t T = nees[0].size();
  for (const auto& v : nees) T = std::min(T, v.size());
  std::FILE* nf = std::fopen((out + "/ekf_nees.csv").c_str(), "w");
  std::fprintf(nf, "t,avg_nees,runs\n");
  for (std::size_t k = 0; k < T; ++k) {
    double a = 0;
    for (const auto& v : nees) a += v[k];
    std::fprintf(nf, "%.1f,%.4f,%zu\n", 0.1 * static_cast<double>(k), a / static_cast<double>(nees.size()), nees.size());
  }
  std::fclose(nf);
  std::FILE* sf = std::fopen((out + "/ekf_nis.csv").c_str(), "w");
  std::fprintf(sf, "nis\n");
  for (double v : nis_all) std::fprintf(sf, "%.4f\n", v);
  std::fclose(sf);
  std::printf("track: done\n");
  return 0;
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 2) {
    std::fprintf(stderr, "usage: %s <search|dstar|hybrid|rrt|track|all> [--out DIR]\n", argv[0]);
    return 2;
  }
  const std::string what = argv[1];
  std::string out = "results";
  if (argc >= 4 && std::string(argv[2]) == "--out") out = argv[3];
  int rc = 0;
  const auto t0 = Clock::now();
  if (what == "search" || what == "all") rc |= exp_search(out);
  if (what == "dstar" || what == "all") rc |= exp_dstar(out);
  if (what == "hybrid" || what == "all") rc |= exp_hybrid(out);
  if (what == "rrt" || what == "all") rc |= exp_rrt(out);
  if (what == "track" || what == "all") rc |= exp_track(out);
  std::printf("[%s] %.1f s\n", what.c_str(), ms_since(t0) / 1000.0);
  return rc;
}
