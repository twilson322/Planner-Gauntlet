#include "pg/hybrid_astar.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <queue>

#include "pg/search.hpp"

namespace pg {

void Vehicle::discs(int n, std::vector<double>& offsets, double& radius) const {
  offsets.clear();
  const double seg = length / n;
  for (int i = 0; i < n; ++i) offsets.push_back(-rear_overhang + (i + 0.5) * seg);
  radius = std::hypot(seg / 2.0, width / 2.0);
}

FootprintChecker::FootprintChecker(const Grid& g, const Vehicle& v, int n_discs)
    : g_(g), edt_(g.distance_transform()) {
  v.discs(n_discs, off_, r_);
  // Disc centre can sit anywhere in its cell, and the obstacle fills its whole cell, so pad
  // by one cell diagonal to make the check conservative.
  r_ += g.resolution() * std::sqrt(2.0);
}

bool FootprintChecker::free(const Pose& p) const {
  const double c = std::cos(p.th), s = std::sin(p.th);
  for (double o : off_) {
    const Cell cell = g_.to_cell(p.x + o * c, p.y + o * s);
    if (!g_.inside(cell.x, cell.y)) return false;
    if (edt_[static_cast<std::size_t>(g_.index(cell.x, cell.y))] < r_) return false;
  }
  return true;
}

namespace {

// Distance-to-goal over the grid, 8-connected, in metres. Computed once per query; a cheap
// lower bound that knows about walls (the Dubins term knows about turning radius).
std::vector<double> holonomic_distance(const Grid& g, const Cell& goal) {
  std::vector<double> d(static_cast<std::size_t>(g.size()), kInf);
  if (g.occupied(goal.x, goal.y)) return d;
  using E = std::pair<double, int>;
  std::priority_queue<E, std::vector<E>, std::greater<E>> pq;
  d[static_cast<std::size_t>(g.index(goal.x, goal.y))] = 0.0;
  pq.push({0.0, g.index(goal.x, goal.y)});
  const GridGraph8 G(g);
  while (!pq.empty()) {
    const auto [du, ui] = pq.top();
    pq.pop();
    if (du > d[static_cast<std::size_t>(ui)]) continue;
    G.for_each_neighbor(g.cell(ui), [&](const Cell& v, double c) {
      const int vi = g.index(v.x, v.y);
      const double nd = du + c * g.resolution();
      if (nd < d[static_cast<std::size_t>(vi)]) {
        d[static_cast<std::size_t>(vi)] = nd;
        pq.push({nd, vi});
      }
    });
  }
  return d;
}

struct Node {
  Pose p;
  double g;
  int parent;
  double steer;
  int dir;
};

// Integrate one motion primitive (kinematic bicycle, constant steering) in substeps of about
// `substep` metres, calling on_pose at each.
template <typename F>
bool roll(const Pose& from, double steer, int dir, double length, double wheelbase, double substep, F&& on_pose,
          Pose& end) {
  const int n = std::max(1, static_cast<int>(std::lround(length / substep)));
  const double ds = length / n;
  Pose q = from;
  for (int i = 0; i < n; ++i) {
    q.x += dir * ds * std::cos(q.th);
    q.y += dir * ds * std::sin(q.th);
    q.th = wrap_pi(q.th + dir * ds * std::tan(steer) / wheelbase);
    if (!on_pose(q)) return false;
  }
  end = q;
  return true;
}

}  // namespace

HybridResult hybrid_astar(const Grid& grid, const Vehicle& veh, const Pose& start, const Pose& goal,
                          const HybridConfig& cfg) {
  const auto t0 = std::chrono::steady_clock::now();
  HybridResult res;
  const FootprintChecker fc(grid, veh, 5);
  if (!fc.free(start) || !fc.free(goal)) return res;

  // Holonomic heuristic on a map inflated by the half width (the car centre line must stay
  // at least that far from walls).
  const Grid infl = grid.inflate(veh.width / 2.0);
  const Cell gc = grid.to_cell(goal.x, goal.y);
  const std::vector<double> holo = cfg.use_holonomic_heuristic ? holonomic_distance(infl, gc)
                                                               : std::vector<double>();
  const double rho = veh.wheelbase / std::tan(veh.max_steer);
  auto heuristic = [&](const Pose& p) {
    double h = std::hypot(goal.x - p.x, goal.y - p.y);
    if (cfg.use_holonomic_heuristic) {
      const Cell c = grid.to_cell(p.x, p.y);
      if (grid.inside(c.x, c.y)) h = std::max(h, holo[static_cast<std::size_t>(grid.index(c.x, c.y))]);
    }
    if (cfg.use_dubins_heuristic) h = std::max(h, dubins_shortest(p, goal, rho).length());
    return h;
  };

  const int nx = static_cast<int>(std::ceil(grid.width() * grid.resolution() / cfg.xy_res));
  const int ny = static_cast<int>(std::ceil(grid.height() * grid.resolution() / cfg.xy_res));
  const int nth = cfg.heading_bins;
  auto key = [&](const Pose& p) {
    const int ix = std::clamp(static_cast<int>(p.x / cfg.xy_res), 0, nx - 1);
    const int iy = std::clamp(static_cast<int>(p.y / cfg.xy_res), 0, ny - 1);
    int it = static_cast<int>(std::floor((wrap_pi(p.th) + M_PI) / (2 * M_PI) * nth));
    it = (it % nth + nth) % nth;
    return (static_cast<std::size_t>(iy) * nx + ix) * nth + it;
  };
  std::vector<double> best_g(static_cast<std::size_t>(nx) * ny * nth, kInf);
  std::vector<char> closed(best_g.size(), 0);
  std::vector<Node> nodes;
  nodes.reserve(1 << 16);
  struct Entry { double f; int id; };
  auto cmp = [](const Entry& a, const Entry& b) { return a.f > b.f; };
  std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)> open(cmp);

  nodes.push_back({start, 0.0, -1, 0.0, 1});
  best_g[key(start)] = 0.0;
  open.push({heuristic(start), 0});

  int goal_node = -1;
  DubinsPath finish;
  std::vector<double> steers;
  for (int i = 0; i < cfg.steer_samples; ++i)
    steers.push_back(-veh.max_steer + 2.0 * veh.max_steer * i / (cfg.steer_samples - 1));

  while (!open.empty() && res.expansions < cfg.max_expansions) {
    const Entry e = open.top();
    open.pop();
    const Node cur = nodes[static_cast<std::size_t>(e.id)];
    const std::size_t k = key(cur.p);
    if (closed[k]) continue;
    closed[k] = 1;
    ++res.expansions;

    if (std::hypot(cur.p.x - goal.x, cur.p.y - goal.y) < cfg.goal_xy_tol &&
        std::abs(wrap_pi(cur.p.th - goal.th)) < cfg.goal_th_tol) {
      goal_node = e.id;
      break;
    }
    // Analytic expansion: if a collision-free Dubins curve reaches the goal, stop here.
    const double dist = std::hypot(cur.p.x - goal.x, cur.p.y - goal.y);
    if (cfg.analytic_expansion && (res.expansions % cfg.analytic_every == 0 || dist < cfg.analytic_always_within)) {
      const DubinsPath dp = dubins_shortest(cur.p, goal, rho);
      bool ok = dp.valid;
      if (ok)
        for (const Pose& q : dp.sample_all(cfg.collision_check_step))
          if (!fc.free(q)) { ok = false; break; }
      if (ok) {
        goal_node = e.id;
        finish = dp;
        res.analytic_finish = true;
        break;
      }
    }
    for (int dir : {1, -1}) {
      if (dir < 0 && !cfg.allow_reverse) continue;
      for (double st : steers) {
        Pose end;
        if (!roll(cur.p, st, dir, cfg.step, veh.wheelbase, cfg.collision_check_step, [&](const Pose& q) { return fc.free(q); }, end))
          continue;
        const std::size_t nk = key(end);
        if (closed[nk]) continue;
        double c = cfg.step * (dir < 0 ? cfg.reverse_penalty : 1.0) + cfg.steer_penalty * std::abs(st) +
                   cfg.steer_change_penalty * std::abs(st - cur.steer);
        if (cur.parent >= 0 && dir != cur.dir) c += cfg.switch_penalty;
        const double ng = cur.g + c;
        if (ng < best_g[nk]) {
          best_g[nk] = ng;
          nodes.push_back({end, ng, e.id, st, dir});
          open.push({ng + heuristic(end), static_cast<int>(nodes.size() - 1)});
        }
      }
    }
  }
  res.time_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  if (goal_node < 0) return res;

  // Rebuild the dense path by re-rolling each primitive from its parent.
  std::vector<int> chain;
  for (int i = goal_node; i != -1; i = nodes[static_cast<std::size_t>(i)].parent) chain.push_back(i);
  std::reverse(chain.begin(), chain.end());
  res.path.push_back(start);
  res.dir.push_back(chain.size() > 1 ? nodes[static_cast<std::size_t>(chain[1])].dir : 1);
  for (std::size_t i = 1; i < chain.size(); ++i) {
    const Node& n = nodes[static_cast<std::size_t>(chain[i])];
    const Node& p = nodes[static_cast<std::size_t>(n.parent)];
    Pose end;
    roll(p.p, n.steer, n.dir, cfg.step, veh.wheelbase, cfg.collision_check_step,
         [&](const Pose& q) { res.path.push_back(q); res.dir.push_back(n.dir); return true; }, end);
  }
  if (res.analytic_finish) {
    const auto pts = finish.sample_all(cfg.collision_check_step);
    for (std::size_t i = 1; i < pts.size(); ++i) {
      res.path.push_back(pts[i]);
      res.dir.push_back(1);
    }
  }
  res.found = true;
  res.cost = nodes[static_cast<std::size_t>(goal_node)].g + (res.analytic_finish ? finish.length() : 0.0);
  for (std::size_t i = 1; i < res.path.size(); ++i) {
    res.length += std::hypot(res.path[i].x - res.path[i - 1].x, res.path[i].y - res.path[i - 1].y);
    if (res.dir[i] != res.dir[i - 1]) ++res.direction_switches;
  }
  return res;
}

}  // namespace pg
