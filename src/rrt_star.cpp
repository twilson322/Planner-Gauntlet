#include "pg/rrt_star.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <random>

namespace pg {

bool segment_free(const Grid& g, const Vec2& a, const Vec2& b) {
  const double L = std::hypot(b.x - a.x, b.y - a.y);
  const int n = std::max(1, static_cast<int>(std::ceil(L / (0.5 * g.resolution()))));
  for (int i = 0; i <= n; ++i) {
    const double t = static_cast<double>(i) / n;
    if (g.occupied_world(a.x + t * (b.x - a.x), a.y + t * (b.y - a.y))) return false;
  }
  return true;
}

namespace {
struct Node {
  Vec2 p;
  int parent;
  double cost;
  std::vector<int> children;
};
double dist(const Vec2& a, const Vec2& b) { return std::hypot(a.x - b.x, a.y - b.y); }
}  // namespace

RrtResult rrt_star(const Grid& cs, const Vec2& start, const Vec2& goal, const RrtConfig& cfg) {
  const auto t0 = std::chrono::steady_clock::now();
  RrtResult res;
  const double W = cs.width() * cs.resolution(), H = cs.height() * cs.resolution();
  std::mt19937_64 rng(cfg.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<Node> T;
  T.reserve(static_cast<std::size_t>(cfg.iterations) + 2);
  T.push_back({start, -1, 0.0, {}});
  std::vector<int> goal_parents;  // nodes that can connect straight to the goal

  // gamma_RRT* > 2 (1 + 1/d)^(1/d) (mu(X_free) / zeta_d)^(1/d), d = 2.
  const double gamma = 2.0 * std::sqrt(1.5) * std::sqrt(W * H / M_PI);
  const double c_min = dist(start, goal);
  const double cth = (goal.x - start.x) / c_min, sth = (goal.y - start.y) / c_min;
  const Vec2 centre{0.5 * (start.x + goal.x), 0.5 * (start.y + goal.y)};

  auto sample = [&]() -> Vec2 {
    if (U(rng) < cfg.goal_bias) return goal;
    if (cfg.informed && res.best_cost < 1e17) {
      // Uniform sample in the prolate ellipse with foci start/goal and major axis c_best.
      // Every point that could shorten the current solution lies inside it.
      const double a = 0.5 * res.best_cost, b = 0.5 * std::sqrt(std::max(0.0, res.best_cost * res.best_cost - c_min * c_min));
      for (;;) {
        const double r = std::sqrt(U(rng)), th = 2 * M_PI * U(rng);
        const double ux = a * r * std::cos(th), uy = b * r * std::sin(th);
        const Vec2 p{centre.x + cth * ux - sth * uy, centre.y + sth * ux + cth * uy};
        if (p.x >= 0 && p.y >= 0 && p.x < W && p.y < H) return p;
      }
    }
    return {U(rng) * W, U(rng) * H};
  };

  // Lower every cost in a subtree by delta after its root was rewired.
  auto propagate = [&](int root, double delta) {
    std::vector<int> stack{root};
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      for (int c : T[static_cast<std::size_t>(u)].children) {
        T[static_cast<std::size_t>(c)].cost -= delta;
        stack.push_back(c);
      }
    }
  };

  std::vector<int> near;
  for (int it = 0; it < cfg.iterations; ++it) {
    const Vec2 xr = sample();
    int ni = 0;
    double nd = dist(T[0].p, xr);
    for (std::size_t i = 1; i < T.size(); ++i) {
      const double d = dist(T[i].p, xr);
      if (d < nd) { nd = d; ni = static_cast<int>(i); }
    }
    const Vec2& xn = T[static_cast<std::size_t>(ni)].p;
    const Vec2 xnew = nd <= cfg.step ? xr : Vec2{xn.x + (xr.x - xn.x) * cfg.step / nd, xn.y + (xr.y - xn.y) * cfg.step / nd};
    if (cs.occupied_world(xnew.x, xnew.y) || !segment_free(cs, xn, xnew)) continue;

    const double n = static_cast<double>(T.size() + 1);
    const double r = std::min(gamma * std::sqrt(std::log(n) / n), cfg.max_rewire_radius);
    near.clear();
    for (std::size_t i = 0; i < T.size(); ++i)
      if (dist(T[i].p, xnew) <= r) near.push_back(static_cast<int>(i));

    // Choose the parent that gives the cheapest collision-free connection.
    int best_parent = ni;
    double best_cost = T[static_cast<std::size_t>(ni)].cost + dist(xn, xnew);
    for (int j : near) {
      const double c = T[static_cast<std::size_t>(j)].cost + dist(T[static_cast<std::size_t>(j)].p, xnew);
      if (c < best_cost && segment_free(cs, T[static_cast<std::size_t>(j)].p, xnew)) {
        best_cost = c;
        best_parent = j;
      }
    }
    const int id = static_cast<int>(T.size());
    T.push_back({xnew, best_parent, best_cost, {}});
    T[static_cast<std::size_t>(best_parent)].children.push_back(id);

    // Rewire: route neighbours through the new node when that is cheaper.
    for (int j : near) {
      if (j == best_parent) continue;
      Node& nj = T[static_cast<std::size_t>(j)];
      const double c = best_cost + dist(xnew, nj.p);
      if (c + 1e-9 < nj.cost && segment_free(cs, xnew, nj.p)) {
        auto& sib = T[static_cast<std::size_t>(nj.parent)].children;
        sib.erase(std::remove(sib.begin(), sib.end(), j), sib.end());
        const double delta = nj.cost - c;
        nj.parent = id;
        nj.cost = c;
        T[static_cast<std::size_t>(id)].children.push_back(j);
        propagate(j, delta);
      }
    }
    if (dist(xnew, goal) <= cfg.goal_radius && segment_free(cs, xnew, goal)) goal_parents.push_back(id);

    double best = res.best_cost;
    int best_gp = -1;
    for (int gp : goal_parents) {
      const double c = T[static_cast<std::size_t>(gp)].cost + dist(T[static_cast<std::size_t>(gp)].p, goal);
      if (c < best - 1e-9) { best = c; best_gp = gp; }
    }
    if (best_gp >= 0) {
      if (res.first_solution_iter < 0) res.first_solution_iter = it;
      res.best_cost = best;
      res.history.push_back({it, best});
      res.path.clear();
      res.path.push_back(goal);
      for (int u = best_gp; u != -1; u = T[static_cast<std::size_t>(u)].parent) res.path.push_back(T[static_cast<std::size_t>(u)].p);
      std::reverse(res.path.begin(), res.path.end());
    }
  }
  res.found = res.best_cost < 1e17;
  res.nodes = static_cast<int>(T.size());
  res.time_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
  return res;
}

}  // namespace pg
