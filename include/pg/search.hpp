#pragma once
#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <vector>

#include "pg/grid.hpp"

namespace pg {

inline constexpr double kInf = std::numeric_limits<double>::infinity();
inline const double kSqrt2 = std::sqrt(2.0);

template <typename Node>
struct SearchResult {
  bool found{false};
  double cost{kInf};
  std::vector<Node> path;
  long expansions{0};
};

// 8-connected grid graph. A diagonal move is only allowed when both orthogonal cells it
// passes are free, otherwise a robot of non-zero size would clip the obstacle corner.
// Edge costs are in cells (1 or sqrt 2).
class GridGraph8 {
 public:
  using Node = Cell;
  explicit GridGraph8(const Grid& g) : g_(&g) {}
  int size() const { return g_->size(); }
  int index(const Cell& c) const { return g_->index(c.x, c.y); }
  Cell node(int idx) const { return g_->cell(idx); }
  bool passable(const Cell& c) const { return g_->free(c.x, c.y); }

  template <typename F>
  void for_each_neighbor(const Cell& c, F&& f) const {
    static constexpr int DX[8] = {1, -1, 0, 0, 1, 1, -1, -1};
    static constexpr int DY[8] = {0, 0, 1, -1, 1, -1, 1, -1};
    for (int k = 0; k < 8; ++k) {
      const int nx = c.x + DX[k], ny = c.y + DY[k];
      if (g_->occupied(nx, ny)) continue;
      if (k >= 4 && (g_->occupied(c.x + DX[k], c.y) || g_->occupied(c.x, c.y + DY[k]))) continue;
      f(Cell{nx, ny}, k < 4 ? 1.0 : kSqrt2);
    }
  }

 private:
  const Grid* g_;
};

struct ZeroHeuristic {
  template <typename N>
  double operator()(const N&, const N&) const { return 0.0; }
};
struct OctileHeuristic {  // exact cost on an empty 8-connected grid, so admissible and consistent
  double operator()(const Cell& a, const Cell& b) const { return octile(a, b); }
};
struct EuclideanHeuristic {
  double operator()(const Cell& a, const Cell& b) const { return std::hypot(a.x - b.x, a.y - b.y); }
};

// Best-first search over any graph type providing size(), index(), node(), passable() and
// for_each_neighbor(). weight = 1 gives A* (Dijkstra with ZeroHeuristic); weight > 1 gives
// weighted A*, whose cost is bounded by weight * optimal for a consistent heuristic.
// Binary heap with lazy deletion; ties on f broken toward larger g (deeper nodes), which
// cuts expansions on open grids without affecting optimality.
template <typename Graph, typename Heuristic>
SearchResult<typename Graph::Node> astar(const Graph& G, const typename Graph::Node& start,
                                         const typename Graph::Node& goal, Heuristic h,
                                         double weight = 1.0) {
  using Node = typename Graph::Node;
  SearchResult<Node> res;
  if (!G.passable(start) || !G.passable(goal)) return res;
  const int n = G.size();
  std::vector<double> g(static_cast<std::size_t>(n), kInf);
  std::vector<int> parent(static_cast<std::size_t>(n), -1);
  std::vector<char> closed(static_cast<std::size_t>(n), 0);
  struct Entry { double f, g; int idx; };
  auto cmp = [](const Entry& a, const Entry& b) { return a.f > b.f || (a.f == b.f && a.g < b.g); };
  std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)> open(cmp);

  const int s = G.index(start), t = G.index(goal);
  g[static_cast<std::size_t>(s)] = 0.0;
  open.push({weight * h(start, goal), 0.0, s});
  while (!open.empty()) {
    const Entry e = open.top();
    open.pop();
    if (closed[static_cast<std::size_t>(e.idx)]) continue;  // stale duplicate
    closed[static_cast<std::size_t>(e.idx)] = 1;
    ++res.expansions;
    if (e.idx == t) break;
    const Node u = G.node(e.idx);
    G.for_each_neighbor(u, [&](const Node& v, double c) {
      const int vi = G.index(v);
      if (closed[static_cast<std::size_t>(vi)]) return;
      const double ng = g[static_cast<std::size_t>(e.idx)] + c;
      if (ng < g[static_cast<std::size_t>(vi)]) {
        g[static_cast<std::size_t>(vi)] = ng;
        parent[static_cast<std::size_t>(vi)] = e.idx;
        open.push({ng + weight * h(v, goal), ng, vi});
      }
    });
  }
  if (g[static_cast<std::size_t>(t)] == kInf) return res;
  res.found = true;
  res.cost = g[static_cast<std::size_t>(t)];
  for (int i = t; i != -1; i = parent[static_cast<std::size_t>(i)]) res.path.push_back(G.node(i));
  std::reverse(res.path.begin(), res.path.end());
  return res;
}

// Jump Point Search (Harabor & Grastien 2011), in the variant that forbids corner cutting so
// it searches exactly the same graph as GridGraph8. Returns the full cell path.
SearchResult<Cell> jps(const Grid& grid, const Cell& start, const Cell& goal);

// Sum of edge costs along a cell path (validates that consecutive cells are legal moves).
double path_cost(const Grid& grid, const std::vector<Cell>& path, bool* valid = nullptr);

}  // namespace pg
