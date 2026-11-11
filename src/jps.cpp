#include <algorithm>
#include <queue>

#include "pg/search.hpp"

namespace pg {

namespace {

inline int sgn(int v) { return (v > 0) - (v < 0); }

class Jps {
 public:
  Jps(const Grid& g, Cell goal) : g_(g), goal_(goal) {}

  // Walk from (x, y) in direction (dx, dy) until reaching a jump point, the goal, or a wall.
  // Iterative so long corridors do not blow the stack.
  bool jump(int x, int y, int dx, int dy, Cell& out) const {
    for (;;) {
      if (g_.occupied(x, y)) return false;
      if (x == goal_.x && y == goal_.y) { out = {x, y}; return true; }
      if (dx != 0 && dy != 0) {
        Cell tmp;
        // A diagonal step is a jump point if either straight scan from it finds one.
        if (jump(x + dx, y, dx, 0, tmp) || jump(x, y + dy, 0, dy, tmp)) { out = {x, y}; return true; }
      } else if (dx != 0) {
        if ((g_.free(x, y - 1) && g_.occupied(x - dx, y - 1)) ||
            (g_.free(x, y + 1) && g_.occupied(x - dx, y + 1))) { out = {x, y}; return true; }
      } else {
        if ((g_.free(x - 1, y) && g_.occupied(x - 1, y - dy)) ||
            (g_.free(x + 1, y) && g_.occupied(x + 1, y - dy))) { out = {x, y}; return true; }
      }
      // No corner cutting: a diagonal continues only if both orthogonal cells are free.
      if (!(g_.free(x + dx, y) && g_.free(x, y + dy))) return false;
      x += dx;
      y += dy;
    }
  }

  // Pruned neighbour set given the direction we arrived from.
  template <typename F>
  void neighbors(const Cell& c, const Cell* parent, F&& f) const {
    const int x = c.x, y = c.y;
    if (!parent) {
      for (int dy = -1; dy <= 1; ++dy)
        for (int dx = -1; dx <= 1; ++dx) {
          if (!dx && !dy) continue;
          if (g_.occupied(x + dx, y + dy)) continue;
          if (dx && dy && (g_.occupied(x + dx, y) || g_.occupied(x, y + dy))) continue;
          f(dx, dy);
        }
      return;
    }
    const int dx = sgn(x - parent->x), dy = sgn(y - parent->y);
    if (dx != 0 && dy != 0) {
      const bool v = g_.free(x, y + dy), h = g_.free(x + dx, y);
      if (v) f(0, dy);
      if (h) f(dx, 0);
      if (v && h) f(dx, dy);
    } else if (dx != 0) {
      const bool next = g_.free(x + dx, y), up = g_.free(x, y + 1), dn = g_.free(x, y - 1);
      if (next) {
        f(dx, 0);
        if (up) f(dx, 1);
        if (dn) f(dx, -1);
      }
      if (up) f(0, 1);
      if (dn) f(0, -1);
    } else {
      const bool next = g_.free(x, y + dy), rt = g_.free(x + 1, y), lt = g_.free(x - 1, y);
      if (next) {
        f(0, dy);
        if (rt) f(1, dy);
        if (lt) f(-1, dy);
      }
      if (rt) f(1, 0);
      if (lt) f(-1, 0);
    }
  }

 private:
  const Grid& g_;
  Cell goal_;
};

}  // namespace

SearchResult<Cell> jps(const Grid& grid, const Cell& start, const Cell& goal) {
  SearchResult<Cell> res;
  if (grid.occupied(start.x, start.y) || grid.occupied(goal.x, goal.y)) return res;
  const Jps J(grid, goal);
  const int n = grid.size();
  std::vector<double> g(static_cast<std::size_t>(n), kInf);
  std::vector<int> parent(static_cast<std::size_t>(n), -1);
  std::vector<char> closed(static_cast<std::size_t>(n), 0);
  struct Entry { double f, g; int idx; };
  auto cmp = [](const Entry& a, const Entry& b) { return a.f > b.f || (a.f == b.f && a.g < b.g); };
  std::priority_queue<Entry, std::vector<Entry>, decltype(cmp)> open(cmp);
  const int s = grid.index(start.x, start.y), t = grid.index(goal.x, goal.y);
  g[static_cast<std::size_t>(s)] = 0.0;
  open.push({octile(start, goal), 0.0, s});
  while (!open.empty()) {
    const Entry e = open.top();
    open.pop();
    if (closed[static_cast<std::size_t>(e.idx)]) continue;
    closed[static_cast<std::size_t>(e.idx)] = 1;
    ++res.expansions;
    if (e.idx == t) break;
    const Cell c = grid.cell(e.idx);
    const int pi = parent[static_cast<std::size_t>(e.idx)];
    const Cell pc = pi >= 0 ? grid.cell(pi) : Cell{};
    J.neighbors(c, pi >= 0 ? &pc : nullptr, [&](int dx, int dy) {
      Cell jp;
      if (!J.jump(c.x + dx, c.y + dy, dx, dy, jp)) return;
      const int ji = grid.index(jp.x, jp.y);
      if (closed[static_cast<std::size_t>(ji)]) return;
      const double ng = e.g + octile(c, jp);  // jump segments are straight or diagonal
      if (ng < g[static_cast<std::size_t>(ji)]) {
        g[static_cast<std::size_t>(ji)] = ng;
        parent[static_cast<std::size_t>(ji)] = e.idx;
        open.push({ng + octile(jp, goal), ng, ji});
      }
    });
  }
  if (g[static_cast<std::size_t>(t)] == kInf) return res;
  res.found = true;
  res.cost = g[static_cast<std::size_t>(t)];
  // Expand jump points back into a full cell path.
  std::vector<Cell> jps_pts;
  for (int i = t; i != -1; i = parent[static_cast<std::size_t>(i)]) jps_pts.push_back(grid.cell(i));
  std::reverse(jps_pts.begin(), jps_pts.end());
  res.path.push_back(jps_pts.front());
  for (std::size_t k = 1; k < jps_pts.size(); ++k) {
    Cell c = jps_pts[k - 1];
    const Cell d = jps_pts[k];
    while (c != d) {
      c.x += sgn(d.x - c.x);
      c.y += sgn(d.y - c.y);
      res.path.push_back(c);
    }
  }
  return res;
}

double path_cost(const Grid& grid, const std::vector<Cell>& path, bool* valid) {
  double c = 0.0;
  bool ok = !path.empty() && grid.free(path.front().x, path.front().y);
  for (std::size_t k = 1; k < path.size(); ++k) {
    const int dx = path[k].x - path[k - 1].x, dy = path[k].y - path[k - 1].y;
    if (std::abs(dx) > 1 || std::abs(dy) > 1 || (dx == 0 && dy == 0) || grid.occupied(path[k].x, path[k].y))
      ok = false;
    if (dx && dy && (grid.occupied(path[k - 1].x + dx, path[k - 1].y) ||
                     grid.occupied(path[k - 1].x, path[k - 1].y + dy)))
      ok = false;
    c += (dx && dy) ? kSqrt2 : 1.0;
  }
  if (valid) *valid = ok;
  return c;
}

}  // namespace pg
