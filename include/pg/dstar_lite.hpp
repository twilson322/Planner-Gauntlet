#pragma once
#include <cmath>
#include <queue>
#include <utility>
#include <vector>

#include "pg/grid.hpp"
#include "pg/search.hpp"

namespace pg {

// D* Lite (Koenig & Likhachev 2002) on the same 8-connected, no-corner-cutting grid graph.
// It searches backward from the goal, so when the robot moves only the key modifier km
// changes, and when the map changes only the vertices whose edges changed are repaired.
// The robot plans on its *belief* map, which starts optimistic (unknown = free).
class DStarLite {
 public:
  DStarLite(Grid belief, Cell start, Cell goal);

  void compute_shortest_path();
  // Move the robot. The key modifier km is only advanced when the map changes.
  void move_to(const Cell& s);
  // Apply sensed occupancy changes and repair affected vertices. Returns true if any changed.
  bool update_cells(const std::vector<std::pair<Cell, bool>>& changes);

  double cost_to_goal() const { return g_[idx(start_)]; }
  std::vector<Cell> extract_path(int max_steps = 1 << 20) const;
  Cell next_step() const;  // greedy successor of the current start
  long expansions() const { return expansions_; }
  const Grid& belief() const { return belief_; }
  const Cell& start() const { return start_; }
  double g_at(const Cell& c) const { return g_[idx(c)]; }
  double rhs_at(const Cell& c) const { return rhs_[idx(c)]; }
  double km() const { return km_; }

 private:
  struct Key {
    double k1, k2;
    // Tolerant comparison. On a grid with sqrt(2) edges, g + h for vertices on the same
    // optimal path agree only up to rounding; an exact compare decided a tie by 1e-15 and
    // stopped the search before the vertex that mattered was expanded.
    bool operator<(const Key& o) const {
      return k1 < o.k1 - 1e-9 || (std::abs(k1 - o.k1) <= 1e-9 && k2 < o.k2 - 1e-9);
    }
  };
  struct Entry {
    Key key;
    int idx;
    bool operator>(const Entry& o) const { return o.key < key; }
  };

  std::size_t idx(const Cell& c) const { return static_cast<std::size_t>(belief_.index(c.x, c.y)); }
  double h(const Cell& a, const Cell& b) const { return octile(a, b); }
  double cost(const Cell& a, const Cell& b) const;
  Key calc_key(const Cell& s) const;
  void update_vertex(const Cell& u);
  template <typename F>
  void for_each_adjacent(const Cell& c, F&& f) const {  // all 8 in-bounds neighbours
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx)
        if ((dx || dy) && belief_.inside(c.x + dx, c.y + dy)) f(Cell{c.x + dx, c.y + dy});
  }
  bool top_key(Key& k);

  Grid belief_;
  Cell start_, goal_, last_;
  double km_{0.0};
  std::vector<double> g_, rhs_;
  std::vector<Key> open_key_;
  std::vector<char> in_open_;
  std::priority_queue<Entry, std::vector<Entry>, std::greater<Entry>> open_;
  long expansions_{0};
};

}  // namespace pg
