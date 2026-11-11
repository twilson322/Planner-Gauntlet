#include <gtest/gtest.h>

#include "pg/dstar_lite.hpp"

using namespace pg;

// Property test: a robot that senses a radius-3 disc each step and replans with D* Lite must
// always hold the same cost-to-goal as A* run from scratch on the same belief, and must reach
// the goal whenever the true map allows it.
TEST(DStarLite, IncrementalEqualsFromScratch) {
  int runs = 0;
  for (int trial = 0; trial < 150; ++trial) {
    const int S = 30, R = 3;
    Grid truth = Grid::random_rects(S, S, 1.0, 12, 1, 6, static_cast<std::uint64_t>(trial));
    const Cell s{0, 0}, g{S - 1, S - 1};
    truth.set(s.x, s.y, false);
    truth.set(g.x, g.y, false);
    if (!astar(GridGraph8(truth), s, g, OctileHeuristic{}).found) continue;
    ++runs;
    DStarLite ds(Grid(S, S, 1.0), s, g);
    ds.compute_shortest_path();
    Cell pos = s;
    for (int step = 0; step < 4 * S * S && pos != g; ++step) {
      std::vector<std::pair<Cell, bool>> ch;
      for (int dy = -R; dy <= R; ++dy)
        for (int dx = -R; dx <= R; ++dx) {
          const Cell c{pos.x + dx, pos.y + dy};
          if (truth.inside(c.x, c.y) && truth.occupied(c.x, c.y) != ds.belief().occupied(c.x, c.y))
            ch.push_back({c, truth.occupied(c.x, c.y)});
        }
      if (!ch.empty()) {
        ds.update_cells(ch);
        ds.compute_shortest_path();
        const auto a = astar(GridGraph8(ds.belief()), pos, g, OctileHeuristic{});
        ASSERT_NEAR(ds.cost_to_goal(), a.cost, 1e-9) << "trial " << trial << " step " << step;
      }
      pos = ds.next_step();
      ds.move_to(pos);
    }
    EXPECT_EQ(pos, g) << "trial " << trial;
  }
  EXPECT_GT(runs, 100);
}
