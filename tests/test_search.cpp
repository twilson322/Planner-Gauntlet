#include <gtest/gtest.h>

#include <random>

#include "pg/search.hpp"

using namespace pg;

namespace {
Cell rand_free(const Grid& g, std::mt19937_64& r) {
  std::uniform_int_distribution<int> X(0, g.width() - 1), Y(0, g.height() - 1);
  for (;;) {
    Cell c{X(r), Y(r)};
    if (g.free(c.x, c.y)) return c;
  }
}
}  // namespace

// Property test: on 300 random maps every optimal algorithm returns the same cost, every
// returned path is a legal sequence of moves whose length equals the reported cost, and
// weighted A* stays inside its suboptimality bound.
TEST(Search, OptimalAlgorithmsAgreeOnRandomMaps) {
  std::mt19937_64 r(11);
  int solved = 0;
  for (int m = 0; m < 300; ++m) {
    const Grid g = Grid::random_rects(48, 48, 1.0, 20, 1, 8, static_cast<std::uint64_t>(m));
    const Cell s = rand_free(g, r), t = rand_free(g, r);
    const GridGraph8 G(g);
    const auto dj = astar(G, s, t, ZeroHeuristic{});
    const auto as = astar(G, s, t, OctileHeuristic{});
    const auto jp = jps(g, s, t);
    ASSERT_EQ(dj.found, as.found);
    ASSERT_EQ(dj.found, jp.found);
    if (!dj.found) continue;
    ++solved;
    EXPECT_NEAR(as.cost, dj.cost, 1e-9);
    EXPECT_NEAR(jp.cost, dj.cost, 1e-9);
    for (const auto* res : {&dj, &as, &jp}) {
      bool valid = false;
      EXPECT_NEAR(path_cost(g, res->path, &valid), res->cost, 1e-9);
      EXPECT_TRUE(valid);
      EXPECT_EQ(res->path.front(), s);
      EXPECT_EQ(res->path.back(), t);
    }
    EXPECT_LE(as.expansions, dj.expansions);
    for (double w : {1.5, 2.0, 5.0}) {
      const auto wa = astar(G, s, t, OctileHeuristic{}, w);
      EXPECT_LE(wa.cost, w * dj.cost + 1e-9);
    }
  }
  EXPECT_GT(solved, 200);
}

TEST(Search, NoCornerCutting) {
  Grid g(3, 3, 1.0);
  g.set(1, 0, true);
  g.set(0, 1, true);  // (0,0) -> (1,1) diagonal passes between two obstacles
  const auto r = astar(GridGraph8(g), Cell{0, 0}, Cell{2, 2}, OctileHeuristic{});
  EXPECT_FALSE(r.found);
  EXPECT_FALSE(jps(g, Cell{0, 0}, Cell{2, 2}).found);
}

TEST(Search, BlockedEndpoints) {
  Grid g(5, 5, 1.0);
  g.set(4, 4, true);
  EXPECT_FALSE(astar(GridGraph8(g), Cell{0, 0}, Cell{4, 4}, OctileHeuristic{}).found);
  EXPECT_FALSE(jps(g, Cell{0, 0}, Cell{4, 4}).found);
}
