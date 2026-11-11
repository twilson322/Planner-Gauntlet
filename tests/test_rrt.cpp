#include <gtest/gtest.h>

#include "pg/rrt_star.hpp"

using namespace pg;

TEST(RrtStar, FindsCollisionFreePathAndCostOnlyImproves) {
  Grid g(200, 200, 0.25);
  g.fill_rect(80, 0, 100, 150);  // a wall with a gap at the top
  for (int informed = 0; informed <= 1; ++informed) {
    RrtConfig c;
    c.iterations = 3000;
    c.informed = informed;
    c.seed = 3;
    const auto r = rrt_star(g.inflate(0.5), {5, 5}, {45, 5}, c);
    ASSERT_TRUE(r.found);
    for (std::size_t k = 1; k < r.path.size(); ++k) EXPECT_TRUE(segment_free(g.inflate(0.5), r.path[k - 1], r.path[k]));
    for (std::size_t k = 1; k < r.history.size(); ++k) EXPECT_LT(r.history[k].second, r.history[k - 1].second);
    // Must go over the wall: the path cannot be shorter than the two straight legs via the gap.
    EXPECT_GT(r.best_cost, 2 * std::hypot(20.0, 33.0) - 1.0);
  }
}
