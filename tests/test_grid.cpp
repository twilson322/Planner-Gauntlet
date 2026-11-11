#include <gtest/gtest.h>

#include <cmath>

#include "pg/grid.hpp"

using namespace pg;

TEST(Grid, DistanceTransformMatchesBruteForce) {
  for (int seed = 0; seed < 20; ++seed) {
    const Grid g = Grid::random_rects(23, 17, 0.5, 4, 1, 4, static_cast<std::uint64_t>(seed));
    const auto dt = g.distance_transform();
    for (int y = 0; y < g.height(); ++y)
      for (int x = 0; x < g.width(); ++x) {
        double best = 1e30;
        for (int v = 0; v < g.height(); ++v)
          for (int u = 0; u < g.width(); ++u)
            if (g.occupied(u, v)) best = std::min(best, std::hypot(u - x, v - y) * 0.5);
        if (best > 1e29) continue;  // no obstacle at all
        ASSERT_NEAR(dt[static_cast<std::size_t>(g.index(x, y))], best, 1e-5) << x << "," << y;
      }
  }
}

TEST(Grid, OutOfBoundsIsOccupied) {
  Grid g(4, 4, 1.0);
  EXPECT_TRUE(g.occupied(-1, 0));
  EXPECT_TRUE(g.occupied(0, 4));
  EXPECT_FALSE(g.occupied(3, 3));
}

TEST(Grid, InflateIsMonotone) {
  const Grid g = Grid::random_rects(40, 40, 0.25, 6, 1, 6, 3);
  const Grid a = g.inflate(0.5), b = g.inflate(1.0);
  for (int y = 0; y < 40; ++y)
    for (int x = 0; x < 40; ++x) {
      if (g.occupied(x, y)) {
        EXPECT_TRUE(a.occupied(x, y));
      }
      if (a.occupied(x, y)) {
        EXPECT_TRUE(b.occupied(x, y));
      }
    }
}
