#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "pg/dubins.hpp"
#include "pg/hybrid_astar.hpp"

using namespace pg;

TEST(Dubins, EndsAtGoalAndRespectsCurvature) {
  std::mt19937_64 r(5);
  std::uniform_real_distribution<double> P(-20, 20), A(-M_PI, M_PI);
  for (int i = 0; i < 2000; ++i) {
    const Pose a{P(r), P(r), A(r)}, b{P(r), P(r), A(r)};
    const double rho = 4.0;
    const DubinsPath d = dubins_shortest(a, b, rho);
    ASSERT_TRUE(d.valid);
    const Pose e = d.sample(d.length());
    EXPECT_NEAR(e.x, b.x, 1e-6);
    EXPECT_NEAR(e.y, b.y, 1e-6);
    EXPECT_NEAR(wrap_pi(e.th - b.th), 0.0, 1e-6);
    EXPECT_GE(d.length(), std::hypot(b.x - a.x, b.y - a.y) - 1e-9);
    // Heading changes by at most ds / rho per step of length ds.
    const auto pts = d.sample_all(0.05);
    for (std::size_t k = 1; k < pts.size(); ++k) {
      const double ds = std::hypot(pts[k].x - pts[k - 1].x, pts[k].y - pts[k - 1].y);
      EXPECT_LE(std::abs(wrap_pi(pts[k].th - pts[k - 1].th)), ds / rho + 1e-6);
    }
  }
}

TEST(Dubins, StraightLineIsExact) {
  const DubinsPath d = dubins_shortest({0, 0, 0}, {10, 0, 0}, 3.0);
  EXPECT_NEAR(d.length(), 10.0, 1e-9);
}

TEST(HybridAStar, ParksWithoutCollision) {
  const Grid lot = Grid::parking_lot(0.1);
  const Vehicle v;
  const Pose start{6.0, 30.0, 0.0}, goal{39.25, 4.85, -M_PI / 2};
  const HybridResult r = hybrid_astar(lot, v, start, goal, HybridConfig{});
  ASSERT_TRUE(r.found);
  const FootprintChecker fc(lot, v, 5);
  for (const Pose& p : r.path) ASSERT_TRUE(fc.free(p));
  EXPECT_NEAR(r.path.back().x, goal.x, 0.3);
  EXPECT_NEAR(r.path.back().y, goal.y, 0.3);
  EXPECT_NEAR(wrap_pi(r.path.back().th - goal.th), 0.0, 0.09);
  // Consecutive samples never jump (continuity of the stitched primitives + Dubins shot).
  for (std::size_t k = 1; k < r.path.size(); ++k)
    EXPECT_LT(std::hypot(r.path[k].x - r.path[k - 1].x, r.path[k].y - r.path[k - 1].y), 0.11);
}
