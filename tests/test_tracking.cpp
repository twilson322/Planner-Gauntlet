#include <gtest/gtest.h>

#include <cmath>

#include "pg/tracking.hpp"

using namespace pg;

TEST(Tracking, RoadIsContinuousAndSpeedRespectsLateralLimit) {
  const RefPath p = RefPath::road(20.0, 2.0, 1.5);
  for (std::size_t i = 1; i < p.s.size(); ++i) {
    EXPECT_NEAR(std::hypot(p.x[i] - p.x[i - 1], p.y[i] - p.y[i - 1]), 0.5, 1e-3);
    EXPECT_LE(p.v[i] * p.v[i] * std::abs(p.kappa[i]), 2.0 + 1e-9);
  }
}

TEST(Tracking, LqrWithTrueStateStaysWithinFiveCentimetres) {
  ScenarioConfig c;
  c.tracker = "lqr";
  c.use_ekf = false;
  const auto r = run_scenario(RefPath::road(20.0, 2.0, 1.5), c);
  EXPECT_TRUE(r.completed);
  EXPECT_LT(r.rms_lat, 0.05);
}

TEST(Tracking, GatedEkfRejectsOutliersAndStaysAccurate) {
  ScenarioConfig c;
  c.tracker = "pure_pursuit";
  c.gps_outlier_prob = 0.05;
  const auto r = run_scenario(RefPath::road(20.0, 2.0, 1.5), c);
  EXPECT_TRUE(r.completed);
  EXPECT_GE(r.gps_rejected, r.gps_outliers * 9 / 10);  // nearly every 8 m jump gated out
  EXPECT_LT(r.rms_pos_err, 0.3);
}

TEST(RefPath, IndexAheadUsesArcLengthNotSampleCount) {
  pg::RefPath p;
  p.s = {0.0, 0.2, 1.0, 1.1, 3.0, 4.0};  // deliberately uneven spacing
  p.x = p.y = p.th = p.kappa = p.v = std::vector<double>(p.s.size(), 0.0);
  EXPECT_EQ(p.index_ahead(0, 1.05), 2);  // last sample within 1.05 m
  EXPECT_EQ(p.index_ahead(1, 0.85), 2);
  EXPECT_EQ(p.index_ahead(3, 0.5), 3);   // nothing new within reach: stay
  EXPECT_EQ(p.index_ahead(2, 100.0), 5);  // clamped to the end
}
