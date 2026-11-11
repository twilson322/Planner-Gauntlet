#pragma once
#include <cstdint>
#include <utility>
#include <vector>

#include "pg/grid.hpp"

namespace pg {

struct Vec2 { double x{0}, y{0}; };

struct RrtConfig {
  int iterations{4000};
  double step{2.0};          // steer distance [m]
  double goal_bias{0.05};
  double goal_radius{1.5};
  double max_rewire_radius{6.0};
  bool informed{false};      // sample the prolate ellipse once a solution exists
  std::uint64_t seed{1};
};

struct RrtResult {
  bool found{false};
  double best_cost{1e18};
  std::vector<Vec2> path;
  std::vector<std::pair<int, double>> history;  // (iteration, best cost) whenever it improves
  int first_solution_iter{-1};
  int nodes{0};
  double time_ms{0};
};

// RRT* (Karaman & Frazzoli 2011) for a point robot on a configuration-space grid (inflate
// obstacles by the robot radius first). Rewiring propagates cost changes to the whole subtree.
RrtResult rrt_star(const Grid& cspace, const Vec2& start, const Vec2& goal, const RrtConfig& cfg);

// Collision check of a straight segment by sampling at half-cell spacing.
bool segment_free(const Grid& g, const Vec2& a, const Vec2& b);

}  // namespace pg
