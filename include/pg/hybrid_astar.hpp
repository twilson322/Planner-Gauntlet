#pragma once
#include <vector>

#include "pg/dubins.hpp"
#include "pg/grid.hpp"

namespace pg {

// Car footprint, referenced to the rear-axle centre.
struct Vehicle {
  double wheelbase{2.8};
  double length{4.6}, width{1.9};
  double rear_overhang{0.9};
  double max_steer{0.6};  // rad
  // Footprint approximated by n discs along the centre line; returns centre offsets + radius.
  void discs(int n, std::vector<double>& offsets, double& radius) const;
};

struct HybridConfig {
  double xy_res{0.5};          // closed-set cell size [m]
  int heading_bins{72};        // 5 degrees
  double step{1.0};            // arc length per expansion [m]
  int steer_samples{5};        // -max .. +max
  bool allow_reverse{true};
  double reverse_penalty{2.0}, switch_penalty{6.0}, steer_penalty{0.3}, steer_change_penalty{0.4};
  bool analytic_expansion{true};
  int analytic_every{8};       // try a Dubins shot every N expansions...
  double analytic_always_within{10.0};  // ...and at every expansion within this distance of the goal [m]
  double collision_check_step{0.1};     // spacing of footprint checks along primitives and shots [m]
  double goal_xy_tol{0.3}, goal_th_tol{0.087};
  long max_expansions{200000};
  bool use_holonomic_heuristic{true};
  bool use_dubins_heuristic{true};  // only admissible-ish when reverse is disallowed
};

struct HybridResult {
  bool found{false};
  bool analytic_finish{false};
  std::vector<Pose> path;   // densely sampled (0.1 m)
  std::vector<int> dir;     // +1 forward, -1 reverse, per path point
  long expansions{0};
  double cost{0};
  double length{0};
  int direction_switches{0};
  double time_ms{0};
};

// Collision checker using the obstacle distance field: a pose is free if every footprint disc
// clears the nearest obstacle cell centre by its radius plus half a cell diagonal.
class FootprintChecker {
 public:
  FootprintChecker(const Grid& g, const Vehicle& v, int n_discs = 3);
  bool free(const Pose& p) const;

 private:
  const Grid& g_;
  std::vector<float> edt_;
  std::vector<double> off_;
  double r_;
};

HybridResult hybrid_astar(const Grid& grid, const Vehicle& veh, const Pose& start, const Pose& goal,
                          const HybridConfig& cfg);

}  // namespace pg
