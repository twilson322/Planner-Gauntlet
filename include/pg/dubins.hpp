#pragma once
#include <array>
#include <vector>

namespace pg {

struct Pose {
  double x{0}, y{0}, th{0};
};

inline double wrap_pi(double a) {
  while (a > 3.14159265358979323846) a -= 2 * 3.14159265358979323846;
  while (a < -3.14159265358979323846) a += 2 * 3.14159265358979323846;
  return a;
}

// Shortest forward-only path for a car with minimum turning radius rho (Dubins 1957).
// Six candidate words (LSL, RSR, LSR, RSL, RLR, LRL); the shortest valid one wins.
struct DubinsPath {
  Pose q0;
  double rho{1.0};
  std::array<double, 3> seg{};  // segment lengths normalised by rho
  std::array<char, 3> type{};   // 'L', 'S', 'R'
  bool valid{false};
  double length() const { return (seg[0] + seg[1] + seg[2]) * rho; }
  Pose sample(double s) const;                         // pose at arc length s
  std::vector<Pose> sample_all(double step) const;     // includes the end pose
};

DubinsPath dubins_shortest(const Pose& q0, const Pose& q1, double rho);

}  // namespace pg
