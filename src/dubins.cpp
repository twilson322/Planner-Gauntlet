#include "pg/dubins.hpp"

#include <cmath>
#include <limits>

namespace pg {

namespace {
constexpr double kTwoPi = 2.0 * 3.14159265358979323846;
double mod2pi(double a) {
  a = std::fmod(a, kTwoPi);
  return a < 0 ? a + kTwoPi : a;
}

// One segment of unit radius starting at pose q, length t (radians or unit lengths).
Pose step(const Pose& q, char type, double t) {
  Pose r;
  if (type == 'L') {
    r.x = q.x + std::sin(q.th + t) - std::sin(q.th);
    r.y = q.y - std::cos(q.th + t) + std::cos(q.th);
    r.th = q.th + t;
  } else if (type == 'R') {
    r.x = q.x - std::sin(q.th - t) + std::sin(q.th);
    r.y = q.y + std::cos(q.th - t) - std::cos(q.th);
    r.th = q.th - t;
  } else {
    r.x = q.x + std::cos(q.th) * t;
    r.y = q.y + std::sin(q.th) * t;
    r.th = q.th;
  }
  return r;
}
}  // namespace

DubinsPath dubins_shortest(const Pose& q0, const Pose& q1, double rho) {
  const double dx = q1.x - q0.x, dy = q1.y - q0.y;
  const double d = std::hypot(dx, dy) / rho;
  const double theta = mod2pi(std::atan2(dy, dx));
  const double a = mod2pi(q0.th - theta), b = mod2pi(q1.th - theta);
  const double sa = std::sin(a), sb = std::sin(b), ca = std::cos(a), cb = std::cos(b), cab = std::cos(a - b);

  DubinsPath best;
  best.q0 = q0;
  best.rho = rho;
  double best_len = std::numeric_limits<double>::infinity();
  auto consider = [&](double t, double p, double q, const char* w) {
    const double L = t + p + q;
    if (L < best_len) {
      best_len = L;
      best.seg = {t, p, q};
      best.type = {w[0], w[1], w[2]};
      best.valid = true;
    }
  };
  {  // LSL
    const double p2 = 2 + d * d - 2 * cab + 2 * d * (sa - sb);
    if (p2 >= 0) {
      const double t1 = std::atan2(cb - ca, d + sa - sb);
      consider(mod2pi(-a + t1), std::sqrt(p2), mod2pi(b - t1), "LSL");
    }
  }
  {  // RSR
    const double p2 = 2 + d * d - 2 * cab + 2 * d * (sb - sa);
    if (p2 >= 0) {
      const double t1 = std::atan2(ca - cb, d - sa + sb);
      consider(mod2pi(a - t1), std::sqrt(p2), mod2pi(-b + t1), "RSR");
    }
  }
  {  // LSR
    const double p2 = -2 + d * d + 2 * cab + 2 * d * (sa + sb);
    if (p2 >= 0) {
      const double p = std::sqrt(p2);
      const double t2 = std::atan2(-ca - cb, d + sa + sb) - std::atan2(-2.0, p);
      consider(mod2pi(-a + t2), p, mod2pi(-mod2pi(b) + t2), "LSR");
    }
  }
  {  // RSL
    const double p2 = d * d - 2 + 2 * cab - 2 * d * (sa + sb);
    if (p2 >= 0) {
      const double p = std::sqrt(p2);
      const double t2 = std::atan2(ca + cb, d - sa - sb) - std::atan2(2.0, p);
      consider(mod2pi(a - t2), p, mod2pi(b - t2), "RSL");
    }
  }
  {  // RLR
    const double tmp = (6.0 - d * d + 2 * cab + 2 * d * (sa - sb)) / 8.0;
    if (std::abs(tmp) <= 1.0) {
      const double p = mod2pi(kTwoPi - std::acos(tmp));
      const double t = mod2pi(a - std::atan2(ca - cb, d - sa + sb) + p / 2.0);
      consider(t, p, mod2pi(a - b - t + p), "RLR");
    }
  }
  {  // LRL
    const double tmp = (6.0 - d * d + 2 * cab + 2 * d * (-sa + sb)) / 8.0;
    if (std::abs(tmp) <= 1.0) {
      const double p = mod2pi(kTwoPi - std::acos(tmp));
      const double t = mod2pi(-a - std::atan2(ca - cb, d + sa - sb) + p / 2.0);
      consider(t, p, mod2pi(mod2pi(b) - a - t + p), "LRL");
    }
  }
  return best;
}

Pose DubinsPath::sample(double s) const {
  // Work in a frame scaled by 1/rho with origin at q0.
  double t = s / rho;
  Pose q{0.0, 0.0, q0.th};
  for (int i = 0; i < 3 && t > 0.0; ++i) {
    const double l = std::min(t, seg[static_cast<std::size_t>(i)]);
    q = step(q, type[static_cast<std::size_t>(i)], l);
    t -= l;
  }
  return {q0.x + q.x * rho, q0.y + q.y * rho, wrap_pi(q.th)};
}

std::vector<Pose> DubinsPath::sample_all(double ds) const {
  std::vector<Pose> out;
  const double L = length();
  for (double s = 0.0; s < L; s += ds) out.push_back(sample(s));
  out.push_back(sample(L));
  return out;
}

}  // namespace pg
