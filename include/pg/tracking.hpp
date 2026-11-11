#pragma once
#include <Eigen/Dense>
#include <cstdint>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "pg/dubins.hpp"

namespace pg {

// ------------------------------------------------------------------ reference path
// Dense path (0.5 m spacing) built by integrating a curvature profile, with a speed profile
// limited by lateral acceleration and longitudinal accel/decel (forward-backward pass).
struct RefPath {
  std::vector<double> s, x, y, th, kappa, v;
  static RefPath road(double v_max, double a_lat_max, double a_lon_max);
  // Closest index to (px, py), searching a window ahead of `hint` to stay on the right lap.
  int closest(double px, double py, int hint, int window = 200) const;
  // Last index whose arc length is within `distance_m` ahead of index i (clamped to the end).
  // Uses s, so callers do not depend on the sample spacing.
  int index_ahead(int i, double distance_m) const;
  double length() const { return s.back(); }
};

// ------------------------------------------------------------------ vehicle
// Kinematic bicycle (rear-axle reference) with a first-order, rate-limited steering actuator
// and a speed loop with acceleration limits. Integrated at 1 kHz.
struct TruckParams {
  double wheelbase{4.0};
  double max_steer{0.5};
  double steer_tau{0.2};       // actuator lag [s]
  double steer_rate{0.6};      // [rad/s]
  double a_max{2.0};           // [m/s^2]
};

struct TruckState { double x{0}, y{0}, th{0}, v{0}, steer{0}; };

class Truck {
 public:
  explicit Truck(const TruckParams& p) : p_(p) {}
  void step(double dt, double steer_cmd, double v_cmd);
  TruckState& state() { return s_; }
  const TruckState& state() const { return s_; }
  double yaw_rate() const { return s_.v * std::tan(s_.steer) / p_.wheelbase; }
  const TruckParams& params() const { return p_; }

 private:
  TruckParams p_;
  TruckState s_;
};

// ------------------------------------------------------------------ EKF
// x = [px, py, heading, speed, gyro bias]. Gyro and wheel speed at 100 Hz drive the
// prediction and a speed update; GPS position at 10 Hz corrects drift. GPS innovations are
// gated on the normalized innovation squared (chi^2, 2 dof).
class Ekf {
 public:
  using V5 = Eigen::Matrix<double, 5, 1>;
  using M5 = Eigen::Matrix<double, 5, 5>;
  struct Noise {
    double gyro{0.005}, accel{0.5}, bias_walk{2e-4};  // continuous-time densities
    double wheel{0.05}, gps{0.5};                       // measurement std
  };
  Ekf(const V5& x0, const M5& P0, const Noise& n, bool gate) : x_(x0), P_(P0), n_(n), gate_(gate) {}
  void predict(double gyro, double dt);
  void update_speed(double v_meas);
  // Returns the NIS; sets *accepted to false if the gate rejected the measurement.
  double update_gps(double px, double py, bool* accepted);
  const V5& x() const { return x_; }
  const M5& P() const { return P_; }
  long rejected() const { return rejected_; }

 private:
  V5 x_;
  M5 P_;
  Noise n_;
  bool gate_;
  long rejected_{0};
};

// ------------------------------------------------------------------ path trackers
struct TrackerInput {
  double x, y, th, v, steer;  // estimated pose and speed, measured steering angle
  const RefPath* path;
  int idx;                    // closest path index
};

class Tracker {
 public:
  virtual ~Tracker() = default;
  virtual std::string name() const = 0;
  virtual double steer(const TrackerInput& in) = 0;
};

std::unique_ptr<Tracker> make_tracker(const std::string& name, const TruckParams& tp, double dt);

// Signed lateral error (positive = vehicle left of the path) and heading error at index i.
void path_errors(const RefPath& p, int i, double x, double y, double th, double& e_lat, double& e_head);

// ------------------------------------------------------------------ closed-loop scenario
struct ScenarioConfig {
  std::string tracker{"lqr"};
  bool use_ekf{true};
  bool gps_gate{true};
  double gps_outlier_prob{0.02};
  double gps_outlier_mag{8.0};
  double dropout_t0{24.0}, dropout_t1{32.0};
  std::uint64_t seed{1};
  std::string log_path;
};

struct ScenarioResult {
  double rms_lat{0}, max_lat{0}, rms_head{0}, steer_rate_rms{0};
  double rms_pos_err{0}, max_pos_err{0};  // estimate vs truth
  bool completed{false};
  long gps_rejected{0}, gps_outliers{0};
  std::vector<double> nees;  // per 10 Hz GPS epoch (3 dof: x, y, heading)
  std::vector<double> nis;   // per accepted GPS update (2 dof)
};

ScenarioResult run_scenario(const RefPath& path, const ScenarioConfig& cfg);

}  // namespace pg
