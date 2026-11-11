#include "pg/tracking.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>

namespace pg {

// ============================================================================ reference path
RefPath RefPath::road(double v_max, double a_lat_max, double a_lon_max) {
  // Piecewise curvature profile: straight, 100 m-radius left curve with clothoid transitions,
  // a 3.7 m lane change, and a 40 m-radius right curve.
  struct Seg { double s0, s1; int kind; double a, b; };  // kind 0 const, 1 ramp a->b, 2 sine
  const std::vector<Seg> segs = {
      {0, 80, 0, 0, 0},          {80, 110, 1, 0, 0.01},     {110, 230, 0, 0.01, 0},
      {230, 260, 1, 0.01, 0},    {260, 320, 0, 0, 0},       {320, 380, 2, 0.00646, 60},
      {380, 430, 0, 0, 0},       {430, 455, 1, 0, -0.025},  {455, 520, 0, -0.025, 0},
      {520, 545, 1, -0.025, 0},  {545, 640, 0, 0, 0}};
  auto kappa = [&](double s) {
    for (const auto& g : segs) {
      if (s < g.s0 || s > g.s1) continue;
      if (g.kind == 0) return g.a;
      if (g.kind == 1) return g.a + (g.b - g.a) * (s - g.s0) / (g.s1 - g.s0);
      return g.a * std::sin(2 * M_PI * (s - g.s0) / g.b);
    }
    return 0.0;
  };
  RefPath p;
  const double ds = 0.5;
  double x = 0, y = 0, th = 0;
  for (double s = 0; s <= 640.0 + 1e-9; s += ds) {
    p.s.push_back(s);
    p.x.push_back(x);
    p.y.push_back(y);
    p.th.push_back(th);
    p.kappa.push_back(kappa(s));
    const double km = kappa(s + 0.5 * ds);  // midpoint rule
    const double thm = th + 0.5 * ds * kappa(s);
    x += ds * std::cos(thm);
    y += ds * std::sin(thm);
    th += ds * km;
  }
  const std::size_t n = p.s.size();
  p.v.resize(n);
  for (std::size_t i = 0; i < n; ++i)
    p.v[i] = std::min(v_max, std::sqrt(a_lat_max / std::max(1e-6, std::abs(p.kappa[i]))));
  for (std::size_t i = 1; i < n; ++i) p.v[i] = std::min(p.v[i], std::sqrt(p.v[i - 1] * p.v[i - 1] + 2 * a_lon_max * ds));
  for (std::size_t i = n - 1; i-- > 0;) p.v[i] = std::min(p.v[i], std::sqrt(p.v[i + 1] * p.v[i + 1] + 2 * a_lon_max * ds));
  return p;
}

int RefPath::index_ahead(int i, double distance_m) const {
  const double target = s[static_cast<std::size_t>(i)] + distance_m;
  const auto after = std::upper_bound(s.begin(), s.end(), target);
  return std::max(i, static_cast<int>(after - s.begin()) - 1);
}

int RefPath::closest(double px, double py, int hint, int window) const {
  const int n = static_cast<int>(s.size());
  int best = std::clamp(hint, 0, n - 1);
  double bd = 1e18;
  for (int i = std::max(0, hint - 20); i < std::min(n, hint + window); ++i) {
    const double d = (x[static_cast<std::size_t>(i)] - px) * (x[static_cast<std::size_t>(i)] - px) +
                     (y[static_cast<std::size_t>(i)] - py) * (y[static_cast<std::size_t>(i)] - py);
    if (d < bd) { bd = d; best = i; }
  }
  return best;
}

void path_errors(const RefPath& p, int i, double x, double y, double th, double& e_lat, double& e_head) {
  const auto k = static_cast<std::size_t>(i);
  const double dx = x - p.x[k], dy = y - p.y[k];
  e_lat = -std::sin(p.th[k]) * dx + std::cos(p.th[k]) * dy;  // left of path is positive
  e_head = wrap_pi(th - p.th[k]);
}

// ============================================================================ vehicle
void Truck::step(double dt, double steer_cmd, double v_cmd) {
  const double target = std::clamp(steer_cmd, -p_.max_steer, p_.max_steer);
  const double rate = std::clamp((target - s_.steer) / p_.steer_tau, -p_.steer_rate, p_.steer_rate);
  s_.steer = std::clamp(s_.steer + rate * dt, -p_.max_steer, p_.max_steer);
  const double a = std::clamp(2.0 * (v_cmd - s_.v), -p_.a_max, p_.a_max);
  s_.v = std::max(0.0, s_.v + a * dt);
  s_.x += s_.v * std::cos(s_.th) * dt;
  s_.y += s_.v * std::sin(s_.th) * dt;
  s_.th = wrap_pi(s_.th + s_.v * std::tan(s_.steer) / p_.wheelbase * dt);
}

// ============================================================================ EKF
void Ekf::predict(double gyro, double dt) {
  const double c = std::cos(x_(2)), s = std::sin(x_(2)), v = x_(3);
  x_(0) += v * c * dt;
  x_(1) += v * s * dt;
  x_(2) = wrap_pi(x_(2) + (gyro - x_(4)) * dt);
  M5 F = M5::Identity();
  F(0, 2) = -v * s * dt;
  F(0, 3) = c * dt;
  F(1, 2) = v * c * dt;
  F(1, 3) = s * dt;
  F(2, 4) = -dt;
  M5 Q = M5::Zero();
  Q(2, 2) = n_.gyro * n_.gyro * dt * dt;  // per-sample gyro noise integrated over one step
  Q(3, 3) = n_.accel * n_.accel * dt;     // speed as a random walk
  Q(4, 4) = n_.bias_walk * n_.bias_walk * dt;
  P_ = F * P_ * F.transpose() + Q;
}

void Ekf::update_speed(double v_meas) {
  const double S = P_(3, 3) + n_.wheel * n_.wheel;
  const V5 K = P_.col(3) / S;
  x_ += K * (v_meas - x_(3));
  M5 IKH = M5::Identity();
  IKH.col(3) -= K;
  P_ = IKH * P_ * IKH.transpose() + K * K.transpose() * n_.wheel * n_.wheel;  // Joseph form
}

double Ekf::update_gps(double px, double py, bool* accepted) {
  Eigen::Matrix<double, 2, 5> H = Eigen::Matrix<double, 2, 5>::Zero();
  H(0, 0) = 1.0;
  H(1, 1) = 1.0;
  const Eigen::Matrix2d R = Eigen::Matrix2d::Identity() * n_.gps * n_.gps;
  const Eigen::Vector2d nu(px - x_(0), py - x_(1));
  const Eigen::Matrix2d S = H * P_ * H.transpose() + R;
  const double nis = nu.dot(S.ldlt().solve(nu));
  // chi^2(2) at 99.9 %: a consistent filter rejects 1 in 1000 good fixes.
  if (gate_ && nis > 13.82) {
    ++rejected_;
    if (accepted) *accepted = false;
    return nis;
  }
  const Eigen::Matrix<double, 5, 2> K = P_ * H.transpose() * S.inverse();
  x_ += K * nu;
  x_(2) = wrap_pi(x_(2));
  const M5 IKH = M5::Identity() - K * H;
  P_ = IKH * P_ * IKH.transpose() + K * R * K.transpose();
  if (accepted) *accepted = true;
  return nis;
}

// ============================================================================ trackers
namespace {

class PurePursuit final : public Tracker {
 public:
  explicit PurePursuit(const TruckParams& p) : L_(p.wheelbase) {}
  std::string name() const override { return "pure_pursuit"; }
  double steer(const TrackerInput& in) override {
    // Lookahead grows with speed (0.9 s of travel) but never below 6 m, where the 4 m wheelbase
    // truck starts to weave. The 1 m floor on the chord only guards the division.
    constexpr double kMinLookaheadM = 6.0, kLookaheadTimeS = 0.9, kMinChordM = 1.0;
    const double Ld = std::max(kMinLookaheadM, kLookaheadTimeS * in.v);
    const int ti = in.path->index_ahead(in.idx, Ld);
    const double dx = in.path->x[static_cast<std::size_t>(ti)] - in.x, dy = in.path->y[static_cast<std::size_t>(ti)] - in.y;
    const double alpha = wrap_pi(std::atan2(dy, dx) - in.th);
    const double ld = std::max(kMinChordM, std::hypot(dx, dy));
    return std::atan(2.0 * L_ * std::sin(alpha) / ld);
  }

 private:
  double L_;
};

class Stanley final : public Tracker {
 public:
  explicit Stanley(const TruckParams& p) : L_(p.wheelbase) {}
  std::string name() const override { return "stanley"; }
  double steer(const TrackerInput& in) override {
    // Stanley uses the front-axle error.
    const double xf = in.x + L_ * std::cos(in.th), yf = in.y + L_ * std::sin(in.th);
    fi_ = in.path->closest(xf, yf, std::max(fi_, in.idx));
    double e, eh;
    path_errors(*in.path, fi_, xf, yf, in.th, e, eh);
    return -eh - std::atan2(k_ * e, in.v + kSoftening_m_s);
  }

 private:
  double L_;
  double k_{0.8};                        // cross-track gain, 1/s
  static constexpr double kSoftening_m_s = 1.0;  // keeps the gain finite at low speed
  int fi_{0};
};

Eigen::MatrixXd expm(const Eigen::MatrixXd& M) {
  const double norm = M.cwiseAbs().rowwise().sum().maxCoeff();
  const int s = norm > 0.5 ? static_cast<int>(std::ceil(std::log2(norm / 0.5))) : 0;
  const Eigen::MatrixXd X = M / std::pow(2.0, s);
  Eigen::MatrixXd E = Eigen::MatrixXd::Identity(M.rows(), M.cols()), term = E;
  for (int k = 1; k <= 14; ++k) {
    term = term * X / static_cast<double>(k);
    E += term;
  }
  for (int i = 0; i < s; ++i) E = E * E;
  return E;
}

// Lateral LQR on [e, heading error, steering angle] with the steering actuator lag in the
// model, plus curvature feedforward. Gains are re-solved per 0.5 m/s speed bin and cached.
// Weights: Q = diag(0.5 per m^2, 4 per rad^2, 0), R = 8 per rad^2 of steering command.
class LqrTracker final : public Tracker {
 public:
  LqrTracker(const TruckParams& p, double dt) : p_(p), dt_(dt) {}
  std::string name() const override { return "lqr"; }
  double steer(const TrackerInput& in) override {
    double e, eh;
    path_errors(*in.path, in.idx, in.x, in.y, in.th, e, eh);
    const double kap = in.path->kappa[static_cast<std::size_t>(in.idx)];
    const double dff = std::atan(p_.wheelbase * kap);
    const Eigen::RowVector3d K = gain(std::max(1.0, in.v));
    return dff - K.dot(Eigen::Vector3d(e, eh, in.steer - dff));
  }

 private:
  Eigen::RowVector3d gain(double v) {
    constexpr double kSpeedBinMPerS = 0.5;
    const int bin = static_cast<int>(std::lround(v / kSpeedBinMPerS));
    auto it = cache_.find(bin);
    if (it != cache_.end()) return it->second;
    const double vb = bin * kSpeedBinMPerS;
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(4, 4);
    M(0, 1) = vb;
    M(1, 2) = vb / p_.wheelbase;
    M(2, 2) = -1.0 / p_.steer_tau;
    M(2, 3) = 1.0 / p_.steer_tau;
    const Eigen::MatrixXd E = expm(M * dt_);
    const Eigen::Matrix3d A = E.topLeftCorner(3, 3);
    const Eigen::Vector3d B = E.topRightCorner(3, 1);
    const Eigen::Matrix3d Q = Eigen::Vector3d(0.5, 4.0, 0.0).asDiagonal();
    const double R = 8.0;
    Eigen::Matrix3d P = Q;
    for (int i = 0; i < 5000; ++i) {  // small system at 50 Hz: plain Riccati iteration converges
      const double S = R + B.dot(P * B);
      const Eigen::Matrix3d Pn = A.transpose() * P * A - (A.transpose() * P * B) * (B.transpose() * P * A) / S + Q;
      if ((Pn - P).norm() < 1e-10 * Pn.norm()) { P = Pn; break; }
      P = Pn;
    }
    const Eigen::RowVector3d K = (B.transpose() * P * A) / (R + B.dot(P * B));
    cache_[bin] = K;
    return K;
  }

  TruckParams p_;
  double dt_;
  std::map<int, Eigen::RowVector3d> cache_;
};

}  // namespace

std::unique_ptr<Tracker> make_tracker(const std::string& name, const TruckParams& tp, double dt) {
  if (name == "pure_pursuit") return std::make_unique<PurePursuit>(tp);
  if (name == "stanley") return std::make_unique<Stanley>(tp);
  if (name == "lqr") return std::make_unique<LqrTracker>(tp, dt);
  return nullptr;
}

// ============================================================================ scenario
ScenarioResult run_scenario(const RefPath& path, const ScenarioConfig& cfg) {
  // Multirate loop on a 1 kHz plant: IMU and wheel speed at 100 Hz, GPS at 10 Hz, control at 50 Hz.
  constexpr double dt = 0.001;
  constexpr int kImuEveryTicks = 10, kGpsEveryTicks = 100, kControlEveryTicks = 20;
  constexpr double kImuPeriodS = kImuEveryTicks * dt, kControlPeriodS = kControlEveryTicks * dt;
  constexpr long kTimeoutTicks = 200000;           // 200 s, far longer than any run
  constexpr double kTrueGyroBiasWalk = 1e-4;        // rad/s per sqrt(s); the EKF assumes 2e-4
  constexpr double kSpeedTargetAheadM = 5.0;
  constexpr double kScoreAfterS = 2.0;               // the initial 0.5 m offset is gone by then
  constexpr double kFinishWithinM = 1.5;
  ScenarioResult R;
  const TruckParams tp;
  Truck truck(tp);
  truck.state() = {path.x[0], path.y[0] + 0.5, path.th[0], path.v[0], 0.0};  // start 0.5 m off
  std::mt19937_64 rng(cfg.seed);
  std::normal_distribution<double> N01(0.0, 1.0);
  std::uniform_real_distribution<double> U(0.0, 1.0);

  const Ekf::Noise nz;
  double bias_true = 0.01;  // rad/s gyro bias (0.57 deg/s, a cheap MEMS part)
  Ekf::V5 x0;
  x0 << truck.state().x + nz.gps * N01(rng), truck.state().y + nz.gps * N01(rng), truck.state().th + 0.02 * N01(rng),
      truck.state().v, 0.0;
  Ekf::M5 P0 = Eigen::Matrix<double, 5, 1>(nz.gps * nz.gps, nz.gps * nz.gps, 4e-4, 0.01, 4e-4).asDiagonal();
  Ekf ekf(x0, P0, nz, cfg.gps_gate);

  auto tracker = make_tracker(cfg.tracker, tp, kControlPeriodS);
  int idx = 0, est_idx = 0;
  double steer_cmd = 0.0, prev_steer = 0.0;
  double sum_lat2 = 0, sum_head2 = 0, sum_sr2 = 0, sum_pos2 = 0;
  long n_m = 0, n_sr = 0;
  std::FILE* log = cfg.log_path.empty() ? nullptr : std::fopen(cfg.log_path.c_str(), "w");
  if (log) std::fprintf(log, "t,x,y,xh,yh,e_lat,steer,v,gps_x,gps_y,gps_used\n");
  double last_gx = NAN, last_gy = NAN;
  int last_used = -1;

  for (long k = 0; k < kTimeoutTicks; ++k) {
    const double t = k * dt;
    const TruckState& s = truck.state();
    idx = path.closest(s.x, s.y, idx);
    if (path.s[static_cast<std::size_t>(idx)] >= path.length() - kFinishWithinM) { R.completed = true; break; }

    if (k % kImuEveryTicks == 0) {
      bias_true += kTrueGyroBiasWalk * std::sqrt(kImuPeriodS) * N01(rng);
      const double gyro = truck.yaw_rate() + bias_true + nz.gyro * N01(rng);
      ekf.predict(gyro, kImuPeriodS);
      ekf.update_speed(s.v + nz.wheel * N01(rng));
    }
    if (k % kGpsEveryTicks == 0) {
      const bool dropout = t >= cfg.dropout_t0 && t < cfg.dropout_t1;
      if (!dropout) {
        double gx = s.x + nz.gps * N01(rng), gy = s.y + nz.gps * N01(rng);
        if (U(rng) < cfg.gps_outlier_prob) {  // multipath jump
          const double a = 2 * M_PI * U(rng);
          gx += cfg.gps_outlier_mag * std::cos(a);
          gy += cfg.gps_outlier_mag * std::sin(a);
          ++R.gps_outliers;
        }
        bool acc = true;
        const double nis = ekf.update_gps(gx, gy, &acc);
        if (acc) R.nis.push_back(nis);
        last_gx = gx;
        last_gy = gy;
        last_used = acc ? 1 : 0;
      }
      // NEES on (x, y, heading)
      const auto& xh = ekf.x();
      const Eigen::Vector3d e(s.x - xh(0), s.y - xh(1), wrap_pi(s.th - xh(2)));
      Eigen::Matrix3d P3;
      P3 << ekf.P()(0, 0), ekf.P()(0, 1), ekf.P()(0, 2), ekf.P()(1, 0), ekf.P()(1, 1), ekf.P()(1, 2),
          ekf.P()(2, 0), ekf.P()(2, 1), ekf.P()(2, 2);
      R.nees.push_back(e.dot(P3.ldlt().solve(e)));
    }
    if (k % kControlEveryTicks == 0) {
      TrackerInput in{};
      if (cfg.use_ekf) {
        const auto& xh = ekf.x();
        est_idx = path.closest(xh(0), xh(1), est_idx);
        in = {xh(0), xh(1), xh(2), xh(3), s.steer, &path, est_idx};
      } else {
        in = {s.x, s.y, s.th, s.v, s.steer, &path, idx};
      }
      steer_cmd = tracker->steer(in);
    }
    const int vi = path.index_ahead(idx, kSpeedTargetAheadM);
    truck.step(dt, steer_cmd, path.v[static_cast<std::size_t>(vi)]);

    if (k % kControlEveryTicks == 0 && t > kScoreAfterS) {
      double e, eh;
      path_errors(path, idx, s.x, s.y, s.th, e, eh);
      sum_lat2 += e * e;
      sum_head2 += eh * eh;
      R.max_lat = std::max(R.max_lat, std::abs(e));
      const double pe = std::hypot(s.x - ekf.x()(0), s.y - ekf.x()(1));
      sum_pos2 += pe * pe;
      R.max_pos_err = std::max(R.max_pos_err, pe);
      ++n_m;
      const double sr = (s.steer - prev_steer) / kControlPeriodS;
      sum_sr2 += sr * sr;
      ++n_sr;
      if (log)
        std::fprintf(log, "%.2f,%.3f,%.3f,%.3f,%.3f,%.4f,%.4f,%.3f,%.3f,%.3f,%d\n", t, s.x, s.y, ekf.x()(0),
                     ekf.x()(1), e, s.steer, s.v, last_gx, last_gy, last_used);
      last_used = -1;
    }
    if (k % kControlEveryTicks == 0) prev_steer = s.steer;
  }
  if (log) std::fclose(log);
  if (n_m) {
    R.rms_lat = std::sqrt(sum_lat2 / n_m);
    R.rms_head = std::sqrt(sum_head2 / n_m);
    R.rms_pos_err = std::sqrt(sum_pos2 / n_m);
  }
  if (n_sr) R.steer_rate_rms = std::sqrt(sum_sr2 / n_sr);
  R.gps_rejected = ekf.rejected();
  return R;
}

}  // namespace pg
