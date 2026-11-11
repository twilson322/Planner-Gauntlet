#include "pg/grid.hpp"

#include <algorithm>
#include <cstdio>
#include <limits>
#include <random>

namespace pg {

void Grid::fill_rect(int x0, int y0, int x1, int y1, bool occ) {
  for (int y = std::max(0, y0); y <= std::min(h_ - 1, y1); ++y)
    for (int x = std::max(0, x0); x <= std::min(w_ - 1, x1); ++x) set(x, y, occ);
}

namespace {
// "Infinity" for the transform. A large finite value keeps the parabola-intersection
// arithmetic well defined; real infinities turn it into inf - inf = NaN.
constexpr double kFar = 1e20;

// 1D squared-distance transform (Felzenszwalb & Huttenlocher 2012).
void dt1d(const std::vector<double>& f, std::vector<double>& d, std::vector<int>& v,
          std::vector<double>& z, int n) {
  int k = 0;
  v[0] = 0;
  z[0] = -kFar;
  z[1] = kFar;
  for (int q = 1; q < n; ++q) {
    double s = 0.0;
    for (;;) {
      const int p = v[static_cast<std::size_t>(k)];
      s = ((f[static_cast<std::size_t>(q)] + double(q) * q) - (f[static_cast<std::size_t>(p)] + double(p) * p)) /
          (2.0 * q - 2.0 * p);
      if (s <= z[static_cast<std::size_t>(k)]) {
        --k;
        continue;
      }
      break;
    }
    ++k;
    v[static_cast<std::size_t>(k)] = q;
    z[static_cast<std::size_t>(k)] = s;
    z[static_cast<std::size_t>(k + 1)] = kFar;
  }
  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[static_cast<std::size_t>(k + 1)] < q) ++k;
    const int p = v[static_cast<std::size_t>(k)];
    d[static_cast<std::size_t>(q)] = double(q - p) * (q - p) + f[static_cast<std::size_t>(p)];
  }
}
}  // namespace

std::vector<float> Grid::distance_transform() const {
  const int n = std::max(w_, h_);
  std::vector<double> f(static_cast<std::size_t>(n)), d(static_cast<std::size_t>(n)), z(static_cast<std::size_t>(n + 1));
  std::vector<int> v(static_cast<std::size_t>(n));
  std::vector<double> tmp(static_cast<std::size_t>(size()));
  // columns
  for (int x = 0; x < w_; ++x) {
    f.assign(static_cast<std::size_t>(h_), 0.0);
    d.assign(static_cast<std::size_t>(h_), 0.0);
    for (int y = 0; y < h_; ++y) f[static_cast<std::size_t>(y)] = occupied(x, y) ? 0.0 : kFar;
    dt1d(f, d, v, z, h_);
    for (int y = 0; y < h_; ++y) tmp[static_cast<std::size_t>(index(x, y))] = d[static_cast<std::size_t>(y)];
  }
  // rows
  std::vector<float> out(static_cast<std::size_t>(size()));
  for (int y = 0; y < h_; ++y) {
    f.assign(static_cast<std::size_t>(w_), 0.0);
    d.assign(static_cast<std::size_t>(w_), 0.0);
    for (int x = 0; x < w_; ++x) f[static_cast<std::size_t>(x)] = tmp[static_cast<std::size_t>(index(x, y))];
    dt1d(f, d, v, z, w_);
    for (int x = 0; x < w_; ++x)
      out[static_cast<std::size_t>(index(x, y))] = static_cast<float>(std::sqrt(d[static_cast<std::size_t>(x)]) * res_);
  }
  return out;
}

Grid Grid::inflate(double radius) const {
  const auto dt = distance_transform();
  Grid g(w_, h_, res_);
  for (int i = 0; i < size(); ++i) {
    const Cell c = cell(i);
    g.set(c.x, c.y, dt[static_cast<std::size_t>(i)] <= radius);
  }
  return g;
}

Grid Grid::random_rects(int w, int h, double res, int n_rects, int min_side, int max_side,
                        std::uint64_t seed) {
  Grid g(w, h, res);
  std::mt19937_64 rng(seed);
  std::uniform_int_distribution<int> X(0, w - 1), Y(0, h - 1), S(min_side, max_side);
  for (int i = 0; i < n_rects; ++i) {
    const int x = X(rng), y = Y(rng);
    g.fill_rect(x, y, x + S(rng), y + S(rng));
  }
  return g;
}

Grid Grid::parking_lot(double res) {
  // 60 m x 40 m lot, 2.5 m x 5 m stalls, two rows of parked cars with one empty stall.
  const int W = static_cast<int>(60.0 / res), H = static_cast<int>(40.0 / res);
  Grid g(W, H, res);
  auto rect_m = [&](double x0, double y0, double x1, double y1) {
    g.fill_rect(static_cast<int>(x0 / res), static_cast<int>(y0 / res),
                static_cast<int>(x1 / res) - 1, static_cast<int>(y1 / res) - 1);
  };
  rect_m(0, 0, 60, 0.4);  // walls
  rect_m(0, 39.6, 60, 40);
  rect_m(0, 0, 0.4, 40);
  rect_m(59.6, 0, 60, 40);
  for (int row = 0; row < 2; ++row) {
    const double y0 = row == 0 ? 0.4 : 34.6;
    for (int s = 0; s < 22; ++s) {
      if (row == 0 && s == 14) continue;  // the target stall
      const double x0 = 3.0 + s * 2.5;
      rect_m(x0 + 0.3, y0 + 0.2, x0 + 2.2, y0 + 4.8);  // a parked car
    }
  }
  rect_m(10, 15, 50, 25);  // central island
  return g;
}

bool Grid::save_pgm(const std::string& path) const {
  std::FILE* f = std::fopen(path.c_str(), "wb");
  if (!f) return false;
  std::fprintf(f, "P5\n%d %d\n255\n", w_, h_);
  for (int y = h_ - 1; y >= 0; --y)
    for (int x = 0; x < w_; ++x) std::fputc(occupied(x, y) ? 0 : 255, f);
  std::fclose(f);
  return true;
}

}  // namespace pg
