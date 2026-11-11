#pragma once
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace pg {

struct Cell {
  int x{0}, y{0};
  bool operator==(const Cell& o) const { return x == o.x && y == o.y; }
  bool operator!=(const Cell& o) const { return !(*this == o); }
};

// Binary occupancy grid. Out-of-bounds cells count as occupied, so planners never need a
// separate bounds check before an occupancy query.
class Grid {
 public:
  Grid() = default;
  Grid(int width, int height, double resolution)
      : w_(width), h_(height), res_(resolution), occ_(static_cast<std::size_t>(width) * height, 0) {}

  int width() const { return w_; }
  int height() const { return h_; }
  double resolution() const { return res_; }
  int size() const { return w_ * h_; }
  bool inside(int x, int y) const { return x >= 0 && y >= 0 && x < w_ && y < h_; }
  int index(int x, int y) const { return y * w_ + x; }
  Cell cell(int idx) const { return {idx % w_, idx / w_}; }

  bool occupied(int x, int y) const { return !inside(x, y) || occ_[static_cast<std::size_t>(index(x, y))] != 0; }
  bool free(int x, int y) const { return !occupied(x, y); }
  void set(int x, int y, bool occ) {
    if (inside(x, y)) occ_[static_cast<std::size_t>(index(x, y))] = occ ? 1 : 0;
  }
  void fill_rect(int x0, int y0, int x1, int y1, bool occ = true);

  // World <-> grid (cell centres at (i + 0.5) * res).
  Cell to_cell(double wx, double wy) const {
    return {static_cast<int>(std::floor(wx / res_)), static_cast<int>(std::floor(wy / res_))};
  }
  double cx(int x) const { return (x + 0.5) * res_; }
  double cy(int y) const { return (y + 0.5) * res_; }
  bool occupied_world(double wx, double wy) const {
    const Cell c = to_cell(wx, wy);
    return occupied(c.x, c.y);
  }

  // Exact Euclidean distance (metres) from each cell centre to the nearest occupied cell
  // centre, Felzenszwalb-Huttenlocher separable transform, O(N).
  std::vector<float> distance_transform() const;
  // Occupy every cell whose centre is within `radius` of an obstacle (configuration space
  // for a disc robot).
  Grid inflate(double radius) const;

  static Grid random_rects(int w, int h, double res, int n_rects, int min_side, int max_side,
                           std::uint64_t seed);
  // A walled parking lot with rows of parked cars, for the kinodynamic planner.
  static Grid parking_lot(double res);
  bool save_pgm(const std::string& path) const;

 private:
  int w_{0}, h_{0};
  double res_{1.0};
  std::vector<std::uint8_t> occ_;
};

inline double octile(const Cell& a, const Cell& b) {
  const double dx = std::abs(a.x - b.x), dy = std::abs(a.y - b.y);
  return (dx + dy) + (std::sqrt(2.0) - 2.0) * std::min(dx, dy);
}

}  // namespace pg
