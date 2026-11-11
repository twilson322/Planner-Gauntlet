#include "pg/dstar_lite.hpp"

#include <algorithm>

namespace pg {

DStarLite::DStarLite(Grid belief, Cell start, Cell goal)
    : belief_(std::move(belief)), start_(start), goal_(goal), last_(start) {
  const auto n = static_cast<std::size_t>(belief_.size());
  g_.assign(n, kInf);
  rhs_.assign(n, kInf);
  open_key_.assign(n, Key{kInf, kInf});
  in_open_.assign(n, 0);
  rhs_[idx(goal_)] = 0.0;
  const Key k{h(start_, goal_), 0.0};
  open_key_[idx(goal_)] = k;
  in_open_[idx(goal_)] = 1;
  open_.push({k, static_cast<int>(idx(goal_))});
}

double DStarLite::cost(const Cell& a, const Cell& b) const {
  if (belief_.occupied(a.x, a.y) || belief_.occupied(b.x, b.y)) return kInf;
  const int dx = b.x - a.x, dy = b.y - a.y;
  if (dx && dy) {
    if (belief_.occupied(a.x + dx, a.y) || belief_.occupied(a.x, a.y + dy)) return kInf;
    return kSqrt2;
  }
  return 1.0;
}

DStarLite::Key DStarLite::calc_key(const Cell& s) const {
  const double m = std::min(g_[idx(s)], rhs_[idx(s)]);
  return {m + h(start_, s) + km_, m};
}

void DStarLite::update_vertex(const Cell& u) {
  const std::size_t ui = idx(u);
  if (u != goal_) {
    double best = kInf;
    for_each_adjacent(u, [&](const Cell& s) {
      const double c = cost(u, s);
      if (c < kInf) best = std::min(best, c + g_[idx(s)]);
    });
    rhs_[ui] = best;
  }
  in_open_[ui] = 0;  // lazy removal: stale heap entries are skipped when popped
  if (g_[ui] != rhs_[ui]) {
    const Key k = calc_key(u);
    open_key_[ui] = k;
    in_open_[ui] = 1;
    open_.push({k, static_cast<int>(ui)});
  }
}

bool DStarLite::top_key(Key& k) {
  while (!open_.empty()) {
    const Entry& e = open_.top();
    const auto i = static_cast<std::size_t>(e.idx);
    if (in_open_[i] && !(e.key < open_key_[i]) && !(open_key_[i] < e.key)) {
      k = e.key;
      return true;
    }
    open_.pop();  // stale
  }
  k = {kInf, kInf};
  return false;
}

void DStarLite::compute_shortest_path() {
  Key kt;
  while (top_key(kt) && (kt < calc_key(start_) || rhs_[idx(start_)] != g_[idx(start_)])) {
    const int ui = open_.top().idx;
    const Cell u = belief_.cell(ui);
    const Key k_old = kt;
    const Key k_new = calc_key(u);
    ++expansions_;
    if (k_old < k_new) {  // key grew because km moved on: re-insert with the current key
      open_.pop();
      open_key_[static_cast<std::size_t>(ui)] = k_new;
      open_.push({k_new, ui});
    } else if (g_[static_cast<std::size_t>(ui)] > rhs_[static_cast<std::size_t>(ui)]) {  // overconsistent
      open_.pop();
      in_open_[static_cast<std::size_t>(ui)] = 0;
      g_[static_cast<std::size_t>(ui)] = rhs_[static_cast<std::size_t>(ui)];
      for_each_adjacent(u, [&](const Cell& s) { update_vertex(s); });
    } else {  // underconsistent
      open_.pop();
      in_open_[static_cast<std::size_t>(ui)] = 0;
      g_[static_cast<std::size_t>(ui)] = kInf;
      update_vertex(u);
      for_each_adjacent(u, [&](const Cell& s) { update_vertex(s); });
    }
  }
}

void DStarLite::move_to(const Cell& s) { start_ = s; }

bool DStarLite::update_cells(const std::vector<std::pair<Cell, bool>>& changes) {
  // km grows by the heuristic distance from where the robot was at the *previous* repair,
  // not by the sum of every step since. Summing per step overestimates (triangle inequality),
  // which breaks the lower-bound property of the old keys and the search can stop early.
  km_ += h(last_, start_);
  last_ = start_;
  bool any = false;
  std::vector<Cell> touched;
  for (const auto& [c, occ] : changes) {
    if (belief_.occupied(c.x, c.y) == occ || !belief_.inside(c.x, c.y)) continue;
    belief_.set(c.x, c.y, occ);
    any = true;
    // Every edge whose cost depends on c (including diagonals that pass its corner) has
    // both endpoints inside the 3x3 block around c.
    for (int dy = -1; dy <= 1; ++dy)
      for (int dx = -1; dx <= 1; ++dx)
        if (belief_.inside(c.x + dx, c.y + dy)) touched.push_back({c.x + dx, c.y + dy});
  }
  for (const Cell& u : touched) update_vertex(u);
  return any;
}

Cell DStarLite::next_step() const {
  Cell best = start_;
  double bv = kInf;
  for_each_adjacent(start_, [&](const Cell& s) {
    const double v = cost(start_, s) + g_[idx(s)];
    if (v < bv) { bv = v; best = s; }
  });
  return best;
}

std::vector<Cell> DStarLite::extract_path(int max_steps) const {
  std::vector<Cell> p{start_};
  if (g_[idx(start_)] == kInf) return {};
  Cell c = start_;
  for (int i = 0; i < max_steps && c != goal_; ++i) {
    Cell best = c;
    double bv = kInf;
    for_each_adjacent(c, [&](const Cell& s) {
      const double v = cost(c, s) + g_[idx(s)];
      if (v < bv) { bv = v; best = s; }
    });
    if (best == c) return {};
    c = best;
    p.push_back(c);
  }
  return p;
}

}  // namespace pg
