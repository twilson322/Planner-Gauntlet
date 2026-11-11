// Google Benchmark micro-benchmarks for the planners.
#include <benchmark/benchmark.h>

#include <cstdio>

#include "pg/dstar_lite.hpp"
#include "pg/hybrid_astar.hpp"
#include "pg/search.hpp"

using namespace pg;

namespace {
// Random-rectangle maps with the corners cleared; the query is corner to corner, and
// main() checks it is solvable so no benchmark silently measures an exhaustive failure.
Grid make_map(int n) {
  for (std::uint64_t seed = 9;; ++seed) {  // first seed whose corners are connected
    Grid g = Grid::random_rects(n, n, 1.0, n / 2, 2, n / 12, seed);
    g.fill_rect(0, 0, 4, 4, false);
    g.fill_rect(n - 5, n - 5, n - 1, n - 1, false);
    if (astar(GridGraph8(g), Cell{1, 1}, Cell{n - 2, n - 2}, OctileHeuristic{}).found) return g;
  }
}
const Grid& map_of(int n) {
  static const Grid g512 = make_map(512);
  static const Grid g1024 = make_map(1024);
  return n == 512 ? g512 : g1024;
}
Cell corner(const Grid& g, bool far) { return far ? Cell{g.width() - 2, g.height() - 2} : Cell{1, 1}; }
}  // namespace

static void BM_Dijkstra(benchmark::State& st) {
  const Grid& g = map_of(static_cast<int>(st.range(0)));
  const GridGraph8 G(g);
  long exp = 0;
  for (auto _ : st) exp = astar(G, corner(g, false), corner(g, true), ZeroHeuristic{}).expansions;
  st.counters["expansions"] = static_cast<double>(exp);
}
static void BM_AStar(benchmark::State& st) {
  const Grid& g = map_of(static_cast<int>(st.range(0)));
  const GridGraph8 G(g);
  long exp = 0;
  for (auto _ : st) exp = astar(G, corner(g, false), corner(g, true), OctileHeuristic{}).expansions;
  st.counters["expansions"] = static_cast<double>(exp);
}
static void BM_WeightedAStar3(benchmark::State& st) {
  const Grid& g = map_of(static_cast<int>(st.range(0)));
  const GridGraph8 G(g);
  long exp = 0;
  for (auto _ : st) exp = astar(G, corner(g, false), corner(g, true), OctileHeuristic{}, 3.0).expansions;
  st.counters["expansions"] = static_cast<double>(exp);
}
static void BM_JPS(benchmark::State& st) {
  const Grid& g = map_of(static_cast<int>(st.range(0)));
  long exp = 0;
  for (auto _ : st) exp = jps(g, corner(g, false), corner(g, true)).expansions;
  st.counters["expansions"] = static_cast<double>(exp);
}
static void BM_DistanceTransform(benchmark::State& st) {
  const Grid& g = map_of(static_cast<int>(st.range(0)));
  for (auto _ : st) benchmark::DoNotOptimize(g.distance_transform());
}
static void BM_HybridAStarParking(benchmark::State& st) {
  const Grid lot = Grid::parking_lot(0.1);
  for (auto _ : st)
    benchmark::DoNotOptimize(hybrid_astar(lot, Vehicle{}, {6.0, 30.0, 0.0}, {39.25, 4.85, -M_PI / 2}, HybridConfig{}));
}
BENCHMARK(BM_Dijkstra)->Arg(512)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_AStar)->Arg(512)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_WeightedAStar3)->Arg(512)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_JPS)->Arg(512)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_DistanceTransform)->Arg(512)->Arg(1024)->Unit(benchmark::kMillisecond);
BENCHMARK(BM_HybridAStarParking)->Unit(benchmark::kMillisecond)->Iterations(5);
int main(int argc, char** argv) {
  for (int n : {512, 1024}) {
    const Grid& g = map_of(n);
    if (!astar(GridGraph8(g), corner(g, false), corner(g, true), OctileHeuristic{}).found) {
      std::fprintf(stderr, "benchmark map %d has no corner-to-corner path\n", n);
      return 1;
    }
  }
  benchmark::Initialize(&argc, argv);
  benchmark::RunSpecifiedBenchmarks();
  benchmark::Shutdown();
  return 0;
}
