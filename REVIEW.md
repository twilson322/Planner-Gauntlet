# Code Review: planner-gauntlet

Reviewed as a staff engineer on an autonomy team would review a pull request that adds grid search, Hybrid A*, RRT*, D* Lite, an EKF and three path trackers. Scope: `include/pg`, `src`, the experiment and benchmark apps, tests and build.

## 1. Overall assessment

Careful, well-tested planning code. The search algorithms are verified against each other the way they should be (JPS and A* on the same corner-cutting-free graph, D* Lite against a from-scratch A* after every map change, weighted A* within its bound), and the EKF reports NIS and NEES rather than just RMS error.

The review found no algorithmic defect. It found a hidden assumption in the tracking scenario: three places (pure-pursuit lookahead, the speed target and the finish test) converted metres to sample indices by assuming the reference path is sampled every 0.5 m. Any other spacing would have silently changed the lookahead distance. Those now use arc length. The remaining changes name constants that were bare literals in the trackers, the multirate loop and Hybrid A*. Every experiment was rerun after the changes: all 24 output files are identical except for wall-clock timing columns.

## 2. Scorecard

| Area | Score | Basis |
|---|---|---|
| Correctness | 8 | Cross-checked search algorithms, consistent heuristics, D* Lite equal to from-scratch replanning, a Joseph-form EKF with a chi-squared gate. |
| Architecture | 8 | Generic `astar` over any graph type, separate graph, heuristic and planner types, tracker interface with a factory. |
| Readability | 7 | Readable, but `run_scenario` mixes simulation, estimation, control and scoring in one 100-line loop. |
| Self-documentation | 8 | After this pass, loop rates, gains and thresholds are named with units and reasons. The map geometry in `grid.cpp` is literal coordinates, which is acceptable for scenario data. |
| Maintainability | 7 | Deterministic outputs make regressions easy to see with `cmp`, which is how this review checked itself. |
| Testing | 8 | 15 tests, including cross-algorithm equivalence, the weighted-A* bound and monotone RRT* cost. The LQR and gated-EKF results are pinned by tests; pure pursuit and Stanley are not. |
| Performance | 8 | JPS and lazy-deletion A* benchmarked with Google Benchmark. Hybrid A* is 167 to 391 ms on its map, reported honestly. |
| Reliability | 7 | Expansion caps on every search, finite timeouts. `fopen` for the optional log is checked. |
| API design | 8 | Config structs with units in comments, `SearchResult` with expansions, `index_ahead` added so callers stop doing spacing arithmetic. |
| Domain-specific correctness | 8 | Corner cutting forbidden on 8-connected grids, ZOH discretization of the steering lag by matrix exponential, curvature feedforward, rate-limited actuator. |
| Safety and robustness | Not assessed | Offline planning and a simulated truck, no safety function. |

## 3. Problems

### Fixed in this review

| Severity | Location | What | Why | Fix and evidence |
|---|---|---|---|---|
| Medium | `PurePursuit::steer`, `run_scenario` | Lookahead `idx + Ld / 0.5`, speed target `idx + 10`, finish test `idx >= n - 4`: metres converted to indices assuming 0.5 m path spacing. | A path built with any other spacing would change the lookahead, the speed preview and the finish line with no error. | `RefPath::index_ahead(i, metres)` on arc length, finish test in metres. Outputs identical on the 0.5 m path. New test with deliberately uneven spacing. |
| Low | `run_scenario` | Loop rates as `k % 10`, `k % 100`, `k % 20`, the control period as `0.02` in two places, the true gyro bias walk as `1e-4 * sqrt(0.01)`. | A rate change in one place would desynchronize the others. | Named rates and periods derived from one `dt`. |
| Low | trackers | Pure-pursuit lookahead (6 m, 0.9 s), Stanley gain and softening, LQR weights and speed bin. | Tuning without units or reasons. | Named constants with comments. |
| Low | `hybrid_astar.cpp` | 0.1 m footprint-check spacing in three places, 10 m analytic-shot radius. | The three spacings could drift apart. | `HybridConfig::collision_check_step`, `analytic_always_within`. |

### Open

| Severity | Location | What | Why | Recommended fix |
|---|---|---|---|---|
| Medium | `run_scenario` | Simulation, estimation, control and scoring in one function. | Hard to test the EKF-tracker interaction separately. Only the LQR and gated-EKF results have regression tests. | Split the loop, and pin the pure-pursuit and Stanley metrics within tolerances too. |
| Low | `RefPath::closest` | Windowed search (20 behind, 200 ahead of a hint). | Correct here, but a large estimate jump (a GPS outlier accepted with the gate off) could lock onto the wrong part of a self-intersecting path. | Document the assumption, or fall back to a global search when the windowed distance is large. |
| Low | RRT* | The rewiring constant uses the whole map area instead of the free area. | Larger than necessary, so it stays above the asymptotic-optimality bound. Only a performance effect. | Estimate free area from the occupancy grid. |
| Low | EKF | The filter assumes a gyro bias walk twice the truth. | Deliberate conservatism, but undocumented in the README. | One line in the README. |

## 4. Readability and self-documentation

```cpp
// before
const double Ld = std::max(6.0, 0.9 * in.v);
const int ti = std::min(n - 1, in.idx + static_cast<int>(Ld / 0.5));
...
if (k % 10 == 0) {  // 100 Hz IMU + wheel speed
  bias_true += 1e-4 * std::sqrt(0.01) * N01(rng);

// after
// Lookahead grows with speed (0.9 s of travel) but never below 6 m, where the 4 m wheelbase
// truck starts to weave. The 1 m floor on the chord only guards the division.
constexpr double kMinLookaheadM = 6.0, kLookaheadTimeS = 0.9, kMinChordM = 1.0;
const double Ld = std::max(kMinLookaheadM, kLookaheadTimeS * in.v);
const int ti = in.path->index_ahead(in.idx, Ld);
...
if (k % kImuEveryTicks == 0) {
  bias_true += kTrueGyroBiasWalk * std::sqrt(kImuPeriodS) * N01(rng);
```

## 5. Architecture and design

`astar` is templated on the graph and heuristic, so Dijkstra, A*, weighted A* and the Hybrid A* holonomic heuristic all share one tested implementation. JPS is written to search exactly the same graph (no corner cutting), which is what makes "same cost, fewer expansions" a meaningful comparison. D* Lite is checked against from-scratch A* after every edit, the strongest possible test for an incremental algorithm.

## 6. Testing and verification

| Claim | Verified by |
|---|---|
| A*, JPS, Dijkstra agree on cost | Tests and `results/search_compare.csv` |
| Weighted A* stays within its bound | Test |
| D* Lite equals replanning from scratch | Test after every map change |
| RRT* cost only improves and paths are collision-free | Test |
| Dubins paths reach the goal within the curvature limit, and a straight shot is exact | Tests (optimality among the six words is not tested) |
| Gated EKF rejects outliers and stays accurate | Test. Consistency (NEES, NIS) is reported in `results/` but not asserted |
| This review changed no behavior | All experiments rerun, every deterministic column identical |
| Memory safety | All 15 tests under ASan and UBSan |
| Any of it on a vehicle | Not verified |

## 7. Prioritized action plan

Must fix before commit: none outstanding.

Before sharing: split `run_scenario`, pin the pure-pursuit and Stanley metrics, and assert EKF consistency (mean NEES within its chi-squared interval).

Recommended: a global fallback in `RefPath::closest`.

Optional: free-area estimate for the RRT* constant.

## 8. Final verdict

| Context | Verdict |
|---|---|
| Personal or portfolio project | STRONG |
| Shared codebase | ACCEPTABLE, after the scenario split and regression tests |
| Simulation | STRONG |
| Vehicle | Not applicable. The kinematic truck and the scenario are not a vehicle validation. |
