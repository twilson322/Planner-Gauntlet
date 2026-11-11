#!/usr/bin/env bash
# Build, test, run every experiment and benchmark, then plot. About 1 minute on 2 cores.
set -euo pipefail
cd "$(dirname "$0")/.."
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
(cd build && ctest --output-on-failure)
mkdir -p results
./build/pg_experiments all --out results | tee results/experiments.txt
./build/pg_bench | tee results/benchmark.txt
python3 scripts/plot.py --results results
