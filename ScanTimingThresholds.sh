#!/bin/bash
# ScanTimingThresholds.sh
# Usage: ./ScanTimingThresholds.sh <time1> [<time2> ...]
# Example: ./ScanTimingThresholds.sh 0 2 3 5

if [ -z "$1" ]; then
    echo "Usage: $0 <time_ns> [<time_ns> ...]"
    exit 1
fi

for TIME_NS in "$@"; do
    echo "========================================"
    echo "  Timing threshold: ${TIME_NS} ns"
    echo "========================================"

    ./SetTimingThreshold.sh ${TIME_NS}
    ./TestAllParticleFlow.sh
    ./NtupleAllParticleFlow.sh
    ./PlotAllParticleFlow.sh
    mv hcal_comparison.pdf hcal_comparison_${TIME_NS}ns.pdf
done
