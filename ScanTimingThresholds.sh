#!/bin/bash
# ScanTimingThresholds.sh
# Usage: ./ScanTimingThresholds.sh <time1> [<time2> ...]
# Example: ./ScanTimingThresholds.sh 0 2 3 5

set -e

if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <time_ns> [<time_ns> ...]"
    exit 1
fi

for TIME_NS in "$@"; do
    echo "========================================"
    echo "  Timing threshold: ${TIME_NS} ns"
    echo "========================================"
    # cd ..
    # ./scan.sh
    # cd PF-Reco-Analysis/
    ./SetTimingThreshold.sh "$TIME_NS"
    ./TestAllParticleFlow_DiPion.sh
    ./NtupleAllParticleFlow_dipi.sh
    ./PlotAllParticleFlow_dipi.sh
    python3 Plotting/dipi_2dheatmap.py \
        -d /eos/user/c/chtong/Public/Rereco/DiPionGun_20GeV_heatmap_2ns \
        -i "hcal_comparison_DR*_DT*.root" \
        -o /eos/user/c/chtong/Public/Rereco/DiPionGun_20GeV_heatmap_2ns
    #python3 Plotting/plot_pfrh_hcal_timealg.py
    echo "Finished timing threshold ${TIME_NS} ns"
done
