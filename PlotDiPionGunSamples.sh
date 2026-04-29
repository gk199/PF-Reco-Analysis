#!/bin/bash
# PlotDiPionGunSamples.sh
# Run compare_hcal_clusters.py for DiPionGun ntuples across timing thresholds.
#
# Usage: ./PlotDiPionGunSamples.sh [<time_ns> ...]
# Example: ./PlotDiPionGunSamples.sh 0 3 5
#
# Defaults to 0 ns if no arguments are given.
#
# Expects ntuples named:
#   pfObjectsNtuple_DiPionGun_DR<X>_DT<Y>_<algo>_<T>ns.root
# Output per DR/DT/timing combination:
#   hcal_comparison_DiPionGun_DR<X>_DT<Y>_<T>ns.{pdf,root}

TIMES="${*:-0}"
PF_DIR="$(cd "$(dirname "$0")" && pwd)"
PLOT_SCRIPT="$PF_DIR/Plotting/compare_hcal_clusters.py"

# DRS="0.1 0.4 0.8"
# DTS="0.0 1.0 5.0"
DRS="0.4"
DTS="0.0"

# ── main loop ─────────────────────────────────────────────────────────────────

for TIME_NS in $TIMES; do
    echo "========================================"
    echo "  Timing threshold: ${TIME_NS} ns"
    echo "========================================"

    for DR in $DRS; do
        for DT in $DTS; do
            LABEL="DiPionGun_DR${DR}_DT${DT}"
            PREFIX="pfObjectsNtuple_${LABEL}_"
            SUFFIX="_${TIME_NS}ns"
            OUT_ROOT="${PF_DIR}/hcal_comparison_${LABEL}_${TIME_NS}ns.root"
            OUT_PDF="${PF_DIR}/hcal_comparison_${LABEL}_${TIME_NS}ns.pdf"

            echo "-----------------------------------"
            echo "  Sample: ${LABEL}  |  threshold: ${TIME_NS} ns"
            echo "-----------------------------------"

            python3 "$PLOT_SCRIPT" \
                --inputdir "$PF_DIR" \
                --prefix   "$PREFIX" \
                --suffix   "$SUFFIX" \
                --output   "$OUT_ROOT" \
                --pdf      "$OUT_PDF"

            if [ $? -ne 0 ]; then
                echo "WARNING: plotting failed for ${LABEL} at ${TIME_NS}ns"
            else
                echo "Saved: ${OUT_PDF}"
            fi
        done
    done
done

echo "========================================"
echo "  Done."
echo "========================================"
