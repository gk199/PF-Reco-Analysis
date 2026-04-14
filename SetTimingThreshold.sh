#!/bin/bash
# SetTimingThreshold.sh
# Usage: ./SetTimingThreshold.sh <time_ns>
#
# Updates the timeThreshold parameter in all timing cfi.py and .cc.edit files
# under PFTestingAlgos/.
#
# In the .cc.edit files the value is a C++ fallback default; the cfi.py value
# takes precedence at runtime.  Both are updated here for consistency.
#
# Run this before TestAllParticleFlow.sh (or before manually copying files)
# so the desired threshold is baked into every variant.
#
# Example: ./SetTimingThreshold.sh 3
#          ./SetTimingThreshold.sh 5.0

if [ -z "$1" ]; then
    echo "Usage: $0 <time_ns>"
    echo "  e.g. $0 3    # sets timeThreshold to 3 ns"
    echo "  e.g. $0 5.0  # sets timeThreshold to 5.0 ns"
    exit 1
fi

TIME_NS="$1"

# Validate: must be a number (integer or decimal, optionally negative)
if ! [[ "$TIME_NS" =~ ^-?[0-9]+(\.[0-9]+)?$ ]]; then
    echo "Error: '$TIME_NS' is not a valid number."
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# cfi.py files: timeThreshold = cms.untracked.double(<value>)
CFI_FILES=(
    "PFTestingAlgos/particleFlowClusterHCAL_timing_cfi.py"
    "PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py"
    "PFTestingAlgos/particleFlowClusterHCAL_depth1Timing_cfi.py"
)

# .cc.edit files: getUntrackedParameter<double>("timeThreshold", <default>)
# The cfi.py value takes precedence; these defaults are updated for consistency.
CC_FILES=(
    "PFTestingAlgos/Basic2DGenericTopoClusterizer_timing.cc.edit"
    "PFTestingAlgos/PFMultiDepthClusterizer_timing.cc.edit"
    "PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit"
    "PFTestingAlgos/PFMultiDepthClusterizer_depth1Timing.cc.edit"
)

echo "Setting timeThreshold = ${TIME_NS} ns in all timing files..."
echo ""

echo "  --- cfi.py files ---"
for CFI in "${CFI_FILES[@]}"; do
    FULL_PATH="${SCRIPT_DIR}/${CFI}"
    if [ ! -f "$FULL_PATH" ]; then
        echo "  WARNING: file not found, skipping: $CFI"
        continue
    fi
    sed -i "s|timeThreshold = cms.untracked.double([^)]*)|timeThreshold = cms.untracked.double(${TIME_NS})|g" "$FULL_PATH"
    echo "  Updated: $CFI"
done

echo ""
echo "  --- .cc.edit files (C++ fallback default) ---"
for CC in "${CC_FILES[@]}"; do
    FULL_PATH="${SCRIPT_DIR}/${CC}"
    if [ ! -f "$FULL_PATH" ]; then
        echo "  WARNING: file not found, skipping: $CC"
        continue
    fi
    sed -i 's|getUntrackedParameter<double>("timeThreshold", [^)]*)|getUntrackedParameter<double>("timeThreshold", '"${TIME_NS}"')|g' "$FULL_PATH"
    echo "  Updated: $CC"
done

echo ""
echo "Done. timeThreshold = ${TIME_NS} ns is now set in all timing files."
echo "Re-run TestAllParticleFlow.sh (or copy the desired files manually + scram b) to apply."
