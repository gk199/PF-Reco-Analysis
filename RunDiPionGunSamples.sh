#!/bin/bash
# RunDiPionGunSamples.sh
# 3D scan: timing threshold x PF algorithm x DiPionGun DR/DT sample.
#
# Usage: ./RunDiPionGunSamples.sh <time_ns> [<time_ns> ...]
# Example: ./RunDiPionGunSamples.sh 3 5
#
# For each timing value and each PF algorithm, compiles once, then runs
# cmsRun + ntupler for every DR/DT input file.
#
# Outputs: pfObjectsNtuple_DiPionGun_DR<X>_DT<Y>_<algo>_<T>ns.root

if [ -z "$1" ]; then
    echo "Usage: $0 <time_ns> [<time_ns> ...]"
    echo "  e.g. $0 3 5"
    exit 1
fi

MC_DIR="/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src"
PF_DIR="$(cd "$(dirname "$0")" && pwd)"
CMSSW_DIR="$(dirname "$PF_DIR")"

cd "$CMSSW_DIR"
cmsenv
cd "$PF_DIR"

# DRS="0.1 0.4 0.8"
# DTS="0.0 1.0 5.0"
DRS="0.1 0.4 0.8"
DTS="0.0 1.0 5.0"

ALGOS="standardPF cellTimingPF seedTimingPF depth1SeedTimingPF"

# Helper: copy the right plugin files for a given algorithm
install_algo() {
    local ALGO="$1"
    cd "$CMSSW_DIR"
    case "$ALGO" in
        standardPF)
            cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit  RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit     RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit         RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py          RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
            ;;
        cellTimingPF)
            cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_timing.cc.edit     RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit       RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_timing.cc.edit           RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_timing_cfi.py            RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
            ;;
        seedTimingPF)
            cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit   RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit       RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit       RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py        RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
            ;;
        depth1SeedTimingPF)
            cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit   RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit       RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_depth1Timing.cc.edit     RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
            cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_depth1Timing_cfi.py      RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
            ;;
    esac
    scram b -j 8
    cd "$PF_DIR"
}

# ── main loop ─────────────────────────────────────────────────────────────────

for TIME_NS in "$@"; do
    echo "========================================"
    echo "  Timing threshold: ${TIME_NS} ns"
    echo "========================================"

    ./SetTimingThreshold.sh "${TIME_NS}"

    for ALGO in $ALGOS; do
        echo "-----------------------------------"
        echo "  Algorithm: ${ALGO}  |  threshold: ${TIME_NS} ns"
        echo "-----------------------------------"

        install_algo "$ALGO"

        for DR in $DRS; do
            for DT in $DTS; do
                LABEL="DiPionGun_DR${DR}_DT${DT}"
                INPUT="${MC_DIR}/${LABEL}_GEN-SIM-RAW.root"
                RECO_OUT="pf_only_reReco_${LABEL}_${ALGO}_${TIME_NS}ns.root"
                NTUPLE_OUT="pfObjectsNtuple_${LABEL}_${ALGO}_${TIME_NS}ns.root"

                if [ ! -f "$INPUT" ]; then
                    echo "WARNING: input not found, skipping: $INPUT"
                    continue
                fi

                cmsRun MyPFStudy_ReReco_DiPionGun_DIGI_RAW2DIGI_L1Reco_RECO.py \
                    inputFiles=file:${INPUT} \
                    outputFile=${RECO_OUT}
                if [ $? -ne 0 ] || [ ! -f "${RECO_OUT}" ]; then
                    echo "ERROR: RECO failed for ${LABEL} ${ALGO} ${TIME_NS}ns — skipping ntupler"
                    continue
                fi

                cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
                    inputFiles=file:${RECO_OUT} \
                    outputFile=${NTUPLE_OUT}
            done
        done
    done
done

# Revert CMSSW to standard PF
cd "$CMSSW_DIR"
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit  RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit     RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit         RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py          RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
scram b -j 8
cd "$PF_DIR"

echo "========================================"
echo "  Done. CMSSW reverted to standard PF."
echo "========================================"
