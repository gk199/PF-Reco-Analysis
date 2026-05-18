#!/bin/bash
# Condor wrapper for HCAL cluster timing fix validation.
# Runs re-reco and ntupling in a single job; only the ntuple is kept.
#
# Arguments:
#   $1 : input ROOT file (e.g. root://cms-xrd-global.cern.ch//store/mc/...)
#   $2 : job index (integer, used to name output file uniquely)
#   $3 : output directory (EOS or AFS path for the ntuple)
#   $4 : CMSSW base directory (full path, e.g. .../CMSSW_15_0_6)
#   $5 : max events per job (default: 200)

INPUT_FILE=$1
JOB_INDEX=$2
OUTPUT_DIR=$3
CMSSW_BASE_DIR=${4:-/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6}
MAX_EVENTS=${5:-200}

if [ -z "$INPUT_FILE" ] || [ -z "$JOB_INDEX" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Usage: run_timingFix_job.sh <inputFile> <jobIndex> <outputDir> [cmssw_dir] [maxEvents]"
    exit 1
fi

SCRATCH=${_CONDOR_SCRATCH_DIR:-/tmp}
LOCAL_RECO="${SCRATCH}/pf_reco_timingFix_job${JOB_INDEX}.root"
LOCAL_NTUPLE="${SCRATCH}/pfObjectsNtuple_timingFix_job${JOB_INDEX}.root"
FINAL_NTUPLE="${OUTPUT_DIR}/pfObjectsNtuple_job${JOB_INDEX}.root"

echo "=== Job ${JOB_INDEX} starting ==="
echo "Input:      ${INPUT_FILE}"
echo "Output:     ${FINAL_NTUPLE}"
echo "CMSSW:      ${CMSSW_BASE_DIR}"
echo "MaxEvents:  ${MAX_EVENTS}"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"

# -----------------------------------------------------------------------
# CMSSW environment
# -----------------------------------------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd "${CMSSW_BASE_DIR}/src"
eval "$(scramv1 runtime -sh)"

cd "${CMSSW_BASE_DIR}/src/PF-Reco-Analysis"

# -----------------------------------------------------------------------
# Step 1: re-reco (RAW2DIGI -> L1Reco -> RECO)
# -----------------------------------------------------------------------
echo "--- Running re-reco ---"
cmsRun MyPFStudy_ReReco_MC_condor.py \
    inputFiles="${INPUT_FILE}" \
    outputFile="${LOCAL_RECO}" \
    maxEvents=${MAX_EVENTS}

if [ $? -ne 0 ]; then
    echo "ERROR: re-reco cmsRun failed"
    exit 1
fi

# -----------------------------------------------------------------------
# Step 2: ntupler
# -----------------------------------------------------------------------
echo "--- Running ntupler ---"
cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles="file:${LOCAL_RECO}" \
    outputFile="${LOCAL_NTUPLE}"

if [ $? -ne 0 ]; then
    echo "ERROR: ntupler cmsRun failed"
    exit 1
fi

# -----------------------------------------------------------------------
# Copy ntuple to final destination; discard the intermediate RECO file
# -----------------------------------------------------------------------
echo "--- Copying ntuple to ${OUTPUT_DIR}/ ---"
xrdcp "${LOCAL_NTUPLE}" "${FINAL_NTUPLE}"

if [ $? -ne 0 ]; then
    echo "ERROR: xrdcp failed"
    exit 1
fi

rm -f "${LOCAL_RECO}" "${LOCAL_NTUPLE}"

echo "=== Job ${JOB_INDEX} finished successfully ==="
echo "Date: $(date)"
