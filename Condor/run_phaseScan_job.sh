#!/bin/bash
# Condor wrapper for PF phase scan re-reco jobs.
# Arguments:
#   $1 : input ROOT file path (e.g. file:/eos/cms/... or root://cms-xrd-global.cern.ch//store/...)
#   $2 : job index (integer, used to name the output file uniquely)
#   $3 : output directory (AFS or EOS path where the output .root is copied after the job)

INPUT_FILE=$1
JOB_INDEX=$2
OUTPUT_DIR=$3

if [ -z "$INPUT_FILE" ] || [ -z "$JOB_INDEX" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "ERROR: Usage: run_phaseScan_job.sh <inputFile> <jobIndex> <outputDir>"
    exit 1
fi

OUTPUT_BASENAME="pf_only_reReco_phaseScan_job${JOB_INDEX}.root"
# Write to local scratch first — EOS/AFS don't support direct ROOT file creation
LOCAL_OUTPUT="${_CONDOR_SCRATCH_DIR:-/tmp}/${OUTPUT_BASENAME}"
FINAL_OUTPUT="${OUTPUT_DIR}/${OUTPUT_BASENAME}"

echo "=== Job ${JOB_INDEX} starting ==="
echo "Input:      ${INPUT_FILE}"
echo "Output:     ${FINAL_OUTPUT}"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"

# -----------------------------------------------------------------------
# Set up CMSSW environment
# -----------------------------------------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh

CMSSW_BASE_DIR=/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6
cd ${CMSSW_BASE_DIR}/src
eval `scramv1 runtime -sh`   # cmsenv

cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis

# -----------------------------------------------------------------------
# Run cmsRun
# -----------------------------------------------------------------------
cmsRun MyPFStudy_ReReco_RAW2DIGI_L1Reco_RECO_phaseScan_condor.py \
    inputFiles="${INPUT_FILE}" \
    outputFile="${LOCAL_OUTPUT}" \
    maxEvents=-1

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: cmsRun failed with exit code ${EXIT_CODE}"
    exit ${EXIT_CODE}
fi

# -----------------------------------------------------------------------
# Copy output from local scratch to final destination
# -----------------------------------------------------------------------
echo "Copying ${OUTPUT_BASENAME} to ${OUTPUT_DIR}/"
xrdcp "${LOCAL_OUTPUT}" "${FINAL_OUTPUT}"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to copy output file"
    exit 1
fi

rm -f "${LOCAL_OUTPUT}"

echo "=== Job ${JOB_INDEX} finished successfully ==="
echo "Date: $(date)"
