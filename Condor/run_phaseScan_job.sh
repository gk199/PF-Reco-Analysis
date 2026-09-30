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
LOCAL_INPUT="${_CONDOR_SCRATCH_DIR:-/tmp}/input_job${JOB_INDEX}.root"

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
# Stage the input file to local scratch first. Streaming events directly
# over xrootd from the shared EOS redirector showed wide, unpredictable
# per-job runtime variance (25 min to 24h+ for the same maxEvents cap) —
# a single upfront copy plus local reads is far more consistent.
# -----------------------------------------------------------------------
echo "Staging input to local scratch..."
xrdcp "${INPUT_FILE}" "${LOCAL_INPUT}"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to stage input file"
    exit 1
fi

# -----------------------------------------------------------------------
# Run cmsRun
# -----------------------------------------------------------------------
cmsRun MyPFStudy_ReReco_RAW2DIGI_L1Reco_RECO_phaseScan_condor.py \
    inputFiles="file:${LOCAL_INPUT}" \
    outputFile="${LOCAL_OUTPUT}" \
    maxEvents=4000

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: cmsRun failed with exit code ${EXIT_CODE}"
    rm -f "${LOCAL_INPUT}"
    exit ${EXIT_CODE}
fi

rm -f "${LOCAL_INPUT}"

# -----------------------------------------------------------------------
# Locate the actual output file. PoolOutputModule appends "_numEventN" to
# the requested filename when maxEvents truncates processing before a
# run/lumi boundary is naturally reached (splits the file to keep run
# bookkeeping self-consistent) — with a finite maxEvents cap the file on
# disk may not exactly match LOCAL_OUTPUT, so find it rather than assume.
# -----------------------------------------------------------------------
ACTUAL_OUTPUT=$(ls "${LOCAL_OUTPUT%.root}"*.root 2>/dev/null | head -1)
if [ -z "${ACTUAL_OUTPUT}" ]; then
    echo "ERROR: no output file found matching ${LOCAL_OUTPUT%.root}*.root"
    exit 1
fi
if [ "${ACTUAL_OUTPUT}" != "${LOCAL_OUTPUT}" ]; then
    echo "Note: cmsRun wrote ${ACTUAL_OUTPUT} (renaming to ${OUTPUT_BASENAME} for upload)"
    mv "${ACTUAL_OUTPUT}" "${LOCAL_OUTPUT}"
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
