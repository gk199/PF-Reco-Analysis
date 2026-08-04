#!/bin/bash
# Condor wrapper for Z+jets PU re-reco jobs (step 4 of ZJets_PU_Resolution_Study_TODO.md).
# Runs HLT:Fake2,RAW2DIGI,L1Reco,RECO on one GEN-SIM-RAW input file, producing
# AODSIM. Which PF variant (nominal / timing3ns) this produces is determined
# entirely by which clusterizer build is active in the CMSSW area when this
# job runs — this script is identical for both, only the submit-time build
# state and output directory differ (same pattern as run_phaseScan_job.sh).
#
# Arguments:
#   $1 : input GEN-SIM-RAW file path (xrootd URL, e.g. from Wisconsin storage)
#   $2 : job index (integer, used to name the output file uniquely)
#   $3 : output directory (EOS path the output .root is copied to after the job)
#   $4 : basename of the grid proxy file (submit script copies it to a
#        persistent, non-/tmp location and lists it in transfer_input_files,
#        so it lands here in the job's own sandbox under this basename —
#        same pattern as gk199/Run3-HCAL-LLP-Analysis. /tmp proxies aren't
#        reliably present by the time the schedd actually transfers them.)

INPUT_FILE=$1
JOB_INDEX=$2
OUTPUT_DIR=$3
PROXY_BASENAME=$4

if [ -z "$INPUT_FILE" ] || [ -z "$JOB_INDEX" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$PROXY_BASENAME" ]; then
    echo "ERROR: Usage: run_zjetsReco_job.sh <inputFile> <jobIndex> <outputDir> <proxyBasename>"
    exit 1
fi

# Capture the job's starting sandbox dir before any `cd` below — the
# transferred proxy lands here.
export X509_USER_PROXY="$(pwd)/${PROXY_BASENAME}"
echo "=== Proxy check ==="
voms-proxy-info -all -file "${X509_USER_PROXY}"

OUTPUT_BASENAME="ZJetsPU_AOD_job${JOB_INDEX}.root"
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
# Stage the input file to local scratch first (the phase-scan study found
# streaming directly over xrootd gives wide, unpredictable per-job runtime
# variance — a single upfront copy plus local reads is far more consistent,
# and this input lives at T2_US_Wisconsin, not CERN-local, so the same
# WAN-read concern applies here too).
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
cmsRun MyPFStudy_ZJetsPU_RECO_condor.py \
    inputFiles="file:${LOCAL_INPUT}" \
    outputFile="${LOCAL_OUTPUT}" \
    maxEvents=-1

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: cmsRun failed with exit code ${EXIT_CODE}"
    rm -f "${LOCAL_INPUT}"
    exit ${EXIT_CODE}
fi

rm -f "${LOCAL_INPUT}"

# -----------------------------------------------------------------------
# Locate the actual output file. Not expected to trigger this time (we're
# not truncating early with maxEvents like the phase-scan study did, so
# PoolOutputModule shouldn't need the "_numEventN" split-renaming) — but
# checking costs nothing and the phase-scan study got burned once already
# by assuming the requested filename would always match reality.
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
# Copy output from local scratch to final destination. -f (force) so a
# stale/partial file from a previous failed attempt doesn't block this.
# -----------------------------------------------------------------------
echo "Copying ${OUTPUT_BASENAME} to ${OUTPUT_DIR}/"
xrdcp -f "${LOCAL_OUTPUT}" "${FINAL_OUTPUT}"

if [ $? -ne 0 ]; then
    echo "ERROR: Failed to copy output file"
    exit 1
fi

rm -f "${LOCAL_OUTPUT}"

echo "=== Job ${JOB_INDEX} finished successfully ==="
echo "Date: $(date)"
