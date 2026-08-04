#!/bin/bash
# Condor wrapper for Z+jets PU NanoAOD jobs (step 5 of ZJets_PU_Resolution_Study_TODO.md).
# Runs PAT,NANO on one AODSIM input file (step 4's output), producing NANOAODSIM.
# Same variant convention as run_zjetsReco_job.sh: this script is identical
# for nominal/timing3ns — point it at the right AOD source directory and
# output directory per variant, the PF variant distinction was already
# baked in upstream at the RECO step.
#
# Arguments:
#   $1 : input AODSIM file path (from step 4's output, e.g. on EOS)
#   $2 : job index (integer, used to name the output file uniquely)
#   $3 : output directory (EOS path the output .root is copied to after the job)
#   $4 : basename of the grid proxy file (see run_zjetsReco_job.sh for why —
#        same fix applied here for consistency even though this input is
#        CERN-local, not the Wisconsin/WAN case that actually needed it)

INPUT_FILE=$1
JOB_INDEX=$2
OUTPUT_DIR=$3
PROXY_BASENAME=$4

if [ -z "$INPUT_FILE" ] || [ -z "$JOB_INDEX" ] || [ -z "$OUTPUT_DIR" ] || [ -z "$PROXY_BASENAME" ]; then
    echo "ERROR: Usage: run_zjetsNano_job.sh <inputFile> <jobIndex> <outputDir> <proxyBasename>"
    exit 1
fi

export X509_USER_PROXY="$(pwd)/${PROXY_BASENAME}"
echo "=== Proxy check ==="
voms-proxy-info -all -file "${X509_USER_PROXY}"

OUTPUT_BASENAME="ZJetsPU_NANO_job${JOB_INDEX}.root"
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
# No local input staging here — unlike step 4 (Wisconsin/WAN), this input
# is step 4's own AODSIM output on our EOS (CERN-local), same situation as
# the phase-scan ntupling step which read directly over xrootd successfully.
# -----------------------------------------------------------------------
cmsRun MyPFStudy_ZJetsPU_NANO_condor.py \
    inputFiles="${INPUT_FILE}" \
    outputFile="${LOCAL_OUTPUT}" \
    maxEvents=-1

EXIT_CODE=$?

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: cmsRun failed with exit code ${EXIT_CODE}"
    exit ${EXIT_CODE}
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
