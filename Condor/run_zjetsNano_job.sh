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
LOCAL_INPUT="${_CONDOR_SCRATCH_DIR:-/tmp}/input_job${JOB_INDEX}.root"
LOCAL_OUTPUT="${_CONDOR_SCRATCH_DIR:-/tmp}/${OUTPUT_BASENAME}"
FINAL_OUTPUT="${OUTPUT_DIR}/${OUTPUT_BASENAME}"

# Retry settings for both EOS transfers. 3 attempts with a growing backoff
# covers the failure actually seen: eosuser refused transfers for a ~35 s
# window while all 50 jobs hit it at once, first on the way out (xrdcp
# "Local error: protocol error (destination)") and, on the resubmission,
# on the way in (cmsRun exit 84, "[3014] ... Network is unreachable").
XRDCP_ATTEMPTS=3
XRDCP_BACKOFF=60          # seconds, multiplied by the attempt number
STAGGER_MAX=120           # start delay is spread over this many seconds

echo "=== Job ${JOB_INDEX} starting ==="
echo "Input:      ${INPUT_FILE}"
echo "Output:     ${FINAL_OUTPUT}"
echo "Hostname:   $(hostname)"
echo "Date:       $(date)"

# -----------------------------------------------------------------------
# Stagger the start. Condor releases the whole cluster at once — all 50 of
# these jobs began within 27 seconds and hammered the same EOS instance,
# which is what triggered the failures above. Spreading the first EOS
# access over a couple of minutes costs nothing here (these jobs are short
# and run on longlunch) and takes the concurrent load down by ~an order of
# magnitude.
# -----------------------------------------------------------------------
# Spread by job index rather than at random: 50 jobs over 120 s puts two per
# 5 s slot, where RANDOM would leave gaps and clumps. The jitter keeps the
# two jobs sharing a slot from starting in lockstep.
STAGGER=$(((JOB_INDEX * 5 + RANDOM % 5) % STAGGER_MAX))
echo "Staggering start by ${STAGGER} s to spread the load on EOS"
sleep ${STAGGER}

# -----------------------------------------------------------------------
# xrdcp with retries. Returns non-zero only if every attempt failed.
# -----------------------------------------------------------------------
xrdcp_retry() {
    local src=$1
    local dst=$2
    local attempt
    for attempt in $(seq 1 ${XRDCP_ATTEMPTS}); do
        if xrdcp -f "${src}" "${dst}"; then
            [ ${attempt} -gt 1 ] && echo "xrdcp succeeded on attempt ${attempt}"
            return 0
        fi
        echo "WARNING: xrdcp attempt ${attempt}/${XRDCP_ATTEMPTS} failed: ${src} -> ${dst}"
        if [ ${attempt} -lt ${XRDCP_ATTEMPTS} ]; then
            local wait=$((XRDCP_BACKOFF * attempt))
            echo "Retrying in ${wait} s..."
            sleep ${wait}
        fi
    done
    return 1
}

# -----------------------------------------------------------------------
# Set up CMSSW environment
# -----------------------------------------------------------------------
source /cvmfs/cms.cern.ch/cmsset_default.sh

CMSSW_BASE_DIR=/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6
cd ${CMSSW_BASE_DIR}/src
eval `scramv1 runtime -sh`   # cmsenv

cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis

# -----------------------------------------------------------------------
# Stage the input to local scratch, same as run_zjetsReco_job.sh. This used
# to read straight over xrootd on the grounds that the input is CERN-local,
# which held until eosuser refused 39 of 50 opens in one submission and
# cmsRun died with exit 84 — an open that fails mid-job is unrecoverable,
# whereas a staged copy can simply be retried. It also drops the xrootd
# connection cmsRun would otherwise hold open for its whole runtime.
# -----------------------------------------------------------------------
echo "Staging input to local scratch..."
if ! xrdcp_retry "${INPUT_FILE}" "${LOCAL_INPUT}"; then
    echo "ERROR: Failed to stage input file after ${XRDCP_ATTEMPTS} attempts"
    rm -f "${LOCAL_INPUT}"
    exit 1
fi

cmsRun MyPFStudy_ZJetsPU_NANO_condor.py \
    inputFiles="file:${LOCAL_INPUT}" \
    outputFile="${LOCAL_OUTPUT}" \
    maxEvents=-1

EXIT_CODE=$?

rm -f "${LOCAL_INPUT}"

if [ ${EXIT_CODE} -ne 0 ]; then
    echo "ERROR: cmsRun failed with exit code ${EXIT_CODE}"
    exit ${EXIT_CODE}
fi

# -----------------------------------------------------------------------
# Copy output from local scratch to final destination. -f (force) so a
# stale/partial file from a previous failed attempt doesn't block this.
#
# On permanent failure the partial destination file is removed: a failed
# xrdcp leaves one behind, and it is worse than no file at all. The last
# round left both 0-byte stubs and one 2.8 MB file that looked the right
# size but was unreadable (uproot: "Input/output error"), so a size check
# alone would not have caught it — nothing downstream should ever see it.
# -----------------------------------------------------------------------
echo "Copying ${OUTPUT_BASENAME} to ${OUTPUT_DIR}/"
if ! xrdcp_retry "${LOCAL_OUTPUT}" "${FINAL_OUTPUT}"; then
    echo "ERROR: Failed to copy output file after ${XRDCP_ATTEMPTS} attempts"
    echo "Removing partial destination file ${FINAL_OUTPUT}"
    # OUTPUT_DIR is a plain /eos path, so try the xrootd removal first and
    # fall back to a plain rm for a node with EOS mounted (or a non-EOS dest).
    xrdfs root://eosuser.cern.ch rm "${FINAL_OUTPUT}" 2>/dev/null || rm -f "${FINAL_OUTPUT}"
    exit 1
fi

rm -f "${LOCAL_OUTPUT}"

echo "=== Job ${JOB_INDEX} finished successfully ==="
echo "Date: $(date)"
