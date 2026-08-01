#!/bin/bash
# HTCondor submission script for Z+jets PU re-reco jobs (step 4 of
# ZJets_PU_Resolution_Study_TODO.md). One job per input GEN-SIM-RAW file.
#
# Which PF variant this produces (nominal / timing3ns) is determined by
# which clusterizer build is active in the CMSSW area at submission time —
# submit once per variant, rebuilding in between, pointing each run at a
# different output_dir (same pattern as the phase-scan nominal_fresh/ vs
# timing_3ns/ directories).
#
# Usage:
#   bash submit_zjetsReco.sh [-t] <input_files.txt> <output_dir>
#
# Options:
#   -t   Test mode: submit only the first job in the input list
#
# IMPORTANT: memory and runtime requirements below are unverified guesses —
# RECO (tracking/vertexing/PF) on high-PU events may be as memory/CPU-heavy
# as generation was (which needed real tuning: 8000MB, 4 cores, 480min
# before it stopped failing). Run with -t first, check the actual memory/
# runtime the test job used (condor_q / job logs), and adjust
# request_memory / +JobFlavour below before submitting the full batch —
# don't assume these defaults are right.

TEST_MODE=0
if [[ "$1" == "-t" ]]; then
    TEST_MODE=1
    shift
fi

INPUT_LIST=$1
OUTPUT_DIR=$2

if [ -z "${INPUT_LIST}" ] || [ -z "${OUTPUT_DIR}" ]; then
    echo "Usage: $0 [-t] <input_files.txt> <output_dir>"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER="${SCRIPT_DIR}/run_zjetsReco_job.sh"
LOG_DIR="${SCRIPT_DIR}/logs"

if [ ! -f "${INPUT_LIST}" ]; then
    echo "ERROR: Input file list not found: ${INPUT_LIST}"
    exit 1
fi

if [ ! -f "${WRAPPER}" ]; then
    echo "ERROR: Wrapper script not found: ${WRAPPER}"
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

SUBMIT_LIST="${INPUT_LIST}"
if [ ${TEST_MODE} -eq 1 ]; then
    SUBMIT_LIST=$(mktemp)
    grep -v '^\s*#' "${INPUT_LIST}" | grep -v '^\s*$' | head -1 > "${SUBMIT_LIST}"
fi

N_JOBS=$(grep -v '^\s*#' "${SUBMIT_LIST}" | grep -v '^\s*$' | wc -l)

echo "Input file list : ${INPUT_LIST}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Logs directory  : ${LOG_DIR}"
[ ${TEST_MODE} -eq 1 ] && echo "TEST MODE: submitting first job only"
echo "Jobs to submit  : ${N_JOBS}"
echo ""

# -----------------------------------------------------------------------
# Submit all jobs as a single cluster: one proc ID per input file.
# $(ProcId) and $(inputFile) are filled by condor for each queued item.
# -----------------------------------------------------------------------
condor_submit -terse - << EOF
universe              = vanilla
executable            = ${WRAPPER}
arguments             = \$(inputFile) \$(ProcId) ${OUTPUT_DIR}
output                = ${LOG_DIR}/zjetsReco_job_\$(ProcId).out
error                 = ${LOG_DIR}/zjetsReco_job_\$(ProcId).err
log                   = ${LOG_DIR}/zjetsReco_condor.log
+JobFlavour           = "tomorrow"
request_memory        = 4000
request_cpus          = 1
request_disk          = 50000000
should_transfer_files = YES
transfer_output_files = ""
queue inputFile from ${SUBMIT_LIST}
EOF

# Clean up temp file used in test mode
[ ${TEST_MODE} -eq 1 ] && rm -f "${SUBMIT_LIST}"

echo ""
echo "Submitted ${N_JOBS} jobs."
echo "Monitor with:  condor_q"
echo "Cancel all:    condor_rm <clusterID>"
echo "Logs in:       ${LOG_DIR}/"
