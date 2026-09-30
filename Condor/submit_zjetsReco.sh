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

# Always strip comments/blank lines — a bug here previously submitted raw
# comment lines as literal (bogus) inputFile values, which also shifted
# every real job's $(ProcId) away from its input file's own job index.
SUBMIT_LIST=$(mktemp)
grep -v '^\s*#' "${INPUT_LIST}" | grep -v '^\s*$' > "${SUBMIT_LIST}"
if [ ${TEST_MODE} -eq 1 ]; then
    head -1 "${SUBMIT_LIST}" > "${SUBMIT_LIST}.tmp"
    mv "${SUBMIT_LIST}.tmp" "${SUBMIT_LIST}"
fi

N_JOBS=$(wc -l < "${SUBMIT_LIST}")

echo "Input file list : ${INPUT_LIST}"
echo "Output directory: ${OUTPUT_DIR}"
echo "Logs directory  : ${LOG_DIR}"
[ ${TEST_MODE} -eq 1 ] && echo "TEST MODE: submitting first job only"
echo "Jobs to submit  : ${N_JOBS}"
echo ""

# -----------------------------------------------------------------------
# Proxy: the input GEN-SIM-RAW files live at T2_US_Wisconsin, a genuine
# remote WLCG site (unlike the phase-scan study's CERN-internal EOS reads)
# — xrdcp there needs a real grid proxy. use_x509userproxy=True failed
# here in practice ("Transfer input files failure ... reading from file
# /tmp/x509up_u<uid>: No such file or directory") — /tmp proxies aren't
# reliably present by the time the schedd actually transfers them (can be
# a while after submission if the job sits idle). Fix, matching
# gk199/Run3-HCAL-LLP-Analysis: copy the current proxy to a persistent,
# non-/tmp location and explicitly ship it via transfer_input_files.
# -----------------------------------------------------------------------
echo "Checking for a valid grid proxy..."
CURRENT_PROXY=$(voms-proxy-info --path 2>/dev/null)
if [ -z "${CURRENT_PROXY}" ] || ! voms-proxy-info --exists 2>/dev/null; then
    echo "ERROR: no valid grid proxy found. Run:"
    echo "  voms-proxy-init --rfc --voms cms --valid 172:00"
    echo "before submitting."
    exit 1
fi

PERSISTENT_PROXY_DIR="${HOME}/private"
mkdir -p "${PERSISTENT_PROXY_DIR}"
chmod 700 "${PERSISTENT_PROXY_DIR}"
PROXY_BASENAME="$(basename "${CURRENT_PROXY}")"
PERSISTENT_PROXY="${PERSISTENT_PROXY_DIR}/${PROXY_BASENAME}"
cp -f "${CURRENT_PROXY}" "${PERSISTENT_PROXY}"
chmod 600 "${PERSISTENT_PROXY}"
echo "Proxy copied to persistent location: ${PERSISTENT_PROXY}"

condor_submit -terse - << EOF
universe              = vanilla
executable            = ${WRAPPER}
arguments             = \$(inputFile) \$(ProcId) ${OUTPUT_DIR} ${PROXY_BASENAME}
output                = ${LOG_DIR}/zjetsReco_job_\$(ProcId).out
error                 = ${LOG_DIR}/zjetsReco_job_\$(ProcId).err
log                   = ${LOG_DIR}/zjetsReco_condor.log
transfer_input_files  = ${PERSISTENT_PROXY}
+JobFlavour           = "tomorrow"
request_memory        = 4000
request_cpus          = 1
request_disk          = 50000000
should_transfer_files = YES
transfer_output_files = ""
queue inputFile from ${SUBMIT_LIST}
EOF

rm -f "${SUBMIT_LIST}"

echo ""
echo "Submitted ${N_JOBS} jobs."
echo "Monitor with:  condor_q"
echo "Cancel all:    condor_rm <clusterID>"
echo "Logs in:       ${LOG_DIR}/"
