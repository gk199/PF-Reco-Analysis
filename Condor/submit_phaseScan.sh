#!/bin/bash
# HTCondor submission script for PF phase scan re-reco jobs.
# All jobs are submitted as a single cluster (one proc ID per file).
#
# Usage:
#   bash submit_phaseScan.sh [-t] [input_files.txt] [output_dir]
#
# Options:
#   -t   Test mode: submit only the first job in the input list
#
# Defaults:
#   input_files.txt  ->  ./input_files.txt
#   output_dir       ->  /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/condor_output/phaseScan

TEST_MODE=0
if [[ "$1" == "-t" ]]; then
    TEST_MODE=1
    shift
fi

INPUT_LIST=${1:-input_files.txt}
OUTPUT_DIR=${2:-/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/condor_output/phaseScan}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER="${SCRIPT_DIR}/run_phaseScan_job.sh"
LOG_DIR="${SCRIPT_DIR}/logs"

# -----------------------------------------------------------------------
# Sanity checks
# -----------------------------------------------------------------------
if [ ! -f "${INPUT_LIST}" ]; then
    echo "ERROR: Input file list not found: ${INPUT_LIST}"
    exit 1
fi

if [ ! -f "${WRAPPER}" ]; then
    echo "ERROR: Wrapper script not found: ${WRAPPER}"
    exit 1
fi

# Create output and log directories if they don't exist
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"

# In test mode, build a temporary single-entry file list
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
output                = ${LOG_DIR}/job_\$(ProcId).out
error                 = ${LOG_DIR}/job_\$(ProcId).err
log                   = ${LOG_DIR}/condor.log
+JobFlavour           = "workday"
request_memory        = 4000
request_cpus          = 1
should_transfer_files = NO
queue inputFile from ${SUBMIT_LIST}
EOF

# Clean up temp file used in test mode
[ ${TEST_MODE} -eq 1 ] && rm -f "${SUBMIT_LIST}"

echo ""
echo "Submitted ${N_JOBS} jobs."
echo "Monitor with:  condor_q"
echo "Cancel all:    condor_rm <clusterID>"
echo "Logs in:       ${LOG_DIR}/"
