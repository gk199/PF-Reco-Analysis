#!/bin/bash
# HTCondor submission for HCAL cluster timing fix validation on physics samples.
# Submits one job per input file; each job runs re-reco + ntupler and writes
# only the ntuple to the output directory.
#
# Usage:
#   bash submit_timingFix.sh [-t] <before|after> <input_files.txt> <output_dir> [cmssw_dir] [max_events]
#
# Options:
#   -t             Test mode: submit only the first job
#   before|after   Label appended to output filenames and used in log directory
#   input_files.txt  One xrootd or file: path per line
#   output_dir     EOS or AFS destination for ntuple files
#   cmssw_dir      CMSSW_BASE directory to source (default: the area containing this script)
#   max_events     Events per job (default: 200)
#
# Typical workflow:
#   # Produce "before" ntuples (revert Basic2DGenericPFlowPositionCalc.cc, scram b)
#   bash submit_timingFix.sh before input_files_ttbar.txt \
#       /eos/user/g/gkopp/timingFix/before /path/to/CMSSW_15_0_6_before
#
#   # Produce "after" ntuples (with timing fix compiled in)
#   bash submit_timingFix.sh after input_files_ttbar.txt \
#       /eos/user/g/gkopp/timingFix/after /path/to/CMSSW_15_0_6

TEST_MODE=0
if [[ "$1" == "-t" ]]; then
    TEST_MODE=1
    shift
fi

VARIANT=${1}          # "before" or "after"
INPUT_LIST=${2}
OUTPUT_DIR=${3}
CMSSW_BASE_DIR=${4:-/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6}
MAX_EVENTS=${5:-200}

if [ -z "$VARIANT" ] || [ -z "$INPUT_LIST" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Usage: submit_timingFix.sh [-t] <before|after> <input_files.txt> <output_dir> [cmssw_dir] [max_events]"
    exit 1
fi

if [[ "$VARIANT" != "before" && "$VARIANT" != "after" ]]; then
    echo "ERROR: variant must be 'before' or 'after', got: $VARIANT"
    exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WRAPPER="${SCRIPT_DIR}/run_timingFix_job.sh"
LOG_DIR="${SCRIPT_DIR}/logs/timingFix_${VARIANT}"

if [ ! -f "${INPUT_LIST}" ]; then
    echo "ERROR: Input file list not found: ${INPUT_LIST}"
    exit 1
fi
if [ ! -f "${WRAPPER}" ]; then
    echo "ERROR: Wrapper not found: ${WRAPPER}"
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

echo "Variant        : ${VARIANT}"
echo "Input list     : ${INPUT_LIST}"
echo "Output dir     : ${OUTPUT_DIR}"
echo "CMSSW dir      : ${CMSSW_BASE_DIR}"
echo "Max events/job : ${MAX_EVENTS}"
echo "Jobs to submit : ${N_JOBS}"
[ ${TEST_MODE} -eq 1 ] && echo "TEST MODE: first job only"
echo ""

condor_submit -terse - << EOF
universe              = vanilla
executable            = ${WRAPPER}
arguments             = \$(inputFile) \$(ProcId) ${OUTPUT_DIR} ${CMSSW_BASE_DIR} ${MAX_EVENTS}
output                = ${LOG_DIR}/job_\$(ProcId).out
error                 = ${LOG_DIR}/job_\$(ProcId).err
log                   = ${LOG_DIR}/condor.log
+JobFlavour           = "tomorrow"
request_memory        = 4000
request_cpus          = 1
request_disk          = 10000000
should_transfer_files = YES
transfer_output_files = ""
queue inputFile from ${SUBMIT_LIST}
EOF

[ ${TEST_MODE} -eq 1 ] && rm -f "${SUBMIT_LIST}"

echo ""
echo "Submitted ${N_JOBS} '${VARIANT}' jobs."
echo "Monitor:  condor_q"
echo "Logs in:  ${LOG_DIR}/"
echo ""
echo "When all jobs finish, merge ntuples with:"
echo "  hadd pfObjectsNtuple_${VARIANT}_merged.root ${OUTPUT_DIR}/pfObjectsNtuple_job*.root"
