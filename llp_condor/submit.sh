#!/usr/bin/env bash
# Run after cmsenv in the CMSSW area that contains your PF modifications.
# voms-proxy-init --voms cms --valid 96:00 
# bash llp_condor/submit.sh
set -euo pipefail
HERE=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
: "${CMSSW_BASE:?Run cmsenv in your existing CMSSW area first}"
: "${CMSSW_VERSION:?Run cmsenv first}"
: "${SCRAM_ARCH:?Run cmsenv first}"

INPUT_LIST=${INPUT_LIST:-"$HERE/input_files_MH125_MS50_CTau3000.txt"}
#INPUT_LIST=${INPUT_LIST:-"$HERE/input_files_MH350_MS160_CTau10000.txt"}
NFILES=${NFILES:-20}                 # 0 selects the whole input list
FILE_OFFSET=${FILE_OFFSET:-0}      # Skip this many unique input entries first
MAX_EVENTS=${MAX_EVENTS:--1}       # -1 processes the entire input file
MERGE_OUTPUTS=${MERGE_OUTPUTS:-1}  # One final ntuple per PF algorithm
SAVE_RECO=${SAVE_RECO:-1}          # Also save separate merged EDM files
KEEP_PARTS=${KEEP_PARTS:-0}        # Keep per-input files after a successful merge
TIMING_THRESHOLD=${TIMING_THRESHOLD:-2.0}
NCPUS=${NCPUS:-4}
MEMORY_MB=${MEMORY_MB:-12000}
DISK_MB=${DISK_MB:-40000}
JOB_FLAVOUR=${JOB_FLAVOUR:-tomorrow}
CONDITIONS=${CONDITIONS:-auto:phase1_2023_realistic_postBPix}
ERA=${ERA:-Run3_2023}
EOS_BASE=${EOS_BASE:-/eos/user/c/chtong/Public/Rereco/LLP_MH125_MS50_CTau3000}
#EOS_BASE=${EOS_BASE:-/eos/user/c/chtong/Public/Rereco/LLP_MH350_MS160_CTau10000}
EOS_HOST=${EOS_HOST:-root://eosuser.cern.ch}
SOURCE_HOST=${SOURCE_HOST:-root://cmseos.fnal.gov}
RUN_ID=${RUN_ID:-"$(date -u +%Y%m%dT%H%M%SZ)_${RANDOM}"}
BATCH_BASE=${BATCH_BASE:-"$(dirname "$CMSSW_BASE")/src/PF-Reco-Analysis/llp_condor_runs"}
PREPARE_ONLY=${PREPARE_ONLY:-0}

die() { echo "ERROR: $*" >&2; exit 1; }
[[ "$CMSSW_VERSION" == CMSSW_15_0_6 ]] || die "This recipe targets CMSSW_15_0_6; active release is $CMSSW_VERSION."
[[ "$SCRAM_ARCH" == el[89]_amd64_* ]] || die "Expected an el8/el9 amd64 CMSSW architecture, got $SCRAM_ARCH."
[[ "$NFILES" =~ ^[0-9]+$ ]] || die 'NFILES must be 0 or a positive integer.'
[[ "$FILE_OFFSET" =~ ^[0-9]+$ ]] || die 'FILE_OFFSET must be a nonnegative integer.'
for flag in MERGE_OUTPUTS SAVE_RECO KEEP_PARTS; do
    [[ "${!flag}" == 0 || "${!flag}" == 1 ]] || die "$flag must be 0 or 1."
done
[[ "$MAX_EVENTS" == -1 || "$MAX_EVENTS" =~ ^[1-9][0-9]*$ ]] || die 'MAX_EVENTS must be -1 or positive.'
[[ "$TIMING_THRESHOLD" =~ ^[0-9]+([.][0-9]+)?$ ]] || die 'Invalid timing threshold.'
[[ "$NCPUS" =~ ^[1-9][0-9]*$ ]] || die 'Invalid NCPUS.'
[[ "$MEMORY_MB" =~ ^[1-9][0-9]*$ && "$DISK_MB" =~ ^[1-9][0-9]*$ ]] || die 'Invalid memory/disk request.'
[[ "$RUN_ID" =~ ^[a-zA-Z0-9_-]+$ ]] || die 'Use only letters, numbers, underscores and hyphens in RUN_ID.'
[[ -n "${INPUT_FILES:-}" || -f "$INPUT_LIST" ]] || die "Cannot read $INPUT_LIST"
for cmd in tar python3 condor_submit voms-proxy-info xrdfs; do command -v "$cmd" >/dev/null || die "Missing $cmd"; done
if [[ "$MERGE_OUTPUTS" == 1 ]]; then command -v condor_submit_dag >/dev/null || die 'Missing condor_submit_dag'; fi

SOURCE_SRC=$(realpath "$CMSSW_BASE/src")
BATCH_BASE=$(realpath -m "$BATCH_BASE")
case "$BATCH_BASE/" in "$SOURCE_SRC/"*) die 'BATCH_BASE must be outside CMSSW_BASE/src, to avoid recursive packaging.';; esac
case "$BATCH_BASE" in /eos/*) die 'Put BATCH_BASE on AFS/work, not EOS; it holds Condor logs and submission files.';; esac
[[ "$BATCH_BASE" != *[[:space:]]* ]] || die 'BATCH_BASE must not contain whitespace.'
ANALYSIS="$SOURCE_SRC/PF-Reco-Analysis"
[[ -f "$ANALYSIS/SetTimingThreshold.sh" ]] || die "Missing $ANALYSIS/SetTimingThreshold.sh"
[[ -d "$SOURCE_SRC/RecoParticleFlow/PFClusterProducer/plugins" ]] || die 'The local PFClusterProducer package is missing.'
for f in Basic2DGenericTopoClusterizer_original.cc.edit PFMultiDepthClusterProducer_original.cc.edit PFMultiDepthClusterizer_original.cc.edit particleFlowClusterHCAL_original_cfi.py PFMultiDepthClusterProducer_timing.cc.edit PFMultiDepthClusterizer_seedTiming.cc.edit particleFlowClusterHCAL_seedTiming_cfi.py; do
    [[ -f "$ANALYSIS/PFTestingAlgos/$f" ]] || die "Missing PFTestingAlgos/$f"
done

# HTCondor delegates the proxy; do not rely on /tmp on the worker.
PROXY=$(voms-proxy-info -path) || die 'Create a CMS proxy first: voms-proxy-init --voms cms --valid 96:00'
voms-proxy-info -exists -valid 30:00 || die 'Need at least 30 hours remaining on the proxy; renew it before submission.'
mkdir -p "$BATCH_BASE"
RUN_DIR="$BATCH_BASE/$RUN_ID"
mkdir "$RUN_DIR" || die "Run directory already exists: $RUN_DIR"
mkdir "$RUN_DIR/logs" "$RUN_DIR/reports"
cp "$PROXY" "$RUN_DIR/x509up"
chmod 600 "$RUN_DIR/x509up"
cp "$HERE/worker.sh" "$HERE/add_ntupler.py" "$HERE/merge.sh" "$HERE/merge_outputs.py" "$HERE/root_counts.py" "$RUN_DIR/"
chmod +x "$RUN_DIR/worker.sh" "$RUN_DIR/merge.sh"

python3 "$HERE/select_inputs.py" "$INPUT_LIST" "$NFILES" "$FILE_OFFSET" "$SOURCE_HOST" "$RUN_DIR"

EOS_RUN="${EOS_BASE%/}/$RUN_ID"
TIMING_TAG="${TIMING_THRESHOLD//./p}ns"
for variant in standardPF "seedTimingPF_${TIMING_TAG}"; do
    xrdfs "$EOS_HOST" mkdir -p "$EOS_RUN/parts/$variant"
done
xrdfs "$EOS_HOST" mkdir -p "$EOS_RUN/completed"

# Snapshot all checked-out source packages (including uncommitted changes).
# ROOT outputs and common generated artifacts are not shipped or compiled.
# Exclude this downloaded bundle if it was unpacked inside src/.
tar --dereference --exclude=.git --exclude=__pycache__ --exclude=.ipynb_checkpoints \
    --exclude='*.root' --exclude='*.root.*' --exclude='*.tgz' --exclude='*.tar.gz' \
    --exclude='*.log' --exclude='*.out' --exclude='*.err' \
    --exclude='*.png' --exclude='*.pdf' --exclude='*.jpg' \
    --exclude='llp_condor' --exclude='llp_condor_runs' \
    -czf "$RUN_DIR/source.tgz" -C "$CMSSW_BASE" src

WANT_OS=${SCRAM_ARCH%%_*}
{
    for key in CMSSW_VERSION SCRAM_ARCH MAX_EVENTS TIMING_THRESHOLD TIMING_TAG NCPUS CONDITIONS ERA EOS_RUN EOS_HOST MERGE_OUTPUTS SAVE_RECO KEEP_PARTS; do
        printf '%s=%q\n' "$key" "${!key}"
    done
} > "$RUN_DIR/settings.sh"

cat > "$RUN_DIR/submit.sub" <<EOF
universe = vanilla
executable = worker.sh
arguments = \$(index) \$(tag)
output = logs/\$(tag).\$(ClusterId).\$(ProcId).out
error = logs/\$(tag).\$(ClusterId).\$(ProcId).err
log = logs/cluster.\$(ClusterId).log
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = source.tgz,settings.sh,selected_inputs.txt,add_ntupler.py,root_counts.py
transfer_output_files = diagnostics
transfer_output_remaps = "diagnostics=reports/\$(tag).\$(ClusterId).\$(ProcId)"
use_x509userproxy = True
x509userproxy = $RUN_DIR/x509up
delegate_job_GSI_credentials_lifetime = 0
request_cpus = $NCPUS
request_memory = ${MEMORY_MB}MB
request_disk = ${DISK_MB}MB
MY.WantOS = "$WANT_OS"
+JobFlavour = "$JOB_FLAVOUR"
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
queue index,tag from jobs.txt
EOF

if [[ "$MERGE_OUTPUTS" == 1 ]]; then
    # One DAG node per input: a failed input cannot trigger an incomplete merge.
    sed 's/^queue index,tag from jobs.txt$/queue/' "$RUN_DIR/submit.sub" > "$RUN_DIR/node.sub"
    cat > "$RUN_DIR/merge.sub" <<EOF
universe = vanilla
executable = merge.sh
output = logs/merge.\$(ClusterId).out
error = logs/merge.\$(ClusterId).err
log = logs/merge.\$(ClusterId).log
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = settings.sh,jobs.txt,selected_inputs.txt,merge_outputs.py,root_counts.py
transfer_output_files = diagnostics
transfer_output_remaps = "diagnostics=reports/merge.\$(ClusterId)"
use_x509userproxy = True
x509userproxy = $RUN_DIR/x509up
delegate_job_GSI_credentials_lifetime = 0
request_cpus = 1
request_memory = 4000MB
request_disk = ${DISK_MB}MB
MY.WantOS = "$WANT_OS"
+JobFlavour = "$JOB_FLAVOUR"
on_exit_hold = (ExitBySignal == True) || (ExitCode != 0)
queue
EOF
    python3 - "$RUN_DIR" <<'PY'
from pathlib import Path
import sys
p = Path(sys.argv[1])
lines, parents = [], []
for row in (p / 'jobs.txt').read_text().splitlines():
    index, tag = row.split()
    name = 'RECO_' + index
    parents.append(name)
    lines += [f'JOB {name} node.sub', f'VARS {name} index="{index}" tag="{tag}"']
lines.append('JOB MERGE merge.sub')
lines += [f'PARENT {name} CHILD MERGE' for name in parents]
(p / 'workflow.dag').write_text('\n'.join(lines) + '\n')
PY
fi

echo "Submission directory: $RUN_DIR"
echo "EOS output: $EOS_RUN"
echo "Input events per file: $MAX_EVENTS; seed threshold: $TIMING_THRESHOLD ns"
echo "The original CMSSW area is only read; each worker compiles its own copy."
cd "$RUN_DIR"
if [[ "$PREPARE_ONLY" == 1 ]]; then
    if [[ "$MERGE_OUTPUTS" == 1 ]]; then
        echo "Prepared. Submit with: cd '$RUN_DIR' && condor_submit_dag workflow.dag"
    else
        echo "Prepared. Submit with: cd '$RUN_DIR' && condor_submit submit.sub"
    fi
elif [[ "$MERGE_OUTPUTS" == 1 ]]; then
    condor_submit_dag workflow.dag
else
    condor_submit submit.sub
fi
