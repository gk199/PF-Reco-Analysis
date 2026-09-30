#!/usr/bin/env bash
set -eo pipefail
JOB_DIR=$PWD
mkdir -p diagnostics
exec > >(tee -a "$JOB_DIR/diagnostics/merge.log") 2>&1
trap 'rc=$?; printf "exit_code=%s\n" "$rc" > "$JOB_DIR/diagnostics/status.txt"' EXIT
set -a
source "$JOB_DIR/settings.sh"
set +a
source /cvmfs/cms.cern.ch/cmsset_default.sh
mkdir runtime
cd runtime
scram project CMSSW "$CMSSW_VERSION"
cd "$CMSSW_VERSION/src"
eval "$(scram runtime -sh)"
cd "$JOB_DIR"
[[ -n "${X509_USER_PROXY:-}" && -r "$X509_USER_PROXY" ]] || { echo 'Delegated proxy missing'; exit 1; }
python3 merge_outputs.py
