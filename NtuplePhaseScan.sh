#!/bin/bash

cmsenv

INPUT_DIR="/eos/user/g/gkopp/PF_PhaseScan"
OUTPUT_DIR="/eos/user/g/gkopp/PF_PhaseScan"
OUTPUT_FILE="${OUTPUT_DIR}/pfObjectsNtuple_phaseScan.root"

# Build comma-separated list of xrootd input files
FILES=$(ls ${INPUT_DIR}/pf_only_reReco_phaseScan_job*.root 2>/dev/null \
    | sed 's|/eos/|root://eosuser.cern.ch//eos/|' \
    | paste -sd,)

if [ -z "$FILES" ]; then
    echo "ERROR: No phase scan files found in ${INPUT_DIR}"
    exit 1
fi

NFILES=$(echo "$FILES" | tr ',' '\n' | wc -l)
echo "Found ${NFILES} phase scan file(s)"
echo "Output: ${OUTPUT_FILE}"

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles="${FILES}" \
    outputFile=${OUTPUT_FILE}
