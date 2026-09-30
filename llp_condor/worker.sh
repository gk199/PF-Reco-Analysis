#!/usr/bin/env bash
set -eo pipefail
JOB_DIR=$PWD
mkdir -p "$JOB_DIR/diagnostics"
exec > >(tee -a "$JOB_DIR/diagnostics/worker.log") 2>&1
trap 'rc=$?; printf "exit_code=%s\n" "$rc" > "$JOB_DIR/diagnostics/status.txt"' EXIT
source "$JOB_DIR/settings.sh"
INDEX=$1
TAG=$2
INPUT_URL=$(sed -n "$((INDEX + 1))p" "$JOB_DIR/selected_inputs.txt")
[[ -n "$INPUT_URL" ]] || { echo "No input at index $INDEX"; exit 1; }
cp "$JOB_DIR/settings.sh" "$JOB_DIR/diagnostics/settings.sh"
printf '%s\n' "$INPUT_URL" > "$JOB_DIR/diagnostics/input.txt"

export SCRAM_ARCH
source /cvmfs/cms.cern.ch/cmsset_default.sh
[[ -n "${X509_USER_PROXY:-}" && -r "$X509_USER_PROXY" ]] || { echo 'Delegated X509 proxy missing'; exit 1; }
voms-proxy-info -timeleft

mkdir "$JOB_DIR/runtime"
cd "$JOB_DIR/runtime"
scram project CMSSW "$CMSSW_VERSION"
cd "$CMSSW_VERSION"
tar -xzf "$JOB_DIR/source.tgz"
cd src
eval "$(scram runtime -sh)"
set -u
BUILD_SRC=$PWD
ANALYSIS="$BUILD_SRC/PF-Reco-Analysis"

# Stage once; both variants see precisely the same local input events.
retry_copy() {
    local source=$1 destination=$2
    local attempt
    for attempt in 1 2 3; do
        if xrdcp --force --cksum adler32 --posc "$source" "$destination"; then return 0; fi
        echo "Transfer attempt $attempt failed: $source -> $destination"
        [[ "$attempt" == 3 ]] || sleep 10
    done
    return 1
}
retry_copy "$INPUT_URL" "$JOB_DIR/input.root"
edmDumpEventContent "$JOB_DIR/input.root" > "$JOB_DIR/diagnostics/input_content.txt"
for label in hbhereco hfreco horeco ecalRecHit ecalPreshowerRecHit generalTracks genParticles; do
    grep -q "\"$label\"" "$JOB_DIR/diagnostics/input_content.txt" || { echo "Required collection absent: $label"; exit 1; }
done

# Match the threshold-setting order in the supplied local script.
cd "$ANALYSIS"
bash ./SetTimingThreshold.sh "$TIMING_THRESHOLD"

for VARIANT in standardPF seedTimingPF; do
    cd "$BUILD_SRC"
    PF_SRC="$ANALYSIS/PFTestingAlgos"
    PF_DST="$BUILD_SRC/RecoParticleFlow/PFClusterProducer"
    cp "$PF_SRC/Basic2DGenericTopoClusterizer_original.cc.edit" "$PF_DST/plugins/Basic2DGenericTopoClusterizer.cc"
    if [[ "$VARIANT" == standardPF ]]; then
        cp "$PF_SRC/PFMultiDepthClusterProducer_original.cc.edit" "$PF_DST/plugins/PFMultiDepthClusterProducer.cc"
        cp "$PF_SRC/PFMultiDepthClusterizer_original.cc.edit" "$PF_DST/plugins/PFMultiDepthClusterizer.cc"
        cp "$PF_SRC/particleFlowClusterHCAL_original_cfi.py" "$PF_DST/python/particleFlowClusterHCAL_cfi.py"
        EOS_VARIANT=standardPF
    else
        cp "$PF_SRC/PFMultiDepthClusterProducer_timing.cc.edit" "$PF_DST/plugins/PFMultiDepthClusterProducer.cc"
        cp "$PF_SRC/PFMultiDepthClusterizer_seedTiming.cc.edit" "$PF_DST/plugins/PFMultiDepthClusterizer.cc"
        cp "$PF_SRC/particleFlowClusterHCAL_seedTiming_cfi.py" "$PF_DST/python/particleFlowClusterHCAL_cfi.py"
        EOS_VARIANT="seedTimingPF_${TIMING_TAG}"
    fi
    echo "Building $VARIANT with $NCPUS CPUs"
    scram b -j "$NCPUS" 2>&1 | tee "$JOB_DIR/diagnostics/build_${VARIANT}.log"
    sha256sum "$PF_DST/plugins/Basic2DGenericTopoClusterizer.cc" \
        "$PF_DST/plugins/PFMultiDepthClusterProducer.cc" \
        "$PF_DST/plugins/PFMultiDepthClusterizer.cc" \
        "$PF_DST/python/particleFlowClusterHCAL_cfi.py" \
        > "$JOB_DIR/diagnostics/sources_${VARIANT}.sha256"
    cp "$PF_DST/python/particleFlowClusterHCAL_cfi.py" "$JOB_DIR/diagnostics/hcal_${VARIANT}.py"

    WORK="$JOB_DIR/$VARIANT"
    mkdir "$WORK"
    cd "$WORK"
    RECO_NAME="pf_reReco_${TAG}_${VARIANT}.root"
    NTUPLE_NAME="pfObjectsNtuple_${TAG}_${VARIANT}.root"

    # The inspected RECOSIM input lacks rawDataCollector. Re-run only PF rechits
    # and PF clusters from the stored RECO rechits; blocks/candidates stay from input RECO.
    cmsDriver.py LLPReReco \
        --mc --conditions "$CONDITIONS" --era "$ERA" \
        --geometry DB:Extended \
        --step RECO:particleFlowCluster \
        --eventcontent RECO --datatier RECO --process ReRECO \
        --filein "file:$JOB_DIR/input.root" --fileout "file:$RECO_NAME" \
        --python_filename rereco_cfg.py --no_exec \
        --nThreads "$NCPUS" -n "$MAX_EVENTS"

    # Add the supplied analyzer in the same process, after reconstruction.
    # It consumes ReRECO explicitly; old input PF products cannot substitute.
    export LLP_NTUPLE_FILE="$NTUPLE_NAME"
    cat "$JOB_DIR/add_ntupler.py" >> rereco_cfg.py
    cp rereco_cfg.py "$JOB_DIR/diagnostics/rereco_${VARIANT}_cfg.py"
    edmConfigDump rereco_cfg.py > "$JOB_DIR/diagnostics/expanded_${VARIANT}_cfg.py"
    cmsRun -j "$JOB_DIR/diagnostics/report_${VARIANT}.xml" rereco_cfg.py \
        2>&1 | tee "$JOB_DIR/diagnostics/cmsRun_${VARIANT}.log"
    test -s "$RECO_NAME"
    test -s "$NTUPLE_NAME"
    edmFileUtil "$RECO_NAME" > "$JOB_DIR/diagnostics/events_${VARIANT}.txt"
    edmDumpEventContent "$RECO_NAME" > "$JOB_DIR/diagnostics/content_${VARIANT}.txt"
    python3 - "$JOB_DIR/diagnostics/content_${VARIANT}.txt" <<'PY'
import pathlib, sys
lines = pathlib.Path(sys.argv[1]).read_text().splitlines()
for label in ('particleFlowRecHitHBHE', 'particleFlowClusterHCAL', 'particleFlowClusterECAL'):
    if not any(f'"{label}"' in s and '"ReRECO"' in s for s in lines):
        raise SystemExit(f'Missing freshly reconstructed collection: {label}:ReRECO')
PY
    python3 "$JOB_DIR/root_counts.py" "$NTUPLE_NAME" > "$JOB_DIR/diagnostics/counts_${VARIANT}.json"
    if [[ "$SAVE_RECO" == 1 ]]; then
        python3 "$JOB_DIR/root_counts.py" "$RECO_NAME" > "$JOB_DIR/diagnostics/edm_counts_${VARIANT}.json"
        retry_copy "$RECO_NAME" "${EOS_HOST}/${EOS_RUN}/parts/$EOS_VARIANT/$RECO_NAME"
    fi
    retry_copy "$NTUPLE_NAME" "${EOS_HOST}/${EOS_RUN}/parts/$EOS_VARIANT/$NTUPLE_NAME"
    # Reclaim only this job's generated files, after both verified transfers.
    rm -- "$RECO_NAME" "$NTUPLE_NAME"
done

python3 - "$JOB_DIR/diagnostics" "$TAG" "$INPUT_URL" "$MAX_EVENTS" "$TIMING_THRESHOLD" "$SAVE_RECO" <<'PY'
import json, pathlib, sys
directory, tag, url, maximum, threshold, save = sys.argv[1:]
p = pathlib.Path(directory)
record = dict(tag=tag, input=url, max_events=int(maximum), threshold=float(threshold), save_reco=bool(int(save)))
record['trees'] = {v: json.loads((p / f'counts_{v}.json').read_text()) for v in ('standardPF', 'seedTimingPF')}
if record['save_reco']:
    record['edm_events'] = {v: json.loads((p / f'edm_counts_{v}.json').read_text())['Events'] for v in ('standardPF', 'seedTimingPF')}
(p / f'{tag}.done.json').write_text(json.dumps(record, indent=2) + '\n')
PY
retry_copy "$JOB_DIR/diagnostics/${TAG}.done.json" "${EOS_HOST}/${EOS_RUN}/completed/${TAG}.done.json"
echo "Both variants completed: $EOS_RUN"
