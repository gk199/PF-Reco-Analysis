#!/bin/bash
set -eo pipefail
NEVENTS=100
R_HCAL=177.7
TIMING_THRESHOLD=4.0

# Etas to scan. Each eta gets its own GEN-SIM folder and its own re-RECO folder.
ETA_LIST=(0.1 0.2 0.4)

# Base output locations; per-eta subfolders are derived from these.
GENSIM_BASE="/eos/user/c/chtong/Public/Public_PFSimulations"
RECO_BASE="/eos/user/c/chtong/Public/Rereco"
RECO_SUFFIX="_${TIMING_THRESHOLD%.0}ns"   # 4.0 -> "_4ns", 2.5 -> "_2.5ns"

gensim_dir() { echo "${GENSIM_BASE}/Dipion_eta${1}"; }
out_dir()    { echo "${RECO_BASE}/Dipion_eta${1}${RECO_SUFFIX}"; }

# This combined script is intended to be run from:
# /afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/
# The generation part is normally run from src/, so we cd .. before cmsDriver GEN/SIM.
if [[ "$(basename "$PWD")" != "PF-Reco-Analysis" ]]; then
  echo "ERROR: Please run this script from:"
  echo "  /afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/"
  exit 1
fi

cd ..
cmsenv
set -u   # enabled after cmsenv, since the CMSSW environment setup may touch unset variables

# Scan grid, run at every eta in ETA_LIST: offsetFirst, DeltaT (ns), DeltaR, pT (GeV)
oF_DT_DR_PT_PAIRS=(
  "0.0 0.0 0.1 20.0"
  "0.0 0.0 0.2 20.0"
  "0.0 0.0 0.3 20.0"
  "0.0 0.0 0.4 20.0"
  "0.0 0.0 0.5 20.0"

  "0.0 1.0 0.1 20.0"
  "0.0 1.0 0.2 20.0"
  "0.0 1.0 0.3 20.0"
  "0.0 1.0 0.4 20.0"
  "0.0 1.0 0.5 20.0"

  "0.0 2.0 0.1 20.0"
  "0.0 2.0 0.2 20.0"
  "0.0 2.0 0.3 20.0"
  "0.0 2.0 0.4 20.0"
  "0.0 2.0 0.5 20.0"

  "0.0 3.0 0.1 20.0"
  "0.0 3.0 0.2 20.0"
  "0.0 3.0 0.3 20.0"
  "0.0 3.0 0.4 20.0"
  "0.0 3.0 0.5 20.0"

  "0.0 4.0 0.1 20.0"
  "0.0 4.0 0.2 20.0"
  "0.0 4.0 0.3 20.0"
  "0.0 4.0 0.4 20.0"
  "0.0 4.0 0.5 20.0"

  "0.0 5.0 0.1 20.0"
  "0.0 5.0 0.2 20.0"
  "0.0 5.0 0.3 20.0"
  "0.0 5.0 0.4 20.0"
  "0.0 5.0 0.5 20.0")

# Keep a list of the exact samples produced so re-RECO/ntupling uses matching
# input/output names. Each entry is "<eta> <sample>" so later loops can find
# the right per-eta folder.
SAMPLES_TO_PROCESS=()

echo
echo "========================================"
echo "  GEN-SIM + DIGI for etas: ${ETA_LIST[*]}"
echo "========================================"

for ETA in "${ETA_LIST[@]}"; do
  GENSIM_DIR=$(gensim_dir "${ETA}")
  OUT_DIR=$(out_dir "${ETA}")
  mkdir -p "${GENSIM_DIR}" "${OUT_DIR}"

  # z = R * sinh(eta), rounded to 0.1 cm so it matches the filename tag exactly
  Z=$(python3 -c "import math; print(f'{${R_HCAL} * math.sinh(${ETA}):.1f}')")

  echo
  echo "######## eta=${ETA}  ->  Z=${Z} cm  (R=${R_HCAL} cm) ########"
  echo "  GEN-SIM dir: ${GENSIM_DIR}"
  echo "  re-RECO dir: ${OUT_DIR}"

  for PAIR in "${oF_DT_DR_PT_PAIRS[@]}"; do
    read -r oF DT DR PT <<< "$PAIR"

    TAG="Z${Z}_offsetFirst${oF}_DT${DT}_DR${DR}_pT${PT}_eta${ETA}"
    SAMPLE="DiPi_${TAG}"
    SAMPLES_TO_PROCESS+=("${ETA} ${SAMPLE}")
    echo "=== Z=${Z}, offsetFirst=${oF}, DT=${DT} ns, DR=${DR}, pT=${PT}, eta=${ETA} ==="

    # Z, ETA, PT, and Delta setup
    ZMIN=$(python3 -c "print(${Z} - 0.000000001)")
    ZMAX=$(python3 -c "print(${Z} + 0.000000001)")
    MINETA=$(python3 -c "print(${ETA} - 0.000000001)")
    MAXETA=$(python3 -c "print(${ETA} + 0.000000001)")

    DELTA=$(python3 -c "import math; print(${DR} * ${Z} / math.sinh(${ETA}))")

    echo "  ZMIN=${ZMIN}, ZMAX=${ZMAX}, DELTA=${DELTA}, PT = ${PT}, ETA = ${ETA}"

    Z_CMD="process.generator.PGunParameters.ZMin=cms.double(${ZMIN});process.generator.PGunParameters.ZMax=cms.double(${ZMAX})"
    ETA_CMD="process.generator.PGunParameters.MinEta=cms.double(${MINETA});process.generator.PGunParameters.MaxEta=cms.double(${MAXETA})"
    PT_CMD="process.generator.PGunParameters.VarMin=cms.double(${PT});process.generator.PGunParameters.VarMax=cms.double(${PT})"
    DELTA_CMD="process.generator.PGunParameters.Delta=cms.double(${DELTA})"

    # DeltaT timing setup
    if python3 -c "import sys; sys.exit(0 if float('${DT}') == 0 else 1)"; then
      TIMING_CMD="process.generator.PGunParameters.UseDeltaT=cms.bool(False)"
    else
      TMIN=$(python3 -c "print(${DT} - 0.000000001)")
      TMAX=$(python3 -c "print(${DT} + 0.000000001)")
      TIMING_CMD="process.generator.PGunParameters.UseDeltaT=cms.bool(True);process.generator.PGunParameters.TMin=cms.double(${TMIN});process.generator.PGunParameters.TMax=cms.double(${TMAX})"
    fi

    # OffsetFirst setup
    OFFSET_CMD="process.generator.PGunParameters.OffsetFirst=cms.double(${oF})"

    cmsDriver.py Configuration/Generator/python/DipionGun_Eta_cfi.py \
      --mc \
      --step GEN,SIM \
      --era Run3_2024 \
      --geometry DB:Extended \
      --conditions auto:phase1_2024_realistic \
      --beamspot Realistic25ns13p6TeVEarly2022Collision \
      --eventcontent RAWSIM \
      --datatier GEN-SIM \
      --customise_commands "${TIMING_CMD};${Z_CMD};${ETA_CMD};${PT_CMD};${DELTA_CMD};${OFFSET_CMD}" \
      --fileout "file:${GENSIM_DIR}/${SAMPLE}_GEN-SIM.root" \
      --python_filename "${SAMPLE}_cfg.py" \
      -n ${NEVENTS}

    cmsDriver.py step1 \
      --python_filename "${SAMPLE}_step1_cfg.py" \
      --filein "file:${GENSIM_DIR}/${SAMPLE}_GEN-SIM.root" \
      --fileout "file:${GENSIM_DIR}/${SAMPLE}_step1_GEN-SIM-RAW.root" \
      --pileup NoPileUp \
      --customise Configuration/DataProcessing/Utils.addMonitoring \
      --eventcontent RAWSIM \
      --datatier GEN-SIM-RAW \
      --conditions auto:phase1_2024_realistic \
      --step DIGI,L1,DIGI2RAW \
      --geometry DB:Extended \
      --era Run3_2024 \
      --mc \
      -n ${NEVENTS}

  done
done

echo "All GEN-SIM and GEN-SIM-RAW files are done (${#SAMPLES_TO_PROCESS[@]} samples)."

# ---------------------------------------------------------------------------
# re-RECO of one sample with whichever PF variant is currently compiled.
# Usage: run_rereco <eta> <sample> <variant-name>
# Must be called from PF-Reco-Analysis/.
# ---------------------------------------------------------------------------
run_rereco() {
  local ETA="$1" SAMPLE="$2" VARIANT="$3"
  local GENSIM_DIR OUT_DIR
  GENSIM_DIR=$(gensim_dir "${ETA}")
  OUT_DIR=$(out_dir "${ETA}")
  local INFILE="${GENSIM_DIR}/${SAMPLE}_step1_GEN-SIM-RAW.root"
  local OUTFILE="${OUT_DIR}/pf_only_reReco_${SAMPLE}_${VARIANT}.root"

  echo
  echo "-----------------------------------"
  echo "Running re-RECO (${VARIANT}): ${SAMPLE}"
  echo "Input:  ${INFILE}"
  echo "Output: ${OUTFILE}"
  echo "-----------------------------------"
  echo

  cmsDriver.py MyPFStudy_ReReco_MC_Sim \
    --mc \
    --conditions auto:phase1_2025_realistic \
    --step RAW2DIGI,L1Reco,RECO \
    --geometry DB:Extended \
    --era Run3 \
    --filein "file:${INFILE}" \
    --fileout file:pf_only_reReco_MC_Sim.root \
    --eventcontent RECO \
    --datatier RECO \
    --process ReRECO \
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco_MC_Sim.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_particleFlowRecHit*_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*'))" \
    --no_exec \
    -n "${NEVENTS}"

  rm -f pf_only_reReco_MC_Sim.root
  cmsRun MyPFStudy_ReReco_MC_Sim_RAW2DIGI_L1Reco_RECO.py
  mv pf_only_reReco_MC_Sim.root "${OUTFILE}"
}

cd PF-Reco-Analysis
./SetTimingThreshold.sh "${TIMING_THRESHOLD}"
echo "========================================"
echo "  Timing threshold: ${TIMING_THRESHOLD} ns"
echo "========================================"
cd ..

echo
echo "Running all ParticleFlow tests on ${NEVENTS} events MC sample"
echo

echo "-----------------------------------"
echo "Standard Particle Flow:"
echo "-----------------------------------"
echo

cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit \
   RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc

cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit \
   RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc

cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit \
   RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc

cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py \
   RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8

cd PF-Reco-Analysis

for ENTRY in "${SAMPLES_TO_PROCESS[@]}"; do
  read -r ETA SAMPLE <<< "$ENTRY"
  run_rereco "${ETA}" "${SAMPLE}" "standardPF"
done


echo " "
echo "-----------------------------------"
echo "Modified Particle Flow: seed level timing cut vs global highest energy seed"
echo "-----------------------------------"
echo " "

cd ..

cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc

echo "Note: only timing cut on seeds, not on all cells, so only the clusterizer and cluster producer are modified"

cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8

cd PF-Reco-Analysis

for ENTRY in "${SAMPLES_TO_PROCESS[@]}"; do
  read -r ETA SAMPLE <<< "$ENTRY"
  run_rereco "${ETA}" "${SAMPLE}" "seedTimingPF"
done


echo " "
echo "All tests completed, output files are in:"
for ETA in "${ETA_LIST[@]}"; do
  echo "  $(out_dir "${ETA}")"
done
echo " "

echo "-----------------------------------"
echo "Reverting to standard Particle Flow configuration"
echo "-----------------------------------"
echo " "

cd ..

cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8

cd PF-Reco-Analysis

echo " "
echo "-----------------------------------"
echo "Reverted to standard Particle Flow configuration"
echo "Done!"
echo "-----------------------------------"


echo
echo "-----------------------------------"
echo "Producing PFObjects ntuples"
echo "-----------------------------------"
echo

PF_VARIANTS=("standardPF" "seedTimingPF")

for PF in "${PF_VARIANTS[@]}"; do
  for ENTRY in "${SAMPLES_TO_PROCESS[@]}"; do
    read -r ETA SAMPLE <<< "$ENTRY"
    OUT_DIR=$(out_dir "${ETA}")
    INFILE="${OUT_DIR}/pf_only_reReco_${SAMPLE}_${PF}.root"
    OUTFILE="pfObjectsNtuple_${PF}_${SAMPLE}.root"

    echo
    echo "PF:     ${PF}"
    echo "Input:  ${INFILE}"
    echo "Output: ${OUT_DIR}/${OUTFILE}"
    echo

    rm -f pfObjectsNtuple.root

    cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
      inputFiles="file:${INFILE}" \
      outputFile="pfObjectsNtuple.root"

    mv pfObjectsNtuple.root "${OUT_DIR}/${OUTFILE}"
  done
done



echo "Done."