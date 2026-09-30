#!/bin/bash
set -e

NEVENTS=100
HCAL_RADIUS_CM=177.7
PION_MASS_GEV=0.13957039

# Run this script from:
# /afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/
#
# It produces single-pion samples with:
#   1. CloseByParticleGunProducer, launched from the HCAL barrel surface.
#   2. Pythia8PtGun, launched from the interaction point.


# The script then runs DIGI, the three PF reconstruction variants, and the
# PFObjects ntupler. Plotting is intentionally not included.
if [[ "$(basename "$PWD")" != "PF-Reco-Analysis" ]]; then
  echo "ERROR: Please run this script from:"
  echo "  /afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/"
  exit 1
fi

# Convert a numeric value into a filename-safe tag.
number_tag() {
  python3 -c "x=float('$1'); print(int(x) if x.is_integer() else str(x).replace('-', 'm').replace('.', 'p'))"
}

# Convert a requested total pion energy E and pseudorapidity eta into pT:
#
#   |p| = sqrt(E^2 - m_pi^2)
#   pT  = |p| / cosh(eta)
#
# Pythia8PtGun requires pT as its input. CloseByParticleGunProducer can instead
# use the requested energy directly when FlatPtGeneration is disabled.
energy_to_pt() {
  local energy="$1"
  local eta="$2"

  python3 - "${energy}" "${eta}" "${PION_MASS_GEV}" <<'PYENERGY'
import math
import sys

energy = float(sys.argv[1])
eta = float(sys.argv[2])
mass = float(sys.argv[3])

if energy <= mass:
    raise SystemExit(
        f"ERROR: requested pion energy {energy} GeV must exceed "
        f"the charged-pion mass {mass} GeV"
    )

momentum = math.sqrt(energy * energy - mass * mass)
pt = momentum / math.cosh(eta)
print(f"{pt:.10f}")
PYENERGY
}

cd ..
cmsenv

GENSIM_DIR="/eos/user/c/chtong/Public/Public_PFSimulations/SinglePion_SmallEta_modified"
OUT_DIR="/eos/user/c/chtong/Public/Rereco/SinglePion_SmallEta_modified_5ns"

mkdir -p "${GENSIM_DIR}"
mkdir -p "${OUT_DIR}"

# Keep the exact sample names so the subsequent re-RECO and ntupling loops
# use matching input and output filenames.
SAMPLES_TO_PROCESS=()
# -----------------------------------------------------------------------------
# CloseByParticleGunProducer single-pion samples
# -----------------------------------------------------------------------------
# The first value in each pair is the desired TOTAL pion energy in GeV.
# It is also used in the sample filename. For CloseByParticleGunProducer,
# VarMin/VarMax are set directly to this energy and FlatPtGeneration is disabled.
#
# For a cylindrical HCAL barrel surface at radius R_HCAL, the matching z
# coordinate for a particle at pseudorapidity eta is:
#
#   z = R_HCAL * sinh(eta)
#
# Setting Z this way makes the generated vertex radius approximately 177 cm for
# every eta below. DeltaR and DeltaT are not scanned because NParticles = 1.
CLOSEBY_ENERGY_ETA_PAIRS=(
  "20.0 0.01"
  "20.0 0.1"
  "20.0 0.2"
  "20.0 0.4"
  "20.0 0.6"
  "20.0 0.8"
  "20.0 1.0"

  "40.0 0.01"
  "40.0 0.1"
  "40.0 0.2"
  "40.0 0.4"
  "40.0 0.6"
  "40.0 0.8"
  "40.0 1.0"

  "60.0 0.01"
  "60.0 0.1"
  "60.0 0.2"
  "60.0 0.4"
  "60.0 0.6"
  "60.0 0.8"
  "60.0 1.0"

  "80.0 0.01"
  "80.0 0.1"
  "80.0 0.2"
  "80.0 0.4"
  "80.0 0.6"
  "80.0 0.8"
  "80.0 1.0"

  "100.0 0.01"
  "100.0 0.1"
  "100.0 0.2"
  "100.0 0.4"
  "100.0 0.6"
  "100.0 0.8"
  "100.0 1.0"

  "120.0 0.01"
  "120.0 0.1"
  "120.0 0.2"
  "120.0 0.4"
  "120.0 0.6"
  "120.0 0.8"
  "120.0 1.0")

for PAIR in "${CLOSEBY_ENERGY_ETA_PAIRS[@]}"; do
  read -r ENERGY ETA <<< "${PAIR}"

  ENERGY_TAG=$(number_tag "${ENERGY}")
  ETA_TAG=$(number_tag "${ETA}")
  SAMPLE="SinglePiCloseByE${ENERGY_TAG}_eta${ETA_TAG}"
  SAMPLES_TO_PROCESS+=("${SAMPLE}")

  # # This pT is calculated only for logging/cross-checking. The CloseBy gun is
  # # configured using ENERGY directly below.
  # PT=$(energy_to_pt "${ENERGY}" "${ETA}")

  # # Choose Z so the generated vertex lies at R = HCAL_RADIUS_CM.
  # Z=$(python3 -c "import math; print(f'{float(${HCAL_RADIUS_CM}) * math.sinh(float(${ETA})):.6f}')")
  # ZMIN=$(python3 -c "print(float(${Z}) - 0.000000000001)")
  # ZMAX=$(python3 -c "print(float(${Z}) + 0.000000000001)")
  # MINETA=$(python3 -c "print(float(${ETA}) - 0.000000000001)")
  # MAXETA=$(python3 -c "print(float(${ETA}) + 0.000000000001)")

  # echo "=============================================="
  # echo "Running ${SAMPLE}"
  # echo "  generator = CloseByParticleGunProducer"
  # echo "  energy    = ${ENERGY} GeV"
  # echo "  eta       = ${ETA}"
  # echo "  derived pT= ${PT} GeV"
  # echo "  HCAL R    = ${HCAL_RADIUS_CM} cm"
  # echo "  vertex Z  = ${Z} cm"
  # echo "  particles = 1 pi+"
  # echo "=============================================="

  # Z_CMD="process.generator.PGunParameters.ZMin=cms.double(${ZMIN});process.generator.PGunParameters.ZMax=cms.double(${ZMAX})"
  # ETA_CMD="process.generator.PGunParameters.MinEta=cms.double(${MINETA});process.generator.PGunParameters.MaxEta=cms.double(${MAXETA})"
  # ENERGY_CMD="process.generator.PGunParameters.VarMin=cms.double(${ENERGY});process.generator.PGunParameters.VarMax=cms.double(${ENERGY});process.generator.PGunParameters.FlatPtGeneration=cms.bool(False)"
  # SINGLE_PION_CMD="process.generator.PGunParameters.NParticles=cms.int32(1);process.generator.PGunParameters.RandomShoot=cms.bool(False);process.generator.PGunParameters.PartID=cms.vint32(211);process.generator.AddAntiParticle=cms.bool(False)"
  # NO_RELATIVE_TIMING_CMD="process.generator.PGunParameters.UseDeltaT=cms.bool(False);process.generator.PGunParameters.OffsetFirst=cms.double(0.0)"

  # cmsDriver.py Configuration/Generator/python/DipionGun_Eta_cfi.py \
  #   --mc \
  #   --step GEN,SIM \
  #   --era Run3_2024 \
  #   --geometry DB:Extended \
  #   --conditions auto:phase1_2024_realistic \
  #   --beamspot Realistic25ns13p6TeVEarly2022Collision \
  #   --eventcontent RAWSIM \
  #   --datatier GEN-SIM \
  #   --customise_commands "${Z_CMD};${ETA_CMD};${ENERGY_CMD};${SINGLE_PION_CMD};${NO_RELATIVE_TIMING_CMD}" \
  #   --fileout "file:${GENSIM_DIR}/${SAMPLE}_GEN-SIM.root" \
  #   --python_filename "${SAMPLE}_cfg.py" \
  #   -n "${NEVENTS}"

  # cmsDriver.py step1 \
  #   --python_filename "${SAMPLE}_step1_cfg.py" \
  #   --filein "file:${GENSIM_DIR}/${SAMPLE}_GEN-SIM.root" \
  #   --fileout "file:${GENSIM_DIR}/${SAMPLE}_step1_GEN-SIM-RAW.root" \
  #   --pileup NoPileUp \
  #   --customise Configuration/DataProcessing/Utils.addMonitoring \
  #   --eventcontent RAWSIM \
  #   --datatier GEN-SIM-RAW \
  #   --conditions auto:phase1_2024_realistic \
  #   --step DIGI,L1,DIGI2RAW \
  #   --geometry DB:Extended \
  #   --era Run3_2024 \
  #   --mc \
  #   -n "${NEVENTS}"

done

echo "All GEN-SIM and GEN-SIM-RAW files are done."


cd PF-Reco-Analysis
./SetTimingThreshold.sh 5.0
echo "========================================"
echo "  Timing threshold: 5 ns"
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

#mkdir -p "${OUT_DIR}"

echo
echo "-----------------------------------"
echo "Running re-RECO for all samples"
echo "-----------------------------------"
echo

for SAMPLE in "${SAMPLES_TO_PROCESS[@]}"; do
  INFILE="${SAMPLE}_step1_GEN-SIM-RAW.root"
  OUTFILE="pf_only_reReco_${SAMPLE}_standardPF.root"

  echo
  echo "-----------------------------------"
  echo "Running re-RECO: ${SAMPLE}"
  echo "Input:  ${GENSIM_DIR}/${INFILE}"
  echo "Output: ${OUT_DIR}/${OUTFILE}"
  echo "-----------------------------------"
  echo

  cmsDriver.py MyPFStudy_ReReco_MC_Sim \
    --mc \
    --conditions auto:phase1_2025_realistic \
    --step RAW2DIGI,L1Reco,RECO \
    --geometry DB:Extended \
    --era Run3 \
    --filein "file:${GENSIM_DIR}/${INFILE}" \
    --fileout file:pf_only_reReco_MC_Sim.root \
    --eventcontent RECO \
    --datatier RECO \
    --process ReRECO \
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco_MC_Sim.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_particleFlowRecHit*_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*'))" \
    --no_exec \
    -n "${NEVENTS}"
  rm -f pf_only_reReco_MC_Sim.root

  cmsRun MyPFStudy_ReReco_MC_Sim_RAW2DIGI_L1Reco_RECO.py
  mv pf_only_reReco_MC_Sim.root "${OUT_DIR}/${OUTFILE}"

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

for SAMPLE in "${SAMPLES_TO_PROCESS[@]}"; do
  INFILE="${SAMPLE}_step1_GEN-SIM-RAW.root"
  OUTFILE="pf_only_reReco_${SAMPLE}_seedTimingPF.root"

  echo "=== Seed timing PF: ${SAMPLE} ==="

  cmsDriver.py MyPFStudy_ReReco_MC_Sim \
      --mc \
      --conditions auto:phase1_2025_realistic \
      --step RAW2DIGI,L1Reco,RECO \
      --geometry DB:Extended \
      --era Run3 \
      --filein "file:${GENSIM_DIR}/${INFILE}" \
      --fileout file:pf_only_reReco_MC_Sim.root \
      --eventcontent RECO \
      --datatier RECO \
      --process ReRECO \
      --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco_MC_Sim.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_particleFlowRecHit*_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*'))" \
      --no_exec -n ${NEVENTS}

  rm -f pf_only_reReco_MC_Sim.root
  cmsRun MyPFStudy_ReReco_MC_Sim_RAW2DIGI_L1Reco_RECO.py

  mv pf_only_reReco_MC_Sim.root "${OUT_DIR}/${OUTFILE}"

done


echo " "
echo "All tests completed, output files are in:"
echo "${OUT_DIR}"
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
  for SAMPLE in "${SAMPLES_TO_PROCESS[@]}"; do
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

python3 onepion_2dheatmap_modified.py

echo
echo "All generation, re-RECO, and ntupling steps are complete."
echo "Output files are in: ${OUT_DIR}"