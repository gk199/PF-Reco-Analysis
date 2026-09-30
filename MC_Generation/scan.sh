#!/bin/bash
set -e

NEVENTS=100
DR_VALUES=(0.1 0.4 0.8)   # pion-pair opening angle (delta-R)
DT_VALUES=(0.0 1.0 5.0)   # timing offset of pi- relative to pi+ (ns)

for DR in "${DR_VALUES[@]}"; do
  for DT in "${DT_VALUES[@]}"; do

    TAG="DR${DR}_DT${DT}"
    echo "=== DR=${DR}  DT=${DT} ns ==="

    # Convert deltaR -> physical separation at HCAL barrel inner face (R ~ 177 cm).
    # Pions shower hadronically so HCAL sets the relevant angular scale.
    DELTA=$(python3 -c "print(${DR} * 177.0)")

    # Timing: the plugin adds ip*fT ns to each particle's vertex time, so particle 1
    # arrives DT ns later than particle 0.  TMin < TMax is required (strict), so for
    # DT=0 we simply disable UseDeltaT; for DT>0 we bracket tightly around DT.
    if python3 -c "import sys; sys.exit(0 if float('${DT}') == 0 else 1)"; then
      TIMING_CMD="process.generator.PGunParameters.UseDeltaT=cms.bool(False)"
    else
      TMIN=$(python3 -c "print(${DT} - 0.001)")
      TMAX=$(python3 -c "print(${DT} + 0.001)")
      TIMING_CMD="process.generator.PGunParameters.UseDeltaT=cms.bool(True);process.generator.PGunParameters.TMin=cms.double(${TMIN});process.generator.PGunParameters.TMax=cms.double(${TMAX})"
    fi

    # Build the CMSSW config (--no_exec skips cmsRun so we can inspect first)
    cmsDriver.py Configuration/Generator/python/DiPionGun_cfi.py \
      --mc \
      --step GEN,SIM,DIGI,L1,DIGI2RAW \
      --era Run3_2024 \
      --geometry DB:Extended \
      --conditions auto:phase1_2024_realistic \
      --beamspot Realistic25ns13p6TeVEarly2022Collision \
      --pileup NoPileUp \
      --eventcontent RAWSIM \
      --datatier GEN-SIM-RAW \
      --customise_commands "${TIMING_CMD};process.generator.PGunParameters.Delta=cms.double(${DELTA})" \
      --fileout "file:DiPionGun_${TAG}_GEN-SIM-RAW.root" \
      --python_filename "DiPionGun_${TAG}_cfg.py" \
      -n ${NEVENTS} \
      --no_exec

    # Run the generated config, capturing stdout+stderr to a log
    if cmsRun "DiPionGun_${TAG}_cfg.py" >& "log_${TAG}.txt"; then
      echo "    -> DiPionGun_${TAG}_GEN-SIM-RAW.root"
    else
      echo "    FAILED (exit $?) — see log_${TAG}.txt"
    fi

  done
done

echo "All done."
