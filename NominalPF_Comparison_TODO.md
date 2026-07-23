# QCD-like Background: Nominal PF vs. Timing PF — TODO

Goal: quantify how much the seed-level timing cut (Option 2, `seedTiming`) disturbs an in-time background sample — cluster multiplicity and cluster energy spectrum, nominal vs. one or more timing thresholds, as ratio plots.

**Sample: reuse the HCAL phase-scan data (`laserType == 0`).** RelVal QCD GEN-SIM-DIGI-RAW isn't on disk for this release cycle (RelVal RAW/DIGI tiers are moved to tape quickly). Rather than sourcing a fresh JetMET dataset via DAS, the phase-scan input files (`Condor/input_files.txt`, real JetMET0 collision data from run 392175, already sitting on local EOS) are a QIE-phase calibration scan over the *same* data — `laserType` (the uMNio phase-delay branch, see [Plotting/plot_phaseScan_timing.py:47](Plotting/plot_phaseScan_timing.py#L47)) lets us select the nominal in-time phase (`laserType == 0`) as the background baseline, with no new files needed at all.

**Important caveat resolved:** the *existing* `pf_only_reReco_phaseScan_job*.root`/ `pfObjectsNtuple_phaseScan.root` outputs were produced by a Condor batch submitted mid-April in the shared CMSSW area, during the same window the pion-gun timing-variant clusterizers were being swapped in and out for other tests — so we can't fully trust that as "definitely nominal PF". Decision: **don't reuse that archive**; produce a fresh nominal pass now (today's build is = original clusterizer + timingFix position-calc, confirmed by diffing against `PFTestingAlgos/*_original.cc.edit`), using the existing `Condor/input_files_first100.txt` subset so nominal and timing-PF outputs come from an identical, known file list.

`Plotting/compare_qcd_ratio.py` supports `--laser-type <int>` to select only events with that `laserType` value when filling histograms (0 = nominal phase) — added and smoke-tested (confirmed it still passes 100/100 events when given `--laser-type -1000` against ntuples where every event defaults to that sentinel, i.e. no filtering regression).

## 1. Fresh nominal-PF pass over the first-100 phase-scan files
- [x] Confirm current build is still original+timingFix (it was as of writing — re-check
      if anything else has been tested in between: `diff RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit`)
- [x] `scram b -j8` if anything needed rebuilding
```bash
cd PF-Reco-Analysis
mkdir -p /eos/user/g/gkopp/PF_PhaseScan/nominal_fresh
cd Condor
./submit_phaseScan.sh -t input_files_first100.txt /eos/user/g/gkopp/PF_PhaseScan/nominal_fresh   # test 1 job first
./submit_phaseScan.sh input_files_first100.txt /eos/user/g/gkopp/PF_PhaseScan/nominal_fresh       # then the full 100
```
- [x] `condor_q` / `tail -f Condor/logs/condor.log` to monitor
- [x] Confirm all 100 jobs completed (`ls /eos/user/g/gkopp/PF_PhaseScan/nominal_fresh | wc -l`)

## 2. Timing-PF pass (seed-level, Option 2) — once per threshold
- [ ] Decide threshold list to scan, e.g. `1 2 3 5` ns (matches values used elsewhere). 3ns is used since that is consistent with the recent dipion tests. 

For each threshold `<ns>`:
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis
./SetTimingThreshold.sh <ns>

cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit \
   RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit \
   RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py \
   RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py
scram b -j8
cd PF-Reco-Analysis

mkdir -p /eos/user/g/gkopp/PF_PhaseScan/timing_<ns>ns
cd Condor
./submit_phaseScan.sh -t input_files_first100.txt /eos/user/g/gkopp/PF_PhaseScan/timing_<ns>ns   # test 1 job first
./submit_phaseScan.sh input_files_first100.txt /eos/user/g/gkopp/PF_PhaseScan/timing_<ns>ns
```
- [ ] Run for each threshold in the list above (same 100 input files each time — only the clusterizer build changes)
- [ ] After the last threshold, revert to the original clusterizer files (see README "Reverting back to the original ones") so the working tree is clean for other studies

## 3. Ntuple each batch
Adapt the `NtuplePhaseScan.sh` pattern (it globs a directory of
`pf_only_reReco_phaseScan_job*.root` and runs them through the ntupler in one
`cmsRun` call) to point at each new output directory:
```bash
cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles="$(ls /eos/user/g/gkopp/PF_PhaseScan/nominal_fresh/pf_only_reReco_phaseScan_job*.root | sed 's|/eos/|root://eosuser.cern.ch//eos/|' | paste -sd,)" \
    outputFile=/tmp/pfObjectsNtuple_bkg_nominal.root
xrdcp -f /tmp/pfObjectsNtuple_bkg_nominal.root root://eosuser.cern.ch//eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_bkg_nominal.root

# repeat per threshold, pointing at timing_<ns>ns/ and naming pfObjectsNtuple_bkg_timing_<ns>ns.root
```
- [ ] Run for nominal + each threshold

## 4. Plot — select laserType == 0 (nominal in-time phase) from both samples
```bash
python3 Plotting/compare_qcd_ratio.py \
    --nominal /eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_bkg_nominal.root \
    --timing  /eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_bkg_timing_1ns.root \
              /eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_bkg_timing_2ns.root ... \
    --labels  1ns 2ns ... \
    --laser-type 0 \
    --output  qcd_timing_ratio.pdf \
    --rootfile qcd_timing_ratio.root
```
Produces (each as an overlay + ratio-panel page): N clusters/event, ΣE/event,
per-cluster energy spectrum, rechits/cluster — computed only over
`laserType == 0` events in each file. Console output also prints
cluster-count-as-%-of-nominal per threshold.

## Open decisions before starting
- [ ] Threshold values to scan
- [ ] Whether 100 files / laserType==0 subset gives enough statistics, or whether to widen to the full `input_files.txt` list once the pipeline is validated on the smaller batch
