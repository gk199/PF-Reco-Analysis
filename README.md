# Particle Flow Reconstruction, Ntuples, and Plotting

To setup, inside of the CMSSW release where the PF development is done (such as [here](https://github.com/gk199/cmssw/blob/PFdevelopment/PF_README.md)):
```
git clone git@github.com:gk199/PF-Reco-Analysis.git
cd PF-Reco-Analysis
git branch <your-dev-branch>
```

To get both this setup as well as the DQM PF monitoring area: 
```
cmsrel CMSSW_15_0_6
cd CMSSW_15_0_6/src
cmsenv
git cms-addpkg HLTrigger/Configuration
git cms-addpkg RecoParticleFlow/PFClusterProducer
git clone git@github.com:gk199/PF-Reco-Analysis.git
scram b -j 8

git remote add origin git@github.com:gk199/cmssw.git
git fetch origin PFdevelopment
git sparse-checkout list
git checkout origin/PFdevelopment -- DQMOffline/ParticleFlow PF_README.md
```

# Testing new PF Algorithms
## Bash Script for New Algos and Re-Reco
```
cmsenv
./SetTimingThreshold.sh <time in ns, such as (5.0)>
./TestAllParticleFlow.sh
./NtupleAllParticleFlow.sh
./PlotAllParticleFlow.sh
mv hcal_comparison.pdf hcal_comparison_<ns>ns.pdf

# all contained in:
cmsenv
./ScanTimingThresholds.sh 0 2 3 5 
```
Since most of the timing checks are done at the depth cluster stacking level, the timing cut can actually be scanned by just changing a parameter in the python file (no re-compiling needed):
```
emacs PFTestingAlgos/particleFlowClusterHCAL_*_cfi.py
```
Edit the line `timeThreshold = cms.untracked.double(3.0),  # ns`. Make sure this is a smaller value than is listed in `Basic2DGenericTopoClusterizer_timing.cc.edit` in the line `_timeThreshold(conf.getUntrackedParameter<double>("timeThreshold", 5.0)) {}`, otherwise this should also be edited. 

## Full Details
The timing algorithm also places a "similar in time" constraint in the gathering step. The threshold can be further optimized as well. To use these files in the re-reco, copy the modified files and remember to recompile.

Option 1 (`timing`): cell level timing cut vs global highest energy seed:
```
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_timing_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8
cd PF-Reco-Analysis
```

Option 2 (`seedTiming`): seed level timing cut vs global highest energy seed:
```
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_seedTiming.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_seedTiming_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8
cd PF-Reco-Analysis
```

Option 3 (`depth1Timing`): seed level timing cut vs first depth cluster seed:
```
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_depth1Timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_timing.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_depth1Timing_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8
cd PF-Reco-Analysis
```

Reverting back to the original ones:
```
cp PF-Reco-Analysis/PFTestingAlgos/Basic2DGenericTopoClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/Basic2DGenericTopoClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterizer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterizer.cc
cp PF-Reco-Analysis/PFTestingAlgos/PFMultiDepthClusterProducer_original.cc.edit RecoParticleFlow/PFClusterProducer/plugins/PFMultiDepthClusterProducer.cc
cp PF-Reco-Analysis/PFTestingAlgos/particleFlowClusterHCAL_original_cfi.py RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

scram b -j 8
```

# Run re-reco
Given a data or MC file at RECO level, re-reco is run to save the PF clusters (ECAL, HCAL), blocks, and candidates. The ECAL (EBEE, saving preshower not working currently) and HCAL (HBHE, HO) raw rechits are also saved. To run, change the path to the data file in `MyPFStudy_ReReco*_RAW2DIGI_L1Reco_RECO.py` and the number of events to run over, and run:
```
cmsenv
voms-proxy-init --rfc --voms cms --valid 172:00

cmsRun MyPFStudy_ReReco*_RAW2DIGI_L1Reco_RECO.py
```
The output will be `pf_only_reReco*.root` depending which file is run. Each one creates an output file at a different datatier. For data, the options are: RECO, AOD, AODfull (with trigger results). AOD with trigger results can be run through the [DQM plotting framework](https://github.com/gk199/cmssw/blob/PFdevelopment/PF_README.md#monitoring-and-plotting-dqmoffline). For MC, use `MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py` to save PF clusters, cands, rechits, g4 sim hits, and gen particles. **Do not use** `MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py` — it re-runs the DIGI step from GEN-SIM-RAW SimHits with potentially mismatched conditions, inflating HCAL cluster counts.

To check the event content, use `edmDumpEventContent` and search for the collection you are interested in:
```
edmDumpEventContent pf_only_reReco_MC.root | grep EcalRecHits 
edmDumpEventContent pf_only_reReco_MC.root | grep particleFlow 
```

An option is also prepared to run over the HCAL phase scan files, to save the uMNio word that encodes the phase delay. To run this:
```
cmsRun MyPFStudy_ReReco_RAW2DIGI_L1Reco_RECO_phaseScanFiles.py
```

## Condor Submission (HCAL phase scan files)
1. Create `input_files.txt` of the phase scan files to run over. These are listed in the format: `file:/eos/cms/store/group/dpg_hcal/comm_hcal/QIEPhaseScan2025/JetMET*/*.root`.

2. Test one job:
```
cmsenv
cd Condor
./submit_phaseScan.sh -t input_files.txt /your/output/dir
```

3. Submit all jobs:
```
cmsenv
cd Condor
./submit_phaseScan.sh input_files.txt /your/output/dir
```

4. Montor:
```
condor_q                          # see running/queued jobs
condor_q -better-analyze <jobid>  # diagnose a held/failing job
tail -f Condor/logs/condor.log    # watch the log
```

## DiPion Gun MC samples

DiPion gun GEN-SIM-RAW files (DR = angular separation, DT = time offset in ns) live in `/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src/`. Generated with `Run3` / `auto:phase1_2024_realistic`. The RECO config uses `auto:phase1_2025_realistic` — the 2025 conditions have tighter HCAL noise thresholds that suppress noise clusters present in 2024-simulated MC, giving a cleaner baseline for PF algorithm comparisons. The RECO config does **not** re-run the DIGI step. A dedicated config accepts `inputFiles` and `outputFile` as command-line arguments:
```
cmsRun MyPFStudy_ReReco_DiPionGun_DIGI_RAW2DIGI_L1Reco_RECO.py \
    inputFiles=file:/path/to/DiPionGun_DR0.1_DT0.0_GEN-SIM-RAW.root \
    outputFile=pf_only_reReco_DiPionGun_DR0.1_DT0.0.root
```

To run RECO + ntupling over all 9 DR/DT combinations in one go:
```
cmsenv
./RunDiPionGunSamples.sh 0 <ns>
```
Outputs are labeled `pf_only_reReco_DiPionGun_DR<X>_DT<Y>.root` and `pfObjectsNtuple_DiPionGun_DR<X>_DT<Y>.root`.

To plot:
```
./PlotDiPionGunSamples.sh 0 <3 5> 
```

## MC cmsDriver command
For MC, a slightly different python config is needed (MC specific GlobalTag, MC flag, no pp scenario). The config below produces `MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py`, starting from RAW2DIGI (no DIGI step) and keeping sim hits and gen particles for truth matching. Use `auto:phase1_2025_realistic` even for 2024-produced MC — the 2025 conditions have tighter HCAL noise thresholds that prevent spurious noise clusters.
```
cmsDriver.py MyPFStudy_ReReco_MC \
    --mc --conditions auto:phase1_2025_realistic \
    --step RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --filein file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/SinglePiPt10_step1_GEN-SIM-RAW.root \
    --fileout file:pf_only_reReco_MC.root \
    --eventcontent RECO --datatier RECO --process ReRECO \
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco_MC.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*'))" \
    --no_exec -n 100
```

## Data cmsDriver command
This keeps a reduced set of output collections: HCAL and ECAL rechits, PF clusters, blocks, and candidates.
```
cmsDriver.py MyPFStudy_ReReco \
    --data --conditions 150X_dataRun3_Prompt_v1 \
    --step RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --scenario pp --filein root://cms-xrd-global.cern.ch//store/data/Run2025E/Muon0/RAW-RECO/MUOJME-PromptReco-v1/000/395/982/00000/01c7900e-0585-4df0-8f2e-23ba45358ed8.root \
    --fileout file:pf_only_reReco.root \
    --eventcontent RECO --datatier RECO --process ReRECO \
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep HcalUMNioDigi_hcalDigis_*_*'))" \
    --no_exec -n 100
```

At AOD level (useful for putting back into the DQM plotting code):
```
cmsDriver.py MyPFStudy_ReRecoAOD \
    --data --conditions 150X_dataRun3_Prompt_v1 \
    --step RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --scenario pp --filein root://cms-xrd-global.cern.ch//store/data/Run2025E/Muon0/RAW-RECO/MUOJME-PromptReco-v1/000/395/982/00000/01c7900e-0585-4df0-8f2e-23ba45358ed8.root \
    --fileout file:pf_only_reRecoAOD.root \
    --eventcontent AOD --datatier AOD --process ReRECOtoAOD \
    --customise_commands="process.AODoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reRecoAOD.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep HcalUMNioDigi_hcalDigis_*_*'))" \
    --no_exec -n 100
```
And saving all the trigger information:
```
cmsDriver.py MyPFStudy_ReRecoAODfull \
    --data --conditions 150X_dataRun3_Prompt_v1 \
    --step RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --scenario pp --filein root://cms-xrd-global.cern.ch//store/data/Run2025E/Muon0/RAW-RECO/MUOJME-PromptReco-v1/000/395/982/00000/01c7900e-0585-4df0-8f2e-23ba45358ed8.root --fileout file:pf_only_reRecoAODfull.root \
    --eventcontent AOD --datatier AOD --process ReRECOtoAOD \
    --customise_commands="process.AODoutput.outputCommands.extend(['keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep HcalUMNioDigi_hcalDigis_*_*'])" \
    --no_exec -n 100 
```

# Ntuple Production
Change the root file listed in the plotting script as needed, and then run: 
```
scram b -j 8
cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_standardPF.root \
    outputFile=pfObjectsNtuple_standardPF.root
```
or with the bash script:
```
./NtupleAllParticleFlow.sh
```
The clusters and rechits are looped over in `PFObjectsNtupler.cc`. Edit this to change the matching and what final output variables are filled in the ntuple. This will produce a file called `pfObjectsNtuple.root` that is now used in the plotting / analysis step.

Note that the Global Tag in `runPFObjectsNtupler_cfg.py` may need to be adjusted between data and MC. 

Note: when you run the ntupler on the re-reco output (not directly on RAW), the process label in the input tag will need to match. If needed, you can override taguMNio in `runPFObjectsNtupler_cfg.py` by adding e.g. `taguMNio = cms.untracked.InputTag("hcalDigis", "", "ReRECO")` to the pfObjectsNtupler PSet.

## Phase Scan Ntupling
To ntuple all phase scan re-reco output files from EOS into a single ntuple:
```
cmsenv
proxy
./NtuplePhaseScan.sh
```
This discovers all `pf_only_reReco_phaseScan_job*.root` files under `/eos/user/g/gkopp/PF_PhaseScan/` via xrootd and runs them through the standard ntupler in one `cmsRun` call, writing the output ntuple directly to `/eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_phaseScan.root` to avoid filling AFS quota. The `laserType` branch (from the HCAL uMNio digi, `valueUserWord(1)`) encodes the phase delay for each event.

# Plotting
```
python3 Plotting/plot_pfObjects.py
```
Outputs a root file with the histograms. This will save the PF candidate pT, eta, and ECAL and HCAL cluster energy.

To compare the clusters across different algorithms:
```
python3 Plotting/compare_hcal_clusters.py --inputdir . --output hcal_comparison.root --pdf hcal_comparison.pdf
```
This will plot the number of HCAL clusters per event, total HCAL cluster energy, and HBHE hits per cluster.

With the bash script:
```
./PlotAllParticleFlow.sh
```

## Testing the HCAL Cluster Timing Fix

`Basic2DGenericPFlowPositionCalc.cc` was edited to skip rechits with `time == -999` (the MAHI invalid sentinel) when computing the E²-weighted cluster time, and to set cluster time to `-999` when all constituent rechits are invalid. This prevents mixed-validity clusters from reporting a bogus average pulled into the −500 to −10 ns range.

Two validation workflows are documented below: a quick pion gun check and a higher-statistics physics sample check via Condor.

`Plotting/compare_timing_fix.py` is used for both. It produces:
- `hcal_time` full range (log y) — shows the population at exactly −999 and any contaminated clusters
- `hcal_time` contaminated range (−998 to −10 ns) — the key diagnostic; should be empty after the fix
- `hcal_time` zoom (valid clusters, −20 to +20 ns) — shows clusters migrating into the physics range after fix
- rechit invalid fraction per cluster and rechits per cluster — informational, identical before/after
- 2D cluster time vs E²-weighted mean valid-rechit time — off-diagonal before fix, on-diagonal after
- Cluster energy and η for contaminated-timing clusters — shows what energy scale and detector region are affected
- 2D cluster time (valid) vs cluster energy — checks no energy-dependent bias is introduced

### Pion Gun Validation (SinglePiPt10, interactive)

Quick check using the SinglePiPt10 GEN-SIM-RAW sample run through `MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py` (100 events, no pileup). Because there is only one pion per event, HCAL occupancy is low and ~80% of clusters have all-invalid rechit timing; the fix is visible as contaminated clusters (−500 to −10 ns) disappearing.

**1. Build the fix**
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src
scram b -j8
cd PF-Reco-Analysis
```

**2. Produce the "before" RECO output** (revert fix temporarily, or use a saved pre-fix file)
```bash
# Revert Basic2DGenericPFlowPositionCalc.cc, then scram b, then:
cmsRun MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC.root pf_only_reReco_MC_standardPF_oldClusterTiming.root
# Restore the fix and scram b again before step 3.
```

**3. Produce the "after" RECO output**
```bash
cmsRun MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py
mv pf_only_reReco_MC.root pf_only_reReco_MC_standardPF.root
```

**4. Run the ntupler on both**
```bash
cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_standardPF_oldClusterTiming.root \
    outputFile=pfObjectsNtuple_standardPF_before.root

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_only_reReco_MC_standardPF.root \
    outputFile=pfObjectsNtuple_standardPF_after.root
```

**5. Plot**
```bash
python3 Plotting/compare_timing_fix.py \
    --before pfObjectsNtuple_standardPF_before.root \
    --after  pfObjectsNtuple_standardPF_after.root \
    --output timing_fix_comparison.pdf
```

### Physics Sample Validation (TTbar or QCD RelVal, interactive)

Physics samples have higher HCAL occupancy: more clusters per event, more clusters with mixed-validity rechits, and realistic jet activity. These give better statistics on the contaminated-timing population than the single-pion gun.

`MyPFStudy_ReReco_MC_condor.py` is a VarParsing-enabled re-reco config that accepts `inputFiles`, `outputFile`, and `maxEvents` on the command line — no file editing needed.

**1. Find GEN-SIM-DIGI-RAW files on DAS**

RelVal samples are the most reliable since they are generated for each CMSSW release:
```bash
voms-proxy-init --rfc --voms cms --valid 172:00
dasgoclient --query="dataset=/RelValTTbar_14TeV/CMSSW_15_0_6*/GEN-SIM-DIGI-RAW"
dasgoclient --query="dataset=/RelValQCD_Pt15To7000_Flat_14TeV/CMSSW_15_0_6*/GEN-SIM-DIGI-RAW"
```

Pick one dataset from the output, then get a file path:
```bash
dasgoclient --query="file dataset=<dataset_from_above>" | head -3
```

**2. Run re-reco and ntupler — before the fix**

Revert `Basic2DGenericPFlowPositionCalc.cc` and `scram b`, then:
```bash
cmsRun MyPFStudy_ReReco_MC_condor.py \
    inputFiles="root://cms-xrd-global.cern.ch//<path_from_DAS>" \
    outputFile=pf_reco_physics_oldClusterTiming.root \
    maxEvents=500

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_reco_physics_oldClusterTiming.root \
    outputFile=pfObjectsNtuple_physics_before.root
```

**3. Run re-reco and ntupler — after the fix**

Restore the fix and `scram b`, then run with the same input file:
```bash
cmsRun MyPFStudy_ReReco_MC_condor.py \
    inputFiles="root://cms-xrd-global.cern.ch//<path_from_DAS>" \
    outputFile=pf_reco_physics.root \
    maxEvents=500

cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py \
    inputFiles=file:pf_reco_physics.root \
    outputFile=pfObjectsNtuple_physics_after.root
```

**4. Plot**
```bash
python3 Plotting/compare_timing_fix.py \
    --before pfObjectsNtuple_physics_before.root \
    --after  pfObjectsNtuple_physics_after.root \
    --output timing_fix_comparison_physics.pdf
```

The contaminated-timing plots (cluster time in −998 to −10 ns, and the energy/η distributions of those clusters) show the fix working across all cluster energies and η values with no detector-region-specific bias.

## Phase Scan Timing Plots
To plot PF cluster time vs phase delay (`laserType`) from the phase scan ntuple:
```
python3 Plotting/plot_phaseScan_timing.py \
    --input /eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_phaseScan.root \
    --pdf phaseScan_timing.pdf \
    --output phaseScan_timing.root
```
Produces 2D histograms and mean-value profiles (TProfile) of HCAL cluster time, HBHE rechit time, and ECAL cluster time vs `laserType`, plus a combined overlay of all three profiles. An energy threshold (`--emin`, default 1 GeV) is applied to suppress noise hits.

# Event Display
The event display is done from the `pfObjectsNtuple.root`.
```
cd EventDisplay
python3 PF_HCAL_rechit_cluster.py
```
This shows 4 depths of HB, with the cluster outlined in red, the top four plots showing HCAL rechits energy in each depth, and the lower four plots showing HCAL rechits time. 

Run from a conda virtual environment (I run this locally on my laptop to interact with the GUI) with `package-list.txt`. This is way more requirements than actually needed, just found from running `conda list --export` in my area. 

ToDo: currently using MAHI time, would like to plot from TDC time! 

## Files from SPVCNN 
The instructions are [here](https://github.com/wpmccormack/spvcnn_instructions/tree/main), however, I am still in the process of adopting this to the newer CMSSW release where current PF studies are done. 

Follow the instructions in the SPVCNN area to put this in `CMSSW_13_0_X`. Get the files with: 
```
git diff --name-only CMSSW_13_3_0_pre5...wpmccormack/add_spvcnn 
git diff CMSSW_13_3_0_pre5...wpmccormack/add_spvcnn <file_path>
```
and try to combine with the `15_0_6` area. 

The files added are:
```
new files:
HLTrigger/Configuration/python/HLT_GRun_SONIC_cff.py (copied over)
RecoParticleFlow/PFClusterProducer/plugins/PFClusterSonicProducer.cc (copied over, added PF cut from DB)

RecoParticleFlow/PFClusterProducer/plugins/PFTruthClusterProducer2.cc (causes many errors, unclear if truth is needed now)
RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_TRUTH2_thresholds2_cfi.py

edited files:
RecoParticleFlow/PFClusterProducer/BuildFile.xml (changes copied over)
RecoParticleFlow/PFClusterProducer/plugins/BuildFile.xml (changes copied over)
RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py (copied over, merged with existing code)

RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_original_cfi.py
```

# Wish List
- Setup re-reco to run with CRAB jobs for re-processing CMS datasets
- Setup re-reco step to run with Condor jobs for re-processing of MC data files
- Add ECAL rechits to PFObjectsNtupler.cc, with matching to clusters
- Extend plotting to include rechits, clusters, blocks, and PF candidates
- Add HCAL PF cuts to rechits kept for analysis: DONE, and also moved to using PF rechits instead of HCAL rechits
- Write event display code for clusters and rechits: DONE
