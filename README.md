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
./SetTimingThreshold.sh <time in ns (5.0)>
./TestAllParticleFlow.sh
./PlotAllParticleFlow.sh
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

cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
```
The output will be `pf_only_reReco*.root` depending which files is run. Each one creates an output file at a different datatier. For data, the options are: RECO, AOD, AODfull (with trigger results). AOD with trigger results can be run through the [DQM plotting framework](https://github.com/gk199/cmssw/blob/PFdevelopment/PF_README.md#monitoring-and-plotting-dqmoffline). For MC, the options are: `MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py` to save PF clusters, cands, and rechits; or `MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py` to also save g4 sim hits and gen particles (the default, as it is compatible with Simon's framework). 

To check the event content, use `edmDumpEventContent` and search for the collection you are interested in:
```
edmDumpEventContent pf_only_reReco_MC_Sim.root | grep EcalRecHits 
edmDumpEventContent pf_only_reReco_MC_Sim.root | grep particleFlow 
```

An option is also prepared to run over the HCAL phase scan files, to save the uMNio word that encodes the phase delay. To run this:
```
cmsRun MyPFStudy_ReReco_RAW2DIGI_L1Reco_RECO_phaseScanFiles.py
```

## Condor Submission (HCAL phase scan files)
1. Create `input_files.txt` of the phase scan files to run over. These are listed in the format: `file:/eos/cms/store/group/dpg_hcal/comm_hcal/QIEPhaseScan2025/JetMET*/*.root`.

2. Test one job:
```
cd Condor
./submit_phaseScan.sh -t input_files.txt /your/output/dir
```

3. Submit all jobs:
```
cd Condor
./submit_phaseScan.sh input_files.txt /your/output/dir
```

4. Montor:
```
condor_q                          # see running/queued jobs
condor_q -better-analyze <jobid>  # diagnose a held/failing job
tail -f Condor/logs/condor.log    # watch the log
```

## MC cmsDriver command
For MC, a slightly different python config is needed (MC specific GlobalTag, MC flag, no pp scenario). The two config generations are given below. For MC, add `outputCommands = cms.untracked.vstring('drop *', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*')` in `process.out` to keep the sim hits and gen particles (for Simon's truth matching studies). 
```
cmsDriver.py MyPFStudy_ReReco_MC_Sim \
    --mc --conditions auto:phase1_2025_realistic \
    --step DIGI,RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --filein file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/SinglePiPt10_step1_GEN-SIM-RAW.root \
    --fileout file:pf_only_reReco_MC_Sim.root \
    --eventcontent RECO --datatier RECO --process ReRECO \
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco_MC_Sim.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*'))" \
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
    inputFiles=file:pf_only_reReco_MC_Sim_standardPF.root \
    outputFile=pfObjectsNtuple_standardPF.root
```
or with the bash script:
```
./NtupleAllParticleFlow.sh
```
The clusters and rechits are looped over in `PFObjectsNtupler.cc`. Edit this to change the matching and what final output variables are filled in the ntuple. This will produce a file called `pfObjectsNtuple.root` that is now used in the plotting / analysis step.

Note that the Global Tag in `runPFObjectsNtupler_cfg.py` may need to be adjusted between data and MC. 

Note: when you run the ntupler on the re-reco output (not directly on RAW), the process label in the input tag will need to match. If needed, you can override taguMNio in `runPFObjectsNtupler_cfg.py` by adding e.g. `taguMNio = cms.untracked.InputTag("hcalDigis", "", "ReRECO")` to the pfObjectsNtupler PSet.

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
