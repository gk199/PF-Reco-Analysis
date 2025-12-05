# Particle Flow Reconstruction, Ntuples, and Plotting

To setup, inside of the CMSSW release where the PF development is done (such as [here](https://github.com/gk199/cmssw/blob/PFdevelopment/PF_README.md)):
```
git clone git@github.com:gk199/PF-Reco-Analysis.git
cd PF-Reco-Analysis
git branch <your-dev-branch>
```

# Run re-reco
Given a data or MC file at RECO level, re-reco is run to save the PF clusters (ECAL, HCAL), blocks, and candidates. The ECAL (EBEE, saving preshower not working currently) and HCAL (HBHE, HO) raw rechits are also saved. To run, change the path to the data file in `MyPFStudy_ReReco*_RAW2DIGI_L1Reco_RECO.py` and the number of events to run over, and run:
```
cmsenv
voms-proxy-init --rfc --voms cms --valid 172:00

cmsRun MyPFStudy_ReReco*_RAW2DIGI_L1Reco_RECO.py

cmsRun MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py
```
The output will be `pf_only_reReco*.root` depending which files is run. Each one creates an output file at a different datatier. For data, the options are: RECO, AOD, AODfull (with trigger results). AOD with trigger results can be run through the [DQM plotting framework](https://github.com/gk199/cmssw/blob/PFdevelopment/PF_README.md#monitoring-and-plotting-dqmoffline). For MC, the options are: `MyPFStudy_ReReco_MC_RAW2DIGI_L1Reco_RECO.py` to save PF clusters, cands, and rechits; or `MyPFStudy_ReReco_MC_Sim_DIGI_RAW2DIGI_L1Reco_RECO.py` to also save g4 sim hits and gen particles (compatible with Simon's framework). 

To check the event content, use `edmDumpEventContent` and search for the collection you are interested in:
```
edmDumpEventContent pf_only_reReco_MC_Sim.root | grep EcalRecHits 
edmDumpEventContent pf_only_reReco_MC_Sim.root | grep particleFlow 
```

## MC cmsDriver command
For MC, a slightly different python config is needed (MC specific GlobalTag, MC flag, no pp scenario). The two config generations are given below. For MC, add `outputCommands = cms.untracked.vstring('drop *', 'keep *_g4SimHits_*_*', 'keep *_genParticles_*_*') in `process.out` to keep the sim hits and gen particles (for Simon's truth matching studies). 
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
    --customise_commands="process.RECOoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reReco.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*'))" \
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
    --customise_commands="process.AODoutput = cms.OutputModule('PoolOutputModule', fileName = cms.untracked.string('pf_only_reRecoAOD.root'), outputCommands = cms.untracked.vstring('drop *', 'keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*'))" \
    --no_exec -n 100
```
And saving all the trigger information:
```
cmsDriver.py MyPFStudy_ReRecoAODfull \
    --data --conditions 150X_dataRun3_Prompt_v1 \
    --step RAW2DIGI,L1Reco,RECO --geometry DB \
    --era Run3 --scenario pp --filein root://cms-xrd-global.cern.ch//store/data/Run2025E/Muon0/RAW-RECO/MUOJME-PromptReco-v1/000/395/982/00000/01c7900e-0585-4df0-8f2e-23ba45358ed8.root --fileout file:pf_only_reRecoAODfull.root \
    --eventcontent AOD --datatier AOD --process ReRECOtoAOD \
    --customise_commands="process.AODoutput.outputCommands.extend(['keep *_particleFlowClusterECAL_*_*', 'keep *_particleFlowClusterHCAL_*_*', 'keep *_particleFlowBlock_*_*', 'keep *_particleFlow_*_*', 'keep *_hbhereco_*_*', 'keep *_horeco_*_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*', 'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*', 'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*'])" \
    --no_exec -n 100 
```

# Ntuple Production
Change the root file listed in the plotting script as needed, and then run: 
```
scram b -j 8
cmsRun PFObjectsNtupler/python/runPFObjectsNtupler_cfg.py
```
The clusters and rechits are looped over in `PFObjectsNtupler.cc`. Edit this to change the matching and what final output variables are filled in the ntuple. This will produce a file called `pfObjectsNtuple.root` that is now used in the plotting / analysis step.

# Plotting
```
python3 Plotting/plot_pfObjects.py
```
Outputs a root file with the histograms.

# Event Display
The event display is done from the `pfObjectsNtuple.root`.
```
cd EventDisplay
python3 PF_HCAL_rechit_cluster.py
```
This shows 4 depths of HB, with the cluster outlined in red, the top four plots showing HCAL rechits energy in each depth, and the lower four plots showing HCAL rechits time. 

Run from a conda virtual environment (I run this locally on my laptop to interact with the GUI) with `package-list.txt`. This is way more requirements than actually needed, just found from running `conda list --export` in my area. 

ToDo: currently using MAHI time, would like to plot from TDC time! 

# Wish List
- Setup re-reco to run with CRAB jobs for re-processing CMS datasets
- Setup re-reco step to run with Condor jobs for re-processing of MC data files
- Add ECAL rechits to PFObjectsNtupler.cc, with matching to clusters
- Extend plotting to include rechits, clusters, blocks, and PF candidates
- Add HCAL PF cuts to rechits kept for analysis 
- Write event display code for clusters and rechits: DONE
