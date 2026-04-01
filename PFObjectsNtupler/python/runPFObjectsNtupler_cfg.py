import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/ttbar/pf_only_reReco_MC_Sim_new_timing.root'
        #'file:/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/SinglePiPt100_1000_n1000/pf_only_reReco_MC_Sim_new_timing.root'
        #'file:/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/SinglePi0E100_1000_n5000/pf_only_reReco_MC_Sim_new_timing.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reRecoAOD.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reRecoAODfull.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('/eos/user/c/chtong/Public/Rereco/ttbar_rereco/pfObjectsNtuple_new_timing.root')
    #fileName = cms.string('/eos/user/c/chtong/Public/Rereco/SinglePiPt100_1000_n1000/pfObjectsNtuple_new_timing.root')                               
    #fileName = cms.string('/eos/user/c/chtong/Public/Rereco/SinglePi0E100_1000_n5000/pfObjectsNtuple_new_timing.root')                              
    #fileName = cms.string('/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/ttbar/pfObjectsNtuple_new.root')
    #fileName = cms.string('/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/SinglePiPt100_1000_n1000/pfObjectsNtuple_new.root')
    #fileName = cms.string('/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root')
)

process.pfObjectsNtupler = cms.EDAnalyzer("PFObjectsNtupler",
    pfCandidates = cms.InputTag("particleFlow"),
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    pfBlocks = cms.InputTag("particleFlowBlock"),
    # adding rechits now
    hbheRechits = cms.InputTag("hbhereco"), 
    ecalRechitsEB = cms.InputTag("ecalRecHit", "EcalRecHitsEB"), 
    ecalRechitsEE = cms.InputTag("ecalRecHit", "EcalRecHitsEE"),      
    ecalRechitsES = cms.InputTag("ecalPreshowerRecHit", "EcalRecHitsES")
)

process.p = cms.Path(process.pfObjectsNtupler)