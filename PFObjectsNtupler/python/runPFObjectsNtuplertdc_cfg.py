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
        'file:/afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reReco_MC.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reRecoAOD.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reRecoAODfull.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pfObjectsNtuple_tdc.root')
)

process.pfObjectsNtuplertdc = cms.EDAnalyzer("PFObjectsNtuplertdc",
    pfCandidates = cms.InputTag("particleFlow"),
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    pfBlocks = cms.InputTag("particleFlowBlock"),
    # adding rechits now
    hbheRechits = cms.InputTag("hbhereco"), 
)

process.p = cms.Path(process.pfObjectsNtuplertdc)