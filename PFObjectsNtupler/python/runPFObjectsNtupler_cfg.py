import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/pf_only_reReco.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/pf_only_reRecoAOD.root'
        # 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/pf_only_reRecoAODfull.root'
    )
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('pfObjectsNtuple.root')
)

process.pfObjectsNtupler = cms.EDAnalyzer("PFObjectsNtupler",
    pfCandidates = cms.InputTag("particleFlow"),
    ecalClusters = cms.InputTag("particleFlowClusterECAL"),
    hcalClusters = cms.InputTag("particleFlowClusterHCAL"),
    pfBlocks = cms.InputTag("particleFlowBlock")
)

process.p = cms.Path(process.pfObjectsNtupler)