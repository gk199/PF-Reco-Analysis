import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.inputFiles = 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/pf_only_reReco_MC_Sim.root'
options.outputFile = 'pfObjectsNtuple.root'
options.parseArguments()

process = cms.Process("NTUPLE")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressWarning = cms.untracked.vstring("pfObjectsNtupler")

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2025_realistic', '')
process.CaloGeometryBuilder.SelectedCalos = [x for x in process.CaloGeometryBuilder.SelectedCalos if x != 'CASTOR']

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
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