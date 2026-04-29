# DIGI -> RAW2DIGI -> L1Reco -> RECO for DiPionGun GEN-SIM-RAW samples.
# Era and GlobalTag match the SinglePiPt10 generation (Run3 / phase1_2025_realistic).
#
# Usage:
#   cmsRun MyPFStudy_ReReco_DiPionGun_DIGI_RAW2DIGI_L1Reco_RECO.py \
#       inputFiles=file:/path/to/DiPionGun_DR0.1_DT0.0_GEN-SIM-RAW.root \
#       outputFile=pf_only_reReco_DiPionGun_DR0.1_DT0.0.root

import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing('analysis')
options.inputFiles = 'file:/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src/DiPionGun_DR0.1_DT0.0_GEN-SIM-RAW.root'
options.outputFile = 'pf_only_reReco_DiPionGun_DR0.1_DT0.0.root'
options.parseArguments()
# maxEvents is set directly on process, not via VarParsing, to avoid VarParsing
# appending _numEventN to the output filename.

from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('ReRECO', Run3)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100),
    output = cms.optional.untracked.allowed(cms.int32, cms.PSet)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(options.inputFiles),
    inputCommands = cms.untracked.vstring(
        'keep *',
        'drop *_genParticles_*_*',
        'drop *_genParticlesForJets_*_*',
        'drop *_kt4GenJets_*_*',
        'drop *_kt6GenJets_*_*',
        'drop *_iterativeCone5GenJets_*_*',
        'drop *_ak4GenJets_*_*',
        'drop *_ak7GenJets_*_*',
        'drop *_ak8GenJets_*_*',
        'drop *_ak4GenJetsNoNu_*_*',
        'drop *_ak8GenJetsNoNu_*_*',
        'drop *_genCandidatesForMET_*_*',
        'drop *_genParticlesForMETAllVisible_*_*',
        'drop *_genMetCalo_*_*',
        'drop *_genMetCaloAndNonPrompt_*_*',
        'drop *_genMetTrue_*_*',
        'drop *_genMetIC5GenJs_*_*'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    TryToContinue = cms.untracked.vstring(),
    accelerators = cms.untracked.vstring('*'),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    deleteNonConsumedUnscheduledModules = cms.untracked.bool(True),
    dumpOptions = cms.untracked.bool(False),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(0)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    holdsReferencesToDeleteEarly = cms.untracked.VPSet(),
    makeTriggerResults = cms.obsolete.untracked.bool,
    modulesToCallForTryToContinue = cms.untracked.vstring(),
    modulesToIgnoreForDeleteEarly = cms.untracked.vstring(),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(0),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)

# GlobalTag matching the Run3 / phase1_2025_realistic conditions (same as SinglePiPt10)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2025_realistic', '')

# Output: slim RECO keeping only PF and calorimeter objects
process.RECOoutput = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('file:' + options.outputFile),
    outputCommands = cms.untracked.vstring(
        'drop *',
        'keep *_particleFlowClusterECAL_*_*',
        'keep *_particleFlowClusterHCAL_*_*',
        'keep *_particleFlowBlock_*_*',
        'keep *_particleFlow_*_*',
        'keep *_hbhereco_*_*',
        'keep *_horeco_*_*',
        'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_*',
        'keep EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_*',
        'keep EcalRecHitsSorted_ecalPreshowerRecHit_EcalRecHitsES_*',
        'keep *_g4SimHits_*_*',
        'keep *_genParticles_*_*',
    )
)

# Path and EndPath definitions
process.raw2digi_step      = cms.Path(process.RawToDigi)
process.L1Reco_step        = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step        = cms.EndPath(process.endOfProcess)
process.RECOoutput_step    = cms.EndPath(process.RECOoutput)

process.schedule = cms.Schedule(
    process.raw2digi_step,
    process.L1Reco_step,
    process.reconstruction_step,
    process.endjob_step,
    process.RECOoutput_step,
)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
process = customiseLogErrorHarvesterUsingOutputCommands(process)

from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
