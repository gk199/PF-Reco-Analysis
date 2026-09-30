
# Appended to the cmsDriver configuration by worker.sh.
import os

# Preserve reference targets in this two-file test. Both old input products
# and newly produced products remain, distinguished by their process names.
# This is intentionally larger than the original PF-only output selection.
process.RECOoutput.outputCommands = cms.untracked.vstring('keep *')
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.options.wantSummary = cms.untracked.bool(True)

process.TFileService = cms.Service(
    'TFileService', fileName=cms.string(os.environ['LLP_NTUPLE_FILE'])
)
process.pfObjectsNtupler = cms.EDAnalyzer(
    'PFObjectsNtupler',
    pfCandidates=cms.InputTag('particleFlow', '', 'RECO'),  # not re-made in PF-cluster-only rereco
    ecalClusters=cms.InputTag('particleFlowClusterECAL', '', 'ReRECO'),
    hcalClusters=cms.InputTag('particleFlowClusterHCAL', '', 'ReRECO'),
    pfBlocks=cms.InputTag('particleFlowBlock', '', 'RECO'),  # not re-made in PF-cluster-only rereco
    hbheRechits=cms.InputTag('hbhereco', '', 'RECO'),
    ecalRechitsEB=cms.InputTag('ecalRecHit', 'EcalRecHitsEB', 'RECO'),
    ecalRechitsEE=cms.InputTag('ecalRecHit', 'EcalRecHitsEE', 'RECO'),
    ecalRechitsES=cms.InputTag('ecalPreshowerRecHit', 'EcalRecHitsES', 'RECO'),
    genParticles=cms.InputTag('genParticles'),
)
process.ntupling_step = cms.EndPath(process.pfObjectsNtupler)
process.schedule.append(process.ntupling_step)
