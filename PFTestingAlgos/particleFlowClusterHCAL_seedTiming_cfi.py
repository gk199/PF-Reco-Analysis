# Configuration for Option 2: timing cut only on seeds when linking depth clusters.
#
# Use together with:
#   Basic2DGenericTopoClusterizer.cc  (original, no per-cell timing cut)
#   PFMultiDepthClusterizer_seedTiming.cc.edit  -> seed-to-seed timing at link stage
#
# A depth cluster is stacked only if its seed is within timeThreshold of the
# leading HB depth seed of its linked group, and it is linked to that seed through
# stacked clusters (no holes); the others are stacked among themselves in later passes.  Depth clusters without a usable
# seed time (-999, or outside HB) follow the stack they are linked to.
# Individual cells within an accepted depth cluster are all included regardless
# of their timing.
# Copy this file to:
#   RecoParticleFlow/PFClusterProducer/python/particleFlowClusterHCAL_cfi.py

import FWCore.ParameterSet.Config as cms

_thresholdsHB = cms.vdouble(0.8, 0.8, 0.8, 0.8)
_thresholdsHE = cms.vdouble(0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
_thresholdsHBphase1 = cms.vdouble(0.1, 0.2, 0.3, 0.3)
_thresholdsHEphase1 = cms.vdouble(0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
#updated HB RecHit threshold for 2023
_thresholdsHBphase1_2023 = cms.vdouble(0.4, 0.3, 0.3, 0.3)

particleFlowClusterHCAL = cms.EDProducer('PFMultiDepthClusterProducer',
       clustersSource = cms.InputTag("particleFlowClusterHBHE"),
       usePFThresholdsFromDB = cms.bool(False),
       pfClusterBuilder = cms.PSet(
           algoName = cms.string("PFMultiDepthClusterizer"),
           nSigmaEta = cms.double(2.),
           nSigmaPhi = cms.double(2.),
           # timing parameters: depth-cluster seed time vs the group's leading HB seed
           useTiming = cms.untracked.bool(True),
           timeThreshold = cms.untracked.double(3),  # ns # PFTesting - change this
           #pf clustering parameters
           minFractionToKeep = cms.double(1e-7),
           allCellsPositionCalc = cms.PSet(
               algoName = cms.string("Basic2DGenericPFlowPositionCalc"),
               minFractionInCalc = cms.double(1e-9),
               posCalcNCrystals = cms.int32(-1),
               logWeightDenominatorByDetector = cms.VPSet(
                cms.PSet( detector = cms.string("HCAL_BARREL1"),
                          depths = cms.vint32(1, 2, 3, 4),
                          logWeightDenominator = _thresholdsHB,
                          ),
                cms.PSet( detector = cms.string("HCAL_ENDCAP"),
                          depths = cms.vint32(1, 2, 3, 4, 5, 6, 7),
                          logWeightDenominator = _thresholdsHE,
                          )
                ),
               minAllowedNormalization = cms.double(1e-9)
           )
       ),
       positionReCalc = cms.PSet(),
       energyCorrector = cms.PSet()
)

# offline 2018 -- uncollapsed
from Configuration.Eras.Modifier_run2_HE_2018_cff import run2_HE_2018
from Configuration.ProcessModifiers.run2_HECollapse_2018_cff import run2_HECollapse_2018
(run2_HE_2018 & ~run2_HECollapse_2018).toModify(particleFlowClusterHCAL,
    pfClusterBuilder = dict(
        allCellsPositionCalc = dict(logWeightDenominatorByDetector = {1 : dict(logWeightDenominator = _thresholdsHEphase1) } ),
    ),
)

# offline 2021
from Configuration.Eras.Modifier_run3_HB_cff import run3_HB
run3_HB.toModify(particleFlowClusterHCAL,
    pfClusterBuilder = dict(
        allCellsPositionCalc = dict(logWeightDenominatorByDetector = {0 : dict(logWeightDenominator = _thresholdsHBphase1) } ),
    ),
)

# offline 2023
from Configuration.Eras.Modifier_run3_egamma_2023_cff import run3_egamma_2023
run3_egamma_2023.toModify(particleFlowClusterHCAL,
    pfClusterBuilder = dict(
        allCellsPositionCalc = dict(logWeightDenominatorByDetector = {0 : dict(logWeightDenominator = _thresholdsHBphase1_2023) } ),
    ),
)

# HCALonly WF
particleFlowClusterHCALOnly = particleFlowClusterHCAL.clone(
    clustersSource = "particleFlowClusterHBHEOnly"
)

#--- Use DB conditions for cuts&seeds for Run3 and Phase2
from Configuration.Eras.Modifier_hcalPfCutsFromDB_cff import hcalPfCutsFromDB
hcalPfCutsFromDB.toModify( particleFlowClusterHCAL,
                           usePFThresholdsFromDB = True)
