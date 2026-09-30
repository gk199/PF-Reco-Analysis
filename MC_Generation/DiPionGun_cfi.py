import FWCore.ParameterSet.Config as cms

# HCAL barrel inner radius used to convert deltaR -> Delta (physical separation).
# Pions shower hadronically, so the relevant surface is the HCAL barrel (~177 cm).
# deltaR = Delta / R_barrel  =>  Delta [cm] = deltaR * 177 cm
# Default below is deltaR=0.4  =>  Delta=70.8 cm

generator = cms.EDProducer("CloseByParticleGunProducer",
    PGunParameters = cms.PSet(
        PartID      = cms.vint32(211, -211),
        NParticles  = cms.int32(2),

        # Fixed pT = 10 GeV per pion.
        # FlatPtGeneration=True means VarMin/VarMax are pT [GeV], not energy.
        FlatPtGeneration = cms.bool(True),
        VarMin           = cms.double(10.0),
        VarMax           = cms.double(10.0),
        MaxVarSpread     = cms.bool(False),
        LogSpacedVar     = cms.bool(False),

        # Barrel geometry: R ~ 177 cm (HCAL barrel inner face), |eta| < 1.4 => Z < ~337 cm.
        # RMin < RMax and ZMin < ZMax are required (strict inequalities).
        ControlledByEta  = cms.bool(False),
        ControlledByREta = cms.bool(False),
        RMin = cms.double(176.9),
        RMax = cms.double(177.1),
        ZMin = cms.double(0.1),
        ZMax = cms.double(337.0),

        MinPhi = cms.double(-3.14159265359),
        MaxPhi = cms.double(3.14159265359),

        # Angular separation between the two pions (cm at the calorimeter surface).
        # deltaR  ~ Delta / R_barrel.  Overridden per point in scan.sh.
        Delta = cms.double(70.8),   # <-> deltaR = 0.4

        Pointing    = cms.bool(True),
        Overlapping = cms.bool(False),
        RandomShoot = cms.bool(False),

        # Timing: particle ip gets vertex time offset = fOffsetFirst + ip*fT (ns).
        # So particle 1 (pi-) arrives fT ns after particle 0 (pi+).
        # UseDeltaT=True: fT sampled from [TMin, TMax].  TMin < TMax required.
        # UseDeltaT=False: fT = 0 (no relative delay).
        # Overridden per scan point in scan.sh.
        UseDeltaT   = cms.bool(False),
        TMin        = cms.double(0.0),
        TMax        = cms.double(0.05),
        OffsetFirst = cms.double(0.0),
    ),
    Verbosity        = cms.untracked.int32(0),
    psethack         = cms.string('di-pion gun'),
    AddAntiParticle  = cms.bool(False),
    firstRun         = cms.untracked.uint32(1),
)

ProductionFilterSequence = cms.Sequence(generator)
