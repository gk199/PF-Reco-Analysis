# Di-Pion Gun

CMSSW generator for pairs of pions (π⁺ π⁻) with configurable opening angle and timing offset. Intended for particle flow studies in the CMS barrel.

## Files

| File | Description |
|---|---|
| `DiPionGun_cfi.py` | Generator fragment using `CloseByParticleGunProducer` |
| `scan.sh` | Runs GEN-SIM-DIGI over a grid of opening angles and timing offsets |

## Setup

```bash
cd CMSSW_15_0_6/src
cmsenv

# Copy the generator fragment
cp DiPionGun_cfi.py Configuration/Generator/python/

scram b -j8
```

## Running a scan

Edit the scan points at the top of `scan.sh`:

```bash
DR_VALUES=(0.1 0.4 0.8)   # opening angle between the two pions in delta-R
DT_VALUES=(0.0 1.0 5.0)   # timing offset of pi- relative to pi+ in ns
```

Then run:

```bash
chmod +x scan.sh
./scan.sh
```

Each scan point produces:
- `DiPionGun_DR<X>_DT<Y>_GEN-SIM-RAW.root` — output
- `DiPionGun_DR<X>_DT<Y>_cfg.py` — full CMSSW config (useful for debugging)
- `log_DR<X>_DT<Y>.txt` — stdout/stderr from cmsRun

## Parameters

**Opening angle** is specified as ΔR and converted to a physical separation at the HCAL barrel inner face (R = 177 cm):

```
Delta [cm] = deltaR × 177 cm
```

For example, ΔR = 0.4 corresponds to a 70.8 cm separation at the calorimeter surface. Pions shower hadronically so the HCAL barrel sets the relevant angular scale.

**Timing** is the arrival time of the π⁻ relative to the π⁺ at the production vertex. At DT = 0 the timing offset is disabled entirely. For DT > 0, `UseDeltaT=True` is set and `TMin`/`TMax` are bracketed tightly around the requested value (±0.001 ns), giving effectively fixed timing.

## Conditions

Configured for Run 3 2024 realistic conditions:

```
--conditions auto:phase1_2024_realistic
--beamspot Realistic25ns13p6TeVEarly2022Collision
--era Run3_2024
--geometry DB:Extended
```

Update these in `scan.sh` to target a different era.