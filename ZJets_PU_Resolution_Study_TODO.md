# Z+jets (PU) — Nominal vs. Timing PF Jet/MET Resolution — Production Plan

Goal: generate 50k DYJetsToLL (Z+jets) events with realistic Run3 pileup, re-reco with nominal PF and 3ns seed-timing PF, produce NanoAOD for both, and compare jet/MET resolution (reco vs. gen) between the two.

## Directory separation

- **Generation** (GEN → SIM → DIGI+PU → RAW) happens in the separate CMSSW area: `/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src` — mirrors the existing `DiPionGun` convention there (generator fragment in `Configuration/Generator/python/`, a driver script, per-job `_cfg.py` + log files). New files go in a new `ZJetsPU/` subdirectory there, same pattern as the top-level `DiPionGun_cfi.py`/`scan.sh`.
- **Re-reco, NanoAOD, ntuplization, and plotting** happen here in `PF-Reco-Analysis`, reusing the existing timing-PF infrastructure (`SetTimingThreshold.sh`, `PFTestingAlgos/`, `Condor/`).
- **Handoff point**: generation goes through CRAB (see step 3), writing to `T2_US_Wisconsin` (no write access to `T2_CH_CERN`), landing at `/store/user/<username>/PF_ZJetsPU/gensim/...` — **not** the personal `/eos/user/g/gkopp/` area used everywhere else in this study, and not CERN storage at all. Re-reco (step 4) reads directly from `/store/user/...` via xrootd, same as any other CMS dataset; no manual copy needed, but don't confuse the two EOS namespaces when looking for files.

Given the EOS quota problems from the phase-scan study, **plan storage deliberately this time**: delete GEN-SIM-RAW files after they're re-reco'd (same "convert to smaller product, then delete the large intermediate" pattern used for `nominal_fresh/` → its ntuple), rather than letting all three tiers (RAW, RECO, NanoAOD) accumulate simultaneously.

## 0. Confirmed inputs (already checked in this conversation)

- No disk-resident DYJetsToLL RAW sample exists — generating from scratch is necessary (RelVal-equivalent and central campaigns are tape-only).
- MinBias PU library **is** disk-resident and CERN-local: `/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter25GS-Winter25PU_correctBS_142X_mcRun3_2025_realistic_v4-v1/GEN-SIM` at `T2_CH_CERN` (plus FNAL, NCHC disk) — no transfer needed.
- Conditions target: `142X_mcRun3_2025_realistic` / `auto:phase1_2025_realistic`,  matching both the PU library and the GlobalTag already used throughout this repo's other RECO configs.

## 1. Generator fragment

- [x] Write a plain Pythia8 Drell-Yan fragment. Base it on CMSSW's standard DY Pythia8 fragment pattern: `Configuration/Generator/python/ZJetsPU/ZJets_pythia8_cfi.py`.
- [x] Copy into place and build: `cp ZJetsPU/ZJets_pythia8_cfi.py Configuration/Generator/python/ && scram b -j8`
- [x] Sanity-check with a tiny interactive run (10 events, no PU, then inspect the generated `_cfg.py`) before touching PU at all
``` bash
cmsDriver.py Configuration/Generator/python/ZJets_pythia8_cfi.py \
    --mc \
    --step GEN,SIM,DIGI,L1,DIGI2RAW \
    --era Run3_2025 \
    --geometry DB:Extended \
    --conditions auto:phase1_2025_realistic \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --pileup NoPileUp \
    --eventcontent RAWSIM \
    --datatier GEN-SIM-RAW \
    --fileout "file:ZJetsPU_job_10eventTest_GEN-SIM-RAW.root" \
    --python_filename "ZJetsPU_job_10eventTest_cfg.py" \
    -n 10 \
    --no_exec

cmsRun ZJetsPU_job_10eventTest_cfg.py
```

## 2. Confirm the PU scenario name

- [x] `--pileup Run3_Flat55To75_PoissonOOTPU` as used for pion gun MC
- [x] `--pileup_input das:/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter25GS-Winter25PU_correctBS_142X_mcRun3_2025_realistic_v4-v1/GEN-SIM`

## 3. Generation config (combined GEN,SIM,DIGI,L1,DIGI2RAW — one step, matching the DiPionGun pattern) + CRAB submission

**MC Generation: CRAB.** CRAB is built for this (`PrivateMC`: generate from scratch, no input dataset, split by event count) and per-job random seeds are handled automatically (CRAB injects a unique seed per job), and failed jobs retry automatically (unlike Condor scripts). CRAB has its own setup (a `crabConfig.py`, `crab submit`/`crab status`).

**Generate the python config** (single cfg — CRAB reuses this across all jobs and injects a per-job seed itself):
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src
cmsenv

cmsDriver.py Configuration/Generator/python/ZJets_pythia8_cfi.py \
    --mc \
    --step GEN,SIM,DIGI,L1,DIGI2RAW \
    --era Run3_2025 \
    --geometry DB:Extended \
    --conditions auto:phase1_2025_realistic \
    --beamspot Realistic25ns13p6TeVEarly2022Collision \
    --pileup Run3_Flat55To75_PoissonOOTPU \
    --pileup_input "das:/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter25GS-Winter25PU_correctBS_142X_mcRun3_2025_realistic_v4-v1/GEN-SIM" \
    --eventcontent RAWSIM \
    --datatier GEN-SIM-RAW \
    --fileout "file:ZJetsPU_GEN-SIM-RAW.root" \
    --python_filename "ZJetsPU_cfg.py" \
    --nThreads 4 \
    -n 1000 \
    --no_exec
```
- [x] Confirmed the 10-event local test (with PU) runs successfully — reuse that validated fragment/config, just regenerate at `-n 1000` (the per-job event count CRAB will use, `unitsPerJob` below) as the actual pset CRAB submits
- [x] `-n` in this pset should equal `config.Data.unitsPerJob` in the CRAB config below — CRAB overrides `maxEvents` at submission using that value, but keep them consistent for clarity when testing locally
- [x] `--nThreads` here must equal `config.JobType.numCores` in the CRAB config below

**`crabConfig.py`** (new file, `ZJetsPU/crabConfig.py` in the `MC_Generation` area).
Follows the same pattern as the earlier `crab_DisplacedHcalJetNTuplizer_MC_cfg.py` wrapper (writing to Wisconsin), and CRAB builds the `/store/user/...` output path itself from `requestName`/`outputPrimaryDataset`/`outputDatasetTag`, reported after submission via `crab status`/`crab report`:
```python
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'ZJetsPU_gensim_v1'
config.General.workArea = 'crab_ZJetsPU'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'ZJetsPU_cfg.py'
# same filename in every job's pset is fine — CRAB appends the job number
# at stage-out (ZJetsPU_GEN-SIM-RAW.root -> ZJetsPU_GEN-SIM-RAW_1.root, _2.root, ...)
config.JobType.outputFiles = ['ZJetsPU_GEN-SIM-RAW.root']
# Jump to 4 cores/8000MB 
config.JobType.numCores = 4
config.JobType.maxMemoryMB = 8000
# runtimes turned out to be 2:26:17 average, up to 4:00:36 —
# right at the 240min/4hr cap, which then wall-clock-killed 24/50 jobs
# (exit 50664). Fixed via `crab resubmit --maxjobruntime=480` for the
# already-failed jobs; bumping the base config here too for future fresh
# submissions.
config.JobType.maxJobRuntimeMin = 480

config.Data.outputPrimaryDataset = 'ZJetsPU'
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 1000
config.Data.totalUnits = 50000        # 50 jobs x 1000 events
config.Data.publication = False
config.Data.outputDatasetTag = 'ZJetsPU_nominal_v1'

config.Site.storageSite = 'T2_US_Wisconsin'
```
- [x] `source /cvmfs/cms.cern.ch/crab3/crab.sh` (or however CRAB is set up in this environment) before `crab submit -c crabConfig.py`
- [x] `crab status -d crab_ZJetsPU/crab_ZJetsPU_gensim_v1` to monitor
- [ ] Once complete, `crab report` / the task's output path gives the exact list of output `.root` files under `/store/user/...` for step 4

## 4. Re-reco: nominal PF and 3ns seed-timing PF

Generation already ran `GEN,SIM,DIGI,L1,DIGI2RAW`, so the output is RAW-tier content — same situation as the real JetMET data re-reco'd elsewhere in this repo. Same `RAW2DIGI,L1Reco,RECO` step list as
`MyPFStudy_ReReco_MC_condor.py` (proven throughout this whole study), but **not** that config directly — this needs standard full AOD content (jets, MET, tracks, vertices) for PAT downstream, not that config's HCAL-cluster-only `outputCommands`. `--eventcontent AODSIM` gives the standard full set.

**`HLT:Fake2` needed** — generation never ran an actual HLT simulation step (`GEN,SIM,DIGI,L1,DIGI2RAW` skips it, matching the DiPionGun convention), so `TriggerResults::HLT` genuinely doesn't exist in these files. The standard NanoAOD sequence (step 5) unconditionally expects it (HLT trigger-bit tables, MET-filter flags) and crashes without it: `Principal::getByLabel: Found zero products... Looking for process: HLT`. `HLT:Fake2` is the fix for this — a pass-through HLT menu for private MC that doesn't need a real trigger decision, just a valid `TriggerResults` collection to exist. HLT reads RAW-level data (same as real HLT during data-taking), so this runs as an extra step on the *existing* GEN-SIM-RAW files — no need to redo the CRAB generation.

```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis

cmsDriver.py \
    --mc \
    --step HLT:Fake2,RAW2DIGI,L1Reco,RECO \
    --era Run3_2025 \
    --geometry DB:Extended \
    --conditions auto:phase1_2025_realistic \
    --eventcontent AODSIM \
    --datatier AODSIM \
    --filein "file:ZJetsPU_GEN-SIM-RAW.root" \
    --fileout "file:ZJetsPU_<variant>_AOD.root" \
    --python_filename "ZJetsPU_<variant>_RECO_cfg.py" \
    -n -1 \
    --no_exec

cmsRun ZJetsPU_<variant>_RECO_cfg.py
```
`<variant>` = `nominal` or `timing3ns`; `<N>` = the per-job index from the CRAB output filenames (`ZJetsPU_GEN-SIM-RAW_1.root`, `_2.root`, ...). Same config for both variants — only the compiled clusterizer plugin differs between the two build states below, not the cmsDriver command itself.

- [x] **Verified on the 10-event local sample (nominal PF build)**: `HLT:Fake2,RAW2DIGI,L1Reco,RECO` runs cleanly, and the resulting AODSIM fed through step 5's PAT,NANO no longer hits the `TriggerResults::HLT` crash. `Jet_pt`/`GenJet_pt`/ `Jet_genJetIdx` all present. Note: use `PFMET_pt`/`PFMET_phi` for the resolution study; `PFMET_*` and `PuppiMET_*` are both fully populated under their real names). Leaning `PFMET_pt` over `PuppiMET_pt` for the actual comparison — PUPPI's own per-particle PU reweighting on top of PF candidates could partially mask the effect we're isolating from the seed-timing PF clustering change itself.
- [ ] Still need: same local test on the timing-PF (3ns seed timing) build to confirm the pipeline works identically there — same input file, just rebuild the clusterizer and rerun the same three cmsRun commands
- [ ] Build state for the nominal pass — confirm current build matches `PFTestingAlgos/*_original.cc.edit` (same diff check as the QCD study) before running any `<N>`
- [ ] Then for the timing pass: `./SetTimingThreshold.sh 3`, copy in `PFMultiDepthClusterizer_seedTiming.cc.edit` + similar, `scram b -j8`, before running the same cmsDriver command again
- [ ] Get the list of CRAB output files/indices once the generation task finishes (`crab report` or listing the `/store/user/...` output path directly) — `<N>` above is a placeholder until then
- [ ] Confirm `142X`/`auto:phase1_2025_realistic` GlobalTag is appropriate for RECO on this sample (should be, matches generation)
- [x] Batch submission tool: **Condor**. Re-reco is deterministic (no per-job seed handling needed, unlike generation), and reusing `run_phaseScan_job.sh` as a template avoids paying the CRAB memory/runtime tuning cost a second time for a different (likely also memory-hungry) workload. 
- [ ] Write `Condor/run_zjetsReco_job.sh` (new, modeled on `run_phaseScan_job.sh`): stage input GEN-SIM-RAW locally via `xrdcp`, `cmsRun` the RAW2DIGI,L1Reco,RECO cfg, `xrdcp -f` the AODSIM output to EOS. 50 files x 2 PF variants = 100 jobs total (submit nominal and timing as two separate batches, matching the two build states)
- [ ] Same lessons as the phase-scan study apply here too: cap nothing on `maxEvents` since `-n -1` processes each job's full 1000 events (no equivalent of the `_numEvent<N>` rename surprise expected since we're not truncating early this time, but worth re-checking output filenames after the first test job regardless), `xrdcp -f` always, watch EOS quota

## 5. NanoAOD

Separate job from step 4, reading the AODSIM output — kept separate rather
than one combined `RAW2DIGI,L1Reco,RECO,PAT,NANO` step, since RECO on
high-PU events already needed real memory tuning during generation (8000MB)
and stacking PAT+NANO on top risks repeating that blind-guessing cycle in a
bigger, more expensive job. If PAT/NANO config needs iteration, this way it
doesn't require redoing the expensive RECO step each time. Delete the AODSIM
intermediate after NanoAOD succeeds (same "convert then delete" pattern as
`nominal_fresh/` in the phase-scan study).

```bash
cmsDriver.py \
    --mc \
    --step PAT,NANO \
    --era Run3_2025 \
    --conditions auto:phase1_2025_realistic \
    --eventcontent NANOAODSIM \
    --datatier NANOAODSIM \
    --filein "file:ZJetsPU_<variant>_AOD.root" \
    --fileout "file:ZJetsPU_<variant>_NANO.root" \
    --python_filename "ZJetsPU_<variant>_NANO_cfg.py" \
    -n -1 \
    --no_exec

cmsRun ZJetsPU_<variant>_NANO_cfg.py
```
No `--geometry` here — PAT/NANO slim already-reconstructed AOD content,
no geometry-dependent reconstruction runs at this stage.

- [ ] Verify PAT's jet/MET sequences run cleanly on the timing-PF variant —
      should work since PAT only depends on standard collection names
      (`particleFlow`, `ak4PFJets`, etc.), not on which clusterizer produced
      the upstream PF candidates, but this is untested in this repo and
      worth confirming on one file before scaling up
- [ ] Confirm `GenJet`/`GenMET` branches are present by default (they are,
      standard NanoAOD content for MC) — needed for resolution truth-matching

## 6. Resolution plotting (new script)

`Plotting/compare_zjets_resolution.py` (doesn't exist yet):
- [ ] Match `Jet` to `GenJet` (NanoAOD's `Jet_genJetIdx` branch does this
      directly — no manual ΔR matching needed)
- [ ] Jet resolution: `(Jet_pt[genJetIdx] - GenJet_pt) / GenJet_pt`, binned
      in `GenJet_pt` (and optionally η), nominal vs. timing PF overlaid —
      same ratio/overlay pattern as `compare_qcd_ratio.py`
- [ ] MET resolution: `(PFMET_pt - GenMET_pt)` (`PFMET_*`, not `MET_pt` —
      see step 4 note) or the standard parallel/perpendicular resolution
      decomposition relative to the hadronic recoil, nominal vs. timing PF
      overlaid
- [ ] Reuse the `dataviz` skill for the actual plot styling/colors when
      building this

## Open decisions before starting
- [x] `--pileup` scenario (`Run3_Flat55To75_PoissonOOTPU`) and beamspot
      (`Realistic25ns13p6TeVEarly2022Collision`) — set in step 3, PU
      confirmed working via the local 10-event test
- [x] Events per generation job — 1000/job × 50 jobs via CRAB `unitsPerJob`/`totalUnits`
- [x] Storage site (`T2_US_Wisconsin`) — matches the existing
      `crab_DisplacedHcalJetNTuplizer_MC_cfg.py` convention; output path
      left to CRAB's default (`crab status`/`crab report` gives the exact
      `/store/user/...` path after submission)
- [ ] RECO `--step`/`--eventcontent` for step 4 (needs enough content for
      PAT/NanoAOD jets+MET, not the narrow phase-scan ntupler subset)
- [ ] EOS storage budget for this study before starting, given the recent
      quota exhaustion — RECO+NanoAOD for 50k PU events × 2 PF variants
      could be substantial; delete each tier after converting to the next
      (RAW itself lives in CRAB's `/store/user/` area, separate from the
      personal-quota pressure hit during the phase-scan study)
