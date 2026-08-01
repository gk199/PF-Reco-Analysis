# Z+jets (PU) — Nominal vs. Timing PF Jet/MET Resolution — Production Plan

Goal: generate 50k Z+jets events with Run 3 pileup, re-reco with nominal PF and 3ns seed-timing PF, produce NanoAOD for both, and compare jet/MET resolution (reco vs. gen) between the two.

## Directory separation

- **Generation** (GEN → SIM → DIGI+PU → RAW) happens in the separate CMSSW area: `/afs/cern.ch/work/g/gkopp/2025_ParticleFlow/MC_Generation/CMSSW_15_0_6/src` — follows the existing `DiPionGun` generation done there (with the generator fragment in `Configuration/Generator/python/`, a driver script, per-job `_cfg.py` + log files). New files go in a new `ZJetsPU/` subdirectory there, same pattern as the `DiPionGun_cfi.py`/`scan.sh`.
- **Re-reco, NanoAOD, ntuplization, and plotting** happen here in `PF-Reco-Analysis`, reusing the existing timing-PF infrastructure (`SetTimingThreshold.sh`, `PFTestingAlgos/`, `Condor/`).
- **File processing**: generation goes through CRAB (see step 3), writing to `T2_US_Wisconsin` (no write access to `T2_CH_CERN`), with output at `/hdfs/store/user/gkopp/ZJetsPU/ZJetsPU_nominal_v1/<date>/0000/*.root`. Re-reco (step 4) reads directly from `/store/user/...` via xrootd.

## 0. MC Checks

- No disk accessible DYJetsToLL (Z+jets when only neutral current process is enabled) RAW sample exists — generating is necessary (RelVal-equivalent and central campaigns are tape-only).
- MinBias PU library **is** disk-resident and at CERN: `/MinBias_TuneCP5_13p6TeV-pythia8/Run3Winter25GS-Winter25PU_correctBS_142X_mcRun3_2025_realistic_v4-v1/GEN-SIM` at `T2_CH_CERN`.
- Conditions target: `142X_mcRun3_2025_realistic` / `auto:phase1_2025_realistic`,  matching the PU library and the GlobalTag already used throughout the other RECO configs.

## 1. Generator fragment

- [x] Write a Pythia8 Drell-Yan fragment. Base it on CMSSW's standard DY Pythia8 fragment pattern: `Configuration/Generator/python/ZJetsPU/ZJets_pythia8_cfi.py`.
- [x] Copy into place and build: `cp ZJetsPU/ZJets_pythia8_cfi.py Configuration/Generator/python/ && scram b -j8`
- [x] Sanity-check with a small interactive run (10 events, no PU, then inspect the generated `_cfg.py`) before adding PU
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

**MC Generation: CRAB.** CRAB is used for this (`PrivateMC`: generate from scratch, no input dataset, split by event count) and per-job random seeds are handled automatically (CRAB injects a unique seed per job), and failed jobs retry automatically. CRAB has its own setup (a `crabConfig.py`, `crab submit`/`crab status`).

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
# at the 240min/4hr cap, which then wall-clock-killed 24/50 jobs
# `crab resubmit --maxjobruntime=480` for the failed jobs
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
- [x] Once complete, `crab report` / the task's output path gives the exact list of output `.root` files under `/store/user/...` for step 4

## 4. Re-reco: nominal PF and 3ns seed-timing PF

Generation ran `GEN,SIM,DIGI,L1,DIGI2RAW`, so the output is RAW-tier content — same as the JetMET data re-reco'd (same `RAW2DIGI,L1Reco,RECO` step list as `MyPFStudy_ReReco_MC_condor.py`). This needs full AOD content (jets, MET, tracks, vertices) so we keep more than just HCAL clusters, using `--eventcontent AODSIM`.

**`HLT:Fake2` needed** — generation never ran an actual HLT simulation step (`GEN,SIM,DIGI,L1,DIGI2RAW` skips it, matching the DiPionGun convention), so `TriggerResults::HLT` doesn't exist in these files. The standard NanoAOD sequence (step 5) expects it (HLT trigger-bit tables, MET-filter flags) and crashes without it: `Principal::getByLabel: Found zero products... Looking for process: HLT`. `HLT:Fake2` fixes this — a pass-through HLT menu for private MC that doesn't need a trigger decision, just a valid `TriggerResults` collection to exist. HLT reads RAW-level data (same as real HLT during data-taking), so this runs as an extra step on the *existing* GEN-SIM-RAW files — no need to redo the CRAB generation.

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

- [x] **Verified on the 10-event local sample (nominal PF build)**: `HLT:Fake2,RAW2DIGI,L1Reco,RECO` runs, and the resulting AODSIM fed through step 5's PAT,NANO. `Jet_pt`/`GenJet_pt`/ `Jet_genJetIdx` all present. Note: use `PFMET_pt`/`PFMET_phi` for the resolution study; `PFMET_*` and `PuppiMET_*` are both fully populated under their real names). Leaning `PFMET_pt` over `PuppiMET_pt` for the actual comparison — PUPPI's own per-particle PU reweighting on top of PF candidates could partially mask the effect we're isolating from the seed-timing PF clustering change itself. (? to check)
- [ ] Still need: same local test on the timing-PF (3ns seed timing) build to confirm the pipeline works identically there — same input file, just rebuild the clusterizer and rerun the same three cmsRun commands
- [ ] Build state for the nominal pass — confirm current build matches `PFTestingAlgos/*_original.cc.edit` (same diff check as the QCD study) before running any `<N>`
- [ ] Then for the timing pass: `./SetTimingThreshold.sh 3`, copy in `PFMultiDepthClusterizer_seedTiming.cc.edit` + similar, `scram b -j8`, before running the same cmsDriver command again
- [ ] Get the list of CRAB output files/indices once the generation task finishes (`crab report` or listing the `/store/user/...` output path directly) — `<N>` above is a placeholder until then
- [x] Batch submission tool: **Condor**. Re-reco is deterministic (no per-job seed handling needed, unlike generation), and reusing `run_phaseScan_job.sh` as a template avoids paying the CRAB memory/runtime tuning cost a second time for a different (likely also memory-hungry) workload. 
- [x] `Condor/run_zjetsReco_job.sh` + `Condor/submit_zjetsReco.sh` written
      (modeled on `run_phaseScan_job.sh`/`submit_phaseScan.sh`) and the
      underlying cfg (`MyPFStudy_ZJetsPU_RECO_condor.py`, VarParsing-enabled
      HLT:Fake2,RAW2DIGI,L1Reco,RECO) verified against the real 10-event
      test GEN-SIM-RAW file — exit 0, valid AODSIM output. Example input
      list format: `Condor/input_files_zjetsReco_example.txt`.
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/Condor
./submit_zjetsReco.sh -t input_files_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_nominal   # test 1 job first
./submit_zjetsReco.sh input_files_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_nominal      # then all 50

# rebuild with the timing-PF clusterizer (SetTimingThreshold.sh 3, etc.), then:
./submit_zjetsReco.sh -t input_files_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_timing3ns
./submit_zjetsReco.sh input_files_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_timing3ns
```
- [ ] **Unverified — memory/runtime in `submit_zjetsReco.sh` are guesses**
      (`request_memory = 4000`, `+JobFlavour = "tomorrow"`/24h). RECO
      (tracking/vertexing/PF) on high-PU events may need real tuning like
      generation did (which needed 8000MB/4 cores/480min before it stopped
      failing). Check the `-t` test job's actual memory/runtime before
      trusting these defaults for the full 100-job batch.
- [ ] Same input list (`input_files_nominal.txt` above) is reused for both
      the nominal and timing submissions — it's the same 50 GEN-SIM-RAW
      files either way, only the compiled clusterizer differs
<!-- - [ ] Confirm `142X`/`auto:phase1_2025_realistic` GlobalTag is appropriate for RECO on this sample (should be, matches generation)
- [ ] Same lessons as the phase-scan study apply here too: cap nothing on `maxEvents` since `-n -1` processes each job's full 1000 events (no equivalent of the `_numEvent<N>` rename surprise expected since we're not truncating early this time, but worth re-checking output filenames after the first test job regardless), `xrdcp -f` always, watch EOS quota -->

## 5. NanoAOD

Separate job from step 4, reading the AODSIM output — kept separate rather
than one combined `RAW2DIGI,L1Reco,RECO,PAT,NANO` step, since RECO on
high-PU events already needed real memory tuning during generation (8000MB)
and stacking PAT+NANO on top risks repeating that blind-guessing cycle in a
bigger, more expensive job. If PAT/NANO config needs iteration, this way it
doesn't require redoing the expensive RECO step each time. 

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
      the upstream PF candidates
- [ ] Confirm `GenJet`/`GenMET` branches are present by default (they are,
      standard NanoAOD content for MC) — needed for resolution truth-matching

## 6. Resolution plotting (new script)

`Plotting/compare_zjets_resolution.py` (doesn't exist yet):
- [ ] Match `Jet` to `GenJet` (NanoAOD's `Jet_genJetIdx` branch does this
      directly — no manual dR matching needed)
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
- [x] RECO `--step`/`--eventcontent` for step 4 (needs enough content for
      PAT/NanoAOD jets+MET)
