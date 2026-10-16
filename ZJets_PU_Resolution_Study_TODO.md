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
- [x] Build state for the nominal pass — confirm current build matches `PFTestingAlgos/*_original.cc.edit` (same check as the QCD study) before running any `<N>`
- [x] Then for the timing pass: `./SetTimingThreshold.sh 3`, copy .edit files over, `scram b -j8`, before running the same cmsDriver command again
- [x] Get the list of CRAB output files/indices once the generation task finishes (`crab report` or listing the `/store/user/...` output path directly) 
- [x] Batch submission tool: **Condor**. Re-reco is deterministic (no per-job seed handling needed, unlike generation), and reusing `run_phaseScan_job.sh` as a template. 
- [x] `Condor/run_zjetsReco_job.sh` + `Condor/submit_zjetsReco.sh` written
      (modeled on `run_phaseScan_job.sh`/`submit_phaseScan.sh`) and the
      underlying cfg (`MyPFStudy_ZJetsPU_RECO_condor.py`, VarParsing-enabled
      HLT:Fake2,RAW2DIGI,L1Reco,RECO) verified against the real 10-event
      test GEN-SIM-RAW file — exit 0, valid AODSIM output. Example input
      list format: `Condor/input_files_zjetsReco_example.txt`.
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/Condor
./submit_zjetsReco.sh -t input_files_zjets.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_nominal   # test 1 job first
./submit_zjetsReco.sh input_files_zjets.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_nominal      # then all 50

# rebuild with the timing-PF clusterizer (SetTimingThreshold.sh 3, copy files, etc.), then:
./submit_zjetsReco.sh -t input_files_zjets.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_timing3ns
./submit_zjetsReco.sh input_files_zjets.txt /eos/user/g/gkopp/PF_ZJetsPU/reco_timing3ns
```

- [x] Same input list (`input_files_zjets.txt` above) is reused for both the nominal and timing submissions — it's the same 50 GEN-SIM-RAW files either way, only the compiled clusterizer differs
<!-- - [ ] Confirm `142X`/`auto:phase1_2025_realistic` GlobalTag is appropriate for RECO on this sample (should be, matches generation)
- [ ] Same lessons as the phase-scan study apply here too: cap nothing on `maxEvents` since `-n -1` processes each job's full 1000 events (no equivalent of the `_numEvent<N>` rename surprise expected since we're not truncating early this time, but worth re-checking output filenames after the first test job regardless), `xrdcp -f` always, watch EOS quota -->

## 5. NanoAOD

From the AODSIM output, run PAT/NANO config.  

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
No `--geometry` here — PAT/NANO slim already-reconstructed AOD content, no geometry-dependent reconstruction runs at this stage.

- [x] Verified on the 10-event local sample (nominal PF build): PAT's jet/MET sequences run on privately-generated AOD. `GenJet_pt`/`GenMET_pt`/`GenMET_phi` all present as expected.

**`Condor/run_zjetsNano_job.sh` + `Condor/submit_zjetsNano.sh`** with the cfg `MyPFStudy_ZJetsPU_NANO_condor.py`. These jobs are very fast, using `longlunch` which is even longer than needed.
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/Condor
./submit_zjetsNano.sh -t input_files_reco_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/nano_nominal   # test 1 first
./submit_zjetsNano.sh input_files_reco_nominal.txt /eos/user/g/gkopp/PF_ZJetsPU/nano_nominal      # then all 50
# repeat with input_files_reco_timing3ns.txt -> nano_timing3ns once step 4's timing batch is done
```
`input_files_reco_<variant>.txt` here lists step 4's AODSIM output files directly (e.g. `ls /eos/user/g/gkopp/PF_ZJetsPU/reco_nominal/*.root`).
- [x] Confirm `GenJet`/`GenMET` branches are present by default — needed for resolution truth-matching

### EOS transfer failures on the timing3ns batch, and the hardening that followed

Two submissions of the timing3ns batch both lost files to `eosuser`, in *different* places, which made the directory confusing to read:
- **First run:** all 50 jobs reconstructed fine, then 15 failed on the output copy with `xrdcp ... [ERROR] Local error: protocol error: (destination)`, leaving 15 zero-byte stubs. The failures were contiguous in time — every copy landing in a ~35 s window failed, everything outside it succeeded.
- **Resubmission (all 50 again):** 39 jobs failed at the *input*, before cmsRun could read anything — `cmsRun failed with exit code 84`, `XrdCl::File::Open ... [3014] Unable to open file ...; Network is unreachable`, problematic data server `eoshome-g.cern.ch:1094`. Only 10 finished cleanly, and they rewrote files that were already good, so the 12 remaining stubs were never touched.

Both are the same underlying thing: all 50 jobs start within ~30 s and hit one EOS instance at once. Reading the same files from lxplus immediately afterwards worked fine.

**A failed copy can leave a file that looks real.** Validating the whole directory with uproot found 37 complete files (1000 events each), 12 empty, and **one 2.8 MB file that was the right size but unreadable** (`OSError: [Errno 5] Input/output error`). A size check alone does not catch that — check `num_entries`:
```python
import uproot, glob
for f in sorted(glob.glob("/eos/user/g/gkopp/PF_ZJetsPU/nano_timing3ns/*.root")):
    try: print(f, uproot.open(f)["Events"].num_entries)
    except Exception as e: print(f, "BAD:", e)
```
Note `plot_jetmet_resolution.py` will not skip these — uproot raises `OSError: expected Chunk of length 403, received 0 bytes` on a zero-size file rather than dropping it.

**`run_zjetsNano_job.sh` now has** (verified against real EOS, including the give-up path):
- an `xrdcp_retry` helper — 3 attempts, 60 s × attempt backoff — used for *both* transfers
- **input pre-staging to local scratch**, as `run_zjetsReco_job.sh` already did. The script previously read straight over xrootd since the input is CERN-local; that assumption broke, and an open that fails mid-job is unrecoverable whereas a staged copy is just retried
- **removal of the partial destination file** when a copy fails permanently, so nothing downstream ever sees a stub or a truncated file
- **a start stagger** of `(JOB_INDEX * 5 + jitter) mod 120` seconds — spread by index rather than at random, which puts two of the 50 jobs in each 5 s slot instead of leaving RANDOM's gaps and clumps

⚠️ **Do not resubmit a subset with `submit_zjetsNano.sh` as it stands.** The output index is `$(ProcId)`, which restarts at 0 for each submission, so a 13-line input list writes `ZJetsPU_NANO_job0..12` from the *wrong* inputs and overwrites good files. Until the index is derived from the input filename instead, only a full 50-job resubmission is safe. `run_zjetsReco_job.sh` (step 4) stages its input with a single un-retried `xrdcp` and would benefit from the same helper.

## 6. Resolution plotting

**`JetMET_Resolution/plot_jetmet_resolution.py`** — new dedicated subdirectory). Reads NanoAOD with `uproot`/`awkward` (handles the jagged per-jet arrays more naturally than a PyROOT TTree loop — same tool `EventDisplay/PF_HCAL_rechit_cluster.py` already uses), fills/draws with PyROOT to match the `Plotting/*.py` convention. Verified against `ZJetsPU_NANO_job0.root`, 1000 events, 1573 gen-matched jets.
```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/JetMET_Resolution
python3 plot_jetmet_resolution.py \
    --input "/eos/user/g/gkopp/PF_ZJetsPU/nano_timing3ns/*.root" \
    --output resolution_timing3ns.root \
    --pdf resolution_timing3ns.pdf
python3 plot_jetmet_resolution.py \
    --input "/eos/user/g/gkopp/PF_ZJetsPU/nano_nominal/*.root" \
    --output resolution_nominal.root \
    --pdf resolution_nominal.pdf
```
- [x] Match `Jet` to `GenJet` via `Jet_genJetIdx` (NanoAOD's own gen-matching, `>=0` = matched) — done
- [x] Jet resolution: `(Jet_pt - GenJet_pt) / GenJet_pt` for gen-matched jets above a 20 GeV `GenJet_pt` floor (avoids soft/noise-dominated jets) — done, one inclusive histogram (not yet binned in `GenJet_pt`/η)
- [x] MET resolution: both `PFMET_pt - GenMET_pt` and `PuppiMET_pt - GenMET_pt` (absolute difference, not relative — GenMET near zero would blow up a relative metric) — done. Both show a large positive bias vs GenMET (PF: mean 42.9 GeV; Puppi: mean 32.6 GeV) — confirmed this isn't missing JEC/Type-1 MET corrections (`Jet_rawFactor` nonzero for 100% of jets, ~26% avg correction; `CorrT1METJet` populated; GT resolves to the real, official `150X_mcRun3_2025_realistic_v2`, same one used throughout this study) — most likely genuine MET degradation from the aggressive 55-75 in-time PU choice combined with this private sample not exactly matching the JEC payload's validation conditions. Since nominal and timing-PF share the same generator/PU library/GT, this bias should mostly cancel in the differential comparison — not treated as a blocker, but absolute MET numbers from this study aren't officially-calibrated performance figures.
- [x] Nominal vs. timing-PF overlay — done, see step 7 below
- [x] Binning jet resolution by `GenJet_pt` and |η| instead of one inclusive histogram (a timing effect concentrated in one pT/η region would be washed out in the inclusive histogram) — done, `--pt-bins` / `--eta-bins`, default η split at 1.3 / 2.5 / 3.0 to follow HB / HE / HE-HF transition / HF since the timing cut acts on HCAL clusters. pT bins are **20-30, 30-40, 40-60, >60** — after lepton cleaning there are only ~22k hadronic jets per variant, so finer or higher bins are empty (the nominal file has 6222 / 5960 / 7442 / 2339 jets in these four)
- [x] Jet **response** `R = pT_reco/pT_gen` distributions, per (pT, |η|) cell, written out alongside the summary curves so a suspicious summary point can be traced back to the distribution it came from. Binned in 0.05-wide bins over 0 < R < 2.5 — fine enough for the fit to have ~10 bins across the core, coarse enough that a single cell is not mostly empty bins
- [x] Per-bin centre/width from a **truncated Gaussian fit** (iterative, `--fit-nsigma` default ±1.5σ, range re-derived from the previous fit) and from the **median + quantile width** (`σ_eff = (q84-q16)/2`, computed unbinned). Both avoid mean/RMS, which the long low-response tail biases badly. Both are written to the ROOT file and appear in the printed table, but **only one is ever drawn** (`--resolution-estimator`, default `fit`) — plotting both on the same axes was more confusing than useful. Where they disagree strongly the core is not Gaussian: in the 20-30 GeV bin the fit gives σ/median = 0.377 against the quantile 0.344, with fit μ = 1.13 vs. median 1.23, i.e. the fit is fighting a skewed distribution there
- [x] **Relative resolution = width / median response.** Note JECs *are* applied — `Jet_pt` is the corrected pT (`Jet_rawFactor` nonzero for 100% of jets, ~21% mean correction) and `--raw-jets` switches to `Jet_pt*(1-Jet_rawFactor)`. What is *not* redone is the PF calibration after the clusterizer is modified for timing PF, and the JEC payload is the central one derived for standard PF on central MC. Measured on this sample: median response **1.19 corrected vs. 1.02 raw**, i.e. the central JEC over-corrects here (worst in 2.5 < |η| < 3.0, where it applies a factor 1.66 and lands at a response of 1.55). Since the response is not 1, an absolute width is not comparable across bins and dividing by the median is the right normalization — valid across variants only while their responses agree, which step 7 checks per bin
- [x] **Seven pages, on purpose.** 1-3 are the original inclusive distributions, 4 is the inclusive response with the fit drawn on it, 5 is R in pT bins (η-inclusive, area-normalized), 6 is median(R) vs pT and 7 is σ_R/median(R) vs pT, each with one curve per |η| slice. The y axis on 6-7 is spelled out on the page itself (what R is, which σ, over what) rather than left as a bare "σ/median". Earlier drafts also had a fit-vs-quantile comparison page and response/resolution vs |η| pages; those are removed until 1-7 are fully understood, and the per-η response distributions stay in the ROOT file without getting their own pages
- [x] **Lepton cleaning (new default).** AK4 gen jets include the `Z->ee/μμ` decay leptons, so most "jets" in this sample were the leptons: **72% of gen-matched jets above 20 GeV are within ΔR < 0.4 of a dressed gen lepton**, leaving 22067 genuine hadronic jets out of 78462. The give-away is a Jacobian peak at `M_Z/2 ≈ 45 GeV` in the `GenJet_pt` spectrum (a *rise* from 20 to 50 GeV, then a cliff) — not a falling jet spectrum at all. Those lepton-jets are ECAL/tracker objects that any HCAL-clustering change leaves alone, so they diluted the earlier inclusive numbers. Cleaning uses `GenDressedLepton` at **gen level on purpose**: the same jets are then dropped in every variant, whereas a reco-level lepton veto could itself differ between nominal and timing PF and fake a resolution difference. `--no-lepton-cleaning` restores the old behaviour, `--lepton-cleaning-dr` sets the cone
<!-- - [ ] Parallel/perpendicular MET resolution decomposition relative to the
      hadronic recoil — deferred, absolute-difference histograms were
      enough to get a first look and surface the JEC question above
- [ ] Reuse the `dataviz` skill for plot styling/colors if this gets
      revisited beyond the current simple PyROOT canvases -->

## 7. Nominal vs. timing overlay

**`JetMET_Resolution/compare_jetmet_resolution.py`** — overlays the variants produced in step 6. It reads the histogram ROOT files that `plot_jetmet_resolution.py` already writes (not the NanoAOD again), so it is fast and the two scripts stay decoupled: re-run step 6 only when the ntuples change, re-run this whenever the plot styling does. Every page is an overlay on top with a timing/nominal ratio panel below — same layout as `Plotting/compare_qcd_ratio.py`. Six pages by default: the three inclusive distributions (jet pT resolution, PF MET, Puppi MET), the inclusive jet response, then median(R) and σ_R/median(R) vs. GenJet pT, η-inclusive. `--eta-slices 0 1 2 3` (or `all`) adds two more pages per |η| slice once the inclusive ones are understood. The η binning is read off the histogram names/titles rather than re-declared, so a custom `--eta-bins` in step 6 carries through automatically.

```bash
cd /afs/cern.ch/work/g/gkopp/2025_ParticleFlow/CMSSW_15_0_6/src/PF-Reco-Analysis/JetMET_Resolution
python3 compare_jetmet_resolution.py \
    --nominal resolution_nominal.root \
    --timing  resolution_timing3ns.root \
    --labels  3ns \
    --output  resolution_nominal_vs_timing.pdf \
    --rootfile resolution_nominal_vs_timing.root
```

Several timing thresholds can be overlaid against the one nominal baseline at once:
```bash
python3 compare_jetmet_resolution.py \
    --nominal resolution_nominal.root \
    --timing  resolution_timing1ns.root resolution_timing3ns.root resolution_timing5ns.root \
    --labels  1ns 3ns 5ns
```

Other options: `--resolution-estimator {fit,quantile}` picks which width to overlay (both are in the input files), `--normalize` scales the *distribution* pages to unit area — needed if the variants were run over different numbers of events, otherwise the ratio panel measures statistics rather than the PF change, and never applied to the response/resolution curves — plus `--logy`, `--ratio-min`/`--ratio-max` (distribution pages), `--response-ratio-range`/`--resolution-ratio-range` (summary pages), and `--ratio-min-entries` (default 10) which blanks ratio points whose nominal bin is too sparse to say anything, so tail noise doesn't swamp the bulk.

**Response-agreement check.** Because the relative resolution divides by each variant's *own* median response, the comparison is only like-for-like where the responses themselves agree; `--response-tolerance` (default 1%) flags every bin where they don't. On the η-inclusive default pages it passes — the medians agree to better than 1% in all four pT bins. Running `--eta-slices all` it fires in 8 bins, and not randomly:
- **1.3 ≤ |η| < 2.5 (HE):** timing PF raises the median response by **+1.2% to +3.1%**, in every pT bin
- **2.5 ≤ |η| < 3.0 (HE-HF transition):** timing PF *lowers* it by **−2.5% to −3.6%**
- barrel (|η| < 1.3) and HF (|η| > 3) agree within the 1% tolerance everywhere

That is a real HCAL-region-dependent energy shift from the seed-timing cut, and it is what the η-inclusive numbers hide (the two regions pull in opposite directions). In those bins the resolution ratio is partly a response shift, so read the two curves together.

**Results with lepton cleaning (50k events per variant, 22067 vs. 22260 hadronic gen-matched jets).** The inclusive resolution RMS ratio is 0.987 (jet), 1.004 (PF MET), 0.987 (Puppi MET) — ~1%. Per pT bin, η-inclusive, timing/nominal σ_R/median(R): **1.005, 1.021, 1.057, 1.044** for 20-30 / 30-40 / 40-60 / >60 GeV, i.e. timing PF is consistently a few % *worse*, growing with pT. Split by |η| the same ratios scatter both ways (0.92 to 1.16 in HE, 0.70 to 1.16 in the transition region) — those cells are statistics-limited.

So: a consistent but small η-inclusive degradation, plus a genuine region-dependent response shift. Neither is established at the per-(pT, η) level from these files. The honest next step is more statistics — and note 72% of the jets were thrown away by lepton cleaning, so a jet-enriched generator setup would buy far more than simply running more of this sample.

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