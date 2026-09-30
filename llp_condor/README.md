# LLP re-reconstruction and automatic merging

The default output is **one analysis ntuple ROOT file per PF algorithm**,
combining all selected input files:

```
/eos/user/c/chtong/Public/Rereco/LLP/<run-id>/pfObjectsNtuple_standardPF.root
/eos/user/c/chtong/Public/Rereco/LLP/<run-id>/pfObjectsNtuple_seedTimingPF.root
```

Each contains the existing PFObjectsNtupler tree structure, with the entries
from the selected input files concatenated. Standard and seed timing remain
separate files. The threshold, selected inputs and run settings are recorded
in the logs and completion metadata.

**By default the reconstructed EDM files are temporary:** they are produced and
checked on the workers, but only the ntuples are saved as final ROOT outputs.
Use `SAVE_RECO=1` to also retain one merged EDM file per algorithm:
`pf_reReco_standardPF.root` and `pf_reReco_seedTimingPF.root` (four final ROOT
files total). EDM files and analysis ntuples are different formats and are
merged separately with the appropriate tools.

## Run on lxplus

Download the updated `llp_condor_test.tar.gz` into `PF-Reco-Analysis`, then:

```bash
cd /afs/cern.ch/user/c/chtong/PF/CMSSW_15_0_6/src
cmsenv
cd PF-Reco-Analysis

tar -xzf llp_condor_test.tar.gz
voms-proxy-init --voms cms --valid 96:00
bash llp_condor/submit.sh
```

This is a new submission; replacing the downloaded scripts does not alter jobs
already submitted. The original CMSSW area is only read. Each input job receives
its own source snapshot, compiles and runs both PF variants sequentially, and
uploads intermediate ntuples. A Condor DAG automatically starts the merge job
only after **every selected input job succeeds**. The default two-input run
therefore has two reconstruction jobs plus one merge job (and a DAG manager).

The source snapshot includes the current local `PFObjectsNtupler` C++ plugin,
including any edits you have installed before submission. This bundle does not
supply or replace that plugin. It also needs your `SetTimingThreshold.sh`, the
seven PF variant source/config files, and the local PFClusterProducer package.

## Choose inputs

Run these commands from `PF-Reco-Analysis`, after `cmsenv`. Settings can also be
changed directly in the defaults near the top of `llp_condor/submit.sh`.

### Number and position of files

```bash
# First two files (the default: step2_177.root and step2_397.root).
bash llp_condor/submit.sh

# First ten files.
NFILES=10 bash llp_condor/submit.sh

# Skip the first ten files and process the next two (entries 11 and 12).
FILE_OFFSET=10 NFILES=2 bash llp_condor/submit.sh

# All files from the input list (only use after the small test works).
NFILES=0 bash llp_condor/submit.sh
```

`FILE_OFFSET` is zero-based. Blank lines and lines beginning with `#` are ignored;
duplicate URLs are removed while preserving first-occurrence order, then the
offset and count are applied. `NFILES=0` means all entries remaining after the
offset. If fewer than NFILES entries remain, all remaining entries are used.
The exact selected URLs are saved in `selected_inputs.txt` in the submission
folder. The supplied list has 9,980 entries; for production at that scale,
replace per-job compilation with two prebuilt runtime archives first.

### Different input list

Create a text file with one `/store/...root` path or `root://...root` URL per line:

```bash
INPUT_LIST=/absolute/path/my_inputs.txt NFILES=2 bash llp_condor/submit.sh

# All entries in the alternative list.
INPUT_LIST=/absolute/path/my_inputs.txt NFILES=0 bash llp_condor/submit.sh
```

Bare `/store/...` paths use the Fermilab endpoint by default. Full XRootD URLs
are kept as written. For a different site, use full URLs or set `SOURCE_HOST`.
Input lists must contain compatible RECO/RECOSIM files with the products needed
by this reconstruction recipe; changing the list does not automatically adapt
the sequence to MiniAOD, NanoAOD, a different era, or other file content.

### Explicit filenames without editing a list

```bash
INPUT_FILES='/store/path/first.root /store/path/second.root' NFILES=0 \
    bash llp_condor/submit.sh
```

Replace the example paths with real files. `INPUT_FILES` takes precedence over
`INPUT_LIST`. It is a whitespace-separated list of `/store/...` paths or full
`root://...` URLs. Offset and count still apply; `NFILES=0` selects all explicit
filenames supplied.

### Events and timing

```bash
# First 100 events of EACH selected input, for EACH PF algorithm.
MAX_EVENTS=100 bash llp_condor/submit.sh

# Combine settings.
INPUT_LIST=/absolute/path/my_inputs.txt NFILES=5 MAX_EVENTS=500 \
    TIMING_THRESHOLD=2.5 bash llp_condor/submit.sh
```

`MAX_EVENTS=-1` processes all events. With two files and `MAX_EVENTS=100`, each
merged algorithm file contains up to 200 events, subject to the inputs' sizes
and the analyzer's filling behavior. It is a per-input limit, not a limit on the
merged file. The same input selection and event limit are used for both PF
algorithms.

## Output options

```bash
# Default: two merged ntuples, remove temporary per-input ROOT files after success.
MERGE_OUTPUTS=1 SAVE_RECO=0 KEEP_PARTS=0 bash llp_condor/submit.sh

# Also save two merged reconstructed EDM files (four final ROOT files).
SAVE_RECO=1 bash llp_condor/submit.sh

# Retain intermediate per-input files as well as the merged results.
KEEP_PARTS=1 bash llp_condor/submit.sh

# Disable merging; save separate per-input ntuples under parts/<variant>/.
MERGE_OUTPUTS=0 bash llp_condor/submit.sh
```

With `MERGE_OUTPUTS=0 SAVE_RECO=1`, the separate per-input EDM files are saved too.
`KEEP_PARTS` has no effect when merging is disabled. Intermediate files are under
`<run-id>/parts/standardPF/` and `<run-id>/parts/seedTimingPF_4p0ns/`. They are
removed only after both final algorithm outputs pass verification and are
uploaded. Completion JSON files remain for audit. A cleanup failure leaves
validated final outputs in place and reports the remaining temporary paths.

The merge reads exactly this submission's input manifest, not a wildcard over
old ROOT files. All completion records must match the selected input URL,
threshold, event limit and output mode. Ntuple entry counts must equal the sum
of input tree entries; optional EDM counts must also match. A missing or failed
input prevents a final completed merge. `merge_complete.json` at the run root
is written only after both algorithms' final files have been validated.

## Other settings

| Variable | Default | Meaning |
|---|---|---|
| `INPUT_LIST` | Bundle's `input_files.txt` | Source file list |
| `NFILES` | `2` | Number of input files; `0` means all remaining |
| `FILE_OFFSET` | `0` | Number of unique input entries to skip |
| `MAX_EVENTS` | `-1` | Events per input file per algorithm |
| `TIMING_THRESHOLD` | `4.0` | Seed timing threshold in ns |
| `MERGE_OUTPUTS` | `1` | Automatically merge all selected inputs |
| `SAVE_RECO` | `0` | Save EDM outputs in addition to ntuples |
| `KEEP_PARTS` | `0` | Retain per-input ROOT files after merging |
| `NCPUS` | `4` | Reconstruction threads and build parallelism |
| `MEMORY_MB` | `12000` | Reconstruction job memory |
| `DISK_MB` | `40000` | Scratch disk requested by each worker/merger |
| `JOB_FLAVOUR` | `tomorrow` | CERN 24-hour limit per job |
| `ERA` | `Run3_2023` | Reconstruction era |
| `CONDITIONS` | `auto:phase1_2023_realistic_postBPix` | Conditions |
| `EOS_BASE` | `/eos/user/c/chtong/Public/Rereco/LLP` | Output parent directory |
| `BATCH_BASE` | Sibling `llp_condor_runs` beside the CMSSW release directory | Submit files and logs |
| `PREPARE_ONLY` | `0` | Use `1` to prepare without submitting |

The merger requests one CPU and 4 GB memory. Increase DISK_MB if a combined
output can exceed the default worker space, especially when SAVE_RECO=1.
A proxy with at least 30 hours remaining is checked at submission; request the
96-hour proxy shown above for the reconstruction-plus-merge workflow and allow
for queue delays. An expired proxy or restricted Fermilab group permissions can
still prevent remote I/O.

## Reconstruction details

The reconstruction recipe is unchanged from the previous bundle:
`RECO:reconstruction_fromRECO`, running in CMSSW_15_0_6 with 2023 BPix conditions.
The inspected files lack `rawDataCollector` but include stored detector rechits
and other RECO inputs. The analyzer runs in the same cmsRun as reconstruction,
reading PF candidates, clusters and blocks explicitly from `ReRECO` and detector
rechits from input `RECO`. The local EDM output uses `keep *` to retain references.
Its new PF collections are checked even when the EDM file is only temporary.

Ntuples are merged with ROOT `hadd`. Optional EDM outputs are merged through
CMSSW PoolSource/PoolOutputModule, never through `hadd`. Event-ID duplicate
checking is enabled for EDM merging; an event-count discrepancy stops the merge
rather than silently accepting a loss of events.

Source packaging excludes ROOT files, plots, archives, logs, Git metadata and
caches. Any custom calibration files in excluded formats must be shipped
separately. Standard release data comes from CVMFS. BATCH_BASE must be outside
CMSSW_BASE/src and on AFS/work, not EOS; use your AFS work area if home quota is
insufficient.

## Monitoring and validation

```bash
condor_q
condor_q -hold
```

The script prints the submission directory. Inspect its `logs/` and `reports/`
folders for reconstruction or merge failures. In merged mode, manual submission
uses `condor_submit_dag workflow.dag` from that directory; submitting
`submit.sub` alone runs only reconstruction and does not schedule the merger.
Failed jobs are held for inspection. There are no automatic endless retries.

Shell/Python syntax, input selection, DAG generation and merge orchestration
were checked outside CERN. **No real CMSSW, HTCondor or ROOT merge has been run
from this chat.** The first two-file test still needs to exercise your custom
plugin and the CMSSW_13_0_13 input compatibility with CMSSW_15_0_6. Send the first
exception and generated config if cmsRun fails.

References:
- https://github.com/cms-sw/cmssw/blob/CMSSW_15_0_6/Configuration/StandardSequences/python/Reconstruction_cff.py
- https://root.cern.ch/doc/v634/hadd_8cxx.html
- https://htcondor.readthedocs.io/en/main/automated-workflows/dagman-introduction.html
- https://batchdocs.web.cern.ch/local/submit.html
