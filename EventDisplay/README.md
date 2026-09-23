# HCAL Event Display

`HCAL_all_collections.py` is an interactive viewer for HCAL clusters in `pfObjectsNtupler` ntuples. It shows one cluster at a time, drawing its matched rechits in η–φ for each depth. The top row colors the rechits by energy and the bottom row by time.

## Rechit collections and timing

The display can draw two families of rechits. For HBHE rechits, there are also two choices of timing variable.

| Collection | `--rechit-type` | Timing | `--hbhe-time` |
|---|---|---|---|
| **PFRecHits**: the rechits PF clustering actually used; PF energy thresholds already applied | `pfrh` *(default)* | Time in ns (`hbhe_pfrh_time`) | n/a |
| **HBHE rechits**: the full collection, before PF selection | `hbhe` | Raw TDC codes (`hbhe_rechit_tdc`) | `tdc` *(default)* |
| | | Reconstructed time in ns (`hbhe_rechit_time`) | `time` |

Rechits with invalid timing measurement are drawn in **gray** and kept off the colorbar:
- **ns timing (`pfrh`, or `hbhe` with `--hbhe-time time`):** time = −999.
- **TDC codes:** code 3 in HB, 62 or 63 in HE.

## Requirements

Python ≥ 3.9 with `uproot`, `awkward`, `numpy` and `matplotlib`. You also need a display: a local machine, `ssh -Y`, or VNC.

```bash
pip install uproot awkward numpy matplotlib
```

## Quick start

```bash
python3 HCAL_all_collections.py pfObjectsNtuple_standardPF.root
```

The filename is the only required argument. It is looked up in the current directory unless you give a full path or set `--input-dir`. At startup the script prints the branches and cuts it is using.

Use **Prev/Next Cluster** and **Prev/Next Event** to navigate; both wrap around at the ends.

## Command-line options

| Option | Default | Applies to | Description |
|---|---|---|---|
| `filename` | *(required)* | all | ROOT file from the ntupler |
| `--input-dir` | current directory | all | Folder containing the file (ignored if `filename` is an absolute path) |
| `--tree-dir` | `pfObjectsNtupler` | all | Directory *inside* the ROOT file that holds `pfTree` |
| `--rechit-type` | `pfrh` | all | `pfrh` or `hbhe` |
| `--hbhe-time` | `tdc` | `hbhe` | `tdc` (raw codes) or `time` (ns) |
| `--pf-energy-cuts` / `--no-pf-energy-cuts` | on | `hbhe` | Apply depth-dependent PF energy thresholds |
| `--cluster-emin` | `NONE` | all | Only show clusters with E > value [GeV] |
| `--start-event` | `5` | all | Event index to open first (wraps if too large) |

**`--hbhe-time`**
- **`tdc`:** drops rechits with negative codes. Codes with no measurement are drawn gray. The colorbar is 0–2 (HB scale) for pure HB clusters, and 0–61 (HE scale) for any cluster with HE rechits, including mixed HB+HE clusters.
- **`time`:** no TDC logic. The colorbar spans the min–max of the valid times in the cluster.

**`--pf-energy-cuts`** drops HBHE rechits below the PF threshold for their subdetector and depth:

| Depth | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---|---|---|---|---|---|---|
| HB [GeV] | 0.6 | 0.4 | 0.5 | 0.5 | – | – | – |
| HE [GeV] | 0.2 | 0.3 | 0.3 | 0.3 | 0.3 | 0.3 | 0.3 |

**Options that only apply to `hbhe`:**
- **`--hbhe-time`:** with `pfrh`, the script prints a note and ignores it.
- **`--pf-energy-cuts` / `--no-pf-energy-cuts`:** with `pfrh`, the script prints `currently using pfrh collections, PF cuts on rechits already applied` and ignores it.

**`--cluster-emin`:** clusters always need E > 0. With a value set, clusters at or below it are hidden from the list and skipped by the cluster buttons.

**`--input-dir` vs `--tree-dir`:**
- **`--input-dir`** is a folder on disk.
- **`--tree-dir`** is the ntupler module name inside the ROOT file. Change it only if the ntupler ran under a different label.

## What the defaults give you

Running with just a filename is equivalent to:

```bash
python3 HCAL_all_collections.py <file> --input-dir . --tree-dir pfObjectsNtupler \
    --rechit-type pfrh --cluster-emin NONE --start-event 5
```

That is:
- **Rechits:** PFRecHits, with time in ns and −999 shown in gray.
- **Energy cuts:** none added by the script, since PF already applied its thresholds.
- **Clusters:** every cluster with E > 0.
- **Start:** event 5.

## Examples

```bash
# HBHE rechits, raw TDC codes, PF thresholds applied
python3 HCAL_all_collections.py ntuple.root --rechit-type hbhe

# HBHE rechits in ns, no energy thresholds
python3 HCAL_all_collections.py ntuple.root --rechit-type hbhe --hbhe-time time --no-pf-energy-cuts

# Only clusters above 30 GeV, from the first event
python3 HCAL_all_collections.py ntuple.root --cluster-emin 30 --start-event 0
```

## Things to keep in mind

- **Mixed clusters in TDC mode.** Mixed HB+HE clusters use the HE scale (0–61), so their HB rechits (codes 0–2) all look dark. `--hbhe-time time` avoids this.
- **Color scales are per cluster.** Don't compare colors between clusters.
- **Zoom and window size.** Each cluster is shown in a square window of at least ±0.27, widened when needed so every selected rechit (up to ±0.4) is fully visible. The panels resize with the window; maximize it on small screens.
