# Timing PF vs standard PF on LLP decays

Compares HCAL PF clusters from `standardPF` and `seedTimingPF` on the same LLP events, using gen truth matching.

| File | Role |
|---|---|
| `TruthInfoHelper.py` | Truth matching and expected delay (port of [TruthInfoHelper.cxx](https://github.com/gk199/Run3-HCAL-LLP-Analysis/blob/main/DisplacedHcalJetAnalyzer/src/TruthInfoHelper.cxx), jets replaced by clusters) |
| `label_llp_clusters.py` | Writes label trees for one ntuple |
| `compare_timing_vs_standard.py` | Reads both label files, writes the PDF |

## Usage

```bash
python3 label_llp_clusters.py --input pfObjectsNtuple_standardPF.root   --output labels_standardPF.root
python3 label_llp_clusters.py --input pfObjectsNtuple_seedTimingPF.root --output labels_seedTimingPF.root

python3 compare_timing_vs_standard.py --standard labels_standardPF.root --timing labels_seedTimingPF.root \
    --output timing_vs_standard_LLP.pdf
```

| Script | Option | Default | Meaning |
|---|---|---|---|
| label | `--deltaR` | 0.4 | Matching cone |
| label | `--clusterE` | 0 | Cluster E cut for the efficiency flag only |
| compare | `--clusterE` | 2 GeV | Cluster E cut for pages 7 to 15 and 18 to 21 |
| compare | `--energyFrac` | 0.4 0.65 0.8 | Correctly matched: ΣE/E_b above this (3 values) |
| compare | `--timeCut` | 1 0.5 0.1 ns | Correctly matched: distance to the time line below this (3 values) |
| compare | `--dRCut` | 0.3 0.15 0.05 | Correctly matched: mean ΔR below this (3 values) |
| compare | `--timeOffset` | median | Offset of the time line (see below) |
| compare | `--hbDepthEdges` | 190.2 214.2 244.8 | Radii [cm] of the HB depth 1/2, 2/3, 3/4 boundaries |

Example: `--energyFrac 0.3 0.5 0.7 --timeCut 2 1 0.5`

## Truth matching

- **LLP:** |pdgId| in `LLP_PDGIDS` with a valid decay vertex. **Decay products:** its non-LLP gen daughters.
- **Match by decay radius R:**

| LLP decay | Cluster is compared to |
|---|---|
| R < 177 cm | Decay product direction, after placing the cluster at R = 177 cm and shifting to the decay vertex frame |
| 177 ≤ R < 295 cm, \|eta_LLP\| ≤ 1.26 | LLP direction, unshifted |
| Other | Not matched |

- **Cluster → truth:** each HB cluster (|eta| < 1.3) goes to its closest candidate if ΔR < cone. Exclusive. Gives matched energy, cluster count, and all per-cluster labels.
- **Truth → cluster:** decay product is matched if any cluster is in the cone. Not exclusive. Gives the efficiency flag.
- **Expected delay** vs a prompt particle at c reaching the same point `x_c` (decay product assumed at c, straight line):

```
dt = L/c (1/beta - 1)  +  (L + |x_c - x_d| - |x_c|)/c
     LLP slowness         path length difference
```

- **Acceptance:** LLP decays before or inside HB, and `x_c` has |eta| < 1.3.

## Plots

"Better" and "worse" refer to timing PF relative to standard PF. Δt is the truth expected delay.

| Page | Plot | Purpose | Better | Worse |
|---|---|---|---|---|
| 1 | Δt distribution, with components | Context: which Δt range is populated | n/a | n/a |
| 2 | min ΔR to an HB cluster, raw vs shifted | Validates the matching (standard PF only) | Expect: shifted curve peaks lower, inside the cone | n/a |
| 3 | Fraction matched vs Δt | LLP cluster efficiency | ≥ standard, flat in Δt | Drops at large Δt |
| 4 | ΣE(matched)/E_b vs Δt | Energy response | ≥ standard, flat in Δt | Falls with Δt |
| 5 | Matched clusters per decay product vs Δt | Splitting or merging | No clear answer; read with page 4 | |
| 6 | (E_timing − E_standard)/E_b vs Δt | Paired energy difference, most sensitive | ≥ 0 | Negative, growing with Δt |
| 7, 8 | Cluster time vs Δt, one per algo | Does cluster time track the delay | Tight band on the dashed line | Band cut off at large Δt, or blob at t ≈ 0 |
| 9 | HB clusters per event | Global clustering change | No clear answer; see 11/12, 14/15 for the source | |
| 10 to 12 | Cluster energy: all, matched, unmatched | Signal vs background effect | Matched unchanged, unmatched reduced | Matched reduced or softer |
| 13 to 15 | Cluster time: all, matched, unmatched | Which times are affected | Matched keeps its delayed tail, unmatched narrows | Matched tail cut off |
| 16, 17 | Decay depth vs cluster depth | Does cluster depth follow the decay depth | Points on the diagonal, \|diff\| peaked at 0 | Flat in decay depth |
| 18 to 21 | Correctly matched per event, 3x3 | Combined energy, time and ΔR quality | Histogram shifted right, higher pooled % | Lower pooled % |

Notes on pages 1 to 15:
- Page 3 saturates near 1 with a 0.4 cone and no energy cut. Label with `--clusterE` > 0 to make it sensitive.
- Page 4: absolute value is not expected to be 1 (non-compensation, ECAL energy, pileup). Only the difference between algos matters.
- Pages 7 and 8: the dashed line is `time = Δt + offset`, slope 1.

### Pages 16 and 17: decay depth vs cluster depth

- One dot per decay product, LLPs decaying inside HB only.
- x: decay R mapped to depth 1 to 4 with `--hbDepthEdges`. y: energy-weighted mean depth of all matched clusters (no energy cut).
- Depth radii (177 / 190.2 / 214.2 / 244.8 / 295 cm) are those used in Run3-HCAL-LLP-Analysis. One set for all eta.
- Bottom panel: |cluster depth − decay depth|.
- Page order: standard, then timing.

### Pages 18 to 21: correctly matched per event

Only matched HB clusters with E > `--clusterE` (2 GeV) are used, including in ΣE. A decay product is correctly matched if it passes all three:

| Cut | Quantity | Default values |
|---|---|---|
| Energy | ΣE(matched clusters)/E_b > cut | 0.4, 0.65, 0.8 |
| Time | mean of \|t − Δt − offset\|/√2 over matched clusters < cut [ns] | 1, 0.5, 0.1 |
| ΔR | mean ΔR of matched clusters < cut | 0.3, 0.15, 0.05 |

- **Layout:** row 1 scans the energy cut, row 2 the time cut, row 3 the ΔR cut. In each row the other two cuts sit at the middle value of their list, so the middle column is the same configuration in all rows.
- **Histogram:** per-event percentage of correctly matched decay products, both algos overlaid.
- **Text in each panel:** pooled percentage per algo (all decay products summed, not averaged over events).
- **Offset:** median of (cluster time − Δt) for standard PF matched clusters (the page 7 selection), used for both algos. Override with `--timeOffset`.
- **Four pages:** denominator (all in HB acceptance, or only those with ≥ 1 matched cluster) × cluster average (energy-weighted, or simple).
- Clusters without a valid time are left out of the time average; a decay product with none fails.

## Caveats

- "Unmatched" means outside the cone of any LLP decay product, not strictly pileup or noise.
- Page 3 uses all HCAL clusters and non-exclusive matching; energy and cluster counts use HB clusters with exclusive matching.
- For LLPs decaying inside HB, ΔR is measured to the LLP direction, not the decay product.
- With the ≥ 1 matched cluster denominator, the two algos have different denominators.
- Pages 3 to 6, 16 and 17 use all matched HB clusters (no 2 GeV cut); pages 18 to 21 apply it.