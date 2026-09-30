#!/usr/bin/env python3
"""
Process one-pion PF-cluster ntuples directly and make eta/energy scan plots.

For every requested 2D quantity, standardPF and seedTimingPF are drawn on the
same landscape PDF page with identical z-axis limits and the same color map,
unchanged from before.  Only the ratio panel's color axis has been switched to
the dipion convention: a diverging magenta-white-green palette centered at
one, with each ratio panel scaled independently from its own values.
The event-level quantities are evaluated both with no HCAL-cluster energy cut
and after requiring E_cluster >= sqrt(0.7 E_pi) GeV, where E_pi is the
generator pion energy encoded in the sample filename.

The script computes, for each cluster selection:

  N_clusters
  N_eff = (sum_i E_i)^2 / sum_i E_i^2
  E_cluster1 / E_pi, with the generator pion energy read event by event from
      gen_energy/gen_pdgId in the ntuple
  E_cluster1 / E_all_PFRecHits, where the denominator is the sum of all HB/HE
      PFRecHits in the event, including unclustered hits
  % events with exactly 1 cluster
  % events with 0 clusters
  % events with >= 2 clusters
  % events with exactly 1 cluster and E_cluster1 / E_all_PFRecHits > 0.8
  % events with >= 2 clusters and E_cluster1 / E_clustered_PFRecHits > 0.8
  % events with >= 2 clusters and E_cluster1 / E_all_PFRecHits > 0.8

The 1D overlay pages include N_clusters, N_eff, the clustered-RecHit
percentage, and reconstructed HCAL energy.  Every 2D scan heatmap includes
the two algorithms plus their ratio.

Finally, density plots combine all pion energies at fixed eta and show
|Delta t| versus E2/E1 for the two leading clusters.  These are produced
twice: once with the PF cluster time, and once with the cluster seed time,
defined as the time of the highest-energy PFRecHit inside each cluster.  Seed
times are whole-nanosecond values, so the seed pages use one bin per integer.
The two leading retained clusters must have finite eta, phi, and the relevant
time, each time must be >= 0 ns (which drops -999 and similar invalid seed
times), and their cluster-cluster separation must satisfy Delta R_12 <= 0.4.
No generator-pion direction enters this cluster-pair selection.

The scan axes are eta and energy, parsed from ntuple names of the form:

  pfObjectsNtuple_<PF>_SinglePiCloseByE20_eta0p1.root
  pfObjectsNtuple_<PF>_SinglePiCloseByE120_eta1.root

The filename-safe tags use ``p`` for the decimal point and ``m`` for a minus
sign. Only CloseByParticleGunProducer samples are processed.
"""

import argparse
import glob
import math
import os
import re
from array import array

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Match the cell-value rounding used in the reference heatmap PDF:
# averages and ratios use two decimal places, while percentages use one.
HEATMAP_VALUE_FORMAT = "4.2f"
HEATMAP_PERCENT_FORMAT = "5.1f"

DEFAULT_PF_NAMES = ["standardPF", "seedTimingPF"]
DEFAULT_SAMPLE_TYPES = ["CloseBy"]
DEFAULT_NCL_MAX = 6.0
DEFAULT_NEFF_MAX = 4.0
DEFAULT_NEFF_BINS = 20
DEFAULT_RECHIT_FRACTION_BINS = 10
# Reconstructed-energy overlays span 0 to 1.5 * max(E_pi), so a modest bin
# count keeps the low-energy samples from being sliced too finely.
DEFAULT_RECO_ENERGY_BINS = 20
DEFAULT_SCATTER_X_BINS = 20
DEFAULT_SCATTER_TIME_BINS = 28
DEFAULT_SCATTER_TIME_MAX = 7.0
# Seed times are integers, so this is the largest integer |Delta t_seed| shown.
DEFAULT_SEED_TIME_MAX = 7.0

# ROOT's ``COLZ TEXT`` labels use marker-size units, whereas the manually
# drawn ratio labels use NDC-like text-size units.  These paired values produce
# the same visual cell-label size in all three heatmaps.
HEATMAP_CELL_MARKER_SIZE = 1.80
HEATMAP_CELL_TEXT_SIZE = 0.036

# Cluster-pair density pages are produced once per time definition.
TIME_KINDS = ["clustertime", "seedtime"]


def cluster_energy_threshold(pion_energy):
    """Return sqrt(0.7 * E_pi), in GeV."""
    pion_energy = float(pion_energy)
    if pion_energy <= 0.0:
        raise ValueError("Generator-pion energy must be positive.")
    return math.sqrt(0.7 * pion_energy)

# -----------------------------------------------------------------------------
# Argument parsing
# -----------------------------------------------------------------------------

def parse_float_list(value, default=None):
    if value is None:
        return default
    if isinstance(value, list):
        return [float(v) for v in value]
    value = str(value).strip()
    if value == "":
        return default
    return [float(x.strip()) for x in value.split(",") if x.strip()]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Loop over one-pion PF-cluster ntuples and make clean eta/energy heatmaps + overlay grids."
    )

    parser.add_argument("--input-dir", "--inputdir", dest="input_dir", default="/eos/user/c/chtong/Public/Rereco/SinglePion_SmallEta_modified_5ns/",
                        help="Directory containing the ntuple ROOT files.")
    parser.add_argument("--output-dir", "--outputdir", dest="output_dir", default="/eos/user/c/chtong/Public/Rereco/SinglePion_SmallEta_modified_5ns/",
                        help="Directory where the output PDF and ROOT file will be saved.")
    parser.add_argument("--prefix", default="pfObjectsNtuple_",
                        help="Filename prefix before the PF algorithm name.")
    parser.add_argument("--tree-name", default="pfObjectsNtupler/pfTree",
                        help="TTree path inside each ntuple.")
    parser.add_argument("--pf-names", nargs="+", default=DEFAULT_PF_NAMES,
                        help="PF variants to compare, e.g. standardPF seedTimingPF.")

    parser.add_argument(
        "--sample-types",
        nargs="+",
        default=DEFAULT_SAMPLE_TYPES,
        help=(
            "Single-pion sample family to process. Only CloseBy "
            "(CloseByParticleGunProducer) is supported."
        ),
    )

    parser.add_argument("--eta-values", default="0.1,0.2,0.4,0.6,0.8,1.0",
                        help="Comma-separated eta scan values for heatmaps and scatter plots. Default: 0.1,0.2,0.4,0.6,0.8,1.0.")
    parser.add_argument("--energy-values", default="20,40,60,80,100,120",
                        help="Comma-separated energy scan values for heatmaps and scatter plots. Default: 20,40,60,80,100,120.")
    parser.add_argument("--eta-dist-values", default="0.1,0.4,0.8",
                        help="Comma-separated eta values for the 3x3 overlay distributions. Default: 0.1,0.4,0.8.")
    parser.add_argument("--energy-dist-values", default="20,40,60",
                        help="Comma-separated energy values for the 3x3 overlay distributions. Default: 20,40,60.")

    parser.add_argument("--energy-branch", default="hcal_energy",
                        help="Per-event vector branch with HCAL cluster energies.")
    parser.add_argument("--cluster-eta-branch", default="hcal_eta",
                        help="Per-event HCAL-cluster eta branch.")
    parser.add_argument("--cluster-phi-branch", default="hcal_phi",
                        help="Per-event HCAL-cluster phi branch.")
    parser.add_argument("--cluster-time-branch", default="hcal_time",
                        help="Per-event HCAL-cluster time branch.")
    parser.add_argument("--all-pfrh-energy-branch", default="all_hbhe_pfrh_energy",
                        help=("Per-event vector branch containing every accepted HB/HE PFRecHit, "
                              "including PFRecHits not assigned to an HCAL cluster."))
    parser.add_argument("--all-pfrh-cluster-index-branch", default="all_hbhe_pfrh_clusterIdx",
                        help=("Cluster-index branch aligned with --all-pfrh-energy-branch; "
                              "-1 denotes an unclustered PFRecHit."))
    parser.add_argument("--clustered-pfrh-energy-branch", default="hbhe_pfrh_energyFracInCluster",
                        help=("PFRecHit energy contribution inside HCAL clusters. This is used "
                              "for the in-cluster-PFRecHit threshold denominator and to pick "
                              "each cluster's seed PFRecHit."))
    parser.add_argument("--clustered-pfrh-cluster-index-branch", default="hbhe_pfrh_clusterIdx",
                        help="Cluster-index branch aligned with --clustered-pfrh-energy-branch.")
    parser.add_argument("--seed-time-branch", default="hbhe_pfrh_time",
                        help=("PFRecHit time branch aligned with "
                              "--clustered-pfrh-energy-branch. The seed time of a cluster is "
                              "the time of its highest-energy PFRecHit. Optional: files "
                              "without this branch are skipped for the seed-time pages only."))
    parser.add_argument("--gen-energy-branch", default="gen_energy",
                        help="Generator-particle energy branch used for event-by-event pion energy.")
    parser.add_argument("--gen-pdgid-branch", default="gen_pdgId",
                        help="Generator-particle PDG-ID branch used to select charged pions.")
    parser.add_argument("--gen-status-branch", default="gen_status",
                        help="Optional generator status branch. When present, status-1 charged pions are preferred.")
    parser.add_argument("--ncl-max", type=float, default=DEFAULT_NCL_MAX,
                        help="Largest integer N_clusters value shown in the overlay plots. Default shows 0 through 6.")
    parser.add_argument("--neff-max", type=float, default=DEFAULT_NEFF_MAX,
                        help="Upper x-axis edge for N_eff overlay histograms. Default is 4.0.")
    parser.add_argument("--neff-bins", type=int, default=DEFAULT_NEFF_BINS,
                        help="Number of bins for N_eff overlay histograms. Default is 20.")
    parser.add_argument("--rechit-fraction-bins", type=int, default=DEFAULT_RECHIT_FRACTION_BINS,
                        help="Number of bins for clustered-RecHit-fraction overlays. Default is 10.")
    parser.add_argument("--reco-energy-bins", type=int, default=DEFAULT_RECO_ENERGY_BINS,
                        help="Number of bins for reconstructed-energy overlays. Default is 20.")
    parser.add_argument("--scatter-x-bins", type=int, default=DEFAULT_SCATTER_X_BINS,
                        help="Number of E2/E1 bins in cluster-pair density plots.")
    parser.add_argument("--scatter-time-bins", type=int, default=DEFAULT_SCATTER_TIME_BINS,
                        help="Number of |Delta t| bins in cluster-time density plots.")
    parser.add_argument("--scatter-time-max", type=float, default=DEFAULT_SCATTER_TIME_MAX,
                        help="Upper |Delta t| edge in ns for cluster-time density plots. Default is 7 ns.")
    parser.add_argument("--seed-time-max", type=float, default=DEFAULT_SEED_TIME_MAX,
                        help=("Largest integer |Delta t_seed| in ns shown in the seed-time "
                              "density plots. Seed times are whole numbers, so the bin count "
                              "is fixed at one bin per integer and is not separately "
                              "configurable. Default is 7 ns."))

    parser.add_argument("--output-pdf", default=None,
                        help="Optional explicit output PDF filename.")
    parser.add_argument("--output-root", default=None,
                        help="Optional explicit output ROOT filename.")
    parser.add_argument("--debug-files", action="store_true",
                        help="Print file matching details.")

    return parser.parse_args()


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def sanitize_name(text):
    return re.sub(r"[^0-9a-zA-Z_]+", "_", str(text)).strip("_")


def key_float(x):
    return round(float(x), 6)


def vector_to_list(v):
    if v is None:
        return []
    if hasattr(v, "size") and hasattr(v, "at"):
        return [v.at(i) for i in range(int(v.size()))]
    try:
        return [v[i] for i in range(len(v))]
    except TypeError:
        try:
            return [v]
        except Exception:
            return []


def clean_positive_energies(values):
    out = []
    for x in values:
        try:
            x = float(x)
        except Exception:
            continue
        if math.isfinite(x) and x > 0.0:
            out.append(x)
    return out


def compute_neff(energies):
    energies = clean_positive_energies(energies)
    if len(energies) == 0:
        return 0.0
    etot = sum(energies)
    sum_e2 = sum(e * e for e in energies)
    if sum_e2 <= 0.0:
        return 0.0
    return (etot * etot) / sum_e2


def branch_exists(tree, branch_name):
    return bool(tree.GetBranch(branch_name))


# -----------------------------------------------------------------------------
# File discovery
# -----------------------------------------------------------------------------

def pf_from_filename(basename, prefix, pf_names):
    if not basename.startswith(prefix):
        return None
    rest = basename[len(prefix):]
    for pf_name in sorted(pf_names, key=len, reverse=True):
        if rest.startswith(pf_name + "_") or rest == pf_name + ".root":
            return pf_name
    return None


def canonical_sample_type(value):
    normalized = str(value).strip().lower()
    aliases = {
        "closeby": "CloseBy",
        "closebyparticlegun": "CloseBy",
        "closebyparticlegunproducer": "CloseBy",
    }
    if normalized not in aliases:
        raise ValueError(
            f"Unknown sample type '{value}'. Supported value: CloseBy (CloseByParticleGunProducer)."
        )
    return aliases[normalized]


def decode_number_tag(tag):
    """Decode the Bash number_tag convention: 0p5 -> 0.5, m0p5 -> -0.5."""
    value = str(tag).strip()
    if value.startswith("m"):
        value = "-" + value[1:]
    value = value.replace("p", ".")
    return float(value)


def parse_singlepi_filename(basename, prefix, pf_names):
    """
    Parse the exact ntuple naming convention produced by the generation script.

    Examples:
      pfObjectsNtuple_standardPF_SinglePiCloseByE20_eta0p1.root
      pfObjectsNtuple_seedTimingPF_SinglePiCloseByE120_eta1.root
    """
    pf_name = pf_from_filename(basename, prefix, pf_names)
    if pf_name is None:
        return None

    expected_prefix = f"{prefix}{pf_name}_"
    payload = basename[len(expected_prefix):]

    match = re.fullmatch(
        r"SinglePi(?P<sample_type>CloseBy)"
        r"E(?P<energy>[-+mp0-9.]+)_eta(?P<eta>[-+mp0-9.]+)\.root",
        payload,
        flags=re.IGNORECASE,
    )
    if match is None:
        return None

    sample_type = canonical_sample_type(match.group("sample_type"))
    try:
        energy = decode_number_tag(match.group("energy"))
        eta = decode_number_tag(match.group("eta"))
    except ValueError:
        return None

    return {
        "pf": pf_name,
        "sample_type": sample_type,
        "eta": eta,
        "energy": energy,
    }


def collect_file_metadata(input_dir, prefix, pf_names, sample_types, debug=False):
    pattern = os.path.join(input_dir, f"{prefix}*.root")
    candidates = sorted(glob.glob(pattern))
    requested_types = set(sample_types)

    if debug:
        print(f"Scanning for ntuples with pattern: {pattern}")
        print(f"Found {len(candidates)} candidate ROOT files")

    metadata = []
    for path in candidates:
        base = os.path.basename(path)
        parsed = parse_singlepi_filename(base, prefix, pf_names)
        if parsed is None:
            if debug:
                print(f"  Could not parse supported single-pion filename: {base}")
            continue
        if parsed["sample_type"] not in requested_types:
            continue

        parsed["path"] = path
        metadata.append(parsed)
        if debug:
            print(
                f"  Matched {base}: pf={parsed['pf']}, "
                f"type={parsed['sample_type']}, eta={parsed['eta']}, "
                f"E={parsed['energy']}"
            )
    return metadata


def select_scan_values(metadata, requested_etas, requested_energies):
    available_etas = sorted({key_float(m["eta"]) for m in metadata})
    available_energies = sorted({key_float(m["energy"]) for m in metadata})

    eta_values = (
        available_etas
        if requested_etas is None
        else [key_float(x) for x in requested_etas]
    )
    energy_values = (
        available_energies
        if requested_energies is None
        else [key_float(x) for x in requested_energies]
    )

    return eta_values, energy_values, available_etas, available_energies


def build_file_lookup(metadata, eta_values, energy_values, debug=False):
    lookup = {}
    eta_set = {key_float(x) for x in eta_values}
    energy_set = {key_float(x) for x in energy_values}

    for item in metadata:
        eta_key = key_float(item["eta"])
        energy_key = key_float(item["energy"])
        if eta_key not in eta_set or energy_key not in energy_set:
            continue

        key = (item["pf"], eta_key, energy_key)
        if key in lookup:
            print(
                "WARNING: duplicate file for "
                f"pf={item['pf']}, type={item['sample_type']}, "
                f"eta={item['eta']}, E={item['energy']}. Keeping first:\n"
                f"  first: {lookup[key]}\n"
                f"  extra: {item['path']}"
            )
            continue
        lookup[key] = item["path"]

    return lookup


# -----------------------------------------------------------------------------
# ROOT drawing helpers
# -----------------------------------------------------------------------------

_ALGORITHM_PALETTE_EXEC_COMMAND = None
_RATIO_PALETTE_EXEC_COMMAND = None
_DECLARED_PALETTE_GLOBALS = set()

# Original algorithm-panel gradient, unchanged.
ALGORITHM_PALETTE_GRADIENT = (
    [0.00, 0.50, 1.00],
    [0.20, 1.00, 0.95],
    [0.55, 1.00, 0.55],
    [0.80, 1.00, 0.20],
)

# Ratio-panel gradient: magenta at low, pure white at one, green at high.
# ROOT's kRedBlue has a dark desaturated midpoint, which turns a map of
# near-unity ratios into flat olive sludge; a genuinely white center keeps
# small deviations readable, and magenta/green stays clear of the blue-white
# gradient the algorithm panels use.
RATIO_PALETTE_GRADIENT = (
    [0.00, 0.25, 0.50, 0.75, 1.00],
    [0.77, 0.96, 1.00, 0.60, 0.30],
    [0.11, 0.55, 1.00, 0.81, 0.57],
    [0.49, 0.75, 1.00, 0.40, 0.13],
)


def _register_gradient_palette(global_name, stops, red, green, blue):
    """Create a gradient color table and return a TExec command restoring it.

    The color table is created exactly once.  Its color indices are captured in
    an interpreter global so a pad's TExec can put the palette back on every
    repaint without allocating a fresh table each time.
    """
    if global_name in _DECLARED_PALETTE_GLOBALS:
        return f"gStyle->SetPalette(255, {global_name});"
    first_index = int(
        ROOT.TColor.CreateGradientColorTable(
            len(stops),
            array("d", stops),
            array("d", red),
            array("d", green),
            array("d", blue),
            255,
        )
    )
    ROOT.gStyle.SetNumberContours(255)
    values = ",".join(str(first_index + offset) for offset in range(255))
    try:
        ROOT.gInterpreter.Declare(f"int {global_name}[255] = {{{values}}};")
    except Exception:
        print(
            f"WARNING: could not register palette '{global_name}' with the ROOT "
            "interpreter; falling back to a stock palette."
        )
        return None
    _DECLARED_PALETTE_GLOBALS.add(global_name)
    return f"gStyle->SetPalette(255, {global_name});"


def install_heatmap_palettes():
    """Build both heatmap palettes once, before any page is drawn.

    The algorithm panels and the ratio panel want different palettes, so gStyle
    cannot simply be set once and left alone: each pad restores its own palette
    whenever ROOT repaints the canvas.
    """
    global _ALGORITHM_PALETTE_EXEC_COMMAND, _RATIO_PALETTE_EXEC_COMMAND
    _RATIO_PALETTE_EXEC_COMMAND = _register_gradient_palette(
        "gOnePiRatioPalette", *RATIO_PALETTE_GRADIENT
    )
    _ALGORITHM_PALETTE_EXEC_COMMAND = _register_gradient_palette(
        "gOnePiAlgorithmPalette", *ALGORITHM_PALETTE_GRADIENT
    )


def palette_exec(name, is_ratio):
    """Apply a pad-local palette whenever ROOT repaints the canvas.

    The two algorithm panels keep the original blue-white-red gradient; the
    ratio panel uses the magenta-white-green gradient above.  Attaching
    both per pad with TExec keeps them correct even when ROOT repaints the
    whole multipad canvas at PDF-output time, and makes the heatmaps immune to
    the gStyle palette the density pages install.
    """
    if is_ratio:
        command = (
            _RATIO_PALETTE_EXEC_COMMAND
            or f"gStyle->SetPalette({int(ROOT.kRedBlue)});"
        )
    else:
        command = (
            _ALGORITHM_PALETTE_EXEC_COMMAND
            or f"gStyle->SetPalette({int(ROOT.kBird)});"
        )
    return ROOT.TExec(name, command)


def set_density_palette(show_zero_as_white=False):
    """Use a perceptually ordered palette for nonnegative event counts.

    The scan heatmaps use a diverging palette on their ratio panels because
    those are centered around one.  Cluster-pair densities are counts, so
    reusing that palette incorrectly suggests negative/positive deviations.
    For sparse, linearly scaled pages, reserve the zero-count end of the
    palette as pure white; every positive count then receives a clearly
    visible Viridis-like color.  Log-scaled pages start directly at dark
    violet because zero lies below the displayed range.
    """
    if show_zero_as_white:
        stops = array("d", [0.00, 0.015, 0.25, 0.50, 0.75, 1.00])
        red = array("d", [1.000, 0.267, 0.230, 0.128, 0.369, 0.993])
        green = array("d", [1.000, 0.005, 0.322, 0.567, 0.789, 0.906])
        blue = array("d", [1.000, 0.329, 0.546, 0.551, 0.383, 0.144])
    else:
        stops = array("d", [0.00, 0.25, 0.50, 0.75, 1.00])
        red = array("d", [0.267, 0.230, 0.128, 0.369, 0.993])
        green = array("d", [0.005, 0.322, 0.567, 0.789, 0.906])
        blue = array("d", [0.329, 0.546, 0.551, 0.383, 0.144])
    ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, 255)
    ROOT.gStyle.SetNumberContours(255)


def style_for_pf(pf_name, index=0):
    if pf_name == "standardPF":
        return {"color": ROOT.kAzure + 1, "style": 1, "width": 1, "fill": True, "alpha": 0.22}
    if pf_name == "seedTimingPF":
        return {"color": ROOT.kBlack, "style": 1, "width": 1, "fill": False, "alpha": 0.0}
    colors = [ROOT.kAzure + 1, ROOT.kBlack, ROOT.kGreen + 3, ROOT.kMagenta + 1]
    return {"color": colors[index % len(colors)], "style": 1 + index, "width": 1, "fill": False, "alpha": 0.0}


def apply_pf_style(h, pf_name, index=0):
    style = style_for_pf(pf_name, index)
    h.SetLineColor(style["color"])
    h.SetLineStyle(style["style"])
    h.SetLineWidth(style["width"])
    if style["fill"]:
        h.SetFillColorAlpha(style["color"], style["alpha"])
    else:
        h.SetFillStyle(0)
    return h


def make_scan_hist(name, title, eta_values, energy_values):
    h = ROOT.TH2F(name, title, len(eta_values), 0, len(eta_values), len(energy_values), 0, len(energy_values))
    for i, eta in enumerate(eta_values, start=1):
        h.GetXaxis().SetBinLabel(i, str(eta))
    for j, energy in enumerate(energy_values, start=1):
        h.GetYaxis().SetBinLabel(j, str(energy))
    return h


def draw_single_heatmap(hist, title, z_title, x_title, y_title, text_format, zmin, zmax,
                        draw_text=True, palette_name=None, is_ratio=False):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat(text_format)

    if zmin is not None:
        hist.SetMinimum(zmin)
    if zmax is not None:
        hist.SetMaximum(zmax)

    hist.SetTitle(title)
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetZaxis().SetTitle(z_title)

    hist.GetXaxis().CenterTitle()
    hist.GetYaxis().CenterTitle()
    hist.GetZaxis().CenterTitle()

    # Match the readable sizing used by the two algorithm heatmaps in the
    # reference pages.  The ratio panel uses the same axis sizing below.
    hist.GetXaxis().SetTitleSize(0.050)
    hist.GetYaxis().SetTitleSize(0.050)
    hist.GetZaxis().SetTitleSize(0.045)

    hist.GetXaxis().SetLabelSize(0.043)
    hist.GetYaxis().SetLabelSize(0.043)
    hist.GetZaxis().SetLabelSize(0.038)

    hist.GetXaxis().SetTitleOffset(1.05)
    hist.GetYaxis().SetTitleOffset(1.10)
    hist.GetZaxis().SetTitleOffset(1.16)
    hist.GetXaxis().SetNdivisions(506)
    hist.GetYaxis().SetNdivisions(506)
    hist.GetZaxis().SetNdivisions(506)

    # SetMarkerSize controls the numeric labels produced by ROOT's TEXT option.
    # Keep the values legible without letting them dominate small heatmap cells.
    hist.SetMarkerSize(HEATMAP_CELL_MARKER_SIZE)

    # ROOT stores this setting globally, so reset it immediately before every
    # heatmap draw. This guarantees identical rounding for no-cut and
    # sqrt(0.7 E_pi)-cut pages in every sample family.
    ROOT.gStyle.SetPaintTextFormat(text_format)

    # Draw the frame first so the TExec palette is registered in the pad's
    # primitive list before the colored contents are painted.
    hist.Draw("AXIS")
    exec_object = None
    if palette_name is not None:
        exec_object = palette_exec(palette_name, is_ratio)
        exec_object.Draw()
    hist.Draw("COLZ TEXT SAME" if draw_text else "COLZ SAME")
    return exec_object


def make_heatmap_ratio(numerator, denominator, name):
    """Return numerator/denominator and the set of bins with valid denominators."""
    ratio = numerator.Clone(name)
    ratio.Reset("ICES")
    ratio.SetDirectory(0)
    valid_bins = set()
    for ix in range(1, denominator.GetNbinsX() + 1):
        for iy in range(1, denominator.GetNbinsY() + 1):
            den = denominator.GetBinContent(ix, iy)
            num = numerator.GetBinContent(ix, iy)
            if not math.isfinite(den) or not math.isfinite(num) or den == 0.0:
                ratio.SetBinContent(ix, iy, float("nan"))
                continue
            ratio.SetBinContent(ix, iy, num / den)
            valid_bins.add((ix, iy))
    return ratio, valid_bins


def build_ratio_heatmaps(results_by_sample_type, pf_names):
    """Build every ratio map once, each with its own range centered at one.

    Every ratio panel is scaled independently from its own values, so a metric
    whose ratios sit within a percent of unity still shows real structure even
    when some other metric on some other page has a near-zero denominator.  The
    cost is that two ratio panels are no longer directly comparable by color,
    so read the printed range (and the cell labels) before comparing pages.
    """
    if len(pf_names) < 2:
        raise RuntimeError("Ratio heatmaps require at least two PF algorithms.")
    if "standardPF" in pf_names and "seedTimingPF" in pf_names:
        pf_left, pf_right = "standardPF", "seedTimingPF"
    else:
        pf_left, pf_right = pf_names[0], pf_names[1]

    saturating = []
    for sample_type, result in results_by_sample_type.items():
        ratios_by_mode = {}
        for mode, heatmaps in result["heatmaps_by_mode"].items():
            ratios = {}
            for metric_key in heatmaps[pf_left]:
                name = "h2_ratio_" + sanitize_name(
                    f"{sample_type}_{metric_key}_{mode}"
                )
                ratio, valid_bins = make_heatmap_ratio(
                    heatmaps[pf_right][metric_key],
                    heatmaps[pf_left][metric_key],
                    name,
                )
                values = []
                for ix, iy in valid_bins:
                    value = ratio.GetBinContent(ix, iy)
                    if math.isfinite(value) and value > 0.0:
                        values.append(value)
                # The 0.05 floor stops an all-identical panel from being
                # scaled down to pure numerical noise.
                raw_deviation = max(
                    [abs(value - 1.0) for value in values] + [0.05]
                )
                # A ratio scale wider than [0, 2] can no longer be symmetric
                # about one while remaining physical.  Saturate rarer extremes
                # at the endpoints.
                deviation = min(1.08 * raw_deviation, 1.0)
                if raw_deviation > 1.0:
                    saturating.append(f"{sample_type}/{mode}/{metric_key}")
                ratios[metric_key] = (
                    ratio,
                    valid_bins,
                    1.0 - deviation,
                    1.0 + deviation,
                )
            ratios_by_mode[mode] = ratios
        result["ratio_heatmaps_by_mode"] = ratios_by_mode

    if saturating:
        print(
            "WARNING: ratio values outside [0, 2] saturate the color scale on "
            f"{len(saturating)} panel(s): " + ", ".join(saturating)
        )


def draw_ratio_heatmap(hist, valid_bins, title, pdf_z_title, text_format, zmin, zmax,
                       palette_name):
    drawn_objects = []
    drawn_objects.append(
        draw_single_heatmap(
            hist,
            title,
            pdf_z_title,
            "#eta",
            "Pion energy [GeV]",
            text_format,
            zmin,
            zmax,
            draw_text=False,
            palette_name=palette_name,
            is_ratio=True,
        )
    )
    # ROOT maps NaN/invalid bins to the low end of the palette on some
    # versions.  Explicitly mask zero-denominator cells in white so they are
    # blank rather than appearing as saturated ratio values.
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            if (ix, iy) in valid_bins:
                continue
            box = ROOT.TBox(
                hist.GetXaxis().GetBinLowEdge(ix),
                hist.GetYaxis().GetBinLowEdge(iy),
                hist.GetXaxis().GetBinUpEdge(ix),
                hist.GetYaxis().GetBinUpEdge(iy),
            )
            box.SetFillColor(ROOT.kWhite)
            box.SetLineColor(ROOT.kWhite)
            box.Draw("same")
            drawn_objects.append(box)

    for ix, iy in sorted(valid_bins):
        value = hist.GetBinContent(ix, iy)
        if not math.isfinite(value):
            continue
        label = ROOT.TLatex()
        label.SetTextAlign(22)
        label.SetTextFont(42)
        # Keep ratio values visually consistent with the values in the two
        # algorithm heatmaps, without returning to the oversized labels.
        label.SetTextSize(HEATMAP_CELL_TEXT_SIZE)
        label.DrawLatex(hist.GetXaxis().GetBinCenter(ix), hist.GetYaxis().GetBinCenter(iy), f"{value:.2f}")
        drawn_objects.append(label)
    hist.Draw("AXIS SAME")
    return drawn_objects


def make_paired_heatmap(
    heatmaps,
    ratio_entry,
    pf_names,
    metric_key,
    title,
    z_title,
    pdf_path,
    cut_label,
    sample_label,
    text_format=HEATMAP_VALUE_FORMAT,
    zmin=None,
    zmax=None,
):
    if len(pf_names) < 2:
        raise RuntimeError("Paired heatmaps require at least two PF algorithms.")

    if "standardPF" in pf_names and "seedTimingPF" in pf_names:
        pf_left, pf_right = "standardPF", "seedTimingPF"
    else:
        pf_left, pf_right = pf_names[0], pf_names[1]
    h_left = heatmaps.get(pf_left, {}).get(metric_key)
    h_right = heatmaps.get(pf_right, {}).get(metric_key)
    if h_left is None or h_right is None or ratio_entry is None:
        print(f"WARNING: missing heatmap '{metric_key}' for paired drawing")
        return

    # Each ratio panel carries its own range, computed from its own values.
    ratio, valid_ratio_bins, ratio_zmin, ratio_zmax = ratio_entry

    object_tag = sanitize_name(f"{sample_label}_{metric_key}_{cut_label}")
    keep = []
    c = ROOT.TCanvas(f"c_pair_{object_tag}", "", 2200, 900)

    # Dedicated three-line header: sample family, metric, and cut selection.
    header = ROOT.TPad(f"header_pair_{object_tag}", "", 0.0, 0.86, 1.0, 1.0)
    header.SetFillStyle(0)
    header.Draw()
    header.cd()

    sample_text = ROOT.TLatex(0.5, 0.78, sample_label)
    sample_text.SetNDC()
    sample_text.SetTextAlign(22)
    sample_text.SetTextFont(62)
    sample_text.SetTextSize(0.22)
    sample_text.Draw()

    title_text = ROOT.TLatex(0.5, 0.48, title)
    title_text.SetNDC()
    title_text.SetTextAlign(22)
    title_text.SetTextFont(62)
    title_text.SetTextSize(0.22)
    title_text.Draw()

    cut_text = ROOT.TLatex(0.5, 0.17, cut_label)
    cut_text.SetNDC()
    cut_text.SetTextAlign(22)
    cut_text.SetTextFont(42)
    cut_text.SetTextSize(0.18)
    cut_text.Draw()
    keep.extend([header, sample_text, title_text, cut_text])

    c.cd()
    body = ROOT.TPad(f"body_pair_{object_tag}", "", 0.0, 0.0, 1.0, 0.86)
    body.SetFillStyle(0)
    body.Draw()
    body.cd()
    body.Divide(3, 1, 0.002, 0.002)
    keep.append(body)

    for ipad, (pf_name, hist) in enumerate(
        [(pf_left, h_left), (pf_right, h_right)], start=1
    ):
        pad = body.cd(ipad)
        pad.SetLeftMargin(0.11)
        pad.SetRightMargin(0.17)
        pad.SetBottomMargin(0.12)
        pad.SetTopMargin(0.10)
        pad.SetTicks(1, 1)
        keep.append(
            draw_single_heatmap(
                hist,
                pf_name,
                z_title,
                "#eta",
                "Pion energy [GeV]",
                text_format,
                zmin,
                zmax,
                palette_name=f"palette_{object_tag}_{ipad}",
                is_ratio=False,
            )
        )

    ratio_pad = body.cd(3)
    ratio_pad.SetLeftMargin(0.11)
    ratio_pad.SetRightMargin(0.17)
    ratio_pad.SetBottomMargin(0.12)
    ratio_pad.SetTopMargin(0.10)
    ratio_pad.SetTicks(1, 1)
    keep.append(
        draw_ratio_heatmap(
            ratio,
            valid_ratio_bins,
            f"{pf_right} / {pf_left}",
            "Ratio",
            HEATMAP_VALUE_FORMAT,
            ratio_zmin,
            ratio_zmax,
            palette_name=f"palette_{object_tag}_ratio",
        )
    )

    c.cd()
    c.Modified()
    c.Update()
    c.Print(pdf_path)
    print(
        f"Added heatmap pair: {sample_label}: {title} [{cut_label}] "
        f"(ratio range [{ratio_zmin:.3f}, {ratio_zmax:.3f}])"
    )


def clone_and_normalize(h, normalize=True):
    hn = h.Clone(h.GetName() + "_draw")
    hn.SetDirectory(0)
    if normalize:
        integral = hn.Integral(0, hn.GetNbinsX() + 1)
        if integral > 0:
            hn.Scale(1.0 / integral)
    return hn


def last_nonzero_x(h):
    last = 1
    for b in range(1, h.GetNbinsX() + 1):
        if h.GetBinContent(b) > 0:
            last = b
    return h.GetXaxis().GetBinUpEdge(last)


def make_overlay_distribution_grid(dist_hists, pf_names, pdf_path, eta_values, energy_values,
                                   xaxis_title, grid_title, tag, normalize=True,
                                   ideal_line=None, x_range=None):
    ncols = len(eta_values)
    nrows = len(energy_values)
    pf_tag = "_".join(sanitize_name(pf) for pf in pf_names)

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas(f"c_grid_{sanitize_name(tag)}_{pf_tag}", "", 680 * ncols, 620 * nrows)
    keep = []

    title_pad = ROOT.TPad(f"title_{sanitize_name(tag)}_{pf_tag}", "", 0.0, 0.955, 1.0, 1.0)
    title_pad.SetFillStyle(0)
    title_pad.Draw()
    title_pad.cd()
    title = ROOT.TLatex(0.5, 0.45, grid_title)
    title.SetNDC()
    title.SetTextAlign(22)
    title.SetTextSize(0.40)
    title.SetTextFont(62)
    title.Draw()
    keep.extend([title_pad, title])

    c.cd()
    main_pad = ROOT.TPad(f"main_{sanitize_name(tag)}_{pf_tag}", "", 0.0, 0.0, 1.0, 0.955)
    main_pad.SetFillStyle(0)
    main_pad.Draw()
    main_pad.cd()
    main_pad.Divide(ncols, nrows, 0.002, 0.002)
    keep.append(main_pad)

    draw_hists = {}
    global_ymax = 0.0
    global_xlo = None
    global_xhi = 0.0
    for energy in energy_values:
        for eta in eta_values:
            point = (key_float(eta), key_float(energy))
            for ipf, pf_name in enumerate(pf_names):
                raw = dist_hists.get((pf_name, point[0], point[1]))
                if raw is None:
                    continue
                hist = clone_and_normalize(raw, normalize=normalize)
                apply_pf_style(hist, pf_name, ipf)
                draw_hists[(pf_name, point[0], point[1])] = hist
                keep.append(hist)
                global_ymax = max(global_ymax, hist.GetMaximum())
                global_xlo = hist.GetXaxis().GetXmin() if global_xlo is None else min(global_xlo, hist.GetXaxis().GetXmin())
                global_xhi = max(global_xhi, last_nonzero_x(hist))

    if x_range is not None:
        global_xlo, global_xhi = x_range
    else:
        global_xlo = 0.0 if global_xlo is None else global_xlo
        global_xhi = global_xlo + 1.0 if global_xhi <= global_xlo else global_xhi
    global_ymax = 1.0 if global_ymax <= 0.0 else global_ymax
    ylabel = "fraction of events" if normalize else "# of events"

    for row, energy in enumerate(energy_values):
        for col, eta in enumerate(eta_values):
            pad_index = row * ncols + col + 1
            cell = main_pad.cd(pad_index)
            cell.SetFillStyle(0)
            cell.SetLeftMargin(0.13)
            cell.SetRightMargin(0.035)
            cell.SetTopMargin(0.10)
            cell.SetBottomMargin(0.13)
            cell.SetTicks(1, 1)
            point = (key_float(eta), key_float(energy))
            present = [
                draw_hists.get((pf_name, point[0], point[1]))
                for pf_name in pf_names
            ]
            if not any(hist is not None for hist in present):
                missing = ROOT.TLatex(0.5, 0.55, f"missing: #eta={eta}, E={energy} GeV")
                missing.SetNDC()
                missing.SetTextAlign(22)
                missing.SetTextColor(ROOT.kGray + 2)
                missing.Draw()
                keep.append(missing)
                continue

            first = True
            for pf_name in pf_names:
                hist = draw_hists.get((pf_name, point[0], point[1]))
                if hist is None:
                    continue
                hist.SetTitle(f"#eta = {eta}, E = {energy} GeV")
                hist.GetYaxis().SetTitle(ylabel)
                hist.GetYaxis().SetTitleSize(0.055)
                hist.GetYaxis().SetLabelSize(0.050)
                hist.GetYaxis().SetTitleOffset(1.05)
                hist.GetXaxis().SetTitle(xaxis_title)
                hist.GetXaxis().SetTitleSize(0.055)
                hist.GetXaxis().SetLabelSize(0.050)
                hist.GetXaxis().SetTitleOffset(1.05)
                hist.GetXaxis().SetRangeUser(global_xlo, global_xhi)
                hist.SetMinimum(0.0)
                hist.SetMaximum(1.28 * global_ymax)
                hist.Draw("HIST" if first else "HIST SAME")
                first = False

            if ideal_line is not None and global_xlo <= ideal_line <= global_xhi:
                line = ROOT.TLine(ideal_line, 0.0, ideal_line, 1.28 * global_ymax)
                line.SetLineColor(ROOT.kGray + 2)
                line.SetLineStyle(2)
                line.Draw()
                keep.append(line)

            means = ROOT.TLatex()
            means.SetNDC()
            means.SetTextSize(0.040)
            y = 0.82
            for ipf, pf_name in enumerate(pf_names):
                raw = dist_hists.get((pf_name, point[0], point[1]))
                if raw is None:
                    continue
                means.SetTextColor(style_for_pf(pf_name, ipf)["color"])
                means.DrawLatex(0.52, y, f"{pf_name}: mean={raw.GetMean():.2f}")
                y -= 0.060
            keep.append(means)

            if pad_index == 1:
                leg = ROOT.TLegend(0.52, 0.56, 0.91, 0.70)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetTextSize(0.040)
                for ipf, pf_name in enumerate(pf_names):
                    hist = draw_hists.get((pf_name, point[0], point[1]))
                    if hist is not None:
                        leg.AddEntry(hist, pf_name, "lf" if style_for_pf(pf_name, ipf)["fill"] else "l")
                leg.Draw()
                keep.append(leg)

    c.Print(pdf_path)
    print(f"Added overlay grid: {grid_title}")


# -----------------------------------------------------------------------------
# Ntuple processing
# -----------------------------------------------------------------------------

def book_distribution_hists(args, pf_names, eta_values, energy_values, sample_type):
    hists = {
        "ncl": {},
        "neff": {},
        "clustered_rh_pct": {},
        "reco_energy": {},
        "ncl_ecut": {},
        "neff_ecut": {},
        "clustered_rh_pct_ecut": {},
        "reco_energy_ecut": {},
    }
    reco_energy_max = 1.5 * max(energy_values) if energy_values else 100.0
    for pf_name in pf_names:
        for eta in eta_values:
            for energy in energy_values:
                eta_key = key_float(eta)
                energy_key = key_float(energy)
                tag = sanitize_name(f"{sample_type}_{pf_name}_eta{eta}_E{energy}")

                ncl_display_max = int(round(args.ncl_max))
                threshold = cluster_energy_threshold(energy)
                for group, suffix, label in [
                    ("ncl", "nocut", "no cluster-energy cut"),
                    ("ncl_ecut", "ecut", f"E_cluster >= sqrt(0.7 E_pi) = {threshold:.2f} GeV"),
                ]:
                    hist = ROOT.TH1D(
                        f"h_ncl_{tag}_{suffix}",
                        f"{pf_name} eta={eta} E={energy} ({label});N_{{clusters}};events",
                        ncl_display_max + 1,
                        -0.5,
                        ncl_display_max + 0.5,
                    )
                    for ibin in range(1, ncl_display_max + 2):
                        hist.GetXaxis().SetBinLabel(ibin, str(ibin - 1))
                    hists[group][(pf_name, eta_key, energy_key)] = hist

                hists["neff"][(pf_name, eta_key, energy_key)] = ROOT.TH1D(
                    f"h_neff_{tag}_nocut",
                    f"{pf_name} eta={eta} E={energy} (no cluster-energy cut);N_{{eff}};events",
                    args.neff_bins,
                    0.0,
                    args.neff_max,
                )
                hists["neff_ecut"][(pf_name, eta_key, energy_key)] = ROOT.TH1D(
                    f"h_neff_{tag}_ecut",
                    f"{pf_name} eta={eta} E={energy} (E_cluster >= sqrt(0.7 E_pi) = {threshold:.2f} GeV);N_{{eff}};events",
                    args.neff_bins,
                    0.0,
                    args.neff_max,
                )

                for group, suffix, label in [
                    ("clustered_rh_pct", "nocut", "no cluster-energy cut"),
                    ("clustered_rh_pct_ecut", "ecut", f"E_cluster >= sqrt(0.7 E_pi) = {threshold:.2f} GeV"),
                ]:
                    hists[group][(pf_name, eta_key, energy_key)] = ROOT.TH1D(
                        f"h_clustered_rh_pct_{tag}_{suffix}",
                        f"{pf_name} eta={eta} E={energy} ({label});clustered PFRecHits [%];events",
                        args.rechit_fraction_bins,
                        0.0,
                        100.001,
                    )

                for group, suffix, label in [
                    ("reco_energy", "nocut", "no cluster-energy cut"),
                    ("reco_energy_ecut", "ecut", f"E_cluster >= sqrt(0.7 E_pi) = {threshold:.2f} GeV"),
                ]:
                    hists[group][(pf_name, eta_key, energy_key)] = ROOT.TH1D(
                        f"h_reco_energy_{tag}_{suffix}",
                        f"{pf_name} eta={eta} E={energy} ({label});sum E_{{cluster}} [GeV];events",
                        args.reco_energy_bins,
                        0.0,
                        reco_energy_max,
                    )
    return hists


def empty_stats():
    return {
        "n_events": 0,
        "sum_ncl": 0.0,
        "sum_neff": 0.0,
        "sum_top1_over_pion": 0.0,
        "sum_top1_over_all_pfrh": 0.0,
        "sum_reco_energy": 0.0,
        "sum_clustered_rh_pct": 0.0,
        "n_valid_clustered_rh_pct": 0,
        "n_zero": 0,
        "n_eq1": 0,
        "n_ge2": 0,
        "n_eq1_top1_all_pfrh_gt80": 0,
        "n_ge2_top1_clustered_pfrh_gt80": 0,
        "n_ge2_top1_all_pfrh_gt80": 0,
        "n_missing_pion": 0,
        "n_zero_all_pfrh": 0,
        "n_zero_clustered_pfrh": 0,
    }


def generator_single_pion_energy(event, args):
    energies = vector_to_list(getattr(event, args.gen_energy_branch))
    pdg_ids = vector_to_list(getattr(event, args.gen_pdgid_branch))
    statuses = []
    if hasattr(event, args.gen_status_branch):
        statuses = vector_to_list(getattr(event, args.gen_status_branch))

    n = min(len(energies), len(pdg_ids))
    final_state = []
    all_pions = []
    for i in range(n):
        try:
            energy = float(energies[i])
            pdg_id = int(pdg_ids[i])
        except Exception:
            continue
        if not math.isfinite(energy) or energy <= 0.0 or abs(pdg_id) != 211:
            continue
        all_pions.append(energy)
        if i < len(statuses):
            try:
                if int(statuses[i]) == 1:
                    final_state.append(energy)
            except Exception:
                pass

    pion_energies = final_state if len(final_state) >= 1 else all_pions
    pion_energies = sorted(pion_energies, reverse=True)
    if len(pion_energies) < 1:
        return 0.0
    return pion_energies[0]


def delta_r(eta1, phi1, eta2, phi2):
    dphi = math.atan2(math.sin(phi1 - phi2), math.cos(phi1 - phi2))
    return math.hypot(eta1 - eta2, dphi)


def selected_cluster_quantities(raw_cluster_energies, threshold=None):
    selected = []
    selected_indices = set()
    for index, value in enumerate(raw_cluster_energies):
        try:
            energy = float(value)
        except Exception:
            continue
        if not math.isfinite(energy) or energy <= 0.0:
            continue
        if threshold is not None and energy < threshold:
            continue
        selected.append(energy)
        selected_indices.add(index)

    multiplicity = len(selected)
    return selected, selected_indices, multiplicity


def sum_clustered_pfrh_for_selected_clusters(pfrh_energy_fractions, pfrh_cluster_indices, selected_cluster_indices):
    total = 0.0
    n = min(len(pfrh_energy_fractions), len(pfrh_cluster_indices))
    for i in range(n):
        try:
            cluster_index = int(pfrh_cluster_indices[i])
            energy = float(pfrh_energy_fractions[i])
        except Exception:
            continue
        if cluster_index not in selected_cluster_indices:
            continue
        if math.isfinite(energy) and energy > 0.0:
            total += energy
    return total


def clustered_rechit_percentage(all_pfrh_cluster_indices, selected_cluster_indices):
    """Percentage of unique event-level HCAL PFRecHits assigned to retained clusters."""
    if len(all_pfrh_cluster_indices) == 0:
        return None
    clustered = 0
    for value in all_pfrh_cluster_indices:
        try:
            cluster_index = int(value)
        except Exception:
            continue
        if cluster_index in selected_cluster_indices:
            clustered += 1
    return 100.0 * clustered / len(all_pfrh_cluster_indices)


def cluster_seed_times(pfrh_energy_fractions, pfrh_cluster_indices, pfrh_times):
    """Map each cluster index to the time of its highest-energy PFRecHit.

    The three input vectors are the in-cluster PFRecHit arrays and must be
    aligned entry by entry.  The seed is chosen purely by energy so that the
    reported value really is the seed PFRecHit's time; invalid sentinel times
    such as -999 are left in place here and rejected later by the same
    ``time >= 0`` requirement used for cluster times.
    """
    best = {}
    n = min(len(pfrh_energy_fractions), len(pfrh_cluster_indices), len(pfrh_times))
    for i in range(n):
        try:
            cluster_index = int(pfrh_cluster_indices[i])
            energy = float(pfrh_energy_fractions[i])
        except Exception:
            continue
        if cluster_index < 0:
            continue
        if not math.isfinite(energy) or energy <= 0.0:
            continue
        try:
            time = float(pfrh_times[i])
        except Exception:
            time = float("nan")
        current = best.get(cluster_index)
        if current is None or energy > current[0]:
            best[cluster_index] = (energy, time)
    return {index: value[1] for index, value in best.items()}


def indexed_time_lookup(values):
    """Return a per-cluster-index time lookup backed by a positional vector."""
    def lookup(index):
        if index < 0 or index >= len(values):
            return None
        return values[index]
    return lookup


def mapped_time_lookup(times_by_index):
    """Return a per-cluster-index time lookup backed by a dict."""
    def lookup(index):
        return times_by_index.get(index)
    return lookup


def fill_cluster_pair_density(hist, raw_energies, raw_etas, raw_phis, time_lookup,
                              selected_cluster_indices, max_pair_delta_r=0.4):
    """
    Fill |t2-t1| versus E2/E1 for the two leading retained clusters.

    The two clusters are chosen by energy before the timing requirement is
    applied.  Both must have finite eta, phi, and time, both times must be
    nonnegative, and their cluster-cluster separation must satisfy
    Delta R_12 <= ``max_pair_delta_r``.  Generator-pion coordinates are not
    used.  ``time_lookup`` maps a cluster index to that cluster's time, so the
    same selection serves both the PF cluster time and the seed time.
    """
    energy_records = []
    for index in selected_cluster_indices:
        if index >= len(raw_energies):
            continue
        try:
            energy = float(raw_energies[index])
        except Exception:
            continue
        if not math.isfinite(energy) or energy <= 0.0:
            continue
        energy_records.append((energy, index))

    if len(energy_records) < 2:
        return False
    # Descending ordering explicitly defines E1 as the higher-energy cluster
    # and E2 as the lower-energy (second-leading) cluster, so E2 / E1 <= 1.
    energy_records.sort(key=lambda item: item[0], reverse=True)
    (e1, index1), (e2, index2) = energy_records[:2]

    cluster_records = []
    for energy, index in ((e1, index1), (e2, index2)):
        if index >= len(raw_etas) or index >= len(raw_phis):
            return False
        if index not in selected_cluster_indices:
            return False
        raw_time = time_lookup(index)
        if raw_time is None:
            return False
        try:
            eta = float(raw_etas[index])
            phi = float(raw_phis[index])
            time = float(raw_time)
        except Exception:
            return False
        if not all(math.isfinite(value) for value in (eta, phi, time)):
            return False
        # Rejects both genuinely negative times and invalid sentinels (-999).
        if time < 0.0:
            return False
        cluster_records.append((energy, eta, phi, time))

    e1, eta1, phi1, time1 = cluster_records[0]
    e2, eta2, phi2, time2 = cluster_records[1]
    if delta_r(eta1, phi1, eta2, phi2) > max_pair_delta_r:
        return False
    hist.Fill(e2 / e1, abs(time2 - time1))
    return True


def update_stats_for_selection(s, cluster_energies, n_cl, pion_energy, all_pfrh_energy,
                               clustered_pfrh_energy, clustered_rh_pct):
    neff = compute_neff(cluster_energies)
    reco_energy = sum(cluster_energies)
    sorted_energies = sorted(cluster_energies, reverse=True)
    top1_energy = sorted_energies[0] if len(sorted_energies) >= 1 else 0.0

    top1_over_pion = top1_energy / pion_energy if pion_energy > 0.0 else 0.0
    top1_over_all_pfrh = top1_energy / all_pfrh_energy if all_pfrh_energy > 0.0 else 0.0

    s["n_events"] += 1
    s["sum_ncl"] += n_cl
    s["sum_neff"] += neff
    s["sum_top1_over_pion"] += top1_over_pion
    s["sum_top1_over_all_pfrh"] += top1_over_all_pfrh
    s["sum_reco_energy"] += reco_energy
    if clustered_rh_pct is not None:
        s["sum_clustered_rh_pct"] += clustered_rh_pct
        s["n_valid_clustered_rh_pct"] += 1

    if pion_energy <= 0.0:
        s["n_missing_pion"] += 1
    if all_pfrh_energy <= 0.0:
        s["n_zero_all_pfrh"] += 1
    if clustered_pfrh_energy <= 0.0:
        s["n_zero_clustered_pfrh"] += 1

    if n_cl == 0:
        s["n_zero"] += 1
    if n_cl == 1:
        s["n_eq1"] += 1
    if n_cl >= 2:
        s["n_ge2"] += 1

    if n_cl == 1 and all_pfrh_energy > 0.0 and (top1_energy / all_pfrh_energy > 0.8):
        s["n_eq1_top1_all_pfrh_gt80"] += 1

    if n_cl >= 2 and clustered_pfrh_energy > 0.0 and (top1_energy / clustered_pfrh_energy > 0.8):
        s["n_ge2_top1_clustered_pfrh_gt80"] += 1

    if n_cl >= 2 and all_pfrh_energy > 0.0 and (top1_energy / all_pfrh_energy > 0.8):
        s["n_ge2_top1_all_pfrh_gt80"] += 1

    return n_cl, neff, reco_energy, clustered_rh_pct


def book_scatter_hists(args, pf_names, eta_values, sample_type, time_kind):
    """Book the |Delta t| versus E2/E1 densities for one time definition.

    Seed times are whole nanosecond values, so the seed pages get one bin per
    integer with the integer at the bin center and an explicit integer label,
    rather than the finer continuous binning used for PF cluster times.
    """
    if time_kind == "seedtime":
        seed_time_max = int(round(args.seed_time_max))
        time_bins = seed_time_max + 1
        time_lo, time_hi = -0.5, seed_time_max + 0.5
    else:
        time_bins = args.scatter_time_bins
        time_lo, time_hi = 0.0, args.scatter_time_max

    hists = {"nocut": {}, "ecut": {}}
    for mode in hists:
        for pf_name in pf_names:
            for eta in eta_values:
                tag = sanitize_name(f"{sample_type}_{time_kind}_{mode}_{pf_name}_eta{eta}")
                hist = ROOT.TH2D(
                    f"h2_cluster_pair_density_{tag}",
                    "",
                    args.scatter_x_bins,
                    0.0,
                    1.0,
                    time_bins,
                    time_lo,
                    time_hi,
                )
                if time_kind == "seedtime":
                    for ibin in range(1, time_bins + 1):
                        hist.GetYaxis().SetBinLabel(ibin, str(ibin - 1))
                hists[mode][(pf_name, key_float(eta))] = hist
    return hists


def process_one_ntuple(path, pf_name, eta, energy, args, hists, scatter_hists_by_kind,
                       stats_by_mode):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"WARNING: could not open {path}")
        return None

    tree = f.Get(args.tree_name)
    if not tree:
        print(f"WARNING: could not find tree {args.tree_name} in {path}")
        f.Close()
        return None

    required_branches = [
        args.energy_branch,
        args.cluster_eta_branch,
        args.cluster_phi_branch,
        args.cluster_time_branch,
        args.all_pfrh_energy_branch,
        args.all_pfrh_cluster_index_branch,
        args.clustered_pfrh_energy_branch,
        args.clustered_pfrh_cluster_index_branch,
        args.gen_energy_branch,
        args.gen_pdgid_branch,
    ]
    missing = [name for name in required_branches if not branch_exists(tree, name)]
    if missing:
        print(f"WARNING: required branches missing in {path}: " + ", ".join(missing) + ". Skipping this file.")
        f.Close()
        return None

    # The PFRecHit time branch is optional: without it only the seed-time
    # density pages lose this file, so the rest of the document still builds.
    has_seed_time = branch_exists(tree, args.seed_time_branch)
    if not has_seed_time:
        print(
            f"WARNING: {args.seed_time_branch} missing in {os.path.basename(path)}; "
            "this file contributes no seed-time cluster pairs."
        )

    eta_key = key_float(eta)
    energy_key = key_float(energy)
    key = (pf_name, eta_key, energy_key)
    for mode in stats_by_mode:
        stats_by_mode[mode].setdefault(key, empty_stats())

    h_ncl = hists["ncl"][key]
    h_neff = hists["neff"][key]
    h_clustered_rh_pct = hists["clustered_rh_pct"][key]
    h_reco_energy = hists["reco_energy"][key]
    h_ncl_ecut = hists["ncl_ecut"][key]
    h_neff_ecut = hists["neff_ecut"][key]
    h_clustered_rh_pct_ecut = hists["clustered_rh_pct_ecut"][key]
    h_reco_energy_ecut = hists["reco_energy_ecut"][key]
    scatter_slots = {
        (time_kind, mode): scatter_hists_by_kind[time_kind][mode][(pf_name, eta_key)]
        for time_kind in TIME_KINDS
        for mode in ("nocut", "ecut")
    }
    energy_threshold = cluster_energy_threshold(energy)

    n_entries = tree.GetEntries()
    print(f"Processing {pf_name:12s} eta={eta:<5} E={energy:<6}  {n_entries:7d} events  {os.path.basename(path)}")

    for event in tree:
        raw_cluster_energies = vector_to_list(getattr(event, args.energy_branch))
        raw_cluster_etas = vector_to_list(getattr(event, args.cluster_eta_branch))
        raw_cluster_phis = vector_to_list(getattr(event, args.cluster_phi_branch))
        raw_cluster_times = vector_to_list(getattr(event, args.cluster_time_branch))
        all_pfrh_energies = clean_positive_energies(vector_to_list(getattr(event, args.all_pfrh_energy_branch)))
        all_pfrh_cluster_indices = vector_to_list(getattr(event, args.all_pfrh_cluster_index_branch))
        clustered_pfrh_energy_fractions = vector_to_list(getattr(event, args.clustered_pfrh_energy_branch))
        clustered_pfrh_cluster_indices = vector_to_list(getattr(event, args.clustered_pfrh_cluster_index_branch))

        pion_energy = generator_single_pion_energy(event, args)
        all_pfrh_energy = sum(all_pfrh_energies)

        cluster_time_lookup = indexed_time_lookup(raw_cluster_times)
        if has_seed_time:
            seed_time_lookup = mapped_time_lookup(
                cluster_seed_times(
                    clustered_pfrh_energy_fractions,
                    clustered_pfrh_cluster_indices,
                    vector_to_list(getattr(event, args.seed_time_branch)),
                )
            )
        else:
            seed_time_lookup = None

        # No-cut selection
        cluster_energies_nocut, selected_indices_nocut, n_cl_nocut = selected_cluster_quantities(raw_cluster_energies, threshold=None)
        clustered_pfrh_energy_nocut = sum_clustered_pfrh_for_selected_clusters(
            clustered_pfrh_energy_fractions,
            clustered_pfrh_cluster_indices,
            selected_indices_nocut,
        )
        clustered_rh_pct_nocut = clustered_rechit_percentage(
            all_pfrh_cluster_indices,
            selected_indices_nocut,
        )
        n_cl_val, neff_val, reco_energy_val, clustered_rh_pct_val = update_stats_for_selection(
            stats_by_mode["nocut"][key],
            cluster_energies_nocut,
            n_cl_nocut,
            pion_energy,
            all_pfrh_energy,
            clustered_pfrh_energy_nocut,
            clustered_rh_pct_nocut,
        )
        h_ncl.Fill(n_cl_val)
        h_neff.Fill(neff_val)
        h_reco_energy.Fill(reco_energy_val)
        if clustered_rh_pct_val is not None:
            h_clustered_rh_pct.Fill(clustered_rh_pct_val)
        fill_cluster_pair_density(
            scatter_slots[("clustertime", "nocut")],
            raw_cluster_energies,
            raw_cluster_etas,
            raw_cluster_phis,
            cluster_time_lookup,
            selected_indices_nocut,
        )
        if seed_time_lookup is not None:
            fill_cluster_pair_density(
                scatter_slots[("seedtime", "nocut")],
                raw_cluster_energies,
                raw_cluster_etas,
                raw_cluster_phis,
                seed_time_lookup,
                selected_indices_nocut,
            )

        # Energy-dependent selection: E_cluster >= sqrt(0.7 E_pi) GeV.
        cluster_energies_cut, selected_indices_cut, n_cl_cut = selected_cluster_quantities(
            raw_cluster_energies,
            threshold=energy_threshold,
        )
        clustered_pfrh_energy_cut = sum_clustered_pfrh_for_selected_clusters(
            clustered_pfrh_energy_fractions,
            clustered_pfrh_cluster_indices,
            selected_indices_cut,
        )
        clustered_rh_pct_cut = clustered_rechit_percentage(
            all_pfrh_cluster_indices,
            selected_indices_cut,
        )
        n_cl_cut_val, neff_cut_val, reco_energy_cut_val, clustered_rh_pct_cut_val = update_stats_for_selection(
            stats_by_mode["ecut"][key],
            cluster_energies_cut,
            n_cl_cut,
            pion_energy,
            all_pfrh_energy,
            clustered_pfrh_energy_cut,
            clustered_rh_pct_cut,
        )
        h_ncl_ecut.Fill(n_cl_cut_val)
        h_neff_ecut.Fill(neff_cut_val)
        h_reco_energy_ecut.Fill(reco_energy_cut_val)
        if clustered_rh_pct_cut_val is not None:
            h_clustered_rh_pct_ecut.Fill(clustered_rh_pct_cut_val)
        fill_cluster_pair_density(
            scatter_slots[("clustertime", "ecut")],
            raw_cluster_energies,
            raw_cluster_etas,
            raw_cluster_phis,
            cluster_time_lookup,
            selected_indices_cut,
        )
        if seed_time_lookup is not None:
            fill_cluster_pair_density(
                scatter_slots[("seedtime", "ecut")],
                raw_cluster_energies,
                raw_cluster_etas,
                raw_cluster_phis,
                seed_time_lookup,
                selected_indices_cut,
            )

    f.Close()
    return True


# -----------------------------------------------------------------------------
# Build heatmaps from stats
# -----------------------------------------------------------------------------

def book_heatmaps(pf_names, eta_values, energy_values, mode, sample_type):
    heatmaps = {}
    mode_tag = sanitize_name(mode)
    for pf_name in pf_names:
        pftag = sanitize_name(pf_name)
        suffix = f"{sanitize_name(sample_type)}_{pftag}_{mode_tag}"
        heatmaps[pf_name] = {
            "mean_ncl": make_scan_hist(f"h2_mean_nclusters_{suffix}", f"Average Number of HCAL Clusters per Event ({pf_name}, {mode})", eta_values, energy_values),
            "mean_neff": make_scan_hist(f"h2_mean_neff_{suffix}", f"Average Effective Number of Energy-Carrying HCAL Clusters ({pf_name}, {mode})", eta_values, energy_values),
            "mean_top1_over_pion": make_scan_hist(f"h2_mean_top1_over_pion_{suffix}", f"Average E1/Epi ({pf_name}, {mode})", eta_values, energy_values),
            "mean_top1_over_all_pfrh": make_scan_hist(f"h2_mean_top1_over_all_pfrh_{suffix}", f"Average E1/EallPFRecHits ({pf_name}, {mode})", eta_values, energy_values),
            "mean_clustered_rh_pct": make_scan_hist(f"h2_mean_clustered_rh_pct_{suffix}", f"Fraction of Clustered HCAL PFRecHits ({pf_name}, {mode})", eta_values, energy_values),
            "frac_eq1": make_scan_hist(f"h2_frac_eq1_{suffix}", f"Percentage of Events with Exactly 1 HCAL Cluster ({pf_name}, {mode})", eta_values, energy_values),
            "frac_zero": make_scan_hist(f"h2_frac_zero_{suffix}", f"Percentage of Events with 0 HCAL Clusters ({pf_name}, {mode})", eta_values, energy_values),
            "frac_ge2": make_scan_hist(f"h2_frac_ge2_{suffix}", f"Percentage of Events with >=2 HCAL Clusters ({pf_name}, {mode})", eta_values, energy_values),
            "frac_eq1_top1_all_pfrh_gt80": make_scan_hist(f"h2_frac_eq1_top1_all_pfrh_gt80_{suffix}", f"Percentage with 1 Cluster and E1/EallPFRecHits>0.8 ({pf_name}, {mode})", eta_values, energy_values),
            "frac_ge2_top1_clustered_pfrh_gt80": make_scan_hist(f"h2_frac_ge2_top1_clustered_pfrh_gt80_{suffix}", f"Percentage with >=2 Clusters and E1/EclusteredPFRecHits>0.8 ({pf_name}, {mode})", eta_values, energy_values),
            "frac_ge2_top1_all_pfrh_gt80": make_scan_hist(f"h2_frac_ge2_top1_all_pfrh_gt80_{suffix}", f"Percentage with >=2 Clusters and E1/EallPFRecHits>0.8 ({pf_name}, {mode})", eta_values, energy_values),
        }
    return heatmaps


def fill_heatmaps_from_stats(heatmaps, stats, pf_names, eta_values, energy_values):
    for pf_name in pf_names:
        for ix, eta in enumerate(eta_values, start=1):
            for iy, energy in enumerate(energy_values, start=1):
                key = (pf_name, key_float(eta), key_float(energy))
                s = stats.get(key)
                if s is None or s["n_events"] <= 0:
                    continue
                n = float(s["n_events"])
                n_valid_rh = float(s["n_valid_clustered_rh_pct"])
                values = {
                    "mean_ncl": s["sum_ncl"] / n,
                    "mean_neff": s["sum_neff"] / n,
                    "mean_top1_over_pion": s["sum_top1_over_pion"] / n,
                    "mean_top1_over_all_pfrh": s["sum_top1_over_all_pfrh"] / n,
                    "mean_clustered_rh_pct": (
                        s["sum_clustered_rh_pct"] / n_valid_rh
                        if n_valid_rh > 0.0 else 0.0
                    ),
                    "frac_eq1": 100.0 * s["n_eq1"] / n,
                    "frac_zero": 100.0 * s["n_zero"] / n,
                    "frac_ge2": 100.0 * s["n_ge2"] / n,
                    "frac_eq1_top1_all_pfrh_gt80": 100.0 * s["n_eq1_top1_all_pfrh_gt80"] / n,
                    "frac_ge2_top1_clustered_pfrh_gt80": 100.0 * s["n_ge2_top1_clustered_pfrh_gt80"] / n,
                    "frac_ge2_top1_all_pfrh_gt80": 100.0 * s["n_ge2_top1_all_pfrh_gt80"] / n,
                }
                for metric_key, value in values.items():
                    heatmaps[pf_name][metric_key].SetBinContent(ix, iy, value)


def make_scatter_density_pages(scatter_hists, pf_names, pdf_path, eta_values,
                               sample_label, sample_type, time_kind):
    """Draw one page per selection and PF algorithm, with every eta together."""
    if len(pf_names) < 2:
        raise RuntimeError("Cluster-pair pages require two PF algorithms.")
    algorithms = (
        ["standardPF", "seedTimingPF"]
        if "standardPF" in pf_names and "seedTimingPF" in pf_names
        else pf_names[:2]
    )
    is_seed = time_kind == "seedtime"
    delta_symbol = "#Deltat_{seed}" if is_seed else "#Deltat"
    y_axis_title = (
        "#Deltat_{seed} = |t^{seed}_{2} - t^{seed}_{1}| [ns]"
        if is_seed
        else "#Deltat = |t_{2} - t_{1}| [ns]"
    )
    timing_requirement = (
        "Two leading clusters: #DeltaR_{12} #leq 0.4 and seed times "
        "t^{seed}_{1}, t^{seed}_{2} #geq 0 ns "
        "(seed = highest-energy PFRecHit in the cluster; invalid times dropped)"
        if is_seed
        else "Two leading clusters: #DeltaR_{12} #leq 0.4 and t_{1}, t_{2} #geq 0 ns "
        "(valid cluster-time information required)"
    )
    ncols = min(3, max(1, len(eta_values)))
    nrows = max(1, int(math.ceil(len(eta_values) / float(ncols))))

    # Required page order: no-cut standardPF, no-cut seedTimingPF,
    # sqrt(0.7 E_pi)-cut standardPF, then seedTimingPF.
    for mode in ("nocut", "ecut"):
        for pf_name in algorithms:
            selection = (
                "No HCAL-cluster energy cut"
                if mode == "nocut"
                else "E_{cluster} #geq #sqrt{0.7 E_{#pi}} GeV (sample-by-sample)"
            )
            tag = sanitize_name(f"{sample_type}_scatter_{time_kind}_{mode}_{pf_name}")
            canvas_height = 650 + 430 * nrows
            c = ROOT.TCanvas(f"c_{tag}", "", 2100, canvas_height)
            keep = []

            header = ROOT.TPad(f"header_{tag}", "", 0.0, 0.88, 1.0, 1.0)
            header.SetFillStyle(0)
            header.Draw()
            header.cd()
            line1 = ROOT.TLatex(
                0.5,
                0.76,
                f"{sample_label}: {delta_symbol} vs. Energy Fraction - {pf_name}",
            )
            line1.SetNDC()
            line1.SetTextAlign(22)
            line1.SetTextFont(62)
            line1.SetTextSize(0.23)
            line1.Draw()
            line2 = ROOT.TLatex(
                0.5,
                0.43,
                f"{selection}; all E_{{#pi}} samples at this #eta combined",
            )
            line2.SetNDC()
            line2.SetTextAlign(22)
            line2.SetTextFont(42)
            line2.SetTextSize(0.15)
            line2.Draw()
            line3 = ROOT.TLatex(0.5, 0.14, timing_requirement)
            line3.SetNDC()
            line3.SetTextAlign(22)
            line3.SetTextFont(42)
            line3.SetTextSize(0.15)
            line3.Draw()
            keep.extend([header, line1, line2, line3])

            c.cd()
            body = ROOT.TPad(f"body_{tag}", "", 0.0, 0.0, 1.0, 0.88)
            body.SetFillStyle(0)
            body.Draw()
            body.cd()
            body.Divide(ncols, nrows, 0.002, 0.002)
            keep.append(body)

            present = [
                scatter_hists[mode].get((pf_name, key_float(eta)))
                for eta in eta_values
            ]
            global_max = max([hist.GetMaximum() for hist in present if hist is not None] + [1.0])
            # Small integer counts are clearer on a linear 0, 1, 2, ... scale.
            # Retain a logarithmic scale only when the page spans at least one
            # full decade.  Each page uses one common scale across all etas.
            use_logz = global_max >= 10.0
            set_density_palette(show_zero_as_white=not use_logz)

            for index, eta in enumerate(eta_values):
                pad = body.cd(index + 1)
                pad.SetLeftMargin(0.13)
                pad.SetRightMargin(0.17)
                pad.SetBottomMargin(0.13)
                pad.SetTopMargin(0.12)
                pad.SetTicks(1, 1)
                hist = scatter_hists[mode].get((pf_name, key_float(eta)))
                if hist is None or hist.GetEntries() <= 0:
                    missing = ROOT.TLatex(0.5, 0.53, f"No selected pairs: #eta={eta}")
                    missing.SetNDC()
                    missing.SetTextAlign(22)
                    missing.SetTextColor(ROOT.kGray + 2)
                    missing.Draw()
                    keep.append(missing)
                    continue
                pad.SetLogz(use_logz)
                hist.SetTitle(f"#eta = {eta}")
                hist.GetXaxis().SetTitle("Energy fraction E_{2} / E_{1}")
                hist.GetYaxis().SetTitle(y_axis_title)
                hist.GetZaxis().SetTitle("Number of events per bin")
                for axis in (hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()):
                    axis.SetTitleSize(0.045)
                    axis.SetLabelSize(0.036)
                hist.GetXaxis().SetTitleOffset(1.08)
                hist.GetYaxis().SetTitleOffset(1.22)
                hist.GetZaxis().SetTitleOffset(1.30)
                if use_logz:
                    # Zero-count bins lie below the displayed range and remain
                    # white.  Use ordinary count labels rather than exponents.
                    hist.SetMinimum(1.0)
                    hist.SetMaximum(global_max)
                    hist.GetZaxis().SetMoreLogLabels(True)
                    hist.GetZaxis().SetNdivisions(510)
                else:
                    # Include zero explicitly: the palette maps it to white,
                    # while positive integer counts use sequential colors.
                    hist.SetMinimum(0.0)
                    hist.SetMaximum(global_max)
                    hist.GetZaxis().SetMoreLogLabels(False)
                    hist.GetZaxis().SetNdivisions(
                        max(1, int(round(global_max))), False
                    )
                hist.GetZaxis().SetNoExponent(True)
                hist.Draw("COLZ")
                pad.Update()

                # ROOT creates the color-palette axis only after the pad is
                # painted, so apply the same formatting directly to it too.
                palette = hist.GetListOfFunctions().FindObject("palette")
                if palette:
                    palette_axis = palette.GetAxis()
                    palette_axis.SetNoExponent(True)
                    palette_axis.SetMoreLogLabels(use_logz)
                    if use_logz:
                        palette_axis.SetNdivisions(510)
                    else:
                        # TGaxis::SetNdivisions accepts only one argument in
                        # PyROOT.  A negative value is the one-argument
                        # equivalent of TAxis::SetNdivisions(ndiv, False): it
                        # keeps the requested integer divisions unoptimized.
                        palette_axis.SetNdivisions(
                            -max(1, int(round(global_max)))
                        )
                    palette_axis.SetLabelSize(0.036)
                    palette_axis.SetTitleSize(0.045)
                    palette_axis.SetTitleOffset(1.30)
                    pad.Modified()
                    pad.Update()

            c.Print(pdf_path)
            print(
                f"Added {time_kind} cluster-pair density page: "
                f"{sample_label}, {pf_name}, all etas, {selection}"
            )


# -----------------------------------------------------------------------------
# Output helpers
# -----------------------------------------------------------------------------

def print_summary(stats_by_mode, pf_names, eta_values, energy_values, sample_label):
    for mode, stats in stats_by_mode.items():
        cut_label = (
            "no cluster-energy cut"
            if mode == "nocut"
            else "Ecluster >= sqrt(0.7 Epi) GeV"
        )
        print()
        print(f"Summary from ntuples: {sample_label}; {cut_label}")
        print("-" * 235)
        print(
            f"{'PF':<14} {'eta':>6} {'E':>7} {'Events':>8} "
            f"{'<Ncl>':>8} {'<N_eff>':>9} "
            f"{'<E1/Epi>':>12} {'<E1/E_allRH>':>15} "
            f"{'%0':>7} {'%==1':>8} {'%>=2':>8} "
            f"{'%1cl & E1/E_allRH>0.8':>23} "
            f"{'%>=2 & E1/E_inClRH>0.8':>25} "
            f"{'%>=2 & E1/E_allRH>0.8':>24}"
        )
        print("-" * 235)
        for pf_name in pf_names:
            for eta in eta_values:
                for energy in energy_values:
                    key = (pf_name, key_float(eta), key_float(energy))
                    s = stats.get(key)
                    if s is None or s["n_events"] <= 0:
                        continue
                    n = float(s["n_events"])
                    print(
                        f"{pf_name:<14} {eta:6.2f} {energy:7.1f} {int(n):8d} "
                        f"{s['sum_ncl']/n:8.3f} {s['sum_neff']/n:9.3f} "
                        f"{s['sum_top1_over_pion']/n:12.3f} "
                        f"{s['sum_top1_over_all_pfrh']/n:15.3f} "
                        f"{100.0*s['n_zero']/n:7.2f} "
                        f"{100.0*s['n_eq1']/n:8.2f} "
                        f"{100.0*s['n_ge2']/n:8.2f} "
                        f"{100.0*s['n_eq1_top1_all_pfrh_gt80']/n:23.2f} "
                        f"{100.0*s['n_ge2_top1_clustered_pfrh_gt80']/n:25.2f} "
                        f"{100.0*s['n_ge2_top1_all_pfrh_gt80']/n:24.2f}"
                    )
        print("-" * 235)


def write_output_root(out_path, results_by_sample_type, pf_names):
    out = ROOT.TFile(out_path, "RECREATE")

    for sample_type, result in results_by_sample_type.items():
        sample_dir = out.mkdir(sanitize_name(sample_type))
        sample_dir.cd()

        for mode, heatmaps in result["heatmaps_by_mode"].items():
            mode_dir = sample_dir.mkdir(f"heatmaps_{mode}")
            mode_dir.cd()
            for pf_name in pf_names:
                if pf_name not in heatmaps:
                    continue
                pf_dir = mode_dir.mkdir(sanitize_name(pf_name))
                pf_dir.cd()
                for hist in heatmaps[pf_name].values():
                    hist.Write()
                mode_dir.cd()
            ratios = result.get("ratio_heatmaps_by_mode", {}).get(mode, {})
            if ratios:
                ratio_dir = mode_dir.mkdir("ratio")
                ratio_dir.cd()
                for entry in ratios.values():
                    entry[0].Write()
                mode_dir.cd()
            sample_dir.cd()

        distributions_dir = sample_dir.mkdir("distributions")
        distributions_dir.cd()
        for group_name, group in result["dist_hists"].items():
            group_dir = distributions_dir.mkdir(group_name)
            group_dir.cd()
            for hist in group.values():
                hist.Write()
            distributions_dir.cd()

        for time_kind in TIME_KINDS:
            scatter_dir = sample_dir.mkdir(f"cluster_pair_density_{time_kind}")
            scatter_dir.cd()
            for mode, group in result["scatter_hists_by_kind"][time_kind].items():
                mode_dir = scatter_dir.mkdir(mode)
                mode_dir.cd()
                for hist in group.values():
                    hist.Write()
                scatter_dir.cd()
            sample_dir.cd()

        out.cd()

    out.Close()


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    args = parse_args()

    pf_names = []
    for pf in args.pf_names:
        if pf not in pf_names:
            pf_names.append(pf)
    if len(pf_names) < 2:
        raise ValueError("Ratio panels require at least two PF algorithms.")

    sample_types = []
    for value in args.sample_types:
        canonical = canonical_sample_type(value)
        if canonical not in sample_types:
            sample_types.append(canonical)

    requested_eta_values = parse_float_list(args.eta_values, None)
    requested_energy_values = parse_float_list(args.energy_values, None)
    requested_eta_dist = parse_float_list(args.eta_dist_values, None)
    requested_energy_dist = parse_float_list(args.energy_dist_values, None)

    metadata = collect_file_metadata(
        args.input_dir,
        args.prefix,
        pf_names,
        sample_types,
        debug=args.debug_files,
    )
    if len(metadata) == 0:
        raise RuntimeError(
            "No matching ntuple files found. Expected names such as "
            "pfObjectsNtuple_standardPF_SinglePiCloseByE20_eta0p1.root "
            "or pfObjectsNtuple_seedTimingPF_SinglePiCloseByE120_eta1.root."
        )

    print("PF variants:", ", ".join(pf_names))
    print("Sample families:", ", ".join(sample_types))
    print(
        "Generator-pion denominator: event-by-event from "
        f"{args.gen_energy_branch}/{args.gen_pdgid_branch}"
    )
    print(
        "Cluster-pair density selection: DeltaR_12 between the two leading "
        "retained clusters <= 0.4, with finite cluster times t1,t2 >= 0 ns"
    )
    print(
        "Seed-time cluster-pair density: same selection, using the time of the "
        f"highest-energy PFRecHit in each cluster from {args.seed_time_branch}"
    )
    print(
        "In-cluster PFRecHit denominator for threshold plots: event-by-event from "
        f"{args.clustered_pfrh_energy_branch}, grouped by "
        f"{args.clustered_pfrh_cluster_index_branch}"
    )
    print(
        "All-HCAL-PFRecHit denominator: event-by-event from "
        f"{args.all_pfrh_energy_branch}"
    )
    print(
        "Cluster-energy cut: E_cluster >= sqrt(0.7 E_pi) GeV, "
        "using E_pi from each sample filename"
    )
    print()

    results_by_sample_type = {}

    for sample_type in sample_types:
        sample_metadata = [
            item for item in metadata if item["sample_type"] == sample_type
        ]
        if len(sample_metadata) == 0:
            print(f"WARNING: no {sample_type} files were found; skipping this family.")
            continue

        eta_values, energy_values, available_etas, available_energies = (
            select_scan_values(
                sample_metadata,
                requested_eta_values,
                requested_energy_values,
            )
        )
        if len(eta_values) == 0 or len(energy_values) == 0:
            print(
                f"WARNING: no eta/energy scan points remain for {sample_type}; "
                "skipping this family."
            )
            continue

        eta_dist_values = (
            [key_float(x) for x in requested_eta_dist]
            if requested_eta_dist is not None
            else eta_values[:3]
        )
        energy_dist_values = (
            [key_float(x) for x in requested_energy_dist]
            if requested_energy_dist is not None
            else energy_values[:3]
        )

        print(f"--- {sample_type} samples ---")
        print("Available eta values:", available_etas)
        print("Available energy values:", available_energies)
        print("Using eta scan values:", eta_values)
        print("Using energy scan values:", energy_values)
        print("Overlay eta values:", eta_dist_values)
        print("Overlay energy values:", energy_dist_values)

        file_lookup = build_file_lookup(
            sample_metadata,
            eta_values,
            energy_values,
            debug=args.debug_files,
        )
        if len(file_lookup) == 0:
            print(f"WARNING: no matching {sample_type} files after filtering; skipping.")
            continue

        for pf_name in pf_names:
            for eta in eta_values:
                for energy in energy_values:
                    key = (pf_name, key_float(eta), key_float(energy))
                    if key not in file_lookup:
                        print(
                            f"Missing ntuple: type={sample_type}, pf={pf_name}, "
                            f"eta={eta}, E={energy}"
                        )
        print()

        dist_hists = book_distribution_hists(
            args,
            pf_names,
            eta_values,
            energy_values,
            sample_type,
        )
        scatter_hists_by_kind = {
            time_kind: book_scatter_hists(
                args, pf_names, eta_values, sample_type, time_kind
            )
            for time_kind in TIME_KINDS
        }
        stats_by_mode = {"nocut": {}, "ecut": {}}

        for pf_name in pf_names:
            for eta in eta_values:
                for energy in energy_values:
                    key = (pf_name, key_float(eta), key_float(energy))
                    path = file_lookup.get(key)
                    if path is None:
                        continue
                    process_one_ntuple(
                        path,
                        pf_name,
                        eta,
                        energy,
                        args,
                        dist_hists,
                        scatter_hists_by_kind,
                        stats_by_mode,
                    )

        sample_label = "CloseByParticleGunProducer (HCAL-surface launch)"
        print_summary(
            stats_by_mode,
            pf_names,
            eta_values,
            energy_values,
            sample_label,
        )

        heatmaps_by_mode = {}
        for mode in stats_by_mode:
            heatmaps_by_mode[mode] = book_heatmaps(
                pf_names,
                eta_values,
                energy_values,
                mode,
                sample_type,
            )
            fill_heatmaps_from_stats(
                heatmaps_by_mode[mode],
                stats_by_mode[mode],
                pf_names,
                eta_values,
                energy_values,
            )

        results_by_sample_type[sample_type] = {
            "sample_label": sample_label,
            "eta_values": eta_values,
            "energy_values": energy_values,
            "eta_dist_values": eta_dist_values,
            "energy_dist_values": energy_dist_values,
            "dist_hists": dist_hists,
            "scatter_hists_by_kind": scatter_hists_by_kind,
            "stats_by_mode": stats_by_mode,
            "heatmaps_by_mode": heatmaps_by_mode,
        }

    if len(results_by_sample_type) == 0:
        raise RuntimeError("No requested single-pion sample family could be processed.")

    # Every ratio panel gets its own range, computed from its own values.
    build_ratio_heatmaps(results_by_sample_type, pf_names)
    install_heatmap_palettes()
    os.makedirs(args.output_dir, exist_ok=True)

    pf_tag = "_".join(sanitize_name(pf) for pf in pf_names)
    sample_tag = "_".join(sanitize_name(x) for x in results_by_sample_type)
    if args.output_pdf is None:
        pdf_name = (
            f"hcal_cluster_onepi_scan_{sample_tag}_{pf_tag}_eta_energy.pdf"
        )
    else:
        pdf_name = args.output_pdf
    if args.output_root is None:
        root_name = (
            f"hcal_cluster_onepi_scan_{sample_tag}_{pf_tag}_eta_energy.root"
        )
    else:
        root_name = args.output_root

    pdf_path = os.path.join(args.output_dir, pdf_name)
    root_path = os.path.join(args.output_dir, root_name)

    c_pdf = ROOT.TCanvas("c_pdf_open", "", 2200, 900)
    c_pdf.SetCanvasSize(2200, 900)
    c_pdf.Print(pdf_path + "[")

    metric_pages = [
        ("mean_ncl", "Average Number of HCAL Clusters per Event", "Average # clusters/event", HEATMAP_VALUE_FORMAT, (0.0, 3.0), (0.0, 3.0)),
        ("mean_neff", "Average Effective Number of Energy-Carrying HCAL Clusters", "#LTN_{eff}#GT = (#Sigma E)^{2}/#Sigma E^{2}", HEATMAP_VALUE_FORMAT, (0.0, 3.0), (0.0, 3.0)),
        ("mean_top1_over_pion", "Average Leading HCAL Cluster Energy / Generator-Pion Energy", "#LTE_{1}/E_{#pi}#GT", HEATMAP_VALUE_FORMAT, (0.0, 1.0), (0.0, 1.0)),
        ("mean_top1_over_all_pfrh", "Average Leading HCAL Cluster Energy / All Event HCAL PFRecHit Energy", "#LTE_{1}/E_{all event PFRecHits}#GT", HEATMAP_VALUE_FORMAT, (0.0, 1.2), (0.0, 1.2)),
        ("mean_clustered_rh_pct", "Fraction of Clustered HCAL PFRecHits", "#LTclustered PFRecHits / all PFRecHits#GT [%]", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_eq1", "Percentage of Events with Exactly 1 HCAL Cluster", "% events with exactly 1 cluster", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_zero", "Percentage of Events with 0 HCAL Clusters", "% events with 0 clusters", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_ge2", "Percentage of Events with #geq 2 HCAL Clusters", "% events with #geq 2 clusters", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_eq1_top1_all_pfrh_gt80", "Percentage with Exactly 1 Cluster and E_{1}/E_{all PFRecHits}>0.8", "% all events passing both conditions", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_ge2_top1_clustered_pfrh_gt80", "Percentage with #geq 2 Clusters and E_{1}/E_{PFRecHits in clusters}>0.8", "% all events passing both conditions", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
        ("frac_ge2_top1_all_pfrh_gt80", "Percentage with #geq 2 Clusters and E_{1}/E_{all event PFRecHits}>0.8", "% all events passing both conditions", HEATMAP_PERCENT_FORMAT, (0.0, 100.0), (0.0, 100.0)),
    ]

    for sample_type, result in results_by_sample_type.items():
        sample_label = result["sample_label"]
        heatmaps_by_mode = result["heatmaps_by_mode"]
        ratio_heatmaps_by_mode = result["ratio_heatmaps_by_mode"]
        eta_dist_values = result["eta_dist_values"]
        energy_dist_values = result["energy_dist_values"]
        dist_hists = result["dist_hists"]
        scatter_hists_by_kind = result["scatter_hists_by_kind"]

        # For each quantity, put the no-cut page immediately before its
        # sqrt(0.7 E_pi)-cut page.
        for (
            metric_key,
            title,
            z_title,
            text_format,
            nocut_range,
            cut_range,
        ) in metric_pages:
            for mode in ["nocut", "ecut"]:
                if mode == "nocut":
                    cut_label = "No HCAL-cluster energy cut"
                    zmin, zmax = nocut_range
                else:
                    cut_label = (
                        "HCAL clusters with E #geq #sqrt{0.7 E_{#pi}} GeV "
                        "(sample-by-sample)"
                    )
                    zmin, zmax = cut_range
                make_paired_heatmap(
                    heatmaps_by_mode[mode],
                    ratio_heatmaps_by_mode[mode].get(metric_key),
                    pf_names,
                    metric_key,
                    title,
                    z_title,
                    pdf_path,
                    cut_label,
                    sample_label,
                    text_format=text_format,
                    zmin=zmin,
                    zmax=zmax,
                )

        reco_xmax = 1.5 * max(result["energy_values"])
        overlay_specs = [
            (
                "ncl",
                "N_{clusters}",
                f"{sample_label}: HCAL cluster multiplicity; no cluster-energy cut",
                f"{sample_type}_nclusters",
                1.0,
                (-0.5, int(round(args.ncl_max)) + 0.5),
            ),
            (
                "ncl_ecut",
                "N_{clusters}",
                f"{sample_label}: HCAL cluster multiplicity; E_{{cluster}} #geq #sqrt{{0.7 E_{{#pi}}}} GeV",
                f"{sample_type}_nclusters_ecut",
                1.0,
                (-0.5, int(round(args.ncl_max)) + 0.5),
            ),
            (
                "neff",
                "N_{eff}",
                f"{sample_label}: Effective energy-carrying clusters; no cluster-energy cut",
                f"{sample_type}_neff",
                1.0,
                (0.0, args.neff_max),
            ),
            (
                "neff_ecut",
                "N_{eff}",
                f"{sample_label}: Effective clusters; E_{{cluster}} #geq #sqrt{{0.7 E_{{#pi}}}} GeV",
                f"{sample_type}_neff_ecut",
                1.0,
                (0.0, args.neff_max),
            ),
            (
                "clustered_rh_pct",
                "clustered HCAL PFRecHits [%]",
                f"{sample_label}: Fraction of Clustered HCAL PFRecHits; no cut",
                f"{sample_type}_clustered_rh_pct",
                None,
                (0.0, 100.0),
            ),
            (
                "clustered_rh_pct_ecut",
                "clustered HCAL PFRecHits [%]",
                f"{sample_label}: Fraction of Clustered HCAL PFRecHits; E_{{cluster}} #geq #sqrt{{0.7 E_{{#pi}}}} GeV",
                f"{sample_type}_clustered_rh_pct_ecut",
                None,
                (0.0, 100.0),
            ),
            (
                "reco_energy",
                "#Sigma E_{cluster} [GeV]",
                f"{sample_label}: Reconstructed HCAL cluster energy; no cluster-energy cut",
                f"{sample_type}_reco_energy",
                None,
                (0.0, reco_xmax),
            ),
            (
                "reco_energy_ecut",
                "#Sigma E_{cluster} [GeV]",
                f"{sample_label}: Reconstructed HCAL energy with E_{{cluster}} #geq #sqrt{{0.7 E_{{#pi}}}} GeV",
                f"{sample_type}_reco_energy_ecut",
                None,
                (0.0, reco_xmax),
            ),
        ]

        for (
            group_name,
            xaxis_title,
            grid_title,
            tag,
            ideal_line,
            x_range,
        ) in overlay_specs:
            # Show the clustered-RecHit-fraction distributions as raw event
            # counts.  The other 1D overlays retain their normalized shapes.
            normalize_overlay = group_name not in {
                "clustered_rh_pct",
                "clustered_rh_pct_ecut",
            }
            make_overlay_distribution_grid(
                dist_hists[group_name],
                pf_names,
                pdf_path,
                eta_dist_values,
                energy_dist_values,
                xaxis_title=xaxis_title,
                grid_title=grid_title,
                tag=tag,
                normalize=normalize_overlay,
                ideal_line=ideal_line,
                x_range=x_range,
            )

        # Cluster-time density pages first, then the seed-time versions.
        for time_kind in TIME_KINDS:
            scatter_hists = scatter_hists_by_kind[time_kind]
            has_entries = any(
                hist.GetEntries() > 0
                for group in scatter_hists.values()
                for hist in group.values()
            )
            if not has_entries:
                print(
                    f"WARNING: no {time_kind} cluster pairs were selected for "
                    f"{sample_label}; skipping those density pages."
                )
                continue
            make_scatter_density_pages(
                scatter_hists,
                pf_names,
                pdf_path,
                result["eta_values"],
                sample_label,
                sample_type,
                time_kind,
            )

    c_pdf.Print(pdf_path + "]")
    write_output_root(root_path, results_by_sample_type, pf_names)

    print()
    print("Saved combined PDF:")
    print(f"  {pdf_path}")
    print("Saved ROOT file:")
    print(f"  {root_path}")


if __name__ == "__main__":
    main()