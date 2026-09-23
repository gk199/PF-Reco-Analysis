#!/usr/bin/env python3
"""Make single-pion PF-cluster and generator diagnostics versus energy and eta.

The input files are expected to follow the naming convention

  pfObjectsNtuple_<PF>_SinglePiCloseByE<energy>_eta<eta>.root

with decimal points encoded as ``p`` (for example, ``eta0p01``).  The script
writes the requested fourteen pages to one multipage PDF.  Heatmaps use all
configured energy/eta points, while the distribution grids use E = 10, 60,
110 GeV and eta = 0.01, 0.2, 0.4 by default.

Definitions
-----------
* PF-cluster multiplicity is the length of the corresponding ``*_energy``
  vector in each event.
* Event-level N_eff is (sum_i E_i)^2 / sum_i E_i^2 using positive, finite
  PF-cluster energies.  Events without clusters have N_eff = 0 and remain in
  the sample average.
* Generator pions satisfy abs(gen_pdgId) == 211 and gen_status == 1.
* Generator-pion R is the transverse production radius sqrt(vx^2 + vy^2).
* The leading HCAL PF cluster is the highest-energy cluster in the event
  (positive, finite energy); events without one are skipped.  Its time is
  only used when it is not -999.
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
ROOT.gStyle.SetTitleFont(42, "XYZ")
ROOT.gStyle.SetLabelFont(42, "XYZ")


DEFAULT_ENERGIES = [20.0, 40.0, 60.0, 80.0, 100.0]
DEFAULT_ETAS = [0.01, 0.1, 0.2, 0.4, 0.6]
DEFAULT_GRID_ENERGIES = [20.0, 60.0, 100.0]
DEFAULT_GRID_ETAS = [0.01, 0.2, 0.4]

# Blue -> white -> orange, matching the supplied reference plot.
HEATMAP_PALETTE = (
    [0.00, 0.50, 1.00],
    [0.20, 1.00, 0.95],
    [0.55, 1.00, 0.55],
    [0.80, 1.00, 0.20],
)


def parse_float_list(value, default):
    if value is None:
        return list(default)
    if isinstance(value, (list, tuple)):
        return [float(x) for x in value]
    return [float(x.strip()) for x in str(value).split(",") if x.strip()]


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot single-pion PF-cluster and generator diagnostics."
    )
    parser.add_argument(
        "--input-dir",
        "--inputdir",
        dest="input_dir",
        default="/eos/user/c/chtong/Public/Rereco/SinglePion_SmallEta_modified",
    )
    parser.add_argument(
        "--output-dir",
        "--outputdir",
        dest="output_dir",
        default="/eos/user/c/chtong/Public/Rereco/SinglePion_SmallEta_modified",
    )
    parser.add_argument("--prefix", default="pfObjectsNtuple_")
    parser.add_argument("--pf-name", default="standardPF")
    parser.add_argument("--tree-name", default="pfObjectsNtupler/pfTree")
    parser.add_argument("--energies", default=",".join(map(str, DEFAULT_ENERGIES)))
    parser.add_argument("--etas", default=",".join(map(str, DEFAULT_ETAS)))
    parser.add_argument(
        "--grid-energies", default=",".join(map(str, DEFAULT_GRID_ENERGIES))
    )
    parser.add_argument("--grid-etas", default=",".join(map(str, DEFAULT_GRID_ETAS)))
    parser.add_argument(
        "--output-pdf", default="single_pion_small_eta_diagnostics_standardPF.pdf"
    )
    parser.add_argument("--cluster-eta-bins", type=int, default=64)
    parser.add_argument("--cluster-eta-min", type=float, default=-3.2)
    parser.add_argument("--cluster-eta-max", type=float, default=3.2)
    parser.add_argument("--cluster-energy-bins", type=int, default=50)
    parser.add_argument("--cluster-time-bins", type=int, default=50)
    parser.add_argument("--gen-energy-bins", type=int, default=50)
    parser.add_argument("--gen-eta-bins", type=int, default=50)
    parser.add_argument("--gen-radius-bins", type=int, default=50)
    parser.add_argument("--debug-files", action="store_true")
    return parser.parse_args()


def key_float(value):
    return round(float(value), 8)


def format_value(value):
    value = float(value)
    return str(int(value)) if value.is_integer() else f"{value:g}"


def number_tag(value):
    return format_value(value).replace("-", "m").replace(".", "p")


def decode_number_tag(tag):
    return float(str(tag).replace("m", "-").replace("p", "."))


def sanitize_name(text):
    return re.sub(r"[^0-9a-zA-Z_]+", "_", str(text)).strip("_")


def vector_to_list(vector):
    if vector is None:
        return []
    if hasattr(vector, "size") and hasattr(vector, "at"):
        return [vector.at(i) for i in range(int(vector.size()))]
    try:
        return [vector[i] for i in range(len(vector))]
    except Exception:
        return []


def finite_float(value):
    try:
        value = float(value)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def positive_finite(values):
    cleaned = []
    for value in values:
        value = finite_float(value)
        if value is not None and value > 0.0:
            cleaned.append(value)
    return cleaned


def leading_index(energies):
    """Index of the highest positive, finite energy, or None."""
    best, best_energy = None, 0.0
    for index, value in enumerate(energies):
        value = finite_float(value)
        if value is not None and value > best_energy:
            best, best_energy = index, value
    return best


def compute_neff(energies):
    energies = positive_finite(energies)
    denominator = sum(energy * energy for energy in energies)
    return (sum(energies) ** 2 / denominator) if denominator > 0.0 else 0.0


def make_empty_sample(path):
    return {
        "path": path,
        "n_events": 0,
        "hcal_multiplicity": [],
        "ecal_multiplicity": [],
        "hcal_neff": [],
        "ecal_neff": [],
        "hcal_eta": [],
        "ecal_eta": [],
        "hcal_energy": [],
        "ecal_energy": [],
        "gen_energy": [],
        "gen_eta": [],
        "gen_radius": [],
        "hcal_lead_time": [],
        "hcal_lead_eta": [],
        "hcal_lead_seed_eta": [],
        "hcal_lead_time_skipped": 0,
    }


def build_file_lookup(input_dir, prefix, pf_name, energies, etas, debug=False):
    lookup = {}
    energy_keys = {key_float(x) for x in energies}
    eta_keys = {key_float(x) for x in etas}
    expression = re.compile(
        rf"^{re.escape(prefix + pf_name)}_SinglePiCloseByE(?P<energy>[^_]+)"
        rf"_eta(?P<eta>[^.]+)\.root$"
    )
    pattern = os.path.join(input_dir, f"{prefix}{pf_name}_SinglePiCloseByE*_eta*.root")
    candidates = sorted(glob.glob(pattern))
    if debug:
        print(f"Scanning {pattern}: found {len(candidates)} candidate files")
    for path in candidates:
        match = expression.match(os.path.basename(path))
        if not match:
            if debug:
                print(f"  Ignoring filename that does not match convention: {path}")
            continue
        try:
            energy = key_float(decode_number_tag(match.group("energy")))
            eta = key_float(decode_number_tag(match.group("eta")))
        except ValueError:
            print(f"WARNING: could not decode energy/eta from {path}")
            continue
        if energy not in energy_keys or eta not in eta_keys:
            continue
        key = (energy, eta)
        if key in lookup:
            print(f"WARNING: duplicate sample {key}; keeping {lookup[key]}, skipping {path}")
            continue
        lookup[key] = path
        if debug:
            print(f"  E={energy:g} GeV, eta={eta:g}: {path}")
    return lookup


def selected_generator_pions(event):
    branches = ["gen_energy", "gen_eta", "gen_pdgId", "gen_status",
                "gen_vx", "gen_vy", "gen_vz"]
    values = [vector_to_list(getattr(event, name)) for name in branches]
    pions = []
    for index in range(min(map(len, values))):
        energy = finite_float(values[0][index])
        eta = finite_float(values[1][index])
        vx = finite_float(values[4][index])
        vy = finite_float(values[5][index])
        vz = finite_float(values[6][index])
        try:
            pdgid = int(values[2][index])
            status = int(values[3][index])
        except Exception:
            continue
        if (
            energy is None
            or eta is None
            or vx is None
            or vy is None
            or vz is None
            or energy <= 0.0
            or abs(pdgid) != 211
            or status != 1
        ):
            continue
        pions.append(
            {
                "energy": energy,
                "eta": eta,
                "radius": math.hypot(vx, vy),
            }
        )
    return pions


def read_sample(path, tree_name):
    sample = make_empty_sample(path)
    source = ROOT.TFile.Open(path)
    if not source or source.IsZombie():
        print(f"WARNING: could not open {path}")
        return None
    tree = source.Get(tree_name)
    if not tree:
        print(f"WARNING: missing tree {tree_name} in {path}")
        source.Close()
        return None

    required = [
        "hcal_energy", "hcal_eta", "ecal_energy", "ecal_eta",
        "hcal_time", "hcal_seed_eta",
        "gen_energy", "gen_eta", "gen_pdgId", "gen_status",
        "gen_vx", "gen_vy", "gen_vz",
    ]
    missing = [name for name in required if not tree.GetBranch(name)]
    if missing:
        print(f"WARNING: skipping {path}; missing branches: {', '.join(missing)}")
        source.Close()
        return None

    for event in tree:
        raw_hcal_energy = vector_to_list(event.hcal_energy)
        raw_ecal_energy = vector_to_list(event.ecal_energy)
        hcal_energy = positive_finite(raw_hcal_energy)
        ecal_energy = positive_finite(raw_ecal_energy)
        hcal_eta = [x for x in map(finite_float, vector_to_list(event.hcal_eta))
                    if x is not None]
        ecal_eta = [x for x in map(finite_float, vector_to_list(event.ecal_eta))
                    if x is not None]

        sample["n_events"] += 1
        # Multiplicity counts every stored PF cluster, matching the no-cut
        # definition used in the earlier plotting scripts.  Energy-weighted
        # quantities below use only positive, finite energies.
        sample["hcal_multiplicity"].append(len(raw_hcal_energy))
        sample["ecal_multiplicity"].append(len(raw_ecal_energy))
        sample["hcal_neff"].append(compute_neff(hcal_energy))
        sample["ecal_neff"].append(compute_neff(ecal_energy))
        sample["hcal_energy"].extend(hcal_energy)
        sample["ecal_energy"].extend(ecal_energy)
        sample["hcal_eta"].extend(hcal_eta)
        sample["ecal_eta"].extend(ecal_eta)

        # Leading (highest-energy) HCAL PF cluster: one entry per event at most.
        lead = leading_index(raw_hcal_energy)
        if lead is not None:
            raw_hcal_eta = vector_to_list(event.hcal_eta)
            raw_hcal_time = vector_to_list(event.hcal_time)
            raw_hcal_seed_eta = vector_to_list(event.hcal_seed_eta)
            time = finite_float(raw_hcal_time[lead]) if lead < len(raw_hcal_time) else None
            if time is None or time == -999.0:
                sample["hcal_lead_time_skipped"] += 1
            else:
                sample["hcal_lead_time"].append(time)
            if lead < len(raw_hcal_eta) and lead < len(raw_hcal_seed_eta):
                cluster_eta = finite_float(raw_hcal_eta[lead])
                seed_eta = finite_float(raw_hcal_seed_eta[lead])
                if cluster_eta is not None and seed_eta is not None:
                    sample["hcal_lead_eta"].append(cluster_eta)
                    sample["hcal_lead_seed_eta"].append(seed_eta)

        for pion in selected_generator_pions(event):
            sample["gen_energy"].append(pion["energy"])
            sample["gen_eta"].append(pion["eta"])
            sample["gen_radius"].append(pion["radius"])

    source.Close()
    print(
        f"  Leading HCAL cluster time: {len(sample['hcal_lead_time'])} events used, "
        f"{sample['hcal_lead_time_skipped']} skipped (time = -999)"
    )
    return sample


def mean_or_zero(values):
    return sum(values) / len(values) if values else 0.0


def install_heatmap_palette():
    stops, red, green, blue = HEATMAP_PALETTE
    ROOT.TColor.CreateGradientColorTable(
        len(stops),
        array("d", stops),
        array("d", red),
        array("d", green),
        array("d", blue),
        255,
    )
    ROOT.gStyle.SetNumberContours(255)


def make_canvas(name, width, height):
    canvas = ROOT.TCanvas(name, "", int(width), int(height))
    canvas.SetCanvasSize(int(width), int(height))
    canvas.Modified()
    canvas.Update()
    return canvas


def print_pdf_page(canvas, pdf_path):
    canvas.cd()
    canvas.Modified()
    canvas.Update()
    canvas.Print(pdf_path)


def draw_page_header(canvas, tag, title):
    header = ROOT.TPad(f"header_{tag}", "", 0.0, 0.93, 1.0, 1.0)
    header.SetFillStyle(0)
    header.Draw()
    header.cd()
    label = ROOT.TLatex()
    label.SetNDC()
    label.SetTextAlign(22)
    label.SetTextFont(62)
    label.SetTextSize(0.47)
    label.DrawLatex(0.5, 0.43, title)
    canvas.cd()
    return header, label


def make_scan_hist(name, title, energies, etas):
    hist = ROOT.TH2D(name, title, len(etas), 0, len(etas),
                     len(energies), 0, len(energies))
    hist.SetDirectory(0)
    for index, eta in enumerate(etas, 1):
        hist.GetXaxis().SetBinLabel(index, format_value(eta))
    for index, energy in enumerate(energies, 1):
        hist.GetYaxis().SetBinLabel(index, format_value(energy))
    return hist


def fill_scan_hist(hist, data, energies, etas, value_getter):
    values = []
    for iy, energy in enumerate(energies, 1):
        for ix, eta in enumerate(etas, 1):
            sample = data.get((key_float(energy), key_float(eta)))
            if sample is None:
                continue
            value = float(value_getter(sample))
            hist.SetBinContent(ix, iy, value)
            values.append(value)
    return values


def heatmap_range(values, start_at_zero=True):
    finite = [float(x) for x in values if math.isfinite(float(x))]
    if not finite:
        return (0.0, 1.0)
    maximum = max(finite)
    minimum = min(finite)
    if start_at_zero:
        # Preserve contrast when the mean multiplicity is far below one.
        # Only fall back to [0, 1] when every populated cell is exactly zero.
        return (0.0, maximum * 1.08 if maximum > 0.0 else 1.0)
    if math.isclose(minimum, maximum, rel_tol=1e-12, abs_tol=1e-12):
        padding = max(abs(maximum) * 0.08, 0.01)
    else:
        padding = 0.06 * (maximum - minimum)
    return (minimum - padding, maximum + padding)


def style_heatmap(hist, panel_title, z_title, z_range, text_format):
    hist.SetTitle(panel_title)
    hist.SetMinimum(z_range[0])
    hist.SetMaximum(z_range[1])
    hist.GetXaxis().SetTitle("Generated pion #eta")
    hist.GetYaxis().SetTitle("Generated pion energy [GeV]")
    hist.GetZaxis().SetTitle(z_title)
    for axis in [hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()]:
        axis.CenterTitle()
    hist.GetXaxis().SetTitleSize(0.050)
    hist.GetYaxis().SetTitleSize(0.050)
    hist.GetZaxis().SetTitleSize(0.042)
    hist.GetXaxis().SetLabelSize(0.043)
    hist.GetYaxis().SetLabelSize(0.043)
    hist.GetZaxis().SetLabelSize(0.035)
    hist.GetXaxis().SetTitleOffset(1.02)
    hist.GetYaxis().SetTitleOffset(1.14)
    hist.GetZaxis().SetTitleOffset(1.35)
    hist.SetMarkerSize(1.55)
    ROOT.gStyle.SetPaintTextFormat(text_format)


def make_two_heatmap_page(data, energies, etas, detector, pdf_path):
    detector_upper = detector.upper()
    key_prefix = detector.lower()
    specs = [
        (
            "multiplicity",
            f"Average Number of {detector_upper} PF Clusters per Event",
            "Average # clusters/event",
            lambda sample: mean_or_zero(sample[f"{key_prefix}_multiplicity"]),
        ),
        (
            "neff",
            f"Average Effective Number of {detector_upper} PF Clusters",
            "#LTN_{eff}#GT = #LT(#Sigma E)^{2}/#Sigma E^{2}#GT",
            lambda sample: mean_or_zero(sample[f"{key_prefix}_neff"]),
        ),
    ]
    canvas = make_canvas(f"c_{key_prefix}_heatmaps", 1900, 900)
    header, page_label = draw_page_header(
        canvas, f"{key_prefix}_heatmaps", f"{detector_upper} PF-cluster summary"
    )
    canvas.cd()
    body = ROOT.TPad(f"body_{key_prefix}_heatmaps", "", 0.0, 0.0, 1.0, 0.93)
    body.SetFillStyle(0)
    body.Draw()
    body.cd()
    body.Divide(2, 1, 0.003, 0.003)
    keep = [header, page_label, body]
    for index, (metric, title, z_title, getter) in enumerate(specs, 1):
        hist = make_scan_hist(
            f"h2_{key_prefix}_{metric}", title, energies, etas
        )
        values = fill_scan_hist(hist, data, energies, etas, getter)
        pad = body.cd(index)
        pad.SetLeftMargin(0.13)
        pad.SetRightMargin(0.18)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.10)
        pad.SetTicks(1, 1)
        style_heatmap(hist, title, z_title, heatmap_range(values), "4.2f")
        hist.Draw("COLZ TEXT")
        keep.append(hist)
    print_pdf_page(canvas, pdf_path)
    print(f"Added page: {detector_upper} multiplicity heatmaps")


def make_generator_heatmap_page(data, energies, etas, pdf_path):
    specs = [
        (
            "count",
            "Total Number of Generator Pions",
            "# generator pions",
            lambda sample: len(sample["gen_energy"]),
            "4.0f",
            True,
        ),
        (
            "eta",
            "Average Generator-Pion #eta",
            "Average pion #eta",
            lambda sample: mean_or_zero(sample["gen_eta"]),
            "5.3f",
            False,
        ),
        (
            "radius",
            "Average Generator-Pion Production R",
            "Average pion R [cm]",
            lambda sample: mean_or_zero(sample["gen_radius"]),
            "5.2f",
            False,
        ),
    ]
    canvas = make_canvas("c_generator_heatmaps", 2600, 850)
    header, page_label = draw_page_header(
        canvas, "generator_heatmaps", "Generator-pion sample validation"
    )
    canvas.cd()
    body = ROOT.TPad("body_generator_heatmaps", "", 0.0, 0.0, 1.0, 0.93)
    body.SetFillStyle(0)
    body.Draw()
    body.cd()
    body.Divide(3, 1, 0.003, 0.003)
    keep = [header, page_label, body]
    for index, (metric, title, z_title, getter, text_format, zero_min) in enumerate(specs, 1):
        hist = make_scan_hist(f"h2_gen_{metric}", title, energies, etas)
        values = fill_scan_hist(hist, data, energies, etas, getter)
        pad = body.cd(index)
        pad.SetLeftMargin(0.14)
        pad.SetRightMargin(0.19)
        pad.SetBottomMargin(0.14)
        pad.SetTopMargin(0.10)
        pad.SetTicks(1, 1)
        style_heatmap(
            hist, title, z_title, heatmap_range(values, zero_min), text_format
        )
        hist.Draw("COLZ TEXT")
        keep.append(hist)
    print_pdf_page(canvas, pdf_path)
    print("Added page: generator-pion validation heatmaps")


def padded_data_range(values, fallback, lower_bound=None, fractional_padding=0.08):
    finite = [float(value) for value in values if math.isfinite(float(value))]
    if not finite:
        return fallback
    minimum, maximum = min(finite), max(finite)
    if math.isclose(minimum, maximum, rel_tol=1e-12, abs_tol=1e-12):
        padding = max(abs(maximum) * fractional_padding, 0.05)
    else:
        padding = fractional_padding * (maximum - minimum)
    xmin, xmax = minimum - padding, maximum + padding
    if lower_bound is not None:
        xmin = max(float(lower_bound), xmin)
    if xmax <= xmin:
        xmax = xmin + 1.0
    return xmin, xmax


def make_hist(name, bins, xmin, xmax):
    hist = ROOT.TH1D(name, "", int(bins), float(xmin), float(xmax))
    hist.SetDirectory(0)
    hist.SetLineColor(ROOT.kAzure + 2)
    hist.SetLineWidth(2)
    hist.SetFillColorAlpha(ROOT.kAzure + 1, 0.28)
    return hist


def style_distribution_hist(hist, title, x_title, y_title, ymax):
    hist.SetTitle(title)
    hist.SetMinimum(0.0)
    hist.SetMaximum(max(1.0, ymax))
    hist.GetXaxis().SetTitle(x_title)
    hist.GetYaxis().SetTitle(y_title)
    hist.GetXaxis().SetTitleSize(0.050)
    hist.GetYaxis().SetTitleSize(0.050)
    hist.GetXaxis().SetLabelSize(0.043)
    hist.GetYaxis().SetLabelSize(0.043)
    hist.GetXaxis().SetTitleOffset(1.08)
    hist.GetYaxis().SetTitleOffset(1.20)


def selected_grid_values(data, grid_energies, grid_etas, field):
    values = []
    for energy in grid_energies:
        for eta in grid_etas:
            sample = data.get((key_float(energy), key_float(eta)))
            if sample:
                values.extend(sample[field])
    return values


def distribution_axis(data, grid_energies, grid_etas, field, expected_eta=None):
    values = selected_grid_values(data, grid_energies, grid_etas, field)
    if field.endswith("_multiplicity"):
        maximum = max([int(value) for value in values] + [1])
        maximum = max(2, int(math.ceil(maximum * 1.10)))
        return maximum + 1, -0.5, maximum + 0.5
    if field in ("hcal_eta", "ecal_eta"):
        return None
    if field in ("hcal_energy", "ecal_energy"):
        xmin, xmax = padded_data_range(values, (0.0, 1.0), lower_bound=0.0, fractional_padding=0.10)
        return None, xmin, xmax
    if field == "gen_energy":
        xmin = max(0.0, min(grid_energies) * 0.75)
        xmax = max(grid_energies) * 1.20
        return None, xmin, xmax
    if field == "gen_radius":
        xmin, xmax = padded_data_range(values, (150.0, 210.0), lower_bound=0.0, fractional_padding=0.08)
        return None, xmin, xmax
    if field == "gen_eta" and expected_eta is not None:
        return None, expected_eta - 0.005, expected_eta + 0.005
    return None, *padded_data_range(values, (0.0, 1.0))


def make_distribution_grid(
    data,
    grid_energies,
    grid_etas,
    field,
    page_title,
    x_title,
    y_title,
    bins,
    pdf_path,
    cluster_eta_range=None,
):
    nrows, ncols = len(grid_energies), len(grid_etas)
    tag = sanitize_name(field)
    canvas = make_canvas(f"c_grid_{tag}", 760 * ncols, 590 * nrows)
    header, page_label = draw_page_header(canvas, tag, page_title)
    canvas.cd()
    grid = ROOT.TPad(f"body_grid_{tag}", "", 0.0, 0.0, 1.0, 0.93)
    grid.SetFillStyle(0)
    grid.Draw()
    grid.cd()
    grid.Divide(ncols, nrows, 0.002, 0.002)
    keep = [header, page_label, grid]

    common_axis = None
    if field != "gen_eta":
        common_axis = distribution_axis(data, grid_energies, grid_etas, field)

    for row, energy in enumerate(grid_energies):
        for col, eta in enumerate(grid_etas):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12)
            pad.SetBottomMargin(0.14)
            pad.SetTicks(1, 1)
            sample = data.get((key_float(energy), key_float(eta)))
            if sample is None:
                missing = ROOT.TLatex()
                missing.SetNDC()
                missing.SetTextAlign(22)
                missing.SetTextSize(0.07)
                missing.DrawLatex(0.5, 0.5, "missing sample")
                keep.append(missing)
                continue

            if field == "gen_eta":
                _, xmin, xmax = distribution_axis(
                    data, grid_energies, grid_etas, field, expected_eta=eta
                )
                nbins = bins
            elif field in ("hcal_eta", "ecal_eta") and cluster_eta_range:
                xmin, xmax = cluster_eta_range
                nbins = bins
            else:
                suggested_bins, xmin, xmax = common_axis
                nbins = suggested_bins if suggested_bins is not None else bins

            hist = make_hist(
                f"h_{tag}_E{number_tag(energy)}_eta{number_tag(eta)}",
                nbins,
                xmin,
                xmax,
            )
            for value in sample[field]:
                hist.Fill(value)
            ymax = 1.25 * hist.GetMaximum()
            title = f"E = {format_value(energy)} GeV, #eta = {format_value(eta)}"
            style_distribution_hist(hist, title, x_title, y_title, ymax)
            hist.Draw("HIST")

            entries = ROOT.TLatex()
            entries.SetNDC()
            entries.SetTextAlign(33)
            entries.SetTextFont(42)
            entries.SetTextSize(0.043)
            entries.DrawLatex(0.94, 0.82, f"Entries = {len(sample[field])}")
            keep.extend([hist, entries])

    print_pdf_page(canvas, pdf_path)
    print(f"Added page: {page_title}")


def make_scatter_grid(
    data,
    grid_energies,
    grid_etas,
    x_field,
    y_field,
    page_title,
    x_title,
    y_title,
    pdf_path,
):
    nrows, ncols = len(grid_energies), len(grid_etas)
    tag = sanitize_name(f"{y_field}_vs_{x_field}")
    canvas = make_canvas(f"c_grid_{tag}", 760 * ncols, 590 * nrows)
    header, page_label = draw_page_header(canvas, tag, page_title)
    canvas.cd()
    grid = ROOT.TPad(f"body_grid_{tag}", "", 0.0, 0.0, 1.0, 0.93)
    grid.SetFillStyle(0)
    grid.Draw()
    grid.cd()
    grid.Divide(ncols, nrows, 0.002, 0.002)
    keep = [header, page_label, grid]

    # Common range on both axes so the y = x reference line is the diagonal.
    values = selected_grid_values(data, grid_energies, grid_etas, x_field)
    values += selected_grid_values(data, grid_energies, grid_etas, y_field)
    vmin, vmax = padded_data_range(values, (-1.0, 1.0))

    for row, energy in enumerate(grid_energies):
        for col, eta in enumerate(grid_etas):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.14)
            pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12)
            pad.SetBottomMargin(0.14)
            pad.SetTicks(1, 1)
            sample = data.get((key_float(energy), key_float(eta)))
            if sample is None:
                missing = ROOT.TLatex()
                missing.SetNDC()
                missing.SetTextAlign(22)
                missing.SetTextSize(0.07)
                missing.DrawLatex(0.5, 0.5, "missing sample")
                keep.append(missing)
                continue

            frame = pad.DrawFrame(vmin, vmin, vmax, vmax)
            title = f"E = {format_value(energy)} GeV, #eta = {format_value(eta)}"
            frame.SetTitle(title)
            frame.GetXaxis().SetTitle(x_title)
            frame.GetYaxis().SetTitle(y_title)
            frame.GetXaxis().SetTitleSize(0.050)
            frame.GetYaxis().SetTitleSize(0.050)
            frame.GetXaxis().SetLabelSize(0.043)
            frame.GetYaxis().SetLabelSize(0.043)
            frame.GetXaxis().SetTitleOffset(1.08)
            frame.GetYaxis().SetTitleOffset(1.20)

            diagonal = ROOT.TLine(vmin, vmin, vmax, vmax)
            diagonal.SetLineStyle(2)
            diagonal.SetLineColor(ROOT.kGray + 2)
            diagonal.Draw()

            x_values, y_values = sample[x_field], sample[y_field]
            n_points = len(x_values)
            if n_points > 0:
                graph = ROOT.TGraph(n_points, array("d", x_values), array("d", y_values))
                graph.SetMarkerStyle(20)
                graph.SetMarkerSize(0.6)
                graph.SetMarkerColor(ROOT.kAzure + 2)
                graph.Draw("P SAME")
                keep.append(graph)

            entries = ROOT.TLatex()
            entries.SetNDC()
            entries.SetTextAlign(33)
            entries.SetTextFont(42)
            entries.SetTextSize(0.043)
            entries.DrawLatex(0.94, 0.82, f"Entries = {n_points}")
            keep.extend([frame, diagonal, entries])

    print_pdf_page(canvas, pdf_path)
    print(f"Added page: {page_title}")


def print_sample_summary(data, energies, etas):
    print("\nPer-sample summary")
    print(
        "E [GeV]   eta     events  gen pi  <N_HCAL>  <N_ECAL>  "
        "sum HCAL E/event  sum ECAL E/event  <pion R> [cm]"
    )
    for energy in energies:
        for eta in etas:
            sample = data.get((key_float(energy), key_float(eta)))
            if sample is None:
                print(f"{energy:7g}  {eta:6g}   MISSING")
                continue
            n_events = sample["n_events"]
            hcal_e_per_event = (
                sum(sample["hcal_energy"]) / n_events if n_events else 0.0
            )
            ecal_e_per_event = (
                sum(sample["ecal_energy"]) / n_events if n_events else 0.0
            )
            print(
                f"{energy:7g}  {eta:6g}  {n_events:7d}  "
                f"{len(sample['gen_energy']):6d}  "
                f"{mean_or_zero(sample['hcal_multiplicity']):8.3f}  "
                f"{mean_or_zero(sample['ecal_multiplicity']):8.3f}  "
                f"{hcal_e_per_event:16.3f}  {ecal_e_per_event:16.3f}  "
                f"{mean_or_zero(sample['gen_radius']):13.3f}"
            )


def validate_requested_grid(all_values, grid_values, label):
    available = {key_float(value) for value in all_values}
    missing = [value for value in grid_values if key_float(value) not in available]
    if missing:
        raise ValueError(
            f"Every --grid-{label} value must also appear in --{label}: {missing}"
        )


def main():
    args = parse_args()
    energies = parse_float_list(args.energies, DEFAULT_ENERGIES)
    etas = parse_float_list(args.etas, DEFAULT_ETAS)
    grid_energies = parse_float_list(args.grid_energies, DEFAULT_GRID_ENERGIES)
    grid_etas = parse_float_list(args.grid_etas, DEFAULT_GRID_ETAS)
    validate_requested_grid(energies, grid_energies, "energies")
    validate_requested_grid(etas, grid_etas, "etas")

    lookup = build_file_lookup(
        args.input_dir,
        args.prefix,
        args.pf_name,
        energies,
        etas,
        args.debug_files,
    )
    if not lookup:
        raise RuntimeError(
            "No matching ntuples found. Check --input-dir, --prefix, --pf-name, "
            "and the filename convention."
        )

    data = {}
    for energy in energies:
        for eta in etas:
            key = (key_float(energy), key_float(eta))
            path = lookup.get(key)
            if path is None:
                print(
                    f"WARNING: missing ntuple for E={energy:g} GeV, eta={eta:g}; "
                    "the corresponding heatmap cell will remain zero."
                )
                continue
            print(f"Reading E={energy:g} GeV, eta={eta:g}: {path}")
            sample = read_sample(path, args.tree_name)
            if sample is not None:
                data[key] = sample

    if not data:
        raise RuntimeError("No readable samples remain after checking the ROOT files.")
    print_sample_summary(data, energies, etas)

    os.makedirs(args.output_dir, exist_ok=True)
    pdf_path = os.path.join(args.output_dir, args.output_pdf)
    install_heatmap_palette()
    opener = make_canvas("pdf_open", 1900, 900)
    opener.Print(pdf_path + "[")

    # 1-2: PF-cluster multiplicity/effective-multiplicity heatmaps.
    make_two_heatmap_page(data, energies, etas, "HCAL", pdf_path)
    make_two_heatmap_page(data, energies, etas, "ECAL", pdf_path)

    # 3-4: PF-cluster multiplicity distributions.
    make_distribution_grid(
        data, grid_energies, grid_etas, "hcal_multiplicity",
        "HCAL PF-cluster multiplicity", "Number of HCAL PF clusters", "Events",
        10, pdf_path,
    )
    make_distribution_grid(
        data, grid_energies, grid_etas, "ecal_multiplicity",
        "ECAL PF-cluster multiplicity", "Number of ECAL PF clusters", "Events",
        10, pdf_path,
    )

    cluster_eta_range = (args.cluster_eta_min, args.cluster_eta_max)

    # 5-6: PF-cluster eta distributions.
    make_distribution_grid(
        data, grid_energies, grid_etas, "hcal_eta",
        "HCAL PF-cluster #eta", "HCAL PF-cluster #eta", "HCAL PF clusters",
        args.cluster_eta_bins, pdf_path, cluster_eta_range,
    )
    make_distribution_grid(
        data, grid_energies, grid_etas, "ecal_eta",
        "ECAL PF-cluster #eta", "ECAL PF-cluster #eta", "ECAL PF clusters",
        args.cluster_eta_bins, pdf_path, cluster_eta_range,
    )

    # 7-8: PF-cluster energy distributions.
    make_distribution_grid(
        data, grid_energies, grid_etas, "hcal_energy",
        "HCAL PF-cluster energy", "HCAL PF-cluster energy [GeV]", "HCAL PF clusters",
        args.cluster_energy_bins, pdf_path,
    )
    make_distribution_grid(
        data, grid_energies, grid_etas, "ecal_energy",
        "ECAL PF-cluster energy", "ECAL PF-cluster energy [GeV]", "ECAL PF clusters",
        args.cluster_energy_bins, pdf_path,
    )

    # 9: Generator-pion energy distributions.
    make_distribution_grid(
        data, grid_energies, grid_etas, "gen_energy",
        "Generator-pion energy", "Generator-pion energy [GeV]", "Generator pions",
        args.gen_energy_bins, pdf_path,
    )

    # 10: Generator-pion count, eta, and transverse production-radius heatmaps.
    make_generator_heatmap_page(data, energies, etas, pdf_path)

    # 11-12: Generator-pion eta and transverse production-radius distributions.
    make_distribution_grid(
        data, grid_energies, grid_etas, "gen_eta",
        "Generator-pion #eta", "Generator-pion #eta", "Generator pions",
        args.gen_eta_bins, pdf_path,
    )
    make_distribution_grid(
        data, grid_energies, grid_etas, "gen_radius",
        "Generator-pion transverse production radius",
        "Generator-pion R = #sqrt{v_{x}^{2}+v_{y}^{2}} [cm]", "Generator pions",
        args.gen_radius_bins, pdf_path,
    )

    # 13: Leading HCAL PF-cluster time distributions (time = -999 skipped).
    make_distribution_grid(
        data, grid_energies, grid_etas, "hcal_lead_time",
        "Leading HCAL PF-cluster time", "Leading HCAL PF-cluster time [ns]", "Events",
        args.cluster_time_bins, pdf_path,
    )

    # 14: Leading HCAL PF-cluster eta versus its seed eta.
    make_scatter_grid(
        data, grid_energies, grid_etas, "hcal_lead_seed_eta", "hcal_lead_eta",
        "Leading HCAL PF-cluster #eta vs seed #eta",
        "Seed #eta", "Leading HCAL PF-cluster #eta", pdf_path,
    )

    opener.Modified()
    opener.Update()
    opener.Print(pdf_path + "]")
    print(f"\nSaved combined PDF: {pdf_path}")


if __name__ == "__main__":
    main()