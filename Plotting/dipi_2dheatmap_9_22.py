#!/usr/bin/env python3
"""Make dipion PF HCAL DR/DT heatmaps and distribution/efficiency comparisons.

Every quantity is written as consecutive no-cut and cluster-energy-cut pages.
Heatmaps contain standardPF, seedTimingPF, and seedTimingPF/standardPF panels.
One-dimensional comparisons overlay the algorithms without ratio subpanels. The
default cut is E_cluster >= sqrt(E_pi) GeV, using the mean energy of the two
selected generator pions in each event; --cluster-energy-cut selects a fixed
threshold instead. Four additional no-cut timing pages combine the times of
the two highest-energy clusters per event, skipping invalid times, for each
PF algorithm at DR=0.2 and 0.3 with DT=0,1,2,3,4,5 ns.
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

HEATMAP_VALUE_FORMAT = "4.2f"
HEATMAP_PERCENT_FORMAT = "5.1f"
CL_1SIGMA = 0.682689492137

DEFAULT_DR_VALUES = [0.1, 0.2, 0.3, 0.4, 0.5]
DEFAULT_DT_VALUES = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]
DEFAULT_DR_DIST_VALUES = [0.1, 0.2, 0.3]
DEFAULT_DT_DIST_VALUES = [1.0, 3.0, 5.0]
DEFAULT_PF_NAMES = ["standardPF", "seedTimingPF"]

SCATTER_DR_VALUES = [0.2, 0.3]
SCATTER_DT_VALUES = [3.0, 4.0, 5.0]

TIMING_DR_VALUES = [0.2, 0.3]
TIMING_DT_VALUES = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]


def parse_float_list(value, default=None):
    if value is None:
        return default
    if isinstance(value, list):
        return [float(v) for v in value]
    value = str(value).strip()
    return default if not value else [float(x) for x in value.split(",") if x.strip()]


def parse_args():
    p = argparse.ArgumentParser(description="Make PF-cluster DR/DT comparison plots.")
    p.add_argument("--input-dir", "--inputdir", dest="input_dir",
                   default="/eos/user/c/chtong/Public/Rereco/Dipion_eta0.4_4ns/")
    p.add_argument("--output-dir", "--outputdir", dest="output_dir",
                   default="/eos/user/c/chtong/Public/Rereco/Dipion_eta0.4_4ns/")
    p.add_argument("--prefix", default="pfObjectsNtuple_")
    p.add_argument("--tree-name", default="pfObjectsNtupler/pfTree")
    p.add_argument("--pf-names", nargs="+", default=DEFAULT_PF_NAMES)
    p.add_argument("--dr-values", default=",".join(map(str, DEFAULT_DR_VALUES)))
    p.add_argument("--dt-values", default=",".join(map(str, DEFAULT_DT_VALUES)))
    p.add_argument("--dr-dist-values", default=",".join(map(str, DEFAULT_DR_DIST_VALUES)))
    p.add_argument("--dt-dist-values", default=",".join(map(str, DEFAULT_DT_DIST_VALUES)))

    p.add_argument("--cluster-time-branch", default="hcal_time")
    p.add_argument("--cluster-time-bin-width", type=float, default=1.0,
                   help="Timing bin width in ns (default: 1, centered on integer times); common range is derived from valid times.")
    p.add_argument("--energy-branch", default="hcal_energy")
    p.add_argument("--cluster-eta-branch", default="hcal_eta")
    p.add_argument("--cluster-phi-branch", default="hcal_phi")
    p.add_argument("--all-pfrh-energy-branch", default="all_hbhe_pfrh_energy")
    p.add_argument("--all-pfrh-cluster-index-branch", default="all_hbhe_pfrh_clusterIdx")
    p.add_argument("--clustered-pfrh-energy-branch", default="hbhe_pfrh_energyFracInCluster")
    p.add_argument("--clustered-pfrh-cluster-index-branch", default="hbhe_pfrh_clusterIdx")
    p.add_argument("--cluster-seed-eta-branch", default="hcal_seed_eta")
    p.add_argument("--cluster-seed-phi-branch", default="hcal_seed_phi")
    p.add_argument("--gen-energy-branch", default="gen_energy")
    p.add_argument("--gen-eta-branch", default="gen_eta")
    p.add_argument("--gen-phi-branch", default="gen_phi")
    p.add_argument("--gen-pdgid-branch", default="gen_pdgId")
    p.add_argument("--gen-status-branch", default="gen_status")

    p.add_argument(
        "--cluster-energy-cut",
        type=float,
        default=None,
        help="Fixed GeV cut. Default: event-by-event sqrt(E_pi) GeV.",
    )
    p.add_argument("--ncl-max", type=float, default=6.0)
    p.add_argument("--neff-max", type=float, default=4.0)
    p.add_argument("--neff-bins", type=int, default=20)
    p.add_argument("--rechit-percent-bins", type=int, default=25)
    p.add_argument("--reco-energy-bins", type=int, default=20)
    p.add_argument("--reco-energy-max", type=float, default=80.0)
    p.add_argument("--match-dr-bins", type=int, default=25)
    p.add_argument("--match-dr-max", type=float, default=0.3)
    # Seed positions sit further from the generator pion than cluster
    # positions do, so the seed page needs a wider axis to avoid piling up
    # in the overflow bin.
    p.add_argument("--seed-match-dr-max", type=float, default=0.5)
    p.add_argument("--output-pdf", default=None)
    p.add_argument("--output-root", default=None)
    p.add_argument("--debug-files", action="store_true")
    return p.parse_args()


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
        return [v]
    except Exception:
        return []


def finite_float(value):
    try:
        value = float(value)
    except Exception:
        return None
    return value if math.isfinite(value) else None


def clean_positive_energies(values):
    return [x for x in (finite_float(v) for v in values) if x is not None and x > 0.0]


def compute_neff(energies):
    energies = clean_positive_energies(energies)
    denom = sum(e * e for e in energies)
    return (sum(energies) ** 2 / denom) if denom > 0.0 else 0.0


def delta_phi(phi1, phi2):
    return math.atan2(math.sin(phi1 - phi2), math.cos(phi1 - phi2))


def delta_r(eta1, phi1, eta2, phi2):
    return math.hypot(eta1 - eta2, delta_phi(phi1, phi2))


def branch_exists(tree, name):
    return bool(tree.GetBranch(name))


def cut_label(args, root_text=False):
    if args.cluster_energy_cut is None:
        return (
            "HCAL clusters with E #geq #sqrt{E_{#pi}} GeV"
            if root_text
            else "HCAL clusters with E >= sqrt(E_pi) GeV"
        )
    ge = "#geq" if root_text else ">="
    return f"HCAL clusters with E {ge} {args.cluster_energy_cut:g} GeV"


def event_cut_threshold(gen_pions, args):
    if args.cluster_energy_cut is not None:
        return args.cluster_energy_cut
    if not gen_pions:
        return float("inf")
    return math.sqrt(sum(p["energy"] for p in gen_pions) / len(gen_pions))


# -----------------------------------------------------------------------------
# File discovery
# -----------------------------------------------------------------------------

def parse_dr_dt_from_filename(basename):
    """Read DT/DR from the new names, retaining support for the old names.

    Example: pfObjectsNtuple_seedTimingPF_DiPi_Z17.8_offsetFirst0.0_DT5.0_DR0.2_pT20.0_eta0.1.root
    The Z, offsetFirst, pT and eta tags do not alter DT/DR extraction.
    """
    match = re.search(r"_DT([-+]?\d+(?:\.\d+)?)_DR([-+]?\d+(?:\.\d+)?)", basename)
    return (float(match.group(2)), float(match.group(1))) if match else (None, None)


def pf_from_filename(basename, prefix, pf_names):
    if not basename.startswith(prefix):
        return None
    rest = basename[len(prefix):]
    for pf_name in sorted(pf_names, key=len, reverse=True):
        if rest.startswith(pf_name + "_") or rest == pf_name + ".root":
            return pf_name
    return None


def build_file_lookup(input_dir, prefix, pf_names, dr_values, dt_values, debug=False):
    lookup = {}
    dr_set, dt_set = set(map(key_float, dr_values)), set(map(key_float, dt_values))
    pattern = os.path.join(input_dir, f"{prefix}*.root")
    candidates = sorted(glob.glob(pattern))
    if debug:
        print(f"Scanning {pattern}: {len(candidates)} candidates")
    for path in candidates:
        base = os.path.basename(path)
        pf_name = pf_from_filename(base, prefix, pf_names)
        dr, dt = parse_dr_dt_from_filename(base)
        if pf_name is None or dr is None:
            continue
        dr_key, dt_key = key_float(dr), key_float(dt)
        if dr_key not in dr_set or dt_key not in dt_set:
            continue
        key = (pf_name, dr_key, dt_key)
        if key in lookup:
            print(f"WARNING: duplicate {key}; keeping {lookup[key]}, skipping {path}")
            continue
        lookup[key] = path
        if debug:
            print(f"  {base} -> {key}")
    return lookup


# -----------------------------------------------------------------------------
# Event calculations and histogram booking
# -----------------------------------------------------------------------------

def generator_pions(event, args):
    branches = [args.gen_energy_branch, args.gen_eta_branch, args.gen_phi_branch,
                args.gen_pdgid_branch]
    values = [vector_to_list(getattr(event, branch)) for branch in branches]
    statuses = (vector_to_list(getattr(event, args.gen_status_branch))
                if hasattr(event, args.gen_status_branch) else [])
    final_state, all_pions = [], []
    for i in range(min(map(len, values))):
        energy, eta, phi = map(finite_float, (values[0][i], values[1][i], values[2][i]))
        try:
            pdgid = int(values[3][i])
        except Exception:
            continue
        if energy is None or eta is None or phi is None or energy <= 0 or abs(pdgid) != 211:
            continue
        pion = {"energy": energy, "eta": eta, "phi": phi}
        all_pions.append(pion)
        if i < len(statuses):
            try:
                if int(statuses[i]) == 1:
                    final_state.append(pion)
            except Exception:
                pass
    candidates = final_state if len(final_state) >= 2 else all_pions
    return sorted(candidates, key=lambda p: p["energy"], reverse=True)[:2]


def selected_clusters(raw_energy, raw_eta, raw_phi, threshold=None):
    clusters = []
    n = min(len(raw_energy), len(raw_eta), len(raw_phi))
    selected_indices = set(range(len(raw_energy))) if threshold is None else set()
    for index in range(n):
        energy, eta, phi = map(finite_float, (raw_energy[index], raw_eta[index], raw_phi[index]))
        if energy is None or eta is None or phi is None or energy <= 0.0:
            continue
        if threshold is not None and energy < threshold:
            continue
        clusters.append({"index": index, "energy": energy, "eta": eta, "phi": phi})
        selected_indices.add(index)
    multiplicity = len(raw_energy) if threshold is None else len(clusters)
    return clusters, selected_indices, multiplicity


def sum_clustered_pfrh(energies, indices, selected_indices):
    total = 0.0
    for energy, index in zip(energies, indices):
        energy = finite_float(energy)
        try:
            index = int(index)
        except Exception:
            continue
        if index in selected_indices and energy is not None and energy > 0.0:
            total += energy
    return total


def clustered_rechit_percent(all_pfrh_energies, all_cluster_indices, selected_indices):
    n_total = min(len(all_pfrh_energies), len(all_cluster_indices))
    if n_total == 0:
        return None
    n_clustered = 0
    for index in all_cluster_indices[:n_total]:
        try:
            index = int(index)
        except Exception:
            continue
        if index >= 0 and index in selected_indices:
            n_clustered += 1
    return 100.0 * n_clustered / n_total


def matched_delta_rs(clusters, pions):
    return [match_dr for _, match_dr in matched_cluster_pairs(clusters, pions)]


def matched_cluster_pairs(clusters, pions):
    """Return the leading two clusters paired with their generator-match DR.

    Cluster order is always descending energy.  For two clusters and two
    generator pions, choose the direct or crossed assignment with the smaller
    total DR.  A one-cluster event contributes only its leading cluster.
    """
    clusters = sorted(clusters, key=lambda c: c["energy"], reverse=True)[:2]
    if not clusters or not pions:
        return []
    if len(clusters) == 1 or len(pions) == 1:
        c = clusters[0]
        match_dr = min(delta_r(c["eta"], c["phi"], p["eta"], p["phi"])
                       for p in pions)
        return [(c, match_dr)]
    c1, c2, p1, p2 = clusters[0], clusters[1], pions[0], pions[1]
    direct = [delta_r(c1["eta"], c1["phi"], p1["eta"], p1["phi"]),
              delta_r(c2["eta"], c2["phi"], p2["eta"], p2["phi"])]
    crossed = [delta_r(c1["eta"], c1["phi"], p2["eta"], p2["phi"]),
               delta_r(c2["eta"], c2["phi"], p1["eta"], p1["phi"])]
    chosen = direct if sum(direct) <= sum(crossed) else crossed
    return [(c1, chosen[0]), (c2, chosen[1])]


def empty_stats():
    return {
        "n_events": 0, "sum_ncl": 0.0, "sum_neff": 0.0,
        "sum_top2_over_pion_pair": 0.0, "sum_top2_over_all_pfrh": 0.0,
        "n_ge2": 0, "n_eq2": 0,
        "n_ge2_top2_clustered_pfrh_gt80": 0,
        "n_ge2_top2_all_pfrh_gt80": 0,
        "sum_clustered_rechit_pct": 0.0, "n_rechit_fraction_events": 0,
        "n_missing_pion_pair": 0, "n_zero_all_pfrh": 0,
        "n_zero_clustered_pfrh": 0,
    }


def update_stats(s, clusters, n_cl, pion_pair_energy, all_pfrh_energy,
                 clustered_pfrh_energy, rechit_pct):
    energies = [c["energy"] for c in clusters]
    top2 = sum(sorted(energies, reverse=True)[:2])
    s["n_events"] += 1
    s["sum_ncl"] += n_cl
    s["sum_neff"] += compute_neff(energies)
    s["sum_top2_over_pion_pair"] += top2 / pion_pair_energy if pion_pair_energy > 0 else 0
    s["sum_top2_over_all_pfrh"] += top2 / all_pfrh_energy if all_pfrh_energy > 0 else 0
    if rechit_pct is not None:
        s["sum_clustered_rechit_pct"] += rechit_pct
        s["n_rechit_fraction_events"] += 1
    s["n_missing_pion_pair"] += pion_pair_energy <= 0
    s["n_zero_all_pfrh"] += all_pfrh_energy <= 0
    s["n_zero_clustered_pfrh"] += clustered_pfrh_energy <= 0
    s["n_ge2"] += n_cl >= 2
    s["n_eq2"] += n_cl == 2
    if n_cl >= 2 and clustered_pfrh_energy > 0 and top2 / clustered_pfrh_energy > 0.8:
        s["n_ge2_top2_clustered_pfrh_gt80"] += 1
    if n_cl >= 2 and all_pfrh_energy > 0 and top2 / all_pfrh_energy > 0.8:
        s["n_ge2_top2_all_pfrh_gt80"] += 1


def make_hist(name, title, bins, xmin, xmax):
    hist = ROOT.TH1D(name, title, bins, xmin, xmax)
    hist.Sumw2()
    hist.SetDirectory(0)
    return hist


def book_distribution_hists(args, pf_names, dr_values, dt_values):
    groups = {f"{metric}_{mode}": {} for metric in
              ["ncl", "neff", "rechit_pct", "reco_energy", "gen_match_dr",
               "gen_match_dr_leading", "gen_match_dr_subleading",
               "gen_match_dr_ncl1", "gen_match_dr_ncl2plus"]
              for mode in ["nocut", "cut"]}
    # Seed-position DeltaR is requested for the no-cut selection only.
    groups["gen_match_dr_seed_nocut"] = {}
    ncl_max = int(round(args.ncl_max))
    for pf in pf_names:
        for dr in dr_values:
            for dt in dt_values:
                key = (pf, key_float(dr), key_float(dt))
                tag = sanitize_name(f"{pf}_DR{dr}_DT{dt}")
                for mode in ["nocut", "cut"]:
                    suffix = f"{tag}_{mode}"
                    groups[f"ncl_{mode}"][key] = make_hist(
                        f"h_ncl_{suffix}", "", ncl_max + 1, -0.5, ncl_max + 0.5)
                    for ibin in range(1, ncl_max + 2):
                        groups[f"ncl_{mode}"][key].GetXaxis().SetBinLabel(
                            ibin, str(ibin - 1)
                        )
                    groups[f"neff_{mode}"][key] = make_hist(
                        f"h_neff_{suffix}", "", args.neff_bins, 0, args.neff_max)
                    groups[f"rechit_pct_{mode}"][key] = make_hist(
                        f"h_rechit_pct_{suffix}", "", args.rechit_percent_bins, 0, 100.0001)
                    groups[f"reco_energy_{mode}"][key] = make_hist(
                        f"h_reco_energy_{suffix}", "", args.reco_energy_bins, 0,
                        args.reco_energy_max)
                    groups[f"gen_match_dr_{mode}"][key] = make_hist(
                        f"h_gen_match_dr_{suffix}", "", args.match_dr_bins, 0,
                        args.match_dr_max)
                    groups[f"gen_match_dr_leading_{mode}"][key] = make_hist(
                        f"h_gen_match_dr_leading_{suffix}", "", args.match_dr_bins,
                        0, args.match_dr_max)
                    groups[f"gen_match_dr_subleading_{mode}"][key] = make_hist(
                        f"h_gen_match_dr_subleading_{suffix}", "", args.match_dr_bins,
                        0, args.match_dr_max)
                    # Split by reconstructed cluster multiplicity.
                    for metric in ["gen_match_dr_ncl1", "gen_match_dr_ncl2plus"]:
                        groups[f"{metric}_{mode}"][key] = make_hist(
                            f"h_{metric}_{suffix}", "", args.match_dr_bins,
                            0, args.match_dr_max)
                groups["gen_match_dr_seed_nocut"][key] = make_hist(
                    f"h_gen_match_dr_seed_{tag}_nocut", "", args.match_dr_bins,
                    0, args.seed_match_dr_max)
    return groups


def book_scatter_points(pf_names):
    points = {mode: {} for mode in ["nocut", "cut"]}
    for mode in points:
        for pf in pf_names:
            for dr in SCATTER_DR_VALUES:
                for dt in SCATTER_DT_VALUES:
                    key = (pf, key_float(dr), key_float(dt))
                    points[mode][key] = {
                        "leading": [], "subleading": [],
                        # Same points, split by reconstructed multiplicity.
                        "ncl1": [],
                        "ncl_ge2_leading": [], "ncl_ge2_subleading": [],
                    }
    return points


def leading_cluster_times(raw_energy, raw_time):
    """Rank by energy first; invalid leading times are not backfilled."""
    ranked = []
    for index, value in enumerate(raw_energy):
        energy = finite_float(value)
        if energy is not None:
            ranked.append((index, energy))
    times = []
    for index, _ in sorted(ranked, key=lambda item: item[1], reverse=True)[:2]:
        if index >= len(raw_time):
            continue
        time = finite_float(raw_time[index])
        if time is not None and time != -999.0:
            times.append(time)
    return times


def collect_timing_values(tree, path, key, args, timing_values):
    """Read only energy/time branches, independent of other plot requirements."""
    if timing_values is None or key not in timing_values:
        return
    missing = [name for name in (args.energy_branch, args.cluster_time_branch)
               if not branch_exists(tree, name)]
    if missing:
        print(f"WARNING: no timing histogram from {path}; missing branches: "
              + ", ".join(missing))
        return
    for event in tree:
        timing_values[key].extend(leading_cluster_times(
            vector_to_list(getattr(event, args.energy_branch)),
            vector_to_list(getattr(event, args.cluster_time_branch)),
        ))


def book_timing_hists(timing_values, bin_width):
    """Use identical bins across the four pages, including valid negative times."""
    if not math.isfinite(bin_width) or bin_width <= 0:
        raise ValueError("--cluster-time-bin-width must be finite and positive.")
    populated = [values for values in timing_values.values() if values]
    # Center bins on multiples of the width. With the 1 ns default,
    # integer times lie at bin centers and edges are half-integers.
    if populated:
        first_center = math.floor(min(min(values) for values in populated) / bin_width + 0.5)
        last_center = math.floor(max(max(values) for values in populated) / bin_width + 0.5)
    else:
        first_center, last_center = 0, int(math.ceil(10.0 / bin_width))
    nbins = last_center - first_center + 1
    low = (first_center - 0.5) * bin_width
    high = low + nbins * bin_width
    hists = {}
    for (pf, dr, dt), values in timing_values.items():
        hist = ROOT.TH1D(sanitize_name(f"h_cluster_time_{pf}_DR{dr}_DT{dt}_nocut"),
                         f"{pf}, DR={dr:g}, DT={dt:g} ns;Cluster time [ns];Clusters",
                         nbins, low, high)
        hist.SetDirectory(0); hist.SetStats(0)
        for time in values:
            hist.Fill(time)
        hists[(pf, dr, dt)] = hist
    return hists


def process_one_ntuple(path, pf_name, dr, dt, args, hists, stats_by_mode,
                       scatter_points, eta_scatter_points, seed_scatter_points=None,
                       timing_values=None):
    source = ROOT.TFile.Open(path)
    if not source or source.IsZombie():
        print(f"WARNING: could not open {path}")
        return
    tree = source.Get(args.tree_name)
    if not tree:
        print(f"WARNING: missing tree {args.tree_name} in {path}")
        source.Close()
        return
    # A separate pass lets timing use events even when branches needed by
    # the existing matching/energy plots are absent. Iterating a TTree again
    # starts from entry zero and leaves the original calculations unchanged.
    collect_timing_values(tree, path, (pf_name, key_float(dr), key_float(dt)),
                          args, timing_values)
    required = [
        args.energy_branch, args.cluster_eta_branch, args.cluster_phi_branch,
        args.all_pfrh_energy_branch, args.all_pfrh_cluster_index_branch,
        args.clustered_pfrh_energy_branch, args.clustered_pfrh_cluster_index_branch,
        args.gen_energy_branch, args.gen_eta_branch, args.gen_phi_branch,
        args.gen_pdgid_branch,
    ]
    missing = [name for name in required if not branch_exists(tree, name)]
    if missing:
        print(f"WARNING: skipping {path}; missing branches: {', '.join(missing)}")
        source.Close()
        return
    # The seed-position branches are optional: without them only the
    # seed-position DeltaR page loses this file.
    has_seed_position = (branch_exists(tree, args.cluster_seed_eta_branch)
                         and branch_exists(tree, args.cluster_seed_phi_branch))
    if not has_seed_position:
        print(f"WARNING: {args.cluster_seed_eta_branch}/"
              f"{args.cluster_seed_phi_branch} missing in "
              f"{os.path.basename(path)}; no seed-position DeltaR from this file")
    has_seed_eta = branch_exists(tree, args.cluster_seed_eta_branch)
    if pf_name == "standardPF" and not has_seed_eta:
        print(f"WARNING: {args.cluster_seed_eta_branch} missing in "
              f"{os.path.basename(path)}; no cluster-eta versus seed-eta points")
    key = (pf_name, key_float(dr), key_float(dt))
    for mode in stats_by_mode:
        stats_by_mode[mode].setdefault(key, empty_stats())
    print(f"Processing {pf_name:12s} DR={dr:<4} DT={dt:<4} {tree.GetEntries():7d} events")
    for event in tree:
        raw_e = vector_to_list(getattr(event, args.energy_branch))
        raw_eta = vector_to_list(getattr(event, args.cluster_eta_branch))
        raw_phi = vector_to_list(getattr(event, args.cluster_phi_branch))
        raw_seed_eta = (vector_to_list(getattr(event, args.cluster_seed_eta_branch))
                        if pf_name == "standardPF" and has_seed_eta else [])
        all_rh_e = vector_to_list(getattr(event, args.all_pfrh_energy_branch))
        all_rh_idx = vector_to_list(getattr(event, args.all_pfrh_cluster_index_branch))
        incl_rh_e = vector_to_list(getattr(event, args.clustered_pfrh_energy_branch))
        incl_rh_idx = vector_to_list(getattr(event, args.clustered_pfrh_cluster_index_branch))
        pions = generator_pions(event, args)
        pion_pair_energy = sum(p["energy"] for p in pions) if len(pions) >= 2 else 0.0
        threshold = event_cut_threshold(pions, args)
        selections = {
            "nocut": selected_clusters(raw_e, raw_eta, raw_phi),
            "cut": selected_clusters(raw_e, raw_eta, raw_phi, threshold),
        }
        all_pfrh_energy = sum(clean_positive_energies(all_rh_e))
        for mode, (clusters, selected_indices, n_cl) in selections.items():
            clustered_energy = sum_clustered_pfrh(incl_rh_e, incl_rh_idx, selected_indices)
            rechit_pct = clustered_rechit_percent(all_rh_e, all_rh_idx, selected_indices)
            update_stats(stats_by_mode[mode][key], clusters, n_cl, pion_pair_energy,
                         all_pfrh_energy, clustered_energy, rechit_pct)
            energies = [c["energy"] for c in clusters]
            hists[f"ncl_{mode}"][key].Fill(n_cl)
            hists[f"neff_{mode}"][key].Fill(compute_neff(energies))
            hists[f"reco_energy_{mode}"][key].Fill(sum(energies))
            if rechit_pct is not None:
                hists[f"rechit_pct_{mode}"][key].Fill(rechit_pct)
            # Associate each cluster with its own seed using its original
            # branch index, never by separately sorting the seed positions.
            # Eta correlations do not require a generator-pion match.
            if key in eta_scatter_points[mode]:
                for rank, cluster in enumerate(sorted(
                        clusters, key=lambda c: c["energy"], reverse=True)[:2]):
                    index = cluster["index"]
                    if index >= len(raw_seed_eta):
                        continue
                    seed_eta = finite_float(raw_seed_eta[index])
                    if seed_eta is not None:
                        rank_name = "leading" if rank == 0 else "subleading"
                        eta_scatter_points[mode][key][rank_name].append(
                            (seed_eta, cluster["eta"]))
            matches = matched_cluster_pairs(clusters, pions)
            for rank, (cluster, match_dr) in enumerate(matches):
                hists[f"gen_match_dr_{mode}"][key].Fill(match_dr)
                rank_name = "leading" if rank == 0 else "subleading"
                hists[f"gen_match_dr_{rank_name}_{mode}"][key].Fill(match_dr)
                # matched_cluster_pairs already caps at the leading two
                # clusters, so the >=2 histogram gets exactly those two.
                if n_cl == 1:
                    hists[f"gen_match_dr_ncl1_{mode}"][key].Fill(match_dr)
                elif n_cl >= 2:
                    hists[f"gen_match_dr_ncl2plus_{mode}"][key].Fill(match_dr)
                if key in scatter_points[mode]:
                    point = (match_dr, cluster["energy"])
                    scatter_points[mode][key][rank_name].append(point)
                    if n_cl == 1:
                        scatter_points[mode][key]["ncl1"].append(point)
                    elif n_cl >= 2:
                        scatter_points[mode][key][
                            f"ncl_ge2_{rank_name}"].append(point)
        # Same no-cut matching, but with each cluster's position replaced by
        # its seed position.  Reusing selected_clusters keeps the energy
        # selection and the descending-energy ordering identical.
        if has_seed_position:
            seed_clusters, _, _ = selected_clusters(
                raw_e,
                vector_to_list(getattr(event, args.cluster_seed_eta_branch)),
                vector_to_list(getattr(event, args.cluster_seed_phi_branch)),
            )
            for _, match_dr in matched_cluster_pairs(seed_clusters, pions):
                hists["gen_match_dr_seed_nocut"][key].Fill(match_dr)
        # Additional no-cut seed scatter: retain the original cluster ranks
        # and multiplicity, replacing only the eta/phi used for matching.
        if has_seed_position and seed_scatter_points is not None and key in seed_scatter_points:
            clusters, _, n_cl = selections["nocut"]
            leading_two = sorted(clusters, key=lambda c: c["energy"], reverse=True)[:2]
            ranks_by_index = {c["index"]: rank for rank, c in enumerate(leading_two)}
            seed_eta = vector_to_list(getattr(event, args.cluster_seed_eta_branch))
            seed_phi = vector_to_list(getattr(event, args.cluster_seed_phi_branch))
            seed_matches = []
            for cluster in leading_two:
                index = cluster["index"]
                if index >= len(seed_eta) or index >= len(seed_phi):
                    continue
                eta, phi = finite_float(seed_eta[index]), finite_float(seed_phi[index])
                if eta is not None and phi is not None:
                    seed_matches.append(dict(cluster, eta=eta, phi=phi))
            for cluster, match_dr in matched_cluster_pairs(seed_matches, pions):
                rank_name = "leading" if ranks_by_index[cluster["index"]] == 0 else "subleading"
                point = (match_dr, cluster["energy"])
                seed_scatter_points[key][rank_name].append(point)
                if n_cl == 1:
                    seed_scatter_points[key]["ncl1"].append(point)
                elif n_cl >= 2:
                    seed_scatter_points[key][f"ncl_ge2_{rank_name}"].append(point)
    source.Close()


# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------

def make_scan_hist(name, title, dr_values, dt_values):
    hist = ROOT.TH2D(name, title, len(dr_values), 0, len(dr_values),
                     len(dt_values), 0, len(dt_values))
    hist.SetDirectory(0)
    for i, value in enumerate(dr_values, 1):
        hist.GetXaxis().SetBinLabel(i, str(value))
    for i, value in enumerate(dt_values, 1):
        hist.GetYaxis().SetBinLabel(i, str(value))
    return hist


HEATMAP_METRICS = [
    ("mean_ncl", "Average Number of HCAL Clusters per Event", "Average # clusters/event",
     HEATMAP_VALUE_FORMAT, (1.0, 3.0), (0.0, 3.0)),
    ("mean_neff", "Average Effective Number of Energy-Carrying HCAL Clusters",
     "#LTN_{eff}#GT = (#Sigma E)^{2}/#Sigma E^{2}", HEATMAP_VALUE_FORMAT,
     (1.0, 3.0), (0.0, 3.0)),
    ("mean_top2_over_pion_pair", "Average Leading-Two Cluster Energy / Generator-Pion Energy",
     "#LT(E_{1}+E_{2})/(E_{#pi1}+E_{#pi2})#GT", HEATMAP_VALUE_FORMAT,
     (0.0, 0.8), (0.0, 0.8)),
    ("mean_top2_over_all_pfrh", "Average Leading-Two Cluster Energy / All HCAL PFRecHit Energy",
     "#LT(E_{1}+E_{2})/E_{all PFRecHits}#GT", HEATMAP_VALUE_FORMAT,
     (0.0, 1.2), (0.0, 1.2)),
    ("mean_clustered_rechit_pct", "Average Percentage of HCAL PFRecHits in Selected Clusters",
     "Average % of PFRecHits clustered", HEATMAP_PERCENT_FORMAT,
     (0.0, 100.0), (0.0, 100.0)),
    ("frac_ge2", "Percentage of Events with #geq 2 HCAL Clusters",
     "% events with #geq 2 clusters", HEATMAP_PERCENT_FORMAT, (0, 100), (0, 100)),
    ("frac_eq2", "Percentage of Events with Exactly 2 HCAL Clusters",
     "% events with exactly 2 clusters", HEATMAP_PERCENT_FORMAT, (0, 100), (0, 100)),
    ("frac_ge2_top2_clustered_pfrh_gt80",
     "Percentage with #geq 2 Clusters and (E_{1}+E_{2})/E_{clustered PFRecHits}>0.8",
     "% all events passing both conditions", HEATMAP_PERCENT_FORMAT, (0, 100), (0, 100)),
    ("frac_ge2_top2_all_pfrh_gt80",
     "Percentage with #geq 2 Clusters and (E_{1}+E_{2})/E_{all PFRecHits}>0.8",
     "% all events passing both conditions", HEATMAP_PERCENT_FORMAT, (0, 100), (0, 100)),
]


def book_heatmaps(pf_names, dr_values, dt_values, mode):
    out = {}
    for pf in pf_names:
        out[pf] = {}
        for metric, title, *_ in HEATMAP_METRICS:
            name = f"h2_{metric}_{sanitize_name(pf)}_{mode}"
            out[pf][metric] = make_scan_hist(name, title, dr_values, dt_values)
    return out


def fill_heatmaps(heatmaps, stats, pf_names, dr_values, dt_values):
    for pf in pf_names:
        for ix, dr in enumerate(dr_values, 1):
            for iy, dt in enumerate(dt_values, 1):
                s = stats.get((pf, key_float(dr), key_float(dt)))
                if not s or s["n_events"] <= 0:
                    continue
                n = float(s["n_events"])
                nrh = float(s["n_rechit_fraction_events"])
                values = {
                    "mean_ncl": s["sum_ncl"] / n,
                    "mean_neff": s["sum_neff"] / n,
                    "mean_top2_over_pion_pair": s["sum_top2_over_pion_pair"] / n,
                    "mean_top2_over_all_pfrh": s["sum_top2_over_all_pfrh"] / n,
                    "mean_clustered_rechit_pct": s["sum_clustered_rechit_pct"] / nrh if nrh else 0,
                    "frac_ge2": 100 * s["n_ge2"] / n,
                    "frac_eq2": 100 * s["n_eq2"] / n,
                    "frac_ge2_top2_clustered_pfrh_gt80": 100 * s["n_ge2_top2_clustered_pfrh_gt80"] / n,
                    "frac_ge2_top2_all_pfrh_gt80": 100 * s["n_ge2_top2_all_pfrh_gt80"] / n,
                }
                for metric, value in values.items():
                    heatmaps[pf][metric].SetBinContent(ix, iy, value)


def make_heatmap_ratio(seed, nominal, name):
    ratio = seed.Clone(name)
    ratio.Reset("ICES")
    ratio.SetDirectory(0)
    for ix in range(1, seed.GetNbinsX() + 1):
        for iy in range(1, seed.GetNbinsY() + 1):
            den = nominal.GetBinContent(ix, iy)
            ratio.SetBinContent(ix, iy, seed.GetBinContent(ix, iy) / den if den != 0 else 0)
    return ratio


def finite_hist2_values(hist, positive_only=False):
    values = []
    for ix in range(1, hist.GetNbinsX() + 1):
        for iy in range(1, hist.GetNbinsY() + 1):
            value = float(hist.GetBinContent(ix, iy))
            if math.isfinite(value) and (not positive_only or value > 0):
                values.append(value)
    return values


_ALGORITHM_PALETTE_EXEC_COMMAND = None
_RATIO_PALETTE_EXEC_COMMAND = None
_DECLARED_PALETTE_GLOBALS = set()

# Algorithm-panel gradient, identical to the single-pion plotting code so the
# two documents share one color scheme.
ALGORITHM_PALETTE_GRADIENT = (
    [0.00, 0.50, 1.00],
    [0.20, 1.00, 0.95],
    [0.55, 1.00, 0.55],
    [0.80, 1.00, 0.20],
)

# Ratio-panel gradient: magenta at low, pure white at one, green at high.
# ROOT's kRedBlue has a dark desaturated midpoint, which turns a map of
# near-unity ratios into flat olive sludge; a genuinely white center keeps
# small deviations readable, and magenta/green stays clear of the gradient the
# algorithm panels use.
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
    """Build both heatmap palettes once, before any page is drawn."""
    global _ALGORITHM_PALETTE_EXEC_COMMAND, _RATIO_PALETTE_EXEC_COMMAND
    _RATIO_PALETTE_EXEC_COMMAND = _register_gradient_palette(
        "gDipionRatioPalette", *RATIO_PALETTE_GRADIENT
    )
    _ALGORITHM_PALETTE_EXEC_COMMAND = _register_gradient_palette(
        "gDipionAlgorithmPalette", *ALGORITHM_PALETTE_GRADIENT
    )


def palette_exec(name, is_ratio):
    """Apply a pad-local palette whenever ROOT repaints the canvas."""
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


def build_ratio_heatmaps(heatmaps_by_mode, ratio_heatmaps_by_mode, pf_names):
    """Build every ratio map, each with its own range centered at one.

    Every ratio panel is scaled independently from its own values, so a metric
    whose ratios sit within a percent of unity still shows real structure even
    when some other metric has a near-zero denominator.  The cost is that two
    ratio panels are no longer comparable by color, so read the printed range
    and the cell labels before comparing pages.

    Returns {mode: {metric: (zmin, zmax)}}.
    """
    nominal_name = "standardPF" if "standardPF" in pf_names else pf_names[0]
    seed_name = ("seedTimingPF" if "seedTimingPF" in pf_names
                 else (pf_names[1] if len(pf_names) > 1 else None))
    ranges_by_mode = {mode: {} for mode in heatmaps_by_mode}
    if seed_name is None:
        return ranges_by_mode

    saturating = []
    for mode, heatmaps in heatmaps_by_mode.items():
        for metric, *_ in HEATMAP_METRICS:
            ratio = make_heatmap_ratio(
                heatmaps[seed_name][metric], heatmaps[nominal_name][metric],
                f"h2_ratio_{metric}_{mode}",
            )
            ratio_heatmaps_by_mode[mode][metric] = ratio
            values = finite_hist2_values(ratio, positive_only=True)
            # The 0.05 floor stops an all-identical panel from being scaled
            # down to pure numerical noise.
            raw_deviation = max([abs(value - 1.0) for value in values] + [0.05])
            # A ratio scale wider than [0, 2] can no longer be symmetric about
            # one while remaining physical.  Saturate rarer extremes.
            deviation = min(1.08 * raw_deviation, 1.0)
            if raw_deviation > 1.0:
                saturating.append(f"{mode}/{metric}")
            ranges_by_mode[mode][metric] = (1.0 - deviation, 1.0 + deviation)

    if saturating:
        print(
            "WARNING: ratio values outside [0, 2] saturate the color scale on "
            f"{len(saturating)} panel(s): " + ", ".join(saturating)
        )
    return ranges_by_mode


def make_canvas(name, width, height):
    """Create a batch-safe ROOT canvas with a nonzero drawable area."""
    width, height = int(width), int(height)
    canvas = ROOT.TCanvas(name, "", width, height)
    # In headless ROOT sessions the constructor's window size can collapse to
    # zero after subtracting GUI decorations. SetCanvasSize controls the actual
    # drawable area used by TPDF and prevents zero-width/zero-height CropBoxes.
    canvas.SetCanvasSize(width, height)
    canvas.Modified()
    canvas.Update()
    return canvas


def print_pdf_page(canvas, pdf_path):
    """Flush every pad before ROOT serializes the canvas into the PDF."""
    canvas.cd()
    canvas.Modified()
    canvas.Update()
    canvas.Print(pdf_path)


def style_heatmap(hist, title, z_title, zmin, zmax, text_format):
    hist.SetTitle(title)
    hist.SetMinimum(zmin)
    hist.SetMaximum(zmax)
    hist.GetXaxis().SetTitle("DR")
    hist.GetYaxis().SetTitle("DT [ns]")
    hist.GetZaxis().SetTitle(z_title)
    for axis in [hist.GetXaxis(), hist.GetYaxis(), hist.GetZaxis()]:
        axis.CenterTitle()
    hist.GetXaxis().SetTitleSize(0.048); hist.GetYaxis().SetTitleSize(0.048)
    hist.GetZaxis().SetTitleSize(0.043); hist.GetXaxis().SetLabelSize(0.043)
    hist.GetYaxis().SetLabelSize(0.043); hist.GetZaxis().SetLabelSize(0.036)
    hist.GetXaxis().SetTitleOffset(1.0); hist.GetYaxis().SetTitleOffset(1.05)
    hist.GetZaxis().SetTitleOffset(1.28); hist.SetMarkerSize(1.55)
    ROOT.gStyle.SetPaintTextFormat(text_format)


def make_heatmap_page(heatmaps, ratio_heatmaps, pf_names, metric, title, z_title,
                      text_format, zrange, ratio_range, pdf_path,
                      selection_label, mode):
    nominal_name = "standardPF" if "standardPF" in heatmaps else pf_names[0]
    seed_name = "seedTimingPF" if "seedTimingPF" in heatmaps else (pf_names[1] if len(pf_names) > 1 else None)
    panels = [(nominal_name, heatmaps[nominal_name][metric], zrange)]
    if seed_name:
        panels.append((seed_name, heatmaps[seed_name][metric], zrange))
        panels.append(("seedTimingPF / standardPF", ratio_heatmaps[metric],
                       ratio_range))
    canvas = make_canvas(f"c_heatmap_{metric}_{mode}", 2350, 820)
    header = ROOT.TPad(f"head_{metric}_{mode}", "", 0, 0.90, 1, 1)
    header.SetFillStyle(0); header.Draw(); header.cd()
    latex = ROOT.TLatex(); latex.SetNDC(); latex.SetTextAlign(22)
    latex.SetTextFont(62); latex.SetTextSize(0.30); latex.DrawLatex(0.5, 0.67, title)
    latex.SetTextFont(42); latex.SetTextSize(0.23); latex.DrawLatex(0.5, 0.22, selection_label)
    canvas.cd(); body = ROOT.TPad(f"body_{metric}_{mode}", "", 0, 0, 1, 0.90)
    body.SetFillStyle(0); body.Draw(); body.cd(); body.Divide(len(panels), 1, 0.002, 0.002)
    keep = [header, body, latex]
    for index, (label, hist, limits) in enumerate(panels, 1):
        pad = body.cd(index); pad.SetLeftMargin(0.115); pad.SetRightMargin(0.17)
        pad.SetBottomMargin(0.13); pad.SetTopMargin(0.10); pad.SetTicks(1, 1)
        is_ratio = " / " in label
        style_heatmap(hist, label, "Ratio" if is_ratio else z_title,
                      limits[0], limits[1],
                      HEATMAP_VALUE_FORMAT if is_ratio else text_format)
        # TExec keeps the shared algorithm gradient on the two algorithm maps
        # and the magenta-white-green diverging palette on the ratio map, even
        # when ROOT repaints the full multipad canvas during PDF output.
        hist.Draw("AXIS")
        palette = palette_exec(f"palette_{metric}_{mode}_{index}", is_ratio)
        palette.Draw()
        hist.Draw("COLZ TEXT SAME")
        keep.append(palette)
    print_pdf_page(canvas, pdf_path)
    print(f"Added heatmap: {title} ({selection_label})")


# -----------------------------------------------------------------------------
# 1D distribution pages
# -----------------------------------------------------------------------------

def style_for_pf(pf_name, index=0):
    if pf_name == "standardPF":
        return ROOT.kAzure + 1, 1, 1, True, 0.22
    if pf_name == "seedTimingPF":
        return ROOT.kBlack, 1, 1, False, 0.0
    colors = [ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 7]
    return colors[index % len(colors)], 1 + index, 1, False, 0.0


def apply_pf_style(hist, pf_name, index=0):
    color, line_style, width, fill, alpha = style_for_pf(pf_name, index)
    hist.SetLineColor(color); hist.SetLineStyle(line_style); hist.SetLineWidth(width)
    if fill:
        hist.SetFillColorAlpha(color, alpha)
    else:
        hist.SetFillStyle(0)


def normalized_clone(hist, suffix="draw", fold_outliers=False):
    """Clone and normalize to unit area.

    The normalization denominator spans underflow through overflow, but a pad
    only draws the visible bins.  With ``fold_outliers`` the out-of-range
    entries are first moved into the first and last visible bins, so every
    entry is both drawn and counted and the visible bars really do sum to one.
    A DeltaR of exactly 0 already lands in bin 1 rather than underflow, so this
    matters at the top end of the axis, not at zero.
    """
    clone = hist.Clone(f"{hist.GetName()}_{suffix}")
    clone.SetDirectory(0)
    if fold_outliers:
        nbins = clone.GetNbinsX()
        for source, target in ((0, 1), (nbins + 1, nbins)):
            extra = clone.GetBinContent(source)
            if extra:
                clone.SetBinContent(target, clone.GetBinContent(target) + extra)
                clone.SetBinError(target, math.hypot(clone.GetBinError(target),
                                                     clone.GetBinError(source)))
                clone.SetBinContent(source, 0.0)
                clone.SetBinError(source, 0.0)
    integral = clone.Integral(0, clone.GetNbinsX() + 1)
    if integral > 0:
        clone.Scale(1.0 / integral)
    return clone


def make_distribution_grid(group, pf_names, pdf_path, dr_values, dt_values,
                           x_title, page_title, tag, x_range=None, ideal_line=None,
                           normalize=True, y_title="fraction of entries"):
    ncols, nrows = len(dr_values), len(dt_values)
    canvas = make_canvas(f"c_grid_{tag}", 700 * ncols, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.955, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(0.5, 0.46, page_title); title.SetNDC(); title.SetTextAlign(22)
    title.SetTextFont(62); title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.955)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]
    for row, dt in enumerate(dt_values):
        for col, dr in enumerate(dr_values):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.03)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
            pad_tag = sanitize_name(f"{tag}_{dr}_{dt}")
            key_base = (key_float(dr), key_float(dt))
            draw = {}
            for ipf, pf in enumerate(pf_names):
                raw = group.get((pf, *key_base))
                if raw is None:
                    continue
                if normalize:
                    draw[pf] = normalized_clone(raw, pad_tag)
                else:
                    draw[pf] = raw.Clone(f"{raw.GetName()}_{pad_tag}")
                    draw[pf].SetDirectory(0)
                apply_pf_style(draw[pf], pf, ipf)
                keep.append(draw[pf])
            if not draw:
                missing = ROOT.TLatex(0.5, 0.5, "missing point"); missing.SetNDC(); missing.SetTextAlign(22); missing.Draw()
                keep.append(missing); continue
            ymax = 1.28 * max(h.GetMaximum() for h in draw.values())
            first = True
            for pf in pf_names:
                hist = draw.get(pf)
                if hist is None:
                    continue
                hist.SetTitle(f"DR = {dr}, DT = {dt} ns")
                hist.SetMinimum(0); hist.SetMaximum(ymax)
                if x_range: hist.GetXaxis().SetRangeUser(*x_range)
                hist.GetXaxis().SetTitle(x_title)
                hist.GetYaxis().SetTitle(y_title)
                hist.GetXaxis().SetTitleSize(0.048); hist.GetXaxis().SetLabelSize(0.042)
                hist.GetXaxis().SetTitleOffset(1.08)
                hist.GetYaxis().SetTitleSize(0.048); hist.GetYaxis().SetLabelSize(0.042)
                hist.GetYaxis().SetTitleOffset(1.18)
                hist.Draw("HIST" if first else "HIST SAME"); first = False
            if ideal_line is not None:
                line = ROOT.TLine(ideal_line, 0, ideal_line, ymax); line.SetLineStyle(2)
                line.SetLineColor(ROOT.kGray + 2); line.Draw(); keep.append(line)
            if row == 0 and col == 0:
                legend = ROOT.TLegend(0.56, 0.62, 0.95, 0.86); legend.SetBorderSize(0); legend.SetFillStyle(0)
                for pf in pf_names:
                    if pf in draw: legend.AddEntry(draw[pf], pf, "lf" if pf == "standardPF" else "l")
                legend.Draw(); keep.append(legend)
    print_pdf_page(canvas, pdf_path)
    print(f"Added distribution grid: {page_title}")


# -----------------------------------------------------------------------------
# Efficiency scans with Clopper-Pearson intervals
# -----------------------------------------------------------------------------

EFFICIENCY_METRICS = [
    ("n_ge2", "Fraction with #geq2 resolved HCAL clusters"),
    ("n_ge2_top2_clustered_pfrh_gt80", "#geq2 and (E1+E2)/Eclustered PFRecHits > 0.8"),
    ("n_ge2_top2_all_pfrh_gt80", "#geq2 and (E1+E2)/Eall PFRecHits > 0.8"),
]


def efficiency_point(stats, stat_key):
    if not stats or stats["n_events"] <= 0:
        return None
    n, k = int(stats["n_events"]), int(stats[stat_key])
    value = k / n
    low = ROOT.TEfficiency.ClopperPearson(n, k, CL_1SIGMA, False)
    high = ROOT.TEfficiency.ClopperPearson(n, k, CL_1SIGMA, True)
    return value, value - low, high - value


def efficiency_style(graph, pf, index):
    color = ROOT.kAzure + 1 if pf == "standardPF" else ROOT.kRed + 1
    graph.SetLineColor(color); graph.SetMarkerColor(color)
    graph.SetMarkerStyle(20 + index); graph.SetMarkerSize(1.0); graph.SetLineWidth(1)


def make_efficiency_page(stats, pf_names, pdf_path, scan_values, fixed_values,
                         scan_axis, selection_label, tag):
    ncols, nrows = len(EFFICIENCY_METRICS), len(fixed_values)
    canvas = make_canvas(f"c_eff_{tag}", 2050, 1120)
    title_pad = ROOT.TPad(f"title_eff_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(0.5, 0.48, f"Efficiency scans: {selection_label}")
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62); title.SetTextSize(0.40); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_eff_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.002, 0.002)
    keep = [title_pad, title, grid]
    for row, fixed in enumerate(fixed_values):
        for col, (stat_key, metric_title) in enumerate(EFFICIENCY_METRICS):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.03)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14)
            pad.SetTicks(1, 1); pad.SetGridy(True)
            mg = ROOT.TMultiGraph()
            if row == 0:
                # Top row: move the legend to the upper-right corner
                legend = ROOT.TLegend(0.53, 0.66, 0.94, 0.83)
            else:
                # Bottom row: keep the current position
                legend = ROOT.TLegend(0.53, 0.20, 0.94, 0.37)
            legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetTextSize(0.043)
            have_graph = False
            for ipf, pf in enumerate(pf_names):
                rows = []
                for scan in scan_values:
                    dr, dt = (fixed, scan) if scan_axis == "dt" else (scan, fixed)
                    point = efficiency_point(stats.get((pf, key_float(dr), key_float(dt))), stat_key)
                    if point: rows.append((float(scan), *point))
                if not rows: continue
                x = array("d", [r[0] for r in rows]); y = array("d", [r[1] for r in rows])
                zero = array("d", [0.0] * len(rows)); elo = array("d", [r[2] for r in rows]); ehi = array("d", [r[3] for r in rows])
                graph = ROOT.TGraphAsymmErrors(len(rows), x, y, zero, zero, elo, ehi)
                efficiency_style(graph, pf, ipf); mg.Add(graph, "LP"); legend.AddEntry(graph, pf, "lp")
                have_graph = True; keep.append(graph)
            fixed_text = f"DR = {fixed}" if scan_axis == "dt" else f"DT = {fixed} ns"
            x_title = "DT [ns]" if scan_axis == "dt" else "DR"
            mg.SetTitle(f"{metric_title}, {fixed_text};{x_title};Efficiency")
            if have_graph:
                mg.Draw("A"); mg.SetMinimum(0); mg.SetMaximum(1.08)
                mg.GetXaxis().SetTitleSize(0.048); mg.GetXaxis().SetLabelSize(0.042)
                mg.GetXaxis().SetTitleOffset(1.08)
                mg.GetYaxis().SetTitleSize(0.048); mg.GetYaxis().SetLabelSize(0.042)
                mg.GetYaxis().SetTitleOffset(1.18)
                legend.Draw(); keep.extend([mg, legend])
            else:
                missing = ROOT.TLatex(0.5, 0.5, "missing point")
                missing.SetNDC(); missing.SetTextAlign(22); missing.Draw(); keep.append(missing)
    print_pdf_page(canvas, pdf_path)
    print(f"Added efficiency page: {selection_label}, {scan_axis}")


# -----------------------------------------------------------------------------
# Generator-matching pages: fixed DR panels with DT overlays
# -----------------------------------------------------------------------------

def dt_style(index):
    """Alternate diagonal-hatched fills and thin solid-step DT curves."""
    fill_colors = [ROOT.kAzure + 1, ROOT.kOrange + 1, ROOT.kMagenta + 1]
    step_colors = [ROOT.kBlack, ROOT.kGreen + 2, ROOT.kRed + 1]
    if index % 2 == 0:
        slot = (index // 2) % len(fill_colors)
        # ROOT hatch patterns: use alternating diagonal hatch directions/shades
        # for DT = 0, 2, 4 ns.
        hatch_styles = [3004, 3005, 3006]
        return fill_colors[slot], 1, True, hatch_styles[slot]
    slot = (index // 2) % len(step_colors)
    # DT = 1, 3, 5 ns are all thin, solid step lines.
    return step_colors[slot], 1, False, 0


def make_gen_matching_page(group, pf_names, display_pf, mode, pdf_path,
                           dr_values, dt_values, selection_label):
    ncols, nrows = 3, int(math.ceil(len(dr_values) / 3.0))
    tag = sanitize_name(f"genmatch_{display_pf}_{mode}")
    canvas = make_canvas(f"c_{tag}", 2100, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(0.5, 0.48, f"Generator-pion / cluster #DeltaR: {display_pf}, {selection_label}")
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62); title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]
    for index, dr in enumerate(dr_values):
        pad = grid.cd(index + 1); panel_tag = sanitize_name(f"{tag}_{dr}")
        pad.SetLeftMargin(0.13); pad.SetRightMargin(0.03)
        pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
        display = []
        for idt, dt in enumerate(dt_values):
            key = (display_pf, key_float(dr), key_float(dt))
            raw = group.get(key)
            if raw:
                hist = normalized_clone(raw, panel_tag, fold_outliers=True)
                color, line_style, is_filled, fill_style = dt_style(idt)
                hist.SetLineColor(color); hist.SetLineStyle(line_style); hist.SetLineWidth(1)
                if is_filled:
                    hist.SetFillColorAlpha(color, 0.45)
                    hist.SetFillStyle(fill_style)
                else:
                    hist.SetFillStyle(0)
                display.append((dt, hist, is_filled)); keep.append(hist)
        if display:
            ymax = 1.28 * max(h.GetMaximum() for _, h, _ in display)
            # Draw the transparent fills first, then the step curves on top so
            # none of the thin outlines gets hidden by a later filled curve.
            draw_order = ([entry for entry in display if entry[2]] +
                          [entry for entry in display if not entry[2]])
            for ih, (dt, hist, _) in enumerate(draw_order):
                hist.SetTitle(f"DR = {dr}"); hist.SetMinimum(0); hist.SetMaximum(ymax)
                hist.GetXaxis().SetRangeUser(0.0, 0.3)
                hist.GetXaxis().SetTitle("#DeltaR(cluster, #pi)")
                hist.GetYaxis().SetTitle("fraction of matches")
                hist.GetXaxis().SetTitleSize(0.048); hist.GetXaxis().SetLabelSize(0.042)
                hist.GetXaxis().SetTitleOffset(1.08)
                hist.GetYaxis().SetTitleSize(0.048); hist.GetYaxis().SetLabelSize(0.042)
                hist.GetYaxis().SetTitleOffset(1.18)
                hist.Draw("HIST" if ih == 0 else "HIST SAME")
            legend = ROOT.TLegend(0.66, 0.55, 0.95, 0.88)
            legend.SetBorderSize(0); legend.SetFillStyle(0); legend.SetNColumns(2)
            legend.SetTextSize(0.036)
            for dt, hist, is_filled in display:
                legend.AddEntry(hist, f"DT={dt:g} ns", "f" if is_filled else "l")
            legend.Draw(); keep.append(legend)
        else:
            missing = ROOT.TLatex(0.5, 0.5, "missing point")
            missing.SetNDC(); missing.SetTextAlign(22); missing.Draw(); keep.append(missing)
    print_pdf_page(canvas, pdf_path)
    print(f"Added gen-matching page: {display_pf}, {selection_label}")


def make_pf_row_matching_page(group, pf_names, pdf_path, dr_values, dt_values,
                              page_title, tag_stub, match_dr_max,
                              x_title="#DeltaR(cluster, #pi)"):
    """Draw one match-DR page with PF algorithms in separate rows.

    Columns are DR points and each panel overlays the DT values, so a page is
    fully described by the histogram group it is handed plus its title.
    """
    ncols, nrows = len(dr_values), len(pf_names)
    tag = sanitize_name(tag_stub)
    canvas = make_canvas(f"c_{tag}", max(2100, 620 * ncols), 520 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(0.5, 0.48, page_title)
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]

    for row, pf in enumerate(pf_names):
        for col, dr in enumerate(dr_values):
            pad = grid.cd(row * ncols + col + 1)
            panel_tag = sanitize_name(f"{tag}_{pf}_{dr}")
            pad.SetLeftMargin(0.14); pad.SetRightMargin(0.03)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.15); pad.SetTicks(1, 1)
            display = []
            for idt, dt in enumerate(dt_values):
                raw = group.get((pf, key_float(dr), key_float(dt)))
                if raw is None:
                    continue
                hist = normalized_clone(raw, panel_tag, fold_outliers=True)
                color, line_style, is_filled, fill_style = dt_style(idt)
                hist.SetLineColor(color); hist.SetLineStyle(line_style)
                hist.SetLineWidth(1)
                if is_filled:
                    hist.SetFillColorAlpha(color, 0.45)
                    hist.SetFillStyle(fill_style)
                else:
                    hist.SetFillStyle(0)
                display.append((dt, hist, is_filled)); keep.append(hist)

            if not display:
                missing = ROOT.TLatex(0.5, 0.5, "missing point")
                missing.SetNDC(); missing.SetTextAlign(22); missing.Draw()
                keep.append(missing)
                continue

            ymax = 1.28 * max(hist.GetMaximum() for _, hist, _ in display)
            draw_order = ([entry for entry in display if entry[2]] +
                          [entry for entry in display if not entry[2]])
            for ih, (_, hist, _) in enumerate(draw_order):
                hist.SetTitle(f"{pf}, DR = {dr}")
                hist.SetMinimum(0); hist.SetMaximum(ymax)
                hist.GetXaxis().SetRangeUser(0.0, match_dr_max)
                hist.GetXaxis().SetTitle(x_title)
                hist.GetYaxis().SetTitle("fraction of matches")
                hist.GetXaxis().SetTitleSize(0.050); hist.GetXaxis().SetLabelSize(0.041)
                hist.GetXaxis().SetTitleOffset(1.08)
                hist.GetYaxis().SetTitleSize(0.050); hist.GetYaxis().SetLabelSize(0.041)
                hist.GetYaxis().SetTitleOffset(1.28)
                hist.Draw("HIST" if ih == 0 else "HIST SAME")
            if row == 0 and col == 0:
                legend = ROOT.TLegend(0.58, 0.54, 0.96, 0.88)
                legend.SetBorderSize(0); legend.SetFillStyle(0)
                legend.SetNColumns(2); legend.SetTextSize(0.034)
                for dt, hist, is_filled in display:
                    legend.AddEntry(hist, f"DT={dt:g} ns", "f" if is_filled else "l")
                legend.Draw(); keep.append(legend)

    print_pdf_page(canvas, pdf_path)
    print(f"Added match-DR page: {page_title}")


def make_ranked_gen_matching_page(group, pf_names, mode, rank_name, pdf_path,
                                  dr_values, dt_values, selection_label,
                                  match_dr_max):
    """Draw one cluster-rank page with PF algorithms in separate rows."""
    rank_title = ("Highest-energy cluster" if rank_name == "leading"
                  else "Second-highest-energy cluster")
    make_pf_row_matching_page(
        group, pf_names, pdf_path, dr_values, dt_values,
        f"{rank_title}: generator-pion / cluster #DeltaR, {selection_label}",
        f"genmatch_{rank_name}_{mode}", match_dr_max,
    )


def make_energy_vs_match_dr_page(points_by_key, pf_names, mode, fixed_dr,
                                 dt_values, pdf_path, selection_label,
                                 match_dr_max):
    """Draw a 2x3 E-cluster versus generator-match-DR scatter page."""
    ncols, nrows = len(dt_values), len(pf_names)
    tag = sanitize_name(f"energy_vs_match_dr_DR{fixed_dr}_{mode}")
    canvas = make_canvas(f"c_{tag}", 2100, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(
        0.5, 0.48,
        f"Leading-two cluster energy vs. generator-matching #DeltaR: DR = {fixed_dr}, {selection_label}",
    )
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]

    page_energies = []
    for pf in pf_names:
        for dt in dt_values:
            rank_points = points_by_key.get(
                (pf, key_float(fixed_dr), key_float(dt)), {}
            )
            for rank_name in ["leading", "subleading"]:
                page_energies.extend(y for _, y in rank_points.get(rank_name, []))
    energy_max = 1.08 * max(page_energies) if page_energies else 1.0

    rank_styles = {
        "leading": (ROOT.kAzure + 1, 20, "highest-energy cluster"),
        "subleading": (ROOT.kRed + 1, 24, "second-highest-energy cluster"),
    }
    for row, pf in enumerate(pf_names):
        for col, dt in enumerate(dt_values):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
            frame = ROOT.TH2D(
                f"frame_{tag}_{row}_{col}", f"{pf}, DT = {dt:g} ns",
                1, 0.0, match_dr_max, 1, 0.0, energy_max,
            )
            frame.SetDirectory(0); frame.SetStats(0)
            frame.GetXaxis().SetTitle("#DeltaR(cluster, #pi)")
            frame.GetYaxis().SetTitle("E_{cluster} [GeV]")
            frame.GetXaxis().SetTitleSize(0.048); frame.GetXaxis().SetLabelSize(0.042)
            frame.GetXaxis().SetTitleOffset(1.08)
            frame.GetYaxis().SetTitleSize(0.048); frame.GetYaxis().SetLabelSize(0.042)
            frame.GetYaxis().SetTitleOffset(1.18)
            frame.Draw(); keep.append(frame)

            rank_points = points_by_key.get(
                (pf, key_float(fixed_dr), key_float(dt)), {}
            )
            graphs = []
            for rank_name in ["leading", "subleading"]:
                values = rank_points.get(rank_name, [])
                if not values:
                    continue
                x = array("d", [value[0] for value in values])
                y = array("d", [value[1] for value in values])
                graph = ROOT.TGraph(len(values), x, y)
                color, marker, _ = rank_styles[rank_name]
                graph.SetMarkerColor(color); graph.SetMarkerStyle(marker)
                graph.SetMarkerSize(0.60)
                graph.Draw("P SAME")
                graphs.append((rank_name, graph)); keep.append(graph)
            if row == 0 and col == 0 and graphs:
                legend = ROOT.TLegend(0.54, 0.72, 0.96, 0.89)
                legend.SetBorderSize(0); legend.SetFillStyle(0)
                legend.SetTextSize(0.035)
                for rank_name, graph in graphs:
                    legend.AddEntry(graph, rank_styles[rank_name][2], "p")
                legend.Draw(); keep.append(legend)

    print_pdf_page(canvas, pdf_path)
    print(f"Added energy-vs-match-DR page: DR={fixed_dr}, {selection_label}")


def make_cluster_time_page(timing_hists, display_pf, fixed_dr, pdf_path):
    """Six DT panels (two rows, three columns), with raw cluster counts."""
    tag = sanitize_name(f"cluster_time_{display_pf}_DR{fixed_dr}_nocut")
    canvas = make_canvas(f"c_{tag}", 2100, 1120)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(
        0.5, 0.48,
        f"Leading-two cluster times: {display_pf}, DR = {fixed_dr:g}, no energy cut or generator matching",
    )
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.36); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(3, 2, 0.001, 0.001)
    keep = [title_pad, title, grid]
    panels = [timing_hists[(display_pf, key_float(fixed_dr), key_float(dt))]
              for dt in TIMING_DT_VALUES]
    ymax = max(1.0, 1.2 * max(hist.GetMaximum() for hist in panels))
    color = ROOT.kAzure + 1 if display_pf == "standardPF" else ROOT.kRed + 1
    for index, (dt, hist) in enumerate(zip(TIMING_DT_VALUES, panels), start=1):
        pad = grid.cd(index)
        pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04)
        pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
        hist.SetTitle(f"{display_pf}, DT = {dt:g} ns")
        hist.GetXaxis().SetTitle("Cluster time [ns]")
        hist.GetYaxis().SetTitle("Clusters")
        hist.GetXaxis().SetTitleSize(0.048); hist.GetXaxis().SetLabelSize(0.042)
        hist.GetXaxis().SetTitleOffset(1.08)
        hist.GetYaxis().SetTitleSize(0.048); hist.GetYaxis().SetLabelSize(0.042)
        hist.GetYaxis().SetTitleOffset(1.18)
        hist.SetLineColor(color); hist.SetLineWidth(2)
        hist.SetMinimum(0.0); hist.SetMaximum(ymax)
        hist.Draw("HIST"); keep.append(hist)
        label = ROOT.TLatex(0.94, 0.82, f"{int(hist.GetEntries())} clusters")
        label.SetNDC(); label.SetTextAlign(32); label.SetTextSize(0.036)
        label.Draw(); keep.append(label)
    print_pdf_page(canvas, pdf_path)
    print(f"Added cluster-time page: {display_pf}, DR={fixed_dr:g}")


def make_multiplicity_energy_vs_match_dr_page(points_by_key, display_pf, fixed_dr,
                                             dt_values, pdf_path,
                                             selection_label, match_dr_max):
    """Cluster energy versus match-DR for one algorithm, split by multiplicity.

    The top row is events with exactly one reconstructed cluster; the bottom
    row is events with two or more, showing the two highest-energy clusters.
    Columns are DT points.
    """
    rows = [
        ("Events with exactly 1 HCAL cluster",
         [("ncl1", ROOT.kAzure + 1, 20, "the single cluster")]),
        ("Events with #geq 2 HCAL clusters",
         [("ncl_ge2_leading", ROOT.kAzure + 1, 20, "highest-energy cluster"),
          ("ncl_ge2_subleading", ROOT.kRed + 1, 24,
           "second-highest-energy cluster")]),
    ]
    ncols, nrows = len(dt_values), len(rows)
    tag = sanitize_name(f"ncl_energy_vs_match_dr_{display_pf}_DR{fixed_dr}")
    canvas = make_canvas(f"c_{tag}", 2100, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(
        0.5, 0.48,
        f"Cluster energy vs. generator-matching #DeltaR by cluster "
        f"multiplicity: {display_pf}, DR = {fixed_dr}, {selection_label}",
    )
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]

    # One common energy axis across the whole page.
    page_energies = []
    for _, series in rows:
        for series_key, *_ in series:
            for dt in dt_values:
                ranks = points_by_key.get(
                    (display_pf, key_float(fixed_dr), key_float(dt)), {}
                )
                page_energies.extend(y for _, y in ranks.get(series_key, []))
    energy_max = 1.08 * max(page_energies) if page_energies else 1.0

    for row, (row_title, series) in enumerate(rows):
        for col, dt in enumerate(dt_values):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
            frame = ROOT.TH2D(
                f"frame_{tag}_{row}_{col}", f"{row_title}, DT = {dt:g} ns",
                1, 0.0, match_dr_max, 1, 0.0, energy_max,
            )
            frame.SetDirectory(0); frame.SetStats(0)
            frame.GetXaxis().SetTitle("#DeltaR(cluster, #pi)")
            frame.GetYaxis().SetTitle("E_{cluster} [GeV]")
            frame.GetXaxis().SetTitleSize(0.048); frame.GetXaxis().SetLabelSize(0.042)
            frame.GetXaxis().SetTitleOffset(1.08)
            frame.GetYaxis().SetTitleSize(0.048); frame.GetYaxis().SetLabelSize(0.042)
            frame.GetYaxis().SetTitleOffset(1.18)
            frame.Draw(); keep.append(frame)

            ranks = points_by_key.get(
                (display_pf, key_float(fixed_dr), key_float(dt)), {}
            )
            graphs = []
            for series_key, color, marker, label in series:
                values = ranks.get(series_key, [])
                if not values:
                    continue
                x = array("d", [value[0] for value in values])
                y = array("d", [value[1] for value in values])
                graph = ROOT.TGraph(len(values), x, y)
                graph.SetMarkerColor(color); graph.SetMarkerStyle(marker)
                graph.SetMarkerSize(0.60)
                graph.Draw("P SAME")
                graphs.append((label, graph)); keep.append(graph)
            if col == 0 and graphs:
                legend = ROOT.TLegend(0.54, 0.72, 0.96, 0.89)
                legend.SetBorderSize(0); legend.SetFillStyle(0)
                legend.SetTextSize(0.035)
                for label, graph in graphs:
                    legend.AddEntry(graph, label, "p")
                legend.Draw(); keep.append(legend)

    print_pdf_page(canvas, pdf_path)
    print(f"Added multiplicity energy-vs-match-DR page: {display_pf}, DR={fixed_dr}")


def make_cluster_eta_vs_seed_eta_page(points_by_key, mode, dr_values,
                                      dt_values, pdf_path, selection_label):
    """StandardPF leading-two clusters; DR columns and DT rows."""
    ncols, nrows = len(dr_values), len(dt_values)
    tag = f"cluster_eta_vs_seed_eta_standardPF_{mode}"
    canvas = make_canvas(f"c_{tag}", 700 * ncols, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.955, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(
        0.5, 0.46,
        f"Cluster #eta vs. seed #eta: standardPF, leading two clusters, {selection_label}",
    )
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.955)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]
    # Shared limits on both axes keep the diagonal meaningful in every panel.
    coordinates = [coordinate
                   for ranks in points_by_key.values()
                   for values in ranks.values()
                   for point in values for coordinate in point]
    if coordinates:
        lo, hi = min(coordinates), max(coordinates)
        padding = max(0.08 * (hi - lo), 0.02)
        lo, hi = lo - padding, hi + padding
    else:
        lo, hi = -1.0, 1.0
    styles = [
        ("leading", ROOT.kAzure + 1, 20, "highest-energy / only cluster"),
        ("subleading", ROOT.kRed + 1, 24, "second-highest-energy cluster"),
    ]
    for row, dt in enumerate(dt_values):
        for col, dr in enumerate(dr_values):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
            frame = ROOT.TH2D(
                f"frame_{tag}_{row}_{col}", f"DR = {dr:g}, DT = {dt:g} ns",
                1, lo, hi, 1, lo, hi,
            )
            frame.SetDirectory(0); frame.SetStats(0)
            frame.GetXaxis().SetTitle("#eta_{seed}")
            frame.GetYaxis().SetTitle("#eta_{cluster}")
            for axis in [frame.GetXaxis(), frame.GetYaxis()]:
                axis.SetTitleSize(0.048); axis.SetLabelSize(0.042)
            frame.GetXaxis().SetTitleOffset(1.08)
            frame.GetYaxis().SetTitleOffset(1.18)
            frame.Draw(); keep.append(frame)
            diagonal = ROOT.TLine(lo, lo, hi, hi)
            diagonal.SetLineColor(ROOT.kGray + 1); diagonal.SetLineStyle(2)
            diagonal.Draw(); keep.append(diagonal)
            ranks = points_by_key.get(
                ("standardPF", key_float(dr), key_float(dt)), {})
            graphs = []
            for rank_name, color, marker, label in styles:
                values = ranks.get(rank_name, [])
                graph = ROOT.TGraph(
                    len(values), array("d", [x for x, _ in values]),
                    array("d", [y for _, y in values]),
                ) if values else ROOT.TGraph()
                graph.SetMarkerColor(color); graph.SetMarkerStyle(marker)
                graph.SetMarkerSize(0.60)
                if values:
                    graph.Draw("P SAME")
                graphs.append((label, graph)); keep.append(graph)
            if row == 0 and col == 0:
                legend = ROOT.TLegend(0.40, 0.73, 0.96, 0.89)
                legend.SetBorderSize(0); legend.SetFillStyle(0)
                legend.SetTextSize(0.030)
                for label, graph in graphs:
                    legend.AddEntry(graph, label, "p")
                legend.Draw(); keep.append(legend)
            if not any(ranks.values()):
                missing = ROOT.TLatex(0.5, 0.5, "No valid cluster / seed pairs")
                missing.SetNDC(); missing.SetTextAlign(22)
                missing.SetTextSize(0.035); missing.Draw(); keep.append(missing)
    print_pdf_page(canvas, pdf_path)
    print(f"Added cluster-eta vs seed-eta page: standardPF, {selection_label}")


def make_merged_energy_vs_match_dr_page(points_by_key, display_pf, dr_values,
                                        dt_values, pdf_path,
                                        selection_label, match_dr_max, use_seed=False):
    """Overlay multiplicities, with DR in rows and DT in columns."""
    series = [
        ("ncl1", ROOT.kBlue + 2, 20, "Cluster (N = 1)"),
        ("ncl_ge2_leading", ROOT.kAzure + 1, 20,
         "Leading Cluster(N #geq 2)"),
        ("ncl_ge2_subleading", ROOT.kRed + 1, 24,
         "Subleading Cluster (N #geq 2)"),
    ]
    position_label = "cluster seed" if use_seed else "cluster"
    rows = list(dr_values)
    ncols, nrows = len(dt_values), len(rows)
    tag = sanitize_name(f"merged_ncl_energy_vs_match_dr_{display_pf}")
    if use_seed:
        tag += "_seed"
    canvas = make_canvas(f"c_{tag}", 2100, 560 * nrows)
    title_pad = ROOT.TPad(f"title_{tag}", "", 0, 0.95, 1, 1)
    title_pad.SetFillStyle(0); title_pad.Draw(); title_pad.cd()
    title = ROOT.TLatex(
        0.5, 0.48,
        f"Cluster energy vs. {'seed-position ' if use_seed else ''}generator-matching #DeltaR by cluster "
        f"multiplicity: {display_pf}, {selection_label}",
    )
    title.SetNDC(); title.SetTextAlign(22); title.SetTextFont(62)
    title.SetTextSize(0.39); title.Draw()
    canvas.cd(); grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.95)
    grid.Draw(); grid.cd(); grid.Divide(ncols, nrows, 0.001, 0.001)
    keep = [title_pad, title, grid]

    # One common energy axis across the whole page.
    page_energies = []
    for fixed_dr in rows:
        for series_key, *_ in series:
            for dt in dt_values:
                ranks = points_by_key.get(
                    (display_pf, key_float(fixed_dr), key_float(dt)), {}
                )
                page_energies.extend(y for _, y in ranks.get(series_key, []))
    energy_max = 1.08 * max(page_energies) if page_energies else 1.0

    for row, fixed_dr in enumerate(rows):
        for col, dt in enumerate(dt_values):
            pad = grid.cd(row * ncols + col + 1)
            pad.SetLeftMargin(0.13); pad.SetRightMargin(0.04)
            pad.SetTopMargin(0.12); pad.SetBottomMargin(0.14); pad.SetTicks(1, 1)
            frame = ROOT.TH2D(
                f"frame_{tag}_{row}_{col}", f"DR = {fixed_dr:g}, DT = {dt:g} ns",
                1, 0.0, match_dr_max, 1, 0.0, energy_max,
            )
            frame.SetDirectory(0); frame.SetStats(0)
            frame.GetXaxis().SetTitle(f"#DeltaR({position_label}, #pi)")
            frame.GetYaxis().SetTitle("E_{cluster} [GeV]")
            frame.GetXaxis().SetTitleSize(0.048); frame.GetXaxis().SetLabelSize(0.042)
            frame.GetXaxis().SetTitleOffset(1.08)
            frame.GetYaxis().SetTitleSize(0.048); frame.GetYaxis().SetLabelSize(0.042)
            frame.GetYaxis().SetTitleOffset(1.18)
            frame.Draw(); keep.append(frame)

            ranks = points_by_key.get(
                (display_pf, key_float(fixed_dr), key_float(dt)), {}
            )
            graphs = []
            for series_key, color, marker, label in series:
                values = ranks.get(series_key, [])
                
                x = array("d", [value[0] for value in values])
                y = array("d", [value[1] for value in values])
                graph = ROOT.TGraph(len(values), x, y) if values else ROOT.TGraph()
                graph.SetMarkerColor(color); graph.SetMarkerStyle(marker)
                graph.SetMarkerSize(0.60)
                if values:
                    graph.Draw("P SAME")
                graphs.append((label, graph)); keep.append(graph)
            if col == 0 and graphs:
                legend = ROOT.TLegend(0.60, 0.60, 0.94, 0.83)
                legend.SetFillColor(ROOT.kWhite); legend.SetFillStyle(1001)
                legend.SetBorderSize(1); legend.SetLineColor(ROOT.kBlack)
                legend.SetLineWidth(1)
                legend.SetTextFont(42); legend.SetTextSize(0.033)
                legend.SetTextColor(ROOT.kBlack)
                legend.SetMargin(0.13); legend.SetEntrySeparation(0.15)
                for label, graph in graphs:
                    legend.AddEntry(graph, label, "p")
                legend.Draw(); keep.append(legend)

    print_pdf_page(canvas, pdf_path)
    print(f"Added merged multiplicity energy-vs-match-DR page: {display_pf}, position={position_label}")


# -----------------------------------------------------------------------------
# Output and main
# -----------------------------------------------------------------------------

def print_summary(stats_by_mode, pf_names, dr_values, dt_values, args):
    for mode, stats in stats_by_mode.items():
        label = "no cut" if mode == "nocut" else cut_label(args)
        print(f"\nSummary: {label}")
        print("PF             DR   DT   events   <Ncl>  <Neff>  <top2/Epipi>  <top2/EallRH>  <RH clustered %>")
        for pf in pf_names:
            for dr in dr_values:
                for dt in dt_values:
                    s = stats.get((pf, key_float(dr), key_float(dt)))
                    if not s or not s["n_events"]: continue
                    n, nrh = float(s["n_events"]), float(s["n_rechit_fraction_events"])
                    rh = s["sum_clustered_rechit_pct"] / nrh if nrh else 0
                    print(f"{pf:<14} {dr:3.1f} {dt:4.1f} {int(n):8d} {s['sum_ncl']/n:7.3f} "
                          f"{s['sum_neff']/n:7.3f} {s['sum_top2_over_pion_pair']/n:13.3f} "
                          f"{s['sum_top2_over_all_pfrh']/n:14.3f} {rh:17.2f}")


def write_output_root(path, heatmaps_by_mode, ratio_heatmaps_by_mode, hists,
                      scatter_points, eta_scatter_points, seed_scatter_points=None):
    out = ROOT.TFile(path, "RECREATE")
    for mode, by_pf in heatmaps_by_mode.items():
        directory = out.mkdir(f"heatmaps_{mode}"); directory.cd()
        for pf, metrics in by_pf.items():
            pfdir = directory.mkdir(sanitize_name(pf)); pfdir.cd()
            for hist in metrics.values(): hist.Write()
            directory.cd()
        ratio_dir = directory.mkdir("seedTimingPF_over_standardPF"); ratio_dir.cd()
        for hist in ratio_heatmaps_by_mode[mode].values(): hist.Write()
        out.cd()
    distdir = out.mkdir("distributions")
    for group_name, group in hists.items():
        groupdir = distdir.mkdir(group_name); groupdir.cd()
        for hist in group.values(): hist.Write()
        out.cd()
    scatter_dir = out.mkdir("scatter_points")
    for mode, points_by_key in scatter_points.items():
        mode_dir = scatter_dir.mkdir(mode); mode_dir.cd()
        for (pf, dr, dt), ranks in points_by_key.items():
            for rank_name, values in ranks.items():
                if not values:
                    continue
                x = array("d", [value[0] for value in values])
                y = array("d", [value[1] for value in values])
                graph = ROOT.TGraph(len(values), x, y)
                graph.SetName(sanitize_name(
                    f"g_energy_vs_match_dr_{pf}_DR{dr}_DT{dt}_{rank_name}_{mode}"
                ))
                graph.SetTitle(
                    f"{pf}, DR={dr:g}, DT={dt:g} ns, {rank_name};"
                    "#DeltaR(cluster, #pi);E_{cluster} [GeV]"
                )
                graph.Write()
        out.cd()
    eta_dir = out.mkdir("cluster_eta_vs_seed_eta")
    for mode, points_by_key in eta_scatter_points.items():
        mode_dir = eta_dir.mkdir(mode); mode_dir.cd()
        for (pf, dr, dt), ranks in points_by_key.items():
            for rank_name, values in ranks.items():
                if not values:
                    continue
                graph = ROOT.TGraph(
                    len(values), array("d", [x for x, _ in values]),
                    array("d", [y for _, y in values]),
                )
                graph.SetName(sanitize_name(
                    f"g_cluster_eta_vs_seed_eta_{pf}_DR{dr}_DT{dt}_{rank_name}_{mode}"))
                graph.SetTitle(
                    f"{pf}, DR={dr:g}, DT={dt:g} ns, {rank_name};"
                    "#eta_{seed};#eta_{cluster}")
                graph.Write()
        out.cd()
    seed_dir = out.mkdir("seed_scatter_points_nocut"); seed_dir.cd()
    for (pf, dr, dt), ranks in (seed_scatter_points or {}).items():
        for rank_name, values in ranks.items():
            if not values:
                continue
            graph = ROOT.TGraph(
                len(values), array("d", [x for x, _ in values]),
                array("d", [y for _, y in values]),
            )
            graph.SetName(sanitize_name(
                f"g_energy_vs_seed_match_dr_{pf}_DR{dr}_DT{dt}_{rank_name}_nocut"))
            graph.SetTitle(
                f"{pf}, DR={dr:g}, DT={dt:g} ns, {rank_name};"
                "#DeltaR(cluster seed, #pi);E_{cluster} [GeV]")
            graph.Write()
    out.Close()


def main():
    args = parse_args()
    if not math.isfinite(args.cluster_time_bin_width) or args.cluster_time_bin_width <= 0:
        raise ValueError("--cluster-time-bin-width must be finite and positive.")
    dr_values = parse_float_list(args.dr_values, DEFAULT_DR_VALUES)
    dt_values = parse_float_list(args.dt_values, DEFAULT_DT_VALUES)
    dr_dist = parse_float_list(args.dr_dist_values, DEFAULT_DR_DIST_VALUES)
    dt_dist = parse_float_list(args.dt_dist_values, DEFAULT_DT_DIST_VALUES)
    pf_names = list(dict.fromkeys(args.pf_names))
    if len(pf_names) < 2:
        raise ValueError("Ratio plots require both standardPF and seedTimingPF inputs.")
    print("PF variants:", ", ".join(pf_names))
    print("Cluster cut:", cut_label(args))
    lookup = build_file_lookup(args.input_dir, args.prefix, pf_names, dr_values,
                               dt_values, args.debug_files)
    if not lookup:
        raise RuntimeError("No matching ntuples found; check input directory, prefix, PF names, DR and DT values.")
    for pf in pf_names:
        for dr in dr_values:
            for dt in dt_values:
                if (pf, key_float(dr), key_float(dt)) not in lookup:
                    print(f"Missing ntuple: pf={pf}, DR={dr}, DT={dt}")
    hists = book_distribution_hists(args, pf_names, dr_values, dt_values)
    scatter_points = book_scatter_points(pf_names)
    seed_scatter_points = book_scatter_points(["standardPF"])["nocut"]
    eta_scatter_points = {
        mode: {("standardPF", key_float(dr), key_float(dt)):
               {"leading": [], "subleading": []}
               for dr in dr_dist for dt in dt_dist}
        for mode in ["nocut", "cut"]
    }
    timing_values = {
        (pf, key_float(dr), key_float(dt)): []
        for pf in ("standardPF", "seedTimingPF") if pf in pf_names
        for dr in TIMING_DR_VALUES for dt in TIMING_DT_VALUES
    }
    stats_by_mode = {"nocut": {}, "cut": {}}
    for pf in pf_names:
        for dr in dr_values:
            for dt in dt_values:
                path = lookup.get((pf, key_float(dr), key_float(dt)))
                if path:
                    process_one_ntuple(path, pf, dr, dt, args, hists,
                                       stats_by_mode, scatter_points, eta_scatter_points,
                                       seed_scatter_points, timing_values)
    hists["cluster_time_leading_two_nocut"] = book_timing_hists(
        timing_values, args.cluster_time_bin_width)
    print_summary(stats_by_mode, pf_names, dr_values, dt_values, args)
    heatmaps_by_mode, ratio_heatmaps_by_mode = {}, {"nocut": {}, "cut": {}}
    for mode in ["nocut", "cut"]:
        heatmaps_by_mode[mode] = book_heatmaps(pf_names, dr_values, dt_values, mode)
        fill_heatmaps(heatmaps_by_mode[mode], stats_by_mode[mode], pf_names, dr_values, dt_values)
    ratio_ranges_by_mode = build_ratio_heatmaps(
        heatmaps_by_mode, ratio_heatmaps_by_mode, pf_names
    )
    install_heatmap_palettes()
    os.makedirs(args.output_dir, exist_ok=True)
    pf_tag = "_".join(map(sanitize_name, pf_names))
    pdf_path = os.path.join(args.output_dir, args.output_pdf or f"hcal_cluster_ntuple_scan_{pf_tag}_extended.pdf")
    root_path = os.path.join(args.output_dir, args.output_root or f"hcal_cluster_ntuple_scan_{pf_tag}_extended.root")
    # Never open a ROOT multipage PDF with a tiny canvas: in a headless
    # session its drawable area can collapse to zero and every page then gets
    # an invalid CropBox. A normal explicitly-sized canvas avoids that ROOT
    # failure mode.
    # Match the opener to the heatmap canvases so the first PDF page inherits
    # the same landscape page geometry as every following heatmap page.
    opener = make_canvas("pdf_open", 2350, 820)
    opener.Print(pdf_path + "[")

    # Consecutive no-cut/cut heatmaps for every quantity.
    for metric, title, z_title, text_format, nocut_range, cut_range in HEATMAP_METRICS:
        for mode, selection, zrange in [
            ("nocut", "No HCAL-cluster energy cut", nocut_range),
            ("cut", cut_label(args, root_text=True), cut_range),
        ]:
            make_heatmap_page(heatmaps_by_mode[mode], ratio_heatmaps_by_mode[mode],
                              pf_names, metric, title, z_title, text_format,
                              zrange, ratio_ranges_by_mode[mode][metric],
                              pdf_path, selection, mode)

    # Consecutive no-cut/cut 1D pages at the script's chosen DR/DT points.
    distribution_specs = [
        ("ncl", "N_{clusters}", "HCAL cluster multiplicity", (-0.5, int(round(args.ncl_max)) + 0.5), 2.0),
        ("neff", "N_{eff}", "Effective number of energy-carrying clusters", (0, args.neff_max), 2.0),
        ("rechit_pct", "% of HCAL PFRecHits clustered", "Percentage of Clustered HCAL PFRecHits", (0, 100.01), None),
        ("reco_energy", "#Sigma E_{cluster} [GeV]", "Reconstructed HCAL energy (sum of selected-cluster energies)", (0, args.reco_energy_max), None),
    ]
    for metric, x_title, title, x_range, ideal in distribution_specs:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            # The clustered-rechit percentage is an event-count distribution;
            # the other 1D overlays remain normalized to unit area.
            use_counts = metric == "rechit_pct"
            make_distribution_grid(hists[f"{metric}_{mode}"], pf_names, pdf_path,
                                   dr_dist, dt_dist, x_title, f"{title}: {selection}",
                                   f"{metric}_{mode}", x_range, ideal,
                                   normalize=not use_counts,
                                   y_title="# of events" if use_counts else "fraction of entries")

    # Consecutive no-cut/cut efficiency pages for each scan direction.
    for scan_axis, scan_values, fixed_values in [
        ("dt", dt_values, [0.1, 0.2]), ("dr", dr_values, [3.0, 5.0])
    ]:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            make_efficiency_page(stats_by_mode[mode], pf_names, pdf_path,
                                 scan_values, fixed_values, scan_axis, selection,
                                 f"{scan_axis}_{mode}")

    # Exactly four gen-matching pages: algorithm x selection.
    for display_pf in ["standardPF", "seedTimingPF"]:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            make_gen_matching_page(hists[f"gen_match_dr_{mode}"], pf_names,
                                   display_pf, mode, pdf_path, dr_values, dt_values,
                                   selection)

    # Exactly four additional ranked-cluster matching pages:
    # leading/subleading x no-cut/cut, with PF algorithms in separate rows.
    for rank_name in ["leading", "subleading"]:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            make_ranked_gen_matching_page(
                hists[f"gen_match_dr_{rank_name}_{mode}"], pf_names, mode,
                rank_name, pdf_path, dr_values, dt_values, selection,
                args.match_dr_max,
            )

    # Two pages split by reconstructed cluster multiplicity, no cut, with the
    # PF algorithms in separate rows.
    for group_name, page_title, tag_stub in [
        ("gen_match_dr_ncl1_nocut",
         "Events with exactly 1 HCAL cluster: generator-pion / cluster "
         "#DeltaR, No HCAL-cluster energy cut",
         "genmatch_ncl1_nocut"),
        ("gen_match_dr_ncl2plus_nocut",
         "Events with #geq 2 HCAL clusters, two highest-energy clusters: "
         "generator-pion / cluster #DeltaR, No HCAL-cluster energy cut",
         "genmatch_ncl2plus_nocut"),
    ]:
        make_pf_row_matching_page(
            hists[group_name], pf_names, pdf_path, dr_values, dt_values,
            page_title, tag_stub, args.match_dr_max,
        )

    # One page: DeltaR from the cluster seed position, no cut, all matches.
    make_pf_row_matching_page(
        hists["gen_match_dr_seed_nocut"], pf_names, pdf_path, dr_values, dt_values,
        "Seed-position generator-pion / cluster #DeltaR, "
        "No HCAL-cluster energy cut",
        "genmatch_seed_nocut", args.seed_match_dr_max,
        x_title="#DeltaR(cluster seed, #pi)",
    )

    # Exactly four E-versus-match-DR pages: DR=0.2/0.3 x no-cut/cut.
    for fixed_dr in SCATTER_DR_VALUES:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            make_energy_vs_match_dr_page(
                scatter_points[mode], pf_names, mode, fixed_dr,
                SCATTER_DT_VALUES, pdf_path, selection, args.match_dr_max,
            )

    # Four timing pages: one algorithm and one fixed DR per page; all six DTs.
    for fixed_dr in TIMING_DR_VALUES:
        for display_pf in ("standardPF", "seedTimingPF"):
            if display_pf in pf_names:
                make_cluster_time_page(hists["cluster_time_leading_two_nocut"],
                                       display_pf, fixed_dr, pdf_path)

    # standardPF only, DT = 3/4/5 ns, split by cluster multiplicity, no cut.
    scatter_pf = "standardPF" if "standardPF" in pf_names else pf_names[0]
    for fixed_dr in SCATTER_DR_VALUES:
        make_multiplicity_energy_vs_match_dr_page(
            scatter_points["nocut"], scatter_pf, fixed_dr, SCATTER_DT_VALUES,
            pdf_path, "No HCAL-cluster energy cut", args.match_dr_max,
        )

    # Additional standardPF pages; all existing pages above are retained.
    if "standardPF" in pf_names:
        for mode, selection in [("nocut", "No HCAL-cluster energy cut"),
                                ("cut", cut_label(args, root_text=True))]:
            make_cluster_eta_vs_seed_eta_page(
                eta_scatter_points[mode], mode, dr_dist, dt_dist,
                pdf_path, selection,
            )
        make_merged_energy_vs_match_dr_page(
            scatter_points["nocut"], "standardPF", SCATTER_DR_VALUES,
            SCATTER_DT_VALUES, pdf_path, "No HCAL-cluster energy cut",
            args.match_dr_max,
        )

    # Final page: the same merged layout using cluster-seed eta and phi.
    if "standardPF" in pf_names:
        make_merged_energy_vs_match_dr_page(
            seed_scatter_points, "standardPF", SCATTER_DR_VALUES,
            SCATTER_DT_VALUES, pdf_path, "No HCAL-cluster energy cut",
            args.seed_match_dr_max, use_seed=True,
        )

    opener.Modified()
    opener.Update()
    opener.Print(pdf_path + "]")
    write_output_root(root_path, heatmaps_by_mode, ratio_heatmaps_by_mode,
                      hists, scatter_points, eta_scatter_points, seed_scatter_points)
    print(f"\nSaved combined PDF: {pdf_path}\nSaved ROOT file: {root_path}")


if __name__ == "__main__":
    main()