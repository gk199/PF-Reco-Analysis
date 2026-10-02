#!/usr/bin/env python3
"""
Compare standard PF and timing PF HCAL clusters matched to LLP decay products, as a function of the expected
time delay of the decay product. Inputs are the label files written by label_llp_clusters.py.

Usage:
    python3 compare_timing_vs_standard.py \
        --standard labels_standardPF.root \
        --timing   labels_seedTimingPF.root \
        --output   timing_vs_standard_LLP.pdf
"""

import argparse
import awkward as ak
import numpy as np
import ROOT
import uproot

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser()
parser.add_argument("--standard", default= " ", help="Labels from the standard PF ntuple")
parser.add_argument("--timing",   default= " ", help="Labels from the timing PF ntuple")
parser.add_argument("--timingLabel", default="seed timing PF (4 ns)")
parser.add_argument("--output",   default="timing_vs_standard_LLP.pdf")
parser.add_argument("--clusterE", type=float, default=2.0, help="Minimum HB cluster energy for the cluster-level plots")
# Correctly matched cuts: three values each. Each one is scanned in its own row, with the other two at their middle value
parser.add_argument("--energyFrac", type=float, nargs=3, default=[0.4, 0.65, 0.8], help="sum E(matched clusters) / E_b above this")
parser.add_argument("--timeCut", type=float, nargs=3, default=[1.0, 0.5, 0.1], help="Mean distance [ns] to the cluster time = delay + offset line below this")
parser.add_argument("--dRCut", type=float, nargs=3, default=[0.3, 0.15, 0.05], help="Mean cluster dR below this")
parser.add_argument("--timeOffset", type=float, default=None, help="Offset [ns] of the cluster time vs delay line (default: median from standard PF)")
parser.add_argument("--hbDepthEdges", type=float, nargs=3, default=[190.2, 214.2, 244.8], metavar=("R12", "R23", "R34"),
                    help="Radii [cm] of the HB depth 1/2, 2/3 and 3/4 boundaries (defaults as in Run3-HCAL-LLP-Analysis)")
args = parser.parse_args()

LABELS = ["standard PF", args.timingLabel]
COLORS = [ROOT.kBlue + 1, ROOT.kRed + 1]
DELAY_BINS = np.array([0., 0.5, 1., 2., 3., 4., 5., 6., 8., 10., 15.])


def load(fname, tree, prefix):
    t = uproot.open(fname)[tree]
    ids = t.arrays(["run", "lumi", "event"])
    obj = t.arrays(filter_name=f"{prefix}_*")
    obj = ak.zip({k[len(prefix) + 1:]: obj[k] for k in obj.fields})
    return ids, obj


def align(ids_ref, ids, obj):
    """Reorders obj so its events follow ids_ref (the files are not in the same event order)."""
    key = lambda i: (int(i.run), int(i.lumi), int(i.event))
    position = {key(i): n for n, i in enumerate(ids)}
    missing = [key(i) for i in ids_ref if key(i) not in position]
    if missing:
        raise RuntimeError(f"{len(missing)} events missing from one of the inputs, e.g. {missing[:3]}")
    return obj[np.array([position[key(i)] for i in ids_ref])]


def style(h, color, ytitle=None):
    h.SetLineColor(color)
    h.SetMarkerColor(color)
    h.SetMarkerStyle(20)
    h.SetLineWidth(2)
    h.SetTitle("")
    if ytitle: h.GetYaxis().SetTitle(ytitle)
    return h


def legend(hists, labels, x1=0.55, y1=0.72, x2=0.88, y2=0.88):
    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.036)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "lp")
    return leg


def fill(h, *cols, weights=None):
    cols = [ak.to_numpy(c).astype(np.float64) for c in cols]
    w = np.ones(len(cols[0])) if weights is None else ak.to_numpy(weights).astype(np.float64)
    if len(cols[0]) == 0: return h
    if len(cols) == 1: h.FillN(len(w), cols[0], w)
    else:
        for x, y, wi in zip(cols[0], cols[1], w): h.Fill(x, y, wi)
    return h


def read_cone(fname):
    f = uproot.open(fname)
    if "config" not in f:
        raise RuntimeError(f"{fname} has no config tree; rerun label_llp_clusters.py to record the matching cone")
    return float(f["config"]["deltaR"].array(library="np")[0])


def match_text():
    """Labels the matching cone on plots that depend on the cluster to LLP (decay product) match."""
    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.032)
    text.SetTextAlign(31)
    text.DrawLatex(0.89, 0.92, f"matched: #DeltaR < {DELTA_R:g}")
    return text


def page(canvas, draw):
    canvas.Clear()
    canvas.cd()
    keep = draw()
    canvas.Print(args.output)
    return keep


def matched_clusters(decay, clus, e_min=None):
    """Mask (event, decay product, cluster) of the HB clusters matched to each decay product, optionally above e_min."""
    use = decay.iGen[:, :, np.newaxis] == clus.iGen[:, np.newaxis, :]
    return use if e_min is None else use & (clus.energy > e_min)[:, np.newaxis, :]


def matched_mean(use, clus, value, weighted, valid=None):
    """Mean of a cluster quantity over the matched clusters of each decay product (per event; -1 if none)."""
    if valid is not None: use = use & valid[:, np.newaxis, :]
    w = use * (clus.energy if weighted else ak.ones_like(clus.energy))[:, np.newaxis, :]
    sum_w = ak.sum(w, axis=-1)
    return ak.where(sum_w > 0, ak.sum(w * value[:, np.newaxis, :], axis=-1) / ak.where(sum_w > 0, sum_w, 1.), -1.)


# ----- Load and align -----
ids_std, decay_std = load(args.standard, "llpDecays", "decay")
ids_tim, decay_tim = load(args.timing, "llpDecays", "decay")
_, clus_std = load(args.standard, "clusterLabels", "clus")
_, clus_tim = load(args.timing, "clusterLabels", "clus")
decay_tim = align(ids_std, ids_tim, decay_tim)
clus_tim  = align(ids_std, ids_tim, clus_tim)
if not ak.all(decay_std.iGen == decay_tim.iGen):
    raise RuntimeError("Gen labels differ between inputs after event alignment")
DELTA_R = read_cone(args.standard)
if read_cone(args.timing) != DELTA_R:
    raise RuntimeError("Inputs were labeled with different matching cones")

decays   = [ak.flatten(d) for d in (decay_std, decay_tim)]
clusters = [ak.flatten(c) for c in (clus_std, clus_tim)]
decays   = [d[d.inHBAcceptance] for d in decays]
hb_clus  = [c[c.isHB & (c.energy > args.clusterE)] for c in clusters]
hb_clus_evt = [c[c.isHB & (c.energy > args.clusterE)] for c in (clus_std, clus_tim)]  # per event, for multiplicities
print(f"{len(ids_std)} events, {len(decays[0])} LLP decay products in HB acceptance, matching cone dR < {DELTA_R:g}")

# Offset of the cluster time = delay + offset line (slope 1), from the standard PF clusters of the time vs delay page
TIME_OFFSET = args.timeOffset
if TIME_OFFSET is None:
    m = hb_clus[0][(hb_clus[0].iGen >= 0) & (hb_clus[0].time > -100)]
    TIME_OFFSET = float(np.median(ak.to_numpy(m.time - m.expDelay))) if len(m) > 0 else 0.
print(f"Cluster time vs delay line: time = delay {TIME_OFFSET:+.2f} ns (used for both PF versions)")

canvas = ROOT.TCanvas("c", "c", 800, 600)
canvas.SetLeftMargin(0.13)
canvas.Print(args.output + "[")

# ----- Truth: expected delay and its components -----
def draw_truth():
    d = decays[0]
    hs = [ROOT.TH1F(f"h_delay_{n}", ";LLP #Deltat [ns];LLP decay products", 40, 0, 20) for n in range(3)]
    fill(hs[0], d.expDelay); fill(hs[1], d.slowness); fill(hs[2], d.pathDelay)
    for h, col in zip(hs, [ROOT.kBlack, ROOT.kBlue + 1, ROOT.kGreen + 2]): style(h, col)
    hs[0].SetMaximum(1.2 * max(h.GetMaximum() for h in hs))
    for n, h in enumerate(hs): h.Draw("hist" if n == 0 else "hist same")
    leg = legend(hs, ["total", "LLP slowness L/c(1/#beta-1)", "path length difference"], 0.45)
    leg.Draw()
    return hs + [leg]
page(canvas, draw_truth)

# ----- Validation: deltaR to nearest HB cluster, with and without the LLP frame shift -----
def draw_dR():
    d = decays[0]
    d = d[~d.llp_decaysInHB]
    hs = [ROOT.TH1F(f"h_dR_{n}", ";min #DeltaR(decay product, HB cluster);LLP decay products", 40, 0, 1) for n in range(2)]
    fill(hs[0], d.minDR_raw); fill(hs[1], d.minDR_shifted)
    style(hs[0], ROOT.kBlack); style(hs[1], ROOT.kMagenta + 1)
    hs[0].SetMaximum(1.2 * max(h.GetMaximum() for h in hs))
    hs[0].Draw("hist"); hs[1].Draw("hist same")
    cone = ROOT.TLine(DELTA_R, 0, DELTA_R, hs[0].GetMaximum())
    cone.SetLineStyle(2)
    cone.Draw()
    leg = legend(hs, ["from (0,0,0)", "shifted to LLP decay vertex"], 0.45)
    leg.AddEntry(cone, "matching cone", "l")
    leg.Draw()
    return hs + [cone, leg, match_text()]
page(canvas, draw_dR)

# ----- Per decay product vs expected delay -----
def draw_vs_delay(name, ytitle, value, ymin=None, ymax=None):
    def draw():
        hs = []
        for n, d in enumerate(decays):
            p = ROOT.TProfile(f"p_{name}_{n}", f";LLP #Deltat [ns];{ytitle}", len(DELAY_BINS) - 1, DELAY_BINS)
            fill(p, d.expDelay, value(d))
            hs.append(style(p, COLORS[n]))
        if ymin is not None: hs[0].SetMinimum(ymin)
        if ymax is not None: hs[0].SetMaximum(ymax)
        hs[0].Draw("e1")
        hs[1].Draw("e1 same")
        leg = legend(hs, LABELS)
        leg.Draw()
        return hs + [leg, match_text()]
    page(canvas, draw)

draw_vs_delay("eff", "fraction matched to an HB cluster", lambda d: ak.values_astype(d.isTruthMatched, np.float64), 0, 1.1)
draw_vs_delay("efrac", "#Sigma E_{matched clusters} / E_{b}", lambda d: d.matchedEnergy / d.energy, 0)
draw_vs_delay("nclus", "matched clusters per decay product", lambda d: d.nMatchedClusters, 0)

# ----- Event-by-event difference (same decay product in both PF versions) -----
def draw_diff():
    p = ROOT.TProfile("p_diff", ";LLP #Deltat [ns];(E_{timing} - E_{standard}) / E_{b}", len(DELAY_BINS) - 1, DELAY_BINS)
    fill(p, decays[0].expDelay, (decays[1].matchedEnergy - decays[0].matchedEnergy) / decays[0].energy)
    style(p, ROOT.kBlack)
    p.SetMinimum(-0.3)
    p.SetMaximum(0.3)
    p.Draw("e1")
    line = ROOT.TLine(DELAY_BINS[0], 0, DELAY_BINS[-1], 0)
    line.SetLineStyle(2)
    line.Draw()
    return [p, line, match_text()]
page(canvas, draw_diff)

# ----- Matched cluster time vs expected delay (clusters with a valid time) -----
for n, c in enumerate(hb_clus):
    def draw_time2d(n=n, c=c):
        m = c[(c.iGen >= 0) & (c.time > -100)]
        h = ROOT.TH2F(f"h_t2d_{n}", ";LLP #Deltat [ns];cluster time [ns]", 30, 0, 15, 40, -15, 25)
        fill(h, m.expDelay, m.time)
        canvas.SetRightMargin(0.14)
        h.Draw("colz")
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.032)
        n_valid = ak.sum(c[c.iGen >= 0].time > -100)
        text.DrawLatex(0.16, 0.85, f"{LABELS[n]}: {n_valid}/{ak.sum(c.iGen >= 0)} matched clusters have a valid time")
        text.DrawLatex(0.16, 0.80, f"dashed: time = #Deltat {TIME_OFFSET:+.2f} ns")
        line = ROOT.TLine(0, TIME_OFFSET, min(15., 25. - TIME_OFFSET), min(15., 25. - TIME_OFFSET) + TIME_OFFSET)
        line.SetLineStyle(2)
        line.Draw()
        return [h, text, line, match_text()]
    page(canvas, draw_time2d)
    canvas.SetRightMargin(0.1)

# ----- Cluster multiplicity per event -----
def draw_multiplicity():
    hs = []
    for n, c in enumerate(hb_clus_evt):
        h = ROOT.TH1F(f"h_mult_{n}", f";HB clusters per event (E > {args.clusterE:g} GeV);events", 50, 0, 100)
        fill(h, ak.num(c))
        hs.append(style(h, COLORS[n]))
    hs[0].SetMaximum(1.2 * max(h.GetMaximum() for h in hs))
    hs[0].Draw("hist")
    hs[1].Draw("hist same")
    leg = legend(hs, LABELS)
    leg.Draw()
    return hs + [leg]
page(canvas, draw_multiplicity)

# ----- Cluster energy and time: all, LLP-matched, and unmatched HB clusters -----
SELECTIONS = {
    "all":       ("all HB clusters",         lambda c: c.energy > 0, False),
    "matched":   ("LLP-matched HB clusters", lambda c: c.iGen >= 0,  True),
    "unmatched": ("unmatched HB clusters",   lambda c: c.iGen < 0,   True),
}

def draw_overlay(var, xtitle, nbins, lo, hi, selection, logy=True):
    sel_label, select, matched = SELECTIONS[selection]
    def draw():
        hs = []
        for n, c in enumerate(hb_clus):
            vals = c[select(c)][var]
            if var == "time": vals = vals[vals > -100]
            h = ROOT.TH1F(f"h_{var}_{selection}_{n}", f";{xtitle} ({sel_label});clusters", nbins, lo, hi)
            fill(h, vals)
            hs.append(style(h, COLORS[n]))
        canvas.SetLogy(logy)
        hs[0].SetMaximum(5 * max(h.GetMaximum() for h in hs) if logy else 1.2 * max(h.GetMaximum() for h in hs))
        hs[0].Draw("hist")
        hs[1].Draw("hist same")
        leg = legend(hs, LABELS)
        leg.Draw()
        return hs + [leg] + ([match_text()] if matched else [])
    page(canvas, draw)
    canvas.SetLogy(False)

for selection in SELECTIONS:
    draw_overlay("energy", "cluster energy [GeV]", 50, 0, 100, selection)
for selection in SELECTIONS:
    draw_overlay("time", "cluster time [ns]", 40, -15, 25, selection)

# ----- Per decay product means over its matched clusters (per event, for both averaging choices) -----
EVT = [(decay_std, clus_std), (decay_tim, clus_tim)]
USE = [matched_clusters(d, c, args.clusterE) for d, c in EVT]  # matched HB clusters with E > --clusterE
AVERAGES = {True: "energy-weighted mean", False: "simple mean"}

# ----- LLP decay depth vs energy-weighted mean depth of all matched clusters (LLPs decaying inside HB) -----
def draw_depth(n, weighted=True):
    def draw():
        d, c = EVT[n]
        depth = matched_mean(matched_clusters(d, c), c, c.depth, weighted)
        sel = d.inHBAcceptance & d.llp_decaysInHB & (depth >= 0)
        x = 1. + np.digitize(ak.to_numpy(ak.flatten(d.llp_decayR[sel])), args.hbDepthEdges)
        y = ak.to_numpy(ak.flatten(depth[sel])).astype(np.float64)

        top = ROOT.TPad(f"depth_top_{n}_{weighted}", "", 0, 0.36, 1, 1)
        bottom = ROOT.TPad(f"depth_bottom_{n}_{weighted}", "", 0, 0, 1, 0.36)
        for pad in (top, bottom):
            pad.SetLeftMargin(0.13)
            pad.Draw()
        bottom.SetBottomMargin(0.2)

        top.cd()
        ymax = max(4.5, y.max() + 0.5) if len(y) > 0 else 4.5
        frame = top.DrawFrame(0.5, 0.5, 4.5, ymax, f";LLP decay depth;matched cluster depth ({AVERAGES[weighted]})")
        frame.GetXaxis().SetNdivisions(4)
        diagonal = ROOT.TLine(0.5, 0.5, 4.5, 4.5)
        diagonal.SetLineStyle(2)
        diagonal.Draw()
        graph = ROOT.TGraph(len(x), x, y) if len(x) > 0 else ROOT.TGraph()
        graph.SetMarkerStyle(20)
        graph.SetMarkerSize(0.5)
        graph.SetMarkerColor(COLORS[n])
        graph.Draw("p same")
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        text.DrawLatex(0.16, 0.84, f"{LABELS[n]}: {len(x)} decay products (LLP decays in HB)")
        cone = match_text()

        bottom.cd()
        h = ROOT.TH1F(f"h_depthdiff_{n}_{weighted}", ";|cluster depth - LLP decay depth|;decay products", 30, 0, 3)
        fill(h, np.abs(y - x))
        style(h, COLORS[n])
        for axis in (h.GetXaxis(), h.GetYaxis()):
            axis.SetTitleSize(0.07)
            axis.SetLabelSize(0.06)
        h.GetYaxis().SetTitleOffset(0.8)
        h.Draw("hist")
        return [top, bottom, frame, diagonal, graph, text, cone, h]
    page(canvas, draw)

for n in range(2):
    draw_depth(n)

# ----- Correctly matched decay products per event, for 9 cut configurations -----
E_MID, T_MID, R_MID = args.energyFrac[1], args.timeCut[1], args.dRCut[1]
CONFIGS = ( [(e, T_MID, R_MID) for e in args.energyFrac]
          + [(E_MID, t, R_MID) for t in args.timeCut]
          + [(E_MID, T_MID, r) for r in args.dRCut] )
DENOMINATORS = {
    "all in HB acceptance":      lambda d, n_matched: d.inHBAcceptance,
    "#geq 1 matched cluster": lambda d, n_matched: d.inHBAcceptance & (n_matched > 0),
}

def draw_correct(den_label, weighted):
    tag = f"{list(DENOMINATORS).index(den_label)}_{weighted}"
    def draw():
        title = ROOT.TLatex()
        title.SetNDC()
        title.SetTextSize(0.02)
        title.DrawLatex(0.02, 0.975, f"Correctly matched decay products per event  |  denominator: {den_label}  |  {AVERAGES[weighted]} over clusters with E > {args.clusterE:g} GeV")
        grid = ROOT.TPad(f"grid_{tag}", "", 0, 0, 1, 0.965)
        grid.Draw()
        grid.Divide(3, 3)
        keep = [title, grid]

        # (sum E / E_b, mean distance to the time line, mean dR, number of clusters) per decay product, for each PF version
        quantities = []
        for (d, c), use in zip(EVT, USE):
            e_frac = ak.sum(use * c.energy[:, np.newaxis, :], axis=-1) / d.energy
            t_dist = matched_mean(use, c, abs(c.time - c.expDelay - TIME_OFFSET) / np.sqrt(2.), weighted, valid=c.time > -100)
            quantities.append((e_frac, t_dist, matched_mean(use, c, c.dR, weighted), ak.sum(use, axis=-1)))

        for i, (e_cut, t_cut, r_cut) in enumerate(CONFIGS):
            grid.cd(i + 1)
            ROOT.gPad.SetLeftMargin(0.14)
            ROOT.gPad.SetBottomMargin(0.13)
            hs, summary = [], []
            for n, (d, _) in enumerate(EVT):
                e_frac, t_dist, dr, n_matched = quantities[n]
                den = DENOMINATORS[den_label](d, n_matched)
                good = den & (e_frac > e_cut) & (t_dist >= 0) & (t_dist < t_cut) & (dr >= 0) & (dr < r_cut)
                n_den, n_good = ak.to_numpy(ak.sum(den, axis=1)), ak.to_numpy(ak.sum(good, axis=1))
                h = ROOT.TH1F(f"h_correct_{tag}_{i}_{n}", ";correctly matched per event [%];events", 21, -2.5, 102.5)
                fill(h, 100. * n_good[n_den > 0] / n_den[n_den > 0])
                hs.append(style(h, COLORS[n]))
                h.SetLineWidth(1)
                percent = 100. * n_good.sum() / n_den.sum() if n_den.sum() > 0 else 0.
                summary.append(f"{LABELS[n]}: {percent:.1f}% ({n_good.sum()}/{n_den.sum()})")
            hs[0].SetMaximum(1.7 * max(1., max(h.GetMaximum() for h in hs)))
            for axis in (hs[0].GetXaxis(), hs[0].GetYaxis()):
                axis.SetTitleSize(0.05)
                axis.SetLabelSize(0.045)
            hs[0].Draw("hist")
            hs[1].Draw("hist same")
            text = ROOT.TLatex()
            text.SetNDC()
            text.SetTextSize(0.045)
            text.DrawLatex(0.17, 0.85, f"#SigmaE/E_{{b}} > {e_cut:g},  d_{{t}} < {t_cut:g} ns,  #LT#DeltaR#GT < {r_cut:g}")
            for n, line in enumerate(summary):
                text.SetTextColor(COLORS[n])
                text.DrawLatex(0.17, 0.79 - 0.06 * n, line)
            keep += hs + [text]
        return keep
    page(canvas, draw)

for den_label in DENOMINATORS:
    for weighted in AVERAGES:
        draw_correct(den_label, weighted)

canvas.Print(args.output + "]")
print(f"Wrote {args.output}")