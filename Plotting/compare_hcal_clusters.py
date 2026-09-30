#!/usr/bin/env python3
"""
Compare HCAL cluster properties across three PF algorithm variants:
  - standardPF
  - cellTimingPF
  - seedTimingPF

Source collection: particleFlowClusterHCAL (post-depth-stacking)

Plots:
  1. Number of HCAL clusters per event
  2. Total HCAL cluster energy per event
  3. Number of PF rechits per HCAL cluster (from hcal_nRecHits branch)

Expected input file naming:
  pfObjectsNtuple_standardPF${suffix}.root
  pfObjectsNtuple_cellTimingPF${suffix}.root
  pfObjectsNtuple_seedTimingPF${suffix}.root

For your current dipion scan, call this with:
  --suffix "_${SAMPLE}"

not:
  --suffix "_${TAG}"

because your ntuples are named with the full SAMPLE string.
"""

import argparse
import ROOT
import os

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser(description="Compare HCAL clusters across PF approaches")
parser.add_argument("--inputdir", default=".", help="Directory containing input ROOT files")
parser.add_argument("--prefix",   default="pfObjectsNtuple_", help="Filename prefix before the algorithm label")
parser.add_argument(
    "--suffix",
    default="",
    help='Filename suffix after the algorithm label before .root, e.g. "_${SAMPLE}"',
)
parser.add_argument("--output",   default="hcal_cluster_comparison.root", help="Output ROOT file")
parser.add_argument("--pdf",      default="hcal_cluster_comparison.pdf",  help="Output PDF with all plots")
args = parser.parse_args()

LABELS = ["standardPF", "cellTimingPF", "seedTimingPF"]
COLORS = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue]
STYLES = [1, 2, 7]  # solid, dashed, dot-dashed

# ── helpers ──────────────────────────────────────────────────────────────────

def open_tree(label):
    fname = os.path.join(args.inputdir, f"{args.prefix}{label}{args.suffix}.root")
    f = ROOT.TFile.Open(fname)
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open {fname}")
    t = f.Get("pfObjectsNtupler/pfTree")
    if not t:
        raise RuntimeError(f"TTree not found in {fname}")
    return f, t


def style_hist(h, color, lstyle):
    h.SetLineColor(color)
    h.SetLineStyle(lstyle)
    h.SetLineWidth(2)
    h.SetTitle("")
    return h


def make_legend(hists, labels):
    leg = ROOT.TLegend(0.55, 0.65, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.035)
    for h, lbl in zip(hists, labels):
        leg.AddEntry(h, lbl, "l")
    return leg


def draw_cms_label():
    """Draw the standard CMS Simulation label in the top margin."""
    cms = ROOT.TLatex()
    cms.SetNDC()
    cms.SetTextSize(0.04)
    cms.DrawLatex(0.12, 0.935, "CMS")

    sim = ROOT.TLatex()
    sim.SetNDC()
    sim.SetTextSize(0.035)
    sim.DrawLatex(0.19, 0.935, "#bf{Simulation}")

    coll = ROOT.TLatex()
    coll.SetTextFont(42)
    coll.SetTextAlign(31)
    coll.SetTextSize(0.035)
    coll.DrawLatexNDC(0.90, 0.935, "particleFlowClusterHCAL")
    return cms, sim, coll  # keep alive


def draw_overlay(canvas, hists, labels, xtitle, logy=False):
    """Normalise to unit area and overlay histograms. Returns clones."""
    canvas.Clear()
    canvas.SetTopMargin(0.08)
    canvas.SetLogy(1 if logy else 0)

    normed = []
    for h in hists:
        hn = h.Clone(h.GetName() + "_norm")
        integral = hn.Integral()
        if integral > 0:
            hn.Scale(1.0 / integral)
        normed.append(hn)

    ymax = max(hn.GetMaximum() for hn in normed) * 1.4
    for i, hn in enumerate(normed):
        hn.GetXaxis().SetTitle(xtitle)
        hn.GetYaxis().SetTitle("Normalised entries")
        hn.SetMaximum(ymax)
        hn.Draw("HIST" if i == 0 else "HIST SAME")

    leg = make_legend(normed, labels)
    leg.Draw()
    canvas.Update()
    return normed, leg  # keep alive


# ── book histograms ───────────────────────────────────────────────────────────

h_ncl   = {}   # N clusters per event
h_etot  = {}   # total cluster energy per event
h_nhits = {}   # hits per cluster

sum_cluster_energy = {}
sum_rechit_energy  = {}
n_events           = {}

for label in LABELS:
    h_ncl[label] = ROOT.TH1F(
        f"h_ncl_{label}",
        f"{label} — HCAL clusters/event;N_{{clusters}};Entries",
        8, 0, 8,
    )

    h_etot[label] = ROOT.TH1F(
        f"h_etot_{label}",
        f"{label} — HCAL total cluster energy/event;#SigmaE [GeV];Entries",
        30, 0, 60,
    )

    h_nhits[label] = ROOT.TH1F(
        f"h_nhits_{label}",
        f"{label} — PF rechits per HCAL cluster;N_{{PF rechits}};Entries",
        10, 0, 20,
    )


# ── fill histograms ───────────────────────────────────────────────────────────

files = {}

for label in LABELS:
    try:
        f, tree = open_tree(label)
    except RuntimeError as e:
        print(f"WARNING: {e} — skipping {label}")
        continue

    files[label] = f

    print(f"Processing {label}: {tree.GetEntries()} events")

    sum_cluster_energy[label] = 0.0
    sum_rechit_energy[label]  = 0.0
    n_events[label]           = 0

    for event in tree:
        n_cl = len(event.hcal_energy)
        h_ncl[label].Fill(n_cl)

        ecl = sum(event.hcal_energy)
        erh = sum(event.hbhe_pfrh_energy)

        h_etot[label].Fill(ecl)

        sum_cluster_energy[label] += ecl
        sum_rechit_energy[label]  += erh
        n_events[label]           += 1

        for n in event.hcal_nRecHits:
            h_nhits[label].Fill(n)


# ── energy conservation cross-check ──────────────────────────────────────────

print()
print(
    f"{'Algorithm':<22} "
    f"{'Events':>7}  "
    f"{'Mean cluster ΣE [GeV]':>22}  "
    f"{'Mean raw rechit ΣE [GeV]':>24}  "
    f"{'Offset (calib) [GeV]':>20}"
)
print("-" * 103)

for label in LABELS:
    if label not in n_events or n_events[label] == 0:
        continue

    n = n_events[label]
    ec = sum_cluster_energy[label] / n
    er = sum_rechit_energy[label] / n

    print(f"{label:<22} {n:>7}  {ec:>22.3f}  {er:>24.3f}  {ec - er:>20.3f}")

print()


# ── apply styles ──────────────────────────────────────────────────────────────

for i, label in enumerate(LABELS):
    for d in [h_ncl, h_etot, h_nhits]:
        style_hist(d[label], COLORS[i], STYLES[i])


# ── draw and save ─────────────────────────────────────────────────────────────

active = [label for label in LABELS if label in files]

if len(active) == 0:
    raise RuntimeError("No valid input files found. Check --inputdir, --prefix, and --suffix.")

out = ROOT.TFile(args.output, "RECREATE")
canvas = ROOT.TCanvas("c", "", 800, 600)
canvas.Print(f"{args.pdf}[")

kept = []

n1, l1 = draw_overlay(
    canvas,
    [h_ncl[l] for l in active],
    active,
    "N_{clusters} per event",
)
kept += [n1, l1, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

n2, l2 = draw_overlay(
    canvas,
    [h_etot[l] for l in active],
    active,
    "#SigmaE per event [GeV]",
)
kept += [n2, l2, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

n3, l3 = draw_overlay(
    canvas,
    [h_nhits[l] for l in active],
    active,
    "N_{PF rechits} per cluster",
    logy=True,
)
kept += [n3, l3, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

canvas.Print(f"{args.pdf}]")

# Write raw, unnormalised histograms to ROOT file
out.cd()
for label in active:
    h_ncl[label].Write()
    h_etot[label].Write()
    h_nhits[label].Write()

out.Close()

print(f"Plots saved to {args.pdf}")
print(f"Histograms saved to {args.output}")