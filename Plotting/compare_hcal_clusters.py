#!/usr/bin/env python3
"""
Compare HCAL cluster properties across four PF algorithm variants:
  - standardPF
  - cellTimingPF
  - seedTimingPF
  - depth1SeedTimingPF

Source collection: particleFlowClusterHCAL (post-depth-stacking)

Plots:
  1. Number of HCAL clusters per event
  2. Total HCAL cluster energy per event
  3. Number of PF rechits per HCAL cluster (from hcal_nRecHits branch)
"""

import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

parser = argparse.ArgumentParser(description="Compare HCAL clusters across PF approaches")
parser.add_argument("--inputdir", default=".", help="Directory containing pfObjectsNtuple_*.root files")
parser.add_argument("--output",   default="hcal_cluster_comparison.root", help="Output ROOT file")
parser.add_argument("--pdf",      default="hcal_cluster_comparison.pdf",  help="Output PDF with all plots")
args = parser.parse_args()

LABELS = ["standardPF", "cellTimingPF", "seedTimingPF", "depth1SeedTimingPF"]
COLORS = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2]
STYLES = [1, 2, 7, 9]  # solid, dashed, dot-dashed, dotted

# ── helpers ──────────────────────────────────────────────────────────────────

def open_tree(label):
    fname = f"{args.inputdir}/pfObjectsNtuple_{label}.root"
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
    # sim.SetTextFont(42) 
    sim.SetTextSize(0.035)
    # sim.DrawLatex(0.19, 0.935, "#it{Simulation}")
    sim.DrawLatex(0.19, 0.935, "#bf{Simulation}")

    coll = ROOT.TLatex()
    coll.SetTextFont(42)
    coll.SetTextAlign(31)
    coll.SetTextSize(0.035)
    coll.DrawLatexNDC(0.90, 0.935, "particleFlowClusterHCAL")
    return cms, sim, coll  # keep alive


def draw_overlay(canvas, hists, labels, xtitle, logy=False):
    """Normalise to unit area and overlay histograms. Returns clones (keep-alive)."""
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

# particleFlowClusterHCAL (post-depth-stacking)
h_ncl   = {}   # N clusters per event
h_etot  = {}   # total cluster energy per event
h_nhits = {}   # hits per cluster (uses hbhe_rechit_clusterIndex)

# accumulators for energy conservation cross-check
sum_cluster_energy  = {}   # sum of hcal_energy per event, accumulated over events
sum_rechit_energy   = {}   # sum of hbhe_rechit_energy per event, accumulated over events
n_events            = {}

for label in LABELS:
    h_ncl[label]   = ROOT.TH1F(f"h_ncl_{label}",
        f"{label} — HCAL clusters/event;N_{{clusters}};Entries",
        30, 0, 30)
    h_etot[label]  = ROOT.TH1F(f"h_etot_{label}",
        f"{label} — HCAL total cluster energy/event;#SigmaE [GeV];Entries",
        50, 0, 150)
    h_nhits[label] = ROOT.TH1F(f"h_nhits_{label}",
        f"{label} — PF rechits per HCAL cluster;N_{{PF rechits}};Entries",
        50, 0, 50)

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
        erh = sum(event.hbhe_rechit_energy)
        h_etot[label].Fill(ecl)

        sum_cluster_energy[label] += ecl
        sum_rechit_energy[label]  += erh
        n_events[label]           += 1

        for n in event.hcal_nRecHits:
            h_nhits[label].Fill(n)

# ── energy conservation cross-check ──────────────────────────────────────────

# Column meanings:
#   Mean cluster ΣE  — mean of sum(hcal_energy) per event; should be equal across algorithms
#   Mean raw rechit ΣE — mean of sum(hbhe_rechit_energy) per event; should be equal across
#                        algorithms (same hbhereco input) and is a sanity check that the same
#                        events/rechits were used; the offset vs cluster energy is calibration
print()
print(f"{'Algorithm':<22} {'Events':>7}  {'Mean cluster ΣE [GeV]':>22}  {'Mean raw rechit ΣE [GeV]':>24}  {'Offset (calib) [GeV]':>20}")
print("-" * 103)
for label in LABELS:
    if label not in n_events or n_events[label] == 0:
        continue
    n  = n_events[label]
    ec = sum_cluster_energy[label] / n
    er = sum_rechit_energy[label]  / n
    print(f"{label:<22} {n:>7}  {ec:>22.3f}  {er:>24.3f}  {ec - er:>20.3f}")
print()

# ── apply styles ──────────────────────────────────────────────────────────────

for i, label in enumerate(LABELS):
    for d in [h_ncl, h_etot, h_nhits]:
        style_hist(d[label], COLORS[i], STYLES[i])

# ── draw and save ─────────────────────────────────────────────────────────────

out = ROOT.TFile(args.output, "RECREATE")
canvas = ROOT.TCanvas("c", "", 800, 600)
canvas.Print(f"{args.pdf}[")

active = [l for l in LABELS if l in files]

kept = []  # keep normalised hists and legends alive for the duration of printing

n1, l1 = draw_overlay(canvas,
    [h_ncl[l] for l in active], active,
    "N_{clusters} per event")
kept += [n1, l1, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

n2, l2 = draw_overlay(canvas,
    [h_etot[l] for l in active], active,
    "#SigmaE per event [GeV]")
kept += [n2, l2, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

n3, l3 = draw_overlay(canvas,
    [h_nhits[l] for l in active], active,
    "N_{PF rechits} per cluster",
    logy=True)
kept += [n3, l3, draw_cms_label()]
canvas.Update()
canvas.Print(args.pdf)

canvas.Print(f"{args.pdf}]")

# Write raw (unnormalised) histograms to ROOT file
out.cd()
for label in active:
    h_ncl[label].Write()
    h_etot[label].Write()
    h_nhits[label].Write()
out.Close()

print(f"Plots saved to {args.pdf}")
print(f"Histograms saved to {args.output}")
