#!/usr/bin/env python3
"""
Plot PF cluster timing vs phase scan delay (laserType from uMNio).

Inputs: pfObjectsNtuple_phaseScan.root produced by NtuplePhaseScan.sh
Outputs: phase scan PDFs + ROOT file with histograms

Plots produced:
  - HCAL cluster time vs laserType  (2D + profile)
  - HBHE rechit time vs laserType   (2D + profile)
  - ECAL cluster time vs laserType  (2D + profile)
  - Mean cluster time vs laserType  (TGraph, one point per delay step)
"""

import argparse
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(ROOT.kBird)

parser = argparse.ArgumentParser(description="PF timing vs phase scan delay")
parser.add_argument("--input",  default="/eos/user/g/gkopp/PF_PhaseScan/pfObjectsNtuple_phaseScan.root",
                    help="Input ntuple ROOT file")
parser.add_argument("--output", default="phaseScan_timing.root",
                    help="Output ROOT file")
parser.add_argument("--pdf",    default="phaseScan_timing.pdf",
                    help="Output PDF")
parser.add_argument("--emin",   type=float, default=1.0,
                    help="Minimum cluster energy [GeV] for timing plots")
args = parser.parse_args()

# ── open input ────────────────────────────────────────────────────────────────

f_in = ROOT.TFile.Open(args.input)
if not f_in or f_in.IsZombie():
    raise RuntimeError(f"Cannot open {args.input}")
tree = f_in.Get("pfObjectsNtupler/pfTree")
if not tree:
    raise RuntimeError("TTree 'pfObjectsNtupler/pfTree' not found")

print(f"Opened {args.input}: {tree.GetEntries()} events")

# ── discover laserType range (fast C++-level scan, excludes sentinel -1000) ───

tree.Draw("laserType>>_hlt_range", "laserType > -999", "goff")
_hlt_range = ROOT.gDirectory.Get("_hlt_range")
lt_min = int(_hlt_range.GetXaxis().GetXmin())
lt_max = int(_hlt_range.GetXaxis().GetXmax())
n_lt_bins = lt_max - lt_min + 1

print(f"laserType range: {lt_min} – {lt_max}  ({n_lt_bins} bins)")

# ── book histograms ───────────────────────────────────────────────────────────

TIME_BINS, TIME_LO, TIME_HI = 100, -25.0, 25.0

h2_hcal = ROOT.TH2F("h2_hcal_time_vs_laser",
    "HCAL cluster time vs phase delay;laserType;Cluster time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_BINS, TIME_LO, TIME_HI)

h2_hbhe = ROOT.TH2F("h2_hbhe_time_vs_laser",
    "HBHE rechit time vs phase delay;laserType;Rechit time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_BINS, TIME_LO, TIME_HI)

h2_ecal = ROOT.TH2F("h2_ecal_time_vs_laser",
    "ECAL cluster time vs phase delay;laserType;Cluster time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_BINS, TIME_LO, TIME_HI)

# Profiles (mean time vs laserType)
prof_hcal = ROOT.TProfile("prof_hcal_time_vs_laser",
    "Mean HCAL cluster time vs phase delay;laserType;Mean cluster time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_LO, TIME_HI)

prof_hbhe = ROOT.TProfile("prof_hbhe_time_vs_laser",
    "Mean HBHE rechit time vs phase delay;laserType;Mean rechit time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_LO, TIME_HI)

prof_ecal = ROOT.TProfile("prof_ecal_time_vs_laser",
    "Mean ECAL cluster time vs phase delay;laserType;Mean cluster time [ns]",
    n_lt_bins, lt_min - 0.5, lt_max + 0.5,
    TIME_LO, TIME_HI)

# ── fill ──────────────────────────────────────────────────────────────────────

n_skip = 0
for ev in tree:
    lt = ev.laserType
    if lt == -1000:
        n_skip += 1
        continue

    for t, e in zip(ev.hcal_time, ev.hcal_energy):
        if e >= args.emin:
            h2_hcal.Fill(lt, t)
            prof_hcal.Fill(lt, t)

    for t, e in zip(ev.hbhe_rechit_time, ev.hbhe_rechit_energy):
        if e >= args.emin:
            h2_hbhe.Fill(lt, t)
            prof_hbhe.Fill(lt, t)

    for t, e in zip(ev.ecal_time, ev.ecal_energy):
        if e >= args.emin:
            h2_ecal.Fill(lt, t)
            prof_ecal.Fill(lt, t)

print(f"Skipped {n_skip} events with no uMNio laserType")

# ── helpers ───────────────────────────────────────────────────────────────────

def draw_cms_label(extra="2025 Phase Scan"):
    cms = ROOT.TLatex(); cms.SetNDC(); cms.SetTextSize(0.04)
    cms.DrawLatex(0.12, 0.935, "CMS")
    sub = ROOT.TLatex(); sub.SetNDC(); sub.SetTextSize(0.035)
    sub.DrawLatex(0.19, 0.935, "#bf{Preliminary}")
    rhs = ROOT.TLatex(); rhs.SetNDC(); rhs.SetTextFont(42)
    rhs.SetTextAlign(31); rhs.SetTextSize(0.035)
    rhs.DrawLatexNDC(0.90, 0.935, extra)
    return cms, sub, rhs


def draw_2d(canvas, h2, pdf):
    canvas.Clear()
    canvas.SetTopMargin(0.08)
    canvas.SetRightMargin(0.13)
    h2.Draw("COLZ")
    kept = draw_cms_label()
    canvas.Update()
    canvas.Print(pdf)
    canvas.SetRightMargin(0.10)
    return kept


def draw_profile(canvas, prof, pdf, color=ROOT.kBlue+1):
    canvas.Clear()
    canvas.SetTopMargin(0.08)
    prof.SetLineColor(color)
    prof.SetMarkerColor(color)
    prof.SetMarkerStyle(20)
    prof.SetMarkerSize(0.8)
    prof.SetLineWidth(2)
    prof.Draw("E1")
    kept = draw_cms_label()
    canvas.Update()
    canvas.Print(pdf)
    return kept

# ── draw ──────────────────────────────────────────────────────────────────────

canvas = ROOT.TCanvas("c", "", 800, 600)
canvas.Print(f"{args.pdf}[")

alive = []  # keep ROOT objects from being garbage-collected

alive += [draw_2d(canvas, h2_hcal, args.pdf)]
alive += [draw_profile(canvas, prof_hcal, args.pdf)]

alive += [draw_2d(canvas, h2_hbhe, args.pdf)]
alive += [draw_profile(canvas, prof_hbhe, args.pdf)]

alive += [draw_2d(canvas, h2_ecal, args.pdf)]
alive += [draw_profile(canvas, prof_ecal, args.pdf)]

# Combined profile overlay: HCAL + HBHE on same canvas
canvas.Clear()
canvas.SetTopMargin(0.08)

ylo = min(prof_hcal.GetMinimum(), prof_hbhe.GetMinimum()) - 2
yhi = max(prof_hcal.GetMaximum(), prof_hbhe.GetMaximum()) + 2

for i, (prof, label, color) in enumerate([
        (prof_hcal, "HCAL cluster", ROOT.kBlue+1),
        (prof_hbhe, "HBHE rechit",  ROOT.kRed+1),
        (prof_ecal, "ECAL cluster", ROOT.kGreen+2)]):
    prof.SetMinimum(ylo)
    prof.SetMaximum(yhi)
    prof.SetLineColor(color)
    prof.SetMarkerColor(color)
    prof.SetMarkerStyle(20 + i)
    prof.Draw("E1" if i == 0 else "E1 SAME")

leg = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.033)
leg.AddEntry(prof_hcal, "HCAL cluster", "lp")
leg.AddEntry(prof_hbhe, "HBHE rechit",  "lp")
leg.AddEntry(prof_ecal, "ECAL cluster", "lp")
leg.Draw()
alive += [draw_cms_label(), leg]
canvas.Update()
canvas.Print(args.pdf)

canvas.Print(f"{args.pdf}]")

# ── save histograms ───────────────────────────────────────────────────────────

f_out = ROOT.TFile(args.output, "RECREATE")
for obj in [h2_hcal, h2_hbhe, h2_ecal, prof_hcal, prof_hbhe, prof_ecal]:
    obj.Write()
f_out.Close()
f_in.Close()

print(f"Plots saved to {args.pdf}")
print(f"Histograms saved to {args.output}")
