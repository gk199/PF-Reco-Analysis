#!/usr/bin/env python3
"""
Overlay nominal PF vs. timing PF jet and MET resolution.

Reads the histogram files written by plot_jetmet_resolution.py (one per PF
variant) and overlays them, one page per observable, with a timing/nominal
ratio panel underneath — same layout as Plotting/compare_qcd_ratio.py.
Mean and RMS per variant are printed in the legend, since for a resolution
study those two numbers are the actual result and the shape is the context.

Usage:
  python3 compare_jetmet_resolution.py \
      --nominal resolution_nominal.root \
      --timing  resolution_timing3ns.root \
      --labels  3ns \
      --output  resolution_nominal_vs_timing.pdf
"""

import argparse
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

COLORS = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen + 2, ROOT.kMagenta + 1, ROOT.kOrange + 1]

# key -> (histogram name in the input files, x-axis title, page title)
OBSERVABLES = [
    ("jet_res", "h_jet_res",
     "(Jet p_{T} - GenJet p_{T}) / GenJet p_{T}", "Jet p_{T} resolution"),
    ("pf_met_res", "h_pf_met_res",
     "PF MET p_{T} - GenMET p_{T} [GeV]", "PF MET resolution"),
    ("puppi_met_res", "h_puppi_met_res",
     "Puppi MET p_{T} - GenMET p_{T} [GeV]", "Puppi MET resolution"),
]


def load_hist(rootfile, name, tag):
    """Read one histogram and detach it from the file, so it survives the
    TFile going out of scope."""
    h = rootfile.Get(name)
    if not h:
        print(f"ERROR: no histogram '{name}' in {rootfile.GetName()}")
        sys.exit(1)
    h = h.Clone(f"{name}_{tag}")
    h.SetDirectory(0)
    # plot_jetmet_resolution.py fills without Sumw2, so the input histograms
    # carry no stored errors and TH1::Divide falls back to a nonsense error
    # (~100% per bin). Enabling it here seeds sumw2 = content, i.e. the sqrt(N)
    # errors these unweighted fills should have, and keeps --normalize honest.
    h.Sumw2()
    return h


def open_file(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"ERROR: cannot open {path}")
        sys.exit(1)
    return f


def draw_overlay_page(pdf, key, xtitle, title, h_nom, timing_hists, labels,
                      logy=False, ratio_range=(0.5, 1.5), ratio_min_entries=10.0,
                      normalized=False):
    """Top pad: nominal + each timing variant overlaid.
    Bottom pad: ratio (timing / nominal) per bin, with a line at 1.0."""
    c = ROOT.TCanvas(f"c_{key}", key, 800, 700)
    pad_top = ROOT.TPad(f"top_{key}", "", 0, 0.32, 1, 1)
    pad_bot = ROOT.TPad(f"bot_{key}", "", 0, 0.0, 1, 0.32)
    pad_top.SetLeftMargin(0.13)
    pad_top.SetBottomMargin(0.02)
    pad_bot.SetLeftMargin(0.13)
    pad_bot.SetTopMargin(0.03)
    pad_bot.SetBottomMargin(0.32)
    pad_top.SetLogy(logy)
    pad_top.Draw()
    pad_bot.Draw()

    pad_top.cd()
    h_nom.SetLineColor(ROOT.kBlack)
    h_nom.SetLineWidth(2)
    h_nom.SetTitle("")
    h_nom.GetXaxis().SetLabelSize(0)
    h_nom.GetXaxis().SetTitle("")
    h_nom.GetYaxis().SetTitle("Normalized entries" if normalized else "Entries")
    h_nom.GetYaxis().SetTitleSize(0.05)
    h_nom.GetYaxis().SetTitleOffset(1.25)
    h_nom.GetYaxis().SetLabelSize(0.04)

    ymax = max([h_nom.GetMaximum()] + [h.GetMaximum() for h in timing_hists])
    h_nom.SetMaximum(ymax * (30 if logy else 1.5))
    if not logy:
        h_nom.SetMinimum(0)
    h_nom.Draw("hist")
    for h, color in zip(timing_hists, COLORS):
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.Draw("hist same")

    leg = ROOT.TLegend(0.45, 0.64, 0.90, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.030)
    leg.AddEntry(h_nom, f"nominal PF  #mu={h_nom.GetMean():.3f}, RMS={h_nom.GetRMS():.3f}", "l")
    for h, lbl in zip(timing_hists, labels):
        leg.AddEntry(h, f"timing PF ({lbl})  #mu={h.GetMean():.3f}, RMS={h.GetRMS():.3f}", "l")
    leg.Draw()

    cms = ROOT.TLatex()
    cms.SetNDC()
    cms.SetTextSize(0.05)
    cms.DrawLatex(0.13, 0.94, "#bf{CMS} #it{Simulation}")
    hdr = ROOT.TLatex()
    hdr.SetNDC()
    hdr.SetTextSize(0.042)
    hdr.DrawLatex(0.55, 0.94, title)

    pad_bot.cd()
    ratios = []
    for i, (h, color) in enumerate(zip(timing_hists, COLORS)):
        r = h.Clone(f"ratio_{key}_{i}")
        r.SetDirectory(0)
        r.Divide(h_nom)
        # Blank bins where the nominal is too sparse to say anything: their
        # ratio errors are ~100% and the tall bars swamp the bulk of the
        # distribution, which is where the real difference lives. Raw counts
        # are recovered as (content/error)^2 so this works after --normalize.
        for b in range(1, r.GetNbinsX() + 1):
            err = h_nom.GetBinError(b)
            raw = (h_nom.GetBinContent(b) / err) ** 2 if err > 0 else 0.0
            if raw < ratio_min_entries:
                r.SetBinContent(b, 0.0)
                r.SetBinError(b, 0.0)
        r.SetTitle("")
        r.SetLineColor(color)
        r.SetMarkerColor(color)
        r.SetMarkerStyle(20)
        r.SetMarkerSize(0.7)
        r.GetYaxis().SetTitle("timing / nominal")
        r.GetYaxis().SetTitleSize(0.11)
        r.GetYaxis().SetTitleOffset(0.42)
        r.GetYaxis().SetLabelSize(0.09)
        r.GetYaxis().SetNdivisions(505)
        r.GetXaxis().SetTitle(xtitle)
        r.GetXaxis().SetTitleSize(0.12)
        r.GetXaxis().SetLabelSize(0.10)
        r.SetMinimum(ratio_range[0])
        r.SetMaximum(ratio_range[1])
        r.Draw("E1" if i == 0 else "E1 same")
        ratios.append(r)

    xmin = h_nom.GetXaxis().GetXmin()
    xmax = h_nom.GetXaxis().GetXmax()
    line = ROOT.TLine(xmin, 1.0, xmax, 1.0)
    line.SetLineStyle(2)
    line.SetLineColor(ROOT.kGray + 2)
    line.Draw()

    c.Print(pdf)
    return c, pad_top, pad_bot, ratios, leg, line, cms, hdr  # keep alive


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--nominal", default="resolution_nominal.root",
                         help="Histogram file for nominal PF (from plot_jetmet_resolution.py)")
    parser.add_argument("--timing", nargs="+", default=["resolution_timing3ns.root"],
                         help="Histogram file(s) for timing PF, one per threshold")
    parser.add_argument("--labels", nargs="+", default=["3ns"],
                         help="Label per --timing file, e.g. 1ns 3ns 5ns")
    parser.add_argument("--output", default="resolution_nominal_vs_timing.pdf", help="Output PDF")
    parser.add_argument("--rootfile", default="resolution_nominal_vs_timing.root",
                         help="Output ROOT file with the overlaid histograms and ratios")
    parser.add_argument("--normalize", action="store_true",
                         help="Scale every histogram to unit area before overlaying. Use when the "
                              "variants were run over different numbers of events — otherwise the "
                              "ratio panel measures statistics, not the PF change.")
    parser.add_argument("--logy", action="store_true", help="Log scale on the overlay pads")
    parser.add_argument("--ratio-min", type=float, default=0.5, help="Ratio panel y-axis minimum")
    parser.add_argument("--ratio-max", type=float, default=1.5, help="Ratio panel y-axis maximum")
    parser.add_argument("--ratio-min-entries", type=float, default=10.0,
                         help="Skip ratio points where the nominal bin has fewer raw entries "
                              "than this (their errors are ~100%% and hide the bulk)")
    args = parser.parse_args()

    if len(args.timing) != len(args.labels):
        print("ERROR: --timing and --labels must have the same length")
        sys.exit(1)

    f_nom = open_file(args.nominal)
    files_timing = [open_file(p) for p in args.timing]

    # Load everything up front so the summary table can be printed before any
    # drawing, and so normalization is applied consistently across variants.
    hists = {}  # key -> (h_nominal, [h_timing...])
    for key, hname, _, _ in OBSERVABLES:
        h_nom = load_hist(f_nom, hname, "nominal")
        h_timing = [load_hist(f, hname, f"timing_{lbl}")
                    for f, lbl in zip(files_timing, args.labels)]
        if args.normalize:
            for h in [h_nom] + h_timing:
                if h.Integral() > 0:
                    h.Scale(1.0 / h.Integral())
        hists[key] = (h_nom, h_timing)

    print(f"{'observable':<16} {'variant':<16} {'entries':>10} {'mean':>10} {'RMS':>10} "
          f"{'RMS/nominal':>12}")
    for key, _, _, title in OBSERVABLES:
        h_nom, h_timing = hists[key]
        print(f"{key:<16} {'nominal':<16} {h_nom.GetEntries():>10.0f} "
              f"{h_nom.GetMean():>10.4f} {h_nom.GetRMS():>10.4f} {1.0:>12.4f}")
        for h, lbl in zip(h_timing, args.labels):
            frac = h.GetRMS() / h_nom.GetRMS() if h_nom.GetRMS() > 0 else float("nan")
            print(f"{'':<16} {'timing ' + lbl:<16} {h.GetEntries():>10.0f} "
                  f"{h.GetMean():>10.4f} {h.GetRMS():>10.4f} {frac:>12.4f}")

    pdf = args.output
    ROOT.TCanvas().Print(pdf + "[")

    kept = []
    for key, _, xtitle, title in OBSERVABLES:
        h_nom, h_timing = hists[key]
        kept.append(draw_overlay_page(
            pdf, key, xtitle, title, h_nom, h_timing, args.labels,
            logy=args.logy, ratio_range=(args.ratio_min, args.ratio_max),
            ratio_min_entries=args.ratio_min_entries, normalized=args.normalize))

    ROOT.TCanvas().Print(pdf + "]")

    out = ROOT.TFile(args.rootfile, "RECREATE")
    for key, _, _, _ in OBSERVABLES:
        h_nom, h_timing = hists[key]
        h_nom.Write()
        for h in h_timing:
            h.Write()
    for entry in kept:
        for r in entry[3]:
            r.Write()
    out.Close()

    print(f"\nPlots saved to {pdf}")
    print(f"Histograms saved to {args.rootfile}")


if __name__ == "__main__":
    main()
