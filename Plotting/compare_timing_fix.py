"""
Compare HCAL cluster timing before and after the fix that skips
rechits with time == -999 in the position-calc weighted average.

Usage:
  python compare_timing_fix.py --before pfObjectsNtuple_before.root \
                               --after  pfObjectsNtuple_after.root  \
                               --output timing_fix_comparison.pdf
"""

import argparse
import sys
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

SENTINEL = -999.0
SENTINEL_TOL = 0.5  # |t - (-999)| < 0.5 ns => treat as invalid


def is_sentinel(t):
    return abs(t - SENTINEL) < SENTINEL_TOL


def load_tree(path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"ERROR: cannot open {path}")
        sys.exit(1)
    tree = f.Get("pfObjectsNtupler/pfTree")
    if not tree:
        print(f"ERROR: no tree in {path}")
        sys.exit(1)
    return f, tree


def make_histos(tag):
    hists = {}
    # cluster time (full range including -999 sentinel)
    hists["cl_time_full"] = ROOT.TH1F(
        f"cl_time_full_{tag}", "HCAL cluster time (all);t [ns];Clusters", 220, -1010, 10
    )
    # cluster time in the contaminated range: between -999 and -10ns.
    # Before the fix, clusters with mixed-validity rechits fall here.
    # After the fix, this range should be empty.
    hists["cl_time_contaminated"] = ROOT.TH1F(
        f"cl_time_contaminated_{tag}",
        "HCAL cluster time (contaminated range);t [ns];Clusters",
        200, -998, -10,
    )
    # cluster time (zoom in on the physics range); 0.5 ns bins
    hists["cl_time_zoom"] = ROOT.TH1F(
        f"cl_time_zoom_{tag}", "HCAL cluster time (valid only);t [ns];Clusters", 40, -20, 20
    )
    # number of rechits per cluster
    hists["cl_nrechits"] = ROOT.TH1F(
        f"cl_nrechits_{tag}", "Rechits per HCAL cluster;N_{rechits};Clusters", 20, 0, 20
    )
    # fraction of rechits with time == -999, per cluster (rechit-level, unaffected by fix)
    hists["rechit_invalid_frac"] = ROOT.TH1F(
        f"rechit_invalid_frac_{tag}",
        "Fraction of rechits with time=-999 per cluster (rechit-level, same before/after);f_{invalid};Clusters",
        20, 0, 1.05,
    )
    # cluster time vs. E²-weighted mean valid-rechit time.
    # The position calc uses E²×fraction weighting (no timeResolutionCalc for standard HCAL).
    # Before fix: contaminated clusters appear off-diagonal (cluster t << E²-mean rechit t)
    # After fix: all points should lie along the diagonal
    hists["cl_vs_rh_time"] = ROOT.TH2F(
        f"cl_vs_rh_time_{tag}",
        f"({tag}) cluster time vs. E^{{2}}-weighted mean rechit time;E^{{2}}-weighted rechit t [ns];cluster t [ns]",
        40, -20, 20, 40, -20, 20,
    )
    # energy distribution of contaminated clusters (−998 to −10 ns).
    # Before fix: shows which energy scale is pulled into the bogus range.
    # After fix: should be empty.
    hists["cl_energy_contaminated"] = ROOT.TH1F(
        f"cl_energy_contaminated_{tag}",
        "HCAL cluster energy (contaminated timing, -998 to -10 ns);E [GeV];Clusters",
        50, 0, 50,
    )
    # energy distribution of clusters with valid timing (normalization reference)
    hists["cl_energy_valid"] = ROOT.TH1F(
        f"cl_energy_valid_{tag}",
        "HCAL cluster energy (valid timing);E [GeV];Clusters",
        50, 0, 50,
    )
    # eta distribution of contaminated clusters — shows which detector region is most affected
    hists["cl_eta_contaminated"] = ROOT.TH1F(
        f"cl_eta_contaminated_{tag}",
        "HCAL cluster #eta (contaminated timing, -998 to -10 ns);#eta;Clusters",
        50, -3.0, 3.0,
    )
    # cluster time (valid) vs cluster energy 2D.
    # Checks for energy-dependent bias after the fix; should be flat across all energies.
    hists["cl_time_vs_energy"] = ROOT.TH2F(
        f"cl_time_vs_energy_{tag}",
        f"({tag}) cluster time (valid) vs. energy;E [GeV];cluster t [ns]",
        50, 0, 50, 30, -15, 15,
    )
    return hists


def fill_histos(tree, hists):
    for ev in tree:
        nclust = len(ev.hcal_energy)
        # group (time, energy) by cluster index
        rh_data = [[] for _ in range(nclust)]
        for rh_t, rh_e, rh_cidx in zip(ev.hbhe_rechit_time, ev.hbhe_rechit_energy, ev.hbhe_rechit_clusterIndex):
            if 0 <= rh_cidx < nclust:
                rh_data[rh_cidx].append((rh_t, rh_e))

        for cidx, (cl_t, cl_e, cl_eta, cl_nrh) in enumerate(
            zip(ev.hcal_time, ev.hcal_energy, ev.hcal_eta, ev.hcal_nRecHits)
        ):
            hists["cl_time_full"].Fill(cl_t)
            hists["cl_nrechits"].Fill(cl_nrh)

            # contaminated range: between sentinel and physics range
            if -998 < cl_t < -10:
                hists["cl_time_contaminated"].Fill(cl_t)
                hists["cl_energy_contaminated"].Fill(cl_e)
                hists["cl_eta_contaminated"].Fill(cl_eta)

            if not is_sentinel(cl_t):
                hists["cl_time_zoom"].Fill(cl_t)
                hists["cl_energy_valid"].Fill(cl_e)
                hists["cl_time_vs_energy"].Fill(cl_e, cl_t)

            # rechit-level stats for this cluster
            data = rh_data[cidx]
            if data:
                n_invalid_rh = sum(1 for t, e in data if is_sentinel(t))
                hists["rechit_invalid_frac"].Fill(n_invalid_rh / len(data))

                # E²-weighted mean of valid rechits — matches the position calc weighting
                valid = [(t, e) for t, e in data if not is_sentinel(t)]
                if valid and not is_sentinel(cl_t):
                    w_sum = sum(e * e for t, e in valid)
                    if w_sum > 0:
                        e2_mean = sum(e * e * t for t, e in valid) / w_sum
                        hists["cl_vs_rh_time"].Fill(e2_mean, cl_t)


def draw_overlay(h_before, h_after, title, output, logy=False, draw_diag=False):
    c = ROOT.TCanvas("c", title, 800, 600)
    if logy:
        c.SetLogy()

    h_before.SetLineColor(ROOT.kBlue)
    h_before.SetLineWidth(2)
    h_after.SetLineColor(ROOT.kRed)
    h_after.SetLineWidth(2)

    mx = max(h_before.GetMaximum(), h_after.GetMaximum())
    h_before.SetMaximum(1.3 * mx)

    h_before.Draw("hist")
    h_after.Draw("hist same")

    leg = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(h_before, "before fix", "l")
    leg.AddEntry(h_after,  "after fix",  "l")
    leg.Draw()

    if draw_diag:
        diag = ROOT.TLine(
            h_before.GetXaxis().GetXmin(), h_before.GetXaxis().GetXmin(),
            h_before.GetXaxis().GetXmax(), h_before.GetXaxis().GetXmax(),
        )
        diag.SetLineStyle(2)
        diag.Draw()

    c.SaveAs(output)
    c.Close()


def draw_2d(h, output):
    c = ROOT.TCanvas("c2d", h.GetTitle(), 800, 700)
    c.SetRightMargin(0.13)
    h.Draw("COLZ")
    # diagonal guide
    mn = min(h.GetXaxis().GetXmin(), h.GetYaxis().GetXmin())
    mx = max(h.GetXaxis().GetXmax(), h.GetYaxis().GetXmax())
    diag = ROOT.TLine(mn, mn, mx, mx)
    diag.SetLineStyle(2)
    diag.Draw()
    c.SaveAs(output)
    c.Close()


def print_summary(tag, hists):
    total = hists["cl_time_full"].GetEntries()
    # exactly -999: all rechits invalid — same before and after the fix
    n_sentinel = hists["cl_time_full"].GetBinContent(
        hists["cl_time_full"].FindBin(SENTINEL)
    )
    # contaminated range (-998 to -10ns): mixed-validity clusters before fix; should be 0 after
    n_contaminated = int(hists["cl_time_contaminated"].GetEntries())
    # valid physics range (not sentinel, not contaminated)
    n_valid = int(total) - int(n_sentinel) - n_contaminated
    pct_s = 100.0 * n_sentinel / total if total > 0 else 0.0
    pct_c = 100.0 * n_contaminated / total if total > 0 else 0.0
    print(f"[{tag}] total HCAL clusters: {int(total)}")
    print(f"  time = -999 (all rechits invalid):   {int(n_sentinel):4d} ({pct_s:.1f}%)  -- expected same before/after")
    print(f"  time in (-998, -10) ns (contaminated): {n_contaminated:4d} ({pct_c:.1f}%)  -- should be ~0 after fix")
    print(f"  valid physics time:                  {n_valid:4d} ({100-pct_s-pct_c:.1f}%)")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", required=True, help="Ntuple before fix")
    parser.add_argument("--after",  required=True, help="Ntuple after fix")
    parser.add_argument("--output", default="timing_fix_comparison.pdf")
    args = parser.parse_args()

    fb, tb = load_tree(args.before)
    fa, ta = load_tree(args.after)

    hb = make_histos("before")
    ha = make_histos("after")

    print("Filling before...")
    fill_histos(tb, hb)
    print("Filling after...")
    fill_histos(ta, ha)

    print_summary("before", hb)
    print_summary("after",  ha)

    pdf = args.output
    ROOT.TCanvas().Print(pdf + "[")  # open multi-page PDF

    def overlay(key, title, logy=False):
        c = ROOT.TCanvas("c_" + key, title, 800, 600)
        if logy:
            c.SetLogy()
        h1 = hb[key]
        h2 = ha[key]
        h1.SetLineColor(ROOT.kBlue); h1.SetLineWidth(2)
        h2.SetLineColor(ROOT.kRed);  h2.SetLineWidth(2)
        mx = max(h1.GetMaximum(), h2.GetMaximum())
        h1.SetMaximum(1.3 * mx if not logy else mx * 10)
        h1.Draw("hist")
        h2.Draw("hist same")
        leg = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.AddEntry(h1, "before fix", "l")
        leg.AddEntry(h2, "after fix",  "l")
        leg.Draw()
        c.Print(pdf)

    overlay("cl_time_full",          "HCAL cluster time (full range, log scale)", logy=True)
    # Key diagnostic: before has entries here, after should be empty
    overlay("cl_time_contaminated",  "HCAL cluster time in contaminated range (-998 to -10 ns)")
    overlay("cl_time_zoom",          "HCAL cluster time (valid clusters only, -20 to 20 ns)")
    # Rechit-level — informational only, identical before/after by design
    overlay("rechit_invalid_frac",   "Fraction of invalid rechits per cluster (rechit-level, same before/after)")
    overlay("cl_nrechits",           "Rechits per cluster")
    # Energy/eta of contaminated clusters: before has entries, after should be empty
    overlay("cl_energy_contaminated", "HCAL cluster energy (contaminated timing, -998 to -10 ns)")
    overlay("cl_energy_valid",        "HCAL cluster energy (valid timing)")
    overlay("cl_eta_contaminated",    "HCAL cluster #eta (contaminated timing, -998 to -10 ns)")

    # 2D cluster time vs mean valid rechit time — drawn separately for before and after
    # Before fix: contaminated clusters appear with cluster t << mean rechit t (off-diagonal)
    # After fix: all points should lie on the diagonal
    for tag, hd in [("before", hb), ("after", ha)]:
        c2d = ROOT.TCanvas(f"c2d_{tag}", hd["cl_vs_rh_time"].GetTitle(), 800, 700)
        c2d.SetRightMargin(0.13)
        hd["cl_vs_rh_time"].Draw("COLZ")
        mn = hd["cl_vs_rh_time"].GetXaxis().GetXmin()
        mx_ax = hd["cl_vs_rh_time"].GetXaxis().GetXmax()
        diag = ROOT.TLine(mn, mn, mx_ax, mx_ax)
        diag.SetLineStyle(2)
        diag.Draw()
        c2d.Print(pdf)

    # 2D cluster time vs cluster energy — drawn separately for before and after
    # Checks that the fix introduces no energy-dependent timing bias
    for tag, hd in [("before", hb), ("after", ha)]:
        c2e = ROOT.TCanvas(f"c2e_{tag}", hd["cl_time_vs_energy"].GetTitle(), 800, 700)
        c2e.SetRightMargin(0.13)
        hd["cl_time_vs_energy"].Draw("COLZ")
        c2e.Print(pdf)

    ROOT.TCanvas().Print(pdf + "]")  # close PDF
    print(f"Saved {pdf}")

    fb.Close()
    fa.Close()


if __name__ == "__main__":
    main()
