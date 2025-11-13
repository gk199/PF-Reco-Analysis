import ROOT

# Load the ROOT file and TTree
file = ROOT.TFile.Open("pfObjectsNtuple.root")
tree = file.Get("pfObjectsNtupler/pfTree")

if not tree:
    print("Error: Could not find TTree 'pfObjectsNtupler/pfTree' in file!")
    exit(1)

# Optional: check entries
print(f"Entries in tree: {tree.GetEntries()}")

# Create output file for histograms
out = ROOT.TFile("pfObjectsHistos.root", "RECREATE")

# Define histograms
h_pf_pt   = ROOT.TH1F("h_pf_pt",   "PF candidate p_{T};p_{T} [GeV];Entries", 100, 0, 200)
h_pf_eta  = ROOT.TH1F("h_pf_eta",  "PF candidate #eta;#eta;Entries", 100, -5, 5)
h_ecal_E  = ROOT.TH1F("h_ecal_E",  "ECAL cluster energy;E [GeV];Entries", 100, 0, 200)
h_hcal_E  = ROOT.TH1F("h_hcal_E",  "HCAL cluster energy;E [GeV];Entries", 100, 0, 200)
# (all depths)
h_hcal_E_all   = ROOT.TH1F("h_hcal_E_all",   "HCAL cluster energy (all depths);E [GeV];Entries", 100, 0, 200)
h_hcal_time_all = ROOT.TH1F("h_hcal_time_all", "HCAL cluster time (all depths);t [ns];Entries", 100, -50, 50)
# Prepare per-depth histograms (depths 1–7 are typical in HCAL, but we can detect dynamically)
depth_hists_E = {}
depth_hists_t = {}

# Loop over entries
for event in tree:
    # PF candidates (vectors)
    for pt, eta in zip(event.pf_pt, event.pf_eta):
        h_pf_pt.Fill(pt)
        h_pf_eta.Fill(eta)
    # ECAL clusters
    for e in event.ecal_energy:
        h_ecal_E.Fill(e)
    # HCAL clusters
    for e in event.hcal_energy:
        h_hcal_E.Fill(e)

    # also do per depth histograms for hcal
    nHcal = len(event.hcal_energy)
    for j in range(nHcal):
        e = event.hcal_energy[j]
        t = event.hcal_time[j]
        d = int(event.hcal_depth[j])

        # Fill overall histos
        h_hcal_E_all.Fill(e)
        h_hcal_time_all.Fill(t)

        # Create depth histos dynamically if not existing
        if d not in depth_hists_E:
            depth_hists_E[d] = ROOT.TH1F(f"h_hcal_E_depth{d}", f"HCAL energy depth {d};E [GeV];Entries", 100, 0, 200)
            depth_hists_t[d] = ROOT.TH1F(f"h_hcal_time_depth{d}", f"HCAL time depth {d};t [ns];Entries", 100, -50, 50)
        
        depth_hists_E[d].Fill(e)
        depth_hists_t[d].Fill(t)

# Save histograms
out.Write()
out.Close()

print("Histograms saved to pfObjectsHistos.root")
