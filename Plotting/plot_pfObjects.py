import argparse
import ROOT

parser = argparse.ArgumentParser(description="Plot PF objects from ntupler output")
parser.add_argument("--input",  default="pfObjectsNtuple.root",  help="Input ROOT file from ntupler")
parser.add_argument("--output", default="pfObjectsHistos.root",  help="Output ROOT file for histograms")
args = parser.parse_args()

# Load the ROOT file and TTree
file = ROOT.TFile.Open(args.input)
tree = file.Get("pfObjectsNtupler/pfTree")

if not tree:
    print("Error: Could not find TTree 'pfObjectsNtupler/pfTree' in file!")
    exit(1)

# Optional: check entries
print(f"Entries in tree: {tree.GetEntries()}")

# Create output file for histograms
out = ROOT.TFile(args.output, "RECREATE")

# Define histograms
h_pf_pt   = ROOT.TH1F("h_pf_pt",   "PF candidate p_{T};p_{T} [GeV];Entries", 100, 0, 200)
h_pf_eta  = ROOT.TH1F("h_pf_eta",  "PF candidate #eta;#eta;Entries", 100, -5, 5)
h_ecal_E  = ROOT.TH1F("h_ecal_E",  "ECAL cluster energy;E [GeV];Entries", 100, 0, 200)
# HCAL clusters (particleFlowClusterHCAL, post-depth-stacking)
h_hcal_E_all    = ROOT.TH1F("h_hcal_E_all",    "HCAL cluster energy;E [GeV];Entries", 100, 0, 200)
h_hcal_time_all = ROOT.TH1F("h_hcal_time_all", "HCAL cluster time;t [ns];Entries",    100, -50, 50)

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
    for e, t in zip(event.hcal_energy, event.hcal_time):
        h_hcal_E_all.Fill(e)
        h_hcal_time_all.Fill(t)

# Save histograms
out.cd()
out.Write()
out.Close()

print(f"Histograms saved to {args.output}")
