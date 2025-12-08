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
h_hcal_E  = ROOT.TH1F("h_hcal_E",  "HCAL cluster energy;E [GeV];Entries", 100, 0, 20)

# Event#, Run #, lumi section histograms
event_number = ROOT.TH1F("event", "Event Number;Event #;Entries", 100, 0, 1000)
run_number   = ROOT.TH1F("run",   "Run Number;Run #;Entries", 100, 0, 1000)
lumi_section = ROOT.TH1F("lumi", "Lumi Section;Lumi Section;Entries", 100, 0, 1e4)

# Ecal 
ecal_energy = ROOT.TH1F("ecal_energy", "ECAL cluster energy;E [GeV];Entries", 100, 0, 40)
ecal_eta    = ROOT.TH1F("ecal_eta",    "ECAL cluster eta;#eta;Entries", 100, -5, 5)
ecal_phi    = ROOT.TH1F("ecal_phi",    "ECAL cluster phi;#phi;Entries", 100, -3.5, 3.5)
ecal_time   = ROOT.TH1F("ecal_time",   "ECAL cluster time;t [ns];Entries", 100, -50, 50)
ecal_clusterIdx = ROOT.TH1F("ecal_clusterIdx", "ECAL cluster index", 100, 0, 100)

#eb rechits
eb_rechit_energy = ROOT.TH1F("eb_rechit_energy", "EB rechit energy;E [GeV];Entries", 50, 0, 15)
eb_rechit_eta    = ROOT.TH1F("eb_rechit_eta",    "EB rechit eta;#eta;Entries", 100, -5, 5)
eb_rechit_phi    = ROOT.TH1F("eb_rechit_phi",    "EB rechit phi;#phi;Entries", 100, -3.5, 3.5)
eb_rechit_time   = ROOT.TH1F("eb_rechit_time",   "EB rechit time;t [ns];Entries", 100, -50, 50)
eb_rechit_clusterIdx = ROOT.TH1F("eb_rechit_clusterIdx", "EB rechit cluster index", 100, 0, 100)
eb_rechit_counts  = ROOT.TH1F("eb_rechit_counts", "EB rechit count per cluster", 100, 0, 1000)

#ee rechits
ee_rechit_energy = ROOT.TH1F("ee_rechit_energy", "EE rechit energy;E [GeV];Entries", 100, 0, 200)
ee_rechit_time   = ROOT.TH1F("ee_rechit_time",   "EE rechit time;t [ns];Entries", 100, -50, 50)
ee_rechit_clusterIdx = ROOT.TH1F("ee_rechit_clusterIdx", "EE rechit cluster index", 100, 0, 100)
ee_rechit_eta    = ROOT.TH1F("ee_rechit_eta",    "EE rechit eta;#eta;Entries", 100, -5, 5)
ee_rechit_phi    = ROOT.TH1F("ee_rechit_phi",    "EE rechit phi;#phi;Entries", 100, -3.5, 3.5)
ee_rechit_counts  = ROOT.TH1F("ee_rechit_counts", "EE rechit count per cluster", 100, 0, 1000)

# HCAL clusters
# (all depths)
h_hcal_E_all   = ROOT.TH1F("h_hcal_E_all",   "HCAL cluster energy (all depths);E [GeV];Entries", 100, 0, 200)
h_hcal_time_all = ROOT.TH1F("h_hcal_time_all", "HCAL cluster time (all depths);t [ns];Entries", 100, -50, 50)
hcal_depth = ROOT.TH1F("hcal_depth", "HCAL cluster depth", 80,0,8)
hcal_eta = ROOT.TH1F("hcal_eta", "HCAL cluster eta;#eta;Entries", 100, -5, 5)
hcal_phi = ROOT.TH1F("hcal_phi", "HCAL cluster phi;#phi;Entries", 100, -3.5, 3.5)

#HBHE Rechits
hbhe_rechit_time = ROOT.TH1F("hbhe_rechit_time","HBHE Rechits time",100,-10,15)
hbhe_rechit_tdc  = ROOT.TH1F("hbhe_rechit_tdc", "HBHE Rechits TDC",128,0,64)
hbhe_rechit_clusteridx = ROOT.TH1F("hbheRechit_clusterIdx", "HBHE Rechits cluster index",100,0,100)
hbhe_clusteridx = ROOT.TH1F("clusterIdx", "HCAL Clusters index",100,0,100)
hbhe_rechit_counts = ROOT.TH1F("hbhe_rechit_counts", "HBHE rechit count per cluster", 100, 0, 1000)
hbhe_rechit_energy = ROOT.TH1F("hbhe_rechit_energy", "HBHE Rechits energy;E [GeV];Entries", 100, 0, 50)
hbhe_rechit_depth = ROOT.TH1F("hbhe_rechit_depth", "HBHE Rechits Depth",80,0,8)

hb_rechit_tdc = ROOT.TH1F("hb_rechit_tdc", "HB Rechits TDC",4,0,4)
he_rechit_tdc = ROOT.TH1F("he_rechit_tdc", "HE Rechits TDC",64,0,64)
hb_rechit_depth = ROOT.TH1F("hb_rechit_depth", "HB Rechits Depth",80,0,8)
he_rechit_depth = ROOT.TH1F("he_rechit_depth", "HE Rechits Depth",80,0,8)
hb_rechit_counts = ROOT.TH1F("hb_rechit_counts", "HB rechit count per cluster", 100, 0, 1200)
he_rechit_counts = ROOT.TH1F("he_rechit_counts", "HE rechit count per cluster", 100, 0, 1200)
hb_rechit_energy = ROOT.TH1F("hb_rechit_energy", "HB Rechits energy;E [GeV];Entries", 100, 0, 50)
he_rechit_energy = ROOT.TH1F("he_rechit_energy", "HE Rechits energy;E [GeV];Entries", 100, 0, 50)

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
        ecal_energy.Fill(e)
    for eta in event.ecal_eta:
        ecal_eta.Fill(eta)
    for phi in event.ecal_phi:
        ecal_phi.Fill(phi)
    for t in event.ecal_time:
        ecal_time.Fill(t)
    for idx in event.ecal_clusterIdx:
        ecal_clusterIdx.Fill(idx)

    #EB rechits
    for e in event.eb_rechit_energy:
        eb_rechit_energy.Fill(e)
    for eta in event.eb_rechit_eta:
        eb_rechit_eta.Fill(eta)
    for phi in event.eb_rechit_phi:
        eb_rechit_phi.Fill(phi)
    for time in event.eb_rechit_time:
        eb_rechit_time.Fill(time)
    for idx in event.eb_rechit_clusterIdx:
        eb_rechit_clusterIdx.Fill(idx)
    for count in event.eb_rechit_counts:
        eb_rechit_counts.Fill(count)
    #EE rechits
    for e in event.ee_rechit_energy:
        ee_rechit_energy.Fill(e)
    for eta in event.ee_rechit_eta:
        ee_rechit_eta.Fill(eta)
    for phi in event.ee_rechit_phi:
        ee_rechit_phi.Fill(phi)
    for time in event.ee_rechit_time:
        ee_rechit_time.Fill(time)
    for idx in event.ee_rechit_clusterIdx:
        ee_rechit_clusterIdx.Fill(idx)
    for count in event.ee_rechit_counts:
        ee_rechit_counts.Fill(count)

    # HCAL clusters
    for e in event.hcal_energy:
        #to-do: energy cuts, then loop through rechits in each cluster
        h_hcal_E.Fill(e)

    for eta in event.hcal_eta:  
        hcal_eta.Fill(eta)

    for phi in event.hcal_phi:
        hcal_phi.Fill(phi)
          
    for depth in event.hcal_depth:
        hcal_depth.Fill(depth)

    for t in event.hbhe_rechit_time:
        hbhe_rechit_time.Fill(t)

    for tdc in event.hbhe_rechit_tdc:
        hbhe_rechit_tdc.Fill(tdc)

    for idx in event.hbheRechit_clusterIdx:
        hbhe_rechit_clusteridx.Fill(idx)

    for idx in event.clusterIdx:
        hbhe_clusteridx.Fill(idx)

    for tdc in event.hb_rechit_tdc:
        hb_rechit_tdc.Fill(idx)
    
    for tdc in event.he_rechit_tdc:
        he_rechit_tdc.Fill(tdc)
    
    for depth in event.hbhe_rechit_depth:
        hbhe_rechit_depth.Fill(depth)

    for depth in event.hb_rechit_depth:
        hb_rechit_depth.Fill(depth)
    
    for depth in event.he_rechit_depth:
        he_rechit_depth.Fill(depth)

    for count in event.hb_rechit_counts:
        hb_rechit_counts.Fill(count)    

    for count in event.he_rechit_counts:
        he_rechit_counts.Fill(count)

    for count in event.hbhe_rechit_counts:
        hbhe_rechit_counts.Fill(count)
    
    for e in event.hb_rechit_energy:
        hb_rechit_energy.Fill(e)
    for e in event.he_rechit_energy:
        he_rechit_energy.Fill(e)
    for e in event.hbhe_rechit_energy:
        hbhe_rechit_energy.Fill(e)
        
    # Event#, Run #, lumi section
    event_number.Fill(event.event)
    run_number.Fill(event.run)
    lumi_section.Fill(event.lumi)

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
