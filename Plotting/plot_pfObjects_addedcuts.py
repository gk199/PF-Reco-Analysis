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
out = ROOT.TFile("pfObjectsHistos_addedcuts.root", "RECREATE")

# Define histograms
eb_rechit_time   = ROOT.TH1F("eb_rechit_time",   "ECAL Barrel, Rechit Time;t [ns];Entries", 50, -25, 25)
eb_rechit_counts  = ROOT.TH1F("eb_rechit_counts", "ECAL Barrel, # of Rechits Per Cluster (After Cuts)", 15, 15, 45)
eb_rechit_energy = ROOT.TH1F("eb_rechit_energy", "ECAL Barrel, Rechit Energy;E [GeV];Entries", 50, 0, 25)

hb_rechit_tdc = ROOT.TH1F("hb_rechit_tdc", "HCAL Barrel Rechits TDC",4,0,4)
hb_rechit_depth = ROOT.TH1F("hb_rechit_depth", "HCAL Barrel Rechits Depth; Depth Layer",7,0,7)
hb_rechit_counts = ROOT.TH1F("hb_rechit_counts", "HCAL Barrel Rechits Counts Per Cluster (After Cuts)", 50, 0, 50)
hb_energy = ROOT.TH1F("hb_energy", "HCAL Barrel Rechits Energy;E [GeV];Entries", 40, 0, 20)
#to-do: do print statements to double check: nEcal, nEcalRechits, nHcal, nHbRechits and check how many passed the cuts

# Loop over entries
for event in tree:
    print(f"\nProcessing new event")
    nEcal = len(event.ecal_energy)
    print(f"Number of ECAL clusters: {nEcal}")
    nEcalRechits = len(event.eb_rechit_energy)
    print(f"Number of EB rechits: {nEcalRechits}")

    eb_rechit_passed = 0
    ecal_cluster_passed = 0

    for cluster_idx in range(nEcal):
        if 8.0 < event.ecal_energy[cluster_idx] < 30.0 and abs(event.ecal_eta[cluster_idx]) < 1.26:
            eb_rechit_passed_percluster = 0 
            for i in range(nEcalRechits):
                if event.eb_rechit_clusterIdx[i] == cluster_idx:
                    eb_rechit_energy.Fill(event.eb_rechit_energy[i])
                    eb_rechit_time.Fill(event.eb_rechit_time[i])
                    eb_rechit_passed_percluster += 1
                    eb_rechit_passed += 1
            eb_rechit_counts.Fill(eb_rechit_passed_percluster)

            ecal_cluster_passed += 1
    print(f"Number of ECAL clusters passing cuts: {ecal_cluster_passed}")
    print(f"Number of EB rechits passing cuts: {eb_rechit_passed}")

    # Loop over all EB rechits:
    """print(f"Number of EB rechit in this event: {nEcalRechits}")
    eb_rechit_passed = 0

    for i in range(nEcalRechits):
        
        #Find the cluster index this EB rechit belongs to
        cluster_idx = event.eb_rechit_clusterIdx[i]

        #Find the energy and eta of the corresponding ECAL cluster, then apply cuts
        if 10.0 < event.ecal_energy[cluster_idx] < 20.0 and abs(event.ecal_eta[cluster_idx]) < 1.26:
            #If passes cuts, fill histograms of this rechit
            eb_rechit_energy.Fill(event.eb_rechit_energy[i])
            eb_rechit_time.Fill(event.eb_rechit_time[i])
            eb_rechit_passed += 1"""

    
    # Loop over all HCAL clusters to fill the rechit counts per cluster histogram
    nHcal = len(event.hcal_energy)
    nHbRechits = len(event.hb_rechit_energy)
    print(f"Number of HB rechits in this event: {nHbRechits}")
    hb_rechit_passed = 0
    print(f"Number of HCAL clusters in this event: {nHcal}")
    hcal_cluster_passed = 0

    for cluster_idx in range(nHcal):
        if 8.0 < event.hcal_energy[cluster_idx] < 30.0 and abs(event.hcal_eta[cluster_idx]) < 1.26:
            #hb_rechit_counts.Fill(event.hb_rechit_counts[cluster_idx])
            hb_rechit_passed_percluster = 0
            for i in range(nHbRechits):
                if event.hb_rechit_clusterIdx[i] == cluster_idx:
                    hb_energy.Fill(event.hb_rechit_energy[i])
                    hb_rechit_depth.Fill(event.hb_rechit_depth[i])
                    hb_rechit_tdc.Fill(event.hb_rechit_tdc[i])
                    hb_rechit_passed_percluster += 1
                    hb_rechit_passed += 1
            hb_rechit_counts.Fill(hb_rechit_passed_percluster)
            hcal_cluster_passed += 1
    print(f"Number of HCAL clusters passing cuts: {hcal_cluster_passed}")
    print(f"Number of HB rechits passing cuts: {hb_rechit_passed}")

    """# Loop over all HB rechits:
    for i in range(nHbRechits):

        #Find the cluster index this HB rechit belongs to
        cluster_idx = event.hb_rechit_clusterIdx[i] 

        #Find the energy and eta of the corresponding HCAL cluster, then apply cuts 
        if 10.0 < event.hcal_energy[cluster_idx] < 20.0 and abs(event.hcal_eta[cluster_idx]) < 1.26:
            #If passes cuts, fill histograms of this rechit
            hb_rechit_tdc.Fill(event.hb_rechit_tdc[i])
            hb_rechit_depth.Fill(event.hb_rechit_depth[i])
            hb_energy.Fill(event.hb_rechit_energy[i])
            hb_rechit_passed += 1"""
    

# Write histograms to output file

# Optionally, create and save plots with CMS header

# Write histograms to the ROOT file
out.cd()
eb_rechit_time.Write()
eb_rechit_counts.Write()
eb_rechit_energy.Write()

hb_rechit_tdc.Write()
hb_rechit_depth.Write()
hb_rechit_counts.Write()
hb_energy.Write()

out.Close()
file.Close()
print("Saved histos + PNGs")

