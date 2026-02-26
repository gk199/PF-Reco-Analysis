import ROOT
import math

def delta_phi(phi1, phi2):
    dphi = phi1 - phi2
    while dphi > math.pi:
        dphi -= 2.0*math.pi
    while dphi <= -math.pi:
        dphi += 2.0*math.pi
    return dphi

def delta_r(eta1, phi1, eta2, phi2):
    deta = eta1 - eta2
    dphi = delta_phi(phi1, phi2)
    return math.sqrt(deta*deta + dphi*dphi)

# Load the ROOT file and TTree
#file = ROOT.TFile.Open("SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root")
inputfile_list =["ttbar/pfObjectsNtuple_new.root", "SinglePiPt100_1000_n1000/pfObjectsNtuple_new.root", "SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root"]
outputfile_list = ["ttbar/ttbar_hcalrh_notdccuts.root", "SinglePiPt100_1000_n1000/PiPt_hcalrh_notdccuts.root", "SinglePi0E100_1000_n5000/Pi0E_hcalrh_notdccuts.root"]

for inputfile, outputfile in zip(inputfile_list, outputfile_list):

    file = ROOT.TFile.Open(inputfile)
    tree = file.Get("pfObjectsNtupler/pfTree")

    if not tree:
        print("Error: Could not find TTree 'pfObjectsNtupler/pfTree' in file!")
        exit(1)

    # Optional: check entries
    print(f"Entries in tree: {tree.GetEntries()}")

    # Create output file for histograms
    out = ROOT.TFile(outputfile, "RECREATE")

    #Cuts applied: cluster energy > 30GeV and rechits with pf cuts

    # Define histograms
    #0) Seeds and rechit cuts
    root_seed_rh_diff_eta = ROOT.TH1F("root_seed_rh_diff_eta", "Difference in Eta Between Highest Energy Rechit and Seed;Delta Eta;Entries", 25, -0.4, 0.4)
    root_seed_rh_diff_phi = ROOT.TH1F("root_seed_rh_diff_phi", "Difference in Phi Between Highest Energy Rechit and Seed;Delta Phi;Entries", 25, -0.4, 0.4)
    root_seed_rh_diff_depth = ROOT.TH1F("root_seed_rh_diff_depth", "Difference in Depth Between Highest Energy Rechit and Seed;Delta Depth;Entries", 20, -2, 2)
    root_he_seed_rh_delta_time = ROOT.TH1F("root_he_seed_rh_delta_time", "Delta TDC between 'Seed' and other Rechits in HE; Delta TDC;Entries", 8, 0, 4)
    root_hb_seed_rh_delta_time= ROOT.TH1F("root_hb_seed_rh_delta_time", "Delta TDC between 'Seed' and other Rechits in HB; Delta TDC;Entries", 64, 0, 64)

    #1) HBHE rechits
    # Note: maxenergy_tdc_percluster only includes clusters with at least 2 rechits
    root_avg_hbhe_energy_percluster = ROOT.TH1F("root_avg_hbhe_energy_percluster", "Average HBHE Rechit Energy per HCAL Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_hbhe_tdc_percluster = ROOT.TH1F("root_avg_hbhe_tdc_percluster", "Average TDC per HBHE Cluster;TDC;Entries", 64, 0, 64)
    root_avg_hbhe_maxenergy_tdc_percluster = ROOT.TH1F("root_avg_hbhe_maxenergy_tdc_percluster", "Average TDC of 2-3 Highest Energy HBHE Rechits per HCAL Cluster; TDC;Entries", 64, 0, 64)

    root_avg_hbhe_tdc_percluster_narrow = ROOT.TH1F("root_avg_hbhe_tdc_percluster_narrow", "Average TDC per HBHE Cluster (narrow cone);TDC;Entries", 64, 0, 64)
    root_avg_hbhe_tdc_percluster_wide = ROOT.TH1F("root_avg_hbhe_tdc_percluster_wide", "Average TDC per HBHE Cluster (wide cone);TDC;Entries", 64, 0, 64)
    root_number_hbhe_percluster = ROOT.TH1F("root_number_hbhe_percluster", "Number of HBHE Rechits per HCAL Cluster; Number of Rechits;Entries", 40, 0, 80) 

    # Note: only cluster with nonzero rechit energy are considered for fraction histograms
    root_hbhe_frac_depth1_energy = ROOT.TH1F("root_hbhe_frac_depth1_energy", "Fraction of Depth 1 HBHE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth2_energy = ROOT.TH1F("root_hbhe_frac_depth2_energy", "Fraction of Depth 2 HBHE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth3_energy = ROOT.TH1F("root_hbhe_frac_depth3_energy", "Fraction of Depth 3 HBHE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth4_energy = ROOT.TH1F("root_hbhe_frac_depth4_energy", "Fraction of Depth 4 HBHE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1) 

    #2) HB Rechits
    root_avg_hb_energy_percluster = ROOT.TH1F("root_avg_hb_energy_percluster", "Average HB Rechit Energy per HCAL Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_hb_tdc_percluster = ROOT.TH1F("root_avg_hb_tdc_percluster", "Average TDC per HCAL Cluster;TDC;Entries", 40, 0, 4)
    root_avg_hb_maxenergy_tdc_percluster = ROOT.TH1F("root_avg_hb_maxenergy_tdc_percluster", "Average TDC of 2-3 Highest Energy HB Rechits per HCAL Cluster; TDC;Entries", 40, 0, 4)
    root_number_hb_percluster = ROOT.TH1F("root_number_hb_percluster", "Number of HB Rechits per HCAL Cluster; Number of Rechits;Entries", 40, 0, 80) 
    root_hb_frac_depth1_energy = ROOT.TH1F("root_hb_frac_depth1_energy", "Fraction of Depth 1 HB Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth2_energy = ROOT.TH1F("root_hb_frac_depth2_energy", "Fraction of Depth 2 HB Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth3_energy = ROOT.TH1F("root_hb_frac_depth3_energy", "Fraction of Depth 3 HB Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth4_energy = ROOT.TH1F("root_hb_frac_depth4_energy", "Fraction of Depth 4 HB Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1) 

    root_avg_hb_tdc_percluster_narrow = ROOT.TH1F("root_avg_hb_tdc_percluster_narrow", "Average TDC per HB Cluster (narrow cone);TDC;Entries", 40, 0, 4)
    root_avg_hb_tdc_percluster_wide = ROOT.TH1F("root_avg_hb_tdc_percluster_wide", "Average TDC per HB Cluster (wide cone);TDC;Entries", 40, 0, 4)

    #3) HE Rechits
    root_avg_he_energy_percluster = ROOT.TH1F("root_avg_he_energy_percluster", "Average HE Rechit Energy per HCAL Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_he_tdc_percluster = ROOT.TH1F("root_avg_he_tdc_percluster", "Average TDC per HCAL Cluster;TDC;Entries", 64, 0, 64)
    root_avg_he_maxenergy_tdc_percluster = ROOT.TH1F("root_avg_he_maxenergy_tdc_percluster", "Average TDC of 2-3 Highest Energy HE Rechits per HCAL Cluster; TDC;Entries", 64, 0, 64)
    root_number_he_percluster = ROOT.TH1F("root_number_he_percluster", "Number of HE Rechits per HCAL Cluster; Number of Rechits;Entries", 40, 0, 80) 
    root_he_frac_depth1_energy = ROOT.TH1F("root_he_frac_depth1_energy", "Fraction of Depth 1 HE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth2_energy = ROOT.TH1F("root_he_frac_depth2_energy", "Fraction of Depth 2 HE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth3_energy = ROOT.TH1F("root_he_frac_depth3_energy", "Fraction of Depth 3 HE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth4_energy = ROOT.TH1F("root_he_frac_depth4_energy", "Fraction of Depth 4 HE Rechit Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1) 

    root_avg_he_tdc_percluster_narrow = ROOT.TH1F("root_avg_he_tdc_percluster_narrow", "Average TDC per HE Cluster (narrow cone);TDC;Entries", 64, 0, 64)
    root_avg_he_tdc_percluster_wide = ROOT.TH1F("root_avg_he_tdc_percluster_wide", "Average TDC per HE Cluster (wide cone);TDC;Entries", 64, 0, 64)

    tree.Print()
    # Loop over entries
    for event in tree:

        print(f"\nProcessing new event")

        nHcal = len(event.hcal_energy)
        nHbHeRechits = len(event.hbhe_rechit_energy)
        nHbRechits = len(event.hb_rechit_energy)
        nHeRechits = len(event.he_rechit_energy)

        print(f"Number of HBHE rechits in this event: {nHbHeRechits}")
        print(f"Number of HCAL clusters in this event: {nHcal}")
        print(f"Number of HB rechits in this event: {nHbRechits}")
        print(f"Number of HE rechits in this event: {nHeRechits}")

        hcal_cluster_passed = 0
        hbhe_rechit_passed = 0
        hb_rechit_passed = 0
        he_rechit_passed = 0

        # Loop over HCAL clusters to apply cuts for both hbhe clusters
        for cluster_idx in range(nHcal):

            if event.hcal_energy[cluster_idx] > 30:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                hcal_cluster_passed += 1
                
                hbhe_rechit_passed_percluster = 0
                hbhe_rechit_passed_percluster_wide = 0
                hbhe_rechit_passed_percluster_narrow = 0

                total_hbhe_rechit_energy_percluster = 0
                hbhe_rechit_energy_tdc_percluster_pairs = []

                total_hbhe_tdc_percluster = 0
                total_hbhe_tdc_percluster_narrow = 0
                total_hbhe_tdc_percluster_wide = 0
                
                hbhe_rechit_depth1_energy_percluster = 0
                hbhe_rechit_depth2_energy_percluster = 0
                hbhe_rechit_depth3_energy_percluster = 0
                hbhe_rechit_depth4_energy_percluster = 0

                rh_list = [] # dictionary to store rh info for seed comparison
                #Cluster seed info from ntuple
                seed_eta, seed_phi, seed_depth =  event.hcal_seed_eta[cluster_idx], event.hcal_seed_phi[cluster_idx], event.hcal_seed_depth[cluster_idx]
            
                for i in range(nHbHeRechits): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.hbheRechit_clusterIdx[i] == cluster_idx:
                        if ((abs(event.hcal_eta[cluster_idx]) < 1.26) and \
                            ((event.hbhe_rechit_depth[i] == 1 and event.hbhe_rechit_energy[i] > 0.6) or \
                            (event.hbhe_rechit_depth[i] == 2 and event.hbhe_rechit_energy[i] > 0.4) or \
                            (event.hbhe_rechit_depth[i] == 3 and event.hbhe_rechit_energy[i] > 0.4) or \
                            (event.hbhe_rechit_depth[i] == 4 and event.hbhe_rechit_energy[i] > 0.5))) or \
                            (abs(event.hcal_eta[cluster_idx]) >= 1.26 and \
                            ((event.hbhe_rechit_depth[i] == 1 and event.hbhe_rechit_energy[i] > 0.2) or \
                            (event.hbhe_rechit_depth[i] == 2 and event.hbhe_rechit_energy[i] > 0.3) or \
                            (event.hbhe_rechit_depth[i] == 3 and event.hbhe_rechit_energy[i] > 0.3) or \
                            (event.hbhe_rechit_depth[i] == 4 and event.hbhe_rechit_energy[i] > 0.3) or \
                            (event.hbhe_rechit_depth[i] == 5 and event.hbhe_rechit_energy[i] > 0.3) or \
                            (event.hbhe_rechit_depth[i] == 6 and event.hbhe_rechit_energy[i] > 0.3) or \
                            (event.hbhe_rechit_depth[i] == 7 and event.hbhe_rechit_energy[i] > 0.3))):
                            
                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]

                            rh_eta = event.hbhe_rechit_eta[i]
                            rh_phi = event.hbhe_rechit_phi[i]

                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                hbhe_rechit_passed_percluster_wide += 1
                                total_hbhe_tdc_percluster_wide += event.hbhe_rechit_tdc[i]
                            if dr <= 0.2:
                                hbhe_rechit_passed_percluster_narrow += 1
                                total_hbhe_tdc_percluster_narrow += event.hbhe_rechit_tdc[i]

                            hbhe_rechit_passed_percluster += 1
                            hbhe_rechit_passed += 1
                            total_hbhe_rechit_energy_percluster += event.hbhe_rechit_energy[i]
                            total_hbhe_tdc_percluster += event.hbhe_rechit_tdc[i]
                            
                            hbhe_rechit_energy_tdc_percluster_pairs.append((event.hbhe_rechit_energy[i], event.hbhe_rechit_tdc[i]))
                            rh_list.append((event.hbhe_rechit_energy[i], event.hbhe_rechit_tdc[i], rh_eta, rh_phi, event.hbhe_rechit_depth[i]))
                            #Fill per depth energy sums
                            if event.hbhe_rechit_depth[i] == 1:
                                hbhe_rechit_depth1_energy_percluster += event.hbhe_rechit_energy[i]
                            elif event.hbhe_rechit_depth[i] == 2:
                                hbhe_rechit_depth2_energy_percluster += event.hbhe_rechit_energy[i]
                            elif event.hbhe_rechit_depth[i] == 3:
                                    hbhe_rechit_depth3_energy_percluster += event.hbhe_rechit_energy[i]
                            elif event.hbhe_rechit_depth[i] == 4:
                                    hbhe_rechit_depth4_energy_percluster += event.hbhe_rechit_energy[i]
                
                #Store seed & pfrh info for clusters with at least 1 pfrh passing cuts
                if len(rh_list) > 0:
                    # pick the highest-energy PFRecHit in this cluster
                    rh_seed_idx = max(range(len(rh_list)), key=lambda j: rh_list[j][0])
                    rh_seed_energy, rh_seed_time, rh_seed_eta, rh_seed_phi, rh_seed_depth = rh_list[rh_seed_idx]

                    # compare PFRecHit seed vs cluster seed
                    root_seed_rh_diff_depth.Fill(rh_seed_depth - seed_depth)
                    root_seed_rh_diff_eta.Fill(rh_seed_eta - seed_eta)
                    root_seed_rh_diff_phi.Fill(rh_seed_phi - seed_phi)

                    # fill delta time for non-seed PFRecHits in the same cluster
                    # remove all invalid time values (-999) (in this section only!) and the seed itself from the delta time histogram
                    if abs(event.hcal_eta[cluster_idx]) < 1.26: #clusters in HB
                        for j, (rh_energy, rh_time, rh_eta, rh_phi, rh_depth) in enumerate(rh_list):
                            if j == rh_seed_idx:
                                continue
                            root_hb_seed_rh_delta_time.Fill(rh_time - rh_seed_time)

                    elif abs(event.hcal_eta[cluster_idx]) >= 1.26: #clusters in HE
                        for j, (rh_energy, rh_time, rh_eta, rh_phi, rh_depth) in enumerate(rh_list):
                            if j == rh_seed_idx:
                                continue
                            root_he_seed_rh_delta_time.Fill(rh_time - rh_seed_time)
                
                root_number_hbhe_percluster.Fill(hbhe_rechit_passed_percluster)

                if hbhe_rechit_passed_percluster > 0:
                    #if hbhe_rechit_passed_percluster > 0 and total_hbhe_maxenergy_tdc_percluster > 0:
                    #Fill in avg number of rechits, avg energy of rechits, avg tdc per cluster                
                    root_avg_hbhe_energy_percluster.Fill(total_hbhe_rechit_energy_percluster / hbhe_rechit_passed_percluster)
                    root_avg_hbhe_tdc_percluster.Fill(total_hbhe_tdc_percluster / hbhe_rechit_passed_percluster) 
                    
                    if hbhe_rechit_passed_percluster_wide > 0:
                        root_avg_hbhe_tdc_percluster_wide.Fill(total_hbhe_tdc_percluster_wide / hbhe_rechit_passed_percluster_wide)
                    
                    if  hbhe_rechit_passed_percluster_narrow > 0:
                        root_avg_hbhe_tdc_percluster_narrow.Fill(total_hbhe_tdc_percluster_narrow / hbhe_rechit_passed_percluster_narrow)
                    
                    #For clusters with nonzero total rechits energy(!), fill fraction of energy per depth
                    if total_hbhe_rechit_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_hbhe_frac_depth1_energy.Fill(hbhe_rechit_depth1_energy_percluster / total_hbhe_rechit_energy_percluster)
                        root_hbhe_frac_depth2_energy.Fill(hbhe_rechit_depth2_energy_percluster / total_hbhe_rechit_energy_percluster)
                        root_hbhe_frac_depth3_energy.Fill(hbhe_rechit_depth3_energy_percluster / total_hbhe_rechit_energy_percluster)
                        root_hbhe_frac_depth4_energy.Fill(hbhe_rechit_depth4_energy_percluster / total_hbhe_rechit_energy_percluster)

                    #For clusters with 2-3+ rechits, find the 2-3 rechits with max energy and fill their corresponding tdc
                    if len(hbhe_rechit_energy_tdc_percluster_pairs) >= 3:
                        top = sorted(hbhe_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_hbhe_maxenergy_tdc_percluster.Fill(avg_tdc) 

                    elif len(hbhe_rechit_energy_tdc_percluster_pairs) == 2:
                        top = sorted(hbhe_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_hbhe_maxenergy_tdc_percluster.Fill(avg_tdc)

                    
            
            # Now process HB clusters   
            if event.hcal_energy[cluster_idx] > 30 and abs(event.hcal_eta[cluster_idx]) < 1.26:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                hb_rechit_passed_percluster = 0
                total_hb_rechit_energy_percluster = 0
                hb_rechit_energy_tdc_percluster_pairs = []
                total_hb_tdc_percluster = 0

                hb_rechit_passed_percluster_wide = 0
                total_hb_tdc_percluster_wide = 0
                hb_rechit_passed_percluster_narrow = 0 
                total_hb_tdc_percluster_narrow = 0

                total_hb_rechit_energy_percluster = 0
                hb_rechit_depth1_energy_percluster = 0
                hb_rechit_depth2_energy_percluster = 0
                hb_rechit_depth3_energy_percluster = 0
                hb_rechit_depth4_energy_percluster = 0

                for i in range(nHbRechits): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.hb_rechit_clusterIdx[i] == cluster_idx:
                        if ((event.hb_rechit_depth[i] == 1 and event.hb_rechit_energy[i] > 0.6) or \
                            (event.hb_rechit_depth[i] == 2 and event.hb_rechit_energy[i] > 0.4) or \
                            (event.hb_rechit_depth[i] == 3 and event.hb_rechit_energy[i] > 0.4) or \
                            (event.hb_rechit_depth[i] == 4 and event.hb_rechit_energy[i] > 0.5)):
                            
                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]
                            rh_eta = event.hb_rechit_eta[i]
                            rh_phi = event.hb_rechit_phi[i]

                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                hb_rechit_passed_percluster_wide += 1
                                total_hb_tdc_percluster_wide += event.hb_rechit_tdc[i]
                            if dr <= 0.2:
                                hb_rechit_passed_percluster_narrow += 1
                                total_hb_tdc_percluster_narrow += event.hb_rechit_tdc[i]


                            hb_rechit_passed_percluster += 1
                            hb_rechit_passed += 1
                            total_hb_rechit_energy_percluster += event.hb_rechit_energy[i]
                            total_hb_tdc_percluster += event.hb_rechit_tdc[i]
                            hb_rechit_energy_tdc_percluster_pairs.append((event.hb_rechit_energy[i], event.hb_rechit_tdc[i]))
                                
                            #Fill per depth energy sums
                            if event.hb_rechit_depth[i] == 1:
                                hb_rechit_depth1_energy_percluster += event.hb_rechit_energy[i]
                            elif event.hb_rechit_depth[i] == 2:
                                hb_rechit_depth2_energy_percluster += event.hb_rechit_energy[i]
                            elif event.hb_rechit_depth[i] == 3:
                                hb_rechit_depth3_energy_percluster += event.hb_rechit_energy[i]
                            elif event.hb_rechit_depth[i] == 4:
                                hb_rechit_depth4_energy_percluster += event.hb_rechit_energy[i]

                root_number_hb_percluster.Fill(hb_rechit_passed_percluster) 
                if hb_rechit_passed_percluster > 0:
                    #Fill in avg number of rechits, avg energy of rechits, avg tdc per cluster                
                    root_avg_hb_energy_percluster.Fill(total_hb_rechit_energy_percluster / hb_rechit_passed_percluster) 
                    root_avg_hb_tdc_percluster.Fill(total_hb_tdc_percluster / hb_rechit_passed_percluster)

                    if hb_rechit_passed_percluster_wide > 0:
                        root_avg_hb_tdc_percluster_wide.Fill(total_hb_tdc_percluster_wide / hb_rechit_passed_percluster_wide)
                    if hb_rechit_passed_percluster_narrow > 0:
                        root_avg_hb_tdc_percluster_narrow.Fill(total_hb_tdc_percluster_narrow / hb_rechit_passed_percluster_narrow)

                    #For clusters with nonzero total rechits energy(!), fill fraction of energy per depth
                    if total_hb_rechit_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_hb_frac_depth1_energy.Fill(hb_rechit_depth1_energy_percluster / total_hb_rechit_energy_percluster)
                        root_hb_frac_depth2_energy.Fill(hb_rechit_depth2_energy_percluster / total_hb_rechit_energy_percluster)
                        root_hb_frac_depth3_energy.Fill(hb_rechit_depth3_energy_percluster / total_hb_rechit_energy_percluster)
                        root_hb_frac_depth4_energy.Fill(hb_rechit_depth4_energy_percluster / total_hb_rechit_energy_percluster)

                    if len(hb_rechit_energy_tdc_percluster_pairs) >= 3:
                        top = sorted(hb_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_hb_maxenergy_tdc_percluster.Fill(avg_tdc) 

                    elif len(hb_rechit_energy_tdc_percluster_pairs) == 2:
                        top = sorted(hb_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_hb_maxenergy_tdc_percluster.Fill(avg_tdc)
                

            
            # Now process HE clusters   
            if event.hcal_energy[cluster_idx] > 30 and abs(event.hcal_eta[cluster_idx]) >= 1.26:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                he_rechit_passed_percluster = 0
                total_he_rechit_energy_percluster = 0
                he_rechit_energy_tdc_percluster_pairs = []
                total_he_tdc_percluster = 0
                total_he_maxenergy_tdc_percluster = 0

                he_rechit_passed_percluster_wide = 0
                total_he_tdc_percluster_wide = 0
                he_rechit_passed_percluster_narrow = 0
                total_he_tdc_percluster_narrow = 0

                total_he_rechit_energy_percluster = 0
                he_rechit_depth1_energy_percluster = 0
                he_rechit_depth2_energy_percluster = 0
                he_rechit_depth3_energy_percluster = 0
                he_rechit_depth4_energy_percluster = 0

                for i in range(nHeRechits): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.he_rechit_clusterIdx[i] == cluster_idx:

                        if ((event.he_rechit_depth[i] == 1 and event.he_rechit_energy[i] > 0.2) or \
                            (event.he_rechit_depth[i] == 2 and event.he_rechit_energy[i] > 0.3) or \
                            (event.he_rechit_depth[i] == 3 and event.he_rechit_energy[i] > 0.3) or \
                            (event.he_rechit_depth[i] == 4 and event.he_rechit_energy[i] > 0.3) or \
                            (event.he_rechit_depth[i] == 5 and event.he_rechit_energy[i] > 0.3) or \
                            (event.he_rechit_depth[i] == 6 and event.he_rechit_energy[i] > 0.3) or \
                            (event.he_rechit_depth[i] == 7 and event.he_rechit_energy[i] > 0.3)):
                            
                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]
                            rh_eta = event.he_rechit_eta[i]
                            rh_phi = event.he_rechit_phi[i] 
                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                he_rechit_passed_percluster_wide += 1
                                total_he_tdc_percluster_wide += event.he_rechit_tdc[i]
                            if dr <= 0.2:   
                                he_rechit_passed_percluster_narrow += 1
                                total_he_tdc_percluster_narrow += event.he_rechit_tdc[i]

                            he_rechit_passed_percluster += 1
                            he_rechit_passed += 1
                            total_he_rechit_energy_percluster += event.he_rechit_energy[i]
                            total_he_tdc_percluster += event.he_rechit_tdc[i]
                            he_rechit_energy_tdc_percluster_pairs.append((event.he_rechit_energy[i], event.he_rechit_tdc[i]))

                            #Fill per depth energy sums
                            if event.he_rechit_depth[i] == 1:
                                he_rechit_depth1_energy_percluster += event.he_rechit_energy[i]
                            elif event.he_rechit_depth[i] == 2:
                                he_rechit_depth2_energy_percluster += event.he_rechit_energy[i]
                            elif event.he_rechit_depth[i] == 3:
                                he_rechit_depth3_energy_percluster += event.he_rechit_energy[i]
                            elif event.he_rechit_depth[i] == 4:
                                he_rechit_depth4_energy_percluster += event.he_rechit_energy[i]
                
                root_number_he_percluster.Fill(he_rechit_passed_percluster)

                if he_rechit_passed_percluster > 0:
                    #Fill in avg number of rechits, avg energy of rechits, avg tdc per cluster                
                    root_avg_he_energy_percluster.Fill(total_he_rechit_energy_percluster / he_rechit_passed_percluster)
                    root_avg_he_tdc_percluster.Fill(total_he_tdc_percluster / he_rechit_passed_percluster) 

                    if he_rechit_passed_percluster_wide > 0:
                        root_avg_he_tdc_percluster_wide.Fill(total_he_tdc_percluster_wide / he_rechit_passed_percluster_wide)
                    if he_rechit_passed_percluster_narrow > 0:
                        root_avg_he_tdc_percluster_narrow.Fill(total_he_tdc_percluster_narrow / he_rechit_passed_percluster_narrow)

                    #For clusters with nonzero total rechits energy(!), fill fraction of energy per depth
                    if total_he_rechit_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_he_frac_depth1_energy.Fill(he_rechit_depth1_energy_percluster / total_he_rechit_energy_percluster)
                        root_he_frac_depth2_energy.Fill(he_rechit_depth2_energy_percluster / total_he_rechit_energy_percluster)
                        root_he_frac_depth3_energy.Fill(he_rechit_depth3_energy_percluster / total_he_rechit_energy_percluster)
                        root_he_frac_depth4_energy.Fill(he_rechit_depth4_energy_percluster / total_he_rechit_energy_percluster)
                    
                    if len(he_rechit_energy_tdc_percluster_pairs) >= 3:
                        top = sorted(he_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_he_maxenergy_tdc_percluster.Fill(avg_tdc)   
                    elif len(he_rechit_energy_tdc_percluster_pairs) == 2:
                        top = sorted(he_rechit_energy_tdc_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_tdc = sum(tdc for _, tdc in top) / len(top)
                        root_avg_he_maxenergy_tdc_percluster.Fill(avg_tdc)
    #Note:

        print(f"Number of HCAL clusters passing cuts: {hcal_cluster_passed}")
        print(f"Number of HBHE rechits passing cuts: {hbhe_rechit_passed}")
        print(f"Number of HB rechits passing cuts: {hb_rechit_passed}")
        print(f"Number of HE rechits passing cuts: {he_rechit_passed}")

    # Save histograms
    out.Write()
    out.Close()

    print("Histograms saved to", outputfile)
