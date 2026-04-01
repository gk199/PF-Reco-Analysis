import ROOT
import math
import numpy as np

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

#file = ROOT.TFile.Open("SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root")
#inputfile_list =["ttbar/pfObjectsNtuple_new.root", "SinglePiPt100_1000_n1000/pfObjectsNtuple_new.root", "SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root"]
#outputfile_list = ["ttbar/ttbar_pfrh_notimecuts.root", "SinglePiPt100_1000_n1000/PiPt_pfrh_notimecuts.root", "SinglePi0E100_1000_n5000/Pi0E_pfrh_notimecuts.root"]
eos_path  = "/eos/user/c/chtong/Public/Rereco/"
inputfile_list =[eos_path + "ttbar_rereco/pfObjectsNtuple_new.root", 
                 eos_path + "SinglePiPt100_1000_n1000/pfObjectsNtuple_new.root", 
                 eos_path + "SinglePi0E100_1000_n5000/pfObjectsNtuple_new.root",
                 eos_path + "ttbar_rereco/pfObjectsNtuple_new_timing.root", 
                 eos_path + "SinglePiPt100_1000_n1000/pfObjectsNtuple_new_timing.root", 
                 eos_path + "SinglePi0E100_1000_n5000/pfObjectsNtuple_new_timing.root"]

outputfile_list = [eos_path + "ttbar_rereco/ttbar_hcalrh_notimecuts.root", 
                   eos_path + "SinglePiPt100_1000_n1000/PiPt_hcalrh_notimecuts.root", 
                   eos_path + "SinglePi0E100_1000_n5000/Pi0E_hcalrh_notimecuts.root",
                   eos_path + "ttbar_rereco/ttbar_hcalrh_notimecuts_timing.root", 
                   eos_path + "SinglePiPt100_1000_n1000/PiPt_hcalrh_notimecuts_timing.root",
                   eos_path + "SinglePi0E100_1000_n5000/Pi0E_hcalrh_notimecuts_timing.root"]
for inputfile, outputfile in zip(inputfile_list, outputfile_list):
    # Load the ROOT file and TTree

    file = ROOT.TFile.Open(inputfile)
    tree = file.Get("pfObjectsNtupler/pfTree")

    if not tree:
        print("Error: Could not find TTree 'pfObjectsNtupler/pfTree' in file!")
        exit(1)

    # Optional: check entries
    print(f"Entries in tree: {tree.GetEntries()}")

    # Create output file for histograms
    out = ROOT.TFile(outputfile, "RECREATE")


    #Cuts applied: cluster energy > 30GeV

    # Define histograms
    #0) seeds and pfrh cuts
    root_seed_pfrh_diff_eta = ROOT.TH1F("root_seed_pfrh_diff_eta", "Difference in Eta Between Highest Energy PFrh and Seed;Delta Eta;Entries", 25, -0.4, 0.4)
    root_seed_pfrh_diff_phi = ROOT.TH1F("root_seed_pfrh_diff_phi", "Difference in Phi Between Highest Energy PFrh and Seed;Delta Phi;Entries", 25, -0.4, 0.4)
    root_seed_pfrh_diff_depth = ROOT.TH1F("root_seed_pfrh_diff_depth", "Difference in Depth Between Highest Energy PFrh and Seed;Delta Depth;Entries", 20, -2, 2)
    root_hbhe_seed_pfrh_delta_time = ROOT.TH1F("root_hbhe_seed_pfrh_delta_time", "Delta Time between 'Seed' and other PFrh;Delta Time;Entries", 80, -8, 8)
    root_delta_pfcuts = ROOT.TH1F("root_delta_pfcuts", "Difference between pfrh counts before and after pfcuts;PFrh Count;Entries", 40, -1, 8)
    #1) HBHE pfrh
    # Note: maxenergy_time_percluster only includes clusters with at least 2 pfrh
    root_avg_hbhe_energy_percluster = ROOT.TH1F("root_avg_pfhbhe_energy_percluster", "Average PFrh Energy per HBHE Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_hbhe_time_percluster = ROOT.TH1F("root_avg_pfhbhe_time_percluster", "Average Time per HBHE Cluster;Time [ns];Entries", 50, -10, 10)
    root_avg_hbhe_maxenergy_time_percluster = ROOT.TH1F("root_avg_pfhbhe_maxenergy_time_percluster", "Average Time of 2-3 Highest Energy PFrh per HBHE Cluster; Time;Entries",  50, -10, 10)
    root_avg_hbhe_time_percluster_narrow = ROOT.TH1F("root_avg_pfhbhe_time_percluster_narrow", "Average Time per HBHE Cluster (PFrh, narrow cone);Time ;Entries", 50, -10, 10)
    root_avg_hbhe_time_percluster_wide = ROOT.TH1F("root_avg_pfhbhe_time_percluster_wide", "Average Time per HBHE Cluster (PFrh,wide cone);Time ;Entries", 50, -10, 10)
    root_number_hbhe_percluster = ROOT.TH1F("root_number_pfhbhe_percluster", "Number of PFrh per HBHE Cluster; Number of PFrh;Entries", 40, 0, 80) 

    # Note: only cluster with nonzero pfrh energy are considered for fraction histograms
    root_hbhe_frac_depth1_energy = ROOT.TH1F("root_pfhbhe_frac_depth1_energy", "Fraction of Depth 1 HBHE PFrh Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth2_energy = ROOT.TH1F("root_pfhbhe_frac_depth2_energy", "Fraction of Depth 2 HBHE PFrh Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth3_energy = ROOT.TH1F("root_pfhbhe_frac_depth3_energy", "Fraction of Depth 3 HBHE PFrh Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1)
    root_hbhe_frac_depth4_energy = ROOT.TH1F("root_pfhbhe_frac_depth4_energy", "Fraction of Depth 4 HBHE PFrh Energy in HCAL Cluster;Fraction;Entries", 20, 0, 1) 

    #Delta time between seed and pfrh in cluster

    #2) HB pfrh
    root_avg_hb_energy_percluster = ROOT.TH1F("root_avg_pfhb_energy_percluster", "Average PFrh Energy per HB Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_hb_time_percluster = ROOT.TH1F("root_avg_pfhb_time_percluster", "Average Time per HB Cluster;Time [ns];Entries", 50, -10, 10)
    root_avg_hb_maxenergy_time_percluster = ROOT.TH1F("root_avg_pfhb_maxenergy_time_percluster", "Average Time of 2-3 Highest Energy PFrh per HB Cluster; Time;Entries", 50, -10, 10)
    root_number_hb_percluster = ROOT.TH1F("root_number_pfhb_percluster", "Number of PFrh per HB Cluster; Number of PFrh;Entries", 40, 0, 80) 
    root_hb_frac_depth1_energy = ROOT.TH1F("root_pfhb_frac_depth1_energy", "Fraction of Depth 1 PFrh Energy in HB Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth2_energy = ROOT.TH1F("root_pfhb_frac_depth2_energy", "Fraction of Depth 2 PFrh Energy in HB Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth3_energy = ROOT.TH1F("root_pfhb_frac_depth3_energy", "Fraction of Depth 3 PFrh Energy in HB Cluster;Fraction;Entries", 20, 0, 1)
    root_hb_frac_depth4_energy = ROOT.TH1F("root_pfhb_frac_depth4_energy", "Fraction of Depth 4 PFrh Energy in HB Cluster;Fraction;Entries", 20, 0, 1) 

    root_avg_hb_time_percluster_narrow = ROOT.TH1F("root_avg_pfhb_time_percluster_narrow", "Average Time per HB Cluster (PFrh, narrow cone);Time;Entries", 50, -10, 10)
    root_avg_hb_time_percluster_wide = ROOT.TH1F("root_avg_pfhb_time_percluster_wide", "Average Time per HB Cluster (PFrh, wide cone);Time;Entries", 50, -10, 10)

    #3) HE pfrh
    root_avg_he_energy_percluster = ROOT.TH1F("root_avg_pfhe_energy_percluster", "Average PFrh Energy per HE Cluster;E [GeV];Entries", 20, 0, 40)
    root_avg_he_time_percluster = ROOT.TH1F("root_avg_pfhe_time_percluster", "Average Time per HE Cluster;Time [ns];Entries", 50, -10, 10)
    root_avg_he_maxenergy_time_percluster = ROOT.TH1F("root_avg_pfhe_maxenergy_time_percluster", "Average Time of 2-3 Highest Energy PFrh per HE Cluster; Time;Entries", 50, -10, 10)
    root_number_he_percluster = ROOT.TH1F("root_number_pfhe_percluster", "Number of PFrh per HE Cluster; Number of PFrh;Entries", 40, 0, 80) 
    root_he_frac_depth1_energy = ROOT.TH1F("root_pfhe_frac_depth1_energy", "Fraction of Depth 1 PFrh Energy in HE Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth2_energy = ROOT.TH1F("root_pfhe_frac_depth2_energy", "Fraction of Depth 2 PFrh Energy in HE Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth3_energy = ROOT.TH1F("root_pfhe_frac_depth3_energy", "Fraction of Depth 3 PFrh Energy in HE Cluster;Fraction;Entries", 20, 0, 1)
    root_he_frac_depth4_energy = ROOT.TH1F("root_pfhe_frac_depth4_energy", "Fraction of Depth 4 PFrh Energy in HE Cluster;Fraction;Entries", 20, 0, 1) 
    root_avg_he_time_percluster_narrow = ROOT.TH1F("root_avg_pfhe_time_percluster_narrow", "Average Time per HE Cluster (PFrh, narrow cone);Time;Entries", 50, -10, 10)
    root_avg_he_time_percluster_wide = ROOT.TH1F("root_avg_pfhe_time_percluster_wide", "Average Time per HE Cluster (PFrh, wide cone);Time;Entries", 50, -10, 10)

    tree.Print()
    seedE_all = []
    seedT_all = []
    seedE_HB, seedT_HB = [], []
    seedE_HE, seedT_HE = [], []

    # Loop over entries
    for event in tree:

        print(f"\nProcessing new event")

        nHcal = len(event.hcal_energy)
        nHbHepfrh = len(event.hbhe_pfrh_energy)
        nHbpfrh = len(event.hb_pfrh_energy)
        nHepfrh = len(event.he_pfrh_energy)

        print(f"Number of HBHE PFrh in this event: {nHbHepfrh}")
        print(f"Number of HCAL clusters in this event: {nHcal}")
        print(f"Number of HB PFrh in this event: {nHbpfrh}")
        print(f"Number of HE PFrh in this event: {nHepfrh}")

        hcal_cluster_passed = 0
        hbhe_pfrh_passed = 0
        hb_pfrh_passed = 0
        he_pfrh_passed = 0

        # Loop over HCAL clusters to apply cuts for both hbhe clusters
        for cluster_idx in range(nHcal):

            if event.hcal_energy[cluster_idx] > 30:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                hcal_cluster_passed += 1
                
                hbhe_pfrh_passed_percluster = 0
                hbhe_pfrh_passed_percluster_wide = 0
                hbhe_pfrh_passed_percluster_narrow = 0

                hbhe_pfrh_passed_percluster_time = 0

                total_hbhe_pfrh_energy_percluster = 0
                hbhe_pfrh_energy_time_percluster_pairs = []

                total_hbhe_time_percluster = 0
                total_hbhe_time_percluster_narrow = 0
                total_hbhe_time_percluster_wide = 0
                
                hbhe_pfrh_depth1_energy_percluster = 0
                hbhe_pfrh_depth2_energy_percluster = 0
                hbhe_pfrh_depth3_energy_percluster = 0
                hbhe_pfrh_depth4_energy_percluster = 0

                rh_list = [] # dictionary to store rh info for seed comparison
                seed_eta, seed_phi, seed_depth =  event.hcal_seed_eta[cluster_idx], event.hcal_seed_phi[cluster_idx], event.hcal_seed_depth[cluster_idx]
            
                rh_count_b4_pfcuts = 0
                rh_count_after_pfcuts = 0
                
                for i in range(nHbHepfrh): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.hbhe_pfrh_clusterIdx[i] == cluster_idx:
                        rh_count_b4_pfcuts += 1
                        if ((abs(event.hcal_eta[cluster_idx]) < 1.26 ) and \
                            ((event.hbhe_pfrh_depth[i] == 1 and event.hbhe_pfrh_energy[i] > 0.6) or \
                            (event.hbhe_pfrh_depth[i] == 2 and event.hbhe_pfrh_energy[i] > 0.4) or \
                            (event.hbhe_pfrh_depth[i] == 3 and event.hbhe_pfrh_energy[i] > 0.4) or \
                            (event.hbhe_pfrh_depth[i] == 4 and event.hbhe_pfrh_energy[i] > 0.5))) or \
                            (abs(event.hcal_eta[cluster_idx]) >= 1.26 and \
                            ((event.hbhe_pfrh_depth[i] == 1 and event.hbhe_pfrh_energy[i] > 0.2) or \
                            (event.hbhe_pfrh_depth[i] == 2 and event.hbhe_pfrh_energy[i] > 0.3) or \
                            (event.hbhe_pfrh_depth[i] == 3 and event.hbhe_pfrh_energy[i] > 0.3) or \
                            (event.hbhe_pfrh_depth[i] == 4 and event.hbhe_pfrh_energy[i] > 0.3) or \
                            (event.hbhe_pfrh_depth[i] == 5 and event.hbhe_pfrh_energy[i] > 0.3) or \
                            (event.hbhe_pfrh_depth[i] == 6 and event.hbhe_pfrh_energy[i] > 0.3) or \
                            (event.hbhe_pfrh_depth[i] == 7 and event.hbhe_pfrh_energy[i] > 0.3))):

                            rh_count_after_pfcuts += 1

                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]

                            rh_eta = event.hbhe_pfrh_eta[i]
                            rh_phi = event.hbhe_pfrh_phi[i]

                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                if event.hbhe_pfrh_time[i] != -999 and event.hbhe_pfrh_time[i] != 100.0:
                                    hbhe_pfrh_passed_percluster_wide += 1
                                    total_hbhe_time_percluster_wide += event.hbhe_pfrh_time[i]

                            if dr <= 0.2:
                                if event.hbhe_pfrh_time[i] != -999 and event.hbhe_pfrh_time[i] != 100.0:
                                    hbhe_pfrh_passed_percluster_narrow += 1
                                    total_hbhe_time_percluster_narrow += event.hbhe_pfrh_time[i]

                            hbhe_pfrh_passed_percluster += 1
                            hbhe_pfrh_passed += 1
                            total_hbhe_pfrh_energy_percluster += event.hbhe_pfrh_energy[i]

                            if event.hbhe_pfrh_time[i] != -999 and event.hbhe_pfrh_time[i] != 100.0:
                                hbhe_pfrh_passed_percluster_time += event.hbhe_pfrh_time[i]
                                total_hbhe_time_percluster += event.hbhe_pfrh_time[i]
                                hbhe_pfrh_energy_time_percluster_pairs.append((event.hbhe_pfrh_energy[i], event.hbhe_pfrh_time[i]))

                            rh_list.append((event.hbhe_pfrh_energy[i],
                                        event.hbhe_pfrh_time[i],
                                        event.hbhe_pfrh_eta[i],
                                        event.hbhe_pfrh_phi[i],
                                        event.hbhe_pfrh_depth[i]))
                            #Fill per depth energy sums
                            if event.hbhe_pfrh_depth[i] == 1:
                                hbhe_pfrh_depth1_energy_percluster += event.hbhe_pfrh_energy[i]
                            elif event.hbhe_pfrh_depth[i] == 2:
                                hbhe_pfrh_depth2_energy_percluster += event.hbhe_pfrh_energy[i]
                            elif event.hbhe_pfrh_depth[i] == 3:
                                hbhe_pfrh_depth3_energy_percluster += event.hbhe_pfrh_energy[i]
                            elif event.hbhe_pfrh_depth[i] == 4:
                                hbhe_pfrh_depth4_energy_percluster += event.hbhe_pfrh_energy[i]
                
                #Check to see if pfcuts removed any pfrh in this cluster
                if rh_count_b4_pfcuts > 0:
                    print(f"\nCluster {cluster_idx} had a difference of {rh_count_b4_pfcuts - rh_count_after_pfcuts} pfrh before vs after cuts")
                    
                    root_delta_pfcuts.Fill(rh_count_b4_pfcuts - rh_count_after_pfcuts)

                #Store seed & pfrh info for clusters with at least 1 pfrh passing cuts
                if len(rh_list) > 0:
                    # pick the highest-energy PFRecHit in this cluster
                    pfrh_seed_idx = max(range(len(rh_list)), key=lambda j: rh_list[j][0])
                    pfrh_seed_energy, pfrh_seed_time, pfrh_seed_eta, pfrh_seed_phi, pfrh_seed_depth = rh_list[pfrh_seed_idx]
                    seedE_all.append(pfrh_seed_energy)
                    seedT_all.append(pfrh_seed_time)

                    if abs(event.hcal_eta[cluster_idx]) < 1.26:
                        seedE_HB.append(pfrh_seed_energy)
                        seedT_HB.append(pfrh_seed_time)

                    elif abs(event.hcal_eta[cluster_idx]) >= 1.26:
                        seedE_HE.append(pfrh_seed_energy)
                        seedT_HE.append(pfrh_seed_time)

                    # compare PFRecHit seed vs cluster seed
                    root_seed_pfrh_diff_depth.Fill(pfrh_seed_depth - seed_depth)
                    root_seed_pfrh_diff_eta.Fill(pfrh_seed_eta - seed_eta)
                    root_seed_pfrh_diff_phi.Fill(pfrh_seed_phi - seed_phi)

                    # fill delta time for non-seed PFRecHits in the same cluster
                    # remove all invalid time values (-999) (in this section only!) and the seed itself from the delta time histogram
                    for j, (rh_energy, rh_time, rh_eta, rh_phi, rh_depth) in enumerate(rh_list):
                        if j == pfrh_seed_idx or rh_time == -999 or rh_time == 100.0:
                            continue
                        root_hbhe_seed_pfrh_delta_time.Fill(rh_time - pfrh_seed_time)

                    if len(rh_list) < 20:
                        print(f"\nCluster {cluster_idx} has {len(rh_list)} PFrh passing cuts:")
                        print(f" Cluster Eta={event.hcal_eta[cluster_idx]:.2f}, Phi={event.hcal_phi[cluster_idx]:.2f}, Depth={event.hcal_depth[cluster_idx]}")
                        print(f" Seed: Eta={seed_eta:.2f}, Phi={seed_phi:.2f}, Depth={seed_depth}")
                        for j, (rh_energy, rh_time, rh_eta, rh_phi, rh_depth) in enumerate(rh_list):
                            if j == pfrh_seed_idx:
                                print(f"  PFrh {j} (Supposed Seed): Eta={rh_eta:.2f}, Phi={rh_phi:.2f}, Depth={rh_depth}, Energy={rh_energy:.2f}")
                            else:
                                print(f"  PFrh {j}: Eta={rh_eta:.2f}, Phi={rh_phi:.2f}, Depth={rh_depth}, Energy={rh_energy:.2f}")

                root_number_hbhe_percluster.Fill(hbhe_pfrh_passed_percluster)

                if hbhe_pfrh_passed_percluster > 0:
                    #if hbhe_pfrh_passed_percluster > 0 and total_hbhe_maxenergy_time_percluster > 0:
                    #Fill in avg number of pfrh, avg energy of pfrh, avg time per cluster                
                    root_avg_hbhe_energy_percluster.Fill(total_hbhe_pfrh_energy_percluster / hbhe_pfrh_passed_percluster)
                    
                    if hbhe_pfrh_passed_percluster_time > 0:
                        root_avg_hbhe_time_percluster.Fill(total_hbhe_time_percluster / hbhe_pfrh_passed_percluster_time) 
                    
                    if hbhe_pfrh_passed_percluster_wide > 0:
                        root_avg_hbhe_time_percluster_wide.Fill(total_hbhe_time_percluster_wide / hbhe_pfrh_passed_percluster_wide)
                    
                    if hbhe_pfrh_passed_percluster_narrow > 0:
                        root_avg_hbhe_time_percluster_narrow.Fill(total_hbhe_time_percluster_narrow / hbhe_pfrh_passed_percluster_narrow)
                    
                    #For clusters with nonzero total pfrh energy(!), fill fraction of energy per depth
                    if total_hbhe_pfrh_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_hbhe_frac_depth1_energy.Fill(hbhe_pfrh_depth1_energy_percluster / total_hbhe_pfrh_energy_percluster)
                        root_hbhe_frac_depth2_energy.Fill(hbhe_pfrh_depth2_energy_percluster / total_hbhe_pfrh_energy_percluster)
                        root_hbhe_frac_depth3_energy.Fill(hbhe_pfrh_depth3_energy_percluster / total_hbhe_pfrh_energy_percluster)
                        root_hbhe_frac_depth4_energy.Fill(hbhe_pfrh_depth4_energy_percluster / total_hbhe_pfrh_energy_percluster)

                    #For clusters with 2-3+ pfrh, find the 2-3 pfrh with max energy and fill their corresponding time
                    if len(hbhe_pfrh_energy_time_percluster_pairs) >= 3:
                        top = sorted(hbhe_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_hbhe_maxenergy_time_percluster.Fill(avg_time) 

                    elif len(hbhe_pfrh_energy_time_percluster_pairs) == 2:
                        top = sorted(hbhe_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_hbhe_maxenergy_time_percluster.Fill(avg_time)

                    
            
            # Now process HB clusters   
            if event.hcal_energy[cluster_idx] > 30 and abs(event.hcal_eta[cluster_idx]) < 1.26:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                hb_pfrh_passed_percluster = 0
                hb_pfrh_passed_percluster_time = 0

                total_hb_pfrh_energy_percluster = 0
                hb_pfrh_energy_time_percluster_pairs = []
                total_hb_time_percluster = 0

                hb_pfrh_passed_percluster_wide = 0
                total_hb_time_percluster_wide = 0
                hb_pfrh_passed_percluster_narrow = 0 
                total_hb_time_percluster_narrow = 0

                total_hb_pfrh_energy_percluster = 0
                hb_pfrh_depth1_energy_percluster = 0
                hb_pfrh_depth2_energy_percluster = 0
                hb_pfrh_depth3_energy_percluster = 0
                hb_pfrh_depth4_energy_percluster = 0

                for i in range(nHbpfrh): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.hb_pfrh_clusterIdx[i] == cluster_idx:
                        if((event.hb_pfrh_depth[i] == 1 and event.hb_pfrh_energy[i] > 0.6) or \
                            (event.hb_pfrh_depth[i] == 2 and event.hb_pfrh_energy[i] > 0.4) or \
                            (event.hb_pfrh_depth[i] == 3 and event.hb_pfrh_energy[i] > 0.4) or \
                            (event.hb_pfrh_depth[i] == 4 and event.hb_pfrh_energy[i] > 0.5)):
                            
                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]
                            rh_eta = event.hb_pfrh_eta[i]
                            rh_phi = event.hb_pfrh_phi[i]

                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                if event.hb_pfrh_time[i] != -999 and event.hb_pfrh_time[i] != 100.0:
                                    hb_pfrh_passed_percluster_wide += 1
                                    total_hb_time_percluster_wide += event.hb_pfrh_time[i]
                            if dr <= 0.2:
                                if event.hb_pfrh_time[i] != -999 and event.hb_pfrh_time[i] != 100.0:
                                    hb_pfrh_passed_percluster_narrow += 1
                                    total_hb_time_percluster_narrow += event.hb_pfrh_time[i]


                            hb_pfrh_passed_percluster += 1
                            hb_pfrh_passed += 1
                            total_hb_pfrh_energy_percluster += event.hb_pfrh_energy[i]
                            if event.hb_pfrh_time[i] != -999 and event.hb_pfrh_time[i] != 100.0:
                                hb_pfrh_passed_percluster_time += event.hb_pfrh_time[i]
                                total_hb_time_percluster += event.hb_pfrh_time[i]
                                hb_pfrh_energy_time_percluster_pairs.append((event.hb_pfrh_energy[i], event.hb_pfrh_time[i]))
                                
                            #Fill per depth energy sums
                            if event.hb_pfrh_depth[i] == 1:
                                hb_pfrh_depth1_energy_percluster += event.hb_pfrh_energy[i]
                            elif event.hb_pfrh_depth[i] == 2:
                                hb_pfrh_depth2_energy_percluster += event.hb_pfrh_energy[i]
                            elif event.hb_pfrh_depth[i] == 3:
                                hb_pfrh_depth3_energy_percluster += event.hb_pfrh_energy[i]
                            elif event.hb_pfrh_depth[i] == 4:
                                hb_pfrh_depth4_energy_percluster += event.hb_pfrh_energy[i]

                root_number_hb_percluster.Fill(hb_pfrh_passed_percluster) 
                if hb_pfrh_passed_percluster > 0:
                    #Fill in avg number of pfrh, avg energy of pfrh, avg time per cluster                
                    root_avg_hb_energy_percluster.Fill(total_hb_pfrh_energy_percluster / hb_pfrh_passed_percluster)

                    if hb_pfrh_passed_percluster_time > 0: 
                        root_avg_hb_time_percluster.Fill(total_hb_time_percluster / hb_pfrh_passed_percluster_time)

                    if hb_pfrh_passed_percluster_wide > 0:
                        root_avg_hb_time_percluster_wide.Fill(total_hb_time_percluster_wide / hb_pfrh_passed_percluster_wide)
                    if hb_pfrh_passed_percluster_narrow > 0:
                        root_avg_hb_time_percluster_narrow.Fill(total_hb_time_percluster_narrow / hb_pfrh_passed_percluster_narrow)

                    #For clusters with nonzero total pfrh energy(!), fill fraction of energy per depth
                    if total_hb_pfrh_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_hb_frac_depth1_energy.Fill(hb_pfrh_depth1_energy_percluster / total_hb_pfrh_energy_percluster)
                        root_hb_frac_depth2_energy.Fill(hb_pfrh_depth2_energy_percluster / total_hb_pfrh_energy_percluster)
                        root_hb_frac_depth3_energy.Fill(hb_pfrh_depth3_energy_percluster / total_hb_pfrh_energy_percluster)
                        root_hb_frac_depth4_energy.Fill(hb_pfrh_depth4_energy_percluster / total_hb_pfrh_energy_percluster)

                    if len(hb_pfrh_energy_time_percluster_pairs) >= 3:
                        top = sorted(hb_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_hb_maxenergy_time_percluster.Fill(avg_time) 

                    elif len(hb_pfrh_energy_time_percluster_pairs) == 2:
                        top = sorted(hb_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_hb_maxenergy_time_percluster.Fill(avg_time)
                

            
            # Now process HE clusters   
            if event.hcal_energy[cluster_idx] > 30 and abs(event.hcal_eta[cluster_idx]) >= 1.26:
                # Define variables to accumulate per-cluster values for clusters that pass cuts
                he_pfrh_passed_percluster = 0
                he_pfrh_passed_percluster_time = 0
                total_he_pfrh_energy_percluster = 0
                he_pfrh_energy_time_percluster_pairs = []
                total_he_time_percluster = 0
                total_he_maxenergy_time_percluster = 0

                he_pfrh_passed_percluster_wide = 0
                total_he_time_percluster_wide = 0
                he_pfrh_passed_percluster_narrow = 0
                total_he_time_percluster_narrow = 0

                total_he_pfrh_energy_percluster = 0
                he_pfrh_depth1_energy_percluster = 0
                he_pfrh_depth2_energy_percluster = 0
                he_pfrh_depth3_energy_percluster = 0
                he_pfrh_depth4_energy_percluster = 0

                for i in range(nHepfrh): 
                    #For all the events in this cluster, apply depth dependent energy cuts
                    #to-do: adjust energy cuts based on depth
                    if event.he_pfrh_clusterIdx[i] == cluster_idx:

                        if ((event.he_pfrh_depth[i] == 1 and event.he_pfrh_energy[i] > 0.2) or \
                            (event.he_pfrh_depth[i] == 2 and event.he_pfrh_energy[i] > 0.3) or \
                            (event.he_pfrh_depth[i] == 3 and event.he_pfrh_energy[i] > 0.3) or \
                            (event.he_pfrh_depth[i] == 4 and event.he_pfrh_energy[i] > 0.3) or \
                            (event.he_pfrh_depth[i] == 5 and event.he_pfrh_energy[i] > 0.3) or \
                            (event.he_pfrh_depth[i] == 6 and event.he_pfrh_energy[i] > 0.3) or \
                            (event.he_pfrh_depth[i] == 7 and event.he_pfrh_energy[i] > 0.3)):
                            
                            cl_eta = event.hcal_eta[cluster_idx]
                            cl_phi = event.hcal_phi[cluster_idx]
                            rh_eta = event.he_pfrh_eta[i]
                            rh_phi = event.he_pfrh_phi[i] 
                            dr = delta_r(rh_eta, rh_phi, cl_eta, cl_phi)
                            if dr <= 0.4:
                                
                                he_pfrh_passed_percluster_wide += 1
                                total_he_time_percluster_wide += event.he_pfrh_time[i]
                            if dr <= 0.2:   
                                he_pfrh_passed_percluster_narrow += 1
                                total_he_time_percluster_narrow += event.he_pfrh_time[i]

                            he_pfrh_passed_percluster += 1
                            he_pfrh_passed += 1
                            total_he_pfrh_energy_percluster += event.he_pfrh_energy[i]

                            if event.he_pfrh_time[i] != -999 and event.he_pfrh_time[i] != 100.0:
                                he_pfrh_passed_percluster_time += event.he_pfrh_time[i]
                                total_he_time_percluster += event.he_pfrh_time[i]
                                he_pfrh_energy_time_percluster_pairs.append((event.he_pfrh_energy[i], event.he_pfrh_time[i]))

                            #Fill per depth energy sums
                            if event.he_pfrh_depth[i] == 1:
                                he_pfrh_depth1_energy_percluster += event.he_pfrh_energy[i]
                            elif event.he_pfrh_depth[i] == 2:
                                he_pfrh_depth2_energy_percluster += event.he_pfrh_energy[i]
                            elif event.he_pfrh_depth[i] == 3:
                                he_pfrh_depth3_energy_percluster += event.he_pfrh_energy[i]
                            elif event.he_pfrh_depth[i] == 4:
                                he_pfrh_depth4_energy_percluster += event.he_pfrh_energy[i]
                
                root_number_he_percluster.Fill(he_pfrh_passed_percluster)

                if he_pfrh_passed_percluster > 0:
                    #Fill in avg number of pfrh, avg energy of pfrh, avg time per cluster                
                    root_avg_he_energy_percluster.Fill(total_he_pfrh_energy_percluster / he_pfrh_passed_percluster)
                    if he_pfrh_passed_percluster_time > 0:
                        root_avg_he_time_percluster.Fill(total_he_time_percluster / he_pfrh_passed_percluster_time)

                    if he_pfrh_passed_percluster_wide > 0:
                        root_avg_he_time_percluster_wide.Fill(total_he_time_percluster_wide / he_pfrh_passed_percluster_wide)
                    if he_pfrh_passed_percluster_narrow > 0:
                        root_avg_he_time_percluster_narrow.Fill(total_he_time_percluster_narrow / he_pfrh_passed_percluster_narrow)

                    #For clusters with nonzero total pfrh energy(!), fill fraction of energy per depth
                    if total_he_pfrh_energy_percluster > 0 :
                        #Fill fraction of energy fractions per depth histograms
                        root_he_frac_depth1_energy.Fill(he_pfrh_depth1_energy_percluster / total_he_pfrh_energy_percluster)
                        root_he_frac_depth2_energy.Fill(he_pfrh_depth2_energy_percluster / total_he_pfrh_energy_percluster)
                        root_he_frac_depth3_energy.Fill(he_pfrh_depth3_energy_percluster / total_he_pfrh_energy_percluster)
                        root_he_frac_depth4_energy.Fill(he_pfrh_depth4_energy_percluster / total_he_pfrh_energy_percluster)
                    
                    if len(he_pfrh_energy_time_percluster_pairs) >= 3:
                        top = sorted(he_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:3]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_he_maxenergy_time_percluster.Fill(avg_time)   
                    elif len(he_pfrh_energy_time_percluster_pairs) == 2:
                        top = sorted(he_pfrh_energy_time_percluster_pairs, key=lambda x: x[0], reverse=True)[:2]  # sort by energy
                        avg_time = sum(time for _, time in top) / len(top)
                        root_avg_he_maxenergy_time_percluster.Fill(avg_time)
    #Note:

        print(f"Number of HCAL clusters passing cuts: {hcal_cluster_passed}")
        print(f"Number of HBHE pfrh passing cuts: {hbhe_pfrh_passed}")
        print(f"Number of HB pfrh passing cuts: {hb_pfrh_passed}")
        print(f"Number of HE pfrh passing cuts: {he_pfrh_passed}")

    # Save histograms
    np.savez(
        outputfile.replace(".root", "_scatter.npz"),
        seedE_all=np.array(seedE_all, dtype=float),
        seedT_all=np.array(seedT_all, dtype=float),
        seedE_HB=np.array(seedE_HB, dtype=float),
        seedT_HB=np.array(seedT_HB, dtype=float),
        seedE_HE=np.array(seedE_HE, dtype=float),
        seedT_HE=np.array(seedT_HE, dtype=float))  

    out.Write()
    out.Close()

    print("Histograms saved to", outputfile)
