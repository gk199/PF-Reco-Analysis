#!/usr/bin/env python3
"""
Label HCAL PF clusters in the LLP ntuples with the LLP decay product they are matched to, its parent LLP, and the
expected time delay relative to a prompt particle. Run separately on each PF version (standardPF, seedTimingPF).

Writes three trees per input:
  clusterLabels: per event, all HCAL clusters with their match (clus_iLLP, clus_iGen, clus_dR, clus_expDelay)
  llpDecays:     per event, all LLP decay products with their expected delay and matched cluster summary
  config:        the matching cone (deltaR) and cluster energy cut used

Usage:
    python3 label_llp_clusters.py \
        --input  /eos/user/c/chtong/Public/Rereco/LLP/20260930T081201Z_19373/pfObjectsNtuple_standardPF.root \
        --output labels_standardPF.root
"""

import argparse
import awkward as ak
import numpy as np
import uproot

from TruthInfoHelper import TruthInfoHelper, GEN_BRANCHES, HB_CLUSTER_ETA_MAX, DeltaR

parser = argparse.ArgumentParser()
parser.add_argument("--input",   required=True, help="PFObjectsNtupler ntuple")
parser.add_argument("--output",  required=True, help="Output ROOT file with the label trees")
parser.add_argument("--deltaR",  type=float, default=0.4, help="Cluster to LLP (decay product) matching cone")
parser.add_argument("--clusterE", type=float, default=0.0, help="Cluster energy cut for LLPIsTruthMatched")
parser.add_argument("--maxEvents", type=int, default=None)
parser.add_argument("--debug", action="store_true")
args = parser.parse_args()

CLUSTER_BRANCHES = ["hcal_energy", "hcal_eta", "hcal_phi", "hcal_time", "hcal_depth", "hcal_nRecHits", "hcal_seed_depth"]

tree = uproot.open(args.input)["pfObjectsNtupler/pfTree"]
arrays = tree.arrays(["run", "lumi", "event"] + CLUSTER_BRANCHES + GEN_BRANCHES, entry_stop=args.maxEvents, library="np")
n_events = len(arrays["event"])

helper = TruthInfoHelper(debug=args.debug)
clus_out  = {k: [] for k in ["energy", "eta", "phi", "time", "depth", "nRecHits", "isHB", "iLLP", "iGen", "dR", "expDelay"]}
decay_out = {k: [] for k in ["iLLP", "iGen", "pdgId", "pt", "eta", "phi", "energy",
                             "llp_decayR", "llp_decayZ", "llp_eta", "llp_beta", "llp_flightTime", "llp_decaysInHB",
                             "expDelay", "slowness", "pathDelay", "hcalEta", "inHBAcceptance",
                             "isTruthMatched", "nMatchedClusters", "matchedEnergy", "matchedTime",
                             "minDR_shifted", "minDR_raw"]}

for i_evt in range(n_events):
    gen = {k: arrays[k][i_evt] for k in GEN_BRANCHES}
    clusters = {"eta": arrays["hcal_eta"][i_evt], "phi": arrays["hcal_phi"][i_evt], "energy": arrays["hcal_energy"][i_evt]}
    helper.LoadEvent(gen, clusters)

    # ----- Clusters -----
    isHB = np.abs(clusters["eta"]) < HB_CLUSTER_ETA_MAX
    iLLP, iGen, dR = helper.ClusterIsMatchedTo(clusters["eta"], clusters["phi"], args.deltaR)
    iLLP[~isHB] = -1
    iGen[~isHB] = -1
    dR[~isHB]   = -1

    delay_cache = {}
    for idx_llp in range(len(helper.gLLP_iGen)):
        for i_tp in helper.map_gLLP_to_gParticle_indices[idx_llp]:
            delay_cache[i_tp] = helper.ExpectedDelay(idx_llp, i_tp)
    expDelay = np.array([delay_cache[g][0] if g >= 0 else -999. for g in iGen])

    clus_out["energy"].append(clusters["energy"])
    clus_out["eta"].append(clusters["eta"])
    clus_out["phi"].append(clusters["phi"])
    clus_out["time"].append(arrays["hcal_time"][i_evt])
    clus_out["depth"].append(arrays["hcal_depth"][i_evt])
    clus_out["nRecHits"].append(arrays["hcal_nRecHits"][i_evt])
    clus_out["isHB"].append(isHB)
    clus_out["iLLP"].append(iLLP)
    clus_out["iGen"].append(iGen)
    clus_out["dR"].append(dR)
    clus_out["expDelay"].append(expDelay)

    # ----- LLP decay products -----
    rows = {k: [] for k in decay_out}
    for idx_gLLPDecay, (idx_llp, i_tp) in enumerate(zip(helper.gLLPDecay_iLLP, helper.gLLPDecay_iParticle)):
        i_llp_gen = helper.gLLP_iGen[idx_llp]
        dt, slowness, path, eta_c = delay_cache[i_tp]
        isTruthMatched, _ = helper.LLPIsTruthMatched(idx_gLLPDecay, args.clusterE, args.deltaR)
        inHB = helper.LLPDecaysInHB(idx_llp)

        mine = (iGen == i_tp)
        e_matched = clusters["energy"][mine].sum()
        t_matched = (clusters["energy"][mine] * arrays["hcal_time"][i_evt][mine]).sum() / e_matched if e_matched > 0 else -999.

        hb_eta, hb_phi = clusters["eta"][isHB], clusters["phi"][isHB]
        if helper.DecayProductsCanReachHB(idx_llp) and len(hb_eta) > 0:
            minDR_shifted = helper.DeltaR_ClusterToDecayProduct(idx_llp, i_tp, hb_eta, hb_phi).min()
        else:
            minDR_shifted = -1.
        minDR_raw = DeltaR(gen["gen_eta"][i_tp], hb_eta, gen["gen_phi"][i_tp], hb_phi).min() if len(hb_eta) > 0 else -1.

        rows["iLLP"].append(idx_llp)
        rows["iGen"].append(i_tp)
        rows["pdgId"].append(gen["gen_pdgId"][i_tp])
        rows["pt"].append(gen["gen_pt"][i_tp])
        rows["eta"].append(gen["gen_eta"][i_tp])
        rows["phi"].append(gen["gen_phi"][i_tp])
        rows["energy"].append(gen["gen_energy"][i_tp])
        rows["llp_decayR"].append(helper.gLLP_DecayVtx_R[idx_llp])
        rows["llp_decayZ"].append(helper.gLLP_DecayVtx[idx_llp][2])
        rows["llp_eta"].append(gen["gen_eta"][i_llp_gen])
        rows["llp_beta"].append(helper.gLLP_Beta[idx_llp])
        rows["llp_flightTime"].append(helper.gLLP_FlightTime[idx_llp])
        rows["llp_decaysInHB"].append(inHB)
        rows["expDelay"].append(dt)
        rows["slowness"].append(slowness)
        rows["pathDelay"].append(path)
        rows["hcalEta"].append(eta_c)
        rows["inHBAcceptance"].append((inHB or helper.DecayProductsCanReachHB(idx_llp)) and abs(eta_c) < HB_CLUSTER_ETA_MAX)
        rows["isTruthMatched"].append(isTruthMatched)
        rows["nMatchedClusters"].append(int(mine.sum()))
        rows["matchedEnergy"].append(e_matched)
        rows["matchedTime"].append(t_matched)
        rows["minDR_shifted"].append(minDR_shifted)
        rows["minDR_raw"].append(minDR_raw)

    for k in decay_out:
        decay_out[k].append(rows[k])

    if args.debug or i_evt % 50 == 0:
        print(f"Event {i_evt}/{n_events}: {len(helper.gLLP_iGen)} LLPs, {len(helper.gLLPDecay_iParticle)} decay products, "
              f"{(iGen >= 0).sum()} matched HB clusters")

event_ids = {"run": arrays["run"], "lumi": arrays["lumi"], "event": arrays["event"]}
with uproot.recreate(args.output) as fout:
    fout["clusterLabels"] = {**event_ids, "clus": ak.zip({k: ak.Array(v) for k, v in clus_out.items()})}
    fout["llpDecays"]     = {**event_ids, "decay": ak.zip({k: ak.Array(v) for k, v in decay_out.items()})}
    fout["config"]        = {"deltaR": np.array([args.deltaR]), "clusterE": np.array([args.clusterE])}
print(f"Wrote {args.output}")
