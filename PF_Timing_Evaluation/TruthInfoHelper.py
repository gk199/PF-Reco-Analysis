"""
Truth matching of HCAL PF clusters to LLPs and LLP decay products, for the PFObjectsNtupler ntuples.
Python port of the jet-matching helpers in Run3-HCAL-LLP-Analysis:
https://github.com/gk199/Run3-HCAL-LLP-Analysis/blob/main/DisplacedHcalJetAnalyzer/src/TruthInfoHelper.cxx
Function names and structure follow that file, with reco jets replaced by HCAL clusters.
"""

import numpy as np

C_CM_PER_NS = 29.9792458

LLP_PDGIDS = [9000006, 9000007, 1023, 1000023, 1000025, 6000113, 9900012, 9900014, 9900016]

HB_R_INNER = 177.   # cm, radius where cluster positions are placed (HB inner radius)
HB_R_OUTER = 295.   # cm
HB_ETA_MAX = 1.26   # LLP |eta| requirement for a decay inside HB (as in the reference)
HB_CLUSTER_ETA_MAX = 1.3  # cluster |eta| for HB clusters; PF timing cuts are HB-only

GEN_BRANCHES = [
    "gen_pt", "gen_eta", "gen_phi", "gen_energy", "gen_pdgId", "gen_status",
    "gen_decayVx", "gen_decayVy", "gen_decayVz", "gen_decayLength3D", "gen_flightTimeEstimate",
    "gen_daughterOffset", "gen_daughterIdx",
]


# ======================================================================================================================
def DeltaPhi( phi1, phi2 ):
    return (np.asarray(phi1) - np.asarray(phi2) + np.pi) % (2*np.pi) - np.pi


# ======================================================================================================================
def DeltaR( eta1, eta2, phi1, phi2 ):
    return np.hypot( np.asarray(eta1) - np.asarray(eta2), DeltaPhi(phi1, phi2) )


# ======================================================================================================================
def UnitVector( eta, phi ):
    theta = 2*np.arctan( np.exp(-eta) )
    return np.array([ np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta) ])


# ======================================================================================================================
def EtaPhi( x, y, z ):
    return np.arcsinh( z / np.hypot(x, y) ), np.arctan2( y, x )


# ======================================================================================================================
class TruthInfoHelper:

    def __init__( self, debug=False ):
        self.debug = debug

    # ==================================================================================================================
    def LoadEvent( self, gen, clusters ):
        """
        Description: Loads one event and sets the LLP variables
        Inputs: gen:         dict of per-event numpy arrays for the GEN_BRANCHES
                clusters:     dict of per-event numpy arrays with keys eta, phi, energy
        """
        self.gen = gen
        self.clusters = clusters
        self.SetLLPVariables()

    # ==================================================================================================================
    def SetLLPVariables( self ):
        """
        Description:
         -- Finds LLPs (|pdgId| in LLP_PDGIDS with a common decay vertex)
         -- Loads indices of LLP decay products (corresponding to which LLP, and which gParticle) (pt-ordered)
         -- Sets LLP decay R, magnitude, beta, and lab-frame flight time
        Decay products are the LLP daughters from the gen daughter links, rather than the production-vertex
        comparison used in the reference.
        """
        if self.debug: print("TruthInfoHelper::SetLLPVariables()")

        g = self.gen
        self.gLLP_iGen = []
        self.gLLP_DecayVtx = []
        self.gLLP_DecayVtx_R = []
        self.gLLP_DecayVtx_Mag = []
        self.gLLP_Beta = []
        self.gLLP_FlightTime = []
        self.map_gLLP_to_gParticle_indices = []

        LLPDecayProducts_temp = []

        for i_gen in np.where( np.isin( np.abs(g["gen_pdgId"]), LLP_PDGIDS ) )[0]:
            if g["gen_decayLength3D"][i_gen] < 0: continue  # no common decay vertex (e.g. intermediate copy)

            daughters = g["gen_daughterIdx"][ g["gen_daughterOffset"][i_gen] : g["gen_daughterOffset"][i_gen+1] ]
            daughters = [ int(j) for j in daughters if j >= 0 and abs(g["gen_pdgId"][j]) not in LLP_PDGIDS ]
            if len(daughters) == 0: continue

            i_llp = len(self.gLLP_iGen)
            vtx = np.array([ g["gen_decayVx"][i_gen], g["gen_decayVy"][i_gen], g["gen_decayVz"][i_gen] ])
            p = g["gen_pt"][i_gen] * np.cosh( g["gen_eta"][i_gen] )

            self.gLLP_iGen.append( int(i_gen) )
            self.gLLP_DecayVtx.append( vtx )
            self.gLLP_DecayVtx_R.append( np.hypot(vtx[0], vtx[1]) )
            self.gLLP_DecayVtx_Mag.append( np.linalg.norm(vtx) )
            self.gLLP_Beta.append( p / g["gen_energy"][i_gen] )
            self.gLLP_FlightTime.append( g["gen_flightTimeEstimate"][i_gen] )  # = L*E/(p*c) = L/(beta*c), ns
            self.map_gLLP_to_gParticle_indices.append( daughters )

            for i_tp in daughters:
                LLPDecayProducts_temp.append( ( g["gen_pt"][i_tp], i_llp, i_tp ) )

        LLPDecayProducts_temp.sort( key=lambda x: -x[0] )
        self.gLLPDecay_iLLP      = [ x[1] for x in LLPDecayProducts_temp ]
        self.gLLPDecay_iParticle = [ x[2] for x in LLPDecayProducts_temp ]

    # ==================================================================================================================
    def LLPDecaysInHB( self, idx_llp ):
        """
        Description: True if the LLP decays inside HB, in which case clusters are matched to the LLP direction
        """
        i_gen = self.gLLP_iGen[idx_llp]
        return ( self.gLLP_DecayVtx_R[idx_llp] >= HB_R_INNER and self.gLLP_DecayVtx_R[idx_llp] < HB_R_OUTER
                 and abs(self.gen["gen_eta"][i_gen]) <= HB_ETA_MAX )

    # ==================================================================================================================
    def DecayProductsCanReachHB( self, idx_llp ):
        """
        Description: True if the LLP decays before the HB inner radius, so its decay products are propagated to
        HB_R_INNER. LLPs decaying beyond HB (or at R >= HB_R_INNER outside the HB eta range) cannot be shifted into
        the LLP frame with a cluster placed at HB_R_INNER, and are not matched.
        """
        return self.gLLP_DecayVtx_R[idx_llp] < HB_R_INNER

    # ==================================================================================================================
    def ClusterPosition( self, cluster_eta, cluster_phi ):
        """
        Description: Delivers cluster coordinates, placing the cluster at the HB inner radius
        Inputs: cluster_eta, cluster_phi: scalars or arrays
        """
        cluster_eta = np.asarray(cluster_eta)
        cluster_phi = np.asarray(cluster_phi)
        return np.stack([ HB_R_INNER*np.cos(cluster_phi), HB_R_INNER*np.sin(cluster_phi), HB_R_INNER*np.sinh(cluster_eta) ], axis=-1)

    # ==================================================================================================================
    def DeltaR_ClusterToDecayProduct( self, idx_llp, idx_gParticle, cluster_eta, cluster_phi ):
        """
        Description: Delivers deltaR between clusters and an LLP decay product, after shifting the clusters into the
        LLP decay frame of reference
        Inputs: idx_llp:         LLP index (generally either 0 or 1)
                idx_gParticle:     gen index of the LLP decay product
                cluster_eta, cluster_phi: arrays
        """
        vec_clus_new = self.ClusterPosition( cluster_eta, cluster_phi ) - self.gLLP_DecayVtx[idx_llp]
        eta_new, phi_new = EtaPhi( vec_clus_new[..., 0], vec_clus_new[..., 1], vec_clus_new[..., 2] )
        return DeltaR( self.gen["gen_eta"][idx_gParticle], eta_new, self.gen["gen_phi"][idx_gParticle], phi_new )

    # ==================================================================================================================
    def DeltaR_ClusterToLLP( self, idx_llp, cluster_eta, cluster_phi ):
        """
        Description: Delivers deltaR between clusters and the LLP direction (used for LLPs decaying in HB)
        """
        i_gen = self.gLLP_iGen[idx_llp]
        return DeltaR( self.gen["gen_eta"][i_gen], cluster_eta, self.gen["gen_phi"][i_gen], cluster_phi )

    # ==================================================================================================================
    def ClusterIsMatchedTo( self, cluster_eta, cluster_phi, deltaR_cut=0.4 ):
        """
        Description: Delivers idx_llp, idx_gParticle and dR for the LLP (decay product) matched to each cluster.
        Cluster version of JetIsMatchedTo, vectorized over clusters:
         -- LLP decays in HB: cluster is matched to the LLP direction. idx_gParticle is the decay product closest
            in (unshifted) deltaR, so every matched cluster carries a decay product label
         -- LLP decays before HB: cluster is placed at HB_R_INNER, shifted into the LLP decay frame, and matched to
            the decay product direction
        Unlike the jet version, the closest candidate is taken rather than the first one passing the cut, since
        one decay product typically makes several clusters. Unmatched clusters get (-1, -1, -1).
        Inputs: cluster_eta, cluster_phi:     arrays of cluster coordinates
                deltaR_cut:                 deltaR between cluster and LLP decay prod (default: 0.4)
        """
        cluster_eta = np.asarray(cluster_eta)
        cluster_phi = np.asarray(cluster_phi)
        n_clus = len(cluster_eta)

        best_dR    = np.full( n_clus, np.inf )
        best_llp   = np.full( n_clus, -1 )
        best_gPart = np.full( n_clus, -1 )

        for idx_llp in range( len(self.gLLP_iGen) ):
            decay_products = self.map_gLLP_to_gParticle_indices[idx_llp]

            if self.LLPDecaysInHB(idx_llp):
                dR_temp = self.DeltaR_ClusterToLLP( idx_llp, cluster_eta, cluster_phi )
                dR_prod = np.stack([ DeltaR( self.gen["gen_eta"][i_tp], cluster_eta, self.gen["gen_phi"][i_tp], cluster_phi ) for i_tp in decay_products ])
                closest_prod = np.asarray(decay_products)[ np.argmin(dR_prod, axis=0) ]
                better = dR_temp < best_dR
                best_dR[better]    = dR_temp[better]
                best_llp[better]   = idx_llp
                best_gPart[better] = closest_prod[better]

            elif self.DecayProductsCanReachHB(idx_llp):
                for i_tp in decay_products:
                    dR_temp = self.DeltaR_ClusterToDecayProduct( idx_llp, i_tp, cluster_eta, cluster_phi )
                    better = dR_temp < best_dR
                    best_dR[better]    = dR_temp[better]
                    best_llp[better]   = idx_llp
                    best_gPart[better] = i_tp

        unmatched = best_dR >= deltaR_cut
        best_dR[unmatched]    = -1
        best_llp[unmatched]   = -1
        best_gPart[unmatched] = -1
        return best_llp, best_gPart, best_dR

    # ==================================================================================================================
    def LLPDecayIsTruthMatched_LLP_b( self, idx_gLLP, idx_gParticle, clusterE_cut, deltaR_cut ):
        """
        Description: Delivers true/false on if an llp decay product is matched to an HCAL cluster, and the eta of the
        first matched cluster (non-exclusive: a cluster can match several decay products)
        clusterE_cut: cut on energy of the cluster that the LLP or decay product is matched to
        deltaR_cut: deltaR between cluster and LLP decay prod (default: 0.4)
        """
        eta = self.clusters["eta"]
        phi = self.clusters["phi"]
        passE = self.clusters["energy"] >= clusterE_cut

        # Check if LLP is directly matched to a cluster
        if self.LLPDecaysInHB(idx_gLLP):
            dR_temp = self.DeltaR_ClusterToLLP( idx_gLLP, eta, phi )
        # Check if LLP decay products are matched to a cluster
        elif self.DecayProductsCanReachHB(idx_gLLP):
            dR_temp = self.DeltaR_ClusterToDecayProduct( idx_gLLP, idx_gParticle, eta, phi )
        else:
            return (False, -99.0)

        matched = np.where( passE & (dR_temp < deltaR_cut) )[0]
        if len(matched) == 0: return (False, -99.0)
        return (True, float(eta[matched[0]]))

    # ==================================================================================================================
    def LLPIsTruthMatched( self, idx_gLLPDecay, clusterE_cut, deltaR_cut ):
        """
        Description: Delivers true/false on if an llp decay product is matched to an HCAL cluster
        deltaR_cut: deltaR between cluster and LLP decay prod (default: 0.4)

        between the two similar LLPIsTruthMatched functions:
        in LLPIsTruthMatched, idx_gLLP          =  idx_llp             in LLPDecayIsTruthMatched_LLP_b
        in LLPIsTruthMatched, idx_gParticle      =  idx_gParticle     in LLPDecayIsTruthMatched_LLP_b
        """
        if idx_gLLPDecay >= len(self.gLLPDecay_iLLP): return (False, -99.0)

        idx_gParticle = self.gLLPDecay_iParticle[idx_gLLPDecay]     # index of this b-quark
        idx_gLLP      = self.gLLPDecay_iLLP[idx_gLLPDecay]         # which LLP is associated to this b-quark

        return self.LLPDecayIsTruthMatched_LLP_b( idx_gLLP, idx_gParticle, clusterE_cut, deltaR_cut )

    # ==================================================================================================================
    def DecayProductHBIntersection( self, idx_llp, idx_gParticle ):
        """
        Description: Delivers the point where the decay product reaches HCAL: the straight-line intersection of the
        decay product direction (from the LLP decay vertex) with r = HB_R_INNER, or the LLP decay vertex itself if the
        LLP decays at or beyond HB_R_INNER
        """
        vtx = self.gLLP_DecayVtx[idx_llp]
        r_vtx = self.gLLP_DecayVtx_R[idx_llp]
        if r_vtx >= HB_R_INNER: return vtx.copy()

        u = UnitVector( self.gen["gen_eta"][idx_gParticle], self.gen["gen_phi"][idx_gParticle] )
        # Solve |vtx_perp + s*u_perp| = HB_R_INNER for s > 0
        a = u[0]**2 + u[1]**2
        b = 2*( vtx[0]*u[0] + vtx[1]*u[1] )
        c = r_vtx**2 - HB_R_INNER**2
        s = ( -b + np.sqrt(b*b - 4*a*c) ) / (2*a)
        return vtx + s*u

    # ==================================================================================================================
    def ExpectedDelay( self, idx_llp, idx_gParticle ):
        """
        Description: Delivers the expected arrival time delay (ns) of an LLP decay product at HCAL, relative to a prompt
        particle from (0,0,0) travelling at c to the same point x_c:
            dt = L/(beta c) + |x_c - x_d|/c - |x_c|/c
               = L/c (1/beta - 1)                  [LLP slowness]
               + (L + |x_c - x_d| - |x_c|)/c       [path length difference, >= 0]
        The decay product is taken to travel at c in a straight line from the decay vertex x_d.
        Returns: (dt, slowness, path length difference, eta of x_c)
        """
        L   = self.gLLP_DecayVtx_Mag[idx_llp]
        x_d = self.gLLP_DecayVtx[idx_llp]
        x_c = self.DecayProductHBIntersection( idx_llp, idx_gParticle )

        t_llp    = self.gLLP_FlightTime[idx_llp]
        slowness = L / C_CM_PER_NS * ( 1./self.gLLP_Beta[idx_llp] - 1. )
        path     = ( L + np.linalg.norm(x_c - x_d) - np.linalg.norm(x_c) ) / C_CM_PER_NS
        dt       = t_llp + ( np.linalg.norm(x_c - x_d) - np.linalg.norm(x_c) ) / C_CM_PER_NS
        eta_c, _ = EtaPhi( x_c[0], x_c[1], x_c[2] )
        return dt, slowness, path, eta_c
