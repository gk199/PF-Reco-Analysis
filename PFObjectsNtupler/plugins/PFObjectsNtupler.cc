#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDigi/interface/HcalUMNioDigi.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" // ECAL rechits
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h" 
#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "FWCore/Utilities/interface/ESGetToken.h"

#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/CaloRecHitAuxSetter.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFraction.h"
#include "DataFormats/Common/interface/Ref.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/PositionVector3D.h"

#include "TTree.h"
#include "Rtypes.h"  // for UInt_t, ULong64_t
#include <vector>
#include <cmath>
#include <utility>


class PFObjectsNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PFObjectsNtupler(const edm::ParameterSet&);
  ~PFObjectsNtupler() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidatesToken_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> ecalClustersToken_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> hcalClustersToken_;     // particleFlowClusterHCAL  (post-depth-stacking)

  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit>> hbheRechitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> ebRechitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRechitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> esRechitsToken_;
  edm::EDGetTokenT<std::vector<reco::PFBlock>> pfBlocksToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;
  edm::EDGetTokenT<HcalUMNioDigi> uMNioToken_;


  void beginRun(const edm::Run&, const edm::EventSetup&);

  // Output tree and variables
  TTree* tree_;
  UInt_t    run_;
  UInt_t    lumi_;
  ULong64_t event_;

  // PF Candidates
  std::vector<float> pf_pt_, pf_eta_, pf_phi_, pf_energy_;
  std::vector<int> pf_charge_, pf_pdgId_;

  // ECAL clusters
  std::vector<float> ecal_energy_, ecal_eta_, ecal_phi_, ecal_time_;
  std::vector<int> ecal_clusterIdx_;
  // ECAL rechits:
  // EB
  std::vector<float> eb_rechit_energy_; std::vector<float> eb_rechit_eta_; std::vector<float> eb_rechit_phi_; 
  std::vector<float> eb_rechit_time_; std::vector<int> eb_rechit_clusterIdx_; 
  std::vector<int> eb_rechit_counts_;
  // EE rechits associated to clusters
  std::vector<float> ee_rechit_energy_; std::vector<float> ee_rechit_eta_; std::vector<float> ee_rechit_phi_; 
  std::vector<float> ee_rechit_time_; std::vector<int> ee_rechit_clusterIdx_; 
  std::vector<int> ee_rechit_counts_;

  // HCAL clusters
  std::vector<float> hcal_energy_, hcal_eta_, hcal_phi_, hcal_time_, hcal_depth_;
  std::vector<float> hcal_seed_eta_;
  std::vector<float> hcal_seed_phi_;
  std::vector<int>   hcal_seed_depth_;
  // HBHE rechits associated to clusters
  std::vector<int> hbhe_rechit_counts_;
  std::vector<float> hbhe_rechit_energy_; std::vector<float> hbhe_rechit_eta_; std::vector<float> hbhe_rechit_phi_; 
  std::vector<float> hbhe_rechit_depth_; std::vector<float> hbhe_rechit_time_; std::vector<int> hbhe_rechit_tdc_; //std::vector<int> hbhe_rechit_clusterIndex_; 
  std::vector<int> hbheRechit_clusterIdx_; std::vector<int> clusterIdx_;
  std::vector<int> hbhe_ietaAbs_;
  // HB and HE rechits separately
  std::vector<int> hb_rechit_counts_;
  std::vector<int> he_rechit_counts_;
  std::vector<float> hb_rechit_eta_; std::vector<float> hb_rechit_phi_;
  std::vector<float> hb_rechit_ieta_; std::vector<float> hb_rechit_iphi_;
  std::vector<float> he_rechit_eta_; std::vector<float> he_rechit_phi_;
  std::vector<float> he_rechit_ieta_; std::vector<float> he_rechit_iphi_;
  std::vector<int> hb_rechit_tdc_; std::vector<int> he_rechit_tdc_;
  std::vector<float> hb_rechit_depth_; std::vector<float> he_rechit_depth_;
  std::vector<float> hb_rechit_energy_;
  std::vector<float> he_rechit_energy_;
  std::vector<int> hb_rechit_clusterIdx_;
  std::vector<int> he_rechit_clusterIdx_;

  // PF HBHE RecHits associated for HCAL
  std::vector<float> hbhe_pfrh_energy_, hbhe_pfrh_energyFracInCluster_;
  std::vector<double> hbhe_pfrh_eta_, hbhe_pfrh_phi_;
  std::vector<float> hbhe_pfrh_depth_;
  std::vector<float> hbhe_pfrh_time_; // if available
  std::vector<int> hbhe_pfrh_clusterIdx_;
  std::vector<int> hbhe_pfrh_counts_;
  std::vector<float> hbhe_pfrh_fracInCluster_;
  std::vector<int> hbhe_pfrh_ieta_;
  std::vector<int>   hbhe_pfrh_iphi_;
  

  // PF HB and HE RecHits for HCAL
  std::vector<float> hb_pfrh_energy_, he_pfrh_energy_;
  std::vector<float> hb_pfrh_energyFracInCluster_, he_pfrh_energyFracInCluster_;
  std::vector<double> hb_pfrh_eta_, he_pfrh_eta_;
  std::vector<double> hb_pfrh_phi_, he_pfrh_phi_;
  std::vector<float> hb_pfrh_depth_, he_pfrh_depth_;
  std::vector<float> hb_pfrh_time_, he_pfrh_time_; // if available
  std::vector<int> hb_pfrh_clusterIdx_, he_pfrh_clusterIdx_;
  std::vector<int> hb_pfrh_counts_, he_pfrh_counts_;
  std::vector<float> hb_pfrh_fracInCluster_, he_pfrh_fracInCluster_; 
  std::vector<int>   hb_pfrh_iphi_, he_pfrh_iphi_;
  std::vector<int> hb_pfrh_ieta_, he_pfrh_ieta_;


  // PF Blocks (just store number of elements for now)
  int num_pfBlocks_;

  // uMNio
  int laserType_;
};

PFObjectsNtupler::PFObjectsNtupler(const edm::ParameterSet& iConfig)
  : caloGeometryToken_{esConsumes<CaloGeometry, CaloGeometryRecord>()}
{
  usesResource("TFileService");

  pfCandidatesToken_ = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"));
  ecalClustersToken_ = consumes<std::vector<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("ecalClusters"));
  hcalClustersToken_ = consumes<std::vector<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("hcalClusters"));
  hbheRechitsToken_ = consumes<edm::SortedCollection<HBHERecHit>>(edm::InputTag("hbhereco", "", "ReRECO")); // make sure this matches the input file! 
  ebRechitsToken_ = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit", "EcalRecHitsEB", "ReRECO")); 
  eeRechitsToken_ = consumes<EcalRecHitCollection>(edm::InputTag("ecalRecHit", "EcalRecHitsEE", "ReRECO")); 
  esRechitsToken_ = consumes<EcalRecHitCollection>(edm::InputTag("ecalPreshowerRecHit", "EcalRecHitsES", "ReRECO")); 
  // hbheRechitsToken_ = consumes<std::vector<reco::PFRecHit>>(edm::InputTag("particleFlowRecHitHBHE", "Cleaned", "ReRECOtoAOD"));
  // hbheRechitsToken_ = consumes<std::vector<reco::PFRecHit>>(edm::InputTag("particleFlowRecHitHBHE", "", "ReRECO"));
  pfBlocksToken_ = consumes<std::vector<reco::PFBlock>>(iConfig.getParameter<edm::InputTag>("pfBlocks"));
  uMNioToken_ = mayConsume<HcalUMNioDigi>(iConfig.getUntrackedParameter<edm::InputTag>("taguMNio", edm::InputTag("hcalDigis")));
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("pfTree", "PF objects");
  
  // Event info branches
  tree_->Branch("run",   &run_,   "run/i");
  tree_->Branch("lumi",  &lumi_,  "lumi/i");
  tree_->Branch("event", &event_, "event/l");

  // PF candidate branches
  tree_->Branch("pf_pt", &pf_pt_);
  tree_->Branch("pf_eta", &pf_eta_);
  tree_->Branch("pf_phi", &pf_phi_);
  tree_->Branch("pf_energy", &pf_energy_);
  tree_->Branch("pf_charge", &pf_charge_);
  tree_->Branch("pf_pdgId", &pf_pdgId_);

  // ECAL cluster branches
  tree_->Branch("ecal_energy", &ecal_energy_);
  tree_->Branch("ecal_eta", &ecal_eta_);
  tree_->Branch("ecal_phi", &ecal_phi_);
  tree_->Branch("ecal_time", &ecal_time_);
  tree_->Branch("ecal_clusterIdx", &ecal_clusterIdx_);

  // ECAL rechit branches
  // EB rechits
  tree_->Branch("eb_rechit_energy", &eb_rechit_energy_);
  tree_->Branch("eb_rechit_eta", &eb_rechit_eta_);
  tree_->Branch("eb_rechit_phi", &eb_rechit_phi_);
  tree_->Branch("eb_rechit_time", &eb_rechit_time_);
  tree_->Branch("eb_rechit_clusterIdx", &eb_rechit_clusterIdx_);
  tree_->Branch("eb_rechit_counts", &eb_rechit_counts_);
  // EE rechits
  tree_->Branch("ee_rechit_energy", &ee_rechit_energy_);
  tree_->Branch("ee_rechit_eta", &ee_rechit_eta_);
  tree_->Branch("ee_rechit_phi", &ee_rechit_phi_);
  tree_->Branch("ee_rechit_time", &ee_rechit_time_);
  tree_->Branch("ee_rechit_clusterIdx", &ee_rechit_clusterIdx_); 
  tree_->Branch("ee_rechit_counts", &ee_rechit_counts_);

  // HCAL cluster branches (post-depth-stacking: particleFlowClusterHCAL)
  tree_->Branch("hcal_energy", &hcal_energy_);
  tree_->Branch("hcal_eta", &hcal_eta_);
  tree_->Branch("hcal_phi", &hcal_phi_);
  tree_->Branch("hcal_time", &hcal_time_);
  tree_->Branch("hcal_depth", &hcal_depth_);
  tree_->Branch("clusterIdx", &clusterIdx_);

  tree_->Branch("hcal_seed_eta",   &hcal_seed_eta_);
  tree_->Branch("hcal_seed_phi",   &hcal_seed_phi_);
  tree_->Branch("hcal_seed_depth", &hcal_seed_depth_);

  // HBHE rechit branches
  tree_->Branch("hbhe_rechit_energy", &hbhe_rechit_energy_);
  tree_->Branch("hbhe_rechit_eta", &hbhe_rechit_eta_);
  tree_->Branch("hbhe_rechit_phi", &hbhe_rechit_phi_);
  tree_->Branch("hbhe_rechit_depth", &hbhe_rechit_depth_);
  tree_->Branch("hbhe_rechit_time", &hbhe_rechit_time_);
  tree_->Branch("hbhe_rechit_tdc", &hbhe_rechit_tdc_);
  tree_->Branch("hbheRechit_clusterIdx", &hbheRechit_clusterIdx_);
  tree_->Branch("hbhe_rechit_ietaAbs", &hbhe_ietaAbs_);
  tree_->Branch("hbhe_rechit_counts", &hbhe_rechit_counts_);

  tree_->Branch("hb_rechit_tdc", &hb_rechit_tdc_);
  tree_->Branch("he_rechit_tdc", &he_rechit_tdc_);
  tree_->Branch("hb_rechit_depth", &hb_rechit_depth_);
  tree_->Branch("he_rechit_depth", &he_rechit_depth_);
  tree_->Branch("hb_rechit_counts", &hb_rechit_counts_);
  tree_->Branch("he_rechit_counts", &he_rechit_counts_);
  tree_->Branch("hb_rechit_energy", &hb_rechit_energy_);  
  tree_->Branch("he_rechit_energy", &he_rechit_energy_);
  tree_->Branch("hb_rechit_clusterIdx", &hb_rechit_clusterIdx_);
  tree_->Branch("he_rechit_clusterIdx", &he_rechit_clusterIdx_);

  tree_->Branch("hb_rechit_eta", &hb_rechit_eta_);
  tree_->Branch("hb_rechit_phi", &hb_rechit_phi_);
  tree_->Branch("hb_rechit_ieta", &hb_rechit_ieta_);
  tree_->Branch("hb_rechit_iphi", &hb_rechit_iphi_);
  tree_->Branch("he_rechit_eta", &he_rechit_eta_);
  tree_->Branch("he_rechit_phi", &he_rechit_phi_);
  tree_->Branch("he_rechit_ieta", &he_rechit_ieta_);
  tree_->Branch("he_rechit_iphi", &he_rechit_iphi_);

  
  // PF HBHE RecHits associated to HCAL clusters
  tree_->Branch("hbhe_pfrh_energy", &hbhe_pfrh_energy_);
  tree_->Branch("hbhe_pfrh_energyFracInCluster", &hbhe_pfrh_energyFracInCluster_);
  tree_->Branch("hbhe_pfrh_eta", &hbhe_pfrh_eta_);
  tree_->Branch("hbhe_pfrh_phi", &hbhe_pfrh_phi_);
  tree_->Branch("hbhe_pfrh_depth", &hbhe_pfrh_depth_);
  tree_->Branch("hbhe_pfrh_time", &hbhe_pfrh_time_);
  tree_->Branch("hbhe_pfrh_clusterIdx", &hbhe_pfrh_clusterIdx_);
  tree_->Branch("hbhe_pfrh_counts", &hbhe_pfrh_counts_);
  tree_->Branch("hbhe_pfrh_ieta", &hbhe_pfrh_ieta_);
  tree_->Branch("hbhe_pfrh_iphi", &hbhe_pfrh_iphi_);
  tree_->Branch("hbhe_pfrh_fracInCluster", &hbhe_pfrh_fracInCluster_);


  // PF HB and HE RecHits for HCAL
  tree_->Branch("hb_pfrh_energy", &hb_pfrh_energy_);
  tree_->Branch("he_pfrh_energy", &he_pfrh_energy_);
  tree_->Branch("hb_pfrh_energyFracInCluster", &hb_pfrh_energyFracInCluster_);
  tree_->Branch("he_pfrh_energyFracInCluster", &he_pfrh_energyFracInCluster_);
  tree_->Branch("hb_pfrh_fracInCluster", &hb_pfrh_fracInCluster_);
  tree_->Branch("he_pfrh_fracInCluster", &he_pfrh_fracInCluster_);
  tree_->Branch("hb_pfrh_eta", &hb_pfrh_eta_);
  tree_->Branch("he_pfrh_eta", &he_pfrh_eta_);
  tree_->Branch("hb_pfrh_phi", &hb_pfrh_phi_);
  tree_->Branch("he_pfrh_phi", &he_pfrh_phi_);
  tree_->Branch("hb_pfrh_depth", &hb_pfrh_depth_);
  tree_->Branch("he_pfrh_depth", &he_pfrh_depth_);
  tree_->Branch("hb_pfrh_time", &hb_pfrh_time_);
  tree_->Branch("he_pfrh_time", &he_pfrh_time_);
  tree_->Branch("hb_pfrh_clusterIdx", &hb_pfrh_clusterIdx_);
  tree_->Branch("he_pfrh_clusterIdx", &he_pfrh_clusterIdx_);
  tree_->Branch("hb_pfrh_counts", &hb_pfrh_counts_);
  tree_->Branch("he_pfrh_counts", &he_pfrh_counts_);
  tree_->Branch("hb_pfrh_ieta", &hb_pfrh_ieta_);
  tree_->Branch("he_pfrh_ieta", &he_pfrh_ieta_);
  tree_->Branch("hb_pfrh_iphi", &hb_pfrh_iphi_);
  tree_->Branch("he_pfrh_iphi", &he_pfrh_iphi_);

  // PF block info
  tree_->Branch("num_pfBlocks", &num_pfBlocks_);

  // uMNio
  tree_->Branch("laserType", &laserType_);
}

// Convert ieta to eta using HCAL mapping

static double hcalEtaFromIeta(int ieta) {
    // HB: |ieta| <= 16, HE: 17 <= |ieta| <= 28
    int sign = (ieta >= 0 ? 1 : -1);
    int absi = std::abs(ieta);

    double eta = 0.0;
    if (absi >= 24){
      eta = 0.1695 * ieta - sign*1.9931; 
    } else {
      eta = 0.0875*ieta - sign*0.0489;
    }
    // if (absi <= 16) { // HB
    //     eta = 0.087 * (absi - 0.5);
    // } else if (absi <= 28) { // HE
    //     eta = 0.087 * (16 + (absi - 16) * 0.9); // approximate
    // } else { // HF etc; nothing should be outside really
    //     eta = 40.0; //placeholder
    // }
    return eta;
}

// Convert iphi to phi (HB/HE have 72 phi bins)
static double hcalPhiFromIphi(int iphi) {
    // HCAL iphi runs from 1..72
    double phi = (iphi - 0.5) * (M_PI / 36.0); // 2pi/72
    // Put phi into -pi, pi
    if (phi > M_PI) phi -= 2.0 * M_PI;

    return phi;
}

// Master helper: return (eta, phi) for HBHE DetId
inline std::pair<double,double> hcalEtaPhiFromDetId(const HcalDetId& detid) {
    int ieta = detid.ieta();  // signed eta index
    int iphi = detid.iphi();  // 1..72

    double eta = hcalEtaFromIeta(ieta);
    double phi = hcalPhiFromIphi(iphi);

    return {eta, phi};
}


void PFObjectsNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Clear all vectors
  pf_pt_.clear(); pf_eta_.clear(); pf_phi_.clear(); pf_energy_.clear(); pf_charge_.clear(); pf_pdgId_.clear();
  
  ecal_energy_.clear(); ecal_eta_.clear(); ecal_phi_.clear(); ecal_time_.clear();
  ecal_clusterIdx_.clear();

  eb_rechit_energy_.clear(); eb_rechit_eta_.clear(); eb_rechit_phi_.clear(); eb_rechit_time_.clear(); eb_rechit_clusterIdx_.clear();
  ee_rechit_energy_.clear(); ee_rechit_eta_.clear(); ee_rechit_phi_.clear(); 
  ee_rechit_time_.clear(); ee_rechit_clusterIdx_.clear();
  ee_rechit_counts_.clear();
  eb_rechit_counts_.clear(); 

  hcal_energy_.clear(); hcal_eta_.clear(); hcal_phi_.clear(); hcal_time_.clear(); hcal_depth_.clear();
  hcal_seed_eta_.clear(); hcal_seed_phi_.clear(); hcal_seed_depth_.clear();

  hbhe_rechit_energy_.clear(); hbhe_rechit_eta_.clear(); 
  hbhe_rechit_phi_.clear(); hbhe_rechit_depth_.clear(); hbhe_rechit_time_.clear(); //hbhe_rechit_clusterIndex_.clear();
  hbhe_rechit_tdc_.clear();
  hbheRechit_clusterIdx_.clear(); clusterIdx_.clear();
  hbhe_ietaAbs_.clear();
  hbhe_rechit_counts_.clear();

  hb_rechit_tdc_.clear(); he_rechit_tdc_.clear();
  hb_rechit_depth_.clear();
  he_rechit_depth_.clear();
  hb_rechit_counts_.clear();  
  he_rechit_counts_.clear(); 
  hb_rechit_energy_.clear();
  he_rechit_energy_.clear();
  hb_rechit_clusterIdx_.clear();
  he_rechit_clusterIdx_.clear();

  hb_rechit_eta_.clear(); hb_rechit_phi_.clear(); hb_rechit_ieta_.clear(); hb_rechit_iphi_.clear();
  he_rechit_eta_.clear(); he_rechit_phi_.clear(); he_rechit_ieta_.clear(); he_rechit_iphi_.clear();

  hbhe_pfrh_energy_.clear(); hbhe_pfrh_energyFracInCluster_.clear();
  hbhe_pfrh_eta_.clear(); hbhe_pfrh_phi_.clear();
  hbhe_pfrh_depth_.clear();
  hbhe_pfrh_time_.clear();
  hbhe_pfrh_clusterIdx_.clear();
  hbhe_pfrh_counts_.clear();
  hbhe_pfrh_ieta_.clear();
  hbhe_pfrh_iphi_.clear();
  hbhe_pfrh_fracInCluster_.clear();

  hb_pfrh_energy_.clear(); he_pfrh_energy_.clear();
  hb_pfrh_energyFracInCluster_.clear(); he_pfrh_energyFracInCluster_.clear();
  hb_pfrh_eta_.clear(); he_pfrh_eta_.clear();
  hb_pfrh_phi_.clear(); he_pfrh_phi_.clear();
  hb_pfrh_depth_.clear(); he_pfrh_depth_.clear();
  hb_pfrh_time_.clear(); he_pfrh_time_.clear();
  hb_pfrh_clusterIdx_.clear(); he_pfrh_clusterIdx_.clear();
  hb_pfrh_counts_.clear(); he_pfrh_counts_.clear();
  hb_pfrh_ieta_.clear(); he_pfrh_ieta_.clear();
  hb_pfrh_iphi_.clear(); he_pfrh_iphi_.clear();
  hb_pfrh_fracInCluster_.clear(); he_pfrh_fracInCluster_.clear();

  auto fill_pfrh = [&](auto& pfrh_energy,
                      auto& pfrh_energyFracInCluster,
                      auto& pfrh_fracInCluster,
                      auto& pfrh_eta,
                      auto& pfrh_phi,
                      auto& pfrh_depth,
                      auto& pfrh_clusterIdx,
                      auto& pfrh_ieta,
                      auto& pfrh_iphi,
                      auto& pfrh_time,
                      float  energy,
                      float  energyFracInCluster,
                      float  fracInCluster,
                      double eta,          // <-- was float
                      double phi,          // <-- was float
                      int    depth,
                      int    clusterIdx,
                      int    ieta,
                      int    iphi,
                      float  time)
  {
    pfrh_energy.push_back(energy);
    pfrh_energyFracInCluster.push_back(energyFracInCluster);
    pfrh_fracInCluster.push_back(fracInCluster);
    pfrh_eta.push_back(eta);
    pfrh_phi.push_back(phi);
    pfrh_depth.push_back(depth);
    pfrh_clusterIdx.push_back(clusterIdx);
    pfrh_ieta.push_back(ieta);
    pfrh_iphi.push_back(iphi);
    pfrh_time.push_back(time);
  };

  num_pfBlocks_ = 0;
  laserType_ = -1000;

  // PF Candidates
  edm::Handle<std::vector<reco::PFCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidatesToken_, pfCandidates);

  run_   = iEvent.id().run();
  lumi_  = iEvent.id().luminosityBlock();
  event_ = iEvent.id().event();

  if (pfCandidates.isValid()) {
    for (const auto& cand : *pfCandidates) {
      pf_pt_.push_back(cand.pt());
      pf_eta_.push_back(cand.eta());
      pf_phi_.push_back(cand.phi());
      pf_energy_.push_back(cand.energy());
      pf_charge_.push_back(cand.charge());
      pf_pdgId_.push_back(cand.pdgId());
    }
  }


  // Get calorimeter geometry

  const CaloGeometry& caloGeom = iSetup.getData(caloGeometryToken_);

  // Get EB subdetector geometry once per event
  const CaloSubdetectorGeometry* ebGeometry =
    caloGeom.getSubdetectorGeometry(DetId::Ecal, EcalBarrel);

  // Get EE subdetector geometry 
  const CaloSubdetectorGeometry* eeGeometry =
    caloGeom.getSubdetectorGeometry(DetId::Ecal, EcalEndcap);

  // ECAL Clusters
  edm::Handle<std::vector<reco::PFCluster>> ecalClusters;
  iEvent.getByToken(ecalClustersToken_, ecalClusters);

  // ECAL RecHits 
  edm::Handle<EcalRecHitCollection> ebRecHits;
  edm::Handle<EcalRecHitCollection> eeRecHits;
  edm::Handle<EcalRecHitCollection> esRecHits;
  iEvent.getByToken(ebRechitsToken_, ebRecHits);
  iEvent.getByToken(eeRechitsToken_, eeRecHits);
  iEvent.getByToken(esRechitsToken_, esRecHits);
  
  int ecal_clusterIndex = 0;

  // Loop over ecal clusters, and for each cluster loop over eb and ee rechits 
  if (ecalClusters.isValid()) {
    for (const auto& cl : *ecalClusters) {
      
      // Save the ecal cluster info
      ecal_energy_.push_back(cl.energy());
      ecal_eta_.push_back(cl.eta());
      ecal_phi_.push_back(cl.phi());
      ecal_time_.push_back(cl.time());
      ecal_clusterIdx_.push_back(ecal_clusterIndex);

      // Loop over EB rechits   
      if (ebRecHits.isValid()) {
        int eb_count = 0;
        for (const auto& rh : *ebRecHits) {

          eb_count++;

          // Raw detid for this rechit
          DetId detid = rh.id();
          // Check to make sure that we are in Ecal Barrel; if not, skip.
          if (detid.subdetId() != EcalBarrel) continue;

          // Convert the generic DetId into an ECAL-barrel-specific ID
          EBDetId ebid(detid);

          // Get position info of the cell from EB rechit geometry
          const CaloCellGeometry* cell = ebGeometry->getGeometry(ebid);
          if (!cell) continue;

          // Ask the geometry object for the center position of this crystal
          GlobalPoint pos = cell->getPosition();
          double rh_eta = pos.eta();
          double rh_phi = pos.phi();

          //Select hits that are close to the cluster
          if (std::abs(rh_eta - cl.eta()) > 0.4) continue;
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.4) continue; //Note: originally 0.2

          // Save the rechit info
          eb_rechit_energy_.push_back(rh.energy());
          eb_rechit_eta_.push_back(rh_eta);
          eb_rechit_phi_.push_back(rh_phi);
          eb_rechit_time_.push_back(rh.time());
          eb_rechit_clusterIdx_.push_back(ecal_clusterIndex); // save cluster index association
        }
        eb_rechit_counts_.push_back(eb_count);
      }
      
      // Loop over EE rechits   
      if (eeRecHits.isValid()) {
        int ee_count = 0;

        for (const auto& rh : *eeRecHits) {
          ee_count++;
         
          // Raw detid for this rechit
          DetId detid = rh.id();
          if (detid.subdetId() != EcalEndcap) continue;  

          //Get Ecal Endcap-specific DetId
          EEDetId eeid(detid);

          //Get geometry for this EE cell
          const CaloCellGeometry* cell = eeGeometry->getGeometry(eeid);
          if (!cell) continue;

          //Get position in global coordinates
          GlobalPoint pos = cell->getPosition();
          double rh_eta = pos.eta();
          double rh_phi = pos.phi();

          if (std::abs(rh_eta - cl.eta()) > 0.4) continue;
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.4) continue;

          // Save the rechit info
          ee_rechit_energy_.push_back(rh.energy());
          ee_rechit_eta_.push_back(rh_eta);
          ee_rechit_phi_.push_back(rh_phi); 
          ee_rechit_time_.push_back(rh.time());
          ee_rechit_clusterIdx_.push_back(ecal_clusterIndex); // save cluster index association
        }
        ee_rechit_counts_.push_back(ee_count);
      }
      ecal_clusterIndex++;
    }

  }

  // HCAL Clusters
  edm::Handle<std::vector<reco::PFCluster>> hcalClusters;
  iEvent.getByToken(hcalClustersToken_, hcalClusters);

  // PF HBHE RecHits (raw rechits)
  edm::Handle<edm::SortedCollection<HBHERecHit>> hbheRechits;
  iEvent.getByToken(hbheRechitsToken_, hbheRechits);
  // PF Rechits cuts for 2025: https://indico.cern.ch/event/1604363/#3-hcal-hcal-scale-summary-of-i
  static const std::vector<float> HB_depth_Emin = {0.6, 0.4, 0.4, 0.5};
  static const std::vector<float> HE_depth_Emin = {0.2, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3};

  int clusterIndex = 0;
  if (hcalClusters.isValid()) {
    for (const auto& cl : *hcalClusters) {

      DetId seedDetId = cl.seed();
      // Defaults (in case seed/cell is invalid)
      
      float seed_eta = -999.f;
      float seed_phi = -999.f;
      int seed_depth = -999;

      // Seed cell indices (HCAL only)
      if (seedDetId.rawId() != 0 && seedDetId.det() == DetId::Hcal) {
          HcalDetId hcalSeed(seedDetId);
          seed_depth = hcalSeed.depth();
      }

      // Seed geometric eta/phi from CaloGeometry
      if (seedDetId.rawId() != 0) {
          const CaloCellGeometry* cell = caloGeom.getGeometry(seedDetId);
          if (cell) {
              const GlobalPoint& p = cell->getPosition();
              seed_eta = p.eta();
              seed_phi = p.phi();
          }
      }

      clusterIdx_.push_back(clusterIndex);
      hcal_energy_.push_back(cl.energy());
      hcal_eta_.push_back(cl.eta());
      hcal_phi_.push_back(cl.phi());
      hcal_time_.push_back(cl.time());
      hcal_depth_.push_back(cl.depth());

      // New seed branches
      hcal_seed_eta_.push_back(seed_eta);
      hcal_seed_phi_.push_back(seed_phi);
      hcal_seed_depth_.push_back(seed_depth);

      // std::cout << "Cluster eta: " << cl.eta() << " phi: " << cl.phi() << std::endl;

      // Loop over HBHE rechits (linked geometrically) associated to this HCAL cluster. Stored as hbhereco, these are the raw ones instead of PF (since that was a transitory collection)
      if (hbheRechits.isValid()) {

        int hbhe_count = 0;
        int hb_count = 0;
        int he_count = 0;

        for (const auto& rh : *hbheRechits) {

          hbhe_count++;
          hbhe_rechit_counts_.push_back(hbhe_count);

          // Get position info from HBHE rechit geometry
          // Takes the HCAL rechit’s detector ID and turns it into an HcalDetId object
          //rh.id() = raw DetId of the cell that the rechit corresponds to
          HcalDetId detid = rh.id();
          // std::cout << "Hit energy: " << rh.energy()
          //     << " detId: " << detid.rawId()
          //     << " depth: " << detid.depth() << std::endl;

          auto [rh_eta, rh_phi] = hcalEtaPhiFromDetId(detid);

          if (std::abs(rh_eta - cl.eta()) > 0.4) continue;
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.4) continue;

          // ---------------- depth-dependent cut (RAW HBHE rechits) ----------------
          float Emin = -1.0f;
          const int depth   = detid.depth();
          const int ietaAbs = detid.ietaAbs();

          // --- current boundary-based HB/HE choice ---
          const bool isHE_boundary = ( (depth == 4 && ietaAbs == 16) || (ietaAbs > 16) );
          const bool isHB_boundary = !isHE_boundary;

          if (isHB_boundary) {
            if (depth >= 1 && depth <= (int)HB_depth_Emin.size())
              Emin = HB_depth_Emin[depth - 1];
          } else { // HE by boundary
            if (depth >= 1 && depth <= (int)HE_depth_Emin.size())
              Emin = HE_depth_Emin[depth - 1];
          }

          /*
          // --- Preferred / cleaner choice: use the official subdetector ---
          // 
          // const int subdet = detid.subdetId();  // HcalBarrel / HcalEndcap
          // if (subdet == HcalBarrel) {
          //   if (depth >= 1 && depth <= (int)HB_depth_Emin.size())
          //     Emin = HB_depth_Emin[depth - 1];
          // } else if (subdet == HcalEndcap) {
          //   if (depth >= 1 && depth <= (int)HE_depth_Emin.size())
          //     Emin = HE_depth_Emin[depth - 1];
          // } else {
          //   continue; // should not happen for HBHERecHit
          // }
          */

          // If depth is out of range, drop it
          if (Emin < 0.0f) continue;

          // Apply the cut on rechit energy
          if (rh.energy() < Emin) continue;
// ------------------------------------------------------------------------

          // Save the rechit info
          hbhe_rechit_energy_.push_back(rh.energy());
          hbhe_rechit_eta_.push_back(rh_eta);
          hbhe_rechit_phi_.push_back(rh_phi);
          hbhe_rechit_depth_.push_back(detid.depth());
          hbhe_ietaAbs_.push_back(ietaAbs);

          // Getting TDC value from auxTDC field
          int six_bits_mask = 0x3f;  // 6-bit mask
          int ts = 3;                // TS3 is SOI
          int SOI_TDC = CaloRecHitAuxSetter::getField(rh.auxTDC(), six_bits_mask, ts * 6);
          hbhe_rechit_tdc_.push_back(SOI_TDC);  //  pushback TDC value
          hbhe_rechit_time_.push_back(rh.time()); // MAHI reconstructed time

          // Saving HB and HE info separately
          if ( (detid.depth() == 4 && detid.ietaAbs() == 16) || (detid.ietaAbs() > 16) ) {
            // HE rechit
            he_rechit_tdc_.push_back(SOI_TDC);
            he_rechit_depth_.push_back(detid.depth());
            he_count++;
            he_rechit_counts_.push_back(he_count);
            he_rechit_energy_.push_back(rh.energy());
            he_rechit_clusterIdx_.push_back(clusterIndex); 
            he_rechit_eta_.push_back(rh_eta);
            he_rechit_phi_.push_back(rh_phi);
            he_rechit_ieta_.push_back(detid.ieta());
            he_rechit_iphi_.push_back(detid.iphi());
          } else {
            // HB rechit
            hb_rechit_tdc_.push_back(SOI_TDC);
            hb_rechit_depth_.push_back(detid.depth());
            hb_count++;
            hb_rechit_counts_.push_back(hb_count);
            hb_rechit_energy_.push_back(rh.energy());
            hb_rechit_clusterIdx_.push_back(clusterIndex);
            hb_rechit_eta_.push_back(rh_eta);
            hb_rechit_phi_.push_back(rh_phi);
            hb_rechit_ieta_.push_back(detid.ieta());
            hb_rechit_iphi_.push_back(detid.iphi());
          }
          
          hbheRechit_clusterIdx_.push_back(clusterIndex); // save cluster index association so it is possible to map backwards to the cluster this rechit was near);
          //hbhe_rechit_clusterIndex_.push_back(hcal_energy_.size() - 1); // save cluster index association so it is possible to map backwards to the cluster this rechit was near
        }
      }

      // Loop over PF RecHits associated to this HCAL cluster
      // loop the PFRecHits that make up this PFCluster

      int hebh_pfrh_count = 0;
      int hb_pfrh_count = 0;
      int he_pfrh_count = 0;

      for (const auto& hitRefAndFrac : cl.recHitFractions()) 
      {
        const auto& pfrh_ref = hitRefAndFrac.recHitRef();
        if (pfrh_ref.isNull()) {
          std::cout << "[PFObjectsNtupler] PFRecHit ref is NULL\n";
          continue;
        }
        if (!pfrh_ref.isAvailable()) {
          std::cout << "[PFObjectsNtupler] PFRecHit ref NOT AVAILABLE. ProductID=" << pfrh_ref.id()
                    << " key=" << pfrh_ref.key()
                    << " run:lumi:event=" << iEvent.id().run() << ":" << iEvent.id().luminosityBlock() << ":" << iEvent.id().event()
                    << "\n";
          continue;
        }

        // std::cout << "[PFObjectsNtupler] PFRecHit ref AVAILABLE. ProductID=" << pfrh_ref.id()
        //           << " key=" << pfrh_ref.key() << "\n";

        const reco::PFRecHit& pfrh = *pfrh_ref;   // <-- keep ONLY THIS ONE


        DetId pfrh_did(pfrh.detId());
        if (pfrh_did.det() != DetId::Hcal) continue;

        HcalDetId pfrh_hid(pfrh.detId());
        const int pfrh_depth = pfrh_hid.depth();
        const int pfrh_ieta  = pfrh_hid.ieta();
        const int pfrh_iphi  = pfrh_hid.iphi();
        const auto pfrh_subdet = pfrh_hid.subdet();

        // --- eta/phi from geometry (cell center) ---
        float pfrh_eta = -999.f;
        float pfrh_phi = -999.f;

        const CaloCellGeometry* cell = caloGeom.getGeometry(pfrh_did);
        if (cell) {
          const GlobalPoint& pos = cell->getPosition();
          pfrh_eta = pos.eta();
          pfrh_phi = pos.phi();
        } else {
          // Fallback: compute eta/phi from ieta/iphi using HCAL mapping
          pfrh_eta = hcalEtaFromIeta(pfrh_ieta);
          pfrh_phi = hcalPhiFromIphi(pfrh_iphi);
        }

        // ---------------- depth-dependent cut ----------------
        // choose the right per-depth threshold
        //float Emin = -1.0f;

        //if (pfrh_subdet == HcalBarrel) {
          //if (pfrh_depth >= 1 && pfrh_depth <= (int)HB_depth_Emin.size())
            //Emin = HB_depth_Emin[pfrh_depth - 1];
        //} else if (pfrh_subdet == HcalEndcap) {
          //if (pfrh_depth >= 1 && pfrh_depth <= (int)HE_depth_Emin.size())
            //Emin = HE_depth_Emin[pfrh_depth - 1];
        //} else {
          //continue;  // ignore HF/HO
        //}

        // if depth is out of range, drop it
        //if (Emin < 0.0f) continue;

        // APPLY the cut on PFRecHit energy
        //if (pfrh.energy() < Emin) continue;
        // ------------------------------------

        const float pfrh_frac = hitRefAndFrac.fraction();
        // returns (eta, phi) from an HcalDetId
        // auto [pfrh_eta, pfrh_phi] = hcalEtaPhiFromDetId(pfrh_hid);

        const bool isHB = (pfrh_subdet == HcalBarrel);
        const bool isHE = (pfrh_subdet == HcalEndcap);

        const float energy = pfrh.energy();
        const float fracInCluster = pfrh_frac;
        const float energyFracInCluster = energy * fracInCluster;
        const float time = pfrh.time();

        // ALWAYS fill HBHE (combined)
        fill_pfrh(hbhe_pfrh_energy_,
                  hbhe_pfrh_energyFracInCluster_,
                  hbhe_pfrh_fracInCluster_,
                  hbhe_pfrh_eta_,
                  hbhe_pfrh_phi_,
                  hbhe_pfrh_depth_,
                  hbhe_pfrh_clusterIdx_,
                  hbhe_pfrh_ieta_,
                  hbhe_pfrh_iphi_,
                  hbhe_pfrh_time_,
                  energy,
                  energyFracInCluster,
                  fracInCluster,
                  pfrh_eta,
                  pfrh_phi,
                  pfrh_depth,
                  clusterIndex,
                  pfrh_ieta,
                  pfrh_iphi,
                  time);
        hebh_pfrh_count++;
        // additionally fill HB or HE
        if (isHB) {
          fill_pfrh(hb_pfrh_energy_,
                    hb_pfrh_energyFracInCluster_,
                    hb_pfrh_fracInCluster_,
                    hb_pfrh_eta_,
                    hb_pfrh_phi_,
                    hb_pfrh_depth_,
                    hb_pfrh_clusterIdx_,
                    hb_pfrh_ieta_,
                    hb_pfrh_iphi_,
                    hb_pfrh_time_,
                    energy,
                    energyFracInCluster,
                    fracInCluster,
                    pfrh_eta,
                    pfrh_phi,
                    pfrh_depth,
                    clusterIndex,
                    pfrh_ieta,
                    pfrh_iphi,
                    time);
          hb_pfrh_count++;
        } else if (isHE) {
          fill_pfrh(he_pfrh_energy_,
                    he_pfrh_energyFracInCluster_,
                    he_pfrh_fracInCluster_,
                    he_pfrh_eta_,
                    he_pfrh_phi_,
                    he_pfrh_depth_,
                    he_pfrh_clusterIdx_,
                    he_pfrh_ieta_,
                    he_pfrh_iphi_,
                    he_pfrh_time_,
                    energy,
                    energyFracInCluster,
                    fracInCluster,
                    pfrh_eta,
                    pfrh_phi,
                    pfrh_depth,
                    clusterIndex,
                    pfrh_ieta,
                    pfrh_iphi,
                    time);
          he_pfrh_count++;
        }
      }

      hbhe_pfrh_counts_.push_back(hebh_pfrh_count);   
      hb_pfrh_counts_.push_back(hb_pfrh_count);
      he_pfrh_counts_.push_back(he_pfrh_count);
      
    clusterIndex++;}
  }

  // PF Blocks
  edm::Handle<std::vector<reco::PFBlock>> pfBlocks;
  iEvent.getByToken(pfBlocksToken_, pfBlocks);
  if (pfBlocks.isValid()) {
    num_pfBlocks_ = pfBlocks->size();
  }
  
  // uMNio (laserType from HCAL uMNio digi; -1000 if not present)
  edm::Handle<HcalUMNioDigi> cumnio;
  iEvent.getByToken(uMNioToken_, cumnio);
  if (cumnio.isValid()) {
    laserType_ = cumnio->valueUserWord(1);
  }

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFObjectsNtupler);