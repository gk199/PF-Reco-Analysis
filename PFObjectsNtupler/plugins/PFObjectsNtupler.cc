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
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/PositionVector3D.h"

#include "TTree.h"
#include <vector>

class PFObjectsNtupler : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PFObjectsNtupler(const edm::ParameterSet&);
  ~PFObjectsNtupler() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  edm::EDGetTokenT<std::vector<reco::PFCandidate>> pfCandidatesToken_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> ecalClustersToken_;
  edm::EDGetTokenT<std::vector<reco::PFCluster>> hcalClustersToken_;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit>> hbheRechitsToken_;
  edm::EDGetTokenT<std::vector<reco::PFBlock>> pfBlocksToken_;

  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geomToken_;
  const CaloGeometry* geometry_ = nullptr;

  void beginRun(const edm::Run&, const edm::EventSetup&);

  // Output tree and variables
  TTree* tree_;

  // PF Candidates
  std::vector<float> pf_pt_, pf_eta_, pf_phi_, pf_energy_;
  std::vector<int> pf_charge_, pf_pdgId_;

  // ECAL clusters
  std::vector<float> ecal_energy_, ecal_eta_, ecal_phi_, ecal_time_;

  // HCAL clusters
  std::vector<float> hcal_energy_, hcal_eta_, hcal_phi_, hcal_time_, hcal_depth_;
  // HBHE rechits associated to clusters
  std::vector<float> hbhe_rechit_energy_; std::vector<float> hbhe_rechit_eta_; std::vector<float> hbhe_rechit_phi_; std::vector<float> hbhe_rechit_depth_; std::vector<float> hbhe_rechit_time_; std::vector<int> hbhe_rechit_clusterIndex_; 

  // PF Blocks (just store number of elements for now)
  int num_pfBlocks_;
};

PFObjectsNtupler::PFObjectsNtupler(const edm::ParameterSet& iConfig)
{
  usesResource("TFileService");
  pfCandidatesToken_ = consumes<std::vector<reco::PFCandidate>>(iConfig.getParameter<edm::InputTag>("pfCandidates"));
  ecalClustersToken_ = consumes<std::vector<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("ecalClusters"));
  hcalClustersToken_ = consumes<std::vector<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("hcalClusters"));
  hbheRechitsToken_ = consumes<edm::SortedCollection<HBHERecHit>>(edm::InputTag("hbhereco", "", "ReRECO")); // make sure this matches the input file! 
  // hbheRechitsToken_ = consumes<std::vector<reco::PFRecHit>>(edm::InputTag("particleFlowRecHitHBHE", "Cleaned", "ReRECOtoAOD"));
  // hbheRechitsToken_ = consumes<std::vector<reco::PFRecHit>>(edm::InputTag("particleFlowRecHitHBHE", "", "ReRECO"));
  pfBlocksToken_ = consumes<std::vector<reco::PFBlock>>(iConfig.getParameter<edm::InputTag>("pfBlocks"));

  geomToken_ = esConsumes<CaloGeometry, CaloGeometryRecord>();

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("pfTree", "PF objects");

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

  // HCAL cluster branches
  tree_->Branch("hcal_energy", &hcal_energy_);
  tree_->Branch("hcal_eta", &hcal_eta_);
  tree_->Branch("hcal_phi", &hcal_phi_);
  tree_->Branch("hcal_time", &hcal_time_);
  tree_->Branch("hcal_depth", &hcal_depth_);
  // HBHE rechit branches
  tree_->Branch("hbhe_rechit_energy", &hbhe_rechit_energy_);
  tree_->Branch("hbhe_rechit_eta", &hbhe_rechit_eta_);
  tree_->Branch("hbhe_rechit_phi", &hbhe_rechit_phi_);
  tree_->Branch("hbhe_rechit_depth", &hbhe_rechit_depth_);
  tree_->Branch("hbhe_rechit_time", &hbhe_rechit_time_);
  tree_->Branch("hbhe_rechit_clusterIndex", &hbhe_rechit_clusterIndex_);

  // PF block info
  tree_->Branch("num_pfBlocks", &num_pfBlocks_);
}

// do this to load the geometry (needed for hbhe rechits)
void PFObjectsNtupler::beginRun(const edm::Run&, const edm::EventSetup& iSetup) {
  geometry_ = &iSetup.getData(geomToken_);
}

// Approximate conversion from ieta/iphi to eta/phi
std::pair<double,double> detidEtaPhiFallback(const HcalDetId& detid) {
    double phi = 0.0;
    double eta = 0.0;

    // iPhi runs from 1..72 (HB/HE)
    int iphi = detid.iphi();
    phi = (iphi - 0.5) * 2.0 * M_PI / 72.0;  // radians

    // iEta runs from 1..16 for HB, 16..29 for HE (example)
    int ieta = detid.ieta();
    bool isNegative = ieta < 0;
    int absIeta = std::abs(ieta);

    // HB: 1..16, HE: 16..29 (simplified central value for eta)
    if (absIeta <= 16) {
        // HB approximate center of tower
        eta = 0.087 * (absIeta - 0.5);
    } else {
        // HE approximate center of tower
        eta = 0.087 * (16 - 0.5) + 0.09 * (absIeta - 16 + 0.5);  // rough spacing
    }

    if (isNegative) eta = -eta;

    return std::make_pair(eta, phi);
}

void PFObjectsNtupler::analyze(const edm::Event& iEvent, const edm::EventSetup&)
{
  // Clear all vectors
  pf_pt_.clear(); pf_eta_.clear(); pf_phi_.clear(); pf_energy_.clear(); pf_charge_.clear(); pf_pdgId_.clear();
  ecal_energy_.clear(); ecal_eta_.clear(); ecal_phi_.clear(); ecal_time_.clear();
  hcal_energy_.clear(); hcal_eta_.clear(); hcal_phi_.clear(); hcal_time_.clear(); hcal_depth_.clear();
  hbhe_rechit_energy_.clear(); hbhe_rechit_eta_.clear(); hbhe_rechit_phi_.clear(); hbhe_rechit_depth_.clear(); hbhe_rechit_time_.clear(); hbhe_rechit_clusterIndex_.clear();
  num_pfBlocks_ = 0;

  // PF Candidates
  edm::Handle<std::vector<reco::PFCandidate>> pfCandidates;
  iEvent.getByToken(pfCandidatesToken_, pfCandidates);
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

  // ECAL Clusters
  edm::Handle<std::vector<reco::PFCluster>> ecalClusters;
  iEvent.getByToken(ecalClustersToken_, ecalClusters);
  if (ecalClusters.isValid()) {
    for (const auto& cl : *ecalClusters) {
      ecal_energy_.push_back(cl.energy());
      ecal_eta_.push_back(cl.eta());
      ecal_phi_.push_back(cl.phi());
      ecal_time_.push_back(cl.time());
    }
  }

  // HCAL Clusters
  edm::Handle<std::vector<reco::PFCluster>> hcalClusters;
  iEvent.getByToken(hcalClustersToken_, hcalClusters);

  // HBHE RecHits (raw rechits)
  edm::Handle<edm::SortedCollection<HBHERecHit>> hbheRechits;
  iEvent.getByToken(hbheRechitsToken_, hbheRechits);

  if (hcalClusters.isValid()) {
    for (const auto& cl : *hcalClusters) {
      hcal_energy_.push_back(cl.energy());
      hcal_eta_.push_back(cl.eta());
      hcal_phi_.push_back(cl.phi());
      hcal_time_.push_back(cl.time());
      hcal_depth_.push_back(cl.depth());

      // std::cout << "Cluster eta: " << cl.eta() << " phi: " << cl.phi() << std::endl;

      // Loop over HBHE rechits associated to this HCAL cluster. Stored as hbhereco, these are the raw ones instead of PF (since that was a transitory collection)
      if (hbheRechits.isValid()) {
        for (const auto& rh : *hbheRechits) {
          // Get position info from HBHE rechit geometry
          HcalDetId detid = rh.id();
          // std::cout << "Hit energy: " << rh.energy()
          //     << " detId: " << detid.rawId()
          //     << " depth: " << detid.depth() << std::endl;
          double rh_eta = 0., rh_phi = 0.;

          if (geometry_) {
              const CaloSubdetectorGeometry* subDetGeom = geometry_->getSubdetectorGeometry(detid);
              if (subDetGeom) {
                  const CaloCellGeometry* cell = subDetGeom->getGeometry(detid);
                  if (cell) {
                      const GlobalPoint pos = cell->getPosition();
                      rh_eta = pos.eta();
                      rh_phi = pos.phi();
                  } else {
                      // fallback: compute eta/phi from detid
                      auto etaPhi = detidEtaPhiFallback(detid);
                      rh_eta = etaPhi.first;
                      rh_phi = etaPhi.second;
                  }
              } else {
                  auto etaPhi = detidEtaPhiFallback(detid);
                  rh_eta = etaPhi.first;
                  rh_phi = etaPhi.second;
              }
          } else {
              auto etaPhi = detidEtaPhiFallback(detid);
              rh_eta = etaPhi.first;
              rh_phi = etaPhi.second;
          }

          // Optional: pre-cut by eta window
          if (std::abs(rh_eta - cl.eta()) > 0.4) continue;

          double dR = reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi);
          if (dR > 0.2) continue;

          // Save the rechit info
          hbhe_rechit_energy_.push_back(rh.energy());
          hbhe_rechit_eta_.push_back(rh_eta);
          hbhe_rechit_phi_.push_back(rh_phi);
          hbhe_rechit_depth_.push_back(detid.depth());
          hbhe_rechit_time_.push_back(rh.time());
          hbhe_rechit_clusterIndex_.push_back(hcal_energy_.size() - 1); // save cluster index association so it is possible to map backwards to the cluster this rechit was near
        }
      }
    }
  }

  // PF Blocks
  edm::Handle<std::vector<reco::PFBlock>> pfBlocks;
  iEvent.getByToken(pfBlocksToken_, pfBlocks);
  if (pfBlocks.isValid()) {
    num_pfBlocks_ = pfBlocks->size();
  }

  tree_->Fill();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PFObjectsNtupler);