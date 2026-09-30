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

#include "DataFormats/DetId/interface/DetId.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/CaloRecHitAuxSetter.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/PositionVector3D.h"

#include "TTree.h"
#include <vector>

class PFObjectsNtuplertdc : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit PFObjectsNtuplertdc(const edm::ParameterSet&);
  ~PFObjectsNtuplertdc() override {}

  void analyze(const edm::Event&, const edm::EventSetup&) override;

private:
  // Tokens
  
  edm::EDGetTokenT<std::vector<reco::PFCluster>> hcalClustersToken_;
  edm::EDGetTokenT<edm::SortedCollection<HBHERecHit>> hbheRechitsToken_;
  edm::EDGetTokenT<std::vector<reco::PFBlock>> pfBlocksToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;

  void beginRun(const edm::Run&, const edm::EventSetup&);
  TTree* tree_;

  // HCAL clusters
  std::vector<float> hcal_energy_, hcal_eta_, hcal_phi_, hcal_depth_;
  // HBHE rechits associated to clusters
  std::vector<float> hbhe_rechit_eta_; std::vector<float> hbhe_rechit_phi_; 
  std::vector<float> hbhe_rechit_tdc_; std::vector<float> hbhe_rechit_energy_; //std::vector<int> hbhe_rechit_clusterIndex_; 
  std::vector<int> hbheRechit_clusterIdx_; std::vector<int> clusterIdx_;
  std::vector<float> hbhe_rechit_depth_;
  std::vector<int> hb_rechit_tdc_; std::vector<int> he_rechit_tdc_;
  std::vector<float> hb_rechit_energy_;
  std::vector<float> he_rechit_energy_;
  std::vector<float> hb_rechit_depth_; std::vector<float> he_rechit_depth_;
  std::vector<int> hb_rechit_clusterIdx_;
  std::vector<int> he_rechit_clusterIdx_;

  // PF Blocks (just store number of elements for now)
  int num_pfBlocks_;
};

PFObjectsNtuplertdc::PFObjectsNtuplertdc(const edm::ParameterSet& iConfig)
  : caloGeometryToken_{esConsumes<CaloGeometry, CaloGeometryRecord>()}
{
  usesResource("TFileService");

  hcalClustersToken_ = consumes<std::vector<reco::PFCluster>>(iConfig.getParameter<edm::InputTag>("hcalClusters"));
  hbheRechitsToken_ = consumes<edm::SortedCollection<HBHERecHit>>(edm::InputTag("hbhereco", "", "ReRECO")); // make sure this matches the input file! 
  pfBlocksToken_ = consumes<std::vector<reco::PFBlock>>(iConfig.getParameter<edm::InputTag>("pfBlocks"));

  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("pfTree", "PF objects");

  // HCAL cluster branches
  tree_->Branch("hcal_energy", &hcal_energy_);
  tree_->Branch("hcal_eta", &hcal_eta_);
  tree_->Branch("hcal_phi", &hcal_phi_);
  tree_->Branch("hcal_depth", &hcal_depth_);
  tree_->Branch("clusterIdx", &clusterIdx_);

  // HBHE rechit branches
  tree_->Branch("hbhe_rechit_energy", &hbhe_rechit_energy_);
  tree_->Branch("hbhe_rechit_eta", &hbhe_rechit_eta_);
  tree_->Branch("hbhe_rechit_phi", &hbhe_rechit_phi_);
  tree_->Branch("hbhe_rechit_tdc", &hbhe_rechit_tdc_);
  tree_->Branch("hbhe_rechit_depth", &hbhe_rechit_depth_);
  tree_->Branch("hbheRechit_clusterIdx", &hbheRechit_clusterIdx_);


  // Separate HB and HE rechit branches
  tree_->Branch("hb_rechit_tdc", &hb_rechit_tdc_);
  tree_->Branch("he_rechit_tdc", &he_rechit_tdc_);
  tree_->Branch("hb_rechit_energy", &hb_rechit_energy_);  
  tree_->Branch("he_rechit_energy", &he_rechit_energy_);
  tree_->Branch("hb_rechit_depth", &hb_rechit_depth_);    
  tree_->Branch("he_rechit_depth", &he_rechit_depth_);
  tree_->Branch("hb_rechit_clusterIdx", &hb_rechit_clusterIdx_);
  tree_->Branch("he_rechit_clusterIdx", &he_rechit_clusterIdx_);

  // PF block info
  tree_->Branch("num_pfBlocks", &num_pfBlocks_);
}

// Convert ieta to eta using HCAL mapping

static double hcalEtaFromIeta(int ieta) {
    // HB: |ieta| <= 16, HE: 17 <= |ieta| <= 28
    int sign = (ieta >= 0 ? 1 : -1);
    int absi = std::abs(ieta);

    double eta = 0.0;

    if (absi <= 16) { // HB
        eta = 0.087 * (absi - 0.5);
    } else if (absi <= 28) { // HE
        eta = 0.087 * (16 + (absi - 16) * 0.9); // approximate
    } else { // HF etc; nothing should be outside really
        eta = 40.0; //placeholder
    }
    return sign * eta;
}

// Convert iphi to phi (HB/HE have 72 phi bins)
static double hcalPhiFromIphi(int iphi) {
    // HCAL iphi runs from 1..72
    double phi = (iphi - 1) * (M_PI / 36.0); // 2pi/72
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

void PFObjectsNtuplertdc::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // Clear all vectors

  hcal_energy_.clear(); hcal_eta_.clear(); hcal_phi_.clear(); hcal_depth_.clear();
  clusterIdx_.clear();

  hbhe_rechit_energy_.clear(); hbhe_rechit_eta_.clear();
  hbhe_rechit_phi_.clear(); hbhe_rechit_depth_.clear(); 
  hbhe_rechit_tdc_.clear();
  hbheRechit_clusterIdx_.clear();

  hb_rechit_tdc_.clear(); he_rechit_tdc_.clear();
  hb_rechit_depth_.clear();
  he_rechit_depth_.clear();
  hb_rechit_energy_.clear();
  he_rechit_energy_.clear();
  hb_rechit_clusterIdx_.clear();  
  he_rechit_clusterIdx_.clear();

  num_pfBlocks_ = 0;

  // HCAL Clusters
  edm::Handle<std::vector<reco::PFCluster>> hcalClusters;
  iEvent.getByToken(hcalClustersToken_, hcalClusters);

  // HBHE RecHits (raw rechits)
  edm::Handle<edm::SortedCollection<HBHERecHit>> hbheRechits;
  iEvent.getByToken(hbheRechitsToken_, hbheRechits);
  //Done: add a counter to count the numbers of clusters we have looped through, use them as cluster indices (0 index)
  int clusterIndex = 0;
  if (hcalClusters.isValid()) {
    for (const auto& cl : *hcalClusters) {
      // Save the cluster info
      // Done: save the cluster index
      clusterIdx_.push_back(clusterIndex);
      hcal_energy_.push_back(cl.energy());
      hcal_eta_.push_back(cl.eta());
      hcal_phi_.push_back(cl.phi());
      hcal_depth_.push_back(cl.depth());

      // std::cout << "Cluster eta: " << cl.eta() << " phi: " << cl.phi() << std::endl;

      // Loop over HBHE rechits associated to this HCAL cluster. Stored as hbhereco, these are the raw ones instead of PF (since that was a transitory collection)
      if (hbheRechits.isValid()) {

        
        for (const auto& rh : *hbheRechits) {

          // Get position info from HBHE rechit geometry
          // Takes the HCAL rechit’s detector ID and turns it into an HcalDetId object
          //rh.id() = raw DetId of the cell that the rechit corresponds to

          HcalDetId detid = rh.id();
          // std::cout << "Hit energy: " << rh.energy()
          //     << " detId: " << detid.rawId()
          //     << " depth: " << detid.depth() << std::endl;

          auto [rh_eta, rh_phi] = hcalEtaPhiFromDetId(detid);

          if (std::abs(rh_eta - cl.eta()) > 0.4) continue;
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.2) continue;

          // Save the rechit info
          hbhe_rechit_energy_.push_back(rh.energy());
          hbhe_rechit_eta_.push_back(rh_eta);
          hbhe_rechit_phi_.push_back(rh_phi);
          hbhe_rechit_depth_.push_back(detid.depth());

          // Getting TDC value from auxTDC field
          int six_bits_mask = 0x3f;  // 6-bit mask
          int ts = 3;                // TS3 is SOI
          int SOI_TDC = CaloRecHitAuxSetter::getField(rh.auxTDC(), six_bits_mask, ts * 6);
          hbhe_rechit_tdc_.push_back(SOI_TDC);  //  pushback TDC value

          // Saving HB and HE info separately
          if ( (detid.depth() == 4 && detid.ietaAbs() == 16) || (detid.ietaAbs() > 16) ) {
            // HE rechit
            he_rechit_tdc_.push_back(SOI_TDC);
            he_rechit_depth_.push_back(detid.depth());
            he_rechit_energy_.push_back(rh.energy());
            he_rechit_clusterIdx_.push_back(clusterIndex);

          } else {
            // HB rechit
            hb_rechit_tdc_.push_back(SOI_TDC);
            hb_rechit_depth_.push_back(detid.depth());
            hb_rechit_energy_.push_back(rh.energy());
            hb_rechit_clusterIdx_.push_back(clusterIndex);
          }
          
          hbheRechit_clusterIdx_.push_back(clusterIndex); // save cluster index association so it is possible to map backwards to the cluster this rechit was near);
          //hbhe_rechit_clusterIndex_.push_back(hcal_energy_.size() - 1); // save cluster index association so it is possible to map backwards to the cluster this rechit was near
        }
      }
    clusterIndex++;}
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
DEFINE_FWK_MODULE(PFObjectsNtuplertdc);