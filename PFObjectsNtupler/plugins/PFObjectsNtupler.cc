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


#include "DataFormats/Math/interface/deltaR.h"
#include "Math/GenVector/PositionVector3D.h"

#include "TTree.h"
#include "Rtypes.h"  // for UInt_t, ULong64_t
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
  edm::EDGetTokenT<EcalRecHitCollection> ebRechitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> eeRechitsToken_;
  edm::EDGetTokenT<EcalRecHitCollection> esRechitsToken_;
  edm::EDGetTokenT<std::vector<reco::PFBlock>> pfBlocksToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> caloGeometryToken_;


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
  // Ecal-barrel rechits associated to clusters
  std::vector<float> eb_rechit_energy_; std::vector<float> eb_rechit_eta_; std::vector<float> eb_rechit_phi_; 
  std::vector<float> eb_rechit_time_; std::vector<int> eb_rechit_clusterIdx_; 
  std::vector<int> eb_rechit_counts_;
  // Ecal-endcap rechits associated to clusters
  std::vector<float> ee_rechit_energy_; std::vector<float> ee_rechit_eta_; std::vector<float> ee_rechit_phi_; 
  std::vector<float> ee_rechit_time_; std::vector<int> ee_rechit_clusterIdx_; 
  std::vector<int> ee_rechit_counts_;

  // HCAL clusters
  std::vector<float> hcal_energy_, hcal_eta_, hcal_phi_, hcal_time_, hcal_depth_;
  // HBHE rechits associated to clusters
  std::vector<int> hbhe_rechit_counts_;
  std::vector<float> hbhe_rechit_energy_; std::vector<float> hbhe_rechit_eta_; std::vector<float> hbhe_rechit_phi_; 
  std::vector<float> hbhe_rechit_depth_; std::vector<float> hbhe_rechit_time_; std::vector<float> hbhe_rechit_tdc_; //std::vector<int> hbhe_rechit_clusterIndex_; 
  std::vector<int> hbheRechit_clusterIdx_; std::vector<int> clusterIdx_;

  std::vector<int> hb_rechit_counts_;
  std::vector<int> he_rechit_counts_;
  std::vector<int> hb_rechit_tdc_; std::vector<int> he_rechit_tdc_;
  std::vector<float> hb_rechit_depth_; std::vector<float> he_rechit_depth_;
  std::vector<int> hbhe_ietaAbs_;

  // PF Blocks (just store number of elements for now)
  int num_pfBlocks_;
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
  // HCAL cluster branches
  tree_->Branch("hcal_energy", &hcal_energy_);
  tree_->Branch("hcal_eta", &hcal_eta_);
  tree_->Branch("hcal_phi", &hcal_phi_);
  tree_->Branch("hcal_time", &hcal_time_);
  tree_->Branch("hcal_depth", &hcal_depth_);
  tree_->Branch("clusterIdx", &clusterIdx_);

  // HBHE rechit branches
  tree_->Branch("hbhe_rechit_energy", &hbhe_rechit_energy_);
  tree_->Branch("hbhe_rechit_eta", &hbhe_rechit_eta_);
  tree_->Branch("hbhe_rechit_phi", &hbhe_rechit_phi_);
  tree_->Branch("hbhe_rechit_depth", &hbhe_rechit_depth_);
  tree_->Branch("hbhe_rechit_time", &hbhe_rechit_time_);
  tree_->Branch("hbhe_rechit_tdc", &hbhe_rechit_tdc_);
  tree_->Branch("hbheRechit_clusterIdx", &hbheRechit_clusterIdx_);
  tree_->Branch("hbhe_rechit_ietaAbs", &hbhe_ietaAbs_);
  tree_->Branch("hb_rechit_tdc", &hb_rechit_tdc_);
  tree_->Branch("he_rechit_tdc", &he_rechit_tdc_);
  tree_->Branch("hb_rechit_depth", &hb_rechit_depth_);
  tree_->Branch("he_rechit_depth", &he_rechit_depth_);
  tree_->Branch("hb_rechit_counts", &hb_rechit_counts_);
  tree_->Branch("he_rechit_counts", &he_rechit_counts_);
  tree_->Branch("hbhe_rechit_counts", &hbhe_rechit_counts_);
  // PF block info
  tree_->Branch("num_pfBlocks", &num_pfBlocks_);
}

// Convert ieta to eta using HCAL mapping

double hcalEtaFromIeta(int ieta) {
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
double hcalPhiFromIphi(int iphi) {
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
  
  hbhe_rechit_energy_.clear(); hbhe_rechit_eta_.clear(); hbhe_rechit_phi_.clear(); hbhe_rechit_depth_.clear(); hbhe_rechit_time_.clear(); //hbhe_rechit_clusterIndex_.clear();
  hbhe_rechit_tdc_.clear();
  hbheRechit_clusterIdx_.clear(); clusterIdx_.clear();
  hb_rechit_tdc_.clear(); he_rechit_tdc_.clear();
  hbhe_ietaAbs_.clear();
  hb_rechit_depth_.clear();
  he_rechit_depth_.clear();
  hbhe_rechit_counts_.clear();
  hb_rechit_counts_.clear();  
  he_rechit_counts_.clear(); 

  num_pfBlocks_ = 0;

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
          eb_rechit_counts_.push_back(eb_count);

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
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.2) continue;

          // Save the rechit info
          eb_rechit_energy_.push_back(rh.energy());
          eb_rechit_eta_.push_back(rh_eta);
          eb_rechit_phi_.push_back(rh_phi);
          eb_rechit_time_.push_back(rh.time());
          eb_rechit_clusterIdx_.push_back(ecal_clusterIndex); // save cluster index association
        }
      }
      
      // Loop over EE rechits   
      if (eeRecHits.isValid()) {
        int ee_count = 0;

        for (const auto& rh : *eeRecHits) {
          ee_count++;
          ee_rechit_counts_.push_back(ee_count);
         
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
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.2) continue;

          // Save the rechit info
          ee_rechit_energy_.push_back(rh.energy());
          ee_rechit_eta_.push_back(rh_eta);
          ee_rechit_phi_.push_back(rh_phi); 
          ee_rechit_time_.push_back(rh.time());
          ee_rechit_clusterIdx_.push_back(ecal_clusterIndex); // save cluster index association
        }
      }
      ecal_clusterIndex++;
    }

  }

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
      hcal_time_.push_back(cl.time());
      hcal_depth_.push_back(cl.depth());

      // std::cout << "Cluster eta: " << cl.eta() << " phi: " << cl.phi() << std::endl;

      // Loop over HBHE rechits associated to this HCAL cluster. Stored as hbhereco, these are the raw ones instead of PF (since that was a transitory collection)
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
          if (reco::deltaR(cl.eta(), cl.phi(), rh_eta, rh_phi) > 0.2) continue;

          // Save the rechit info
          hbhe_rechit_energy_.push_back(rh.energy());
          hbhe_rechit_eta_.push_back(rh_eta);
          hbhe_rechit_phi_.push_back(rh_phi);
          hbhe_rechit_depth_.push_back(detid.depth());

          // Done: get and save the absolute ieta value
          int ietaAbs = detid.ietaAbs();
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

          } else {
            // HB rechit
            hb_rechit_tdc_.push_back(SOI_TDC);
            hb_rechit_depth_.push_back(detid.depth());
            hb_count++;
            hb_rechit_counts_.push_back(hb_count);
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
DEFINE_FWK_MODULE(PFObjectsNtupler);