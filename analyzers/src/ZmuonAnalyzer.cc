// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// new includes
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TTree.h>

// for vertexing                                                                                              
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"  // edm::Handle
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"                                                        
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
     
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"                                                                        
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonQuality.h"
             
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"                                                                
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/PrimaryVertexProducer/interface/VertexHigherPtSquared.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"

#include <math.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <algorithm>  // std::sort, std::swap
#include <iostream>  // std::cout, std::endl
#include <string>

class ZmuonAnalyzer : public edm::EDAnalyzer {
public:
	explicit ZmuonAnalyzer(const edm::ParameterSet&);
private:
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
   edm::EDGetTokenT<std::vector<reco::Vertex>> token_vertices;
   edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
//   edm::EDGetTokenT<pat::TriggerEvent> triggerEventToken_;
//   edm::EDGetTokenT<std::vector<std::string>> muonMatch_;
   edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
   edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

   TTree * tree;
   TTree * treemc;
 
   int truth_nMuPos, truth_nMuNeg, truth_nPhoton, truth_nBs;
   std::vector<double> truth_Bsmuon_pt, truth_Bsmuon_eta, truth_Bsmuon_phi;
   std::vector<double> truth_Bsphoton_pt, truth_Bsphoton_eta, truth_Bsphoton_phi;
   std::vector<double> truth_Bs_pt, truth_Bs_eta, truth_Bs_phi, truth_Bs_mass;

   std::vector<double> truth_Bs_pdgid, truth_Bsmuon_pdgid;

};

ZmuonAnalyzer::ZmuonAnalyzer(const edm::ParameterSet& iConfig)
{
   muonsToken_        = consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection"));
   token_vertices     = consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   triggerBits_       = consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"));
//   triggerEventToken_ = consumes<pat::TriggerEvent>(iConfig.getParameter<edm::InputTag>("triggerEvent"));
//   muonMatch_         = consumes<std::vector<std::string>>(iConfig.getParameter< std::string >( "muonMatch" ) );
   triggerObjects_ = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("objects"));
   genParticlesToken_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "tree");

   treemc = fs->make<TTree>("treemc", "treemc");
   treemc->Branch("truth_nMuPos"   ,  &truth_nMuPos);
   treemc->Branch("truth_nMuNeg"   ,  &truth_nMuNeg);
   treemc->Branch("truth_nPhoton"  ,  &truth_nPhoton);
   treemc->Branch("truth_nBs"      ,  &truth_nBs);
   treemc->Branch("truth_Bsmuon_pt",  &truth_Bsmuon_pt);
   treemc->Branch("truth_Bsmuon_eta", &truth_Bsmuon_eta);
   treemc->Branch("truth_Bsmuon_phi", &truth_Bsmuon_phi);
   treemc->Branch("truth_Bsmuon_pdgid", &truth_Bsmuon_pdgid);
   treemc->Branch("truth_Bsphoton_pt",  &truth_Bsphoton_pt);
   treemc->Branch("truth_Bsphoton_eta", &truth_Bsphoton_eta);
   treemc->Branch("truth_Bsphoton_phi", &truth_Bsphoton_phi);
   treemc->Branch("truth_Bs_pt",   &truth_Bs_pt);
   treemc->Branch("truth_Bs_eta",  &truth_Bs_eta);
   treemc->Branch("truth_Bs_phi",  &truth_Bs_phi);
   treemc->Branch("truth_Bs_mass", &truth_Bs_mass);
   treemc->Branch("truth_Bs_pdgid", &truth_Bs_pdgid);

}

void ZmuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

   bool mc = true;
   truth_nMuPos = 0; truth_nMuNeg = 0; truth_nPhoton = 0; truth_nBs = 0;
   truth_Bsmuon_pt.clear(); 
   truth_Bsmuon_eta.clear();
   truth_Bsmuon_phi.clear();
   truth_Bsmuon_pdgid.clear();
   truth_Bs_pt.clear(); truth_Bs_eta.clear(); truth_Bs_phi.clear(); truth_Bs_mass.clear();
   truth_Bsphoton_pt.clear(); 
   truth_Bsphoton_eta.clear();
   truth_Bsphoton_phi.clear();
 
   truth_Bs_pdgid.clear();  

   edm::Handle<std::vector<reco::Vertex>> reco_vertices;
   iEvent.getByToken(token_vertices, reco_vertices);

// **************
// MC starts here
// **************
   if (mc) {
    edm::Handle<reco::GenParticleCollection> mc_particles;
    iEvent.getByToken(genParticlesToken_, mc_particles);

    int BS = 531;
    int MUON = 13;
    int ANTIMUON = -13;
    int PHOTON = 22;

    bool found_mu = false;
    bool found_anti_mu = false;
    bool found_photon = false;

    for(unsigned int i = 0; i < mc_particles->size(); ++i) {
      const reco::GenParticle & gen_particle = (*mc_particles)[i];                                             
      if ( TMath::Abs(gen_particle.pdgId()) != BS ) continue; 
      int imum(-1), imup(-1), ipho(-1);
      
      // loop over all Bs daughters
      for (size_t j = 0; j < gen_particle.numberOfDaughters(); ++j) {
	const reco::Candidate  &dau = *(gen_particle.daughter(j));
	if (dau.pdgId() == MUON) imup = j;
	if (dau.pdgId() == ANTIMUON) imum = j;
	if (dau.pdgId() == PHOTON)   ipho = j;
      }

      if ( ipho == -1 ) continue;
      const reco::Candidate & pho = *(gen_particle.daughter(ipho));
      const reco::Candidate *mum = NULL;
      const reco::Candidate *mup = NULL;
      
      if (imum != -1 && imup != -1) {
	mum = gen_particle.daughter(imum);
	mup = gen_particle.daughter(imup);
	
      }
      
      if ( mum == NULL || mup == NULL) continue;
    
      //std::cout << "positive pt:" << mup->pt() << std::endl;
      truth_Bsmuon_pt .push_back(mup->pt());
      truth_Bsmuon_eta.push_back(mup->eta());
      truth_Bsmuon_phi.push_back(mup->phi());
      truth_Bsmuon_pdgid.push_back(mup->charge());

      truth_Bsmuon_pt .push_back(mum->pt());
      truth_Bsmuon_eta.push_back(mum->eta());
      truth_Bsmuon_phi.push_back(mum->phi());
      truth_Bsmuon_pdgid.push_back(mum->charge());

      truth_Bsphoton_pt .push_back(pho.pt());
      truth_Bsphoton_eta.push_back(pho.eta());
      truth_Bsphoton_phi.push_back(pho.phi());

      truth_Bs_pt  .push_back(gen_particle.pt());
      truth_Bs_eta .push_back(gen_particle.eta());
      truth_Bs_phi .push_back(gen_particle.phi());
      truth_Bs_mass.push_back(gen_particle.mass());
      truth_nBs++;
      treemc->Fill();
    } // mc particle size

    //treemc->Fill();
   } // end of mc

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZmuonAnalyzer);

