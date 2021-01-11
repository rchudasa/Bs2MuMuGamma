#ifndef _MiniAnalyzer_h
#define _MiniAnalyzer_h

// system include files
#include <memory>

// user include files
#include <map>
#include <string>

#include "TH1.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/Handle.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // muy importante para MiniAOD

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/V0Candidate/interface/V0Candidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerObjectMapRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"

#include "RecoVertex/VertexPrimitives/interface/BasicSingleVertexState.h"
#include "RecoVertex/VertexPrimitives/interface/VertexState.h"

#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

//
// class decleration
//

class MiniAnalyzer : public edm::EDAnalyzer {
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  void fillPsi(const reco::Candidate& genpsi);
  void fillV0(const reco::Candidate& genv0);
  //int const getMuCat(reco::Muon const& muon) const;
  bool IsTheSame(const pat::GenericParticle& tk, const pat::Muon& mu);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  void printout(const RefCountedKinematicVertex& myVertex) const;
  void printout(const RefCountedKinematicParticle& myParticle) const;
  void printout(const RefCountedKinematicTree& myTree) const;

  void clearVariables();
  void saveGenInfo(const edm::Event&);
  bool hasBeamSpot(const edm::Event&);
  bool hasPrimaryVertex(const edm::Event &);

  //void MatchMuonWithTriggers(const pat::Muon &iMuon, const std::vector<std::string>& TrigList, std::string &TrigListNameTmp);
  void CheckHLTTriggers(const std::vector<std::string>& TrigList);


  // ----------member data ---------------------------
  
  edm::EDGetTokenT<edm::View<pat::Muon>> dimuon_label;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trakCollection_label;
  edm::EDGetTokenT<reco::VertexCollection> primaryVertices_label;
  edm::EDGetTokenT<edm::TriggerResults> triggerResults_label;
  edm::EDGetTokenT<reco::BeamSpot> BS_label;

  edm::EDGetTokenT<edm::View<reco::GenParticle>>      pruned_label;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> genColl_label;
  edm::EDGetTokenT<pat::PhotonCollection>  photon_label;

  std::map<std::string,TH1F*> histContainer_; 


  //--------------------                                                                                                                    
  // Across the event                                                                 
  //--------------------                                                                                               
  reco::BeamSpot beamSpot_;
  reco::Vertex primaryVertex_;

  TTree*      tree_;

 
  // variables associated with tree branches
  UInt_t         run_;
  ULong64_t      event_;
  UInt_t         lumis_;
  Bool_t         isData_;
  Int_t          nPV_;
  Int_t          nSlimmedSecondV_;
  Int_t          nprivtx;
  
  std::vector<float>  pvX_;
  std::vector<float>  pvY_;
  std::vector<float>  pvZ_;
  
  // gen variables
  Int_t dimuon_pdgId;
  TLorentzVector gen_dimuon_p4;
  TLorentzVector gen_muonP_p4;
  TLorentzVector gen_muonM_p4;
  TLorentzVector gen_photon_p4;
  
  std::vector<float> genbPt_;
  std::vector<float> genbEta_;
  std::vector<float> genbPhi_;
  std::vector<float> genbMass_;
  std::vector<int>   genbPID_;
  
  std::vector<float> genbPx_;
  std::vector<float> genbPy_;
  std::vector<float> genbPz_;
  

  std::vector<float> genphoPt_;
  std::vector<float> genphoEta_;
  std::vector<float> genphoPhi_;
  std::vector<float> genphoPx_;
  std::vector<float> genphoPy_;
  std::vector<float> genphoPz_;
  std::vector<float> genphoVtxx_;
  std::vector<float> genphoVtxy_;
  std::vector<float> genphoVtxz_;

  std::vector<float> genmumPt_;
  std::vector<float> genmumEta_;
  std::vector<float> genmumPhi_;
  std::vector<float> genmupPt_;
  std::vector<float> genmupEta_;
  std::vector<float> genmupPhi_;
  

  // reco::Photon
  Int_t          nPho_;
  std::vector<float>  phoE_;
  std::vector<float>  phoEt_;
  std::vector<float>  phoEta_;
  std::vector<float>  phoPhi_;
  
  
  // reco::Muon
  Int_t          nMu_;
  std::vector<float>  muPt_;
  std::vector<float>  muEta_;
  std::vector<float>  muPhi_;
  std::vector<int>    muCharge_;
  std::vector<int>    muType_;
  std::vector<int>    muIsGood_;

  std::vector<float>  mumuPt_;
  std::vector<float>  mumuRapidity_;
  std::vector<float>  mumuMass_; 

  int mupCategory;
  int mumCategory;
  int mupME1Clean;
  int mumME1Clean;
  
  std::vector<float>       *mumC2;
  std::vector<int>         *mumNHits, *mumNPHits; 
  std::vector<float>       *mupC2;
  std::vector<int>         *mupNHits, *mupNPHits;
  std::vector<float>       *mumdxy, *mupdxy, *mumdz, *mupdz;
  std::vector<float>       *muon_dca;

  std::vector<int>         *tri_Dim25, *tri_Dim20, *tri_JpsiTk; 
  
  std::vector<bool>        *mu1soft, *mu2soft, *mu1tight, *mu2tight;  
  std::vector<bool>        *mu1PF, *mu2PF, *mu1loose, *mu2loose;  
 
  int                      muAcc, muTrig, weight;
  // *************************************
  unsigned int             nB;
  unsigned int             nMu;
    
  std::vector<float>       *B_mass, *B_px, *B_py, *B_pz, *B_charge;
  std::vector<float>       *B_k_px, *B_k_py, *B_k_pz,  *B_k_charge1; 
  std::vector<float>       *B_k_px_track, *B_k_py_track, *B_k_pz_track;

  std::vector<float>       *B_J_mass, *B_J_px, *B_J_py, *B_J_pz;

  std::vector<float>       *B_J_pt1, *B_J_px1, *B_J_py1, *B_J_pz1;
  std::vector<float>       *B_J_pt2, *B_J_px2, *B_J_py2, *B_J_pz2;
  std::vector<int>         *B_J_charge1, *B_J_charge2;

  // Primary Vertex (PV)
  unsigned int             nVtx;
  float                    priVtxX, priVtxY, priVtxZ, priVtxXE, priVtxYE, priVtxZE, priVtxCL;
  float                    priVtxXYE, priVtxXZE, priVtxYZE;
  
  // ********************************** ************************************************************************
 
  std::vector<float>       *B_chi2, *B_J_chi2;
  std::vector<float>       *B_Prob, *B_J_Prob;

  std::vector<float>       *B_DecayVtxX,  *B_DecayVtxY,  *B_DecayVtxZ;
  std::vector<double>      *B_DecayVtxXE, *B_DecayVtxYE, *B_DecayVtxZE;
  std::vector<double>      *B_DecayVtxXYE, *B_DecayVtxXZE, *B_DecayVtxYZE;



};
#endif
