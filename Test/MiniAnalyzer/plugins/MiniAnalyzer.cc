#include <map>
#include <string>

#include "TH1.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistance3D.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"

//For kinematic fit:
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTSFactory.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h" // most importante para MiniAOD

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



#include "TLorentzVector.h"
#include <TTree.h>



const int MUONMINUS_PDG_ID = 13;
const int KAONPLUS_PDG_ID = 321;
const int PHI_PDG_ID = 333;      // phi(1020)
const int BS_PDG_ID = 531;
const int JPSI_PDG_ID = 443;
const int PSI2S_PDG_ID = 100443;

const double PI = 3.141592653589793;


class MiniAnalyzer : public edm::EDAnalyzer {
  
public:
  explicit MiniAnalyzer(const edm::ParameterSet&);
  ~MiniAnalyzer();
  bool isAncestor(const reco::Candidate * ancestor, const reco::Candidate * particle);
  
private:
  
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void clearVariables();
  void calLS (double, double, double, double, double, double, double,
              double, double,  double, double, double, double, double,
              double, double, double, double, double*, double*);

  void calCosAlpha (double, double, double, double, double,
                    double, double, double, double, double,
                    double, double, double, double,
                    double, double, double, double,
                    double*, double*);
  void saveGenInfo(const edm::Event&);
  bool hasBeamSpot(const edm::Event&);
  bool hasPrimaryVertex(const edm::Event &);
  
  // simple map to contain all histograms; 
  // histograms are booked in the beginJob() 
  // method
  std::map<std::string,TH1F*> histContainer_; 
  
  
  //----------------------                              
  // particle properties                                                                                                
  //----------------------                                                                                                      
  ParticleMass MuonMass_;
  float MuonMassErr_;
  double BsMass_;

  
  // ----------member data ---------------------------    
  edm::EDGetTokenT<reco::BeamSpot> bsCollToken; 
  edm::EDGetTokenT<edm::View<reco::GenParticle>>      prunedCollToken;
  edm::EDGetTokenT<edm::View<pat::PackedGenParticle>> genCollToken;
  edm::EDGetTokenT<reco::VertexCollection> vtxCollToken;
  edm::EDGetTokenT<edm::View<pat::Muon>> muonCollToken; 
  //edm::EDGetTokenT<edm::View<pat::MuonCollection>>    muonCollToken;
  edm::EDGetTokenT<pat::PhotonCollection>  photonCollToken;
  edm::EDGetTokenT<edm::View<pat::PackedCandidate>> trkCollToken;
  
  // input tags 
  edm::InputTag bsSrc_;
  edm::InputTag prunedSrc_;
  edm::InputTag genSrc_;
  edm::InputTag vtxSrc_;
  edm::InputTag muonSrc_;
  edm::InputTag photonSrc_;
  edm::InputTag trkSrc_;

  //---------------------                                                                                                   
  // pre-selection cuts                                                                                                  
  //---------------------                                                                                                 
  double MuonMinPt_;
  double MuonMaxEta_;
  double MuonMaxDcaBs_;

  double MuMuMaxDca_;
  double MuMuMinVtxCl_;
  double MuMuMinLxySigmaBs_;
  double MuMuMinCosAlphaBs_;

  double BsMinVtxCl_;
  double BsMinMass_;
  double BsMaxMass_;


  //--------------------                                                                                                                    
  // Across the event                                                                 
  //--------------------                                                                                               
  reco::BeamSpot beamSpot_;
  edm::ESHandle<MagneticField> bFieldHandle_;
  reco::Vertex primaryVertex_;

  TTree*         tree_;
  
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
 
  Int_t              ngenBs_; 
  std::vector<float> genbPt_;
  std::vector<float> genbEta_;
  std::vector<float> genbRap_;
  std::vector<float> genbPhi_;
  std::vector<float> genbMass_;
  std::vector<int>   genbPID_;
  
  std::vector<float> genbPx_;
  std::vector<float> genbPy_;
  std::vector<float> genbPz_;
  
  Int_t              ngenPhoton_; 
  std::vector<float> genphoPt_;
  std::vector<float> genphoEta_;
  std::vector<float> genphoRap_;
  std::vector<float> genphoPhi_;
  std::vector<float> genphoPx_;
  std::vector<float> genphoPy_;
  std::vector<float> genphoPz_;
  std::vector<float> genphoVtxx_;
  std::vector<float> genphoVtxy_;
  std::vector<float> genphoVtxz_;

  Int_t              ngenMMu_; 
  std::vector<float> genmumPt_;
  std::vector<float> genmumEta_;
  std::vector<float> genmumRap_;
  std::vector<float> genmumPhi_;

  Int_t              ngenPMu_; 
  std::vector<float> genmupPt_;
  std::vector<float> genmupEta_;
  std::vector<float> genmupRap_;
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

  std::vector<float>  mumuLSBS_;
  std::vector<float>  mumuLSBSErr_;
  std::vector<float>  mumuDCA_;
  std::vector<float>  mumuVtxCl_;
};


MiniAnalyzer::MiniAnalyzer(const edm::ParameterSet& iConfig):

  histContainer_(),

  //IsMonteCarlo_(iConfig.getUntrackedParameter<bool>("IsMonteCarlo")),
 // particle properties                            
  
            
  MuonMass_(iConfig.getUntrackedParameter<double>("MuonMass")),
  MuonMassErr_(iConfig.getUntrackedParameter<double>("MuonMassErr")),

  BsMass_(iConfig.getUntrackedParameter<double>("BsMass")),
  
  bsSrc_(iConfig.getUntrackedParameter<edm::InputTag>("bsSrc")),
  prunedSrc_(iConfig.getUntrackedParameter<edm::InputTag>("prunedSrc")),
  genSrc_(iConfig.getUntrackedParameter<edm::InputTag>("genSrc")),
  vtxSrc_(iConfig.getUntrackedParameter<edm::InputTag>("vtxSrc")),
  muonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("muonSrc")),
  photonSrc_(iConfig.getUntrackedParameter<edm::InputTag>("photonSrc")),
  trkSrc_(iConfig.getUntrackedParameter<edm::InputTag>("trkSrc")),

 // pre-selection cuts                                                                                                                 
  MuonMinPt_(iConfig.getUntrackedParameter<double>("MuonMinPt")),
  MuonMaxEta_(iConfig.getUntrackedParameter<double>("MuonMaxEta")),
  MuonMaxDcaBs_(iConfig.getUntrackedParameter<double>("MuonMaxDcaBs")),


  MuMuMaxDca_(iConfig.getUntrackedParameter<double>("MuMuMaxDca")),
  MuMuMinVtxCl_(iConfig.getUntrackedParameter<double>("MuMuMinVtxCl")),
  MuMuMinLxySigmaBs_(iConfig.getUntrackedParameter<double>("MuMuMinLxySigmaBs")),
  MuMuMinCosAlphaBs_(iConfig.getUntrackedParameter<double>("MuMuMinCosAlphaBs")),


  BsMinVtxCl_(iConfig.getUntrackedParameter<double>("BsMinVtxCl")),
  BsMinMass_(iConfig.getUntrackedParameter<double>("BsMinMass")),
  BsMaxMass_(iConfig.getUntrackedParameter<double>("BsMaxMass")){
  
  bsCollToken     = consumes<reco::BeamSpot>(bsSrc_);
  prunedCollToken = consumes<edm::View<reco::GenParticle>>(prunedSrc_);
  genCollToken    = consumes<edm::View<pat::PackedGenParticle>>(genSrc_);
  vtxCollToken    = consumes<reco::VertexCollection>(vtxSrc_);
  muonCollToken   = consumes<edm::View<pat::Muon>>(muonSrc_);
  //muonCollToken   = consumes<edm::View<pat::MuonCollection>>(muonSrc_);
  photonCollToken = consumes<pat::PhotonCollection>(photonSrc_);
  trkCollToken    = consumes<edm::View<pat::PackedCandidate>>(trkSrc_);
  

  // initialize output TTree
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("EventTree", "Event data");
  
  tree_->Branch("run",    &run_);
  tree_->Branch("event",  &event_);
  tree_->Branch("lumis",  &lumis_);
  tree_->Branch("isData", &isData_);
  
  tree_->Branch("nPV",             &nPV_);
  tree_->Branch("nSlimmedSecondV", &nSlimmedSecondV_);
  tree_->Branch("nprivtx",         &nprivtx);
  tree_->Branch("pvX",             &pvX_);
  tree_->Branch("pvY",             &pvY_);
  tree_->Branch("pvZ",             &pvZ_);
  

  tree_->Branch("ngenBs",       &ngenBs_);
  tree_->Branch("genbPt",       &genbPt_);
  tree_->Branch("genbEta",      &genbEta_);
  tree_->Branch("genbRap",      &genbRap_);
  tree_->Branch("genbPhi",      &genbPhi_);
  tree_->Branch("genbMass",     &genbMass_);
  tree_->Branch("genbPID",      &genbPID_);

  tree_->Branch("genbPx",      &genbPx_);
  tree_->Branch("genbPy",      &genbPy_);
  tree_->Branch("genbPz",      &genbPz_);

  tree_->Branch("ngenPhoton",   &ngenPhoton_);
  tree_->Branch("genphoPt",     &genphoPt_);
  tree_->Branch("genphoEta",    &genphoEta_);
  tree_->Branch("genphoRap",    &genphoRap_);
  tree_->Branch("genphoPhi",    &genphoPhi_);

  tree_->Branch("genphoPx",    &genphoPx_);
  tree_->Branch("genphoPy",    &genphoPy_);
  tree_->Branch("genphoPz",    &genphoPz_);
  tree_->Branch("genphoVtxx",  &genphoVtxx_);
  tree_->Branch("genphoVtxy",  &genphoVtxy_);
  tree_->Branch("genphoVtxz",  &genphoVtxz_);

  tree_->Branch("ngenMMu",     &ngenMMu_);
  tree_->Branch("genmumPt",    &genmumPt_);
  tree_->Branch("genmumEta",   &genmumEta_);
  tree_->Branch("genmumRap",   &genmumRap_);
  tree_->Branch("genmumPhi",   &genmumPhi_);

  tree_->Branch("ngenPMu",     &ngenPMu_);
  tree_->Branch("genmupPt",    &genmupPt_);
  tree_->Branch("genmupEta",   &genmupEta_);
  tree_->Branch("genmupRap",   &genmupRap_);
  tree_->Branch("genmupPhi",   &genmupPhi_);
  
 
  tree_->Branch("nPho",                  &nPho_);
  tree_->Branch("phoE",                  &phoE_);
  tree_->Branch("phoEt",                 &phoEt_);
  tree_->Branch("phoEta",                &phoEta_);
  tree_->Branch("phoPhi",                &phoPhi_);
  
  tree_->Branch("nMu",                   &nMu_);
  tree_->Branch("muPt",                  &muPt_);
  tree_->Branch("muEta",                 &muEta_);
  tree_->Branch("muPhi",                 &muPhi_);
  tree_->Branch("muCharge",              &muCharge_);
  tree_->Branch("muType",                &muType_);
  tree_->Branch("muIsGood",              &muIsGood_);
  
  tree_->Branch("dimuon_pdgId",  &dimuon_pdgId,     "dimuon_pdgId/I");
  tree_->Branch("mumuPt",                &mumuPt_);
  tree_->Branch("mumuRapidity",          &mumuRapidity_);
  tree_->Branch("mumuMass",              &mumuMass_);
  tree_->Branch("mumuLSBS",              &mumuLSBS_);
  tree_->Branch("mumuLSBSErr",           &mumuLSBSErr_);
  tree_->Branch("mumuDCA",               &mumuDCA_);
  tree_->Branch("mumuVtxCl",             &mumuVtxCl_);
  
}

MiniAnalyzer::~MiniAnalyzer(){
}

void MiniAnalyzer::clearVariables(){
  run_    = 0;
  event_  = 0; 
  lumis_  = 0;
  isData_ = false; 

  nPV_ = 0;
  nSlimmedSecondV_ = 0;
  nprivtx = 0;
  pvX_                  .clear();
  pvY_                  .clear();
  pvZ_                  .clear();

  ngenBs_ = 0;
  genbPt_      .clear();
  genbEta_     .clear();
  genbRap_     .clear();
  genbPhi_     .clear();
  genbMass_    .clear();
  genbPID_     .clear();

  genbPx_      .clear();
  genbPy_      .clear();
  genbPz_      .clear();
 
  ngenPhoton_ = 0; 
  genphoPt_    .clear();
  genphoEta_   .clear();
  genphoRap_   .clear();
  genphoPhi_   .clear();

  genphoPx_    .clear();
  genphoPy_    .clear();
  genphoPz_    .clear();
  genphoVtxx_  .clear();
  genphoVtxy_  .clear();
  genphoVtxz_  .clear();

  ngenMMu_ = 0;   
  genmumPt_      .clear();
  genmumEta_     .clear();
  genmumRap_     .clear();
  genmumPhi_     .clear();

  ngenPMu_ = 0;   
  genmupPt_      .clear();
  genmupEta_     .clear();
  genmupRap_     .clear();
  genmupPhi_     .clear();
  


  nPho_ = 0;
  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
  
  nMu_ = 0;
  muPt_                 .clear();
  muEta_                .clear();
  muPhi_                .clear();
  muCharge_             .clear();
  muType_               .clear();
  muIsGood_             .clear();

  mumuPt_               .clear();
  mumuRapidity_         .clear();
  mumuMass_             .clear();

  mumuLSBS_             .clear();
  mumuLSBSErr_          .clear();
  mumuDCA_              .clear();
  mumuVtxCl_            .clear();

}

void
MiniAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup){
 
  using std::vector;
  using namespace edm;
  using namespace reco;
  using namespace std;

  clearVariables();
  double MuMuLSBS, MuMuLSBSErr;

  run_    = iEvent.id().run();
  event_  = iEvent.id().event();
  lumis_  = iEvent.luminosityBlock();
  isData_ = iEvent.isRealData();

  // save gen level information
  saveGenInfo(iEvent);
  
  // get primary vertex collection 
  edm::Handle< std::vector<reco::Vertex>> vtxs;
  iEvent.getByToken(vtxCollToken, vtxs);
  
  for (auto vt = vtxs->cbegin(); vt != vtxs->cend(); ++vt) {
    if(!vt->isFake()){
      pvX_        .push_back(vt->x());
      pvY_        .push_back(vt->y());
      pvZ_        .push_back(vt->z());
      nPV_++;
    }
  }
  
  
  // get pat muon collection 
  // Kinematic fit
  edm::ESHandle<TransientTrackBuilder> theB; 
  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",theB); 

  edm::Handle< View<pat::Muon> > muons;
  iEvent.getByToken(muonCollToken, muons);
  
  
  
  // fill pat muon- histograms
  for(View<pat::Muon>::const_iterator iMuon1 = muons->begin(); iMuon1 != muons->end(); ++iMuon1){
    for(View<pat::Muon>::const_iterator iMuon2 = iMuon1+1; iMuon2 != muons->end(); ++iMuon2) {

      if(iMuon1==iMuon2) continue;
      if( muons->size() < 2 || iMuon1->pt() < MuonMinPt_ || iMuon2->pt() < MuonMinPt_ ) continue;
      //opposite charge 
      if( (iMuon1->charge())*(iMuon2->charge()) == 1) continue;
      
      TrackRef glbTrackP;	  
      TrackRef glbTrackM;	  
      
      if(iMuon1->charge() == 1){ glbTrackP = iMuon1->track();}
      if(iMuon1->charge() == -1){ glbTrackM = iMuon1->track();}
      
      if(iMuon2->charge() == 1) { glbTrackP = iMuon2->track();}
      if(iMuon2->charge() == -1){ glbTrackM = iMuon2->track();}
      
      if( glbTrackP.isNull() || glbTrackM.isNull() ) 
	{
	  //std::cout << "continue due to no track ref" << endl;
	  continue;
	}
      
      
      if(!(glbTrackM->quality(reco::TrackBase::highPurity))) continue;
      if(!(glbTrackP->quality(reco::TrackBase::highPurity))) continue;	 
      
      reco::TransientTrack muon1TT((*theB).build(glbTrackP));
      reco::TransientTrack muon2TT((*theB).build(glbTrackM));
      
      // *****  Trajectory states to calculate DCA for the 2 muons *********************
      FreeTrajectoryState mu1State = muon1TT.impactPointTSCP().theState();
      FreeTrajectoryState mu2State = muon2TT.impactPointTSCP().theState();
      
      if( !muon1TT.impactPointTSCP().isValid() || !muon2TT.impactPointTSCP().isValid() ) continue;
      
      // Measure distance between tracks at their closest approach
      ClosestApproachInRPhi cApp;
      cApp.calculate(mu1State, mu2State);
      if( !cApp.status() ) continue;
      float dca = fabs( cApp.distance() );	
      mumuDCA_ .push_back(dca);  
      
      // *****  end DCA for the 2 muons *********************
      
      //Let's check the vertex and mass
      
      //The mass of a muon and the insignificant mass sigma 
      //to avoid singularities in the covariance matrix.
      ParticleMass muon_mass = 0.10565837; //pdg mass
      ParticleMass psi_mass = 3.096916;
      float muon_sigma = muon_mass*1.e-6;
      
      //Creating a KinematicParticleFactory
      KinematicParticleFactoryFromTransientTrack pFactory;
      
      //initial chi2 and ndf before kinematic fits.
      float chi = 0.;
      float ndf = 0.;
      vector<RefCountedKinematicParticle> muonParticles;
      try {
	muonParticles.push_back(pFactory.particle(muon1TT,muon_mass,chi,ndf,muon_sigma));
	muonParticles.push_back(pFactory.particle(muon2TT,muon_mass,chi,ndf,muon_sigma));
      }
      catch(...) { 
	std::cout<<" Exception caught ... continuing 1 "<<std::endl; 
	continue;
      }
      
      KinematicParticleVertexFitter fitter;   
      
      RefCountedKinematicTree psiVertexFitTree;
      try {
	psiVertexFitTree = fitter.fit(muonParticles); 
      }
      catch (...) { 
	std::cout<<" Exception caught ... continuing 2 "<<std::endl; 
	continue;
      }
      
      if (!psiVertexFitTree->isValid()) 
	{
	  //std::cout << "caught an exception in the psi vertex fit" << std::endl;
	  continue; 
	}
      
      psiVertexFitTree->movePointerToTheTop();
      
      RefCountedKinematicParticle psi_vFit_noMC = psiVertexFitTree->currentParticle();
      RefCountedKinematicVertex psi_vFit_vertex_noMC = psiVertexFitTree->currentDecayVertex();
      
      if( psi_vFit_vertex_noMC->chiSquared() < 0 )
	{
	  //std::cout << "negative chisq from psi fit" << endl;
	  continue;
	}
      
      //some loose cuts go here
      
      if(psi_vFit_vertex_noMC->chiSquared()>50.) continue;
      //if(psi_vFit_noMC->currentState().mass()<2.9 || psi_vFit_noMC->currentState().mass()>3.3) continue;
      
      double J_Prob_tmp   = TMath::Prob(psi_vFit_vertex_noMC->chiSquared(),(int)psi_vFit_vertex_noMC->degreesOfFreedom());
      mumuVtxCl_ .push_back(J_Prob_tmp);
      //if(J_Prob_tmp<0.01)
      //	{
      //	  continue;
      //	}
      
      

  // compute the distance between mumu vtx and beam spot
  calLS (psi_vFit_vertex_noMC->position().x(),psi_vFit_vertex_noMC->position().y(),0.0,
	 beamSpot_.position().x(),beamSpot_.position().y(),0.0,
	 psi_vFit_vertex_noMC->error().cxx(),psi_vFit_vertex_noMC->error().cyy(),0.0,
	 psi_vFit_vertex_noMC->error().matrix()(0,1),0.0,0.0,
	 beamSpot_.covariance()(0,0),beamSpot_.covariance()(1,1),0.0,
	 beamSpot_.covariance()(0,1),0.0,0.0,
	 &MuMuLSBS,&MuMuLSBSErr);

      mumuLSBS_.push_back(MuMuLSBS);      
      mumuLSBSErr_.push_back(MuMuLSBSErr);      
      

      reco::TrackRef muTrackm = iMuon1->innerTrack();
      if ( muTrackm.isNull() ) continue;
      

      
      muPt_    .push_back(iMuon1->pt());
      muEta_   .push_back(iMuon1->eta());
      muPhi_   .push_back(iMuon1->phi());
      muCharge_.push_back(iMuon1->charge());
      muType_  .push_back(iMuon1->type());
      
      histContainer_["muonPt"] ->Fill(iMuon1->pt());
      histContainer_["muonEta"]->Fill(iMuon1->eta());
      histContainer_["muonPhi"]->Fill(iMuon1->phi());
      
      
      
      reco::TrackRef muTrackp = iMuon2->innerTrack();
      if ( muTrackp.isNull()) continue;
      muPt_    .push_back(iMuon2->pt());
      muEta_   .push_back(iMuon2->eta());
      muPhi_   .push_back(iMuon2->phi());
      muCharge_.push_back(iMuon2->charge());
      muType_  .push_back(iMuon2->type());
      
      TLorentzVector mup, mum, dimu;
      mup.SetPtEtaPhiM(iMuon2->pt(),iMuon2->eta(),iMuon2->phi(),iMuon2->mass());
      mum.SetPtEtaPhiM(iMuon1->pt(),iMuon1->eta(),iMuon1->phi(),iMuon1->mass());
      
      dimu = mup+mum;
      mumuPt_.push_back(dimu.Pt());      
      mumuRapidity_.push_back(dimu.Rapidity());      
      mumuMass_.push_back(dimu.M());      

      histContainer_["muonPt"] ->Fill(iMuon2->pt());
      histContainer_["muonEta"]->Fill(iMuon2->eta());
      histContainer_["muonPhi"]->Fill(iMuon2->phi());
      
      histContainer_["mumuMass"]->Fill(dimu.M());
      histContainer_["mumupt"]->Fill(dimu.Pt());



    } // +mu
  } //-mu
  
  
  
  
  // get pat photon collection 
  edm::Handle< std::vector<pat::Photon>> photons;
  iEvent.getByToken(photonCollToken, photons);
  
  // loop over photons
  for (auto pho = photons->cbegin(); pho != photons->cend(); ++pho) {
    phoE_             .push_back(pho->energy());
    phoEt_            .push_back(pho->et());
    phoEta_           .push_back(pho->eta());
    phoPhi_           .push_back(pho->phi()); 
    nPho_++;
    
    histContainer_["phoPt"] ->Fill(pho->pt());
    histContainer_["phoEta"]->Fill(pho->eta());
    histContainer_["phoPhi"]->Fill(pho->phi());    
  }
  
  // Multiplicity
  histContainer_["phoMult" ]->Fill(photons->size());
  histContainer_["muonMult"]->Fill(muons->size() );
 
  //if (!hasBeamSpot(iEvent))continue;
  //if (!hasPrimaryVertex(iEvent)) continue;
   
  tree_->Fill();
  clearVariables();
}

void 
MiniAnalyzer::beginJob()
{
  // register to the TFileService
  edm::Service<TFileService> fs;
  
 


  histContainer_["mumuMass"]=fs->make<TH1F>("mumuMass", "mass",    90,   0, 120.);
  histContainer_["mumupt"]=fs->make<TH1F>( "mumupt", "#mu^{+}#mu^{-} pT ; pT [GeV]", 100, 0, 50);
  
  // book histograms for Multiplicity:
  
  histContainer_["phoMult"]=fs->make<TH1F>("phoMult",   "photon multiplicity", 100, 0,  50);
  histContainer_["muonMult"]=fs->make<TH1F>("muonMult",   "muon multiplicity",     100, 0,  50);
  
  // book histograms for Pt:
  
  histContainer_["phoPt"]=fs->make<TH1F>("phoPt",   "photon Pt", 100, 0,  200);
  histContainer_["muonPt"]=fs->make<TH1F>("muonPt",   "muon Pt", 100, 0, 200);

  // book histograms for Eta:
  histContainer_["phoEta"]=fs->make<TH1F>("phoEta",   "photon Eta",100, -5,  5);
  histContainer_["muonEta"]=fs->make<TH1F>("muonEta",   "muon Eta",  100, -5,  5);


  // book histograms for Phi:
  histContainer_["phoPhi"]=fs->make<TH1F>("phoPhi",   "photon Phi", 100, -3.5, 3.5);
  histContainer_["muonPhi"]=fs->make<TH1F>("muonPhi",   "muon Phi",     100, -3.5, 3.5);
    
}

bool
MiniAnalyzer::hasBeamSpot(const edm::Event& iEvent)
{
  edm::Handle<reco::BeamSpot> bs;
  iEvent.getByToken(bsCollToken,bs);

  if ( ! bs.isValid() ) {
    edm::LogError("myBeam") << "No beam spot available from EventSetup" ;
    return false;
  }

  beamSpot_ = *bs;
  return true;
}


bool
MiniAnalyzer::hasPrimaryVertex(const edm::Event& iEvent)
{

 // get primary vertex collection 
  edm::Handle< std::vector<reco::Vertex>> vtxs;
  iEvent.getByToken(vtxCollToken, vtxs);

  nprivtx = vtxs->size();

  for (std::vector<reco::Vertex>::const_iterator iVertex = vtxs->begin();
       iVertex != vtxs->end(); iVertex++) {
    primaryVertex_ = *(iVertex);
    if (primaryVertex_.isValid()) break;
  }

  if (!primaryVertex_.isValid()) return false;

  return true;
}


//Check recursively if any ancestor of particle is the given one
bool MiniAnalyzer::isAncestor(const reco::Candidate* ancestor, const reco::Candidate * particle) {
   if (ancestor == particle && ancestor->pdgId()==particle->pdgId()) return true;
   for (size_t i=0; i< particle->numberOfMothers(); i++) {
      if (isAncestor(ancestor, particle->mother(i))) return true;
   }
   return false;
}


void
MiniAnalyzer::calCosAlpha (double Vx, double Vy, double Vz,
			   double Wx, double Wy, double Wz,
			   double VxErr2, double VyErr2, double VzErr2,
			   double VxyCov, double VxzCov, double VyzCov,
			   double WxErr2, double WyErr2, double WzErr2,
			   double WxyCov, double WxzCov, double WyzCov,
			   double* cosAlpha, double* cosAlphaErr)
{
  double Vnorm = sqrt(Vx*Vx + Vy*Vy + Vz*Vz);
  double Wnorm = sqrt(Wx*Wx + Wy*Wy + Wz*Wz);
  double VdotW = Vx*Wx + Vy*Wy + Vz*Wz;

  if ((Vnorm > 0.) && (Wnorm > 0.)) {
    *cosAlpha = VdotW / (Vnorm * Wnorm);
    *cosAlphaErr = sqrt( (
			  (Vx*Wnorm - VdotW*Wx) * (Vx*Wnorm - VdotW*Wx) * WxErr2 +
			  (Vy*Wnorm - VdotW*Wy) * (Vy*Wnorm - VdotW*Wy) * WyErr2 +
			  (Vz*Wnorm - VdotW*Wz) * (Vz*Wnorm - VdotW*Wz) * WzErr2 +

			  (Vx*Wnorm - VdotW*Wx) * (Vy*Wnorm - VdotW*Wy) * 2.*WxyCov +
			  (Vx*Wnorm - VdotW*Wx) * (Vz*Wnorm - VdotW*Wz) * 2.*WxzCov +
			  (Vy*Wnorm - VdotW*Wy) * (Vz*Wnorm - VdotW*Wz) * 2.*WyzCov) /
			 (Wnorm*Wnorm*Wnorm*Wnorm) +
			 
			 ((Wx*Vnorm - VdotW*Vx) * (Wx*Vnorm - VdotW*Vx) * VxErr2 +
			  (Wy*Vnorm - VdotW*Vy) * (Wy*Vnorm - VdotW*Vy) * VyErr2 +
			  (Wz*Vnorm - VdotW*Vz) * (Wz*Vnorm - VdotW*Vz) * VzErr2 +
			  
			  (Wx*Vnorm - VdotW*Vx) * (Wy*Vnorm - VdotW*Vy) * 2.*VxyCov +
			  (Wx*Vnorm - VdotW*Vx) * (Wz*Vnorm - VdotW*Vz) * 2.*VxzCov +
			  (Wy*Vnorm - VdotW*Vy) * (Wz*Vnorm - VdotW*Vz) * 2.*VyzCov) /
			 (Vnorm*Vnorm*Vnorm*Vnorm) ) / (Wnorm*Vnorm);
  }  else {
    *cosAlpha = 0.;
    *cosAlphaErr = 0.;
  }

}


void
MiniAnalyzer::calLS (double Vx, double Vy, double Vz,
		     double Wx, double Wy, double Wz,
		     double VxErr2, double VyErr2, double VzErr2,
		     double VxyCov, double VxzCov, double VyzCov,
		     double WxErr2, double WyErr2, double WzErr2,
		     double WxyCov, double WxzCov, double WyzCov,
		     double* deltaD, double* deltaDErr)
{
  *deltaD = sqrt((Vx-Wx) * (Vx-Wx) + (Vy-Wy) * (Vy-Wy) + (Vz-Wz) * (Vz-Wz));
  if (*deltaD > 0.)
    *deltaDErr = sqrt((Vx-Wx) * (Vx-Wx) * VxErr2 +
		      (Vy-Wy) * (Vy-Wy) * VyErr2 +
		      (Vz-Wz) * (Vz-Wz) * VzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*VxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*VxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*VyzCov +
		      
		      (Vx-Wx) * (Vx-Wx) * WxErr2 +
		      (Vy-Wy) * (Vy-Wy) * WyErr2 +
		      (Vz-Wz) * (Vz-Wz) * WzErr2 +
		      
		      (Vx-Wx) * (Vy-Wy) * 2.*WxyCov +
		      (Vx-Wx) * (Vz-Wz) * 2.*WxzCov +
		      (Vy-Wy) * (Vz-Wz) * 2.*WyzCov) / *deltaD;
  else *deltaDErr = 0.;
}


void MiniAnalyzer::saveGenInfo(const edm::Event& iEvent){

  // Pruned particles are the one containing "important" stuff
  edm::Handle<edm::View<reco::GenParticle>> pruned;
  iEvent.getByToken(prunedCollToken,pruned);
  
  // Packed particles are all the status 1, so usable to remake jets
  // The navigation from status 1 to pruned is possible (the other direction should be made by hand)
  edm::Handle<edm::View<pat::PackedGenParticle>> packed;
  iEvent.getByToken(genCollToken, packed);
  
  
  if (packed.isValid() && pruned.isValid()) {
    dimuon_pdgId  = 0;
    gen_dimuon_p4.SetPtEtaPhiM(0.,0.,0.,0.);
    int foundit   = 0;
    int ngenPosMu=0;      int ngenNegMu=0;     int ngenPhoton=0; int ngenBs=0;

    for (size_t i=0; i<pruned->size(); i++) {
      int p_id = 0;
      p_id = abs((*pruned)[i].pdgId());
      const reco::Candidate *aonia = &(*pruned)[i];
      if ( p_id == 531 && aonia->status() == 2) {
	dimuon_pdgId = p_id;
        foundit++;

	int imum(-1), imup(-1), ipho(-1);
	//int found_daugter = 0;	
	for (size_t j=0; j<packed->size(); j++) { //get the pointer to the first survied ancestor of a given packed GenParticle in the prunedCollection
	  const reco::Candidate * motherInPrunedCollection = (*packed)[j].mother(0);
	  const reco::Candidate * d = &(*packed)[j];
	  
          if(aonia->pdgId() != motherInPrunedCollection->pdgId()) continue;
	  if ( motherInPrunedCollection != nullptr && (d->pdgId() == 13 ) && isAncestor(aonia , motherInPrunedCollection) && d->status()==1 ){
             imum = j;  foundit++;
            //std::cout << " coming in muon M loop:" << j << " onia pdgID:" << aonia->pdgId() << " mother ID:";
            //std::cout << motherInPrunedCollection->pdgId() << std::endl;
	    ngenMMu_++;
	    gen_muonM_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
	    genmumPt_ .push_back(d->pt());
	    genmumEta_ .push_back(d->eta());
	    genmumRap_ .push_back(d->rapidity());
	    genmumPhi_ .push_back(d->phi()); 
	  } 
	  if ( motherInPrunedCollection != nullptr && (d->pdgId() == -13 ) && isAncestor(aonia , motherInPrunedCollection) && d->status()==1 ){
             imup = j;    foundit++;
            //std::cout << " coming in muon ++ loop:" << j << " onia pdgID:" << aonia->pdgId() << " mother ID:";
            //std::cout << motherInPrunedCollection->pdgId() << std::endl;
	    ngenPMu_++;
	    gen_muonP_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
	    genmupPt_ .push_back(d->pt());
	    genmupEta_ .push_back(d->eta());
	    genmupRap_ .push_back(d->rapidity());
	    genmupPhi_ .push_back(d->phi());    
	  }
	  
	  if ( motherInPrunedCollection != nullptr && (d->pdgId() == 22 ) && isAncestor(aonia , motherInPrunedCollection) && d->status()==1 ){
             ipho = j;   foundit++;
            //std::cout << " coming in photon loop:" << j << " pt:"<< d->pt() << " eta:"<< d->eta() << " phi:"<< d->phi();
            //std::cout <<  " onia pdgID:" << aonia->pdgId() << " mother ID:" << motherInPrunedCollection->pdgId() << std::endl;
	    ngenPhoton_++;
	    gen_photon_p4.SetPtEtaPhiM(d->pt(),d->eta(),d->phi(),d->mass());
	    genphoPt_  .push_back(d->pt());
	    genphoEta_ .push_back(d->eta());
	    genphoRap_ .push_back(d->rapidity());
	    genphoPhi_ .push_back(d->phi());

	    genphoPx_ .push_back(d->px());
	    genphoPy_ .push_back(d->py());
	    genphoPz_ .push_back(d->pz());    
  	    genphoVtxx_.push_back(d->vx());    
   	    genphoVtxy_.push_back(d->vy());    
    	    genphoVtxz_.push_back(d->vz());    
	    
	  }
	} // packed MC size*/
 	
       if ( ipho == -1 || imum ==-1 || imup == -1) continue;
        ngenBs_++;
        //std::cout << " loop i:" << i << "  PDG ID:" << (*pruned)[i].pdgId() << " status:" << aonia->pdgId() ; 
        //std::cout << " bs px:" << aonia->px() << " py:" << aonia->py() << " pz:" << aonia->pz() << std::endl;
	genbPx_ .push_back(aonia->px());
	genbPy_ .push_back(aonia->py());
	genbPz_ .push_back(aonia->pz());

	genbPt_   .push_back(aonia->pt());
	genbEta_  .push_back(aonia->eta());
	genbRap_  .push_back(aonia->rapidity());
	genbPhi_  .push_back(aonia->phi());
	genbMass_ .push_back(aonia->mass());
	genbPID_  .push_back(aonia->pdgId());

     } // if Bs nd status 2     
    } // for pruned size
  }// if MC
}


void 
MiniAnalyzer::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAnalyzer);
