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
  
  int truth_nMuon, truth_nMuPos, truth_nMuNeg, truth_nPhoton, truth_nBs;
  std::vector<double> truth_Bsmuon_pt, truth_Bsmuon_eta, truth_Bsmuon_phi;
  std::vector<double> truth_Bsphoton_pt, truth_Bsphoton_eta, truth_Bsphoton_phi;
  std::vector<double> truth_Bs_pt, truth_Bs_eta, truth_Bs_phi, truth_Bs_mass;
  std::vector<double> truth_Bs_pdgid, truth_Bsmuon_pdgid;
  
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
  
  tree_->Branch("truth_nMuon"    ,  &truth_nMuon);
  tree_->Branch("truth_nMuPos"   ,  &truth_nMuPos);
  tree_->Branch("truth_nMuNeg"   ,  &truth_nMuNeg);
  tree_->Branch("truth_nPhoton"  ,  &truth_nPhoton);
  tree_->Branch("truth_nBs"      ,  &truth_nBs);
  tree_->Branch("truth_Bsmuon_pt",  &truth_Bsmuon_pt);
  tree_->Branch("truth_Bsmuon_eta", &truth_Bsmuon_eta);
  tree_->Branch("truth_Bsmuon_phi", &truth_Bsmuon_phi);
  tree_->Branch("truth_Bsmuon_pdgid", &truth_Bsmuon_pdgid);
  tree_->Branch("truth_Bsphoton_pt",  &truth_Bsphoton_pt);
  tree_->Branch("truth_Bsphoton_eta", &truth_Bsphoton_eta);
  tree_->Branch("truth_Bsphoton_phi", &truth_Bsphoton_phi);
  tree_->Branch("truth_Bs_pt",   &truth_Bs_pt);
  tree_->Branch("truth_Bs_eta",  &truth_Bs_eta);
  tree_->Branch("truth_Bs_phi",  &truth_Bs_phi);
  tree_->Branch("truth_Bs_mass", &truth_Bs_mass);
  tree_->Branch("truth_Bs_pdgid", &truth_Bs_pdgid);
  
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

  truth_nMuon=0; truth_nMuPos = 0; truth_nMuNeg = 0; truth_nPhoton = 0; truth_nBs = 0;
  truth_Bsmuon_pt.clear(); 
  truth_Bsmuon_eta.clear();
  truth_Bsmuon_phi.clear();
  truth_Bsmuon_pdgid.clear();
  truth_Bs_pt.clear(); truth_Bs_eta.clear(); truth_Bs_phi.clear(); truth_Bs_mass.clear();
  truth_Bsphoton_pt.clear(); 
  truth_Bsphoton_eta.clear();
  truth_Bsphoton_phi.clear();
  
  truth_Bs_pdgid.clear();  
  
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
      nMu_++;
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
  
  int BS = 531;
  int MUON = 13;
  int ANTIMUON = -13;
  int PHOTON = 22;
  
  bool found_mu = false;
  bool found_anti_mu = false;
  bool found_photon = false;
  
  for(unsigned int i = 0; i < pruned->size(); ++i) {
    const reco::GenParticle & gen_particle = (*pruned)[i];                                             
    if ( TMath::Abs(gen_particle.pdgId()) != BS ) continue; 
    int imum(-1), imup(-1), ipho(-1);
    
    // loop over all Bs daughters
    for (size_t j = 0; j < gen_particle.numberOfDaughters(); ++j) {
      const reco::Candidate  &dau = *(gen_particle.daughter(j));
      //std::cout << "Event:" << iEvent.id().event() << "  j:" << j << "   PID:" << dau.pdgId() << "  pt:" << dau.pt() << "  Eta:" << dau.eta() << std::endl;
      //if (dau.pdgId() == MUON  && dau.status()==1){
      if (dau.pdgId() == MUON  ){
	imup = j;
	truth_Bsmuon_pt .push_back(dau.pt());
	truth_Bsmuon_eta.push_back(dau.eta());
	truth_Bsmuon_phi.push_back(dau.phi());
	truth_Bsmuon_pdgid.push_back(dau.charge());
	truth_nMuon++;
	truth_nMuPos++;
	//std::cout << "Event:" << iEvent.id().event() << "  Run:"<< iEvent.id().run() << "  lumi:" << iEvent.id().lumis() << ++ muon daughter:" << dau.numberOfDaughters() << std::endl;
      }
      if (dau.pdgId() == ANTIMUON ){
	//if (dau.pdgId() == ANTIMUON && dau.status()==1){
	imum = j;
	truth_Bsmuon_pt .push_back(dau.pt());
	truth_Bsmuon_eta.push_back(dau.eta());
	truth_Bsmuon_phi.push_back(dau.phi());
	truth_Bsmuon_pdgid.push_back(dau.charge());	
	truth_nMuon++;
	truth_nMuNeg++;
	//std::cout << "Event:" << iEvent.id().event() << "  -- muon daughter:" << dau.numberOfDaughters() << std::endl;
      }
      if (dau.pdgId() == PHOTON   && dau.status()==1)  { ipho = j;
	truth_Bsphoton_pt .push_back(dau.pt());
	truth_Bsphoton_eta.push_back(dau.eta());
	truth_Bsphoton_phi.push_back(dau.phi());
	truth_nPhoton++;
	//std::cout << "Event:" << iEvent.id().event() << "  photon multiplicity:" << truth_nPhoton << "   pt:" << dau.pt() << "  Eta:" << dau.eta();
	//std::cout << "  mother:" << (dau.mother(0))->pdgId()  << std::endl;
      }
    }
    //std::cout << " Event:" << iEvent.id().event() << "   photon multiplicity:" << truth_nPhoton << std::endl << std::endl;
    if ( ipho == -1 ) continue;
    
    const reco::Candidate & pho = *(gen_particle.daughter(ipho));
    //std::cout << "event:" << ievent.id().event() << "  photon multiplicity:" << ipho << "   pt:" << pho.pt() << "  eta:" << pho.eta() << std::endl;
    const reco::Candidate *mum = NULL;
    const reco::Candidate *mup = NULL;
    
    if (imum != -1 && imup != -1) {
      mum = gen_particle.daughter(imum);
      mup = gen_particle.daughter(imup);
    }
    
    if ( mum == NULL || mup == NULL) continue;
    
    // std::cout << "event:" << iEvent.id().event() << "  photon multiplicity:" << ipho << "   pt:" << pho.pt() << "  eta:" << pho.eta() <<  "  muon number:" << imum << "   negative:" << imup << std::endl << std::endl;
    //std::cout << "positive pt:" << mup->pt() << std::endl;
    truth_Bs_pt  .push_back(gen_particle.pt());
    truth_Bs_eta .push_back(gen_particle.eta());
    truth_Bs_phi .push_back(gen_particle.phi());
    truth_Bs_mass.push_back(gen_particle.mass());
    truth_nBs++;

  } // mc particle size
}


void 
MiniAnalyzer::endJob() 
{
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MiniAnalyzer);
