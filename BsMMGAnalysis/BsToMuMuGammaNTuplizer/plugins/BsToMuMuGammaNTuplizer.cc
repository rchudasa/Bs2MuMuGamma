// -*- C++ -*-
//

// Package:    BsMMGAnalysis/BsToMuMuGammaNTuplizer
// Class:      BsToMuMuGammaNTuplizer
//
/**\class BsToMuMuGammaNTuplizer BsToMuMuGammaNTuplizer.cc BsMMGAnalysis/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.cc

   Description: Takes in AOD and makes NTuples of Muon/DiMuons/Photons

   Implementation:
   [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Wed, 12 May 2021 18:51:51 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"



// user include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "CondFormats/DataRecord/interface/L1TUtmTriggerMenuRcd.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "MagneticField/Engine/interface/MagneticField.h"

#include "TMath.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/TreeContent.h"
#include "BsMMGAnalysis/BsToMuMuGammaNTuplizer/interface/Utils.h"




//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.


using reco::TrackCollection;

class BsToMuMuGammaNTuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
public:
  explicit BsToMuMuGammaNTuplizer(const edm::ParameterSet&);
  ~BsToMuMuGammaNTuplizer();
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  bool printMsg;
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool isMC;
  bool doBsToMuMuGamma;
  TTree* theTree;
  TreeContent* NTuple;
  Utils* Utility;
 
 
  //l1t::L1TGlobalUtil *fGtUtil;

  // selection cuts;
  double pTMinMuons;
  double cl_dimuon_vtx;
  double ls_max_dimuonBS          ;
  double dcaMax_dimuon_mumu           ;
  double dcaMax_muon_bs           ;
  double cosAlphaMax_dimuonBs    ;
  double MUMINPT           ;
  double etaMax_muon          ;
  double minDimuon_pt         ;
  double minDimuonInvariantMass    ;
  double maxDimuonInvariantMass    ;
  double trackIP_zMax_muon 		;
  double trackIP_rMax_muon	        ;
    
  // variables
  double chi2,ndof;
  float muonMass,muonMassErr;


      
  //edm::EDGetTokenT<edm::TriggerResults>                    triggerBits_;
  //edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjects_;
  //edm::EDGetTokenT<pat::PackedTriggerPrescales>            triggerPrescales_;
  //std::vector<std::string> trigTable_;
  //std::vector<std::string> l1Table_;
      
  edm::EDGetTokenT<reco::BeamSpot>                  beamSpotToken_;
  edm::EDGetTokenT<reco::VertexCollection>          primaryVtxToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>>         muonToken_;
  edm::EDGetTokenT<reco::GenParticleCollection>     simGenTocken_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCBarrelCollection_;
  edm::EDGetTokenT<std::vector<reco::SuperCluster>>   MustacheSCEndcapCollection_;



  // input tags 
  edm::InputTag MustacheSCBarrelSrc_;
  edm::InputTag MustacheSCEndcapSrc_;


  // Dimuon Reco vars
  TrajectoryStateClosestToPoint theDCAXBS;
  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;


  RefCountedKinematicTree mumuVertexFitTree;
	 
  TLorentzVector bsDimuon_lv;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
BsToMuMuGammaNTuplizer::BsToMuMuGammaNTuplizer(const edm::ParameterSet& iConfig) :
  muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons"))),
  MustacheSCBarrelSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCBarrelSrc")),
  MustacheSCEndcapSrc_(iConfig.getParameter<edm::InputTag>("MustacheSCEndcapSrc"))
{
  //now do what ever initialization is needed
 
  beamSpotToken_  	   =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
  primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
  simGenTocken_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
  //trackToken_          =consumes<reco::TrackCollection>             (iConfig.getParameter<edm::InputTag>("tracks"));
  //triggerBits_         =consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("bits"));
  //triggerPrescales_    =consumes<pat::PackedTriggerPrescales>            (iConfig.getParameter<edm::InputTag>("prescales"));
  //trigTable_           =iConfig.getParameter<std::vector<std::string> >("TriggerNames");   
  //triggerObjects_   (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("objects"));
    
  //l1Table_           = iConfig.getParameter<std::vector<std::string> >("L1Names");   
  //mcGenToken_        = consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>(""));

  MustacheSCBarrelCollection_             = consumes<std::vector<reco::SuperCluster>>(MustacheSCBarrelSrc_);
  MustacheSCEndcapCollection_             = consumes<std::vector<reco::SuperCluster>>(MustacheSCEndcapSrc_);

  etaMax_muon               =  iConfig.getUntrackedParameter<double>("muon_EtaMax")        ;
  dcaMax_muon_bs            =  iConfig.getUntrackedParameter<double>("muon_dcaMAX")        ;
  pTMinMuons 		      =  iConfig.getUntrackedParameter<double>("muon_minPt");
  trackIP_zMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_zIPMax")        ;
  trackIP_rMax_muon	      =  iConfig.getUntrackedParameter<double>("muon_rIPMax")        ;

  minDimuon_pt              =  iConfig.getUntrackedParameter<double>("dimuon_minPt")      ;
  cl_dimuon_vtx             =  iConfig.getUntrackedParameter<double>("dimuon_minVtxCL")      ;
  ls_max_dimuonBS           =  iConfig.getUntrackedParameter<double>("dimuon_maxLStoBS")       ;
  dcaMax_dimuon_mumu        =  iConfig.getUntrackedParameter<double>("dimuon_maxDCAMuMu")        ;
  cosAlphaMax_dimuonBs      =  iConfig.getUntrackedParameter<double>("dimuon_maxCosAlphaToBS") ;
  minDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_minInvMass")    ;
  maxDimuonInvariantMass    =  iConfig.getUntrackedParameter<double>("dimuon_maxInvMass")    ;
    
  printMsg=iConfig.getParameter<bool>("verbose");
  isMC=iConfig.getParameter<bool>("isMC");
  doBsToMuMuGamma=iConfig.getParameter<bool>("doBsToMuMuGamma");
    
  NTuple = new TreeContent;
  Utility= new Utils();
  muonMass= Utility->muonMass;
  muonMassErr= Utility->muonMassErr;

}


BsToMuMuGammaNTuplizer::~BsToMuMuGammaNTuplizer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
  //delete NTuple;
  //delete Utility;


}


//
// member functions
//

// ------------ method called for each event  ------------
void
BsToMuMuGammaNTuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
 
  // Get magnetic field
    
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
  //  Get BeamSpot
  edm::Handle<reco::BeamSpot> beamSpotH;
  iEvent.getByToken(beamSpotToken_, beamSpotH);
  reco::BeamSpot beamSpot = *beamSpotH;

  edm::Handle<std::vector<reco::Muon>> muons;
  iEvent.getByToken(muonToken_, muons);
 
  edm::Handle<reco::GenParticleCollection> genParticleCollection;
  if(isMC) iEvent.getByToken(simGenTocken_, genParticleCollection);
 
  edm::Handle<std::vector<reco::Vertex>> primaryVertexCollection;
  iEvent.getByToken(primaryVtxToken_, primaryVertexCollection);
 
  // adding  BEAMSOPT 
  NTuple->beamspot_x			= beamSpot.x0();  ;
  NTuple->beamspot_y			= beamSpot.y0();  ;
  NTuple->beamspot_z			= beamSpot.z0();  ;
  NTuple->beamspot_x_error		= beamSpot.x0Error();  ;
  NTuple->beamspot_y_error		= beamSpot.y0Error();  ;
  NTuple->beamspot_z_error		= beamSpot.z0Error();  ;
  NTuple->beamspot_dxdz   		= beamSpot.dxdz();  ;
  NTuple->beamspot_dydz	         	= beamSpot.dydz();  ;
  NTuple->beamspot_sigmaZ		= beamSpot.sigmaZ();  ;
  NTuple->beamspot_dxdz_error		= beamSpot.dxdzError();  ;
  NTuple->beamspot_dydz_error		= beamSpot.dydzError();  ;
  NTuple->beamspot_sigmaZError		= beamSpot.sigmaZ0Error();  ;
  NTuple->beamspot_beamWidthX		= beamSpot.BeamWidthX();  ;
  NTuple->beamspot_beamWidthY		= beamSpot.BeamWidthY();  ;
  NTuple->beamspot_beamWidthX_error	= beamSpot.BeamWidthXError();  ;
  NTuple->beamspot_beamWidthY_error	= beamSpot.BeamWidthXError();  ;
 
  for(auto&  aVertex : *primaryVertexCollection){

    if( not aVertex.isValid() ) continue;
    
    // # offlinePrimaryVertices # //
    (NTuple->primaryVertex_isFake ).push_back(   aVertex.isFake() );
    (NTuple->primaryVertex_x ).push_back(   aVertex.x() );
    (NTuple->primaryVertex_y ).push_back(   aVertex.y()  );
    (NTuple->primaryVertex_z ).push_back(   aVertex.z()  );
    (NTuple->primaryVertex_t ).push_back(   aVertex.t()  );
    (NTuple->primaryVertex_x_error ).push_back(   aVertex.xError()  );
    (NTuple->primaryVertex_y_error ).push_back(   aVertex.yError() );
    (NTuple->primaryVertex_z_error ).push_back(   aVertex.zError()  );
    (NTuple->primaryVertex_t_error ).push_back(   aVertex.tError()  );
    (NTuple->primaryVertex_ntracks ).push_back(   aVertex.nTracks() );
    (NTuple->primaryVertex_ndof ).push_back(   aVertex.ndof() 	 	  );
    (NTuple->primaryVertex_chi2 ).push_back(   aVertex.chi2()  );
    (NTuple->primaryVertex_normalizedChi2 ).push_back(   aVertex.normalizedChi2()  );
  } // loop over primary vertex collection
  
  int phoMul(0),muMMul(0),muPMul(0);
  
  if(isMC){
    
    for(auto& aBsMeson : *genParticleCollection){
      
      if(abs(aBsMeson.pdgId())!=531) continue;
      
      for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
	
	auto& bsDaughter = *(aBsMeson.daughter(j));
	
	if(bsDaughter.pdgId() == -13) muMMul++;
	if(bsDaughter.pdgId() ==  13) muPMul++;
	if(bsDaughter.pdgId() ==  22) phoMul++;
      }
      if(muMMul!=1 or muPMul!=1 or phoMul!=1 ){
	
	muMMul=0;
	muPMul=0;
	phoMul=0;
	continue;
      }
      
      for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++){
	
	auto& bsDaughter = *(aBsMeson.daughter(j));
	if(bsDaughter.pdgId() == -13) {
	  (NTuple->gen_BsMuonM_pt).push_back(bsDaughter.pt());
	  (NTuple->gen_BsMuonM_eta).push_back(bsDaughter.pt());
	  (NTuple->gen_BsMuonM_phi).push_back(bsDaughter.pt());
	}
	
	if(bsDaughter.pdgId() ==  13){
	  (NTuple->gen_BsMuonP_pt).push_back(bsDaughter.pt());
	  (NTuple->gen_BsMuonP_eta).push_back(bsDaughter.pt());
	  (NTuple->gen_BsMuonP_phi).push_back(bsDaughter.pt());
	}

	if(bsDaughter.pdgId() ==  22){
	  (NTuple->gen_BsPhoton_pt).push_back(bsDaughter.pt());
	  (NTuple->gen_BsPhoton_eta).push_back(bsDaughter.pt());
	  (NTuple->gen_BsPhoton_phi).push_back(bsDaughter.pt());
	  
	}
	
      } // loop over number of daughters
      
      (NTuple->gen_Bs_pt).push_back(aBsMeson.pt());
      (NTuple->gen_Bs_eta).push_back(aBsMeson.eta());
      (NTuple->gen_Bs_phi).push_back(aBsMeson.phi());
      (NTuple->gen_Bs_pz).push_back(aBsMeson.pz());
      (NTuple->gen_Bs_pdgId).push_back(aBsMeson.pdgId());
      
      break;
      
    }
    
    bool hasAValidMCCandidate= true;
    if(muMMul!=1 or muPMul!=1 ) {
      
      if(printMsg) std::cout<<" Ghost event found !! Mu+ Mu- from any of the Bs not found to == 1 "<<std::endl;
      hasAValidMCCandidate = false;
    }
    else if(doBsToMuMuGamma and phoMul!=1){
      if(printMsg) std::cout<<" Ghost event found !! gamma multiplicity from any of the Bs not found to == 1 "<<std::endl;
      hasAValidMCCandidate = false;
    }
    
    (NTuple->gen_hasAValid_candidate).push_back(hasAValidMCCandidate);
  } // If is MC

  //  Muon Ntuplizing
  //  TODO : Add details to closest PV
  //         Add details with BS
  //         Add distance of dimuon vertex to BS and its err
  //
  
  

  reco::TrackRef muTrackm,muTrackp;
  reco::TransientTrack refitMupTT, refitMumTT;
  double mu_mu_vtx_cl, mu_mu_pt, mu_mu_mass, mu_mu_mass_err;
  int num_SC = 0;
 
  //start loop for positive muons 
  for (uint32_t i=0;i<muons->size();i++){
    
    auto &mum=muons->at(i);
    if(mum.pt()  < pTMinMuons) continue;
    muTrackm= mum.innerTrack();
    if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1))   continue;
    
    const reco::TransientTrack muTrackmTT( muTrackm, &(*bFieldHandle));
    if (!muTrackmTT.isValid()) continue;
    
    // # Compute mu- DCA to BeamSpot #
    theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
    double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
    double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
    
    
    //start loop for positive muons 
    for ( uint32_t j=0;j<muons->size();j++){
      auto &mup=muons->at(j);
      muTrackp = mup.innerTrack();
      if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1)) continue;
      
      const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));
      if (!muTrackpTT.isValid()) continue;
      
      // # Compute mu+ DCA to BeamSpot #
      theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
      double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
      double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
      
      
      // # Check goodness of muons closest approach #
      ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
      XingPoint = ClosestApp.crossingPoint();
      
      double mumuDCA = ClosestApp.distance();
      
      // # dimuon inviariant mass and pT before vertex fitting #
      bsDimuon_lv.SetPxPyPzE( muTrackmTT.track().px() + muTrackpTT.track().px(), 
			      muTrackmTT.track().py() + muTrackpTT.track().py(),
			      muTrackmTT.track().pz() + muTrackpTT.track().pz(),
			      sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ) + sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) ) );
      
      
      if (true or (bsDimuon_lv.Pt() < minDimuon_pt)  || (bsDimuon_lv.M() < minDimuonInvariantMass) || (bsDimuon_lv.M() > maxDimuonInvariantMass))
	{
	  if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << bsDimuon_lv.Pt() << "\tinv. mass: " << bsDimuon_lv.M() << std::endl;
	  //        continue;
	}
      
      
      chi2 = 0.;
      ndof  = 0.;
      
      // ####################################################
      // # Try to vertex the two muons to get dimuon vertex #
      // ####################################################
      KinematicParticleFactoryFromTransientTrack partFactory;
      KinematicParticleVertexFitter PartVtxFitter;
      
      std::vector<RefCountedKinematicParticle> muonParticles;
      muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi2,ndof,muonMassErr));
      muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass,chi2,ndof,muonMassErr));
      
      RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles); 
      if ( !mumuVertexFitTree->isValid()) continue;
      
      if (mumuVertexFitTree->isValid() == false){
	if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
	continue; 
      }
      
      mumuVertexFitTree->movePointerToTheTop();
      RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
      RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
      
      if ( !mumu_KV->vertexIsValid()) continue;
      
      mu_mu_vtx_cl = TMath::Prob((double)mumu_KV->chiSquared(),
				 int(rint(mumu_KV->degreesOfFreedom())));
      
      
      // extract the re-fitted tracks
      mumuVertexFitTree->movePointerToTheTop();
      
      mumuVertexFitTree->movePointerToTheFirstChild();
      RefCountedKinematicParticle refitMum = mumuVertexFitTree->currentParticle();
      refitMumTT = refitMum->refittedTransientTrack();
      
      mumuVertexFitTree->movePointerToTheNextChild();
      RefCountedKinematicParticle refitMup = mumuVertexFitTree->currentParticle();
      refitMupTT = refitMup->refittedTransientTrack();
      
      TLorentzVector mymum, mymup, mydimu;
      
      mymum.SetXYZM(refitMumTT.track().momentum().x(),
		    refitMumTT.track().momentum().y(),
		    refitMumTT.track().momentum().z(), muonMass);
      
      mymup.SetXYZM(refitMupTT.track().momentum().x(),
		    refitMupTT.track().momentum().y(),
		    refitMupTT.track().momentum().z(), muonMass);
      
      mydimu = mymum + mymup;
      mu_mu_pt = mydimu.Perp();
      
      mu_mu_mass = mumu_KP->currentState().mass();
      mu_mu_mass_err = sqrt(mumu_KP->currentState().kinematicParametersError().
			    matrix()(6,6));
      
      
      
      // ######################################################
      // # Compute the distance between mumu vtx and BeamSpot #
      // ######################################################
      
      double MuMuLSBS;
      double MuMuLSBSErr;
      Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
			  beamSpot.position().x(),beamSpot.position().y(),0.0,
			  mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
			  mumu_KV->error().matrix()(0,1),0.0,0.0,
			  beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
			  beamSpot.covariance()(0,1),0.0,0.0,
			  &MuMuLSBS,&MuMuLSBSErr);
      
      
      
      // ###################################################################
      // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
      // ###################################################################
      double MuMuCosAlphaBS;
      double MuMuCosAlphaBSErr;
      Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),
                                mumu_KP->currentState().globalMomentum().y(),
				0.0,
				mumu_KV->position().x() - beamSpot.position().x(),
				mumu_KV->position().y() - beamSpot.position().y(),
				0.0,
				mumu_KP->currentState().kinematicParametersError().matrix()(3,3),
				mumu_KP->currentState().kinematicParametersError().matrix()(4,4),
				0.0,
				mumu_KP->currentState().kinematicParametersError().matrix()(3,4),
				0.0,
				0.0,
				mumu_KV->error().cxx() + beamSpot.covariance()(0,0),
				mumu_KV->error().cyy() + beamSpot.covariance()(1,1),
				0.0,
				mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),
				0.0,
				0.0,
				&MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
    
      
      // add loop for photons here...

       
    
      //// #################
      //// # Save: mu+ mu- #
      //// #################
      (NTuple->mumuPt).push_back(mu_mu_pt);
      (NTuple->mumuMass).push_back(mu_mu_mass);
      (NTuple->mumuMassE).push_back(mu_mu_mass_err);
      
      (NTuple->mumuPx).push_back(mumu_KP->currentState().globalMomentum().x());
      (NTuple->mumuPy).push_back(mumu_KP->currentState().globalMomentum().y());
      (NTuple->mumuPz).push_back(mumu_KP->currentState().globalMomentum().z());
      
      (NTuple->mumuVtxCL).push_back(mu_mu_vtx_cl);
      (NTuple->mumuVtxX).push_back(mumu_KV->position().x());
      (NTuple->mumuVtxY).push_back(mumu_KV->position().y());
      (NTuple->mumuVtxZ).push_back(mumu_KV->position().z());
      
      (NTuple->mumuCosAlphaBS).push_back(MuMuCosAlphaBS);
      (NTuple->mumuCosAlphaBSE).push_back(MuMuCosAlphaBSErr);
      (NTuple->mumuLBS).push_back(MuMuLSBS);
      (NTuple->mumuLBSE).push_back(MuMuLSBSErr);
      (NTuple->mumuDCA).push_back(mumuDCA);
      
      
      //// #############
      //// # Save: mu- #
      //// #############
      (NTuple->mumHighPurity).push_back( (int)muTrackm->quality(reco::Track::highPurity));
      (NTuple->mumPt).push_back(muTrackm->pt());
      (NTuple->mumEta).push_back(muTrackm->eta());
      (NTuple->mumPhi).push_back(muTrackm->phi());
      (NTuple->mumCL).push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
      (NTuple->mumNormChi2).push_back(muTrackm->normalizedChi2());
      (NTuple->mumPx).push_back(refitMumTT.track().momentum().x());
      (NTuple->mumPy).push_back(refitMumTT.track().momentum().y());
      (NTuple->mumPz).push_back(refitMumTT.track().momentum().z());
      
      //                     (NTuple->mumDCAVtx).push_back(DCAmumVtx);
      //                     (NTuple->mumDCAVtxE).push_back(DCAmumVtxErr);
      (NTuple->mumDCABS).push_back(DCAmumBS);
      (NTuple->mumDCABSE).push_back(DCAmumBSErr);
      
      //                     (NTuple->mumKinkChi2).push_back(iMuonM->combinedQuality().trkKink);
      (NTuple->mumFracHits).push_back(static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
															   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
      //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
      //                     (NTuple->mumdxyVtx).push_back(theDCAXVtx.second.value());
      //                     (NTuple->mumdzVtx).push_back(muTrackmTT.track().dz( ));
      (NTuple->mumdxyBS).push_back(muTrackmTT.track().dxy( (beamSpot.position() )));
      (NTuple->mumdzBS ).push_back(muTrackmTT.track().dz(  (beamSpot.position() )));
      
      //(NTuple->mumCat).push_back(getMuCat(mum));
      
      (NTuple->mumNPixHits).push_back(muTrackm->hitPattern().numberOfValidPixelHits());
      (NTuple->mumNPixLayers).push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());  
      (NTuple->mumNTrkHits).push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
      (NTuple->mumNTrkLayers).push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
      if (mum.isGlobalMuon() == true) (NTuple->mumNMuonHits).push_back(mum.globalTrack()->hitPattern().numberOfValidMuonHits());
      else (NTuple->mumNMuonHits).push_back(0);
      (NTuple->mumNMatchStation).push_back(mum.numberOfMatchedStations());
      
      
      //// #############
      //// # Save: mu+ #
      //// #############
      (NTuple->mupHighPurity).push_back( (int) muTrackp->quality(reco::Track::highPurity));
      (NTuple->mupPt).push_back(muTrackp->pt());
      (NTuple->mupEta).push_back(muTrackp->eta());
      (NTuple->mupPhi).push_back(muTrackp->phi());
      //(NTuple->mupCL).push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
      (NTuple->mupNormChi2).push_back(muTrackp->normalizedChi2());
      (NTuple->mupPx).push_back(refitMupTT.track().momentum().x());
      (NTuple->mupPy).push_back(refitMupTT.track().momentum().y());
      (NTuple->mupPz).push_back(refitMupTT.track().momentum().z());
      
      (NTuple->mupDCABS).push_back(DCAmupBS);
      (NTuple->mupDCABSE).push_back(DCAmupBSErr);
      (NTuple->mupdxyBS).push_back(muTrackpTT.track().dxy( (beamSpot.position() )));
      (NTuple->mupdzBS ).push_back(muTrackpTT.track().dz(  (beamSpot.position() )));
      
      //                     (NTuple->mupKinkChi2).push_back(iMuonP->combinedQuality().trkKink);
      (NTuple->mupFracHits).push_back(static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackp->hitPattern().numberOfValidHits() +
															   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
															   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
															   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
      //                     theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
      //                     (NTuple->mupdxyVtx).push_back(theDCAXVtx.second.value());
      //                     (NTuple->mupdzVtx).push_back(muTrackpTT.track().dz(bestVtxReFit.position()));
      
      //(NTuple->mupCat).push_back(getMuCat(mup));
      
      (NTuple->mupNPixHits).push_back(muTrackp->hitPattern().numberOfValidPixelHits());
      (NTuple->mupNPixLayers).push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());  
      (NTuple->mupNTrkHits).push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
      (NTuple->mupNTrkLayers).push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
      if (mup.isGlobalMuon() == true) (NTuple->mupNMuonHits).push_back(mup.globalTrack()->hitPattern().numberOfValidMuonHits());
      else (NTuple->mupNMuonHits).push_back(0);
      (NTuple->mupNMatchStation).push_back(mup.numberOfMatchedStations());
    
      // loop for photons from Musctache SC starts here
      edm::Handle<reco::SuperClusterCollection> barrelSCHandle;
      iEvent.getByToken(MustacheSCBarrelCollection_, barrelSCHandle);
      
      edm::Handle<reco::SuperClusterCollection> endcapSCHandle;
      iEvent.getByToken(MustacheSCEndcapCollection_, endcapSCHandle);
      
      for (auto const& scs : { *barrelSCHandle, *endcapSCHandle }) {
        for (uint32_t kk = 0; kk <scs.size(); kk++){
 	    auto &sc = scs.at(kk);
 //	for (reco::SuperClusterCollection::const_iterator sc = scs.begin(); sc != scs.end(); ++sc) {
	//for (auto const& sc : scs) {
	  (NTuple->scE)           .push_back(sc.energy());
	  /*(NTuple->scEt)          .push_back(sc.energy()/cosh(sc.eta()));
	  (NTuple->scEta)         .push_back(sc.eta());
	  (NTuple->scPhi)         .push_back(sc.phi());
	  (NTuple->scX)           .push_back(sc.x());
	  (NTuple->scY)           .push_back(sc.y());
	  (NTuple->scZ)           .push_back(sc.z());
	  (NTuple->scEtaWidth)    .push_back(sc.etaWidth());
	  (NTuple->scPhiWidth)    .push_back(sc.phiWidth());        
	  (NTuple->scRawE)        .push_back(sc.rawEnergy());         
	  (NTuple->scRawEt)       .push_back(sc.rawEnergy()/cosh(sc.eta()));    
	  */++num_SC;

          reco::CompositeCandidate mmg;
          mmg.addDaughter(muons->at(i), "negative muon");
          mmg.addDaughter(muons->at(j), "pos muon");
          mmg.addDaughter(scs.at(kk), "photon");
          //mmg.addDaughter(*sc, "photon2");
  	  AddFourMomenta addP4;
          addP4.set(mmg);
	}
      } // loop over barrel and endcap SC
     NTuple->nSC = num_SC;
      
      muonParticles.clear(); 
    } // muon positive
  } //muon negative
  
  theTree->Fill();
  NTuple->ClearNTuple();

}


// ------------ method called once each job just before starting event loop  ------------
void
BsToMuMuGammaNTuplizer::beginJob()
{
  edm::Service<TFileService> outfile_;
  theTree = outfile_->make<TTree>("Event","Event Tree");
  NTuple->MakeTreeBranches(theTree);
}

// ------------ method called once each job just after ending the event loop  ------------
void
BsToMuMuGammaNTuplizer::endJob()
{

  //theTree->GetDirectory()->cd();
  //theTree->Write();

}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BsToMuMuGammaNTuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BsToMuMuGammaNTuplizer);
