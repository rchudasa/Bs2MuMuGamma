// -*- C++ -*-
//

// Package:    Bs2MuMuGamma/BsToMuMuGammaNTuplizer
// Class:      BsToMuMuGammaNTuplizer
//
/**\class BsToMuMuGammaNTuplizer BsToMuMuGammaNTuplizer.cc Bs2MuMuGamma/BsToMuMuGammaNTuplizer/plugins/BsToMuMuGammaNTuplizer.cc

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



#include "Bs2MuMuGamma/BsToMuMuGammaNtuplizer/interface/TreeContent.h"
#include "Bs2MuMuGamma/BsToMuMuGammaNtuplizer/interface/Utils.h"




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
	double cl_dimuon_vtx         ;
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
      edm::EDGetTokenT<reco::GenParticleCollection>   simGenTocken_;



     // Dimuon Reco vars
         TrajectoryStateClosestToPoint theDCAXBS;
	 ClosestApproachInRPhi ClosestApp;
	 GlobalPoint XingPoint;
         KinematicParticleVertexFitter PartVtxFitter;
	 KinematicParticleFactoryFromTransientTrack partFactory;
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
BsToMuMuGammaNTuplizer::BsToMuMuGammaNTuplizer(const edm::ParameterSet& iConfig) 
: muonToken_(consumes<reco::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
{
   //now do what ever initialization is needed
 
    beamSpotToken_  	   =consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"));
    primaryVtxToken_       =consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"));
    simGenTocken_          =consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"));
    //triggerBits_           =consumes<edm::TriggerResults>                    (iConfig.getParameter<edm::InputTag>("HLTBits"));
    //triggerPrescales_      =consumes<pat::PackedTriggerPrescales>            (iConfig.getParameter<edm::InputTag>("prescales"));
    //trigTable_             =iConfig.getParameter<std::vector<std::string> >("TriggerNames");   
    //triggerObjects_   (consumes<pat::TriggerObjectStandAloneCollection> (iConfig.getParameter<edm::InputTag>("objects"));
    
    //l1Table_           = iConfig.getParameter<std::vector<std::string> >("L1Names");   
    //mcGenToken_        = consumes<reco::GenParticleCollection >(iConfig.getParameter<edm::InputTag>(""));

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
 /*
    edm::Handle<edm::TriggerResults>                    triggerBits;
    //edm::Handle<pat::TriggerObjectStandAloneCollection> triggerObjects;
    //edm::Handle<pat::PackedTriggerPrescales>            triggerPrescales;

    iEvent.getByToken(triggerBits_,      triggerBits);
    //iEvent.getByToken(triggerObjects_,   triggerObjects);
    //iEvent.getByToken(triggerPrescales_, triggerPrescales);

    bool foundOneTrig = false;
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) {
        for (unsigned int it = 0; it < trigTable_.size(); it++){
            if (names.triggerName(i).find(trigTable_[it]) != std::string::npos)
            {
                NTuple->TrigTable->push_back(names.triggerName(i) );
                NTuple->TrigResult->push_back(triggerBits->accept(i));
		//NTuple->TrigPrescales->push_back(triggerPrescales->getPrescaleForIndex(i));
                foundOneTrig = true;
            }
        }
    }
    if ( iEvent.isRealData() && !foundOneTrig) return;
    if (NTuple->TrigTable->size() == 0)
    {
      NTuple->TrigTable->push_back("FINAL");
      NTuple->TrigResult->push_back(false);
    }
    else
    {
      NTuple->TrigTable->push_back( "FINAL" );
      NTuple->TrigResult->push_back( true );
    }
    
   */ 
    
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
   NTuple->beamspot_dxdz		= beamSpot.dxdz();  ;
   NTuple->beamspot_dydz		= beamSpot.dydz();  ;
   NTuple->beamspot_sigmaZ		= beamSpot.sigmaZ();  ;
   NTuple->beamspot_dxdz_error		= beamSpot.dxdzError();  ;
   NTuple->beamspot_dydz_error		= beamSpot.dydzError();  ;
   NTuple->beamspot_sigmaZError		= beamSpot.sigmaZ0Error();  ;
   NTuple->beamspot_beamWidthX		= beamSpot.BeamWidthX();  ;
   NTuple->beamspot_beamWidthY		= beamSpot.BeamWidthY();  ;
   NTuple->beamspot_beamWidthX_error	= beamSpot.BeamWidthXError();  ;
   NTuple->beamspot_beamWidthY_error	= beamSpot.BeamWidthXError();  ;
 
   for(auto&  aVertex : *primaryVertexCollection)
   {

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
	

   }
   
   int phoMul(0),muMMul(0),muPMul(0);
   
   if(isMC)
   {
	   for(auto& aBsMeson : *genParticleCollection)
	   {
		if(abs(aBsMeson.pdgId())!=531) continue;
		
		for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++)
	      	{
		    auto& bsDaughter = *(aBsMeson.daughter(j));
	 	    
		    if(bsDaughter.pdgId() == -13) muMMul++;
		    if(bsDaughter.pdgId() ==  13) muPMul++;
		    if(bsDaughter.pdgId() ==  22) phoMul++;
		}
		if(muMMul!=1 or muPMul!=1)
		{
			muMMul=0;
			muPMul=0;
			phoMul=0;
			continue;
		}
	
	       if(doBsToMuMuGamma and phoMul!=1)
		{
			muMMul=0;
			muPMul=0;
			phoMul=0;
			continue;
		}
	
		for(unsigned int j=0; j<aBsMeson.numberOfDaughters(); j++)
	      	{
	            auto& bsDaughter = *(aBsMeson.daughter(j));
		    if(bsDaughter.pdgId() == -13) 
		    {
			(NTuple->gen_BsMuonM_pt).push_back(bsDaughter.pt());
			(NTuple->gen_BsMuonM_eta).push_back(bsDaughter.eta());
			(NTuple->gen_BsMuonM_phi).push_back(bsDaughter.phi());
		    }
		    if(bsDaughter.pdgId() ==  13)
		    {
			(NTuple->gen_BsMuonP_pt).push_back(bsDaughter.pt());
			(NTuple->gen_BsMuonP_eta).push_back(bsDaughter.eta());
			(NTuple->gen_BsMuonP_phi).push_back(bsDaughter.phi());
		    }
		    if(bsDaughter.pdgId() ==  22)
		    {
			(NTuple->gen_BsPhoton_pt).push_back(bsDaughter.pt());
			(NTuple->gen_BsPhoton_eta).push_back(bsDaughter.eta());
			(NTuple->gen_BsPhoton_phi).push_back(bsDaughter.phi());
	
		    }
		}
	
		(NTuple->gen_Bs_pt).push_back(aBsMeson.pt());
		(NTuple->gen_Bs_eta).push_back(aBsMeson.eta());
		(NTuple->gen_Bs_phi).push_back(aBsMeson.phi());
		(NTuple->gen_Bs_pz).push_back(aBsMeson.pz());
		(NTuple->gen_Bs_pdgId).push_back(aBsMeson.pdgId());
	
		break;
	
        }
	(NTuple->gen_BsMuonMMultiplicity).push_back(muMMul);
	(NTuple->gen_BsMuonPMultiplicity).push_back(muPMul);
	(NTuple->gen_BsPhotonMultiplicity).push_back(phoMul);

	bool hasAValidMCCandidate= true;
	if(muMMul!=1 or muPMul!=1 ) 
	{
		if(printMsg) std::cout<<" Ghost event found !! Mu+ Mu- from any of the Bs not found to == 1 "<<std::endl;
		hasAValidMCCandidate = false;
	}
	else if(doBsToMuMuGamma and phoMul!=1)
	{
		if(printMsg) std::cout<<" Ghost event found !! gamma multiplicity from any of the Bs not found to == 1 "<<std::endl;
		hasAValidMCCandidate = false;
	}

        (NTuple->gen_hasAValid_candidate).push_back(hasAValidMCCandidate);
    }
    //  Muon Ntuplizing
    //  TODO : Add details to closest PV
    //         Add details with BS
    //         Add distance of dimuon vertex to BS and its err
    //

    
    int ii=0,jj=0;
    reco::TrackRef aMuonInnerTrack,bMuonInnerTrack;
    NTuple->nMuons = muons->size();
    for ( uint32_t i=0;i<muons->size();i++) 
    {
    	 auto &aMuon=muons->at(i);
         
         if(aMuon.pt()  < pTMinMuons) continue;
	 aMuonInnerTrack= aMuon.innerTrack();
	 
	 (NTuple->muon_pt		      ).push_back(aMuon.pt());
         (NTuple->muon_eta		      ).push_back(aMuon.eta());
         (NTuple->muon_phi		      ).push_back(aMuon.phi());
	 if(not aMuonInnerTrack.isNull())
	 {
         	(NTuple->mum_dz		              ).push_back(aMuonInnerTrack->dz());
         	(NTuple->muon_dxy		      ).push_back(aMuonInnerTrack->dxy());
         	(NTuple->mum_dz_error	              ).push_back(aMuonInnerTrack->dzError() );
         	(NTuple->muon_dxy_error	              ).push_back(aMuonInnerTrack->dxyError());
         }
	 else
	 {
	         (NTuple->mum_dz		              ).push_back(1e5);
	         (NTuple->muon_dxy		      ).push_back(1e5);
	         (NTuple->mum_dz_error	              ).push_back(1e9);
	         (NTuple->muon_dxy_error	              ).push_back(1e9);
	 }

	 (NTuple->muon_vx		      ).push_back(aMuon.vx() );
         (NTuple->muon_vy		      ).push_back(aMuon.vy());
         (NTuple->muon_vz		      ).push_back(aMuon.vz());
         (NTuple->muon_vertexChi2	      ).push_back(aMuon.vertexChi2());
         (NTuple->muon_vertexNDoF	      ).push_back(aMuon.vertexNdof());
         (NTuple->muon_charge	              ).push_back(aMuon.charge());
         (NTuple->muon_isGlobalMuon	      ).push_back(aMuon.isGlobalMuon());
         (NTuple->muon_isTrackerMuon	      ).push_back(aMuon.isTrackerMuon());
         (NTuple->muon_StandAloneMuon          ).push_back(aMuon.isStandAloneMuon());
         (NTuple->muon_isCaloMuon	      ).push_back(aMuon.isCaloMuon());
         (NTuple->muon_isPFMuon	              ).push_back(aMuon.isPFMuon());

         (NTuple->muon_selector            ).push_back(aMuon.selectors()); 
         (NTuple->muon_isIsolationValid    ).push_back(aMuon.isIsolationValid());
         (NTuple->muon_isPFIsolationValid  ).push_back(aMuon.isIsolationValid());

    auto &MuIsol03 = aMuon.isolationR03();
    (NTuple->muon_isolationR03_trackSumPt).push_back(MuIsol03.sumPt);
    (NTuple->muon_isolationR03_trackEcalSumEt).push_back(MuIsol03.emEt);
    (NTuple->muon_isolationR03_trackHcalSumEt).push_back(MuIsol03.hadEt);
    (NTuple->muon_isolationR03_trackHoSumEt).push_back(MuIsol03.hoEt);
    (NTuple->muon_isolationR03_trackNTracks).push_back(MuIsol03.nTracks);
    (NTuple->muon_isolationR03_trackNJets).push_back(MuIsol03.nJets);
    (NTuple->muon_isolationR03_trackerVetoSumPt).push_back(MuIsol03.trackerVetoPt);
    (NTuple->muon_isolationR03_emVetoSumEt).push_back(MuIsol03.emVetoEt);
    (NTuple->muon_isolationR03_hadVetoSumEt).push_back(MuIsol03.hadVetoEt);
    (NTuple->muon_isolationR03_hoVetoEt).push_back(MuIsol03.hoVetoEt);
  
    auto &MuIsol05 = aMuon.isolationR05();
    (NTuple->muon_isolationR05_trackSumPt).push_back(MuIsol05.sumPt);
    (NTuple->muon_isolationR05_trackEcalSumEt).push_back(MuIsol05.emEt);
    (NTuple->muon_isolationR05_trackHcalSumEt).push_back(MuIsol05.hadEt);
    (NTuple->muon_isolationR05_trackHoSumEt).push_back(MuIsol05.hoEt);
    (NTuple->muon_isolationR05_trackNTracks).push_back(MuIsol05.nTracks);
    (NTuple->muon_isolationR05_trackNJets).push_back(MuIsol05.nJets);
    (NTuple->muon_isolationR05_trackerVetoSumPt).push_back(MuIsol05.trackerVetoPt);
    (NTuple->muon_isolationR05_emVetoSumEt).push_back(MuIsol05.emVetoEt);
    (NTuple->muon_isolationR05_hadVetoSumEt).push_back(MuIsol05.hadVetoEt);
    (NTuple->muon_isolationR05_hoVetoEt).push_back(MuIsol05.hoVetoEt);
  
    auto &MuPFIsol = aMuon.pfIsolationR03();
    (NTuple->muon_PFIsolationR03_sumChargedHadronPt).push_back(MuPFIsol.sumChargedHadronPt);
    (NTuple->muon_PFIsolationR03_sumChargedParticlePt).push_back(MuPFIsol.sumChargedParticlePt);
    (NTuple->muon_PFIsolationR03_sumNeutralHadronEt).push_back(MuPFIsol.sumNeutralHadronEt);
    (NTuple->muon_PFIsolationR03_sumPhotonEt).push_back(MuPFIsol.sumPhotonEt);
    (NTuple->muon_PFIsolationR03_sumNeutralHadronEtHighThreshold).push_back(MuPFIsol.sumNeutralHadronEtHighThreshold);
    (NTuple->muon_PFIsolationR03_sumPhotonEtHighThreshold).push_back(MuPFIsol.sumPhotonEtHighThreshold);
    (NTuple->muon_PFIsolationR03_sumPUPt).push_back(MuPFIsol.sumPUPt);
    
    auto &MuPFIsol04 = aMuon.pfIsolationR04();
    (NTuple->muon_PFIsolationR04_sumChargedHadronPt).push_back(MuPFIsol04.sumChargedHadronPt);
    (NTuple->muon_PFIsolationR04_sumChargedParticlePt).push_back(MuPFIsol04.sumChargedParticlePt);
    (NTuple->muon_PFIsolationR04_sumNeutralHadronEt).push_back(MuPFIsol04.sumNeutralHadronEt);
    (NTuple->muon_PFIsolationR04_sumPhotonEt).push_back(MuPFIsol04.sumPhotonEt);
    (NTuple->muon_PFIsolationR04_sumNeutralHadronEtHighThreshold).push_back(MuPFIsol04.sumNeutralHadronEtHighThreshold);
    (NTuple->muon_PFIsolationR04_sumPhotonEtHighThreshold).push_back(MuPFIsol04.sumPhotonEtHighThreshold);
    (NTuple->muon_PFIsolationR04_sumPUPt).push_back(MuPFIsol04.sumPUPt);
    

  // ### mu- ###
  //mumHighPurity    = nullptr;
  //mumCL            = nullptr;
	 
	 if ((aMuonInnerTrack.isNull() == true))
	 {
		(NTuple->muon_dcaToBS).push_back(0.0);
		(NTuple->muon_dcaToBS_error).push_back(1e4);
	    continue;
	 }
	 const reco::TransientTrack muTrackmTT( aMuonInnerTrack, &(*bFieldHandle));
  	 if (!muTrackmTT.isValid()) continue;
	
	 // # Compute mu- DCA to BeamSpot #
	 theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
	 if (theDCAXBS.isValid() == false)
	 {
	      if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
	            continue;
	 }

        double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
        double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
        if (true /* or fabs(DCAmumBS) > dcaMax_muon_bs*/)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
    //      continue;
        }

	(NTuple->muon_dcaToBS).push_back(DCAmumBS);
	(NTuple->muon_dcaToBS_error).push_back(DCAmumBSErr);
	 ii++;
	 if ((aMuonInnerTrack.isNull() == true)) continue;
	 if ((aMuonInnerTrack->charge() != -1)) continue;
	
	jj=0;
    	for ( uint32_t j=0;j<muons->size();j++) 
    	{
    	 	auto &bMuon=muons->at(j);
               
  	       if(bMuon.pt()  < pTMinMuons) continue;
	       jj++;
	       
	       bMuonInnerTrack = bMuon.innerTrack();
               if ((bMuonInnerTrack.isNull() == true) || (bMuonInnerTrack->charge() != 1)) continue;
	       
		const reco::TransientTrack muTrackpTT(bMuonInnerTrack, &(*bFieldHandle));
                if (!muTrackpTT.isValid()) continue;
	       
	       // # Compute mu+ DCA to BeamSpot #
               theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
               
		if (theDCAXBS.isValid() == false)
               {
                 if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
                 continue;
               }

               double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
               double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
               if (fabs(DCAmupBS) > dcaMax_muon_bs)
               {
                 if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS << std::endl;
                 continue;
               }
	          
	    // # Check goodness of muons closest approach #
               ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
	       XingPoint = ClosestApp.crossingPoint();

	       (NTuple->dimuon_vx).push_back(XingPoint.x());
	       (NTuple->dimuon_vy).push_back(XingPoint.y());
	       (NTuple->dimuon_vz).push_back(XingPoint.z());

            if (  (sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > trackIP_rMax_muon) || (fabs(XingPoint.z()) > trackIP_zMax_muon) )
             {

                if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
          //      continue ;
             }

	    double mumuDCA = ClosestApp.distance();
            
	    if (mumuDCA > dcaMax_dimuon_mumu)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA << std::endl;
        //      continue;
            }
            
	    // # Cut on the dimuon inviariant mass and pT #
            bsDimuon_lv.SetPxPyPzE( muTrackmTT.track().px() + muTrackpTT.track().px(), 
                                muTrackmTT.track().py() + muTrackpTT.track().py(),
                                muTrackmTT.track().pz() + muTrackpTT.track().pz(),
                                sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ) + sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) )
                              );

	    if (true or (bsDimuon_lv.Pt() < minDimuon_pt)  || (bsDimuon_lv.M() < minDimuonInvariantMass) || (bsDimuon_lv.M() > maxDimuonInvariantMass))
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << bsDimuon_lv.Pt() << "\tinv. mass: " << bsDimuon_lv.M() << std::endl;
      //        continue;
            }
	    
	    (NTuple->dimuon_MuPIdx).push_back(jj-1);
	    (NTuple->dimuon_MuMIdx).push_back(ii-1);

	    (NTuple->dimuon_invMass).push_back(bsDimuon_lv.Mag());
	    (NTuple->dimuon_dcaMuMu).push_back(mumuDCA);
	    (NTuple->dimuon_deltaRMuMu).push_back(deltaR(aMuon,bMuon));
	    (NTuple->dimuon_px).push_back(bsDimuon_lv.Px());
	    (NTuple->dimuon_py).push_back(bsDimuon_lv.Py());
	    (NTuple->dimuon_pz).push_back(bsDimuon_lv.Pz());
	    (NTuple->dimuon_energy).push_back(bsDimuon_lv.Energy());
	    (NTuple->dimuon_pt).push_back(bsDimuon_lv.Perp());
	    (NTuple->dimuon_eta).push_back(bsDimuon_lv.Eta());
	    (NTuple->dimuon_phi).push_back(bsDimuon_lv.Phi());

  	    chi2 = 0.;
            ndof  = 0.;
            // ####################################################
            // # Try to vertex the two muons to get dimuon vertex #
            // ####################################################
            std::vector<RefCountedKinematicParticle> muonParticles;
            muonParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi2,ndof,muonMassErr));
            muonParticles.push_back(partFactory.particle(muTrackpTT, muonMass,chi2,ndof,muonMassErr));
        
            mumuVertexFitTree = PartVtxFitter.fit(muonParticles); 
	    if (mumuVertexFitTree->isValid() == false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
             continue; 
            }
        
            mumuVertexFitTree->movePointerToTheTop();
            RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
	    
	   ( NTuple->dimuon_isGoodVertexFit).push_back(mumu_KV->vertexIsValid() );

            if (mumu_KV->vertexIsValid()== false)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
            (NTuple->dimuon_vertex_chi2).push_back(1e7);
	    (NTuple->dimuon_vertex_ndof).push_back(0);
	    (NTuple->dimuon_vertex_proba).push_back(0.0);
	    (NTuple->dimuon_ls).push_back(1e5);
	    (NTuple->dimuon_ls_error).push_back(1e7);
            (NTuple->dimuon_cosAlphaBS).push_back(2);
            (NTuple->dimuon_cosAlphaBS_error).push_back(1e7);
	    
	      continue;
            }
              

            (NTuple->dimuon_vertex_chi2).push_back(mumu_KV->chiSquared());
	    (NTuple->dimuon_vertex_ndof).push_back(mumu_KV->degreesOfFreedom());
	    (NTuple->dimuon_vertex_proba).push_back( TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
            if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) < cl_dimuon_vtx)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: " << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) << std::endl;
             // continue;
            }

            RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();

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

	    (NTuple->dimuon_ls).push_back(MuMuLSBS);
	    (NTuple->dimuon_ls_error).push_back(MuMuLSBSErr);

            if (MuMuLSBS/MuMuLSBSErr < ls_max_dimuonBS)     
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-" << MuMuLSBSErr << std::endl;
             // continue;
            }
            // ###################################################################
            // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
            // ###################################################################
            double MuMuCosAlphaBS;
            double MuMuCosAlphaBSErr;
            Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),mumu_KP->currentState().globalMomentum().y(),0.0,
                                      mumu_KV->position().x() - beamSpot.position().x(),mumu_KV->position().y() - beamSpot.position().y(),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,3),mumu_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                      mumu_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                      mumu_KV->error().cxx() + beamSpot.covariance()(0,0),mumu_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                      mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                      &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
            (NTuple->dimuon_cosAlphaBS).push_back(MuMuCosAlphaBS);
            (NTuple->dimuon_cosAlphaBS_error).push_back(MuMuCosAlphaBSErr);
	    
	    if (MuMuCosAlphaBS < cosAlphaMax_dimuonBs)
            {
              if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS << "+/-" << MuMuCosAlphaBSErr << std::endl;
              //continue;
            }	     
	 
	}
	



   }
 
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
