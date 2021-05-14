#include "Bs2MuMuGamma/BsToMuMuGammaNtuplizer/interface/TreeContent.h"
#include <iostream>

TreeContent::TreeContent ()
{
  ClearScalars();

  // ### Trigger ###
  TrigTable     = nullptr;
  TrigPrescales = nullptr;
  L1Table       = nullptr;
  L1Prescales   = nullptr;
 // hltObjs       = nullptr;
  
    nMuons=0;
 //   muon_pt		=nullptr	;
 //   muon_eta		=nullptr	;
 //   muon_phi		=nullptr	;
 //   mum_dz		=nullptr	;
 //   muon_dxy		=nullptr	;
 //   mum_dz_error	=nullptr	;
 //   muon_dxy_error	=nullptr	;
 //   muon_vx		=nullptr	;
 //   muon_vy		=nullptr	;
 //   muon_vz		=nullptr	;
 //   muon_vertexChi2	=nullptr	;
 //   muon_vertexNDoF	=nullptr	;
 //   muon_charge	        =nullptr	;
 //   muon_isGlobalMuon	=nullptr	;
 //   muon_isTrackerMuon	=nullptr	;
 //   muon_StandAloneMuon =nullptr	;
 //   muon_isCaloMuon	=nullptr	;
 //   muon_isPFMuon	=nullptr	;

 //   muon_selector         	= nullptr; 
 //   muon_isIsolationValid 	= nullptr;
 //   muon_isPFIsolationValid 	= nullptr;

 //  muon_selector=nullptr; 
 //  muon_isIsolationValid=nullptr;
 //  muon_isPFIsolationValid=nullptr;
 // 
 //  muon_isolationR03_trackSumPt=nullptr;
 //  muon_isolationR03_trackEcalSumEt=nullptr;
 //  muon_isolationR03_trackHcalSumEt=nullptr;
 //  muon_isolationR03_trackHoSumEt=nullptr;
 //  muon_isolationR03_trackNTracks=nullptr;
 //  muon_isolationR03_trackNJets=nullptr;
 //  muon_isolationR03_trackerVetoSumPt=nullptr;
 //  muon_isolationR03_emVetoSumEt=nullptr;
 //  muon_isolationR03_hadVetoSumEt=nullptr;
 //  muon_isolationR03_hoVetoEt=nullptr;
 // 
 //  muon_isolationR05_trackSumPt=nullptr;
 //  muon_isolationR05_trackEcalSumEt=nullptr;
 //  muon_isolationR05_trackHcalSumEt=nullptr;
 //  muon_isolationR05_trackHoSumEt=nullptr;
 //  muon_isolationR05_trackNTracks=nullptr;
 //  muon_isolationR05_trackNJets=nullptr;
 //  muon_isolationR05_trackerVetoSumPt=nullptr;
 //  muon_isolationR05_emVetoSumEt=nullptr;
 //  muon_isolationR05_hadVetoSumEt=nullptr;
 //  muon_isolationR05_hoVetoEt=nullptr;
 // 
 //  muon_PFIsolationR03_sumChargedHadronPt=nullptr;
 //  muon_PFIsolationR03_sumChargedParticlePt=nullptr;
 //  muon_PFIsolationR03_sumNeutralHadronEt=nullptr;
 //  muon_PFIsolationR03_sumPhotonEt=nullptr;
 //  muon_PFIsolationR03_sumNeutralHadronEtHighThreshold=nullptr;
 //  muon_PFIsolationR03_sumPhotonEtHighThreshold=nullptr;
 //  muon_PFIsolationR03_sumPUPt=nullptr;
 //  
 //  muon_PFIsolationR04_sumChargedHadronPt=nullptr;
 //  muon_PFIsolationR04_sumChargedParticlePt=nullptr;
 //  muon_PFIsolationR04_sumNeutralHadronEt=nullptr;
 //  muon_PFIsolationR04_sumPhotonEt=nullptr;
 //  muon_PFIsolationR04_sumNeutralHadronEtHighThreshold=nullptr;
 //  muon_PFIsolationR04_sumPhotonEtHighThreshold=nullptr;
 //  muon_PFIsolationR04_sumPUPt=nullptr;
 

  // ### mu- ###
  //mumHighPurity    = nullptr;
  //mumCL            = nullptr;
  //mumNormChi2      = nullptr;
  //mumPx            = nullptr;
  //mumPy            = nullptr;
  //mumPz            = nullptr;
  //mumDCAVtx        = nullptr;
  //mumDCAVtxE       = nullptr;
  //mumDCABS         = nullptr;
  //mumDCABSE        = nullptr;
  //mumKinkChi2      = nullptr;
  //mumFracHits      = nullptr;
  //mumdxyBS         = nullptr;
  //mumdzBS          = nullptr;
  //mumMinIP2D       = nullptr;
  //mumMinIP2DE      = nullptr;
  //mumMinIP         = nullptr;
  //mumMinIPS        = nullptr;
  //mumDeltaRwithMC  = nullptr;
  //mumCat           = nullptr;
  //mumNPixHits      = nullptr;
  //mumNPixLayers    = nullptr;
  //mumNTrkHits      = nullptr;
  //mumNTrkLayers    = nullptr;
  //mumNMuonHits     = nullptr;
  //mumNMatchStation = nullptr;
  //mumIso           = nullptr;
  //mumIsoPt         = nullptr;
  //mumIsodR         = nullptr;

  }

void TreeContent::Init ()
{
  // ### Trigger ###
 // TrigTable     = new std::vector<std::string>;
 // TrigPrescales = new std::vector<int>;
 // L1Table       = new std::vector<std::string>;
 // L1Prescales   = new std::vector<int>;
 // //hltObjs       = new std::vector<miniHLTObj>;
 // 
 //   nMuons= new int;
 //   muon_pt		= new std::vector<double>	;
 //   muon_eta		= new std::vector<double>	;
 //   muon_phi		= new std::vector<double>	;
 //   mum_dz		= new std::vector<double>	;
 //   muon_dxy		= new std::vector<double>	;
 //   mum_dz_error	= new std::vector<double>	;
 //   muon_dxy_error	= new std::vector<double>	;
 //   muon_vx		= new std::vector<double>	;
 //   muon_vy		= new std::vector<double>	;
 //   muon_vz		= new std::vector<double>	;
 //   muon_vertexChi2	= new std::vector<double>	;
 //   muon_vertexNDoF	= new std::vector<double>	;
 //   muon_charge	        = new std::vector<int>	;
 //   muon_isGlobalMuon	= new std::vector<bool>	;
 //   muon_isTrackerMuon	= new std::vector<bool>	;
 //   muon_StandAloneMuon = new std::vector<bool>	;
 //   muon_isCaloMuon	= new std::vector<bool>	;
 //   muon_isPFMuon	= new std::vector<bool>	;

 //   muon_selector         =  new std::vector<uint64_t>; 
 //
 //   muon_isIsolationValid =  new std::vector<bool>;;
 //   muon_isPFIsolationValid =  new std::vector<bool>;;

 //  muon_isolationR03_trackSumPt=new std::vector<double>;
 //  muon_isolationR03_trackEcalSumEt=new std::vector<double>;
 //  muon_isolationR03_trackHcalSumEt=new std::vector<double>;
 //  muon_isolationR03_trackHoSumEt=new std::vector<double>;
 //  muon_isolationR03_trackNTracks=new std::vector<int>;
 //  muon_isolationR03_trackNJets=new std::vector<int>;
 //  muon_isolationR03_trackerVetoSumPt=new std::vector<double>;
 //  muon_isolationR03_emVetoSumEt=new std::vector<double>;
 //  muon_isolationR03_hadVetoSumEt=new std::vector<double>;
 //  muon_isolationR03_hoVetoEt=new std::vector<double>;
 // 
 //  muon_isolationR05_trackSumPt=new std::vector<double>;
 //  muon_isolationR05_trackEcalSumEt=new std::vector<double>;
 //  muon_isolationR05_trackHcalSumEt=new std::vector<double>;
 //  muon_isolationR05_trackHoSumEt=new std::vector<double>;
 //  muon_isolationR05_trackNTracks=new std::vector<int>;
 //  muon_isolationR05_trackNJets=new std::vector<int>;
 //  muon_isolationR05_trackerVetoSumPt=new std::vector<double>;
 //  muon_isolationR05_emVetoSumEt=new std::vector<double>;
 //  muon_isolationR05_hadVetoSumEt=new std::vector<double>;
 //  muon_isolationR05_hoVetoEt=new std::vector<double>;
 // 
 //  muon_PFIsolationR03_sumChargedHadronPt=new std::vector<double>;
 //  muon_PFIsolationR03_sumChargedParticlePt=new std::vector<double>;
 //  muon_PFIsolationR03_sumNeutralHadronEt=new std::vector<double>;
 //  muon_PFIsolationR03_sumPhotonEt=new std::vector<double>;
 //  muon_PFIsolationR03_sumNeutralHadronEtHighThreshold=new std::vector<double>;
 //  muon_PFIsolationR03_sumPhotonEtHighThreshold=new std::vector<double>;
 //  muon_PFIsolationR03_sumPUPt=new std::vector<double>;
 //  
 //  muon_PFIsolationR04_sumChargedHadronPt=new std::vector<double>;
 //  muon_PFIsolationR04_sumChargedParticlePt=new std::vector<double>;
 //  muon_PFIsolationR04_sumNeutralHadronEt=new std::vector<double>;
 //  muon_PFIsolationR04_sumPhotonEt=new std::vector<double>;
 //  muon_PFIsolationR04_sumNeutralHadronEtHighThreshold=new std::vector<double>;
 //  muon_PFIsolationR04_sumPhotonEtHighThreshold=new std::vector<double>;
 //  muon_PFIsolationR04_sumPUPt=new std::vector<double>;
 

  // ### mu- ###
  //mumHighPurity    = new std::vector<double>;
  //mumCL            = new std::vector<double>;
}

TreeContent::~TreeContent ()
{
   //delete  TrigTable     	;	
   //delete  TrigPrescales 	;
   //delete  L1Table       	;
   //delete  L1Prescales   	;
   ////delete  hltObjs       	;
  
   //delete  nMuons		;
   //delete  muon_pt		;
   //delete  muon_eta		;
   //delete  muon_phi		;
   //delete  mum_dz		;
   //delete  muon_dxy		;
   //delete  mum_dz_error	;
   //delete  muon_dxy_error	;
   //delete  muon_vx		;
   //delete  muon_vy		;
   //delete  muon_vz		;
   //delete  muon_vertexChi2	;
   //delete  muon_vertexNDoF	;
   //delete  muon_charge	        ;
   //delete  muon_isGlobalMuon	;
   //delete  muon_isTrackerMuon	;
   //delete  muon_StandAloneMuon ;
   //delete  muon_isCaloMuon	;
   //delete  muon_isPFMuon	;

   //delete  muon_selector       ;
 
   //delete  muon_isIsolationValid;
   //delete  muon_isPFIsolationValid;

   //delete muon_isolationR03_trackSumPt;
   //delete muon_isolationR03_trackEcalSumEt;
   //delete muon_isolationR03_trackHcalSumEt;
   //delete muon_isolationR03_trackHoSumEt;
   //delete muon_isolationR03_trackNTracks;
   //delete muon_isolationR03_trackNJets;
   //delete muon_isolationR03_trackerVetoSumPt;
   //delete muon_isolationR03_emVetoSumEt;
   //delete muon_isolationR03_hadVetoSumEt;
   //delete muon_isolationR03_hoVetoEt;
  
   //delete muon_isolationR05_trackSumPt;
   //delete muon_isolationR05_trackEcalSumEt;
   //delete muon_isolationR05_trackHcalSumEt;
   //delete muon_isolationR05_trackHoSumEt;
   //delete muon_isolationR05_trackNTracks;
   //delete muon_isolationR05_trackNJets;
   //delete muon_isolationR05_trackerVetoSumPt;
   //delete muon_isolationR05_emVetoSumEt;
   //delete muon_isolationR05_hadVetoSumEt;
   //delete muon_isolationR05_hoVetoEt;
  
   //delete muon_PFIsolationR03_sumChargedHadronPt;
   //delete muon_PFIsolationR03_sumChargedParticlePt;
   //delete muon_PFIsolationR03_sumNeutralHadronEt;
   //delete muon_PFIsolationR03_sumPhotonEt;
   //delete muon_PFIsolationR03_sumNeutralHadronEtHighThreshold;
   //delete muon_PFIsolationR03_sumPhotonEtHighThreshold;
   //delete muon_PFIsolationR03_sumPUPt;
   // 
   //delete muon_PFIsolationR04_sumChargedHadronPt;
   //delete muon_PFIsolationR04_sumChargedParticlePt;
   //delete muon_PFIsolationR04_sumNeutralHadronEt;
   //delete muon_PFIsolationR04_sumPhotonEt;
   //delete muon_PFIsolationR04_sumNeutralHadronEtHighThreshold;
   //delete muon_PFIsolationR04_sumPhotonEtHighThreshold;
   //delete muon_PFIsolationR04_sumPUPt;
 

  // ### mu- ###
  //mumHighPurity    = nullptr;
  //mumCL            = nullptr;

}

void TreeContent::ClearScalars ()
{
  ClearScalarsMonteCarlo();
}

void TreeContent::ClearScalarsMonteCarlo ()
{

}

void TreeContent::ClearVectors ()
{
  // ### Trigger ###
  TrigTable->clear();
  TrigPrescales->clear();
  L1Table->clear();
  L1Prescales->clear();
   //hltObjs -> clear();
 
   muon_dcaToBS.clear();
   muon_dcaToBS_error.clear();
   
  
    dimuon_pt.clear();
    dimuon_eta.clear();
    dimuon_phi.clear();
    dimuon_energy.clear();
    dimuon_px.clear();
    dimuon_py.clear();
    dimuon_pz.clear();
    dimuon_vx.clear();
    dimuon_vy.clear();
    dimuon_vz.clear();
    dimuon_vertex_chi2.clear();
    dimuon_vertex_ndof.clear();
    dimuon_vertex_proba.clear();
    dimuon_invMass.clear();
    dimuon_dcaMuMu.clear();
    dimuon_deltaRMuMu.clear();
    dimuon_ls.clear();
    dimuon_ls_error.clear();
    dimuon_cosAlphaBS.clear();
    dimuon_cosAlphaBS_error.clear();
    dimuon_MuPIdx.clear();
    dimuon_MuMIdx.clear();
    dimuon_isGoodVertexFit.clear();
   
   
   nMuons=0;
   muon_pt		.clear();
   muon_eta		.clear();
   muon_phi		.clear();
   mum_dz		.clear();
   muon_dxy		.clear();
   mum_dz_error	.clear();
   muon_dxy_error	.clear();
   muon_vx		.clear();
   muon_vy		.clear();
   muon_vz		.clear();
   muon_vertexChi2	.clear();
   muon_vertexNDoF	.clear();
   muon_charge	        .clear();
   muon_isGlobalMuon	.clear();
   muon_isTrackerMuon	.clear();
   muon_StandAloneMuon .clear();
   muon_isCaloMuon	.clear();
   muon_isPFMuon	.clear();

   muon_selector       .clear();
 
   muon_isIsolationValid.clear();
   muon_isPFIsolationValid.clear();
   
    muon_isolationR03_trackSumPt.clear();
    muon_isolationR03_trackEcalSumEt.clear();
    muon_isolationR03_trackHcalSumEt.clear();
    muon_isolationR03_trackHoSumEt.clear();
    muon_isolationR03_trackNTracks.clear();
    muon_isolationR03_trackNJets.clear();
    muon_isolationR03_trackerVetoSumPt.clear();
    muon_isolationR03_emVetoSumEt.clear();
    muon_isolationR03_hadVetoSumEt.clear();
    muon_isolationR03_hoVetoEt.clear();
  
    muon_isolationR05_trackSumPt.clear();
    muon_isolationR05_trackEcalSumEt.clear();
    muon_isolationR05_trackHcalSumEt.clear();
    muon_isolationR05_trackHoSumEt.clear();
    muon_isolationR05_trackNTracks.clear();
    muon_isolationR05_trackNJets.clear();
    muon_isolationR05_trackerVetoSumPt.clear();
    muon_isolationR05_emVetoSumEt.clear();
    muon_isolationR05_hadVetoSumEt.clear();
    muon_isolationR05_hoVetoEt.clear();
  
    muon_PFIsolationR03_sumChargedHadronPt.clear();
    muon_PFIsolationR03_sumChargedParticlePt.clear();
    muon_PFIsolationR03_sumNeutralHadronEt.clear();
    muon_PFIsolationR03_sumPhotonEt.clear();
    muon_PFIsolationR03_sumNeutralHadronEtHighThreshold.clear();
    muon_PFIsolationR03_sumPhotonEtHighThreshold.clear();
    muon_PFIsolationR03_sumPUPt.clear();
    
    muon_PFIsolationR04_sumChargedHadronPt.clear();
    muon_PFIsolationR04_sumChargedParticlePt.clear();
    muon_PFIsolationR04_sumNeutralHadronEt.clear();
    muon_PFIsolationR04_sumPhotonEt.clear();
    muon_PFIsolationR04_sumNeutralHadronEtHighThreshold.clear();
    muon_PFIsolationR04_sumPhotonEtHighThreshold.clear();
    muon_PFIsolationR04_sumPUPt.clear();
 

  // ### mu- ###
  //mumHighPurity    = nullptr.clear();
  //mumCL            = nullptr.clear();


  ClearVectorsMonteCarlo();
}

void TreeContent::ClearVectorsMonteCarlo ()
{
  // ### Matching Between Reconstructed and Generated ###
}

void TreeContent::ClearNTuple ()
{
  ClearScalars();
  ClearVectors();
}

void TreeContent::ClearMonteCarlo ()
{
  ClearScalarsMonteCarlo();
  ClearVectorsMonteCarlo();
}

void TreeContent::MakeTreeBranches (TTree* theTree)
{
  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  // theTree->Branch("runN",            &runN,            "runN/i");
  // theTree->Branch("eventN",          &eventN,          "eventN/i");
  // theTree->Branch("recoVtxN",        &recoVtxN,        "recoVtxN/i");
  // theTree->Branch("evWeight",        &evWeight,        "evWeight/D");
  // theTree->Branch("evWeightE2",      &evWeightE2,      "evWeightE2/D");
  // theTree->Branch("numEventsTried",  &numEventsTried,  "numEventsTried/i");
  // theTree->Branch("numEventsPassed", &numEventsPassed, "numEventsPassed/i");

  // ### Trigger ###
  theTree->Branch("TrigTable",     &TrigTable);
  theTree->Branch("TrigPrescales", &TrigPrescales);
  theTree->Branch("L1Table",       &L1Table);
  theTree->Branch("L1Prescales",   &L1Prescales);
  //theTree->Branch("hltObjs"    ,   &hltObjs);

   theTree->Branch("nMuons"			,  nMuons			);
   theTree->Branch("muon_pt"			,  &muon_pt			);
   theTree->Branch("muon_eta"		   	,  &muon_eta			);
   theTree->Branch("muon_phi"		   	,  &muon_phi			);
   theTree->Branch("mum_dz"		   	,  &mum_dz			);
   theTree->Branch("muon_dxy"		   	,  &muon_dxy			);
   theTree->Branch("mum_dz_error"	   	,  &mum_dz_error   		);
   theTree->Branch("muon_dxy_error"	   	,  &muon_dxy_error		);
   theTree->Branch("muon_vx"		   	,  &muon_vx			);
   theTree->Branch("muon_vy"		   	,  &muon_vy			);
   theTree->Branch("muon_vz"		   	,  &muon_vz			);
   theTree->Branch("muon_vertexChi2"	   	,  &muon_vertexChi2		);
   theTree->Branch("muon_vertexNDoF"	   	,  &muon_vertexNDoF		);
   theTree->Branch("muon_charge"	        ,  &muon_charge	        	);
   theTree->Branch("muon_isGlobalMuon"	   	,  &muon_isGlobalMuon		);
   theTree->Branch("muon_isTrackerMuon"	   	,  &muon_isTrackerMuon		);
   theTree->Branch("muon_StandAloneMuon"    	,  &muon_StandAloneMuon 	);
   theTree->Branch("muon_isCaloMuon"	   	,  &muon_isCaloMuon		);
   theTree->Branch("muon_isPFMuon"	   	,  &muon_isPFMuon		);
   theTree->Branch("muon_selector"          	,  &muon_selector         	); 
   theTree->Branch("muon_isIsolationValid"  	,  &muon_isIsolationValid 	);
   theTree->Branch("muon_isPFIsolationValid"	,  &muon_isPFIsolationValid 	);
   
   theTree->Branch("muon_isolationR03_trackSumPt"			,&muon_isolationR03_trackSumPt						);
   theTree->Branch("muon_isolationR03_trackEcalSumEt"			,&muon_isolationR03_trackEcalSumEt					);
   theTree->Branch("muon_isolationR03_trackHcalSumEt"			,&muon_isolationR03_trackHcalSumEt					);
   theTree->Branch("muon_isolationR03_trackHoSumEt"			,&muon_isolationR03_trackHoSumEt					);
   theTree->Branch("muon_isolationR03_trackNTracks"			,&muon_isolationR03_trackNTracks					);
   theTree->Branch("muon_isolationR03_trackNJets"			,&muon_isolationR03_trackNJets						);
   theTree->Branch("muon_isolationR03_trackerVetoSumPt"			,&muon_isolationR03_trackerVetoSumPt					);
   theTree->Branch("muon_isolationR03_emVetoSumEt"			,&muon_isolationR03_emVetoSumEt						);
   theTree->Branch("muon_isolationR03_hadVetoSumEt"			,&muon_isolationR03_hadVetoSumEt					);
   theTree->Branch("muon_isolationR03_hoVetoEt"				,&muon_isolationR03_hoVetoEt						);
                                                                                                                             
   theTree->Branch("muon_isolationR05_trackSumPt"			,&muon_isolationR05_trackSumPt						);
   theTree->Branch("muon_isolationR05_trackEcalSumEt"			,&muon_isolationR05_trackEcalSumEt					);
   theTree->Branch("muon_isolationR05_trackHcalSumEt"			,&muon_isolationR05_trackHcalSumEt					);
   theTree->Branch("muon_isolationR05_trackHoSumEt"			,&muon_isolationR05_trackHoSumEt					);
   theTree->Branch("muon_isolationR05_trackNTracks"			,&muon_isolationR05_trackNTracks					);
   theTree->Branch("muon_isolationR05_trackNJets"			,&muon_isolationR05_trackNJets						);
   theTree->Branch("muon_isolationR05_trackerVetoSumPt"			,&muon_isolationR05_trackerVetoSumPt					);
   theTree->Branch("muon_isolationR05_emVetoSumEt"			,&muon_isolationR05_emVetoSumEt						);
   theTree->Branch("muon_isolationR05_hadVetoSumEt"			,&muon_isolationR05_hadVetoSumEt					);
   theTree->Branch("muon_isolationR05_hoVetoEt"				,&muon_isolationR05_hoVetoEt						);
                                                                                                                             
   theTree->Branch("muon_PFIsolationR03_sumChargedHadronPt"		,&muon_PFIsolationR03_sumChargedHadronPt				);
   theTree->Branch("muon_PFIsolationR03_sumChargedParticlePt"		,&muon_PFIsolationR03_sumChargedParticlePt				);
   theTree->Branch("muon_PFIsolationR03_sumNeutralHadronEt"		,&muon_PFIsolationR03_sumNeutralHadronEt				);
   theTree->Branch("muon_PFIsolationR03_sumPhotonEt"			,&muon_PFIsolationR03_sumPhotonEt					);
   theTree->Branch("muon_PFIsolationR03_sumNeutralHadronEtHighThreshold",&muon_PFIsolationR03_sumNeutralHadronEtHighThreshold			);
   theTree->Branch("muon_PFIsolationR03_sumPhotonEtHighThreshold"	,&muon_PFIsolationR03_sumPhotonEtHighThreshold				);
   theTree->Branch("muon_PFIsolationR03_sumPUPt"			,&muon_PFIsolationR03_sumPUPt						);
                                                                                                                             
   theTree->Branch("muon_PFIsolationR04_sumChargedHadronPt"		,&muon_PFIsolationR04_sumChargedHadronPt				);
   theTree->Branch("muon_PFIsolationR04_sumChargedParticlePt"		,&muon_PFIsolationR04_sumChargedParticlePt				);
   theTree->Branch("muon_PFIsolationR04_sumNeutralHadronEt"		,&muon_PFIsolationR04_sumNeutralHadronEt				);
   theTree->Branch("muon_PFIsolationR04_sumPhotonEt"			,&muon_PFIsolationR04_sumPhotonEt					);
   theTree->Branch("muon_PFIsolationR04_sumNeutralHadronEtHighThreshold",&muon_PFIsolationR04_sumNeutralHadronEtHighThreshold			);
   theTree->Branch("muon_PFIsolationR04_sumPhotonEtHighThreshold"	,&muon_PFIsolationR04_sumPhotonEtHighThreshold				);
   theTree->Branch("muon_PFIsolationR04_sumPUPt"			,&muon_PFIsolationR04_sumPUPt						);
   theTree->Branch("muon_dcaToBS"		     			,&muon_dcaToBS 								);
   theTree->Branch("muon_dcaToBS_error"		     			,&muon_dcaToBS_error 							);
 
  
   theTree->Branch("dimuon_pt"				,&dimuon_pt				);
   theTree->Branch("dimuon_eta"				,&dimuon_eta				);
   theTree->Branch("dimuon_phi"				,&dimuon_phi				);
   theTree->Branch("dimuon_energy"			,&dimuon_energy				);
   theTree->Branch("dimuon_px"				,&dimuon_px				);
   theTree->Branch("dimuon_py"				,&dimuon_py				);
   theTree->Branch("dimuon_pz"				,&dimuon_pz				);
   theTree->Branch("dimuon_vx"				,&dimuon_vx				);
   theTree->Branch("dimuon_vy"				,&dimuon_vy				);
   theTree->Branch("dimuon_vz"				,&dimuon_vz				);
   theTree->Branch("dimuon_vertex_chi2"			,&dimuon_vertex_chi2			);
   theTree->Branch("dimuon_vertex_ndof"			,&dimuon_vertex_ndof			);
   theTree->Branch("dimuon_vertex_proba"		,&dimuon_vertex_proba			);
   theTree->Branch("dimuon_invMass"			,&dimuon_invMass			);
   theTree->Branch("dimuon_deltaRMuMu"			,&dimuon_deltaRMuMu			);
   theTree->Branch("dimuon_ls"				,&dimuon_ls				);
   theTree->Branch("dimuon_ls_error"			,&dimuon_ls_error			);
   theTree->Branch("dimuon_cosAlphaBS"			,&dimuon_cosAlphaBS			);
   theTree->Branch("dimuon_cosAlphaBS_error"		,&dimuon_cosAlphaBS_error            	);
   theTree->Branch("dimuon_MuPIdx"			,&dimuon_MuPIdx				);
   theTree->Branch("dimuon_MuMIdx"			,&dimuon_MuMIdx				);
   theTree->Branch("dimuon_isGoodVertexFit"		,&dimuon_isGoodVertexFit		);
   
   

}




