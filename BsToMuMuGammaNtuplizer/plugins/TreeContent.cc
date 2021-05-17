#include "Bs2MuMuGamma/BsToMuMuGammaNtuplizer/interface/TreeContent.h"
#include <iostream>

TreeContent::TreeContent ()
{
  ClearScalars();
}

TreeContent::~TreeContent ()
{
}

void TreeContent::ClearScalars ()
{
   nMuons=0;
      
      // ## BEAMSOPT STUFF  ## //
   beamspot_x  = 0.0   ;
   beamspot_y  = 0.0   ;
   beamspot_z  = 0.0   ;
   beamspot_x_error  = 0.0   ;
   beamspot_y_error  = 0.0   ;
   beamspot_z_error  = 0.0   ;
   beamspot_dxdz  = 0.0   ;
   beamspot_dydz  = 0.0   ;
   beamspot_sigmaZ  = 0.0   ;
   beamspot_dxdz_error  = 0.0   ;
   beamspot_dydz_error  = 0.0   ;
   beamspot_sigmaZError  = 0.0   ;
   beamspot_beamWidthX  = 0.0   ;
   beamspot_beamWidthY  = 0.0   ;
   beamspot_beamWidthX_error  = 0.0   ;
   beamspot_beamWidthY_error  = 0.0   ;
}


void TreeContent::ClearVectors ()
{
  // ### Trigger ###
  TrigTable.clear();
  TrigPrescales.clear();
  L1Table.clear();
  L1Prescales.clear();
 
   muon_dcaToBS.clear();
   muon_dcaToBS_error.clear();
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
 	
	// ## DIMUON STUFF ## //
	
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

  
  // # offlinePrimaryVertices # //
  
  primaryVertex_isFake.clear();
  primaryVertex_x.clear();
  primaryVertex_y.clear();
  primaryVertex_z.clear();
  primaryVertex_t.clear();
  primaryVertex_x_error.clear();
  primaryVertex_y_error.clear();
  primaryVertex_z_error.clear();
  primaryVertex_t_error.clear();
  primaryVertex_ntracks.clear();
  primaryVertex_ndof.clear();
  primaryVertex_chi2.clear();
  primaryVertex_normalizedChi2.clear();
  
 // # simGen Particles //
 
  gen_hasAValid_candidate.clear();
  gen_Bs_pt.clear() ;
  gen_Bs_eta.clear() ;
  gen_Bs_phi.clear() ;
  gen_Bs_pz.clear() ;
  gen_Bs_pdgId.clear();
  gen_BsMuonM_pt.clear() ;
  gen_BsMuonM_eta.clear() ;
  gen_BsMuonM_phi.clear();
  gen_BsMuonP_pt.clear() ;
  gen_BsMuonP_eta.clear() ;
  gen_BsMuonP_phi.clear();
  gen_BsPhoton_pt.clear() ;
  gen_BsPhoton_eta.clear() ;
  gen_BsPhoton_phi.clear();
  gen_BsPhotonMultiplicity.clear() ;
  gen_BsMuonMMultiplicity.clear() ;
  gen_BsMuonMPultiplicity.clear();


  ClearVectorsMonteCarlo();
}

void TreeContent::ClearScalarsMonteCarlo ()
{
  // ### Matching Between Reconstructed and Generated ###

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

   theTree->Branch("beamspot_x"		 	,&beamspot_x		    			);
   theTree->Branch("beamspot_y"			,&beamspot_y		    			);
   theTree->Branch("beamspot_z"			,&beamspot_z		    			);
   theTree->Branch("beamspot_x_error"		,&beamspot_x_error	    			);
   theTree->Branch("beamspot_y_error"		,&beamspot_y_error	    			);
   theTree->Branch("beamspot_z_error"		,&beamspot_z_error	    			);
   theTree->Branch("beamspot_dxdz"		,&beamspot_dxdz		    			);
   theTree->Branch("beamspot_dydz"		,&beamspot_dydz		    			);
   theTree->Branch("beamspot_sigmaZ"		,&beamspot_sigmaZ	    			);
   theTree->Branch("beamspot_dxdz_error"	,&beamspot_dxdz_error	    			);
   theTree->Branch("beamspot_dydz_error"	,&beamspot_dydz_error	    			);
   theTree->Branch("beamspot_sigmaZError"	,&beamspot_sigmaZError	    			);
   theTree->Branch("beamspot_beamWidthX"	,&beamspot_beamWidthX	    			);
   theTree->Branch("beamspot_beamWidthY"	,&beamspot_beamWidthY	    			);
   theTree->Branch("beamspot_beamWidthX_error"	,&beamspot_beamWidthX_error			);
   theTree->Branch("beamspot_beamWidthY_error"	,&beamspot_beamWidthY_error			);


  // # offlinePrimaryVertex  # //
  theTree->Branch("primaryVertex_isFake"	,&primaryVertex_isFake		   );
  theTree->Branch("primaryVertex_x"		,&primaryVertex_x		   );
  theTree->Branch("primaryVertex_y"		,&primaryVertex_y		   );
  theTree->Branch("primaryVertex_z"		,&primaryVertex_z		   );
  theTree->Branch("primaryVertex_t"		,&primaryVertex_t		   );
  theTree->Branch("primaryVertex_x_error"	,&primaryVertex_x_error		   );
  theTree->Branch("primaryVertex_y_error"	,&primaryVertex_y_error		   );
  theTree->Branch("primaryVertex_z_error"	,&primaryVertex_z_error		   );
  theTree->Branch("primaryVertex_t_error"	,&primaryVertex_t_error		   );
  theTree->Branch("primaryVertex_ntracks"	,&primaryVertex_ntracks		   );
  theTree->Branch("primaryVertex_ndof"		,&primaryVertex_ndof		   );
  theTree->Branch("primaryVertex_chi2"		,&primaryVertex_chi2		   );
  theTree->Branch("primaryVertex_normalizedChi2",&primaryVertex_normalizedChi2	   );
 
  theTree->Branch("gen_hasAValid_candidate"	,&gen_hasAValid_candidate	   );
  theTree->Branch("gen_Bs_pt"			,&gen_Bs_pt			   );
  theTree->Branch("gen_Bs_eta"			,&gen_Bs_eta			   );
  theTree->Branch("gen_Bs_phi"			,&gen_Bs_phi			   );
  theTree->Branch("gen_Bs_pz"			,&gen_Bs_pz			   );
  theTree->Branch("gen_Bs_pdgId"		,&gen_Bs_pdgId			   );
  theTree->Branch("gen_BsMuonM_pt"		,&gen_BsMuonM_pt		   );
  theTree->Branch("gen_BsMuonM_eta"		,&gen_BsMuonM_eta		   );
  theTree->Branch("gen_BsMuonM_phi"		,&gen_BsMuonM_phi		   );
  theTree->Branch("gen_BsMuonP_pt"		,&gen_BsMuonP_pt		   );
  theTree->Branch("gen_BsMuonP_eta"		,&gen_BsMuonP_eta		   );
  theTree->Branch("gen_BsMuonP_phi"		,&gen_BsMuonP_phi		   );
  theTree->Branch("gen_BsPhoton_pt"		,&gen_BsPhoton_pt		   );
  theTree->Branch("gen_BsPhoton_eta"		,&gen_BsPhoton_eta		   );
  theTree->Branch("gen_BsPhoton_phi"		,&gen_BsPhoton_phi		   );
  theTree->Branch("gen_BsPhotonMultiplicity"	,&gen_BsPhotonMultiplicity	   );
  theTree->Branch("gen_BsMuonMMultiplicity"	,&gen_BsMuonMMultiplicity	   );
  theTree->Branch("gen_BsMuonMPultiplicity"	,&gen_BsMuonMPultiplicity	   );

  

}




