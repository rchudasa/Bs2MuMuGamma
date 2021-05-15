#ifndef TREECONTENT_H
#define TREECONTENT_H

#include <vector>
#include <string>
#include "TTree.h"

#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"


class TreeContent
{
 public:
  
   TreeContent ();
  ~TreeContent ();

  void Init ();
  void ClearNTuple ();
  void ClearScalars();
  void ClearScalarsMonteCarlo();
  void ClearVectors();
  void ClearVectorsMonteCarlo();
  void ClearMonteCarlo ();
  void MakeTreeBranches (TTree* theTree);


  // ########################################################
  // # Run Number, event number, #reco vtx and event weight #
  // ########################################################
  unsigned int              runN;
  unsigned int              eventN;
  unsigned int              recoVtxN;
  double                    evWeight;
  double                    evWeightE2;
  unsigned int              numEventsTried;
  unsigned int              numEventsPassed;


  // ###########
  // # Trigger #
  // ###########
  std::vector<std::string>  TrigTable;
  std::vector<int>          TrigPrescales;
  std::vector<std::string>  L1Table;
  std::vector<int>          L1Prescales;

  
  
  // ###########################
  // # Number of B0 candidates #
  // ###########################
  unsigned int              nB;
  
  // ############################
  // # Pileup information in MC #
  // ############################
  std::vector<double>       bunchXingMC, numInteractionsMC, trueNumInteractionsMC;
  // Comment:
  // - PileupSummaryInfo::getTrueNumInteractions() gives the distribution of the mean number of interactions per crossing.
  // Since this is the mean value of the poisson distribution from which the number of interactions in- and out-of-time are
  // generated, no additional information should be required for reweighting if these values are matched in data and Monte Carlo.
  // - PileupSummaryInfo::getPU_NumInteractions() gives the expected mean number of interactions per crossing for each LumiSection.
  // Therefore the pileup histogram will contain the distribution of the number of interactions one would actually observe given
  // a poisson of that mean. So, this distribution is what one would see if one counted the number of events seen in a given beam
  // crossing (by looking at the number of vertices in data, for example. This would be appropriate for pileup reweighting based
  // on in-time-only distributions.

  // ################################
  // # Primary Vertex and Beam Spot #
  // ################################
  double                    bsX, bsY;

  
  // # reco::Muons #
  
  std::vector<double>   muon_pt, muon_eta, muon_phi, mum_dz, muon_dxy;
  std::vector<double>   mum_dz_error, muon_dxy_error;
  std::vector<double>   muon_vx,muon_vy,muon_vz,muon_vertexChi2,muon_vertexNDoF;

  std::vector<int>	muon_charge;
  std::vector<bool>	muon_isGlobalMuon,muon_isTrackerMuon,muon_StandAloneMuon,muon_isCaloMuon,muon_isPFMuon;
  int 	            	nMuons;

  std::vector<uint64_t> muon_selector; 
  std::vector<bool>	muon_isIsolationValid;
  std::vector<bool>	muon_isPFIsolationValid;
  
  std::vector<double>  muon_isolationR03_trackSumPt;
  std::vector<double>  muon_isolationR03_trackEcalSumEt;
  std::vector<double>  muon_isolationR03_trackHcalSumEt;
  std::vector<double>  muon_isolationR03_trackHoSumEt;
  std::vector<int>     muon_isolationR03_trackNTracks;
  std::vector<int>     muon_isolationR03_trackNJets;
  std::vector<double>  muon_isolationR03_trackerVetoSumPt;
  std::vector<double>  muon_isolationR03_emVetoSumEt;
  std::vector<double>  muon_isolationR03_hadVetoSumEt;
  std::vector<double>  muon_isolationR03_hoVetoEt;
  
  std::vector<double>  muon_isolationR05_trackSumPt;
  std::vector<double>  muon_isolationR05_trackEcalSumEt;
  std::vector<double>  muon_isolationR05_trackHcalSumEt;
  std::vector<double>  muon_isolationR05_trackHoSumEt;
  std::vector<int>     muon_isolationR05_trackNTracks;
  std::vector<int>     muon_isolationR05_trackNJets;
  std::vector<double>  muon_isolationR05_trackerVetoSumPt;
  std::vector<double>  muon_isolationR05_emVetoSumEt;
  std::vector<double>  muon_isolationR05_hadVetoSumEt;
  std::vector<double>  muon_isolationR05_hoVetoEt;
  
  std::vector<double>  muon_PFIsolationR03_sumChargedHadronPt;
  std::vector<double>  muon_PFIsolationR03_sumChargedParticlePt;
  std::vector<double>  muon_PFIsolationR03_sumNeutralHadronEt;
  std::vector<double>  muon_PFIsolationR03_sumPhotonEt;
  std::vector<double>  muon_PFIsolationR03_sumNeutralHadronEtHighThreshold;
  std::vector<double>  muon_PFIsolationR03_sumPhotonEtHighThreshold;
  std::vector<double>  muon_PFIsolationR03_sumPUPt;
  
  std::vector<double>  muon_PFIsolationR04_sumChargedHadronPt;
  std::vector<double>  muon_PFIsolationR04_sumChargedParticlePt;
  std::vector<double>  muon_PFIsolationR04_sumNeutralHadronEt;
  std::vector<double>  muon_PFIsolationR04_sumPhotonEt;
  std::vector<double>  muon_PFIsolationR04_sumNeutralHadronEtHighThreshold;
  std::vector<double>  muon_PFIsolationR04_sumPhotonEtHighThreshold;
  std::vector<double>  muon_PFIsolationR04_sumPUPt;
  
  std::vector<double> muon_dcaToBS;
  std::vector<double> muon_dcaToBS_error;
  
 // # Dimuon # // 
  
  std::vector<double>   dimuon_pt, dimuon_eta, dimuon_phi, dimuon_energy, dimuon_px,dimuon_py,dimuon_pz;
  std::vector<double>   dimuon_vx,dimuon_vy,dimuon_vz,dimuon_vertex_chi2,dimuon_vertex_ndof,dimuon_vertex_proba;
  std::vector<double>   dimuon_invMass,dimuon_dcaMuMu,dimuon_deltaRMuMu,dimuon_ls,dimuon_ls_error,dimuon_cosAlphaBS,dimuon_cosAlphaBS_error;
  std::vector<int>	dimuon_MuPIdx,dimuon_MuMIdx;
  std::vector<bool>     dimuon_isGoodVertexFit;
  
  // # BeamSpot # //

  double beamspot_x,beamspot_y,beamspot_z,beamspot_x_error,beamspot_y_error,beamspot_z_error;
  double beamspot_dxdz,beamspot_dydz,beamspot_sigmaZ,beamspot_dxdz_error,beamspot_dydz_error,beamspot_sigmaZError;
  double beamspot_beamWidthX,beamspot_beamWidthY,beamspot_beamWidthX_error,beamspot_beamWidthY_error;
 
  // # offlinePrimaryVertices # //
  std::vector<bool> primaryVertex_isFake;
  std::vector<double> primaryVertex_x, primaryVertex_y,primaryVertex_z,primaryVertex_t;
  std::vector<double> primaryVertex_x_error, primaryVertex_y_error,primaryVertex_z_error,primaryVertex_t_error;
  std::vector<double> primaryVertex_ntracks,primaryVertex_ndof,primaryVertex_chi2,primaryVertex_normalizedChi2;
};

#endif
