import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Phase2C9_cff import Phase2C9

process = cms.Process('Demo')

#process = cms.Process("Demo")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.Geometry.GeometryExtended2026D49Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v3', '')

'''
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '92X_dataRun2_Prompt_v4', '')
'''



process.MessageLogger.cerr.FwkReport.reportEvery = 500 

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/CMSSW_11_1_0/src/BsToMuMuGamma_Phase2HLTTDR_MiniAOD_PU140-0B6B07E6-E0DA-7847-87C5-BE84A4E263D5.root'
        'file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/CMSSW_11_1_0/src/BsToMuMuGamma_Phase2HLTTDR_MiniAOD_PU200-844A4C17-49D4-C64F-A0EB-1004E7DF7524.root'
    )
)
process.TFileService = cms.Service("TFileService",
    fileName = cms.string("bsmumuGamma_ntuple_PU200.root"))

# Additional output definition

# Other statements
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '110X_mcRun4_realistic_v3', '')

process.demo = cms.EDAnalyzer("MiniAnalyzer",
    MuonMass = cms.untracked.double(0.10565837), 
    MuonMassErr = cms.untracked.double(3.5e-9),   
    BsMass = cms.untracked.double(5.36677),          ## put the Bs Mass (pdg value)

    MuonMinPt = cms.untracked.double(3.0), # 3.0 [GeV]
    MuonMaxEta = cms.untracked.double(999),  
    MuonMaxDcaBs = cms.untracked.double(2.0), # [cm]

    MuMuMinVtxCl = cms.untracked.double(0.10), # 0.05
    MuMuMinLxySigmaBs = cms.untracked.double(3.0), 
    MuMuMaxDca = cms.untracked.double(0.5), # [cm]
    MuMuMinCosAlphaBs = cms.untracked.double(0.9),

    BsMinVtxCl = cms.untracked.double(0.01), 
    BsMinMass = cms.untracked.double(4.5), # [GeV/c2] 
    BsMaxMass = cms.untracked.double(6.5), # [GeV/c2]  

    bsSrc     = cms.untracked.InputTag("offlineBeamSpot"),
    prunedSrc = cms.untracked.InputTag("prunedGenParticles"),
    genSrc    = cms.untracked.InputTag("packedGenParticles"),
    vtxSrc    = cms.untracked.InputTag("offlineSlimmedPrimaryVertices"),
    photonSrc = cms.untracked.InputTag("slimmedPhotons"),
    muonSrc   = cms.untracked.InputTag("slimmedMuons"),
    trkSrc    = cms.untracked.InputTag("packedPFCandidates"),
    #elecSrc = cms.untracked.InputTag("slimmedElectrons"),
)

process.p = cms.Path(process.demo)
