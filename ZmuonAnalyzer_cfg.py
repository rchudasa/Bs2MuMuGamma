import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")

process.ZmuonAnalyzer = cms.EDAnalyzer("ZmuonAnalyzer",
   muonCollection = cms.InputTag("slimmedMuons"),
   electronCollection = cms.InputTag("slimmedElectrons"),
   bits = cms.InputTag("TriggerResults","", "HLT"),
   objects = cms.InputTag("selectedPatTrigger"),
   genParticles = cms.InputTag("prunedGenParticles"),
   vertexCollection = cms.InputTag("offlineSlimmedPrimaryVertices"),
   metTag = cms.InputTag("slimmedMETs"),
   pfCands = cms.InputTag("packedPFCandidates")

)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
                   "file:/afs/cern.ch/work/r/rchudasa/private/bsmumu/CMSSW_11_1_0/src/BsToMuMuGamma_Phase2HLTTDR_MiniAOD_PU140-0B6B07E6-E0DA-7847-87C5-BE84A4E263D5.root"
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU140_110X_mcRun4_realistic_v3-v2/100000/0B6B07E6-E0DA-7847-87C5-BE84A4E263D5.root",
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU140_110X_mcRun4_realistic_v3-v2/100000/51509B9F-0D61-2749-B610-5ADF28F685C6.root",
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU140_110X_mcRun4_realistic_v3-v2/100000/A9DCEFF7-4E7A-3240-994E-931EA0323F39.root", 
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU140_110X_mcRun4_realistic_v3-v2/100000/F85BFDFF-2FAB-4A4A-B0E9-ED5686D1B8CE.root" 

#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU200_110X_mcRun4_realistic_v3-v2/10000/0B378AB3-5632-5449-ABCD-4C7A45AD1377.root",
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU200_110X_mcRun4_realistic_v3-v2/10000/844A4C17-49D4-C64F-A0EB-1004E7DF7524.root", 
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU200_110X_mcRun4_realistic_v3-v2/10000/8CFF9BD5-57D8-5140-B9BA-61236E594D10.root", 
#                  "/store/mc/Phase2HLTTDRWinter20RECOMiniAOD/BsToMuMuGamma_BMuonFilter_SoftQCDnonD_TuneCP5_14TeV-pythia8-evtgen/MINIAODSIM/PU200_110X_mcRun4_realistic_v3-v2/10000/AFDF2012-6A4E-564E-95AB-D0C414AA3299.root" 
                                     )
)

process.TFileService = cms.Service("TFileService",
   fileName = cms.string("muons_tree_filled_iffound3dau_v2.root")
)

process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(10000)
#   input = cms.untracked.int32(-1)
)

process.options = cms.untracked.PSet(
   wantSummary = cms.untracked.bool(True),
)

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

#process.MessageLogger = cms.Service("MessageLogger",
#                    destinations   = cms.untracked.vstring('messages.txt')
#)

process.p = cms.Path(
   process.ZmuonAnalyzer
)

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data')

