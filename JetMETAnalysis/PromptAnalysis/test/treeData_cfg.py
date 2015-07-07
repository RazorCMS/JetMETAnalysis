import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process('TEST')
process.load('JetMETAnalysis.PromptAnalysis.ntuple_cff')

process.load("Configuration/StandardSequences/Geometry_cff")
process.load("Configuration/StandardSequences/MagneticField_cff")

###########
process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag ='GR_P_V56::All'##make sure to check the GT from DBS
###########

process.load("RecoMuon/Configuration/RecoMuon_cff")

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
fileName = cms.string("reco_test.root"), ##give a name to the output file
                           closeFileFast = cms.untracked.bool(True)  ) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring(
    "file:eos/cms/store/data/Run2015A/Jet/RECO/PromptReco-v1/000/246/960/00000/4CED9FE9-DB0B-E511-A069-02163E012432.root"
    ),
    
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    
    secondaryFileNames = cms.untracked.vstring())

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100

# jet corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')


# summary
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )


process.promptanaTree = cms.EDAnalyzer("PromptAnaTree",
    outputCommands = cms.untracked.vstring(
    'drop *',
    'keep *_promptanaevent_*_*',
    'keep *_promptanamet_*_*',
    'keep *_promptanapfmet_*_*',
    'keep *_promptanaak5pfjet_*_*',
    'keep *_promptanacalotowers_*_*',
    'keep *_promptanavtx_*_*',
    'keep *_promptanatrack_*_*',
    'keep *_promptanaecalspikes_*_*'
    ))

process.theBigNtuple = cms.Path(
    (
    process.promptanaevent +
    process.promptanamet   +
    process.promptanapfmet   +
    process.promptanaak5pfjet +
    process.promptanacalotowers +
    process.promptanavtx +
    process.promptanatrack +
    process.promptanacleanup +
    process.promptanaecalspikes
    )
    * process.promptanaTree )

