import FWCore.ParameterSet.Config as cms
import HLTrigger.HLTfilters.hltHighLevel_cfi as hlt

process = cms.Process('USER')

process.load('JetMETAnalysis.PromptAnalysis.ntuple_cff')

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")

#########################
### MET cleaning
process.load('RecoMET.METAnalyzers.CSCHaloFilter_cfi')


#########################


process.GlobalTag.globaltag = 'FT_R_42_V13A::All'

process.load('Configuration.StandardSequences.Services_cff')
process.add_( cms.Service( "TFileService",
fileName = cms.string("RelValTTbar_CMSSW_3_8_4_patch2-START38_V12_APDHits_ECALCleaningUndone.root"),
                           closeFileFast = cms.untracked.bool(True)  ) )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(200) )
process.source = cms.Source (
    "PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:/uscms_data/d1/apresyan/TTbar_7TeV_GEN-SIM-RECO_99_3_Tng.root'
        # 'file:/uscms_data/d1/apresyan/pickevents_spikesTime.root'
         '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FC5E445A-397C-E011-908E-0024E87663FB.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FC575FAF-1B7C-E011-A689-0015178C4B78.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FC5685A9-427C-E011-A69F-00266CF279F8.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FAF468E6-3F7C-E011-96E7-0015178C69B8.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FAD88D60-037C-E011-955C-001D0967CFCC.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FA6487B2-C77B-E011-B499-0015178C0218.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FA524CFC-F27B-E011-894A-00266CF253C4.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/FA279DC0-037C-E011-B14C-00A0D1EEF364.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F8F78004-2C7C-E011-A128-0024E87683F8.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F859F3EE-F47B-E011-B54A-00A0D1EE968C.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F824813A-F17B-E011-B249-0015178C495C.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F675FF4D-D47B-E011-802B-0015178C4D14.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F62804FB-3C7C-E011-BB95-00A0D1EE8A14.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F626C17D-137C-E011-BCF2-0015178C4CF8.root',
        # '/store/data/Run2011A/Jet/RECO/May10ReReco-v1/0000/F4F9C7E3-F37B-E011-8D20-00A0D1EEE654.root'
        'file:/uscms_data/d1/apresyan/pickevents_HT_PromptReco_MR400_R04_merged004.root'
        ),
    
    duplicateCheckMode = cms.untracked.string("noDuplicateCheck"),
    
    secondaryFileNames = cms.untracked.vstring())

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
process.MessageLogger.cerr.default.limit = 100

# summary
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# jet corrections
process.load('JetMETCorrections.Configuration.DefaultJEC_cff')

##################################################
############ SELECTIONS START HERE ###############
##################################################

# Primary vertex filter
process.primaryVertexFilter = cms.EDFilter("GoodVertexFilter",
                                           vertexCollection = cms.InputTag('offlinePrimaryVertices'),
                                           minimumNDOF = cms.uint32(4) ,
                                           maxAbsZ = cms.double(24),	
                                           maxd0 = cms.double(2)	
                                           )

# Scraping filter
process.scrapingVeto = cms.EDFilter("FilterOutScraping",
                                    applyfilter = cms.untracked.bool(True),
                                    debugOn = cms.untracked.bool(False),
                                    numtrack = cms.untracked.uint32(10),
                                    thresh = cms.untracked.double(0.25)
                                    )

# Instead of rejecting the event, add a flag indicating the HBHE noise 
process.load('CommonTools/RecoAlgos/HBHENoiseFilterResultProducer_cfi')
process.hbheflag = cms.Path(process.HBHENoiseFilterResultProducer)

# HBHE isolation + R45 filter
process.load('CommonTools/RecoAlgos/HBHENoiseFilter_cfi')

##################################################
############## SELECTIONS END HERE ###############
##################################################

process.promptanaTree = cms.EDAnalyzer("PromptAnaTree",
                                       outputCommands = cms.untracked.vstring(
                                           'drop *',
                                           'keep *_promptanaevent_*_*',
                                           'keep *_promptanamet_*_*',
                                           'keep *_promptanametdefault_*_*',
                                           'keep *_promptanatcmet_*_*',
                                           'keep *_promptanapfmet_*_*',
                                           'keep *_promptanacalotowers_*_*',
                                           'keep *_promptanaak5calojet_*_*',
                                           'keep *_promptanatrigger_*_*',
                                           'keep *_promptanavtx_*_*',
                                           'keep *_promptanatrack_*_*',
                                           'keep *_promptanaecalspikes_*_*',
                                           'keep *_HBHENoiseFilterResultProducer_*_*'
                                           ))

# Filter sequence
process.filterSequence = cms.Sequence(
    process.HBHENoiseFilter*
    process.CSCTightHaloFilter*
    process.primaryVertexFilter*
    process.scrapingVeto
)

from Configuration.DataProcessing.RecoTLR import customisePPData 
#call to customisation function customisePPData imported from Configuration.DataProcessing.RecoTLR
process = customisePPData(process)

process.rereco_step = cms.Path(process.filterSequence*process.caloTowersRec*(process.recoJets*process.recoJetIds+process.recoTrackJets)*process.recoJetAssociations*process.met) # re-reco jets and met

process.ntuple_step=cms.Path(
    process.filterSequence*
    (
        process.promptanaevent +
        process.promptanamet   +
        process.promptanametdefault   +
        process.promptanatcmet   +
        process.promptanapfmet   +
        process.promptanacalotowers +
        process.promptanaak5calojet +
        process.promptanatrigger +
        process.promptanavtx +
        process.promptanatrack +
        process.promptanaecalspikes
        # process.promptanacleanup
        )
    *process.promptanaTree
    )

process.schedule = cms.Schedule(process.rereco_step, process.hbheflag, process.ntuple_step)
#process.schedule = cms.Schedule(process.hbheflag, process.ntuple_step)
