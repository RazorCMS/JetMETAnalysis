import FWCore.ParameterSet.Config as cms

promptanaak5pfjet = cms.EDProducer("PromptAna_PFJet",
                            InputTag = cms.InputTag('ak4PFJets'),
                            JetCorrectionService = cms.string('ak4PFL2L3'),       
                            Prefix = cms.string('ak4PFJet'),
                            Suffix = cms.string('')
                            )

