
import FWCore.ParameterSet.Config as cms

promptanamet = cms.EDProducer( "PromptAna_MET"
#                            , InputTag  = cms.InputTag('met::USER')
                            , InputTag  = cms.InputTag('caloMet')
                            , Prefix    = cms.string('calomet')
                            , Suffix    = cms.string('')
                            )

#promptanametdefault = cms.EDProducer("PromptAna_MET",
#                           InputTag  = cms.InputTag('met::RECO'),
#                           Prefix = cms.string('calometDflt'),
#                           Suffix = cms.string('')
#                           )
