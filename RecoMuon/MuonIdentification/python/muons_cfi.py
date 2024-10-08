# This file name is temporary and ment for development only.
# The content of this file will be moved to muons_cfi as soon as the complete work flow is in place.

import FWCore.ParameterSet.Config as cms

muons = cms.EDProducer("MuonProducer",
                       ActivateDebug = cms.untracked.bool(False),
                       InputMuons = cms.InputTag("muons1stStep"),

                       FillPFMomentumAndAssociation = cms.bool(True),
                       PFCandidates = cms.InputTag("particleFlowTmp"),

                       FillTimingInfo = cms.bool(True),
                       
                       FillDetectorBasedIsolation = cms.bool(True),
                       EcalIsoDeposits  = cms.InputTag("muIsoDepositCalByAssociatorHits","ecal"),
                       HcalIsoDeposits  = cms.InputTag("muIsoDepositCalByAssociatorTowers","hcal"),
                       HoIsoDeposits    = cms.InputTag("muIsoDepositCalByAssociatorTowers","ho"),
                       TrackIsoDeposits = cms.InputTag("muIsoDepositTk"),
                       JetIsoDeposits   = cms.InputTag("muIsoDepositJets"),

                       FillPFIsolation = cms.bool(True),                     
                       PFIsolation = cms.PSet(
                        pfIsolationR03 = cms.PSet(chargedParticle            = cms.InputTag("muPFIsoValueChargedAll03"),
                                                  chargedHadron              = cms.InputTag("muPFIsoValueCharged03"),
                                                  neutralHadron              = cms.InputTag("muPFIsoValueNeutral03"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFIsoValueNeutralHighThreshold03"),
                                                  photon                     = cms.InputTag("muPFIsoValueGamma03"),
                                                  photonHighThreshold        = cms.InputTag("muPFIsoValueGammaHighThreshold03"),
                                                  pu                         = cms.InputTag("muPFIsoValuePU03")), 
                        pfIsolationR04 = cms.PSet(chargedParticle            = cms.InputTag("muPFIsoValueChargedAll04"),
                                                  chargedHadron              = cms.InputTag("muPFIsoValueCharged04"),
                                                  neutralHadron              = cms.InputTag("muPFIsoValueNeutral04"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFIsoValueNeutralHighThreshold04"),
                                                  photon                     = cms.InputTag("muPFIsoValueGamma04"),
                                                  photonHighThreshold        = cms.InputTag("muPFIsoValueGammaHighThreshold04"),
                                                  pu                         = cms.InputTag("muPFIsoValuePU04")),
                        pfIsoMeanDRProfileR03 = cms.PSet(chargedParticle     = cms.InputTag("muPFMeanDRIsoValueChargedAll03"),
                                                  chargedHadron              = cms.InputTag("muPFMeanDRIsoValueCharged03"),
                                                  neutralHadron              = cms.InputTag("muPFMeanDRIsoValueNeutral03"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFMeanDRIsoValueNeutralHighThreshold03"),
                                                  photon                     = cms.InputTag("muPFMeanDRIsoValueGamma03"),
                                                  photonHighThreshold        = cms.InputTag("muPFMeanDRIsoValueGammaHighThreshold03"),
                                                  pu                         = cms.InputTag("muPFMeanDRIsoValuePU03")), 
                        pfIsoMeanDRProfileR04 = cms.PSet(chargedParticle     = cms.InputTag("muPFMeanDRIsoValueChargedAll04"),
                                                  chargedHadron              = cms.InputTag("muPFMeanDRIsoValueCharged04"),
                                                  neutralHadron              = cms.InputTag("muPFMeanDRIsoValueNeutral04"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFMeanDRIsoValueNeutralHighThreshold04"),
                                                  photon                     = cms.InputTag("muPFMeanDRIsoValueGamma04"),
                                                  photonHighThreshold        = cms.InputTag("muPFMeanDRIsoValueGammaHighThreshold04"),
                                                  pu                         = cms.InputTag("muPFMeanDRIsoValuePU04")),
                       pfIsoSumDRProfileR03 = cms.PSet(chargedParticle       = cms.InputTag("muPFSumDRIsoValueChargedAll03"),
                                                  chargedHadron              = cms.InputTag("muPFSumDRIsoValueCharged03"),
                                                  neutralHadron              = cms.InputTag("muPFSumDRIsoValueNeutral03"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFSumDRIsoValueNeutralHighThreshold03"),
                                                  photon                     = cms.InputTag("muPFSumDRIsoValueGamma03"),
                                                  photonHighThreshold        = cms.InputTag("muPFSumDRIsoValueGammaHighThreshold03"),
                                                  pu                         = cms.InputTag("muPFSumDRIsoValuePU03")), 
                        pfIsoSumDRProfileR04 = cms.PSet(chargedParticle      = cms.InputTag("muPFSumDRIsoValueChargedAll04"),
                                                  chargedHadron              = cms.InputTag("muPFSumDRIsoValueCharged04"),
                                                  neutralHadron              = cms.InputTag("muPFSumDRIsoValueNeutral04"),
                                                  neutralHadronHighThreshold = cms.InputTag("muPFSumDRIsoValueNeutralHighThreshold04"),
                                                  photon                     = cms.InputTag("muPFSumDRIsoValueGamma04"),
                                                  photonHighThreshold        = cms.InputTag("muPFSumDRIsoValueGammaHighThreshold04"),
                                                  pu                         = cms.InputTag("muPFSumDRIsoValuePU04"))
                       ),

                       FillSelectorMaps = cms.bool(True),
                       SelectorMaps = cms.VInputTag(cms.InputTag("muidTMLastStationOptimizedLowPtLoose"),
                                                    cms.InputTag("muidTMLastStationOptimizedLowPtTight"),
                                                    cms.InputTag("muidTM2DCompatibilityLoose"),
                                                    cms.InputTag("muidTM2DCompatibilityTight"),
                                                    cms.InputTag("muidTrackerMuonArbitrated"),
                                                    cms.InputTag("muidTMLastStationAngLoose"),
                                                    cms.InputTag("muidGlobalMuonPromptTight"),
                                                    cms.InputTag("muidGMStaChiCompatibility"),
                                                    cms.InputTag("muidTMLastStationAngTight"),
                                                    cms.InputTag("muidGMTkChiCompatibility"),
                                                    cms.InputTag("muidTMOneStationAngTight"),
                                                    cms.InputTag("muidTMOneStationAngLoose"),
                                                    cms.InputTag("muidTMLastStationLoose"),
                                                    cms.InputTag("muidTMLastStationTight"),
                                                    cms.InputTag("muidTMOneStationTight"),
                                                    cms.InputTag("muidTMOneStationLoose"),
                                                    cms.InputTag("muidAllArbitrated"),
                                                    cms.InputTag("muidGMTkKinkTight"),
                                                    cms.InputTag("muidRPCMuLoose")
                                                    ),

                       FillShoweringInfo = cms.bool(True),
                       ShowerInfoMap = cms.InputTag("muonShowerInformation"),

                       FillCosmicsIdMap = cms.bool(True),
                       CosmicIdMap = cms.InputTag("cosmicsVeto"),

                       ComputeStandardSelectors = cms.bool(True),
                       vertices = cms.InputTag("offlinePrimaryVertices")
                       
                       )

# not commisoned and not relevant in FastSim (?):
from Configuration.Eras.Modifier_fastSim_cff import fastSim
fastSim.toModify(muons, FillCosmicsIdMap = False, FillSelectorMaps = False)
