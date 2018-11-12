import FWCore.ParameterSet.Config as cms

TTTracksFromTracklet = cms.EDProducer("L1TrackProducer",
                                      SimTrackSource = cms.InputTag("g4SimHits"),
                                      SimVertexSource = cms.InputTag("g4SimHits"),
                                      TTStubSource = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                      MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                      MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                      TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                      TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                      BeamSpotSource = cms.InputTag("offlineBeamSpot"),
                                      asciiFileName = cms.untracked.string(""),
                                      failscenario = cms.untracked.int32(0),
                                      #GenJetInputTag = cms.InputTag("ak4GenJets", ""),
                                      trackerGeometryType  = cms.untracked.string("")  #tilted barrel is assumed, use "flat" if running on flat
    )

TTTracksFromTrackletEmulation = cms.EDProducer("L1FPGATrackProducer",
                                               # general L1 tracking inputs
                                               SimTrackSource = cms.InputTag("g4SimHits"),
                                               SimVertexSource = cms.InputTag("g4SimHits"),
                                               TTStubSource = cms.InputTag("TTStubsFromPhase2TrackerDigis","StubAccepted"),
                                               MCTruthClusterInputTag = cms.InputTag("TTClusterAssociatorFromPixelDigis", "ClusterAccepted"),
                                               MCTruthStubInputTag = cms.InputTag("TTStubAssociatorFromPixelDigis", "StubAccepted"),
                                               TrackingParticleInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                               TrackingVertexInputTag = cms.InputTag("mix", "MergedTrackTruth"),
                                               BeamSpotSource = cms.InputTag("offlineBeamSpot"),
                                               asciiFileName = cms.untracked.string(""),
                                               failscenario = cms.untracked.int32(0),
                                               trackerGeometryType  = cms.untracked.string(""),  #tilted barrel is assumed, use "flat" if running on flat
                                               # specific emulation inputs 
                                               # (if running on CRAB use "../../fitpattern.txt" etc instead)
                                               fitPatternFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/fitpattern.txt'),
#                                               memoryModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/memorymodules_new.dat'),
#                                               processingModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/processingmodules_new.dat'),
#                                               wiresFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/wires_new.dat'),
                                               memoryModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/memorymodules_hourglass.dat'),
                                               processingModulesFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/processingmodules_hourglass.dat'),
                                               wiresFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/wires_hourglass.dat'),
                                               DTCLinkFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/calcNumDTCLinks.txt'),
                                               moduleCablingFile  = cms.FileInPath('L1Trigger/TrackFindingTracklet/test/modules_T5v3_27SP_nonant_tracklet.dat')
    )
