import FWCore.ParameterSet.Config as cms

process = cms.Process("PickEvent")

process.source = cms.Source("PoolSource",
    eventsToProcess = cms.untracked.VEventRange("5511205:16761"),
    fileNames = cms.untracked.vstring()
)
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.Out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('pickevents.root')
)


process.end = cms.EndPath(process.Out)


