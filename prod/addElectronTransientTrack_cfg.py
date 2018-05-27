# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms

process = cms.Process('TTK')

process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.GlobalTag.globaltag = cms.string( "101X_dataRun2_Prompt_v9" )

# Configure the object that reads the input file
process.source = cms.Source('PoolSource', 
    fileNames = cms.untracked.vstring(
        'root://cms-xrd-global.cern.ch//store/data/Run2018A/ParkingBPH1/MINIAOD/14May2018-v1/30000/C870BB22-1D5A-E811-AE46-0025905D1D7A.root',
    ),
)

process.ttk = cms.EDProducer(
    'AddElectronTransientTrack',
    patEleSrc = cms.InputTag('slimmedElectrons'),
)

# Configure the object that writes an output file
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('output.root'),
)

process.ttkPath = cms.Path(process.ttk)

process.prunedOutput = cms.EndPath( process.out )

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

## logger
process.MessageLogger.cerr.FwkReport.reportEvery = 1