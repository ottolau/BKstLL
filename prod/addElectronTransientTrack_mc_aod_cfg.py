# Import CMS python class definitions such as Process, Source, and EDProducer
import FWCore.ParameterSet.Config as cms
from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
import numpy as np

process = cms.Process('TTK')
patAlgosToolsTask = getPatAlgosToolsTask(process)

#process.load( "PhysicsTools.PatAlgos.triggerLayer1.triggerProducer_cff" )
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.load('PhysicsTools.PatAlgos.producersLayer1.patCandidates_cff')
patAlgosToolsTask.add(process.patCandidatesTask)
process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
patAlgosToolsTask.add(process.selectedPatCandidatesTask)
process.load('PhysicsTools.PatAlgos.cleaningLayer1.cleanPatCandidates_cff')
patAlgosToolsTask.add(process.cleanPatCandidatesTask)

### add trigger information to the configuration
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, None, None, None, None, '' )

from PhysicsTools.PatAlgos.tools.coreTools import *
runOnData( process,  names=['Photons', 'Electrons','Muons','Taus','Jets'], outputModules = [] )
removeMCMatching(process, names=['All'], outputModules=[])

#process.TFileService = cms.Service("TFileService", fileName = cms.string('ttk_debug.root'))

############ THOMAS' MC ##############
#process.GlobalTag.globaltag = cms.string( "94X_mc2017_realistic_v12" )  
process.GlobalTag.globaltag = cms.string( "102X_upgrade2018_realistic_v15" )  

inputdata = np.loadtxt('BdKstEE_AODSIM_filename.dat', dtype='str')
inputdata = ['root://cms-xrd-global.cern.ch/'+st for st in inputdata]
readFiles = cms.untracked.vstring()
readFiles.extend(inputdata)

# Configure the object that reads the input file
process.source = cms.Source('PoolSource', 
#     fileNames = readFiles,
    fileNames = cms.untracked.vstring(
#         'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18DR/BdToKstarJpsi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15_ext1-v1/00000/9C4D3D80-F714-1F4E-99E1-C38EC0E3B096.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_9.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_8.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_7.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_6.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_5.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_4.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_3.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_2.root',
#        'file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/Bu_KJPsi_ee_MINIAODSIM_1.root',

#         'file:/afs/cern.ch/work/m/manzoni/hbb/CMSSW_9_4_6_patch1/src/CMGTools/BKstLL/prod/Bu_KJPsi_ee_MINIAODSIM_9.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_9.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_4.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_7.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_3.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_1.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_6.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_5.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_8.root',
#         'root://cms-xrd-global.cern.ch//store/user/tstreble/Bu_KJPsi_ee_Pythia/BuToKJPsiee_Pythia_MINIAODSIM_18_06_05/180605_092537/0000/Bu_KJPsi_ee_MINIAODSIM_2.root',
          '/store/mc/RunIIAutumn18DR/BsToPhi_Toee_MuFilter_SoftQCDnonD_TuneCP5_13TeV-pythia8-evtgen/AODSIM/PUPoissonAve20_102X_upgrade2018_realistic_v15-v3/280000/92AAD7E0-3512-A640-9035-190B8EAEFF0B.root',
    ),
)

process.selectedElectrons = cms.EDFilter(
    "PATElectronSelector",
    src = cms.InputTag( 'slimmedElectrons' ),
    cut = cms.string("abs(eta)<2.5")
)

process.diEleCandProd = cms.EDProducer(
    "CandViewShallowCloneCombiner",
    decay = cms.string("selectedElectrons@+ selectedElectrons@-"),
    cut   = cms.string("mass < 6."),
)

process.countdiEleCand = cms.EDFilter(
    "CandViewCountFilter",
    src = cms.InputTag("diEleCandProd"),
    minNumber = cms.uint32(1)
)

process.ttk = cms.EDProducer(
    'AddElectronTransientTrack',
    patEleSrc = cms.InputTag('selectedPatElectrons'),
)

process.ttkPath = cms.Path(
    #process.selectedElectrons +
    #process.diEleCandProd +
    #process.countdiEleCand +
    process.ttk,
    patAlgosToolsTask
)

# Configure the object that writes an output file
process.out = cms.OutputModule('PoolOutputModule',
    fileName = cms.untracked.string('ttk_debug.root'),
    SelectEvents = cms.untracked.PSet( 
        SelectEvents = cms.vstring("ttkPath")
    )
)


process.prunedOutput = cms.EndPath( process.out )

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 1000 )
)

## logger
process.MessageLogger.cerr.FwkReport.reportEvery = 100
