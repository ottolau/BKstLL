# import dill # needed in order to serialise lambda functions, need to be installed by the user. See http://stackoverflow.com/questions/25348532/can-python-pickle-lambda-functions
from collections import OrderedDict

import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.framework.config     import printComps
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
from CMGTools.RootTools.utils.splitFactor import splitFactor

# import all analysers:
# Heppy analyzers
from PhysicsTools.Heppy.analyzers.core.JSONAnalyzer         import JSONAnalyzer
from PhysicsTools.Heppy.analyzers.core.SkimAnalyzerCount    import SkimAnalyzerCount
from PhysicsTools.Heppy.analyzers.core.EventSelector        import EventSelector
from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer    import VertexAnalyzer
from PhysicsTools.Heppy.analyzers.core.PileUpAnalyzer       import PileUpAnalyzer
from PhysicsTools.Heppy.analyzers.gen.GeneratorAnalyzer     import GeneratorAnalyzer
from PhysicsTools.Heppy.analyzers.gen.LHEWeightAnalyzer     import LHEWeightAnalyzer
        
# Tau-tau analysers        
from CMGTools.H2TauTau.proto.analyzers.JetAnalyzer          import JetAnalyzer

# WTau3Mu analysers
from CMGTools.BKstLL.analyzers.BsJPsiPhiGenAnalyzer            import BsJPsiPhiGenAnalyzer    
from CMGTools.BKstLL.analyzers.BsJPsiPhiGenTreeProducer        import BsJPsiPhiGenTreeProducer

import numpy as np

# import samples, signal
#from CMGTools.BKstLL.samples.signal import BdKstMM

puFileMC   = '$CMSSW_BASE/src/CMGTools/H2TauTau/data/MC_Moriond17_PU25ns_V1.root'
#puFileData = '/afs/cern.ch/user/a/anehrkor/public/Data_Pileup_2016_271036-284044_80bins.root'
puFileData = 'Data_Pileup_2016_271036-284044_80bins.root'

###################################################
###                   OPTIONS                   ###
###################################################
# Get all heppy options; set via "-o production" or "-o production=True"
# production = True run on batch, production = False (or unset) run locally
production         = getHeppyOption('production' , False)
pick_events        = getHeppyOption('pick_events', False)
###################################################
###               HANDLE SAMPLES                ###
###################################################

### Dataset: /BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17MiniAOD-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/MINIAODSIM

inputdata = np.loadtxt('BsToJpsiPhi_MINIAODSIM_filename.dat', dtype='str')
inputdata = ['root://cms-xrd-global.cern.ch/'+st for st in inputdata]

BsJPsiMMPhi = cfg.MCComponent(
    'BsJPsiMMPhi',
    #files = '/store/mc/RunIISummer16MiniAODv2/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/288EFE7A-60B9-E611-B1DB-0CC47A1DF7F2.root', 
    #files = ['root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6E47188D-3BAD-E711-8B3F-0CC47A6C1874.root',
#            ],
#    files = ['root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/6E47188D-3BAD-E711-8B3F-0CC47A6C1874.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A60BFA1C-3CAD-E711-AE31-008CFAC91B60.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/3A0A3362-36AC-E711-BBED-FA163ED8E79E.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/7A22593E-4AAC-E711-84F7-02163E0176A5.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/B8DB27C5-57AC-E711-9112-FA163E7B5756.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A0126530-5EAC-E711-98A1-02163E0164E3.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/8E5FA699-64AC-E711-9772-FA163EE771A0.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/5AC2A21E-64AC-E711-B29B-FA163E7625E2.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/949F3A31-6BAC-E711-B38C-FA163EF027B6.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/80125A7C-7EAC-E711-8CDA-FA163ECC515A.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17MiniAOD/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/4469C355-B2AC-E711-BE8E-02163E012D69.root',
#        ],
    files = inputdata

    # a list of local or xrootd files can be specified by hand.
    )

multi_thread = True

samples = [BsJPsiMMPhi]

for sample in samples:
    sample.triggers  = ['HLT_DoubleMu3_Trk_Tau3mu_v%d' %i for i in range(3, 12)]
    sample.splitFactor = splitFactor(sample, 1e5)
    sample.puFileData = puFileData
    sample.puFileMC   = puFileMC

selectedComponents = samples

###################################################
###                  ANALYSERS                  ###
###################################################
eventSelector = cfg.Analyzer(
    EventSelector,
    name='EventSelector',
    toSelect=[4148011548]
)

lheWeightAna = cfg.Analyzer(
    LHEWeightAnalyzer, name="LHEWeightAnalyzer",
    useLumiInfo=False
)

jsonAna = cfg.Analyzer(
    JSONAnalyzer,
    name='JSONAnalyzer',
)

skimAna = cfg.Analyzer(
    SkimAnalyzerCount,
    name='SkimAnalyzerCount'
)

vertexAna = cfg.Analyzer(
    VertexAnalyzer,
    name='VertexAnalyzer',
    fixedWeight=1,
    keepFailingEvents=True,
    verbose=False
)

pileUpAna = cfg.Analyzer(
    PileUpAnalyzer,
    name='PileUpAnalyzer',
    true=True
)

genAna = GeneratorAnalyzer.defaultConfig
genAna.allGenTaus = True # save in event.gentaus *ALL* taus, regardless whether hadronic / leptonic decay

# see SM HTT TWiki
# https://twiki.cern.ch/twiki/bin/viewauth/CMS/SMTauTau2016#Jet_Energy_Corrections
jetAna = cfg.Analyzer(
    JetAnalyzer,
    name              = 'JetAnalyzer',
    jetCol            = 'slimmedJets',
    jetPt             = 20.,
    jetEta            = 4.7,
    relaxJetId        = False, # relax = do not apply jet ID
    relaxPuJetId      = True, # relax = do not apply pileup jet ID
    jerCorr           = False,
    puJetIDDisc       = 'pileupJetId:fullDiscriminant',
    recalibrateJets   = True,
    applyL2L3Residual = 'MC',
    mcGT              = '80X_mcRun2_asymptotic_2016_TrancheIV_v8',
    dataGT            = '80X_dataRun2_2016SeptRepro_v7',
    #jesCorr = 1., # Shift jet energy scale in terms of uncertainties (1 = +1 sigma)
)

mainAna = cfg.Analyzer(
    BsJPsiPhiGenAnalyzer,
    name = 'BsJPsiPhiGenAnalyzer',
    flavour = 13,
)

treeProducer = cfg.Analyzer(
    BsJPsiPhiGenTreeProducer,
    name = 'BsJPsiPhiGenTreeProducer',
)

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = cfg.Sequence([
#     eventSelector,
    jsonAna,
    skimAna,
    genAna,
    vertexAna,
    pileUpAna,
    mainAna,
#     jetAna,
    treeProducer,
])

###################################################
###            SET BATCH OR LOCAL               ###
###################################################
if not production:
    comp                 = BsJPsiMMPhi
    selectedComponents   = [comp]
    if multi_thread:
        comp.splitFactor = 8
    else:
        comp.splitFactor     = 1
    comp.fineSplitFactor = 1
    #comp.nEvents         = 100
    #comp.files           = comp.files[:3]
#    comp.files           = [
#         'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/001CD385-A2B9-E611-8F02-FA163EC1154A.root',
#     ]

#     selectedComponents   = [QCD_Pt_15to30, QCD_Pt_80to120, QCD_Pt_3200toInf]
#     for comp in selectedComponents:
#         comp.splitFactor     = 1
#         comp.fineSplitFactor = 1
#         comp.files           = comp.files[:3]

#     comp.files = [
#        'file:/afs/cern.ch/work/m/manzoni/diTau2015/CMSSW_9_2_2_minimal_recipe/src/RecoMET/METPUSubtraction/test/output.root',
#        'root://xrootd.unl.edu//store/data/Run2016B/SingleMuon/MINIAOD/PromptReco-v1/000/272/760/00000/68B88794-7015-E611-8A92-02163E01366C.root',
#     ]

preprocessor = None

# the following is declared in case this cfg is used in input to the
# heppy.py script
from PhysicsTools.HeppyCore.framework.eventsfwlite import Events
config = cfg.Config(
    components   = selectedComponents,
    sequence     = sequence,
    services     = [],
    preprocessor = preprocessor,
    events_class = Events
)

printComps(config.components, True)
