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
from CMGTools.H2TauTau.proto.analyzers.TriggerAnalyzer      import TriggerAnalyzer
from CMGTools.H2TauTau.proto.analyzers.JetAnalyzer          import JetAnalyzer

# WTau3Mu analysers
from CMGTools.BKstLL.analyzers.L1RateAnalyzer               import L1RateAnalyzer    
from CMGTools.BKstLL.analyzers.L1RateTreeProducer           import L1RateTreeProducer

# import samples, signal
from CMGTools.BKstLL.samples.hltphysics import EphemeralHLTPhysics_2017F

puFileMC   = '$CMSSW_BASE/src/CMGTools/H2TauTau/data/MC_Moriond17_PU25ns_V1.root'
puFileData = '/afs/cern.ch/user/a/anehrkor/public/Data_Pileup_2016_271036-284044_80bins.root'

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
samples = [EphemeralHLTPhysics_2017F]

for sample in samples:
    sample.triggers  = ['HLT_Physics_part%d_v7' %i for i in range(0, 8)]
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
    toSelect=[]
)

triggerAna = cfg.Analyzer(
    TriggerAnalyzer,
    name='TriggerAnalyzer',
    addTriggerObjects=True,
    requireTrigger=True,
    usePrescaled=True,
    triggerResultsHandle=('TriggerResults', '', 'HLT'),
    triggerObjectsHandle=('slimmedPatTrigger', '', 'RECO'),
)

jsonAna = cfg.Analyzer(
    JSONAnalyzer,
    name='JSONAnalyzer',
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

mainAna = cfg.Analyzer(
    L1RateAnalyzer,
    name = 'L1RateAnalyzer',
    onlyBX0 = True, # BE CAREFUL!
)

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
    recalibrateJets   = False, # don't recalibrate if you take the latest & greatest samples
    applyL2L3Residual = 'MC',
    mcGT              = '80X_mcRun2_asymptotic_2016_TrancheIV_v8',
    dataGT            = '80X_dataRun2_2016SeptRepro_v7',
    selectedLeptons   = [],
    #jesCorr = 1., # Shift jet energy scale in terms of uncertainties (1 = +1 sigma)
)

treeProducer = cfg.Analyzer(
    L1RateTreeProducer,
    name = 'L1RateTreeProducer',
)

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = cfg.Sequence([
#     eventSelector,
    jsonAna,
    triggerAna,
    vertexAna,
    pileUpAna,
    jetAna,
    mainAna,
    treeProducer,
])

###################################################
###            SET BATCH OR LOCAL               ###
###################################################
if not production:
    comp                 = EphemeralHLTPhysics_2017F
    selectedComponents   = [comp]
    comp.splitFactor     = 1
    comp.fineSplitFactor = 1
    comp.files           = comp.files[:10]

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