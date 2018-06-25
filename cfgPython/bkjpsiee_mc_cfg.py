from collections import OrderedDict

import PhysicsTools.HeppyCore.framework.config as cfg
from PhysicsTools.HeppyCore.framework.config     import printComps
from PhysicsTools.HeppyCore.framework.heppy_loop import getHeppyOption
from PhysicsTools.Heppy.utils.cmsswPreprocessor import CmsswPreprocessor
from CMGTools.RootTools.utils.splitFactor import splitFactor

# import all analysers:
# Heppy analyzers
from PhysicsTools.Heppy.analyzers.core.JSONAnalyzer      import JSONAnalyzer
from PhysicsTools.Heppy.analyzers.core.SkimAnalyzerCount import SkimAnalyzerCount
from PhysicsTools.Heppy.analyzers.core.EventSelector     import EventSelector
from PhysicsTools.Heppy.analyzers.objects.VertexAnalyzer import VertexAnalyzer
from PhysicsTools.Heppy.analyzers.core.PileUpAnalyzer    import PileUpAnalyzer
from PhysicsTools.Heppy.analyzers.gen.GeneratorAnalyzer  import GeneratorAnalyzer
from PhysicsTools.Heppy.analyzers.gen.LHEWeightAnalyzer  import LHEWeightAnalyzer
        
from CMGTools.H2TauTau.proto.analyzers.FileCleaner       import FileCleaner

# WTau3Mu analysers
from CMGTools.BKstLL.analyzers.BKJPsiEEAnalyzer          import BKJPsiEEAnalyzer    
from CMGTools.BKstLL.analyzers.BKJPsiTreeProducer        import BKJPsiTreeProducer

# import samples, signal
from CMGTools.BKstLL.samples.signal                      import BdKJPsiEEMuGenFilter

###################################################
###                   OPTIONS                   ###
###################################################
# Get all heppy options; set via "-o production" or "-o production=True"
# production = True run on batch, production = False (or unset) run locally
production         = getHeppyOption('production' , False)
pick_events        = getHeppyOption('pick_events', False)
run_skim           = getHeppyOption('run_skim'   , False)
###################################################
###               HANDLE SAMPLES                ###
###################################################
# samples = [BPHParking6_AOD_2018A]

# FIXME!!
BdKJPsiEEMuGenFilter.dataset_entries = 19221
samples = [BdKJPsiEEMuGenFilter]

for sample in samples:
    sample.triggers  = ['HLT_Mu9_IP6_v%d' %i for i in range(3, 12)]
    sample.splitFactor = splitFactor(sample, 1e3)

selectedComponents = samples

###################################################
###                  ANALYSERS                  ###
###################################################
eventSelector = cfg.Analyzer(
    EventSelector,
    name='EventSelector',
    toSelect=[4148011548],
)

jsonAna = cfg.Analyzer(
    JSONAnalyzer,
    name='JSONAnalyzer',
)

skimAna = cfg.Analyzer(
    SkimAnalyzerCount,
    name='SkimAnalyzerCount',
)

vertexAna = cfg.Analyzer(
    VertexAnalyzer,
    name='VertexAnalyzer',
    fixedWeight=1,
    keepFailingEvents=True,
    verbose=False,
)

pileUpAna = cfg.Analyzer(
    PileUpAnalyzer,
    name='PileUpAnalyzer',
    true=True,
)

mainAna = cfg.Analyzer(
    BKJPsiEEAnalyzer,
    name = 'BKJPsiEEAnalyzer',
)

treeProducer = cfg.Analyzer(
    BKJPsiTreeProducer,
    name = 'BKJPsiTreeProducer',
    addTagMu = True,
)

fileCleaner = cfg.Analyzer(
    FileCleaner,
    name='FileCleaner'
)

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = cfg.Sequence([
#     eventSelector,
    jsonAna,
    skimAna,
    vertexAna,
    pileUpAna,
    mainAna,
    treeProducer,
])

###################################################
###            SET BATCH OR LOCAL               ###
###################################################
if not production:
    comp                 = BdKJPsiEEMuGenFilter
    selectedComponents   = [comp]
    comp.splitFactor     = 1
    comp.files           = ['file:/eos/user/m/manzoni/BKstLL/MC/BKPsiEEMuFilter/output_v2.root']
#     comp.files           = comp.files[:1]
#     comp.files           = ['file:/afs/cern.ch/work/m/manzoni/hbb/CMSSW_10_1_4/src/CMGTools/BKstLL/prod/output_v2.root']
#     comp.fineSplitFactor = 4


## PREPROCESSOR
if run_skim:
    fname = '$CMSSW_BASE/src/CMGTools/BKstLL/prod/addElectronTransientTrack_mc_cfg.py'
    sequence.append(fileCleaner)
    preprocessor = CmsswPreprocessor(fname, addOrigAsSecondary=False)
else:
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