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
        
# WTau3Mu analysers
from CMGTools.BKstLL.analyzers.BKJPsiEEAnalyzer          import BKJPsiEEAnalyzer    
from CMGTools.BKstLL.analyzers.BKJPsiTreeProducer        import BKJPsiTreeProducer

# import samples, signal
from CMGTools.BKstLL.samples.bph_parked_data_2018        import bph_parking_2018A, BPHParking6_AOD_2018A, bph_parking_2018A_skimmed, BPHParking3_2018A_skimmed

###################################################
###                   OPTIONS                   ###
###################################################
# Get all heppy options; set via "-o production" or "-o production=True"
# production = True run on batch, production = False (or unset) run locally
production         = getHeppyOption('production' , True )
pick_events        = getHeppyOption('pick_events', False)
###################################################
###               HANDLE SAMPLES                ###
###################################################
# samples = [BPHParking6_AOD_2018A]
samples = bph_parking_2018A_skimmed

for sample in samples:
    sample.triggers  = ['HLT_DoubleMu3_Trk_Tau3mu_v%d' %i for i in range(3, 12)]
    sample.splitFactor = splitFactor(sample, 2e5)

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
    comp                 = BPHParking3_2018A_skimmed
    selectedComponents   = [comp]
    comp.splitFactor     = 1
    comp.files           = comp.files[:1]
#     comp.files           = ['file:/afs/cern.ch/work/m/manzoni/hbb/CMSSW_10_1_4/src/CMGTools/BKstLL/prod/output_v2.root']
#     comp.fineSplitFactor = 4

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