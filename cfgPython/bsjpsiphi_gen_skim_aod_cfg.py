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
from CMGTools.BKstLL.analyzers.BsJPsiPhiGenAnalyzer_Skim_AOD        import BsJPsiPhiGenAnalyzer_Skim_AOD
from CMGTools.BKstLL.analyzers.BsJPsiPhiGenTreeProducer_Skim_AOD        import BsJPsiPhiGenTreeProducer_Skim_AOD

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

### Dataset: /BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer17DRStdmix-NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/AODSIM

inputdata = np.loadtxt('BsToJpsiPhi_AODSIM_filename.dat', dtype='str')
inputdata = ['root://cms-xrd-global.cern.ch/'+st for st in inputdata]

BsJPsiMMPhi = cfg.MCComponent(
    'BsJPsiMMPhi',
    #files = '/store/mc/RunIISummer16MiniAODv2/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/288EFE7A-60B9-E611-B1DB-0CC47A1DF7F2.root', 
    files = ['file:46BFC0D9-D3AC-E711-850E-008CFAC913F8.root',],
    #files = ['file:46BFC0D9-D3AC-E711-850E-008CFAC913F8.root','file:A6221481-FDAB-E711-B8A0-02163E01765E.root',],

#    files = ['root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/46BFC0D9-D3AC-E711-850E-008CFAC913F8.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/80E3CBD8-F2AC-E711-9191-008CFAC93F0C.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/629DACDD-E0AB-E711-8792-02163E016491.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/8619CBC2-FBAB-E711-8754-02163E01643C.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/A6221481-FDAB-E711-B8A0-02163E01765E.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/569DA21D-02AC-E711-9F91-FA163E3A1163.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/44E127E3-02AC-E711-8EE2-02163E014FDA.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/C6401406-03AC-E711-9CEA-FA163E5DA5D0.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/DAE60915-EFAB-E711-B27C-FA163E74C5BB.root',
#        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer17DRStdmix/BsToJpsiPhi_BMuonFilter_TuneCUEP8M1_13TeV-pythia8-evtgen/AODSIM/NZSFlatPU28to62_92X_upgrade2017_realistic_v10-v1/150000/DC42DBA8-08AC-E711-A101-FA163E3111A9.root',
#    ]
    #files = inputdata

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
    BsJPsiPhiGenAnalyzer_Skim_AOD,
    name = 'BsJPsiPhiGenAnalyzer_Skim_AOD',
    flavour = 11,
)

treeProducer = cfg.Analyzer(
    BsJPsiPhiGenTreeProducer_Skim_AOD,
    name = 'BsJPsiPhiGenTreeProducer_Skim_AOD',
)

###################################################
###                  SEQUENCE                   ###
###################################################
sequence = cfg.Sequence([
#     eventSelector,
#    jsonAna,
#    skimAna,
#    genAna,
#    vertexAna,
#    pileUpAna,
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
        #comp.splitFactor = len(comp.files)
        comp.splitFactor = 8
    else:
        comp.splitFactor     = 1
    comp.fineSplitFactor = 1
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
