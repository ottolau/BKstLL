import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

#BdKstMM = creator.makeMCComponent(
BdKstMM = creator.makeMCComponentFromEOS(
    'BdKstMM', 
    '/BdToKstarMuMu_BMuonFilter_SoftQCDnonD_TuneCUEP8M1_13TeV-pythia8-evtgen/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM', 
    'XXX', 
    '.*root', 
    1.0, 
    #useAAA=True,
)


BdKstEE = creator.makeMCComponentFromLocal(
    'BdKstEE', 
    'XXX', 
    path = '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/', 
#     path = '/afs/cern.ch/user/s/slacapra/public/Kstee/',
    pattern='.*MiniAOD-00038.*root', 
    xSec=1.0, 
)

# BdKstEE.files += [
#     '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/BPH-RunIIFall17MiniAOD-00038_MINIAOD_firstbunch.root',
#     '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/BPH-RunIIFall17MiniAOD-00038_MINIAOD_fourthbunch.root',
#     '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/BPH-RunIIFall17MiniAOD-00038_MINIAOD_secondbunch.root',
#     '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/BPH-RunIIFall17MiniAOD-00038_MINIAOD_thirdbunch.root',
#     '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/BPH-RunIIFall17MiniAOD-00038_MINIAOD_fifthbunch.root',
# ]


BdKstEEOnlyMuGenFilter = creator.makeMCComponentFromLocal(
    'BdKstEEOnlyMuGenFilter', 
    'XXX', 
    path = '/eos/cms/store/group/phys_bphys/fiorendi/13TeV/MC/B0KStarEE/OnlyMUGENfilter/', 
    pattern='.*MiniAOD-00038.*root', 
    xSec=1.0, 
)

BdKJPsiEEMuGenFilter = creator.makeMyPrivateMCComponent(
    name        = 'BdKJPsiEEMuGenFilter', 
    dataset     = '/Bu_KJPsi_ee_Pythia/tstreble-BuToKJPsiee_Pythia_MINIAODSIM_18_06_05-393c5b9eb90ffa47a3c9a0f5562ec979/USER', 
    user        = 'PRIVATE',
    dbsInstance = 'phys03',
    pattern     = '.*root', 
    xSec        = 1.0,
    useAAA      = True, 
)
