import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

QCD_Pt_15to30     = creator.makeMCComponent("QCD_Pt_15to30"    , "/QCD_Pt_15to30_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"    , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_30to50     = creator.makeMCComponent("QCD_Pt_30to50"    , "/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"    , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_50to80     = creator.makeMCComponent("QCD_Pt_50to80"    , "/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM"    , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_80to120    = creator.makeMCComponent("QCD_Pt_80to120"   , "/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v3/MINIAODSIM"   , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_120to170   = creator.makeMCComponent("QCD_Pt_120to170"  , "/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"  , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_170to300   = creator.makeMCComponent("QCD_Pt_170to300"  , "/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"  , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_300to470   = creator.makeMCComponent("QCD_Pt_300to470"  , "/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"  , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_470to600   = creator.makeMCComponent("QCD_Pt_470to600"  , "/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM"  , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_600to800   = creator.makeMCComponent("QCD_Pt_600to800"  , "/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM"  , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_800to1000  = creator.makeMCComponent("QCD_Pt_800to1000" , "/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM" , "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_1000to1400 = creator.makeMCComponent("QCD_Pt_1000to1400", "/QCD_Pt_1000to1400_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM", "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_1400to1800 = creator.makeMCComponent("QCD_Pt_1400to1800", "/QCD_Pt_1400to1800_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM", "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_1800to2400 = creator.makeMCComponent("QCD_Pt_1800to2400", "/QCD_Pt_1800to2400_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM", "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_2400to3200 = creator.makeMCComponent("QCD_Pt_2400to3200", "/QCD_Pt_2400to3200_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v4/MINIAODSIM", "CMS", ".*root", 1.0, useAAA=True,)
QCD_Pt_3200toInf  = creator.makeMCComponent("QCD_Pt_3200toInf" , "/QCD_Pt_3200toInf_TuneCUETP8M1_13TeV_pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v5/MINIAODSIM" , "CMS", ".*root", 1.0, useAAA=True,)

qcd_samples = [
    QCD_Pt_15to30    ,
    QCD_Pt_30to50    ,
    QCD_Pt_50to80    ,
    QCD_Pt_80to120   ,
    QCD_Pt_120to170  ,
    QCD_Pt_170to300  ,
    QCD_Pt_300to470  ,
    QCD_Pt_470to600  ,
    QCD_Pt_600to800  ,
    QCD_Pt_800to1000 ,
    QCD_Pt_1000to1400,
    QCD_Pt_1400to1800,
    QCD_Pt_1800to2400,
    QCD_Pt_2400to3200,
    QCD_Pt_3200toInf ,
]





