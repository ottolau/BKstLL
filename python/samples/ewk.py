import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

creator = ComponentCreator()

DYJetsM10to50 = creator.makeMCComponent("DYJetsM10to50"    , "/DYJetsToLL_M-10to50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1/MINIAODSIM", "CMS", ".*root", 18610.   , useAAA=True,)
DYJetsM50     = creator.makeMCComponent("DYJetsM50"        , "/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v7-v1/MINIAODSIM"     , "CMS", ".*root",  2008.*3 , useAAA=True,)
WJetsToLNu    = creator.makeMCComponent("WJetsToLNu"       , "/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer17MiniAOD-92X_upgrade2017_realistic_v10-v1/MINIAODSIM"         , "CMS", ".*root", 3*20508.9, useAAA=True,)

ewk_samples = [
    DYJetsM10to50,
    DYJetsM50    ,
    WJetsToLNu   ,
]





