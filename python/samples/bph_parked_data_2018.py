import PhysicsTools.HeppyCore.framework.config as cfg
import os

#####COMPONENT CREATOR

from CMGTools.RootTools.samples.ComponentCreator import ComponentCreator

json = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification//Collisions18/13TeV/DCSOnly/json_DCSONLY.txt'

creator = ComponentCreator()

BPHParking1_2018A = creator.makeDataComponent("BPHParking1_2018A", "/ParkingBPH1/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)
BPHParking2_2018A = creator.makeDataComponent("BPHParking2_2018A", "/ParkingBPH2/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)
BPHParking3_2018A = creator.makeDataComponent("BPHParking3_2018A", "/ParkingBPH3/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)
BPHParking4_2018A = creator.makeDataComponent("BPHParking4_2018A", "/ParkingBPH4/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)
BPHParking5_2018A = creator.makeDataComponent("BPHParking5_2018A", "/ParkingBPH5/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)
BPHParking6_2018A = creator.makeDataComponent("BPHParking6_2018A", "/ParkingBPH6/Run2018A-14May2018-v1/MINIAOD", "CMS", ".*root", json, useAAA=True)

bph_parking_2018A = [
    BPHParking1_2018A,
    BPHParking2_2018A,
    BPHParking3_2018A,
    BPHParking4_2018A,
    BPHParking5_2018A,
    BPHParking6_2018A,
]

BPHParking1_AOD_2018A = creator.makeDataComponent("BPHParking1_AOD_2018A", "/ParkingBPH1/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)
BPHParking2_AOD_2018A = creator.makeDataComponent("BPHParking2_AOD_2018A", "/ParkingBPH2/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)
BPHParking3_AOD_2018A = creator.makeDataComponent("BPHParking3_AOD_2018A", "/ParkingBPH3/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)
BPHParking4_AOD_2018A = creator.makeDataComponent("BPHParking4_AOD_2018A", "/ParkingBPH4/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)
BPHParking5_AOD_2018A = creator.makeDataComponent("BPHParking5_AOD_2018A", "/ParkingBPH5/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)
BPHParking6_AOD_2018A = creator.makeDataComponent("BPHParking6_AOD_2018A", "/ParkingBPH6/Run2018A-14May2018-v1/AOD", "CMS", ".*root", json, useAAA=True)

bph_parking_AOD_2018A = [
    BPHParking1_AOD_2018A,
    BPHParking2_AOD_2018A,
    BPHParking3_AOD_2018A,
    BPHParking4_AOD_2018A,
    BPHParking5_AOD_2018A,
    BPHParking6_AOD_2018A,
]



##########################################################################################
from glob import glob

BPHParking1_2018A_skimmed = BPHParking1_2018A
BPHParking2_2018A_skimmed = BPHParking2_2018A
BPHParking3_2018A_skimmed = BPHParking3_2018A
BPHParking4_2018A_skimmed = BPHParking4_2018A
BPHParking5_2018A_skimmed = BPHParking5_2018A
BPHParking6_2018A_skimmed = BPHParking6_2018A


BPHParking1_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV2/ParkingBPH1/skimParkingBPHToEE/180601_111629/0000/*.root')
BPHParking2_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV2/ParkingBPH2/skimParkingBPHToEE/180601_111818/0000/*.root')
BPHParking3_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV1/ParkingBPH3/skimParkingBPHToEE/180529_202857/0000/*.root')
BPHParking4_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV2/ParkingBPH4/skimParkingBPHToEE/180601_112005/0000/*.root')
BPHParking5_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV1/ParkingBPH5/skimParkingBPHToEE/180529_203042/0000/*.root')
BPHParking6_2018A_skimmed.files = glob('/eos/cms/store/group/phys_tau/BKstLL/skimParkingBPHToEEV1/ParkingBPH6/skimParkingBPHToEE/180529_203228/0000/*.root')

bph_parking_2018A_skimmed = [
    BPHParking1_2018A_skimmed,
    BPHParking2_2018A_skimmed,
    BPHParking3_2018A_skimmed,
    BPHParking4_2018A_skimmed,
    BPHParking5_2018A_skimmed,
    BPHParking6_2018A_skimmed,
]



