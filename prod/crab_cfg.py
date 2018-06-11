'''
Effectively a multiple crab job submitter.
Submit with python crab_cfg.py
It creates a cfg and a working directory for each dataset you want to run on.
'''
from CRABClient.UserUtilities import config
from collections import OrderedDict

config = config()

config.General.transferOutputs = True
config.General.transferLogs    = True

config.JobType.psetName        = 'addElectronTransientTrack_cfg.py'
config.JobType.pluginName      = 'Analysis'
config.JobType.outputFiles     = ['output.root']
config.JobType.maxMemoryMB     = 2500
config.JobType.priority        = 99999

config.Data.unitsPerJob        = 100000
config.Data.splitting          = 'EventAwareLumiBased'

# JSON files:
# /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/
config.Data.publication        = True
config.Data.outputDatasetTag   = 'skimParkingBPHToEE'

config.Site.storageSite        = 'T2_CH_CERN'
# config.Site.blacklist          = ['T1_US_FNAL']
# config.Site.whitelist          = ['T2_CH_CERN']

if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    tag = 'skimParkingBPHToEEV2'

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea   = 'crab_bph_parking_' + tag
    config.Data.outLFNDirBase = '/store/group/phys_tau/BKstLL/' + tag 
    
    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    datasets = OrderedDict()

    datasets['ParkingBPH1'] = '/ParkingBPH1/Run2018A-14May2018-v1/MINIAOD'
    datasets['ParkingBPH2'] = '/ParkingBPH2/Run2018A-14May2018-v1/MINIAOD'
    datasets['ParkingBPH4'] = '/ParkingBPH4/Run2018A-14May2018-v1/MINIAOD'
    
    # datasets['ParkingBPH3'] = '/ParkingBPH3/Run2018A-14May2018-v1/MINIAOD'
    # datasets['ParkingBPH5'] = '/ParkingBPH5/Run2018A-14May2018-v1/MINIAOD'
    # datasets['ParkingBPH6'] = '/ParkingBPH6/Run2018A-14May2018-v1/MINIAOD'

    for k, v in datasets.iteritems():
        config.General.requestName = k
        config.Data.inputDataset   = v
        print 'submitting config:'
        print config
        submit(config)
