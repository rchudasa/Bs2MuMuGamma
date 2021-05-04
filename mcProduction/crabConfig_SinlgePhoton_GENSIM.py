from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'SinglePhotonFlatPt1To20_GENSIM_Run2018'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'PrivateMC'
config.JobType.psetName = 'GENToSIM_2018_cfg.py'
config.Data.outputPrimaryDataset = 'SinglePhotonFlatPt1To20_GENSIM_Run2018'

config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.splitting = 'EventBased'
config.Data.unitsPerJob = 5000
NJOBS = 250
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_cep/ntuples'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
