from CRABClient.UserUtilities import config
config = config()

#config.section_('General')
config.General.requestName = 'SinglePhotonFlatPt1To20_DIGIRAW_Run2018'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

#config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'DIGIToRaw_noPU2018_cfg.py'
#config.Data.outputPrimaryDataset = 'SinglePhotonFlatPt1To20_DIGIRAW_Run2018'
config.Data.inputDBS = 'phys03'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.numCores = 8
config.Data.inputDataset ='/SinglePhotonFlatPt1To20_GENSIM_Run2018/rchudasa-crab_SinglePhotonFlatPt1To20_GENSIM_Run2018-ebb6bae7693268003fa4895b6ae4f74c/USER'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1 

#config.Data.outLFNDirBase = '/store/group/phys_diffraction/lbyl_2018/mc_cep/ntuples'
config.Data.outLFNDirBase = '/store/group/phys_heavyions/rchudasa/lowPT_photonReco'
config.Data.publication = True 
config.Site.storageSite = 'T2_CH_CERN'
