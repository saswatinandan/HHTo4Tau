from WMCore.Configuration import Configuration
config = Configuration()
config.section_('General')
config.General.transferOutputs = True
config.General.workArea = 'MC3rdJune'
config.section_('JobType')
config.JobType.psetName = 'run_mc_80X.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['ggtree_mc.root']
config.JobType.inputFiles = ['Summer16_23Sep2016AllV4_DATA.db','Summer16_23Sep2016V4_MC.db']
config.JobType.maxMemoryMB = 3000
config.section_('Data')
config.Data.unitsPerJob = 8
config.Data.splitting = 'FileBased'
config.Data.outLFNDirBase = '/store/user/snandan/Moriond17/MC/WnbinnedJets3rdJune'
config.Data.allowNonValidInputDataset = True
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'


if __name__ == '__main__':

    from CRABAPI.RawCommand import crabCommand
    from CRABClient.ClientExceptions import ClientException
    from httplib import HTTPException

    # We want to put all the CRAB project directories from the tasks we submit here into one common directory.
    # That's why we need to set this parameter (here or above in the configuration file, it does not matter, we will not overwrite it).
    config.General.workArea = 'MC3rdJune'

    def submit(config):
        try:
            crabCommand('submit', config = config)
        except HTTPException as hte:
            print "Failed submitting task: %s" % (hte.headers)
        except ClientException as cle:
            print "Failed submitting task: %s" % (cle)

    #############################################################################################
    ## From now on that's what users should modify: this is the a-la-CRAB2 configuration part. ##
    #############################################################################################



    config.General.requestName = "W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W1JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W2JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
    submit(config)


    config.General.requestName = "W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W3JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
    submit(config)


    config.General.requestName = "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/MINIAODSIM"
    submit(config)

    config.General.requestName = "W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    config.Data.inputDataset = "/W4JetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext2-v1/MINIAODSIM"
    submit(config)



#    config.General.requestName = "DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
 #   config.Data.inputDataset = "/DY1JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
  #  submit(config)

   # config.General.requestName = "DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
    #config.Data.inputDataset = "/DY2JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
#    submit(config)
    
 #   config.General.requestName = "DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
  #  config.Data.inputDataset = "/DY3JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
   # submit(config)

#    config.General.requestName = "DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8"
 #   config.Data.inputDataset = "/DY4JetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16MiniAODv2-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/MINIAODSIM"
  #  submit(config)
