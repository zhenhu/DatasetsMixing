import FWCore.ParameterSet.Config as cms
process = cms.Process("Rootuple")

process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '92X_upgrade2017_realistic_v7', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '94X_dataRun2_ReReco_EOY17_v2', '') #2017 ReReco
#process.GlobalTag.globaltag = '74X_dataRun2_Prompt_v1'  # or some other global tag depending on your CMSSW release and sample. 
#process.GlobalTag = GlobalTag(process.GlobalTag, '102X_dataRun2_Sep2018Rereco_v1', '')	
process.load("FWCore.MessageLogger.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#'file:/eos/uscms/store/user/lpcmuon/fourmuonMC/H0ToUps1SMuMu_m18p5_TuneCUEP8M1_13TeV-pythia8/BPHSkim-v10/180308_232212/0000/BPHSkim_1.root',
#'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/MuOnia/BPHSkim-v4-Run2017F-17Nov2017-v1/180314_061742/0000/BPHSkim_101.root')
#'file:/eos/uscms/store/user/lpcbphy/hui/MuOnia/BPHSkim-2018B_Rereco-Run2018B-17Sep2018-v1/190104_025327/0000/BPHSkim_426.root'),
#'root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/MuOnia_v3/2017_v2/MuOnia/BPHSkim--Run2017F-17Nov2017-v1/190611_021142/0000/BPHSkim_2017_99.root')
#'file:/eos/uscms/store/user/muahmad/FourMuon_Analysis/MuOnia_v3/2017_v2/MuOnia/BPHSkim--Run2017F-17Nov2017-v1/190611_021142/0000/BPHSkim_2017_99.root')
'file:/eos/uscms/store/user/muahmad/FourMuon_Analysis/MuOnia_v3/2017_v2/MuOnia/BPHSkim--Run2017D-17Nov2017-v1/190611_021109/0000/BPHSkim_2017_40.root')
#eventsToProcess = cms.untracked.VEventRange('317641:1352331850-317641:1352331850')
)


process.TFileService = cms.Service("TFileService",
        fileName = cms.string('Rootuple.root'),
)

process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True))

process.load('DatasetsMixing.Onia.MixingRootupler_cfi')
process.p = cms.Path(process.rootuple)

process.rootuple.upsilon_mass = cms.double(9.4603)
process.rootuple.triggerCuts = cms.uint32(36)
process.rootuple.isMC = cms.bool(False)                 # is mc?
process.rootuple.onia_mass_cuts = cms.vdouble(0.1,500)    # you may need to adjust this
process.rootuple.SecondSource.fileNames = cms.untracked.vstring(
'file:/eos/uscms/store/user/muahmad/FourMuon_Analysis/MiniBias/2017_v2/ZeroBias/BPHSkim--Run2017D-17Nov2017-v1/190609_091210/0000/BPHSkim_2017_MiniBias_203.root')
#'root://cmsxrootd.fnal.gov//store/user/muahmad/FourMuon_Analysis/MiniBias/2017_v2/ZeroBias/BPHSkim--Run2017B-17Nov2017-v1/190609_091136/0000/BPHSkim_2017_MiniBias_97.root')
#'file:/eos/uscms/store/user/muahmad/FourMuon_Analysis/MiniBias/2017_v2/ZeroBias/BPHSkim--Run2017F-17Nov2017-v1/190609_091305/0000/BPHSkim_2017_MiniBias_97.root')
#'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/ZeroBias/BPHSkim-v6-Run2017F-17Nov2017-v1/180321_073042/0000/BPHSkim_250.root',
#'file:/eos/uscms/store/user/l1upgrades/Run2017/fourmuon/ZeroBias/BPHSkim-v6-Run2017F-17Nov2017-v1/180321_073042/0000/BPHSkim_250.root'
#    'file:/eos/cms/store/user/jblee/SpyFEDemulated234824.root'
	  #)



