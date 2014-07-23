import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("L1CustomNtupleProc")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

##process.source = cms.Source("PoolSource",                                                                                                                                   
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

#import glob
#fileList = glob.glob('/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/*/*.root')
#for f in fileList:
#    readFiles.append('file:'+f)

readFiles.extend( [
"file:/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_PU40bx25_saveRegsTPs_emulator/85f8a81e12128d811e24c1b783faddb0/output_25_2_JB1.root",
#"file:/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/01c5ccb61a90a376764306e94cf36814/output_90_3_UQN.root",
#"file:/hadoop/cms/store/user/olivito/DYJetsToLL_M-50_13TeV-madgraph-pythia8/L1Trig_Flat20to50_bx25_saveRegsTPs_oldGT/e8988c01989bba6871868d8dd7ab7e67/output_34_1_17g.root",
] )

print readFiles

#from input_800_100 import input_files
#signal_point = '800_100'

#from input_500_250 import input_files
#signal_point = '500_250'

#readFiles.extend( input_files )


# Make the framework shut up.
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

######
######
######

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
## 62X samples, 25ns
process.GlobalTag.globaltag = 'POSTLS162_V2::All'
## 62X samples, 50ns
#process.GlobalTag.globaltag = 'POSTLS162_V1::All'
## 70X samples, 25ns flat PU, running in 710pre8 or higher
#process.GlobalTag.globaltag = 'PRE_LS171V9A::All' ### potentially bad
#process.GlobalTag.globaltag = 'POSTLS170_V5::All'

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('ntuple.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_l1ntuple_*_L1CustomNtupleProc',
    ) 
)

##process.ep = cms.EndPath(process.output)

process.l1ntuple = cms.EDProducer( 'L1CustomNtupleProducer' ,
                              EcalTriggerPrimitiveInputTag = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                              HcalTriggerPrimitiveInputTag = cms.InputTag("hackHCALMIPs"),
                              L1JetsCentInputTag = cms.InputTag("l1extraParticles","Central","L1TEMULATION"),
                              L1JetsFwdInputTag = cms.InputTag("l1extraParticles","Forward","L1TEMULATION"),
                              L1EtMissInputTag = cms.InputTag("l1extraParticles","MET","L1TEMULATION"),
                              L1MHTInputTag = cms.InputTag("l1extraParticles","MHT"),
                              ###
                              L1JetsCent2015InputTag = cms.InputTag("l1extraParticlesEMU","Central","L1TEMULATION"),
                              L1JetsFwd2015InputTag = cms.InputTag("l1extraParticlesEMU","Forward","L1TEMULATION"),
                              L1EtMiss2015InputTag = cms.InputTag("l1extraParticlesEMU","MET","L1TEMULATION"),
                              L1MHT2015InputTag = cms.InputTag("l1extraParticlesEMU","MHT","L1TEMULATION"),
                              Regions2015InputTag = cms.InputTag("uctDigis","","L1TEMULATION"),
                              CorRegions2015InputTag = cms.InputTag("CorrectedDigis","CorrectedRegions","L1TEMULATION"),
                              ###
                              GenJetsInputTag = cms.InputTag("ak4GenJetsNoNu"),
                              ###
                              PUSummaryInfoInputTag = cms.InputTag("addPileupInfo"),
                              ###
                              applyGenWeights = cms.bool(False),
                              doGenKinematics = cms.bool(False),
                              doHTVariations = cms.bool(False),
                              )

process.ntuplePath = cms.Path( process.l1ntuple )
process.outPath = cms.EndPath( process.output )


#file = open('Config_T2tt_800_100_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()
