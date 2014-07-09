'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("L1CustomNtupleProc")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

##process.source = cms.Source("PoolSource",                                                                                                                                   
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

#readFiles.extend( [
#"file:/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/01c5ccb61a90a376764306e94cf36814/output_90_3_UQN.root",
#"file:/hadoop/cms/store/user/olivito/DYJetsToLL_M-50_13TeV-madgraph-pythia8/L1Trig_Flat20to50_bx25_saveRegsTPs_oldGT/e8988c01989bba6871868d8dd7ab7e67/output_34_1_17g.root",
#] )

from input_T2tt_800_100 import input_files
signal_model = 'T2tt'
signal_point = '800_100'

#from input_T2qq_600_100 import input_files
#signal_model = 'T2qq'
#signal_point = '600_100'

#from input_T2qq_450_400 import input_files
#signal_model = 'T2qq'
#signal_point = '450_400'

#from input_T2qq_400_150 import input_files
#signal_model = 'T2qq'
#signal_point = '400_150'

readFiles.extend( input_files )

print readFiles

# Make the framework shut up.
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 10000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

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

# Load emulation and RECO sequences
isMC = True
if not isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
    print "Running on data!"     
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

# Load sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")

###process.CorrectedDigis.puMultCorrect = cms.bool(False) # regions with PU corrections
###process.UCT2015Producer.puMultCorrect = cms.bool(False) # this would do the same trick, applying both to make sure

###process.UCT2015Producer.jetSeed = cms.uint32(10) # study jetSeed 10 --> 5
###process.UCT2015Producer.jetSeed = cms.uint32(5) # study jetSeed 10 --> 5            

## configuration for sums: 7 GeV by default
process.UCT2015Producer.regionETCutForHT = cms.uint32(7)
## eta range:
# 0-21 corresponds to all eta
# 4-17 corresponds to |eta| < 3.0 (default)
# 5-16 corresponds to |eta| < 2.2
process.UCT2015Producer.minGctEtaForSums = cms.uint32(4)
process.UCT2015Producer.maxGctEtaForSums = cms.uint32(17)

process.p1 = cms.Path(
    process.emulationSequence *
    process.uct2015L1Extra
       #  *process.YourFavoritePlottingRoutine  --> This ends at l1extra production, anything after is up to the analyst 
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("ntuple_%s_%s_test.root"%(signal_model,signal_point)),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_l1ntuple_*_L1CustomNtupleProc',
    ) 
)

##process.ep = cms.EndPath(process.output)

process.l1ntuple = cms.EDProducer( 'L1CustomNtupleProducer' ,
                              EcalTriggerPrimitiveInputTag = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                              HcalTriggerPrimitiveInputTag = cms.InputTag("hackHCALMIPs"),
                              L1JetsCentInputTag = cms.InputTag("l1extraParticles","Central","L1CustomNtupleProc"),
                              L1JetsFwdInputTag = cms.InputTag("l1extraParticles","Forward","L1CustomNtupleProc"),
                              L1EtMissInputTag = cms.InputTag("l1extraParticles","MET","L1CustomNtupleProc"),
                              L1MHTInputTag = cms.InputTag("l1extraParticles","MHT"),
                              ###
                              L1JetsCent2015InputTag = cms.InputTag("l1extraParticlesUCT","Central","L1CustomNtupleProc"),
                              L1JetsFwd2015InputTag = cms.InputTag("l1extraParticlesUCT","Forward","L1CustomNtupleProc"),
                              L1EtMiss2015InputTag = cms.InputTag("l1extraParticlesUCT","MET","L1CustomNtupleProc"),
                              L1MHT2015InputTag = cms.InputTag("l1extraParticlesUCT","MHT","L1CustomNtupleProc"),
                              Regions2015InputTag = cms.InputTag("uctDigis","","L1CustomNtupleProc"),
                              CorRegions2015InputTag = cms.InputTag("CorrectedDigis","CorrectedRegions","L1CustomNtupleProc"),
                              ###
                              PUSummaryInfoInputTag = cms.InputTag("addPileupInfo"),
                              ###
                              applyGenWeights = cms.bool(False),
                              doGenKinematics = cms.bool(True),
                              )

process.ntuplePath = cms.Path( process.l1ntuple )
process.outPath = cms.EndPath( process.output )


#file = open('Config_T2tt_800_100_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()
