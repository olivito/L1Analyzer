'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("RateAnalysis")

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

##process.source = cms.Source("PoolSource",                                                                                                                                   
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

#import glob
#fileList = glob.glob('/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/*/*.root')
#for f in fileList:
#    readFiles.append('file:'+f)

readFiles.extend( [
"file:/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/01c5ccb61a90a376764306e94cf36814/output_90_3_UQN.root",
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

process.uctEmulatorStepAgain = cms.Sequence(
#    hackHCALMIPs
    # Now make UCT and L1 objects
#    * uctDigis
     process.CorrectedDigis 
    * process.UCT2015Producer    
)

process.uct2015L1ExtraAgain = cms.Sequence(
#    gtEvmDigis
#    *dttfDigis
#    *csctfDigis
    process.uctGctDigis
    *process.gtUCTDigis
#    *process.scalersRawToDigi
#    *process.l1extraParticles
    *process.l1extraParticlesUCT  
)

process.p2 = cms.Path(
   process.uctEmulatorStepAgain *
   process.uct2015L1ExtraAgain
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_*_*_ReRunningL1',
          'keep *_l1extraParticles*_*_*') 
)

##process.ep = cms.EndPath(process.output)

process.ana = cms.EDAnalyzer( 'L1RateAnalyzer' ,
                              EcalTriggerPrimitiveInputTag = cms.InputTag("ecalDigis","EcalTriggerPrimitives"),
                              HcalTriggerPrimitiveInputTag = cms.InputTag("hackHCALMIPs"),
                              L1JetsCentInputTag = cms.InputTag("l1extraParticles","Central","ReRunningL1"),
                              L1JetsFwdInputTag = cms.InputTag("l1extraParticles","Forward","ReRunningL1"),
                              L1EtMissInputTag = cms.InputTag("l1extraParticles","MET","ReRunningL1"),
                              L1MHTInputTag = cms.InputTag("l1extraParticles","MHT"),
                              ###
                              L1JetsCent2015InputTag = cms.InputTag("l1extraParticlesUCT","Central","RateAnalysis"),
                              L1JetsFwd2015InputTag = cms.InputTag("l1extraParticlesUCT","Forward","RateAnalysis"),
                              L1EtMiss2015InputTag = cms.InputTag("l1extraParticlesUCT","MET","RateAnalysis"),
                              L1MHT2015InputTag = cms.InputTag("l1extraParticlesUCT","MHT","RateAnalysis"),
                              Regions2015InputTag = cms.InputTag("uctDigis","","ReRunningL1"),
                              CorRegions2015InputTag = cms.InputTag("CorrectedDigis","CorrectedRegions","RateAnalysis"),
                              ###
                              PUSummaryInfoInputTag = cms.InputTag("addPileupInfo"),
                              ###
                              regionETCutForHT = process.UCT2015Producer.regionETCutForHT,
                              ## eta range:
                              # 0-21 corresponds to all eta
                              # 4-17 corresponds to |eta| < 3.0 (default)
                              # 5-16 corresponds to |eta| < 2.2
                              minGctEtaForSums = process.UCT2015Producer.minGctEtaForSums,
                              maxGctEtaForSums = process.UCT2015Producer.maxGctEtaForSums,
                              ###
                              doPUPlots = cms.bool(True),
                              doTPPlots = cms.bool(True),
                              applyGenWeights = cms.bool(False),
                              histoname = cms.string("histos_minbias_rate_pu40bx25_mettest.root")
###                              histoname = cms.string("T2tt_800_100_noPU_histos.root")
###                              histoname = cms.string("T2tt_800_100_seed5_histos.root")
                              )

process.pAna = cms.Path( process.ana )


#file = open('Config_T2tt_800_100_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()
