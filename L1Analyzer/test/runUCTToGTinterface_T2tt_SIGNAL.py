'''

Creates L1ExtraNtuples (L1 Style) using a UCT->GT jump

Authors: L. Dodd, N. Woods, T. Perry, A. Levine,, S. Dasu, M. Cepeda, E. Friis (UW Madison)

'''

import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("ReRunningL1")

# Get command line options
from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing ('analysis')

options.register(
    'isMC',
    1,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    'Set to 1 for simulated samples - updates GT, emulates HCAL TPGs.')

options.parseArguments()

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

##process.source = cms.Source("PoolSource",                                                                                                                                   
process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

from input_T2tt_800_100 import input_files
signal_model = 'T2tt'
signal_point = '800_100'

#from input_500_250 import input_files
#signal_model = 'T2tt'
#signal_point = '500_250'

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

#readFiles.extend( ['file:/home/users/olivito/l1trig/CMSSW_7_1_0_pre8/src/L1Analyzer/L1Analyzer/test/T2qq_600_100_Fall13_event1050_step1.root'] )


# Make the framework shut up.
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

######
######
######

# Tested on Monte Carlo, for a test with data edit ahead
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V12::All'

# Load emulation and RECO sequences
if not options.isMC:
    process.load("L1Trigger.UCT2015.emulation_cfi")
    print "Running on data!"     
else:
    process.load("L1Trigger.UCT2015.emulationMC_cfi")

# Load sequences
process.load("Configuration.Geometry.GeometryIdeal_cff")
process.load("L1Trigger.UCT2015.uctl1extraparticles_cfi")
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')

process.antiktGenJets = cms.Path(
    process.genParticlesForJetsNoNu*
    process.ak4GenJetsNoNu*
    process.genParticlesForJetsNoMuNoNu*
    process.ak4GenJetsNoMuNoNu)

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
    fileName = cms.untracked.string('out.root'),
    outputCommands = cms.untracked.vstring('drop *',
          'keep *_*_*_ReRunningL1',
          'keep *_l1extraParticles*_*_*') 
)

##process.ep = cms.EndPath(process.output)

process.ana = cms.EDAnalyzer( 'L1Analyzer' ,
                              L1JetsCentInputTag = cms.InputTag("l1extraParticles","Central","ReRunningL1"),
                              L1JetsFwdInputTag = cms.InputTag("l1extraParticles","Forward","ReRunningL1"),
                              L1MHTInputTag = cms.InputTag("l1extraParticles","MHT"),
                              ###
                              L1JetsCent2015InputTag = cms.InputTag("l1extraParticlesUCT","Central","ReRunningL1"),
                              L1JetsFwd2015InputTag = cms.InputTag("l1extraParticlesUCT","Forward","ReRunningL1"),
                              L1MHT2015InputTag = cms.InputTag("l1extraParticlesUCT","MHT","ReRunningL1"),
                              Regions2015InputTag = cms.InputTag("uctDigis","","ReRunningL1"),
                              CorRegions2015InputTag = cms.InputTag("CorrectedDigis","CorrectedRegions","ReRunningL1"),
                              ###
                              GenJetsInputTag = cms.InputTag("ak4GenJetsNoNu"),
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
                              histoname = cms.string("histos_%s_%s_turnons_all.root"%(signal_model,signal_point))
###                              histoname = cms.string("T2tt_800_100_noPU_histos.root")
###                              histoname = cms.string("T2tt_800_100_seed5_histos.root")
                              )

process.pAna = cms.Path( process.ana )


#file = open('Config_T2tt_800_100_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()
