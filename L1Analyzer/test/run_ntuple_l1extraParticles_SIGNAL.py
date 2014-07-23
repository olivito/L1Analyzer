import FWCore.ParameterSet.Config as cms
import os

process = cms.Process("L1CustomNtupleProc")

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryIdeal_cff')

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()

process.source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles,skipEvents = cms.untracked.uint32(0))

#readFiles.extend( [
#"file:/hadoop/cms/store/user/olivito/Neutrino_Pt-2to20_gun/L1Trig_saveRegsTPs/01c5ccb61a90a376764306e94cf36814/output_90_3_UQN.root",
#"file:/hadoop/cms/store/user/olivito/DYJetsToLL_M-50_13TeV-madgraph-pythia8/L1Trig_Flat20to50_bx25_saveRegsTPs_oldGT/e8988c01989bba6871868d8dd7ab7e67/output_34_1_17g.root",
#] )

#from input_T2tt_800_100 import input_files
#signal_model = 'T2tt'
#signal_point = '800_100'

#from input_T2tt_500_250 import input_files
#signal_model = 'T2tt'
#signal_point = '500_250'

#from input_T2tt_300_100 import input_files
#signal_model = 'T2tt'
#signal_point = '300_100'

#from input_T2tt_300_200 import input_files
#signal_model = 'T2tt'
#signal_point = '300_200'

#from input_T2qq_600_100 import input_files
#signal_model = 'T2qq'
#signal_point = '600_100'

from input_T2qq_450_400 import input_files
signal_model = 'T2qq'
signal_point = '450_400'

#from input_T2qq_450_350 import input_files
#signal_model = 'T2qq'
#signal_point = '450_350'

#from input_T2qq_400_150 import input_files
#signal_model = 'T2qq'
#signal_point = '400_150'

#from input_T2cc_250_190 import input_files
#signal_model = 'T2cc'
#signal_point = '250_190'

#from input_T2cc_250_210 import input_files
#signal_model = 'T2cc'
#signal_point = '250_210'

#from input_T1tttt_1025_625 import input_files
#signal_model = 'T1tttt'
#signal_point = '1025_625'

#from input_T1tttt_825_525 import input_files
#signal_model = 'T1tttt'
#signal_point = '825_525'

#from input_T1qqqq_1200_425 import input_files
#signal_model = 'T1qqqq'
#signal_point = '1200_425'

#from input_T1qqqq_800_575 import input_files
#signal_model = 'T1qqqq'
#signal_point = '800_575'


readFiles.extend( input_files )

print readFiles

# Make the framework shut up.
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

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

# Load sequences
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load('RecoJets.Configuration.RecoGenJets_cff')

process.antiktGenJets = cms.Path(
    process.genParticlesForJetsNoNu*
    process.ak4GenJetsNoNu*
    process.genParticlesForJetsNoMuNoNu*
    process.ak4GenJetsNoMuNoNu)

process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')
process.load('L1Trigger.L1TCalorimeter.caloStage1Params_cfi')

### configuration for sums eta and region Et
#process.caloStage1Params.etSumEtaMin             = cms.vint32(4, 4) #ET, HT
#process.caloStage1Params.etSumEtaMax             = cms.vint32(17, 17) #ET, HT
#process.caloStage1Params.etSumEtThreshold        = cms.vdouble(0., 7.) #ET, HT


# L1TCaloStage1_PPFromRaw = cms.Sequence(
#     L1TRerunHCALTP_FromRAW
#     +ecalDigis
#     +simRctDigis
#     +L1TCaloStage1
#     +simGtDigis
# )

# L1TCaloStage1 = cms.Sequence(
#     rctUpgradeFormatDigis +
#     caloStage1Digis +
#     caloStage1FinalDigis +
#     caloStage1LegacyFormatDigis
# )

process.l1extraParticlesEMU = cms.EDProducer("L1ExtraParticlesProd",
    muonSource = cms.InputTag("simGtDigis"),
    etTotalSource = cms.InputTag("caloStage1LegacyFormatDigis"),
    nonIsolatedEmSource = cms.InputTag("caloStage1LegacyFormatDigis","nonIsoEm"),
    etMissSource = cms.InputTag("caloStage1LegacyFormatDigis"),
    htMissSource = cms.InputTag("caloStage1LegacyFormatDigis"),
    produceMuonParticles = cms.bool(True),
    forwardJetSource = cms.InputTag("caloStage1LegacyFormatDigis","forJets"),
    centralJetSource = cms.InputTag("caloStage1LegacyFormatDigis","cenJets"),
    produceCaloParticles = cms.bool(True),
    tauJetSource = cms.InputTag("caloStage1LegacyFormatDigis","tauJets"),
    isolatedEmSource = cms.InputTag("caloStage1LegacyFormatDigis","isoEm"),
    etHadSource = cms.InputTag("caloStage1LegacyFormatDigis"),
    hfRingEtSumsSource = cms.InputTag("caloStage1LegacyFormatDigis"), # these are empty
    hfRingBitCountsSource = cms.InputTag("caloStage1LegacyFormatDigis"), # these are empty
    centralBxOnly = cms.bool(True),
    ignoreHtMiss = cms.bool(False)
)


process.p1 = cms.Path(
    process.L1TCaloStage1_PPFromRaw*
    process.l1extraParticlesEMU
    )

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("ntuple_%s_%s_emulator_all.root"%(signal_model,signal_point)),
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
                              L1MHTInputTag = cms.InputTag("l1extraParticles","MHT","L1CustomNtupleProc"),
                              ###
                              L1JetsCent2015InputTag = cms.InputTag("l1extraParticlesEMU","Central","L1CustomNtupleProc"),
                              L1JetsFwd2015InputTag = cms.InputTag("l1extraParticlesEMU","Forward","L1CustomNtupleProc"),
                              L1EtMiss2015InputTag = cms.InputTag("l1extraParticlesEMU","MET","L1CustomNtupleProc"),
                              L1MHT2015InputTag = cms.InputTag("l1extraParticlesEMU","MHT","L1CustomNtupleProc"),
                              Regions2015InputTag = cms.InputTag("uctDigis","","L1CustomNtupleProc"),
                              CorRegions2015InputTag = cms.InputTag("CorrectedDigis","CorrectedRegions","L1CustomNtupleProc"),
                              ###
                              GenJetsInputTag = cms.InputTag("ak4GenJetsNoNu"),
                              ###
                              PUSummaryInfoInputTag = cms.InputTag("addPileupInfo"),
                              ###
                              applyGenWeights = cms.bool(False),
                              doGenKinematics = cms.bool(True),
                              doHTVariations = cms.bool(False),
                              )

process.ntuplePath = cms.Path( process.l1ntuple )
process.outPath = cms.EndPath( process.output )


#file = open('Config_T2tt_800_100_cfg.py','w')
#file.write(str(process.dumpPython()))
#file.close()
