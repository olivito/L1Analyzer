// -*- C++ -*-
//
// Package:    L1TrackTriggerObjectsAnalyzer 
// Class:      L1TrackTriggerObjectsAnalyzer
//
/**\class L1TrackTriggerObjectsAnalyzer L1TrackTriggerObjectsAnalyzer.cc SLHCUpgradeSimulations/L1TrackTriggerObjectsAnalyzer/src/L1TrackTriggerObjectsAnalyzer.cc 
 Description: [one line class summary] 

 Implementation: 
     [Notes on implementation]
*/
//
// Original Author:  Maria
//         Created:  Thu Nov 14 11:22:13 CET 2013 
// $Id$ 
//
//

// system include files 
#include <memory>
//#include <string>

// user include files 
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"

#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1CaloTrigger/interface/L1CaloCollections.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegion.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloEmCand.h"
#include "DataFormats/L1CaloTrigger/interface/L1CaloRegionDetId.h"

#include "TFile.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// MT2 stuff
#include "L1Analyzer/L1Analyzer/interface/Davismt2.h"
#include "L1Analyzer/L1Analyzer/interface/Hemisphere.hh"

using namespace l1extra;

using namespace reco;
using namespace std;

typedef math::XYZTLorentzVector LorentzVector;

//////// 
// class declaration
////

/*
class cmpPtGT {
public:
  bool operator() (const L1TkJetParticle &r1, const L1TkJetParticle &r2) {
    return r1.pt()>r2.pt(); 
  }
};
*/

class cmpPtLorentz {
public:
  bool operator() (const LorentzVector &r1, const LorentzVector &r2) {
    return r1.pt()>r2.pt();
  }
};

/////  

class L1CustomNtupleProducer : public edm::EDProducer {
public:

  explicit L1CustomNtupleProducer(const edm::ParameterSet&);
  ~L1CustomNtupleProducer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void produce(edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  Int_t phiINjetCoord(Double_t phi);    
  Int_t etaINjetCoord(Double_t eta);
  inline Double_t degree(Double_t radian);
  float dphi_normal(float phi1, float phi2);
  float dphi(float phi1, float phi2);
  float deltaR_normal(float eta1, float phi1, float eta2, float phi2);
  Float_t CalcHemisphereAndMT2(float testmass, bool massive, std::vector<LorentzVector> jets, LorentzVector MET );
  Float_t CalcMT2(float testmass, bool massive, LorentzVector visible1, LorentzVector visible2, LorentzVector MET );

  double regionPhysicalEt(const L1CaloRegion& cand) const {
    double regionLSB_=0.5;
    //    double regionLSB_=1.0;
    return std::max(0.,regionLSB_*cand.et());
  }


  // ----------member data ---------------------------

  const std::string aliasprefix_ = "l1";
                                                                                                                                          
  edm::InputTag L1JetsCentInputTag;
  edm::InputTag L1JetsFwdInputTag;
  edm::InputTag L1MHTInputTag;
  //
  edm::InputTag L1JetsCent2015InputTag;
  edm::InputTag L1JetsFwd2015InputTag;
  edm::InputTag L1EtMiss2015InputTag;
  edm::InputTag L1MHT2015InputTag;
  edm::InputTag Regions2015InputTag;
  edm::InputTag CorRegions2015InputTag;

  edm::InputTag GenJetsInputTag;

  edm::InputTag PUSummaryInfoInputTag;

  bool doGenKinematics_;
  bool applyGenWeights_;
  bool doHTVariations_;

  bool doDebug=false; 

  const static size_t PHIBINS = 18;
  const Double_t PHIBIN[PHIBINS] = {10,30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350};
  const static size_t ETABINS = 23;
  const Double_t ETABIN[ETABINS] = {-5.,-4.5,-4.,-3.5,-3.,-2.172,-1.74,-1.392,-1.044,-0.696,-0.348,0,0.348,0.696,1.044,1.392,1.74,2.172,3.,3.5,4.,4.5,5.};

};

//
// constants, enums and typedefs 
//

//
// static data member definitions
//

//
// constructors and destructor
//
L1CustomNtupleProducer::L1CustomNtupleProducer(const edm::ParameterSet& iConfig)

{
  L1JetsCentInputTag = iConfig.getParameter<edm::InputTag>("L1JetsCentInputTag");
  L1JetsFwdInputTag = iConfig.getParameter<edm::InputTag>("L1JetsFwdInputTag");
  L1MHTInputTag = iConfig.getParameter<edm::InputTag>("L1MHTInputTag");
  //
  L1JetsCent2015InputTag = iConfig.getParameter<edm::InputTag>("L1JetsCent2015InputTag");
  L1JetsFwd2015InputTag = iConfig.getParameter<edm::InputTag>("L1JetsFwd2015InputTag");
  L1EtMiss2015InputTag = iConfig.getParameter<edm::InputTag>("L1EtMiss2015InputTag");
  L1MHT2015InputTag= iConfig.getParameter<edm::InputTag>("L1MHT2015InputTag");
  Regions2015InputTag= iConfig.getParameter<edm::InputTag>("Regions2015InputTag");
  CorRegions2015InputTag= iConfig.getParameter<edm::InputTag>("CorRegions2015InputTag");
  //
  GenJetsInputTag = iConfig.getParameter<edm::InputTag>("GenJetsInputTag");
  //
  PUSummaryInfoInputTag = iConfig.getParameter<edm::InputTag>("PUSummaryInfoInputTag");

  doGenKinematics_ = iConfig.getParameter<bool>("doGenKinematics");
  applyGenWeights_= iConfig.getParameter<bool>("applyGenWeights"); 
  doHTVariations_ = iConfig.getParameter<bool>("doHTVariations");

  // output branches
  produces<float>  (aliasprefix_+"jet1ptEta30").setBranchAlias(aliasprefix_+"_jet1ptEta30");
  produces<float>  (aliasprefix_+"jet2ptEta30").setBranchAlias(aliasprefix_+"_jet2ptEta30");
  produces<float>  (aliasprefix_+"jet3ptEta30").setBranchAlias(aliasprefix_+"_jet3ptEta30");
  produces<float>  (aliasprefix_+"jet4ptEta30").setBranchAlias(aliasprefix_+"_jet4ptEta30");
  produces<float>  (aliasprefix_+"dphiJet12Eta30").setBranchAlias(aliasprefix_+"_dphiJet12Eta30");

  produces<float>  (aliasprefix_+"jet1ptEta22").setBranchAlias(aliasprefix_+"_jet1ptEta22");
  produces<float>  (aliasprefix_+"jet2ptEta22").setBranchAlias(aliasprefix_+"_jet2ptEta22");
  produces<float>  (aliasprefix_+"jet3ptEta22").setBranchAlias(aliasprefix_+"_jet3ptEta22");
  produces<float>  (aliasprefix_+"jet4ptEta22").setBranchAlias(aliasprefix_+"_jet4ptEta22");
  produces<float>  (aliasprefix_+"dphiJet12Eta22").setBranchAlias(aliasprefix_+"_dphiJet12Eta22");

  if (doHTVariations_) {
    produces<float>  (aliasprefix_+"htEta30RegEt7").setBranchAlias(aliasprefix_+"_htEta30RegEt7");
    produces<float>  (aliasprefix_+"htEta30RegEt10").setBranchAlias(aliasprefix_+"_htEta30RegEt10");
    produces<float>  (aliasprefix_+"htEta30RegEt15").setBranchAlias(aliasprefix_+"_htEta30RegEt15");
    produces<float>  (aliasprefix_+"htEta30RegEt20").setBranchAlias(aliasprefix_+"_htEta30RegEt20");

    produces<float>  (aliasprefix_+"htEta22RegEt7").setBranchAlias(aliasprefix_+"_htEta22RegEt7");
    produces<float>  (aliasprefix_+"htEta22RegEt10").setBranchAlias(aliasprefix_+"_htEta22RegEt10");
    produces<float>  (aliasprefix_+"htEta22RegEt15").setBranchAlias(aliasprefix_+"_htEta22RegEt15");
    produces<float>  (aliasprefix_+"htEta22RegEt20").setBranchAlias(aliasprefix_+"_htEta22RegEt20");
  }

  produces<float>  (aliasprefix_+"met").setBranchAlias(aliasprefix_+"_met");
  produces<float>  (aliasprefix_+"ht").setBranchAlias(aliasprefix_+"_ht");
  produces<float>  (aliasprefix_+"mht").setBranchAlias(aliasprefix_+"_mht");

  produces<float>  (aliasprefix_+"weight").setBranchAlias(aliasprefix_+"_weight");
  produces<float>  (aliasprefix_+"ntrueint").setBranchAlias(aliasprefix_+"_ntrueint");

  if (doGenKinematics_) {
    produces<unsigned int>  (aliasprefix_+"ngenlep").setBranchAlias(aliasprefix_+"_ngenlep");

    produces<float>  (aliasprefix_+"genjet1ptEta30").setBranchAlias(aliasprefix_+"_genjet1ptEta30");
    produces<float>  (aliasprefix_+"genjet2ptEta30").setBranchAlias(aliasprefix_+"_genjet2ptEta30");
    produces<float>  (aliasprefix_+"genjet3ptEta30").setBranchAlias(aliasprefix_+"_genjet3ptEta30");
    produces<float>  (aliasprefix_+"genjet4ptEta30").setBranchAlias(aliasprefix_+"_genjet4ptEta30");
    produces<float>  (aliasprefix_+"genDphiJet12Eta30").setBranchAlias(aliasprefix_+"_genDphiJet12Eta30");

    produces<unsigned int>  (aliasprefix_+"ngenjetsPt40Eta30").setBranchAlias(aliasprefix_+"_ngenjetsPt40Eta30");

    produces<float>  (aliasprefix_+"genjet1ptEta22").setBranchAlias(aliasprefix_+"_genjet1ptEta22");
    produces<float>  (aliasprefix_+"genjet2ptEta22").setBranchAlias(aliasprefix_+"_genjet2ptEta22");
    produces<float>  (aliasprefix_+"genjet3ptEta22").setBranchAlias(aliasprefix_+"_genjet3ptEta22");
    produces<float>  (aliasprefix_+"genjet4ptEta22").setBranchAlias(aliasprefix_+"_genjet4ptEta22");
    produces<float>  (aliasprefix_+"genDphiJet12Eta22").setBranchAlias(aliasprefix_+"_genDphiJet12Eta22");

    produces<unsigned int>  (aliasprefix_+"ngenjetsPt40Eta22").setBranchAlias(aliasprefix_+"_ngenjetsPt40Eta22");

    produces<float>  (aliasprefix_+"genhtEta30").setBranchAlias(aliasprefix_+"_genhtEta30");
    produces<float>  (aliasprefix_+"genhtEta22").setBranchAlias(aliasprefix_+"_genhtEta22");

    produces<float>  (aliasprefix_+"genmet").setBranchAlias(aliasprefix_+"_genmet");
    produces<float>  (aliasprefix_+"genmetEta30").setBranchAlias(aliasprefix_+"_genmetEta30");

    produces<float>  (aliasprefix_+"genmhtEta30").setBranchAlias(aliasprefix_+"_genmhtEta30");
    produces<float>  (aliasprefix_+"genmhtEta22").setBranchAlias(aliasprefix_+"_genmhtEta22");

    produces<float>  (aliasprefix_+"genmt2Eta30").setBranchAlias(aliasprefix_+"_genmt2Eta30");
    produces<float>  (aliasprefix_+"genmt2Eta22").setBranchAlias(aliasprefix_+"_genmt2Eta22");
  }

}


L1CustomNtupleProducer::~L1CustomNtupleProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1CustomNtupleProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  if(doDebug) std::cout << " ----  a new event ----- " << std::endl;
 
  // ----------------------------
  // branch variable defs
  auto_ptr<float>   jet1ptEta30    (new float);
  auto_ptr<float>   jet2ptEta30    (new float);
  auto_ptr<float>   jet3ptEta30    (new float);
  auto_ptr<float>   jet4ptEta30    (new float);
  auto_ptr<float>   dphiJet12Eta30    (new float);

  auto_ptr<float>   jet1ptEta22    (new float);
  auto_ptr<float>   jet2ptEta22    (new float);
  auto_ptr<float>   jet3ptEta22    (new float);
  auto_ptr<float>   jet4ptEta22    (new float);
  auto_ptr<float>   dphiJet12Eta22    (new float);

  auto_ptr<float>   htEta30RegEt7     (new float);
  auto_ptr<float>   htEta30RegEt10    (new float);
  auto_ptr<float>   htEta30RegEt15    (new float);
  auto_ptr<float>   htEta30RegEt20    (new float);

  auto_ptr<float>   htEta22RegEt7     (new float);
  auto_ptr<float>   htEta22RegEt10    (new float);
  auto_ptr<float>   htEta22RegEt15    (new float);
  auto_ptr<float>   htEta22RegEt20    (new float);

  auto_ptr<float>   met           (new float);
  auto_ptr<float>   ht           (new float);
  auto_ptr<float>   mht           (new float);

  auto_ptr<float>   weight    (new float);
  auto_ptr<float>   ntrueint    (new float);

  // if doGenKinematics
  auto_ptr<unsigned int>   ngenlep    (new unsigned int);

  auto_ptr<float>   genjet1ptEta30    (new float);
  auto_ptr<float>   genjet2ptEta30    (new float);
  auto_ptr<float>   genjet3ptEta30    (new float);
  auto_ptr<float>   genjet4ptEta30    (new float);
  auto_ptr<float>   genDphiJet12Eta30    (new float);

  auto_ptr<unsigned int>   ngenjetsPt40Eta30    (new unsigned int);

  auto_ptr<float>   genjet1ptEta22    (new float);
  auto_ptr<float>   genjet2ptEta22    (new float);
  auto_ptr<float>   genjet3ptEta22    (new float);
  auto_ptr<float>   genjet4ptEta22    (new float);
  auto_ptr<float>   genDphiJet12Eta22    (new float);

  auto_ptr<unsigned int>   ngenjetsPt40Eta22    (new unsigned int);

  auto_ptr<float>   genhtEta30     (new float);
  auto_ptr<float>   genhtEta22     (new float);

  auto_ptr<float>   genmet           (new float);
  auto_ptr<float>   genmetEta30     (new float);

  auto_ptr<float>   genmhtEta30     (new float);
  auto_ptr<float>   genmhtEta22     (new float);

  auto_ptr<float>   genmt2Eta30     (new float);
  auto_ptr<float>   genmt2Eta22     (new float);
  // ----------------------------

  *weight = 1.0;  

  // gen level weights
  if (applyGenWeights_) {
    Handle<GenEventInfoProduct> evt_info;
    iEvent.getByLabel("generator", evt_info);
    if(evt_info->weight()>0) *weight = evt_info->weight();       
  }
  
  //
  // ----------------------------------------------------------------------
  // retrieve the PU summary info
  //

  edm::Handle<std::vector<PileupSummaryInfo> > PUSummaryInfoHandle;
  iEvent.getByLabel(PUSummaryInfoInputTag, PUSummaryInfoHandle);

  *ntrueint = -1.;
  if (PUSummaryInfoHandle.isValid()) {
    *ntrueint = PUSummaryInfoHandle->at(0).getTrueNumInteractions();
  } else {
    std::cout << "WARNING: no PU info found! tag was: " << PUSummaryInfoInputTag.encode() << std::endl;
  }


  // gen level quantities
  if (doGenKinematics_) {

  //
  // ----------------------------------------------------------------------
  // retrieve the genParticle objects
  //

    edm::Handle<reco::GenParticleCollection> genparts;
    iEvent.getByLabel("genParticles", genparts);
    // look for prompt lepton from W 
    // recompute the genMET with particle |eta|<3                           
                       
    //  cout << "------- " << endl;

    std::vector<LorentzVector> lsps;
    LorentzVector genMETvec;
    LorentzVector genMETvec_eta30;
    LorentzVector genMETvec_eta22;

    // this is generator-specific and current works for Pythia8!
    // See pythia8xxx/htmldoc/ParticleProperties.html for status codes
    // 21: incoming parton
    // 22: intermediate particle in hard process
    // 23: outgoing parton
                                                                                        
    *ngenlep = 0;
    for (reco::GenParticleCollection::const_iterator genp = genparts->begin(); genp != genparts->end(); ++ genp ) {

      //    if(genp->status() == 23 || genp->status() == 22 ) cout << "genp->pdgId() " << genp->pdgId() << " genp->status() " << genp->status() << endl;

      //    if(abs(genp->pdgId()) == 15) cout << "genp->pdgId() " << genp->pdgId() << " genp->status() " << genp->status() << endl;

      //    if ((abs(genp->pdgId()) == 11 || abs(genp->pdgId()) == 12 || abs(genp->pdgId()) == 13 || abs(genp->pdgId()) == 14 || abs(genp->pdgId()) == 15 || abs(genp->pdgId()) == 16 ) && genp->status() == 23) nGenLep++;

      if ((genp->status() == 3 || genp->status() == 23) && (abs(genp->pdgId()) == 11 || abs(genp->pdgId()) == 13 || abs(genp->pdgId()) == 15)) (*ngenlep)++;

      if((abs(genp->pdgId()) == 1000022) && (genp->status() == 1)) {
	lsps.push_back(genp->p4());
      }

      if(genp->status()!=1) continue;
      if(abs(genp->pdgId()) == 13) continue; // remove muons 
      if(abs(genp->pdgId()) == 12 || abs(genp->pdgId()) == 14 || abs(genp->pdgId()) == 16) continue; // remove neutrinos
      if(abs(genp->pdgId()) == 1000022) continue; // remove LSPs
      if(abs(genp->pdgId()) == 1000012 || abs(genp->pdgId())==1000014 || abs(genp->pdgId())==1000016 ) continue;
      if(abs(genp->pdgId()) == 2000012 || abs(genp->pdgId())==2000014 || abs(genp->pdgId())==2000016 ) continue;
      if(abs(genp->pdgId()) == 1000039 || abs(genp->pdgId())==5100039 ) continue; // remove gravitinos
      if(abs(genp->pdgId()) == 4000012 || abs(genp->pdgId())==4000014 || abs(genp->pdgId())==4000016 ) continue;
      if(abs(genp->pdgId()) == 9900012 || abs(genp->pdgId())==9900014 || abs(genp->pdgId())==9900016 ) continue;
      if(abs(genp->pdgId()) == 39) continue;

      /// https://cmssdt.cern.ch/SDT/lxr/source/RecoMET/Configuration/python/GenMETParticles_cff.py#056
      LorentzVector p3(genp->px(),genp->py(),genp->pz(),genp->energy());
      genMETvec -= p3;

      if(fabs(genp->eta()) < 3.0) genMETvec_eta30 -= p3;
      if(fabs(genp->eta()) < 2.172) genMETvec_eta22 -= p3;

    }

    *genmet = genMETvec.pt();
    *genmetEta30 = genMETvec_eta30.pt();
    //  float genmet_eta22 = genMETvec_eta22.pt();

    // require hadronic events
    //  if (ngenlep > 0) return;

    // --------------------------------------------------------------------------
    // retrieve the genJET objects
    //

    edm::Handle<reco::GenJetCollection> genJetsHandle;
    iEvent.getByLabel(GenJetsInputTag, genJetsHandle);

    *genhtEta30 = 0.;
    *genhtEta22 = 0.;
    *ngenjetsPt40Eta30 = 0;
    *ngenjetsPt40Eta22 = 0;
    *genmhtEta30 = 0.;
    *genmhtEta22 = 0.;

    std::vector<LorentzVector> genJets;
    std::vector<LorentzVector> genJetsCent30;
    std::vector<LorentzVector> genJetsCent22;

    LorentzVector genMHTvec_eta30;
    LorentzVector genMHTvec_eta22;

    for (unsigned iGenJet = 0; iGenJet < genJetsHandle->size(); ++iGenJet) {
      const reco::GenJet& genJet = (*genJetsHandle) [iGenJet];
    
      genJets.push_back(genJet.p4());
      if(fabs(genJet.eta())<3.0) genJetsCent30.push_back(genJet.p4());
      if(fabs(genJet.eta())<2.172) genJetsCent22.push_back(genJet.p4());

      float genJetPt  = genJet.pt();
      float genJetEta = genJet.eta();

      // should apply some smearing to get jet pt to mock up reco?

      LorentzVector p3(genJet.px(),genJet.py(),genJet.pz(),genJet.energy());

      if( genJetPt > 40.0 && fabs(genJetEta)<3.0 ) {
	*genhtEta30 += genJetPt;
	++(*ngenjetsPt40Eta30);
	genMHTvec_eta30 -= p3;
      }
      if( genJetPt > 40.0 && fabs(genJetEta)<2.172 ) {
	*genhtEta22 += genJetPt;
	++(*ngenjetsPt40Eta22);
	genMHTvec_eta22 -= p3;
      }
    } // loop on genjets

    *genmhtEta30 = genMHTvec_eta30.pt();
    *genmhtEta22 = genMHTvec_eta22.pt();

    sort(genJets.begin(),genJets.end(),cmpPtLorentz());
    sort(genJetsCent30.begin(),genJetsCent30.end(),cmpPtLorentz());
    sort(genJetsCent22.begin(),genJetsCent22.end(),cmpPtLorentz());
  
    *genjet1ptEta30 = 0.;
    *genjet2ptEta30 = 0.;
    *genjet3ptEta30 = 0.;
    *genjet4ptEta30 = 0.;
    *genDphiJet12Eta30 = -99.;

    if (genJetsCent30.size() > 0) {
      *genjet1ptEta30 = genJetsCent30.at(0).pt();
    }
    if (genJetsCent30.size() > 1) {
      *genjet2ptEta30 = genJetsCent30.at(1).pt();
      *genDphiJet12Eta30 = dphi_normal(genJetsCent30.at(0).phi(),genJetsCent30.at(1).phi());
    }
    if (genJetsCent30.size() > 2) {
      *genjet3ptEta30 = genJetsCent30.at(2).pt();
    }
    if (genJetsCent30.size() > 3) {
      *genjet4ptEta30 = genJetsCent30.at(3).pt();
    }

    *genjet1ptEta22 = 0.;
    *genjet2ptEta22 = 0.;
    *genjet3ptEta22 = 0.;
    *genjet4ptEta22 = 0.;
    *genDphiJet12Eta22 = -99.;

    if (genJetsCent22.size() > 0) {
      *genjet1ptEta22 = genJetsCent22.at(0).pt();
    }
    if (genJetsCent22.size() > 1) {
      *genjet2ptEta22 = genJetsCent22.at(1).pt();
      *genDphiJet12Eta22 = dphi_normal(genJetsCent22.at(0).phi(),genJetsCent22.at(1).phi());
    }
    if (genJetsCent22.size() > 2) {
      *genjet3ptEta22 = genJetsCent22.at(2).pt();
    }
    if (genJetsCent22.size() > 3) {
      *genjet4ptEta22 = genJetsCent22.at(3).pt();
    }

    // calculate genMT2
    *genmt2Eta30 = CalcHemisphereAndMT2(0.,false,genJetsCent30,genMETvec);
    *genmt2Eta22 = CalcHemisphereAndMT2(0.,false,genJetsCent22,genMETvec);
  }

  //
  // ----------------------------------------------------------------------
  // retrieve the jets objects
  //

  edm::Handle<L1JetParticleCollection> L1JetsCent2015Handle;
  iEvent.getByLabel(L1JetsCent2015InputTag, L1JetsCent2015Handle);
  edm::Handle<L1JetParticleCollection> L1JetsFwd2015Handle;
  iEvent.getByLabel(L1JetsFwd2015InputTag, L1JetsFwd2015Handle);

  std::vector<LorentzVector> newJets2015;
  std::vector<LorentzVector> newJets2015Cent30;
  std::vector<LorentzVector> newJets2015Cent22;

  if ( L1JetsCent2015Handle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsCent2015Handle -> begin(); jetIter != L1JetsCent2015Handle->end(); ++jetIter) {
      // float et = jetIter -> pt();
      float eta = jetIter -> eta();
      // float phi = jetIter -> phi();

      // int etabin = etaINjetCoord(eta);
      // int phibin = phiINjetCoord(phi);

      newJets2015.push_back(jetIter->p4());
      if (fabs(eta) < 3.0) {
	newJets2015Cent30.push_back(jetIter->p4());
      }
      if (fabs(eta) < 2.172) {
	newJets2015Cent22.push_back(jetIter->p4());
      }

    } // loop over jets
  }

  // loop also over forward jets
  if ( L1JetsFwd2015Handle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsFwd2015Handle -> begin(); jetIter != L1JetsFwd2015Handle->end(); ++jetIter) {
      // float et = jetIter -> pt();
      // float eta = jetIter -> eta();
      // float phi = jetIter -> phi();

      // int etabin = etaINjetCoord(eta);
      // int phibin = phiINjetCoord(phi);

      // combine with central jets
      newJets2015.push_back(jetIter->p4());
    }
  }

  sort(newJets2015.begin(),newJets2015.end(),cmpPtLorentz());
  sort(newJets2015Cent30.begin(),newJets2015Cent30.end(),cmpPtLorentz());
  sort(newJets2015Cent22.begin(),newJets2015Cent22.end(),cmpPtLorentz());

  *jet1ptEta30 = 0.;
  *jet2ptEta30 = 0.;
  *jet3ptEta30 = 0.;
  *jet4ptEta30 = 0.;
  *dphiJet12Eta30 = -99.;

  if (newJets2015Cent30.size() > 0) {
    *jet1ptEta30 = newJets2015Cent30.at(0).pt();
  }
  if (newJets2015Cent30.size() > 1) {
    *jet2ptEta30 = newJets2015Cent30.at(1).pt();
    *dphiJet12Eta30 = dphi(newJets2015Cent30.at(0).phi(),newJets2015Cent30.at(1).phi());
  }
  if (newJets2015Cent30.size() > 2) {
    *jet3ptEta30 = newJets2015Cent30.at(2).pt();
  }
  if (newJets2015Cent30.size() > 3) {
    *jet4ptEta30 = newJets2015Cent30.at(3).pt();
  }

  *jet1ptEta22 = 0.;
  *jet2ptEta22 = 0.;
  *jet3ptEta22 = 0.;
  *jet4ptEta22 = 0.;
  *dphiJet12Eta22 = -99.;

  if (newJets2015Cent22.size() > 0) {
    *jet1ptEta22 = newJets2015Cent22.at(0).pt();
  }
  if (newJets2015Cent22.size() > 1) {
    *jet2ptEta22 = newJets2015Cent22.at(1).pt();
    *dphiJet12Eta22 = dphi(newJets2015Cent22.at(0).phi(),newJets2015Cent22.at(1).phi());
  }
  if (newJets2015Cent22.size() > 2) {
    *jet3ptEta22 = newJets2015Cent22.at(2).pt();
  }
  if (newJets2015Cent22.size() > 3) {
    *jet4ptEta22 = newJets2015Cent22.at(3).pt();
  }


  // //
  // // ----------------------------------------------------------------------
  // // retrieve the caloMHT objects
  // //

  // edm::Handle<L1EtMissParticleCollection> L1MHTHandle;
  // iEvent.getByLabel(L1MHTInputTag, L1MHTHandle);
  
  // float caloMHt_=0;
  // float caloHt_=0;
  
  // if (L1MHTHandle.isValid() ) {
  //   for (L1EtMissParticleCollection::const_iterator mhtIter = L1MHTHandle->begin(); mhtIter != L1MHTHandle->end(); ++mhtIter) {
  //     caloMHt_ = mhtIter -> etMiss();
  //     caloHt_ = mhtIter -> etTotal();

  //     //      cout << "overflowHT " << mhtIter->gctEtMiss()->overFlow() << endl;
  //     //      if(mhtIter->gctEtMiss()->overFlow()>0) cout << "!!!!!!!!!!!!! HT BAD OVERFLOW "<< endl;
  //   }
  // }
  
  // plot1D("h_caloMHT",caloMHt_,1,h_1d, nBinMet,0.,maxMet);
  // plot1D("h_caloHT",caloHt_,1,h_1d, nBinHt,0.,maxHt);

  //
  // ----------------------------------------------------------------------
  // retrieve the caloMHT objects 2015
  //

  *mht = 0.;
  *ht = 0.;

  edm::Handle<L1EtMissParticleCollection> L1MHT2015Handle;
  iEvent.getByLabel(L1MHT2015InputTag, L1MHT2015Handle);
  
  if (L1MHT2015Handle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator mhtIter = L1MHT2015Handle->begin(); mhtIter != L1MHT2015Handle->end(); ++mhtIter) {
      *mht = mhtIter -> etMiss();
      *ht = mhtIter -> etTotal();
    }
  }
  
  //
  // ----------------------------------------------------------------------
  // retrieve the caloMET objects 2015
  //

  *met=0;

  edm::Handle<L1EtMissParticleCollection> L1EtMiss2015Handle;
  iEvent.getByLabel(L1EtMiss2015InputTag, L1EtMiss2015Handle);

  if (L1EtMiss2015Handle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator caloEtmIter = L1EtMiss2015Handle->begin(); caloEtmIter != L1EtMiss2015Handle->end(); ++caloEtmIter) {
      *met=caloEtmIter -> et();
    }
  }

  // ----------------------------------------------------------------------
  // input to the sums, after PU correction
  //

  if (doHTVariations_) {

    Handle<L1CaloRegionCollection> corRegions;
    iEvent.getByLabel(CorRegions2015InputTag, corRegions);

    float htEta30RegEt0 = 0.;
    float htEta30RegEt3 = 0.;
    *htEta30RegEt7 = 0.;
    *htEta30RegEt10 = 0.;
    *htEta30RegEt15 = 0.;
    *htEta30RegEt20 = 0.;
    float ht_alletaRegEt7 = 0.;
    float htEta22RegEt0 = 0.;
    float htEta22RegEt3 = 0.;
    *htEta22RegEt7 = 0.;
    *htEta22RegEt10 = 0.;
    *htEta22RegEt15 = 0.;
    *htEta22RegEt20 = 0.;

    for(L1CaloRegionCollection::const_iterator corRegion = corRegions->begin();
	corRegion != corRegions->end(); corRegion++){

      float regionET =  regionPhysicalEt(*corRegion);

      if (regionET >= 7.) {
	ht_alletaRegEt7 += regionET;
      }

      // eta range:
      // 0-21 corresponds to all eta
      // 4-17 corresponds to |eta| < 3.0 (default)
      // 5-16 corresponds to |eta| < 2.2
      if (corRegion->gctEta() >= 4 && corRegion->gctEta() <= 17) {
	htEta30RegEt0 += regionET;
	if (regionET >= 3.) {
	  htEta30RegEt3 += regionET;
	}
	if (regionET >= 7.) {
	  *htEta30RegEt7 += regionET;
	}
	if (regionET >= 10.) {
	  *htEta30RegEt10 += regionET;
	}
	if (regionET >= 15.) {
	  *htEta30RegEt15 += regionET;
	}
	if (regionET >= 20.) {
	  *htEta30RegEt20 += regionET;
	}
      }

      if (corRegion->gctEta() >= 5 && corRegion->gctEta() <= 16) {
	htEta22RegEt0 += regionET;
	if (regionET >= 3.) {
	  htEta22RegEt3 += regionET;
	}
	if (regionET >= 7.) {
	  *htEta22RegEt7 += regionET;
	}
	if (regionET >= 10.) {
	  *htEta22RegEt10 += regionET;
	}
	if (regionET >= 15.) {
	  *htEta22RegEt15 += regionET;
	}
	if (regionET >= 20.) {
	  *htEta22RegEt20 += regionET;
	}
      }

    } // loop over regions

  } // if (doHTVariations_)

  // place branches into the event
  iEvent.put(jet1ptEta30,aliasprefix_+"jet1ptEta30");
  iEvent.put(jet2ptEta30,aliasprefix_+"jet2ptEta30");
  iEvent.put(jet3ptEta30,aliasprefix_+"jet3ptEta30");
  iEvent.put(jet4ptEta30,aliasprefix_+"jet4ptEta30");
  iEvent.put(dphiJet12Eta30,aliasprefix_+"dphiJet12Eta30");

  iEvent.put(jet1ptEta22,aliasprefix_+"jet1ptEta22");
  iEvent.put(jet2ptEta22,aliasprefix_+"jet2ptEta22");
  iEvent.put(jet3ptEta22,aliasprefix_+"jet3ptEta22");
  iEvent.put(jet4ptEta22,aliasprefix_+"jet4ptEta22");
  iEvent.put(dphiJet12Eta22,aliasprefix_+"dphiJet12Eta22");

  if (doHTVariations_) {
    iEvent.put(htEta30RegEt7,aliasprefix_+"htEta30RegEt7");
    iEvent.put(htEta30RegEt10,aliasprefix_+"htEta30RegEt10");
    iEvent.put(htEta30RegEt15,aliasprefix_+"htEta30RegEt15");
    iEvent.put(htEta30RegEt20,aliasprefix_+"htEta30RegEt20");

    iEvent.put(htEta22RegEt7,aliasprefix_+"htEta22RegEt7");
    iEvent.put(htEta22RegEt10,aliasprefix_+"htEta22RegEt10");
    iEvent.put(htEta22RegEt15,aliasprefix_+"htEta22RegEt15");
    iEvent.put(htEta22RegEt20,aliasprefix_+"htEta22RegEt20");
  }

  iEvent.put(met,aliasprefix_+"met");
  iEvent.put(ht,aliasprefix_+"ht");
  iEvent.put(mht,aliasprefix_+"mht");

  iEvent.put(weight,aliasprefix_+"weight");
  iEvent.put(ntrueint,aliasprefix_+"ntrueint");

  if (doGenKinematics_) {
    iEvent.put(ngenlep,aliasprefix_+"ngenlep");

    iEvent.put(genjet1ptEta30,aliasprefix_+"genjet1ptEta30");
    iEvent.put(genjet2ptEta30,aliasprefix_+"genjet2ptEta30");
    iEvent.put(genjet3ptEta30,aliasprefix_+"genjet3ptEta30");
    iEvent.put(genjet4ptEta30,aliasprefix_+"genjet4ptEta30");
    iEvent.put(genDphiJet12Eta30,aliasprefix_+"genDphiJet12Eta30");

    iEvent.put(ngenjetsPt40Eta30,aliasprefix_+"ngenjetsPt40Eta30");

    iEvent.put(genjet1ptEta22,aliasprefix_+"genjet1ptEta22");
    iEvent.put(genjet2ptEta22,aliasprefix_+"genjet2ptEta22");
    iEvent.put(genjet3ptEta22,aliasprefix_+"genjet3ptEta22");
    iEvent.put(genjet4ptEta22,aliasprefix_+"genjet4ptEta22");
    iEvent.put(genDphiJet12Eta22,aliasprefix_+"genDphiJet12Eta22");

    iEvent.put(ngenjetsPt40Eta22,aliasprefix_+"ngenjetsPt40Eta22");

    iEvent.put(genhtEta30,aliasprefix_+"genhtEta30");
    iEvent.put(genhtEta22,aliasprefix_+"genhtEta22");

    iEvent.put(genmet,aliasprefix_+"genmet");
    iEvent.put(genmetEta30,aliasprefix_+"genmetEta30");

    iEvent.put(genmhtEta30,aliasprefix_+"genmhtEta30");
    iEvent.put(genmhtEta22,aliasprefix_+"genmhtEta22");

    iEvent.put(genmt2Eta30,aliasprefix_+"genmt2Eta30");
    iEvent.put(genmt2Eta22,aliasprefix_+"genmt2Eta22");
  }


}


//____________________________________________________________
Int_t L1CustomNtupleProducer::phiINjetCoord(Double_t phi) {

  size_t phiIdx = 0;
  Double_t phidegree = degree(phi);
  for (size_t idx=0; idx<PHIBINS; idx++) {
    if (phidegree>=PHIBIN[idx] and phidegree<PHIBIN[idx+1])
      phiIdx = idx;
    else if (phidegree>=PHIBIN[PHIBINS-1] || phidegree<=PHIBIN[0])
      phiIdx = idx;
  }
  phiIdx = phiIdx + 1;
  if (phiIdx == 18)  phiIdx = 0;

  return int(phiIdx);
}

//____________________________________________________________
Int_t L1CustomNtupleProducer::etaINjetCoord(Double_t eta) {

  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETABINS; idx++) {
    if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
      etaIdx = idx;
  }

  return int(etaIdx);
}

//____________________________________________________________
inline Double_t L1CustomNtupleProducer::degree(Double_t radian) {
  if (radian<0)
    return 360.+(radian/TMath::Pi()*180.);
  else
    return radian/TMath::Pi()*180.;
}

//____________________________________________________________
float L1CustomNtupleProducer::dphi_normal(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  else if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();

  return dphi;
}

//____________________________________________________________
float L1CustomNtupleProducer::dphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  //  if (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  if (dphi < -3.0) dphi += 2*TMath::Pi();
  else if (dphi > 3.2) dphi -= 2*TMath::Pi();
  //  std::cout << "dphi calc: phi1: " << phi1 << ", phi2: " << phi2 << ", diff: " << phi1-phi2 << ", dphi: " << dphi << std::endl;

  return dphi;
}

//____________________________________________________________
float L1CustomNtupleProducer::deltaR_normal(float eta1, float phi1, float eta2, float phi2) {
  float diffeta = eta1 - eta2;
  float diffphi = dphi_normal(phi1,phi2);
  return sqrt(diffeta*diffeta + diffphi*diffphi);
}

//____________________________________________________________
// calculate hemispheres and input them to MT2 calc
Float_t L1CustomNtupleProducer::CalcHemisphereAndMT2(float testmass, bool massive, std::vector<LorentzVector> jets, LorentzVector MET ) {

  if (doDebug) std::cout << "-- debug for MT2 calculation:" << std::endl;

  // fill Pseudojets with selected objects
  vector<float> px, py, pz, E;
  for(unsigned int ijet=0; ijet<jets.size(); ++ijet){
    px.push_back(jets[ijet].Px());
    py.push_back(jets[ijet].Py());
    pz.push_back(jets[ijet].Pz());
    E .push_back(jets[ijet].E ());
    if (doDebug) std::cout << "  - jet: pt: " << jets[ijet].pt() << ", phi: " << jets[ijet].phi() << std::endl;
  }
  if (px.size()<2) return -999.99;

  // get hemispheres (seed 2: max inv mass, association method: default 3 = minimal lund distance)
  const int hemi_seed = 2;
  const int hemi_association = 3;
  //  Hemisphere* hemisp = new Hemisphere(px, py, pz, E, hemi_seed, hemi_association);
  //  vector<int> grouping = hemisp->getGrouping();
  Hemisphere hemisp(px, py, pz, E, hemi_seed, hemi_association);
  vector<int> grouping = hemisp.getGrouping();

  LorentzVector pseudojet1(0.,0.,0.,0.);
  LorentzVector pseudojet2(0.,0.,0.,0.);

  for(unsigned int i=0; i<px.size(); ++i){
    if(grouping[i]==1){
        pseudojet1.SetPx(pseudojet1.Px() + px[i]);
        pseudojet1.SetPy(pseudojet1.Py() + py[i]);
        pseudojet1.SetPz(pseudojet1.Pz() + pz[i]);
        pseudojet1.SetE( pseudojet1.E()  + E[i]);
    }else if(grouping[i] == 2){
        pseudojet2.SetPx(pseudojet2.Px() + px[i]);
        pseudojet2.SetPy(pseudojet2.Py() + py[i]);
        pseudojet2.SetPz(pseudojet2.Pz() + pz[i]);
        pseudojet2.SetE( pseudojet2.E()  + E[i]);
    }
  }
  //  delete hemisp;

  if (doDebug) std::cout << "  - pseudojet1: pt: " << pseudojet1.pt() << ", phi: " << pseudojet1.phi() << std::endl;
  if (doDebug) std::cout << "  - pseudojet2: pt: " << pseudojet2.pt() << ", phi: " << pseudojet2.phi() << std::endl;

  return CalcMT2(testmass, massive, pseudojet1, pseudojet2, MET);

}


//____________________________________________________________
// calculate MT2
Float_t L1CustomNtupleProducer::CalcMT2(float testmass, bool massive, LorentzVector visible1, LorentzVector visible2, LorentzVector MET ){

  double pa[3];
  double pb[3];
  double pmiss[3];

  pmiss[0] = 0;
  pmiss[1] = static_cast<double> (MET.Px());
  pmiss[2] = static_cast<double> (MET.Py());

  pa[0] = static_cast<double> (massive ? visible1.M() : 0);
  pa[1] = static_cast<double> (visible1.Px());
  pa[2] = static_cast<double> (visible1.Py());

  pb[0] = static_cast<double> (massive ? visible2.M() : 0);
  pb[1] = static_cast<double> (visible2.Px());
  pb[2] = static_cast<double> (visible2.Py());

  Davismt2 mt2;
  mt2.set_momenta(pa, pb, pmiss);
  mt2.set_mn(testmass);
  Float_t MT2=mt2.get_mt2();

  if (doDebug) std::cout << "  - MT2 value: " << MT2 << std::endl;
  return MT2;
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1CustomNtupleProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1CustomNtupleProducer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
L1TrackTriggerObjectsAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
L1CustomNtupleProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1CustomNtupleProducer);
