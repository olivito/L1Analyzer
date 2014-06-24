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
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"


// Gen-level stuff:
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveDigi.h"
#include "DataFormats/HcalDigi/interface/HcalTriggerPrimitiveDigi.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"

#include "TFile.h"
#include "TMath.h"
#include "TH1D.h"
#include "TH2D.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

using namespace l1extra;

using namespace reco;
using namespace std;


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
  bool operator() (const math::XYZTLorentzVector &r1, const math::XYZTLorentzVector &r2) {
    return r1.pt()>r2.pt();
  }
};

/////  

class L1RateAnalyzer : public edm::EDAnalyzer {
public:

  explicit L1RateAnalyzer(const edm::ParameterSet&);
  ~L1RateAnalyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  virtual void plot1D(std::string title, float xval, double weight, std::map<std::string, TH1D*> &allhistos, int numbinsx, float xmin, float xmax);
  virtual void plot2D(std::string title, float xval, float yval, std::string titleX, std::string titleY , double weight, std::map<std::string, TH2D*> &allhistos, int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax);

  void plotHT(float ht, float npu, std::string var, double weight, int nbins, float max);

  Int_t phiINjetCoord(Double_t phi);    
  Int_t etaINjetCoord(Double_t eta);
  inline Double_t degree(Double_t radian);
  float dphi(float phi1, float phi2);
  float deltaR(float eta1, float phi1, float eta2, float phi2);
  double regionPhysicalEt(const L1CaloRegion& cand) const {
    double regionLSB_=0.5;
    //    double regionLSB_=1.0;
    return std::max(0.,regionLSB_*cand.et());
  }


  // ----------member data ---------------------------
                                                                                                                                          
  std::map<std::string, TH1D*> h_1d;
  std::map<std::string, TH2D*> h_2d;

  std::string histoname_;
  edm::InputTag EcalTriggerPrimitiveInputTag;
  edm::InputTag HcalTriggerPrimitiveInputTag;
  edm::InputTag L1JetsCentInputTag;
  edm::InputTag L1JetsFwdInputTag;
  edm::InputTag L1MHTInputTag;
  //
  edm::InputTag L1JetsCent2015InputTag;
  edm::InputTag L1JetsFwd2015InputTag;
  edm::InputTag L1MHT2015InputTag;
  edm::InputTag Regions2015InputTag;
  edm::InputTag CorRegions2015InputTag;

  edm::InputTag PUSummaryInfoInputTag;

  int nBinHt=2000;
  double maxHt=2000;

  int nBinJetPt=600;
  double maxJetPt=600;

  int nBinPU = 60;

  unsigned int regionETCutForHT_;
  unsigned int minGctEtaForSums_;
  unsigned int maxGctEtaForSums_;

  bool doDebug=false; 
  bool doPUPlots_;
  bool doTPPlots_;
  bool applyGenWeights_;

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
L1RateAnalyzer::L1RateAnalyzer(const edm::ParameterSet& iConfig)

{
  EcalTriggerPrimitiveInputTag= iConfig.getParameter<edm::InputTag>("EcalTriggerPrimitiveInputTag");
  HcalTriggerPrimitiveInputTag= iConfig.getParameter<edm::InputTag>("HcalTriggerPrimitiveInputTag");
  L1JetsCentInputTag = iConfig.getParameter<edm::InputTag>("L1JetsCentInputTag");
  L1JetsFwdInputTag = iConfig.getParameter<edm::InputTag>("L1JetsFwdInputTag");
  L1MHTInputTag = iConfig.getParameter<edm::InputTag>("L1MHTInputTag");
  //
  L1JetsCent2015InputTag = iConfig.getParameter<edm::InputTag>("L1JetsCent2015InputTag");
  L1JetsFwd2015InputTag = iConfig.getParameter<edm::InputTag>("L1JetsFwd2015InputTag");
  L1MHT2015InputTag= iConfig.getParameter<edm::InputTag>("L1MHT2015InputTag");
  Regions2015InputTag= iConfig.getParameter<edm::InputTag>("Regions2015InputTag");
  CorRegions2015InputTag= iConfig.getParameter<edm::InputTag>("CorRegions2015InputTag");
  //
  PUSummaryInfoInputTag = iConfig.getParameter<edm::InputTag>("PUSummaryInfoInputTag");

  regionETCutForHT_ = iConfig.getParameter<unsigned int>("regionETCutForHT");
  minGctEtaForSums_ = iConfig.getParameter<unsigned int>("minGctEtaForSums");
  maxGctEtaForSums_ = iConfig.getParameter<unsigned int>("maxGctEtaForSums");

  //now do what ever initialization is needed
  histoname_= iConfig.getParameter<std::string>("histoname");  

  doTPPlots_= iConfig.getParameter<bool>("doTPPlots");  
  doPUPlots_= iConfig.getParameter<bool>("doPUPlots");  
  applyGenWeights_= iConfig.getParameter<bool>("applyGenWeights"); 

}


L1RateAnalyzer::~L1RateAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void L1RateAnalyzer::plot1D(std::string title, float xval, double weight, std::map<std::string, TH1D*> &allhistos,
                                    int numbinsx, float xmin, float xmax)
{

  std::map<std::string, TH1D*>::iterator iter= allhistos.find(title);
  if(iter == allhistos.end()) //no histo for this yet, so make a new one 
                                                                                                                     
    {
      TH1D* currentHisto= new TH1D(title.c_str(), title.c_str(), numbinsx, xmin, xmax);
      currentHisto->Sumw2();
      currentHisto->Fill(xval, weight);
      allhistos.insert(std::pair<std::string, TH1D*> (title,currentHisto) );
    }
  else // exists already, so just fill it

    {
      (*iter).second->Fill(xval, weight);
    }
}

void L1RateAnalyzer::plot2D(std::string title, float xval, float yval, std::string titleX, std::string titleY , double weight, std::map<std::string, TH2D*> &allhistos, int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax)
{


  std::map<std::string, TH2D*>::iterator iter= allhistos.find(title);
  if(iter == allhistos.end()) //no histo for this yet, so make a new one 
                                                                                                                     
    {
      TH2D* currentHisto= new TH2D(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax);
      currentHisto->Fill(xval, yval, weight);
      currentHisto->GetXaxis()->SetTitle(titleX.c_str());
      currentHisto->GetYaxis()->SetTitle(titleY.c_str());
      allhistos.insert(std::pair<std::string, TH2D*> (title,currentHisto) );
    }
  else // exists already, so just fill it
                                                                                                                                                     
    {
      (*iter).second->Fill(xval, yval, weight);
      //      (*iter).second->SetTitleX(titleX);
      //      (*iter).second->SetTitleY(titleY);
    }

  return;

}


//
// member functions
//

// ------------ method called for each event  ------------
void
L1RateAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  

  if(doDebug) std::cout << " ----  a new event ----- " << std::endl;

  double weight = 1.;

  // gen level weights
  if (applyGenWeights_) {
    Handle<GenEventInfoProduct> evt_info;
    iEvent.getByLabel("generator", evt_info);
    if(weight>0) weight = evt_info->weight();       
  }
  
  //
  // ----------------------------------------------------------------------
  // retrieve the PU summary info
  //

  edm::Handle<std::vector<PileupSummaryInfo> > PUSummaryInfoHandle;
  iEvent.getByLabel(PUSummaryInfoInputTag, PUSummaryInfoHandle);

  float nPU = -1;
  if (PUSummaryInfoHandle.isValid()) {
    nPU = PUSummaryInfoHandle->at(0).getTrueNumInteractions();
  } else {
    std::cout << "WARNING: no PU info found! tag was: " << PUSummaryInfoInputTag.encode() << std::endl;
  }

  plot1D("h_nPU",nPU,weight,h_1d, nBinPU,0.,nBinPU);

  //
  // ----------------------------------------------------------------------
  // retrieve the jets objects
  //

  edm::Handle<L1JetParticleCollection> L1JetsCentHandle;
  iEvent.getByLabel(L1JetsCentInputTag, L1JetsCentHandle);
  edm::Handle<L1JetParticleCollection> L1JetsFwdHandle;
  iEvent.getByLabel(L1JetsFwdInputTag, L1JetsFwdHandle);

  // std::vector<math::XYZTLorentzVector> newJets;
  // std::vector<math::XYZTLorentzVector> newJetsCent30;
  // std::vector<math::XYZTLorentzVector> newJetsCent22;

  if ( L1JetsCentHandle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsCentHandle -> begin(); jetIter != L1JetsCentHandle->end(); ++jetIter) {
      float et = jetIter -> pt();
      float eta = jetIter -> eta();
      float phi = jetIter -> phi();

      plot1D("h_jets2012_pt", et, weight, h_1d, nBinJetPt, 0., maxJetPt);
      plot1D("h_jets2012_eta", eta, weight, h_1d, 100, -5., 5.);
      plot1D("h_jets2012_phi", phi, weight, h_1d, 100, -3.5,3.5);

      int etabin = etaINjetCoord(eta);
      int phibin = phiINjetCoord(phi);

      plot1D("h_jets2012_etabin", etabin, weight, h_1d, ETABINS, 0., ETABINS);
      plot1D("h_jets2012_phibin", phibin, weight, h_1d, PHIBINS, 0., PHIBINS);

      // newJets.push_back(jetIter->p4());
      // if (fabs(eta) < 3.0) newJetsCent30.push_back(jetIter->p4());
      // if (fabs(eta) < 2.172) newJetsCent22.push_back(jetIter->p4());
    }
  }

  // loop also over forward jets
  if ( L1JetsFwdHandle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsFwdHandle -> begin(); jetIter != L1JetsFwdHandle->end(); ++jetIter) {
      float et = jetIter -> pt();
      float eta = jetIter -> eta();
      float phi = jetIter -> phi();

      plot1D("h_jets2012_fwd_pt", et, weight, h_1d, nBinJetPt, 0., maxJetPt);
      plot1D("h_jets2012_fwd_eta", eta, weight, h_1d, 100, -5., 5.);
      plot1D("h_jets2012_fwd_phi", phi, weight, h_1d, 100, -3.5,3.5);

      int etabin = etaINjetCoord(eta);
      int phibin = phiINjetCoord(phi);

      plot1D("h_jets2012_fwd_etabin", etabin, weight, h_1d, ETABINS, 0., ETABINS);
      plot1D("h_jets2012_fwd_phibin", phibin, weight, h_1d, PHIBINS, 0., PHIBINS);

      // combine with central jets
      //      newJets.push_back(jetIter->p4());
    }
  }

  //
  // ----------------------------------------------------------------------
  // retrieve the jets objects
  //

  edm::Handle<L1JetParticleCollection> L1JetsCent2015Handle;
  iEvent.getByLabel(L1JetsCent2015InputTag, L1JetsCent2015Handle);
  edm::Handle<L1JetParticleCollection> L1JetsFwd2015Handle;
  iEvent.getByLabel(L1JetsFwd2015InputTag, L1JetsFwd2015Handle);

  std::vector<math::XYZTLorentzVector> newJets2015;
  std::vector<math::XYZTLorentzVector> newJets2015Cent30;
  std::vector<math::XYZTLorentzVector> newJets2015Cent22;

  float HTjets = 0.;
  float HTjets_eta30 = 0.;
  float HTjets_eta22 = 0.;

  if ( L1JetsCent2015Handle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsCent2015Handle -> begin(); jetIter != L1JetsCent2015Handle->end(); ++jetIter) {
      float et = jetIter -> pt();
      float eta = jetIter -> eta();
      float phi = jetIter -> phi();

      plot1D("h_jets_pt", et, weight, h_1d, nBinJetPt, 0., maxJetPt);
      plot1D("h_jets_eta", eta, weight, h_1d, 100, -5., 5.);
      plot1D("h_jets_phi", phi, weight, h_1d, 100, -3.5,3.5);

      int etabin = etaINjetCoord(eta);
      int phibin = phiINjetCoord(phi);

      plot1D("h_jets_etabin", etabin, weight, h_1d, ETABINS, 0., ETABINS);
      plot1D("h_jets_phibin", phibin, weight, h_1d, PHIBINS, 0., PHIBINS);

      bool doOverlapCheck = true;
      if (doOverlapCheck) {
	float mindR = 99.;
	for (std::vector<L1JetParticle>::const_iterator jetIter2 = L1JetsCent2015Handle -> begin(); jetIter2 != L1JetsCent2015Handle->end(); ++jetIter2) {
	  if (jetIter2 == jetIter) continue;
	  float eta2 = jetIter2 -> eta();
	  float phi2 = jetIter2 -> phi();
	  float dR = deltaR(eta, phi, eta2, phi2);
	  if (dR < mindR) mindR = dR;
	}
	plot1D("h_jets_mindR", mindR, weight, h_1d, 100, 0., 6.5);

	// for (unsigned int ijet = 0; ijet < newJets2015.size(); ++ijet) {
	//   if ( (fabs(eta - newJets2015.at(ijet).eta()) < 0.1) && (fabs(dphi(phi,newJets2015.at(ijet).phi())) < 0.1) ) {
	//     std::cout << "WARNING: jet overlap.  jet1 pt: " << newJets2015.at(ijet).pt() 
	// 	      << ", eta: " << newJets2015.at(ijet).eta() << ", phi: " << newJets2015.at(ijet).phi()
	// 	      << ",   jet2 pt: " << et << ", eta: " << eta << ", phi: " << phi << std::endl;
	//   } 
	// }
      } // doOverlapCheck

      newJets2015.push_back(jetIter->p4());
      HTjets += et;
      if (fabs(eta) < 3.0) {
	newJets2015Cent30.push_back(jetIter->p4());
	HTjets_eta30 += et;
      }
      if (fabs(eta) < 2.172) {
	newJets2015Cent22.push_back(jetIter->p4());
	HTjets_eta22 += et;
      }

    }
  }

  // loop also over forward jets
  if ( L1JetsFwd2015Handle.isValid() ) {
    for (std::vector<L1JetParticle>::const_iterator jetIter = L1JetsFwd2015Handle -> begin(); jetIter != L1JetsFwd2015Handle->end(); ++jetIter) {
      float et = jetIter -> pt();
      float eta = jetIter -> eta();
      float phi = jetIter -> phi();

      plot1D("h_jets_fwd_pt", et, weight, h_1d, nBinJetPt, 0., maxJetPt);
      plot1D("h_jets_fwd_eta", eta, weight, h_1d, 100, -5., 5.);
      plot1D("h_jets_fwd_phi", phi, weight, h_1d, 100, -3.5,3.5);

      int etabin = etaINjetCoord(eta);
      int phibin = phiINjetCoord(phi);

      plot1D("h_jets_fwd_etabin", etabin, weight, h_1d, ETABINS, 0., ETABINS);
      plot1D("h_jets_fwd_phibin", phibin, weight, h_1d, PHIBINS, 0., PHIBINS);

      // combine with central jets
      newJets2015.push_back(jetIter->p4());
      HTjets += et;
    }
  }

  sort(newJets2015.begin(),newJets2015.end(),cmpPtLorentz());
  sort(newJets2015Cent30.begin(),newJets2015Cent30.end(),cmpPtLorentz());
  sort(newJets2015Cent22.begin(),newJets2015Cent22.end(),cmpPtLorentz());

  float pt_jet1 = 0.;
  float pt_jet2 = 0.;
  float pt_jet3 = 0.;
  float pt_jet4 = 0.;
  float dphi_jet12 = -99.;

  if (newJets2015.size() > 0) {
    pt_jet1 = newJets2015.at(0).pt();
  }
  if (newJets2015.size() > 1) {
    pt_jet2 = newJets2015.at(1).pt();
    dphi_jet12 = dphi(newJets2015.at(0).phi(),newJets2015.at(1).phi());
  }
  if (newJets2015.size() > 2) {
    pt_jet3 = newJets2015.at(2).pt();
  }
  if (newJets2015.size() > 3) {
    pt_jet4 = newJets2015.at(3).pt();
  }

  plot1D("h_jet1_pt", pt_jet1, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet2_pt", pt_jet2, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet3_pt", pt_jet3, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet4_pt", pt_jet4, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet12_dphi", dphi_jet12, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_jet12_dphi_vs_jet2_pt", pt_jet2, dphi_jet12, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  float pt_jet1_eta30 = 0.;
  float pt_jet2_eta30 = 0.;
  float pt_jet3_eta30 = 0.;
  float pt_jet4_eta30 = 0.;
  float dphi_jet12_eta30 = -99.;

  if (newJets2015Cent30.size() > 0) {
    pt_jet1_eta30 = newJets2015Cent30.at(0).pt();
  }
  if (newJets2015Cent30.size() > 1) {
    pt_jet2_eta30 = newJets2015Cent30.at(1).pt();
    dphi_jet12_eta30 = dphi(newJets2015Cent30.at(0).phi(),newJets2015Cent30.at(1).phi());
  }
  if (newJets2015Cent30.size() > 2) {
    pt_jet3_eta30 = newJets2015Cent30.at(2).pt();
  }
  if (newJets2015Cent30.size() > 3) {
    pt_jet4_eta30 = newJets2015Cent30.at(3).pt();
  }

  plot1D("h_jet1_pt_eta30", pt_jet1_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet2_pt_eta30", pt_jet2_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet3_pt_eta30", pt_jet3_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet4_pt_eta30", pt_jet4_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet12_dphi_eta30", dphi_jet12_eta30, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_jet12_dphi_vs_jet2_pt_eta30", pt_jet2_eta30, dphi_jet12_eta30, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  float pt_jet1_eta22 = 0.;
  float pt_jet2_eta22 = 0.;
  float pt_jet3_eta22 = 0.;
  float pt_jet4_eta22 = 0.;
  float dphi_jet12_eta22 = -99.;

  if (newJets2015Cent22.size() > 0) {
    pt_jet1_eta22 = newJets2015Cent22.at(0).pt();
  }
  if (newJets2015Cent22.size() > 1) {
    pt_jet2_eta22 = newJets2015Cent22.at(1).pt();
    dphi_jet12_eta22 = dphi(newJets2015Cent22.at(0).phi(),newJets2015Cent22.at(1).phi());
  }
  if (newJets2015Cent22.size() > 2) {
    pt_jet3_eta22 = newJets2015Cent22.at(2).pt();
  }
  if (newJets2015Cent22.size() > 3) {
    pt_jet4_eta22 = newJets2015Cent22.at(3).pt();
  }

  plot1D("h_jet1_pt_eta22", pt_jet1_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet2_pt_eta22", pt_jet2_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet3_pt_eta22", pt_jet3_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet4_pt_eta22", pt_jet4_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_jet12_dphi_eta22", dphi_jet12_eta22, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_jet12_dphi_vs_jet2_pt_eta22", pt_jet2_eta22, dphi_jet12_eta22, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  // --- HT with jets

  plotHT(HTjets,nPU,"HTjets",weight,nBinHt,maxHt);
  plotHT(HTjets_eta30,nPU,"HTjets_eta30",weight,nBinHt,maxHt);
  plotHT(HTjets_eta22,nPU,"HTjets_eta22",weight,nBinHt,maxHt);

  plot2D("h_HTjets_vs_jet1pt",pt_jet1,HTjets,"Leading Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_vs_jet2pt",pt_jet2,HTjets,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_vs_jet3pt",pt_jet3,HTjets,"3rd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_vs_jet4pt",pt_jet4,HTjets,"4th Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HTjets_eta30_vs_jet1pt_eta30",pt_jet1_eta30,HTjets_eta30,"Leading Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta30_vs_jet2pt_eta30",pt_jet2_eta30,HTjets_eta30,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta30_vs_jet3pt_eta30",pt_jet3_eta30,HTjets_eta30,"3rd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta30_vs_jet4pt_eta30",pt_jet4_eta30,HTjets_eta30,"4th Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HTjets_eta22_vs_jet1pt_eta22",pt_jet1_eta22,HTjets_eta22,"Leading Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta22_vs_jet2pt_eta22",pt_jet2_eta22,HTjets_eta22,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta22_vs_jet3pt_eta22",pt_jet3_eta22,HTjets_eta22,"3rd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HTjets_eta22_vs_jet4pt_eta22",pt_jet4_eta22,HTjets_eta22,"4th Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  // tight dphi cut: excludes back-to-back and one more bin in phi (dphi < 140-160 degrees)
  if (fabs(dphi_jet12_eta30) < 2.6) {
    plot2D("h_HTjets_eta30_vs_jet2pt_eta30_dphiT",pt_jet2_eta30,HTjets_eta30,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // loose dphi cut: excludes back-to-back bin (dphi < 160-180 degrees)
  if (fabs(dphi_jet12_eta30) < 3.0) {
    plot2D("h_HTjets_eta30_vs_jet2pt_eta30_dphiL",pt_jet2_eta30,HTjets_eta30,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // tight dphi cut: excludes back-to-back and one more bin in phi (dphi < 140-160 degrees)
  if (fabs(dphi_jet12_eta22) < 2.6) {
    plot2D("h_HTjets_eta22_vs_jet2pt_eta22_dphiT",pt_jet2_eta22,HTjets_eta22,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // loose dphi cut: excludes back-to-back bin (dphi < 160-180 degrees)
  if (fabs(dphi_jet12_eta22) < 3.0) {
    plot2D("h_HTjets_eta22_vs_jet2pt_eta22_dphiL",pt_jet2_eta22,HTjets_eta22,"2nd Jet PT [GeV]","HT from jets [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }


  //
  // ----------------------------------------------------------------------
  // retrieve the caloMHT objects
  //

  edm::Handle<L1EtMissParticleCollection> L1MHTHandle;
  iEvent.getByLabel(L1MHTInputTag, L1MHTHandle);
  
  double caloMHt_=0;
  double caloHt_=0;
  
  if (L1MHTHandle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator mhtIter = L1MHTHandle->begin(); mhtIter != L1MHTHandle->end(); ++mhtIter) {
      caloMHt_ = mhtIter -> etMiss();
      caloHt_ = mhtIter -> etTotal();

      //      cout << "overflowHT " << mhtIter->gctEtMiss()->overFlow() << endl;
      //      if(mhtIter->gctEtMiss()->overFlow()>0) cout << "!!!!!!!!!!!!! HT BAD OVERFLOW "<< endl;
    }
  }
  
  // plot1D("h_caloMHT",caloMHt_,1,h_1d, nBinMet,0.,maxMet);
  plot1D("h_caloHT2012",caloHt_,weight,h_1d, nBinHt,0.,maxHt);

  //
  // ----------------------------------------------------------------------
  // retrieve the caloMHT objects 2015
  //

  edm::Handle<L1EtMissParticleCollection> L1MHT2015Handle;
  iEvent.getByLabel(L1MHT2015InputTag, L1MHT2015Handle);
  
  double caloMHt2015_=0;
  double caloHt2015_=0;
  
  if (L1MHT2015Handle.isValid() ) {
    for (L1EtMissParticleCollection::const_iterator mhtIter = L1MHT2015Handle->begin(); mhtIter != L1MHT2015Handle->end(); ++mhtIter) {
      caloMHt2015_ = mhtIter -> etMiss();
      caloHt2015_ = mhtIter -> etTotal();

      //      cout << "overflowHT " << mhtIter->gctEtMiss()->overFlow() << endl;
//      if(mhtIter->gctEtMiss()->overFlow()>0) cout << "!!!!!!!!!!!!! HT2015 BAD OVERFLOW "<< endl;

      //      cout << "overflowHT2015 " << mhtIter->gctEtTotalRef()->overFlow() << endl;
    }
  }
  
  plotHT(caloHt_,nPU,"caloHT2012",weight,nBinHt,maxHt);
  plotHT(caloHt2015_,nPU,"caloHT2015",weight,nBinHt,maxHt);

  plot2D("h_HTcorr",caloHt_,caloHt2015_,"2012 HT [GeV]","2015 HT [GeV]",weight,h_2d,nBinHt,0.,maxHt,nBinHt,0.,maxHt);

  plot2D("h_HT_vs_MHT2012",caloMHt_,caloHt_,"2012 MHT [GeV]","2012 HT [GeV]",weight,h_2d,300,0.,3000.,50,0.,500.);
  plot2D("h_HT_vs_MHT2015",caloMHt2015_,caloHt2015_,"MHT [GeV]","HT [GeV]",weight,h_2d,300,0.,3000.,50,0.,500.);

  //
  // ----------------------------------------------------------------------
  // input to the sums
  //
  Handle<L1CaloRegionCollection> newRegions;
  iEvent.getByLabel(Regions2015InputTag, newRegions);

  /*
  if(puMultCorrect) {
    edm::Handle<int> puweightHandle;
    iEvent.getByLabel("CorrectedDigis","CorrectedRegions", newRegions);
    iEvent.getByLabel("CorrectedDigis","PUM0Level",puweightHandle);
    puLevelPUM0=(*puweightHandle);
  }
  else {iEvent.getByLabel("uctDigis", newRegions);}
  */

  double rawHT = 0;

  float regionET1 = 0.;
  float regionET2 = 0.;
  float regionET3 = 0.;

  for(L1CaloRegionCollection::const_iterator newRegion = newRegions->begin();
      newRegion != newRegions->end(); newRegion++){
    // Remove forward stuff 
    if (newRegion->gctEta() < minGctEtaForSums_ || newRegion->gctEta() > maxGctEtaForSums_) {
      continue;
    }

    double regionET =  regionPhysicalEt(*newRegion);
    if (regionET > regionET1) {
      regionET3 = regionET2;
      regionET2 = regionET1;
      regionET1 = regionET;
    } else if (regionET > regionET2) {
      regionET3 = regionET2;
      regionET2 = regionET;
    } else if (regionET > regionET3) {
      regionET3 = regionET;
    }

    if(regionET==256) cout << "found saturation 256 " << endl;
    if(regionET==512) cout << "found saturation 512 " << endl;

    plot1D("h_regionET", regionET ,weight,h_1d, 1000,0.,1000);
    if (newRegion->gctEta() >6 || newRegion->gctEta() < 15) plot1D("h_regionET_barrel", regionET ,weight,h_1d, 1000,0.,1000);
    if (newRegion->gctEta() == 4 || 
	newRegion->gctEta() == 5 ||
	newRegion->gctEta() == 6 ||
	newRegion->gctEta() == 15 ||
	newRegion->gctEta() == 16 ||
	newRegion->gctEta() == 17 
	) plot1D("h_regionET_endcap", regionET ,weight,h_1d, 1000,0.,1000);

    if (regionET >= double(regionETCutForHT_)) {
      rawHT += regionET;
    }
  }

  plot2D("h_regionET1_vs_HT2015",caloHt2015_,regionET1,"HT [GeV]","Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
  plot2D("h_regionET2_vs_HT2015",caloHt2015_,regionET2,"HT [GeV]","2nd Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
  plot2D("h_regionET3_vs_HT2015",caloHt2015_,regionET3,"HT [GeV]","3rd Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);

  plotHT(rawHT,nPU,"rawHT",weight,nBinHt,maxHt);
  plot2D("h_recompRawHTcorr",rawHT,caloHt2015_,"recomputed Raw HT [GeV]","2015 HT [GeV]",weight,h_2d,nBinHt,0.,maxHt,nBinHt,0.,maxHt);

  //
  // ----------------------------------------------------------------------
  // input to the sums, after PU correction
  //
  Handle<L1CaloRegionCollection> corRegions;
  iEvent.getByLabel(CorRegions2015InputTag, corRegions);

  double HT = 0.;
  double HT_reget0 = 0.;
  double HT_reget3 = 0.;
  double HT_reget7 = 0.;
  double HT_reget10 = 0.;
  double HT_reget15 = 0.;
  double HT_reget20 = 0.;
  double HT_alleta_reget7 = 0.;
  double HT_eta22_reget0 = 0.;
  double HT_eta22_reget3 = 0.;
  double HT_eta22_reget7 = 0.;
  double HT_eta22_reget10 = 0.;
  double HT_eta22_reget15 = 0.;
  double HT_eta22_reget20 = 0.;

  float corregionET1 = 0.;
  float corregionET2 = 0.;
  float corregionET3 = 0.;

  for(L1CaloRegionCollection::const_iterator corRegion = corRegions->begin();
      corRegion != corRegions->end(); corRegion++){

    double regionET =  regionPhysicalEt(*corRegion);

    if (regionET >= 7.) {
      HT_alleta_reget7 += regionET;
    }

    // eta range:
    // 0-21 corresponds to all eta
    // 4-17 corresponds to |eta| < 3.0 (default)
    // 5-16 corresponds to |eta| < 2.2
    if (corRegion->gctEta() >= 4 && corRegion->gctEta() <= 17) {
      HT_reget0 += regionET;
      if (regionET >= 3.) {
	HT_reget3 += regionET;
      }
      if (regionET >= 7.) {
	HT_reget7 += regionET;
      }
      if (regionET >= 10.) {
	HT_reget10 += regionET;
      }
      if (regionET >= 15.) {
	HT_reget15 += regionET;
      }
      if (regionET >= 20.) {
	HT_reget20 += regionET;
      }
    }

    if (corRegion->gctEta() >= 5 && corRegion->gctEta() <= 16) {
      HT_eta22_reget0 += regionET;
      if (regionET >= 3.) {
	HT_eta22_reget3 += regionET;
      }
      if (regionET >= 7.) {
	HT_eta22_reget7 += regionET;
      }
      if (regionET >= 10.) {
	HT_eta22_reget10 += regionET;
      }
      if (regionET >= 15.) {
	HT_eta22_reget15 += regionET;
      }
      if (regionET >= 20.) {
	HT_eta22_reget20 += regionET;
      }
    }

    // Remove forward stuff
    if (corRegion->gctEta() < minGctEtaForSums_ || corRegion->gctEta() > maxGctEtaForSums_) {
      continue;
    }

    if (regionET > corregionET1) {
      corregionET3 = corregionET2;
      corregionET2 = corregionET1;
      corregionET1 = regionET;
    } else if (regionET > corregionET2) {
      corregionET3 = corregionET2;
      corregionET2 = regionET;
    } else if (regionET > corregionET3) {
      corregionET3 = regionET;
    }

    if(regionET==256) cout << "found saturation 256 " << endl;
    if(regionET==512) cout << "found saturation 512 " << endl;

    plot1D("h_corregionET", regionET ,weight,h_1d, 1000,0.,1000);
    if (corRegion->gctEta() >6 || corRegion->gctEta() < 15) plot1D("h_corregionET_barrel", regionET ,weight,h_1d, 1000,0.,1000);
    if (corRegion->gctEta() == 4 || 
	corRegion->gctEta() == 5 ||
	corRegion->gctEta() == 6 ||
	corRegion->gctEta() == 15 ||
	corRegion->gctEta() == 16 ||
	corRegion->gctEta() == 17 
	) plot1D("h_corregionET_endcap", regionET ,weight,h_1d, 1000,0.,1000);

    if (regionET >= double(regionETCutForHT_)) {
      // default, configurable HT
      HT += regionET;
    }
  }

  plot2D("h_corregionET1_vs_HT2015",caloHt2015_,corregionET1,"HT [GeV]","Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
  plot2D("h_corregionET2_vs_HT2015",caloHt2015_,corregionET2,"HT [GeV]","2nd Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
  plot2D("h_corregionET3_vs_HT2015",caloHt2015_,corregionET3,"HT [GeV]","3rd Highest Region Et [GeV]",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);

  plot2D("h_recompHTcorr",HT,caloHt2015_,"recomputed HT [GeV]","2015 HT [GeV]",weight,h_2d,nBinHt,0.,maxHt,nBinHt,0.,maxHt);

  plotHT(HT,nPU,"HT",weight,nBinHt,maxHt);
  plotHT(HT_reget0,nPU,"HT_reget0",weight,nBinHt,maxHt);
  plotHT(HT_reget3,nPU,"HT_reget3",weight,nBinHt,maxHt);
  plotHT(HT_reget7,nPU,"HT_reget7",weight,nBinHt,maxHt);
  plotHT(HT_reget10,nPU,"HT_reget10",weight,nBinHt,maxHt);
  plotHT(HT_reget15,nPU,"HT_reget15",weight,nBinHt,maxHt);
  plotHT(HT_reget20,nPU,"HT_reget20",weight,nBinHt,maxHt);
  plotHT(HT_alleta_reget7,nPU,"HT_alleta_reget7",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget0,nPU,"HT_eta22_reget0",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget3,nPU,"HT_eta22_reget3",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget7,nPU,"HT_eta22_reget7",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget10,nPU,"HT_eta22_reget10",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget15,nPU,"HT_eta22_reget15",weight,nBinHt,maxHt);
  plotHT(HT_eta22_reget20,nPU,"HT_eta22_reget20",weight,nBinHt,maxHt);

  plot2D("h_HT_vs_jet1pt",pt_jet1,HT_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet2pt",pt_jet2,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet3pt",pt_jet3,HT_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet4pt",pt_jet4,HT_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_vs_jet1pt_eta30",pt_jet1_eta30,HT_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet2pt_eta30",pt_jet2_eta30,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet3pt_eta30",pt_jet3_eta30,HT_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet4pt_eta30",pt_jet4_eta30,HT_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_vs_jet1pt_eta22",pt_jet1_eta22,HT_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet2pt_eta22",pt_jet2_eta22,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet3pt_eta22",pt_jet3_eta22,HT_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_vs_jet4pt_eta22",pt_jet4_eta22,HT_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_reget0_vs_jet1pt_eta30",pt_jet1_eta30,HT_reget0,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget0_vs_jet2pt_eta30",pt_jet2_eta30,HT_reget0,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget0_vs_jet3pt_eta30",pt_jet3_eta30,HT_reget0,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget0_vs_jet4pt_eta30",pt_jet4_eta30,HT_reget0,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_reget3_vs_jet1pt_eta30",pt_jet1_eta30,HT_reget3,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget3_vs_jet2pt_eta30",pt_jet2_eta30,HT_reget3,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget3_vs_jet3pt_eta30",pt_jet3_eta30,HT_reget3,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_reget3_vs_jet4pt_eta30",pt_jet4_eta30,HT_reget3,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_eta22_vs_jet1pt",pt_jet1,HT_eta22_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet2pt",pt_jet2,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet3pt",pt_jet3,HT_eta22_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet4pt",pt_jet4,HT_eta22_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_eta22_vs_jet1pt_eta30",pt_jet1_eta30,HT_eta22_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet2pt_eta30",pt_jet2_eta30,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet3pt_eta30",pt_jet3_eta30,HT_eta22_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet4pt_eta30",pt_jet4_eta30,HT_eta22_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_eta22_vs_jet1pt_eta22",pt_jet1_eta22,HT_eta22_reget7,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet2pt_eta22",pt_jet2_eta22,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet3pt_eta22",pt_jet3_eta22,HT_eta22_reget7,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_vs_jet4pt_eta22",pt_jet4_eta22,HT_eta22_reget7,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_eta22_reget20_vs_jet1pt_eta30",pt_jet1_eta30,HT_eta22_reget20,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet2pt_eta30",pt_jet2_eta30,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet3pt_eta30",pt_jet3_eta30,HT_eta22_reget20,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet4pt_eta30",pt_jet4_eta30,HT_eta22_reget20,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  plot2D("h_HT_eta22_reget20_vs_jet1pt_eta22",pt_jet1_eta22,HT_eta22_reget20,"1st Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet2pt_eta22",pt_jet2_eta22,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet3pt_eta22",pt_jet3_eta22,HT_eta22_reget20,"3rd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  plot2D("h_HT_eta22_reget20_vs_jet4pt_eta22",pt_jet4_eta22,HT_eta22_reget20,"4th Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);

  // tight dphi cut: excludes back-to-back and one more bin in phi (dphi < 140-160 degrees)
  if (fabs(dphi_jet12_eta30) < 2.6) {
    plot2D("h_HT_vs_jet2pt_eta30_dphiT",pt_jet2_eta30,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_vs_jet2pt_eta30_dphiT",pt_jet2_eta30,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_reget20_vs_jet2pt_eta30_dphiT",pt_jet2_eta30,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // loose dphi cut: excludes back-to-back bin (dphi < 160-180 degrees)
  if (fabs(dphi_jet12_eta30) < 3.0) {
    plot2D("h_HT_vs_jet2pt_eta30_dphiL",pt_jet2_eta30,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_vs_jet2pt_eta30_dphiL",pt_jet2_eta30,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_reget20_vs_jet2pt_eta30_dphiL",pt_jet2_eta30,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // tight dphi cut: excludes back-to-back and one more bin in phi (dphi < 140-160 degrees)
  if (fabs(dphi_jet12_eta22) < 2.6) {
    plot2D("h_HT_vs_jet2pt_eta22_dphiT",pt_jet2_eta22,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_vs_jet2pt_eta22_dphiT",pt_jet2_eta22,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_reget20_vs_jet2pt_eta22_dphiT",pt_jet2_eta22,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  // loose dphi cut: excludes back-to-back bin (dphi < 160-180 degrees)
  if (fabs(dphi_jet12_eta22) < 3.0) {
    plot2D("h_HT_vs_jet2pt_eta22_dphiL",pt_jet2_eta22,HT_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_vs_jet2pt_eta22_dphiL",pt_jet2_eta22,HT_eta22_reget7,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
    plot2D("h_HT_eta22_reget20_vs_jet2pt_eta22_dphiL",pt_jet2_eta22,HT_eta22_reget20,"2nd Jet PT [GeV]","HT [GeV]",weight,h_2d,nBinJetPt,0.,maxJetPt,nBinHt,0.,maxHt);
  }

  //
  // ----------------------------------------------------------------------
  // trigger primitives
  //

  if (doTPPlots_) {
    edm::Handle<EcalTrigPrimDigiCollection> ecal;
    edm::Handle<HcalTrigPrimDigiCollection> hcal;

    iEvent.getByLabel(EcalTriggerPrimitiveInputTag, ecal);
    iEvent.getByLabel(HcalTriggerPrimitiveInputTag, hcal);

    int ecaltp1 = -1;
    int ecaltp2 = -1;
    int ecaltp3 = -1;

    for (unsigned int i=0;i<ecal.product()->size();i++) {
      EcalTriggerPrimitiveDigi d = (*(ecal.product()))[i];
      //      cout << "energy " << d.compressedEt() << endl;
      int energy = d.compressedEt();
      // for a more complicated definition see:
      // https://cmssdt.cern.ch/SDT/lxr/source/L1Trigger/RegionalCaloTrigger/src/L1RCT.cc?v=CMSSW_7_1_0_pre8#150
      plot1D("h_ecalTriggerPrimitive", energy ,weight,h_1d, 1000,0.,1000);

      if (energy > ecaltp1) {
	ecaltp3 = ecaltp2;
	ecaltp2 = ecaltp1;
	ecaltp1 = energy;
      } else if (energy > ecaltp2) {
	ecaltp3 = ecaltp2;
	ecaltp2 = energy;
      } else if (energy > ecaltp3) {
	ecaltp3 = energy;
      }
      //      if(d[0].compressedEt()==255) foundSatRegion=true;
    }
 
    plot2D("h_ecaltp1_vs_HT2015",caloHt2015_,ecaltp1,"HT [GeV]","Highest Et ECAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
    plot2D("h_ecaltp2_vs_HT2015",caloHt2015_,ecaltp2,"HT [GeV]","2nd Highest Et ECAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
    plot2D("h_ecaltp3_vs_HT2015",caloHt2015_,ecaltp3,"HT [GeV]","3rd Highest Et ECAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);

    int hcaltp1 = -1;
    int hcaltp2 = -1;
    int hcaltp3 = -1;

    for (unsigned int i=0;i<hcal.product()->size();i++) {
      HcalTriggerPrimitiveDigi d = (*(hcal.product()))[i];
      //      cout << "energy " << d.SOI_compressedEt() << endl;
      int energy = d.SOI_compressedEt();
      plot1D("h_hcalTriggerPrimitive",  energy,weight,h_1d, 1000,0.,1000);

      if (energy > hcaltp1) {
	hcaltp3 = hcaltp2;
	hcaltp2 = hcaltp1;
	hcaltp1 = energy;
      } else if (energy > hcaltp2) {
	hcaltp3 = hcaltp2;
	hcaltp2 = energy;
      } else if (energy > hcaltp3) {
	hcaltp3 = energy;
      }
    }

    plot2D("h_hcaltp1_vs_HT2015",caloHt2015_,hcaltp1,"HT [GeV]","Highest Et HCAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
    plot2D("h_hcaltp2_vs_HT2015",caloHt2015_,hcaltp2,"HT [GeV]","2nd Highest Et HCAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);
    plot2D("h_hcaltp3_vs_HT2015",caloHt2015_,hcaltp3,"HT [GeV]","3rd Highest Et HCAL TP",weight,h_2d,nBinHt,0.,maxHt,1000,0.,1000.);

  } // if doTPPlots

}

void L1RateAnalyzer::plotHT(float ht, float npu, std::string var, double weight, int nbins, float max) {

  plot1D("h_"+var,ht,weight,h_1d, nbins,0.,max);
  if (doPUPlots_) plot2D("h_"+var+"_vs_nPU",npu,ht,"nPU","HT [GeV]",weight,h_2d,nBinPU,0.,nBinPU,nBinHt,0.,maxHt);

  return;
}

Int_t L1RateAnalyzer::phiINjetCoord(Double_t phi) {

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

Int_t L1RateAnalyzer::etaINjetCoord(Double_t eta) {

  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETABINS; idx++) {
    if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
      etaIdx = idx;
  }

  return int(etaIdx);
}

inline Double_t L1RateAnalyzer::degree(Double_t radian) {
  if (radian<0)
    return 360.+(radian/TMath::Pi()*180.);
  else
    return radian/TMath::Pi()*180.;
}

float L1RateAnalyzer::dphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  //  if (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  if (dphi < -3.0) dphi += 2*TMath::Pi();
  else if (dphi > 3.2) dphi -= 2*TMath::Pi();
  //  std::cout << "dphi calc: phi1: " << phi1 << ", phi2: " << phi2 << ", diff: " << phi1-phi2 << ", dphi: " << dphi << std::endl;

  return dphi;
}

float L1RateAnalyzer::deltaR(float eta1, float phi1, float eta2, float phi2) {
  float diffeta = eta1 - eta2;
  float diffphi = dphi(phi1,phi2);
  return sqrt(diffeta*diffeta + diffphi*diffphi);
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1RateAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1RateAnalyzer::endJob() 
{

  //  TFile*fout = new TFile("ttbar_controlplots.root","RECREATE");
  //  TFile*fout = new TFile("T2tt_controlplots.root","RECREATE");
  TFile*fout = new TFile(histoname_.data(),"RECREATE");
  
  std::map<std::string, TH1D*>::iterator it1d;
  for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
    it1d->second->Write();
    delete it1d->second;
  }

  std::map<std::string, TH2D*>::iterator it2d;
  for(it2d=h_2d.begin(); it2d!=h_2d.end(); it2d++) {
    it2d->second->Write();
    delete it2d->second;
  }

  //  fout->mkdir("Varie");
  //  fout->cd("Varie");

  fout->Write();
  fout->Close();

  //  gROOT->cd();

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
L1RateAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1RateAnalyzer);
