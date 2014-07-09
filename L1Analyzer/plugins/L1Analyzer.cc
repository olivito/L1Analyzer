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

class L1Analyzer : public edm::EDAnalyzer {
public:

  explicit L1Analyzer(const edm::ParameterSet&);
  ~L1Analyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;

  void plot1D(std::string title, float xval, double weight, std::map<std::string, TH1F*> &allhistos,
	 int numbinsx, float xmin, float xmax);
  void plot2D(std::string title, float xval, float yval, std::string titleX, std::string titleY , double weight, std::map<std::string, TH2F*> &allhistos,
		      int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax);
    
  void makeTurnOn(double pt, string var, double l1pt, string l1var, std::vector<int> thresh, double weight, int nBin, int max);

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
                                                                                                                                          
  std::map<std::string, TH1F*> h_1d;
  std::map<std::string, TH2F*> h_2d;

  std::string histoname_;
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

  int nBinMet=500;
  double maxMet=500;

  int nBinJetPt=600;
  double maxJetPt=600;

  int nBinMT2=1000;
  double maxMT2=1000;

  int nBinNJets=10;
  double maxNJets=10;

  int nBinPU = 60;

  unsigned int regionETCutForHT_;
  unsigned int minGctEtaForSums_;
  unsigned int maxGctEtaForSums_;

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
L1Analyzer::L1Analyzer(const edm::ParameterSet& iConfig)

{
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

}


L1Analyzer::~L1Analyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


void L1Analyzer::plot1D(std::string title, float xval, double weight, std::map<std::string, TH1F*> &allhistos,
                                    int numbinsx, float xmin, float xmax)
{

  std::map<std::string, TH1F*>::iterator iter= allhistos.find(title);
  if(iter == allhistos.end()) //no histo for this yet, so make a new one 
                                                                                                                     
    {
      TH1F* currentHisto= new TH1F(title.c_str(), title.c_str(), numbinsx, xmin, xmax);
      currentHisto->Sumw2();
      currentHisto->Fill(xval, weight);
      allhistos.insert(std::pair<std::string, TH1F*> (title,currentHisto) );
    }
  else // exists already, so just fill it

    {
      (*iter).second->Fill(xval, weight);
    }
}

void L1Analyzer::plot2D(std::string title, float xval, float yval, std::string titleX, std::string titleY , double weight, std::map<std::string, TH2F*> &allhistos, int numbinsx, float xmin, float xmax, int numbinsy, float ymin, float ymax)
{


  std::map<std::string, TH2F*>::iterator iter= allhistos.find(title);
  if(iter == allhistos.end()) //no histo for this yet, so make a new one 
                                                                                                                     
    {
      TH2F* currentHisto= new TH2F(title.c_str(), title.c_str(), numbinsx, xmin, xmax, numbinsy, ymin, ymax);
      currentHisto->Fill(xval, yval, weight);
      currentHisto->GetXaxis()->SetTitle(titleX.c_str());
      currentHisto->GetYaxis()->SetTitle(titleY.c_str());
      allhistos.insert(std::pair<std::string, TH2F*> (title,currentHisto) );
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
L1Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;

  double weight = 1.0;  

  if(doDebug) std::cout << " ----  a new event ----- " << std::endl;
 
  // First, retrieve the generated primary vertex
  
  edm::Handle<edm::HepMCProduct> HepMCEvt;
  iEvent.getByLabel("generator",HepMCEvt);
  
  const HepMC::GenEvent* MCEvt = HepMCEvt->GetEvent();
  const double mm=0.1;
  
  float zvtx_gen = -999;
  for ( HepMC::GenEvent::vertex_const_iterator ivertex = MCEvt->vertices_begin(); ivertex != MCEvt->vertices_end(); ++ivertex )
    {
      bool hasParentVertex = false;
      
      // Loop over the parents looking to see if they are coming from a production vertex
      for (
	   HepMC::GenVertex::particle_iterator iparent = (*ivertex)->particles_begin(HepMC::parents);
	   iparent != (*ivertex)->particles_end(HepMC::parents);
	   ++iparent
	   )
	if ( (*iparent)->production_vertex() )
	  {
	    hasParentVertex = true;
	    break;
	  }
      
      // Reject those vertices with parent vertices
      if (hasParentVertex) continue;
      
             // Get the position of the vertex
      HepMC::FourVector pos = (*ivertex)->position();
      zvtx_gen = pos.z()*mm; 
      
      break;  // there should be one single primary vertex
      
    }  // end loop over gen vertices
  
  //  h_zgen -> Fill( zvtx_gen );
  if(doDebug) std::cout << " Generated zvertex : " << zvtx_gen << std::endl;

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
                                                                                        
  int nGenLep = 0;
  for (reco::GenParticleCollection::const_iterator genp = genparts->begin(); genp != genparts->end(); ++ genp ) {

    //    if(genp->status() == 23 || genp->status() == 22 ) cout << "genp->pdgId() " << genp->pdgId() << " genp->status() " << genp->status() << endl;

    //    if(abs(genp->pdgId()) == 15) cout << "genp->pdgId() " << genp->pdgId() << " genp->status() " << genp->status() << endl;

    //    if ((abs(genp->pdgId()) == 11 || abs(genp->pdgId()) == 12 || abs(genp->pdgId()) == 13 || abs(genp->pdgId()) == 14 || abs(genp->pdgId()) == 15 || abs(genp->pdgId()) == 16 ) && genp->status() == 23) nGenLep++;

    if ((genp->status() == 23) && (abs(genp->pdgId()) == 11 || abs(genp->pdgId()) == 13 || abs(genp->pdgId()) == 15)) nGenLep++;

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

  //  cout << "nGenLep " << nGenLep << endl;
  plot1D("h_nGenLep",nGenLep,weight,h_1d, 3,0.,3.);

  for (unsigned int ilsp = 0; ilsp < lsps.size(); ++ilsp) {
    plot1D("h_genlsp_pt",lsps.at(ilsp).pt(),weight,h_1d,nBinJetPt,0.,maxJetPt);
  }

  float genMET = genMETvec.pt();
  float genMET_eta30 = genMETvec_eta30.pt();
  float genMET_eta22 = genMETvec_eta22.pt();
  // //  float genUESumEt_ = genMETeta3.sumEt();

  plot1D("h_genMET",genMET,1,h_1d, nBinMet,0.,maxMet);  
  plot1D("h_genMET_eta30",genMET_eta30,1,h_1d, nBinMet,0.,maxMet);  
  plot1D("h_genMET_eta22",genMET_eta22,1,h_1d, nBinMet,0.,maxMet);  
  // //  plot1D("h_genSumET",genUESumEt_,1,h_1d, nBinHt, 0., maxHt);  


  // //  if(doDebug) cout << "nGenLep " << nGenLep << endl;

  // double genMETeta=genMETeta3.pt();

  // plot1D("h_genMETeta3",genMETeta,1,h_1d, nBinMet,0.,maxMet);  
  // plot2D("METcorrEta",genMETeta,genMEt_,"","",1.,h_2d,nBinMet,0.,maxMet,nBinMet,0.,maxMet);

  // require hadronic events
  if (nGenLep > 0) return;

  // --------------------------------------------------------------------------
  // retrieve the genJET objects
  //

  edm::Handle<reco::GenJetCollection> genJetsHandle;
  iEvent.getByLabel("ak5GenJets", genJetsHandle);

  int ngenjets_pt30 = 0;
  int ngenjets_pt30_eta30 = 0;
  int ngenjets_pt30_eta22 = 0;
  int ngenjets_pt40 = 0;
  int ngenjets_pt40_eta30 = 0;
  int ngenjets_pt40_eta22 = 0;
  int ngenjets_pt50 = 0;
  int ngenjets_pt50_eta30 = 0;
  int ngenjets_pt50_eta22 = 0;
  int ngenjets_pt80 = 0;
  int ngenjets_pt80_eta30 = 0;
  int ngenjets_pt80_eta22 = 0;
  int ngenjets_pt100 = 0;
  int ngenjets_pt100_eta30 = 0;
  int ngenjets_pt100_eta22 = 0;
  double genHT30=0.;
  double genHT30_eta30=0.;
  double genHT30_eta22=0.;
  double genHT40=0.;
  double genHT40_eta30=0.;
  double genHT40_eta22=0.;

  std::vector<LorentzVector> genJets;
  std::vector<LorentzVector> genJetsCent30;
  std::vector<LorentzVector> genJetsCent22;

  for (unsigned iGenJet = 0; iGenJet < genJetsHandle->size(); ++iGenJet) {
    const reco::GenJet& genJet = (*genJetsHandle) [iGenJet];
    
    genJets.push_back(genJet.p4());
    if(fabs(genJet.eta())<3.0) genJetsCent30.push_back(genJet.p4());
    if(fabs(genJet.eta())<2.172) genJetsCent22.push_back(genJet.p4());

    double genJetPt  = genJet.pt();
    double genJetEta = genJet.eta();

    // should apply some smearing to get jet pt to mock up reco..
    
    if( genJetPt > 30.0 && std::abs(genJetEta)<4.7 ) {
      ++ngenjets_pt30;
      genHT30 += genJetPt;
      if (genJetPt > 40.) {
	++ngenjets_pt40;
	genHT40 += genJetPt;
      }
      if (genJetPt > 50.) ++ngenjets_pt50;
      if (genJetPt > 80.) ++ngenjets_pt80;
      if (genJetPt > 100.) ++ngenjets_pt100;
    }
    if( genJetPt > 30.0 && std::abs(genJetEta)<3.0 ) {
      ++ngenjets_pt30_eta30;
      genHT30_eta30 += genJetPt;
      if (genJetPt > 40.) {
	++ngenjets_pt40_eta30;
	genHT40_eta30 += genJetPt;
      }
      if (genJetPt > 50.) ++ngenjets_pt50_eta30;
      if (genJetPt > 80.) ++ngenjets_pt80_eta30;
      if (genJetPt > 100.) ++ngenjets_pt100_eta30;
    }
    if( genJetPt > 30.0 && std::abs(genJetEta)<2.172 ) {
      ++ngenjets_pt30_eta22;
      genHT30_eta22 += genJetPt;
      if (genJetPt > 40.) {
	++ngenjets_pt40_eta22;
	genHT40_eta22 += genJetPt;
      }
      if (genJetPt > 50.) ++ngenjets_pt50_eta22;
      if (genJetPt > 80.) ++ngenjets_pt80_eta22;
      if (genJetPt > 100.) ++ngenjets_pt100_eta22;
    }
  }

  sort(genJets.begin(),genJets.end(),cmpPtLorentz());
  sort(genJetsCent30.begin(),genJetsCent30.end(),cmpPtLorentz());
  sort(genJetsCent22.begin(),genJetsCent22.end(),cmpPtLorentz());
  
  plot1D("h_ngenjets_pt30", ngenjets_pt30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt30_eta30", ngenjets_pt30_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt30_eta22", ngenjets_pt30_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt40", ngenjets_pt40, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt40_eta30", ngenjets_pt40_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt40_eta22", ngenjets_pt40_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt50", ngenjets_pt50, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt50_eta30", ngenjets_pt50_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt50_eta22", ngenjets_pt50_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt80", ngenjets_pt80, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt80_eta30", ngenjets_pt80_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt80_eta22", ngenjets_pt80_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt100", ngenjets_pt100, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt100_eta30", ngenjets_pt100_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_ngenjets_pt100_eta22", ngenjets_pt100_eta22, 1, h_1d, 20, 0., 20);  

  plot1D("h_genHT30", genHT30, 1, h_1d, nBinHt, 0., maxHt);  
  plot1D("h_genHT30_eta30", genHT30_eta30, 1, h_1d, nBinHt, 0., maxHt);  
  plot1D("h_genHT30_eta22", genHT30_eta22, 1, h_1d, nBinHt, 0., maxHt);  
  plot1D("h_genHT40", genHT40, 1, h_1d, nBinHt, 0., maxHt);  
  plot1D("h_genHT40_eta30", genHT40_eta30, 1, h_1d, nBinHt, 0., maxHt);  
  plot1D("h_genHT40_eta22", genHT40_eta22, 1, h_1d, nBinHt, 0., maxHt);  

  float pt_genjet1 = 0.;
  float pt_genjet2 = 0.;
  float pt_genjet3 = 0.;
  float pt_genjet4 = 0.;
  float dphi_genjet12 = -99.;

  if (genJets.size() > 0) {
    pt_genjet1 = genJets.at(0).pt();
  }
  if (genJets.size() > 1) {
    pt_genjet2 = genJets.at(1).pt();
    dphi_genjet12 = dphi_normal(genJets.at(0).phi(),genJets.at(1).phi());
  }
  if (genJets.size() > 2) {
    pt_genjet3 = genJets.at(2).pt();
  }
  if (genJets.size() > 3) {
    pt_genjet4 = genJets.at(3).pt();
  }

  plot1D("h_genjet1_pt", pt_genjet1, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet2_pt", pt_genjet2, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet3_pt", pt_genjet3, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet4_pt", pt_genjet4, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet12_dphi", dphi_genjet12, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_genjet12_dphi_vs_jet2_pt", pt_genjet2, dphi_genjet12, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  float pt_genjet1_eta30 = 0.;
  float pt_genjet2_eta30 = 0.;
  float pt_genjet3_eta30 = 0.;
  float pt_genjet4_eta30 = 0.;
  float dphi_genjet12_eta30 = -99.;

  if (genJetsCent30.size() > 0) {
    pt_genjet1_eta30 = genJetsCent30.at(0).pt();
  }
  if (genJetsCent30.size() > 1) {
    pt_genjet2_eta30 = genJetsCent30.at(1).pt();
    dphi_genjet12_eta30 = dphi_normal(genJetsCent30.at(0).phi(),genJetsCent30.at(1).phi());
  }
  if (genJetsCent30.size() > 2) {
    pt_genjet3_eta30 = genJetsCent30.at(2).pt();
  }
  if (genJetsCent30.size() > 3) {
    pt_genjet4_eta30 = genJetsCent30.at(3).pt();
  }

  plot1D("h_genjet1_pt_eta30", pt_genjet1_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet2_pt_eta30", pt_genjet2_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet3_pt_eta30", pt_genjet3_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet4_pt_eta30", pt_genjet4_eta30, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet12_dphi_eta30", dphi_genjet12_eta30, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_genjet12_dphi_vs_jet2_pt_eta30", pt_genjet2_eta30, dphi_genjet12_eta30, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  float pt_genjet1_eta22 = 0.;
  float pt_genjet2_eta22 = 0.;
  float pt_genjet3_eta22 = 0.;
  float pt_genjet4_eta22 = 0.;
  float dphi_genjet12_eta22 = -99.;

  if (genJetsCent22.size() > 0) {
    pt_genjet1_eta22 = genJetsCent22.at(0).pt();
  }
  if (genJetsCent22.size() > 1) {
    pt_genjet2_eta22 = genJetsCent22.at(1).pt();
    dphi_genjet12_eta22 = dphi_normal(genJetsCent22.at(0).phi(),genJetsCent22.at(1).phi());
  }
  if (genJetsCent22.size() > 2) {
    pt_genjet3_eta22 = genJetsCent22.at(2).pt();
  }
  if (genJetsCent22.size() > 3) {
    pt_genjet4_eta22 = genJetsCent22.at(3).pt();
  }

  plot1D("h_genjet1_pt_eta22", pt_genjet1_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet2_pt_eta22", pt_genjet2_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet3_pt_eta22", pt_genjet3_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet4_pt_eta22", pt_genjet4_eta22, weight, h_1d, nBinJetPt, 0., maxJetPt);
  plot1D("h_genjet12_dphi_eta22", dphi_genjet12_eta22, weight, h_1d, 100, -3.5, 3.5);
  plot2D("h_genjet12_dphi_vs_jet2_pt_eta22", pt_genjet2_eta22, dphi_genjet12_eta22, "2nd jet PT [GeV]", "#Delta#phi(j1,j2)", weight, h_2d, nBinJetPt, 0., maxJetPt, 100, -3.5, 3.5);

  // calculate genMT2
  float genmt2_jeteta30 = CalcHemisphereAndMT2(0.,false,genJetsCent30,genMETvec);
  //  float genmt2_jeteta30_meteta30 = CalcHemisphereAndMT2(0.,false,genJetsCent30,genMETvec_eta30);
  float genmt2_jeteta22 = CalcHemisphereAndMT2(0.,false,genJetsCent22,genMETvec);
  //  float genmt2_jeteta22_meteta30 = CalcHemisphereAndMT2(0.,false,genJetsCent22,genMETvec_eta30);

  plot1D("h_genmt2_jeteta30", genmt2_jeteta30, weight, h_1d, nBinMT2, 0., maxMT2);
  //  plot1D("h_genmt2_jeteta30_meteta30", genmt2_jeteta30_meteta30, weight, h_1d, nBinMT2, 0., maxMT2);
  plot2D("h_genjet12_dphi_vs_genmt2_jeteta30", genmt2_jeteta30, dphi_genjet12_eta30, "Gen M_{T2} [GeV]", "Gen #Delta#phi(j1,j2)", weight, h_2d, nBinMT2, 0., maxMT2, 100, -3.5, 3.5);

  plot1D("h_genmt2_jeteta22", genmt2_jeteta22, weight, h_1d, nBinMT2, 0., maxMT2);
  //  plot1D("h_genmt2_jeteta22_meteta30", genmt2_jeteta22_meteta30, weight, h_1d, nBinMT2, 0., maxMT2);
  plot2D("h_genjet12_dphi_vs_genmt2_jeteta22", genmt2_jeteta22, dphi_genjet12_eta22, "Gen M_{T2} [GeV]", "Gen #Delta#phi(j1,j2)", weight, h_2d, nBinMT2, 0., maxMT2, 100, -3.5, 3.5);



  // //
  // // ----------------------------------------------------------------------
  // // loop back over gen particles to compare to gen jets..
  // //

  // // this is generator-specific and current works for Pythia8! 
  // // See pythia8xxx/htmldoc/ParticleProperties.html for status codes
  // // 21: incoming parton
  // // 22: intermediate particle in hard process
  // // 23: outgoing parton

  // std::vector<double> genjet_sumgenpartpt(genJetsCent30.size(),0.);
  // std::vector<double> genjet_sumgenpartpt_wmu(genJetsCent30.size(),0.);
  // std::vector<double> genjet_sumgenpartpt_wlsp(genJetsCent30.size(),0.);
  // std::vector<double> genjet_sumgenpartpt_wmulsp(genJetsCent30.size(),0.);
                                                                                        
  // for (reco::GenParticleCollection::const_iterator genp = genparts->begin(); genp != genparts->end(); ++ genp ) {

  //   if(genp->status()!=1) continue;
  //   if(abs(genp->pdgId()) == 12 || abs(genp->pdgId()) == 14 || abs(genp->pdgId()) == 16) continue; // remove neutrinos
  //   if(abs(genp->pdgId()) == 1000012 || abs(genp->pdgId())==1000014 || abs(genp->pdgId())==1000016 ) continue;
  //   if(abs(genp->pdgId()) == 2000012 || abs(genp->pdgId())==2000014 || abs(genp->pdgId())==2000016 ) continue;
  //   if(abs(genp->pdgId()) == 1000039 || abs(genp->pdgId())==5100039 ) continue;
  //   if(abs(genp->pdgId()) == 4000012 || abs(genp->pdgId())==4000014 || abs(genp->pdgId())==4000016 ) continue;
  //   if(abs(genp->pdgId()) == 9900012 || abs(genp->pdgId())==9900014 || abs(genp->pdgId())==9900016 ) continue;
  //   if(abs(genp->pdgId()) == 39) continue;

  //   // loop over genjets and check for dR match
  //   for (unsigned int igenjet = 0; igenjet < genJetsCent30.size(); ++igenjet) {
  //     float dR = deltaR_normal(genp->eta(),genp->phi(),genJetsCent30.at(igenjet).eta(),genJetsCent30.at(igenjet).phi());
  //     if (dR > 0.5) continue;
  //     genjet_sumgenpartpt_wmulsp.at(igenjet) += genp->pt();
  //     // muons
  //     if(abs(genp->pdgId()) == 13) {
  //       genjet_sumgenpartpt_wmu.at(igenjet) += genp->pt();
  //     }
  //     // LSPs
  //     else if(abs(genp->pdgId()) == 1000022) {
  //       genjet_sumgenpartpt_wlsp.at(igenjet) += genp->pt();
  //     } else {
  //       genjet_sumgenpartpt_wlsp.at(igenjet) += genp->pt();
  //       genjet_sumgenpartpt_wmu.at(igenjet) += genp->pt();
  // 	genjet_sumgenpartpt.at(igenjet) += genp->pt();
  //     }
  //   } // loop over genjets

  // }

  // std::cout << "-- comparing genjets pt to gen particles" << std::endl;
  // // loop back over genjets and compare pt to sumpt of dR 0.5 cone
  // for (unsigned int igenjet = 0; igenjet < genJetsCent30.size(); ++igenjet) {
  //   std::cout << " " << igenjet << ", genjet pt: " << genJetsCent30.at(igenjet).pt()
  // 	      << ", sum genps: " << genjet_sumgenpartpt.at(igenjet)
  // 	      << ", sum genps w mus: " << genjet_sumgenpartpt_wmu.at(igenjet)
  // 	      << ", sum genps w lsps: " << genjet_sumgenpartpt_wlsp.at(igenjet)
  // 	      << ", sum genps w both: " << genjet_sumgenpartpt_wmulsp.at(igenjet)
  // 	      << std::endl;
  // }

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

  int njets_pt30 = 0;
  int njets_pt30_eta30 = 0;
  int njets_pt30_eta22 = 0;
  int njets_pt40 = 0;
  int njets_pt40_eta30 = 0;
  int njets_pt40_eta22 = 0;
  int njets_pt50 = 0;
  int njets_pt50_eta30 = 0;
  int njets_pt50_eta22 = 0;
  int njets_pt80 = 0;
  int njets_pt80_eta30 = 0;
  int njets_pt80_eta22 = 0;
  int njets_pt100 = 0;
  int njets_pt100_eta30 = 0;
  int njets_pt100_eta22 = 0;

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

      newJets2015.push_back(jetIter->p4());
      if (et > 30.) ++njets_pt30;
      if (et > 40.) ++njets_pt40;
      if (et > 50.) ++njets_pt50;
      if (et > 80.) ++njets_pt80;
      if (et > 100.) ++njets_pt100;
      if (fabs(eta) < 3.0) {
	newJets2015Cent30.push_back(jetIter->p4());
	if (et > 30.) ++njets_pt30_eta30;
	if (et > 40.) ++njets_pt40_eta30;
	if (et > 50.) ++njets_pt50_eta30;
	if (et > 80.) ++njets_pt80_eta30;
	if (et > 100.) ++njets_pt100_eta30;
      }
      if (fabs(eta) < 2.172) {
	newJets2015Cent22.push_back(jetIter->p4());
 	if (et > 30.) ++njets_pt30_eta22;
	if (et > 40.) ++njets_pt40_eta22;
	if (et > 50.) ++njets_pt50_eta22;
	if (et > 80.) ++njets_pt80_eta22;
	if (et > 100.) ++njets_pt100_eta22;
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
      if (et > 30.) ++njets_pt30;
      if (et > 40.) ++njets_pt40;
      if (et > 50.) ++njets_pt50;
    }
  }

  sort(newJets2015.begin(),newJets2015.end(),cmpPtLorentz());
  sort(newJets2015Cent30.begin(),newJets2015Cent30.end(),cmpPtLorentz());
  sort(newJets2015Cent22.begin(),newJets2015Cent22.end(),cmpPtLorentz());

  plot1D("h_njets_pt30", njets_pt30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt30_eta30", njets_pt30_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt30_eta22", njets_pt30_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt40", njets_pt40, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt40_eta30", njets_pt40_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt40_eta22", njets_pt40_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt50", njets_pt50, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt50_eta30", njets_pt50_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt50_eta22", njets_pt50_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt80", njets_pt80, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt80_eta30", njets_pt80_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt80_eta22", njets_pt80_eta22, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt100", njets_pt100, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt100_eta30", njets_pt100_eta30, 1, h_1d, 20, 0., 20);  
  plot1D("h_njets_pt100_eta22", njets_pt100_eta22, 1, h_1d, 20, 0., 20);  

  plot2D("h_njets_pt30_vs_ngenjets_pt30", ngenjets_pt30, njets_pt30, "N(gen jets, pt30)", "N(L1 jets, pt30)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt30_eta30_vs_ngenjets_pt30_eta30", ngenjets_pt30_eta30, njets_pt30_eta30, "N(gen jets, pt30)", "N(L1 jets, pt30)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt30_eta22_vs_ngenjets_pt30_eta22", ngenjets_pt30_eta22, njets_pt30_eta22, "N(gen jets, pt30)", "N(L1 jets, pt30)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt40_vs_ngenjets_pt40", ngenjets_pt40, njets_pt40, "N(gen jets, pt40)", "N(L1 jets, pt40)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt40_eta30_vs_ngenjets_pt40_eta30", ngenjets_pt40_eta30, njets_pt40_eta30, "N(gen jets, pt40)", "N(L1 jets, pt40)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt40_eta22_vs_ngenjets_pt40_eta22", ngenjets_pt40_eta22, njets_pt40_eta22, "N(gen jets, pt40)", "N(L1 jets, pt40)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt50_vs_ngenjets_pt50", ngenjets_pt50, njets_pt50, "N(gen jets, pt50)", "N(L1 jets, pt50)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt50_eta30_vs_ngenjets_pt50_eta30", ngenjets_pt50_eta30, njets_pt50_eta30, "N(gen jets, pt50)", "N(L1 jets, pt50)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt50_eta22_vs_ngenjets_pt50_eta22", ngenjets_pt50_eta22, njets_pt50_eta22, "N(gen jets, pt50)", "N(L1 jets, pt50)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt80_vs_ngenjets_pt80", ngenjets_pt80, njets_pt80, "N(gen jets, pt80)", "N(L1 jets, pt80)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt80_eta30_vs_ngenjets_pt80_eta30", ngenjets_pt80_eta30, njets_pt80_eta30, "N(gen jets, pt80)", "N(L1 jets, pt80)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt80_eta22_vs_ngenjets_pt80_eta22", ngenjets_pt80_eta22, njets_pt80_eta22, "N(gen jets, pt80)", "N(L1 jets, pt80)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt100_vs_ngenjets_pt100", ngenjets_pt100, njets_pt100, "N(gen jets, pt100)", "N(L1 jets, pt100)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt100_eta30_vs_ngenjets_pt100_eta30", ngenjets_pt100_eta30, njets_pt100_eta30, "N(gen jets, pt100)", "N(L1 jets, pt100)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);
  plot2D("h_njets_pt100_eta22_vs_ngenjets_pt100_eta22", ngenjets_pt100_eta22, njets_pt100_eta22, "N(gen jets, pt100)", "N(L1 jets, pt100)", weight, h_2d, nBinNJets,0.,maxNJets,nBinNJets,0.,maxNJets);


  std::vector<int> jetThresholds;
  jetThresholds.push_back(40);
  jetThresholds.push_back(50);
  jetThresholds.push_back(60);
  jetThresholds.push_back(70);
  jetThresholds.push_back(80);
  jetThresholds.push_back(90);
  jetThresholds.push_back(100);

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

  // turn on curves, all eta
  makeTurnOn(pt_genjet1, "genjet1_pt", pt_jet1, "jet1_pt", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet2, "genjet2_pt", pt_jet2, "jet2_pt", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12) < 3.0) makeTurnOn(pt_genjet2, "genjet2_pt", pt_jet2, "jet2_pt_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12) < 2.6) makeTurnOn(pt_genjet2, "genjet2_pt", pt_jet2, "jet2_pt_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet3, "genjet3_pt", pt_jet3, "jet3_pt", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet4, "genjet4_pt", pt_jet4, "jet4_pt", jetThresholds, weight, nBinJetPt, maxJetPt);

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

  // turn on curves, eta 30
  makeTurnOn(pt_genjet1_eta30, "genjet1_pt_eta30", pt_jet1_eta30, "jet1_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30", pt_jet2_eta30, "jet2_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12_eta30) < 3.0) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30", pt_jet2_eta30, "jet2_pt_eta30_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12_eta30) < 2.6) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30", pt_jet2_eta30, "jet2_pt_eta30_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet3_eta30, "genjet3_pt_eta30", pt_jet3_eta30, "jet3_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet4_eta30, "genjet4_pt_eta30", pt_jet4_eta30, "jet4_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);

  // turn-ons with gen mt2 cut
  if (genmt2_jeteta30 > 100.) {
    makeTurnOn(pt_genjet1_eta30, "genjet1_pt_eta30_genmt2_jeteta30_cut100", pt_jet1_eta30, "jet1_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut100", pt_jet2_eta30, "jet2_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta30) < 3.0) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut100", pt_jet2_eta30, "jet2_pt_eta30_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta30) < 2.6) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut100", pt_jet2_eta30, "jet2_pt_eta30_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet3_eta30, "genjet3_pt_eta30_genmt2_jeteta30_cut100", pt_jet3_eta30, "jet3_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet4_eta30, "genjet4_pt_eta30_genmt2_jeteta30_cut100", pt_jet4_eta30, "jet4_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
  }
  if (genmt2_jeteta30 > 200.) {
    makeTurnOn(pt_genjet1_eta30, "genjet1_pt_eta30_genmt2_jeteta30_cut200", pt_jet1_eta30, "jet1_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut200", pt_jet2_eta30, "jet2_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta30) < 3.0) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut200", pt_jet2_eta30, "jet2_pt_eta30_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta30) < 2.6) makeTurnOn(pt_genjet2_eta30, "genjet2_pt_eta30_genmt2_jeteta30_cut200", pt_jet2_eta30, "jet2_pt_eta30_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet3_eta30, "genjet3_pt_eta30_genmt2_jeteta30_cut200", pt_jet3_eta30, "jet3_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet4_eta30, "genjet4_pt_eta30_genmt2_jeteta30_cut200", pt_jet4_eta30, "jet4_pt_eta30", jetThresholds, weight, nBinJetPt, maxJetPt);
  }



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

  // turn on curves, eta 22
  makeTurnOn(pt_genjet1_eta22, "genjet1_pt_eta22", pt_jet1_eta22, "jet1_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22", pt_jet2_eta22, "jet2_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12_eta22) < 3.0) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22", pt_jet2_eta22, "jet2_pt_eta22_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
  if (fabs(dphi_jet12_eta22) < 2.6) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22", pt_jet2_eta22, "jet2_pt_eta22_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet3_eta22, "genjet3_pt_eta22", pt_jet3_eta22, "jet3_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
  makeTurnOn(pt_genjet4_eta22, "genjet4_pt_eta22", pt_jet4_eta22, "jet4_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);

  // turn-ons with gen mt2 cut
  if (genmt2_jeteta22 > 100.) {
    makeTurnOn(pt_genjet1_eta22, "genjet1_pt_eta22_genmt2_jeteta22_cut100", pt_jet1_eta22, "jet1_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut100", pt_jet2_eta22, "jet2_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta22) < 3.0) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut100", pt_jet2_eta22, "jet2_pt_eta22_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta22) < 2.6) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut100", pt_jet2_eta22, "jet2_pt_eta22_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet3_eta22, "genjet3_pt_eta22_genmt2_jeteta22_cut100", pt_jet3_eta22, "jet3_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet4_eta22, "genjet4_pt_eta22_genmt2_jeteta22_cut100", pt_jet4_eta22, "jet4_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
  }
  if (genmt2_jeteta22 > 200.) {
    makeTurnOn(pt_genjet1_eta22, "genjet1_pt_eta22_genmt2_jeteta22_cut200", pt_jet1_eta22, "jet1_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut200", pt_jet2_eta22, "jet2_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta22) < 3.0) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut200", pt_jet2_eta22, "jet2_pt_eta22_dphiL", jetThresholds, weight, nBinJetPt, maxJetPt);
    if (fabs(dphi_jet12_eta22) < 2.6) makeTurnOn(pt_genjet2_eta22, "genjet2_pt_eta22_genmt2_jeteta22_cut200", pt_jet2_eta22, "jet2_pt_eta22_dphiT", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet3_eta22, "genjet3_pt_eta22_genmt2_jeteta22_cut200", pt_jet3_eta22, "jet3_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
    makeTurnOn(pt_genjet4_eta22, "genjet4_pt_eta22_genmt2_jeteta22_cut200", pt_jet4_eta22, "jet4_pt_eta22", jetThresholds, weight, nBinJetPt, maxJetPt);
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
  
  plot1D("h_caloMHT",caloMHt_,1,h_1d, nBinMet,0.,maxMet);
  plot1D("h_caloHT",caloHt_,1,h_1d, nBinHt,0.,maxHt);

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
  
  plot1D("h_caloMHT2015",caloMHt2015_,1,h_1d, nBinMet,0.,maxMet);
  plot1D("h_caloHT2015",caloHt2015_,1,h_1d, nBinHt,0.,maxHt);

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

  // HT thresholds for turn-ons
  std::vector<int> htThresholds;
  htThresholds.push_back(110);
  htThresholds.push_back(125);
  htThresholds.push_back(150);
  htThresholds.push_back(175);
  htThresholds.push_back(200);
  htThresholds.push_back(225);

  // HT plots and turn on curves
  makeTurnOn(genHT40_eta30, "genHT40_eta30", HT_reget7, "HT_reget7", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta30, "genHT40_eta30", HT_reget10, "HT_reget10", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta30, "genHT40_eta30", HT_reget15, "HT_reget15", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta30, "genHT40_eta30", HT_reget20, "HT_reget20", htThresholds, weight, nBinHt, maxHt);

  if (genmt2_jeteta30 > 100.) {
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut100", HT_reget7, "HT_reget7", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut100", HT_reget10, "HT_reget10", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut100", HT_reget15, "HT_reget15", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut100", HT_reget20, "HT_reget20", htThresholds, weight, nBinHt, maxHt);
  }
  if (genmt2_jeteta30 > 200.) {
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut200", HT_reget7, "HT_reget7", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut200", HT_reget10, "HT_reget10", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut200", HT_reget15, "HT_reget15", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta30, "genHT40_eta30_genmt2_jeteta30_cut200", HT_reget20, "HT_reget20", htThresholds, weight, nBinHt, maxHt);
  }


  makeTurnOn(genHT40_eta22, "genHT40_eta22", HT_eta22_reget7, "HT_eta22_reget7", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta22, "genHT40_eta22", HT_eta22_reget10, "HT_eta22_reget10", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta22, "genHT40_eta22", HT_eta22_reget15, "HT_eta22_reget15", htThresholds, weight, nBinHt, maxHt);
  makeTurnOn(genHT40_eta22, "genHT40_eta22", HT_eta22_reget20, "HT_eta22_reget20", htThresholds, weight, nBinHt, maxHt);

  if (genmt2_jeteta22 > 100.) {
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut100", HT_eta22_reget7, "HT_eta22_reget7", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut100", HT_eta22_reget10, "HT_eta22_reget10", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut100", HT_eta22_reget15, "HT_eta22_reget15", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut100", HT_eta22_reget20, "HT_eta22_reget20", htThresholds, weight, nBinHt, maxHt);
  }
  if (genmt2_jeteta22 > 200.) {
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut200", HT_eta22_reget7, "HT_eta22_reget7", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut200", HT_eta22_reget10, "HT_eta22_reget10", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut200", HT_eta22_reget15, "HT_eta22_reget15", htThresholds, weight, nBinHt, maxHt);
    makeTurnOn(genHT40_eta22, "genHT40_eta22_genmt2_jeteta22_cut200", HT_eta22_reget20, "HT_eta22_reget20", htThresholds, weight, nBinHt, maxHt);
  }


}


void L1Analyzer::makeTurnOn(double pt, string var, double l1pt, string l1var, std::vector<int> thresh, double weight, int nBin, int max) {

  plot1D("h_"+var+"_"+l1var+"_denom", pt, weight, h_1d,  nBin, 0.,max); 

  for (unsigned int ithresh = 0; ithresh < thresh.size(); ++ithresh) {
    if (l1pt >= thresh.at(ithresh))   plot1D(Form("h_%s_%s_%d",var.c_str(),l1var.c_str(),thresh.at(ithresh)), pt, weight, h_1d, nBin, 0.,max); 
  }

  return;
}
  
Int_t L1Analyzer::phiINjetCoord(Double_t phi) {

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

Int_t L1Analyzer::etaINjetCoord(Double_t eta) {

  size_t etaIdx = 0.;
  for (size_t idx=0; idx<ETABINS; idx++) {
    if (eta>=ETABIN[idx] and eta<ETABIN[idx+1])
      etaIdx = idx;
  }

  return int(etaIdx);
}

inline Double_t L1Analyzer::degree(Double_t radian) {
  if (radian<0)
    return 360.+(radian/TMath::Pi()*180.);
  else
    return radian/TMath::Pi()*180.;
}

float L1Analyzer::dphi_normal(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  if (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  else if (dphi > TMath::Pi()) dphi -= 2*TMath::Pi();

  return dphi;
}

float L1Analyzer::dphi(float phi1, float phi2) {
  float dphi = phi1 - phi2;
  //  if (dphi <= -TMath::Pi()) dphi += 2*TMath::Pi();
  if (dphi < -3.0) dphi += 2*TMath::Pi();
  else if (dphi > 3.2) dphi -= 2*TMath::Pi();
  //  std::cout << "dphi calc: phi1: " << phi1 << ", phi2: " << phi2 << ", diff: " << phi1-phi2 << ", dphi: " << dphi << std::endl;

  return dphi;
}

float L1Analyzer::deltaR_normal(float eta1, float phi1, float eta2, float phi2) {
  float diffeta = eta1 - eta2;
  float diffphi = dphi_normal(phi1,phi2);
  return sqrt(diffeta*diffeta + diffphi*diffphi);
}

//____________________________________________________________
// calculate hemispheres and input them to MT2 calc
Float_t L1Analyzer::CalcHemisphereAndMT2(float testmass, bool massive, std::vector<LorentzVector> jets, LorentzVector MET ) {

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
Float_t L1Analyzer::CalcMT2(float testmass, bool massive, LorentzVector visible1, LorentzVector visible2, LorentzVector MET ){

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

  // Davismt2 *mt2 = new Davismt2();
  // mt2->set_momenta(pa, pb, pmiss);
  // mt2->set_mn(testmass);
  // Float_t MT2=mt2->get_mt2();
  // delete mt2;

  Davismt2 mt2;
  mt2.set_momenta(pa, pb, pmiss);
  mt2.set_mn(testmass);
  Float_t MT2=mt2.get_mt2();

  if (doDebug) std::cout << "  - MT2 value: " << MT2 << std::endl;
  return MT2;
}

// ------------ method called once each job just before starting event loop  ------------
void 
L1Analyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
L1Analyzer::endJob() 
{

  //  TFile*fout = new TFile("ttbar_controlplots.root","RECREATE");
  //  TFile*fout = new TFile("T2tt_controlplots.root","RECREATE");
  TFile*fout = new TFile(histoname_.data(),"RECREATE");
  
  std::map<std::string, TH1F*>::iterator it1d;
  for(it1d=h_1d.begin(); it1d!=h_1d.end(); it1d++) {
    it1d->second->Write();
    delete it1d->second;
  }

  std::map<std::string, TH2F*>::iterator it2d;
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
L1Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(L1Analyzer);
