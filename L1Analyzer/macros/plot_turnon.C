#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>

TH1F* turnon_plot(TH1* h_num, TH1* h_denom, int rebin = 1);
TCanvas* compare_turnons(std::vector<TH1D*> hists_num, TH1* h_denom, std::vector<std::string> labels, std::string name, int rebin = 1);

bool savePlots_ = false;
TString outdir_ = "";
TFile* fin = 0;
TFile* outfile = 0;

//______________________________________________________
void plot_turnon(TString input_file = "histos_T2qq_400_150_turnons_all.root", bool savePlots = false, TString outdir = "./" ) {

  savePlots_ = savePlots;
  outdir_ = outdir;

  outfile = new TFile("turnons.root","RECREATE");

  fin = new TFile(input_file);

  // ---------------------------------------
  // ----------- Dijets --------------------
  // ---------------------------------------

  int rebin = 4;

  std::vector<std::string> jet2_labels;
  jet2_labels.push_back("jet2 pt > 60");
  jet2_labels.push_back("jet2 pt > 70");
  jet2_labels.push_back("jet2 pt > 80");
  jet2_labels.push_back("jet2 pt > 90");
  jet2_labels.push_back("jet2 pt > 100");

  // --- 2 jets turn on
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_60 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_60");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_70 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_70");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_80 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_80");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_90 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_90");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_100 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_100");
  TH1F* h_denom_genjet2_pt_eta30_jet2_pt_eta30 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_denom");
  if (rebin > 1) h_denom_genjet2_pt_eta30_jet2_pt_eta30->Rebin(rebin);

  std::vector<TH1F*> hists_num_jet2;
  hists_num_jet2.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_60);
  hists_num_jet2.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_70);
  hists_num_jet2.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_80);
  hists_num_jet2.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_90);
  hists_num_jet2.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_100);

  compare_turnons(hists_num_jet2,h_denom_genjet2_pt_eta30_jet2_pt_eta30,jet2_labels,"jet2_eta30",rebin);

  // --- 2 jets turn on, dphiL, no mt2 cut
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_60 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiL_60");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_70 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiL_70");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_80 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiL_80");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_90 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiL_90");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_100 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiL_100");

  std::vector<TH1F*> hists_num_jet2_dphiL;
  hists_num_jet2_dphiL.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_60);
  hists_num_jet2_dphiL.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_70);
  hists_num_jet2_dphiL.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_80);
  hists_num_jet2_dphiL.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_90);
  hists_num_jet2_dphiL.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiL_100);

  // denom without dphiL cut
  compare_turnons(hists_num_jet2_dphiL,h_denom_genjet2_pt_eta30_jet2_pt_eta30,jet2_labels,"jet2_eta30_dphiL",rebin);

  // --- 2 jets turn on, dphiL, mt2 > 100
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_60 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_60");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_70 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_70");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_80 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_80");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_90 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_90");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_100 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_100");
  TH1F* h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_denom");
  if (rebin > 1) h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30->Rebin(rebin);

  std::vector<TH1F*> hists_num_jet2_mt2_cut100_dphiL;
  hists_num_jet2_mt2_cut100_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_60);
  hists_num_jet2_mt2_cut100_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_70);
  hists_num_jet2_mt2_cut100_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_80);
  hists_num_jet2_mt2_cut100_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_90);
  hists_num_jet2_mt2_cut100_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiL_100);

  // denom without dphiL cut
  compare_turnons(hists_num_jet2_mt2_cut100_dphiL,h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30,jet2_labels,"jet2_eta30_mt2_cut100_dphiL",rebin);

  // --- 2 jets turn on, dphiL, mt2 > 200
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_60 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_60");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_70 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_70");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_80 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_80");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_90 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_90");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_100 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_100");
  TH1F* h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_denom");
  if (rebin > 1) h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30->Rebin(rebin);

  std::vector<TH1F*> hists_num_jet2_mt2_cut200_dphiL;
  hists_num_jet2_mt2_cut200_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_60);
  hists_num_jet2_mt2_cut200_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_70);
  hists_num_jet2_mt2_cut200_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_80);
  hists_num_jet2_mt2_cut200_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_90);
  hists_num_jet2_mt2_cut200_dphiL.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiL_100);

  // denom without dphiL cut
  compare_turnons(hists_num_jet2_mt2_cut200_dphiL,h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30,jet2_labels,"jet2_eta30_mt2_cut200_dphiL",rebin);

  // --- 2 jets turn on, dphiT, no mt2 cut
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_60 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiT_60");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_70 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiT_70");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_80 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiT_80");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_90 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiT_90");
  TH1F* h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_100 = fin->Get("h_genjet2_pt_eta30_jet2_pt_eta30_dphiT_100");

  std::vector<TH1F*> hists_num_jet2_dphiT;
  hists_num_jet2_dphiT.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_60);
  hists_num_jet2_dphiT.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_70);
  hists_num_jet2_dphiT.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_80);
  hists_num_jet2_dphiT.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_90);
  hists_num_jet2_dphiT.push_back(h_num_genjet2_pt_eta30_jet2_pt_eta30_dphiT_100);

  // denom without dphiT cut
  compare_turnons(hists_num_jet2_dphiT,h_denom_genjet2_pt_eta30_jet2_pt_eta30,jet2_labels,"jet2_eta30_dphiT",rebin);

  // --- 2 jets turn on, dphiT, mt2 > 100
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_60 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_60");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_70 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_70");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_80 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_80");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_90 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_90");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_100 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_100");

  std::vector<TH1F*> hists_num_jet2_mt2_cut100_dphiT;
  hists_num_jet2_mt2_cut100_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_60);
  hists_num_jet2_mt2_cut100_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_70);
  hists_num_jet2_mt2_cut100_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_80);
  hists_num_jet2_mt2_cut100_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_90);
  hists_num_jet2_mt2_cut100_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30_dphiT_100);

  // denom without dphiT cut
  compare_turnons(hists_num_jet2_mt2_cut100_dphiT,h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut100_jet2_pt_eta30,jet2_labels,"jet2_eta30_mt2_cut100_dphiT",rebin);

  // --- 2 jets turn on, dphiT, mt2 > 200
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_60 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_60");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_70 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_70");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_80 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_80");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_90 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_90");
  TH1F* h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_100 = fin->Get("h_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_100");

  std::vector<TH1F*> hists_num_jet2_mt2_cut200_dphiT;
  hists_num_jet2_mt2_cut200_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_60);
  hists_num_jet2_mt2_cut200_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_70);
  hists_num_jet2_mt2_cut200_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_80);
  hists_num_jet2_mt2_cut200_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_90);
  hists_num_jet2_mt2_cut200_dphiT.push_back(h_num_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30_dphiT_100);

  // denom without dphiT cut
  compare_turnons(hists_num_jet2_mt2_cut200_dphiT,h_denom_genjet2_pt_eta30_genmt2_jeteta30_cut200_jet2_pt_eta30,jet2_labels,"jet2_eta30_mt2_cut200_dphiT",rebin);

  // ---------------------------------------
  // ----------- HT ------------------------
  // ---------------------------------------

  int rebin_ht = 10;

  std::vector<std::string> ht_labels;
  ht_labels.push_back("HT > 110");
  ht_labels.push_back("HT > 125");
  ht_labels.push_back("HT > 150");
  ht_labels.push_back("HT > 175");
  ht_labels.push_back("HT > 200");
  ht_labels.push_back("HT > 225");

  // --- HT turn on, reget7
  TH1F* h_num_genHT40_eta30_HT_reget7_110 = fin->Get("h_genHT40_eta30_HT_reget7_110");
  TH1F* h_num_genHT40_eta30_HT_reget7_125 = fin->Get("h_genHT40_eta30_HT_reget7_125");
  TH1F* h_num_genHT40_eta30_HT_reget7_150 = fin->Get("h_genHT40_eta30_HT_reget7_150");
  TH1F* h_num_genHT40_eta30_HT_reget7_175 = fin->Get("h_genHT40_eta30_HT_reget7_175");
  TH1F* h_num_genHT40_eta30_HT_reget7_200 = fin->Get("h_genHT40_eta30_HT_reget7_200");
  TH1F* h_num_genHT40_eta30_HT_reget7_225 = fin->Get("h_genHT40_eta30_HT_reget7_225");
  TH1F* h_denom_genHT40_eta30_HT_reget7 = fin->Get("h_genHT40_eta30_HT_reget7_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta30_HT_reget7->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_reget7;
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_110);
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_125);
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_150);
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_175);
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_200);
  hists_num_HT_reget7.push_back(h_num_genHT40_eta30_HT_reget7_225);

  compare_turnons(hists_num_HT_reget7,h_denom_genHT40_eta30_HT_reget7,ht_labels,"HT_reget7",rebin_ht);

  // --- HT turn on, reget10
  TH1F* h_num_genHT40_eta30_HT_reget10_110 = fin->Get("h_genHT40_eta30_HT_reget10_110");
  TH1F* h_num_genHT40_eta30_HT_reget10_125 = fin->Get("h_genHT40_eta30_HT_reget10_125");
  TH1F* h_num_genHT40_eta30_HT_reget10_150 = fin->Get("h_genHT40_eta30_HT_reget10_150");
  TH1F* h_num_genHT40_eta30_HT_reget10_175 = fin->Get("h_genHT40_eta30_HT_reget10_175");
  TH1F* h_num_genHT40_eta30_HT_reget10_200 = fin->Get("h_genHT40_eta30_HT_reget10_200");
  TH1F* h_num_genHT40_eta30_HT_reget10_225 = fin->Get("h_genHT40_eta30_HT_reget10_225");
  TH1F* h_denom_genHT40_eta30_HT_reget10 = fin->Get("h_genHT40_eta30_HT_reget10_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta30_HT_reget10->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_reget10;
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_110);
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_125);
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_150);
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_175);
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_200);
  hists_num_HT_reget10.push_back(h_num_genHT40_eta30_HT_reget10_225);

  compare_turnons(hists_num_HT_reget10,h_denom_genHT40_eta30_HT_reget10,ht_labels,"HT_reget10",rebin_ht);

  // --- HT turn on, reget20
  TH1F* h_num_genHT40_eta30_HT_reget20_110 = fin->Get("h_genHT40_eta30_HT_reget20_110");
  TH1F* h_num_genHT40_eta30_HT_reget20_125 = fin->Get("h_genHT40_eta30_HT_reget20_125");
  TH1F* h_num_genHT40_eta30_HT_reget20_150 = fin->Get("h_genHT40_eta30_HT_reget20_150");
  TH1F* h_num_genHT40_eta30_HT_reget20_175 = fin->Get("h_genHT40_eta30_HT_reget20_175");
  TH1F* h_num_genHT40_eta30_HT_reget20_200 = fin->Get("h_genHT40_eta30_HT_reget20_200");
  TH1F* h_num_genHT40_eta30_HT_reget20_225 = fin->Get("h_genHT40_eta30_HT_reget20_225");
  TH1F* h_denom_genHT40_eta30_HT_reget20 = fin->Get("h_genHT40_eta30_HT_reget20_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta30_HT_reget20->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_reget20;
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_110);
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_125);
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_150);
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_175);
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_200);
  hists_num_HT_reget20.push_back(h_num_genHT40_eta30_HT_reget20_225);

  compare_turnons(hists_num_HT_reget20,h_denom_genHT40_eta30_HT_reget20,ht_labels,"HT_reget20",rebin_ht);

  // --- HT turn on, reget7, eta22
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_110 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_110");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_125 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_125");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_150 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_150");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_175 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_175");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_200 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_200");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget7_225 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_225");
  TH1F* h_denom_genHT40_eta22_HT_eta22_reget7 = fin->Get("h_genHT40_eta22_HT_eta22_reget7_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta22_HT_eta22_reget7->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_eta22_reget7;
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_110);
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_125);
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_150);
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_175);
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_200);
  hists_num_HT_eta22_reget7.push_back(h_num_genHT40_eta22_HT_eta22_reget7_225);

  compare_turnons(hists_num_HT_eta22_reget7,h_denom_genHT40_eta22_HT_eta22_reget7,ht_labels,"HT_eta22_reget7",rebin_ht);

  // --- HT turn on, reget10, eta22
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_110 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_110");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_125 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_125");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_150 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_150");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_175 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_175");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_200 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_200");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget10_225 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_225");
  TH1F* h_denom_genHT40_eta22_HT_eta22_reget10 = fin->Get("h_genHT40_eta22_HT_eta22_reget10_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta22_HT_eta22_reget10->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_eta22_reget10;
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_110);
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_125);
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_150);
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_175);
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_200);
  hists_num_HT_eta22_reget10.push_back(h_num_genHT40_eta22_HT_eta22_reget10_225);

  compare_turnons(hists_num_HT_eta22_reget10,h_denom_genHT40_eta22_HT_eta22_reget10,ht_labels,"HT_eta22_reget10",rebin_ht);

  // --- HT turn on, reget20, eta22
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_110 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_110");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_125 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_125");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_150 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_150");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_175 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_175");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_200 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_200");
  TH1F* h_num_genHT40_eta22_HT_eta22_reget20_225 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_225");
  TH1F* h_denom_genHT40_eta22_HT_eta22_reget20 = fin->Get("h_genHT40_eta22_HT_eta22_reget20_denom");
  if (rebin_ht > 1) h_denom_genHT40_eta22_HT_eta22_reget20->Rebin(rebin_ht);

  std::vector<TH1F*> hists_num_HT_eta22_reget20;
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_110);
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_125);
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_150);
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_175);
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_200);
  hists_num_HT_eta22_reget20.push_back(h_num_genHT40_eta22_HT_eta22_reget20_225);

  compare_turnons(hists_num_HT_eta22_reget20,h_denom_genHT40_eta22_HT_eta22_reget20,ht_labels,"HT_eta22_reget20",rebin_ht);

  if (outfile) outfile->Close();

}

//______________________________________________________
TH1F* turnon_plot(TH1* h_num, TH1* h_denom, int rebin) {

  if (rebin > 1) {
    h_num->Rebin(rebin);
  }

  TH1F* h_turnon = h_num->Clone(Form("%s_turnon",h_num->GetName()));
  h_turnon->Divide(h_turnon,h_denom,1.0,1.0,"B");
  if (TString(h_turnon->GetName()).Contains("genjet2")) {
    h_turnon->GetXaxis()->SetRangeUser(0,299);
    h_turnon->GetXaxis()->SetTitle("Subleading GenJet p_{T} [GeV]");
  }
  else if (TString(h_turnon->GetName()).Contains("genHT")) {
    h_turnon->GetXaxis()->SetRangeUser(0,599);
    h_turnon->GetXaxis()->SetTitle("Gen H_{T} [GeV]");
  }
  h_turnon->GetYaxis()->SetTitle("L1 Efficiency wrt Gen");
  h_turnon->GetYaxis()->SetTitleOffset(1.3);

  return h_turnon;
}

//______________________________________________________
TCanvas* compare_turnons(std::vector<TH1F*> hists_num, TH1* h_denom, std::vector<std::string> labels, std::string name, int rebin) {

  // if (rebin > 1) {
  //   h_denom->Rebin(rebin);
  // }

  std::vector<int> colors;
  colors.push_back(kBlack);
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kGreen+1);
  colors.push_back(kMagenta);
  colors.push_back(kOrange);
  colors.push_back(kBlue-4);

  TCanvas* c = new TCanvas(Form("c_%s",name.c_str()),Form("c_%s",name.c_str()));
  c->SetGrid(1,1);

  // for HT (0.71,0.21,0.96,0.52)
  // for jets (0.58,0.21,0.89,0.44)

  float x1(0.71), y1(0.21), x2(0.96), y2(0.52);
  if (TString(h_denom->GetName()).Contains("genjet2")) {
    x1 = 0.58; y1 = 0.21; x2 = 0.89; y2 = 0.44;
  }

  TLegend* leg = new TLegend(x1,y1,x2,y2);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);

  for (unsigned int ihist = 0; ihist < hists_num.size(); ++ihist) {
    TH1F* h_turnon = turnon_plot(hists_num.at(ihist),h_denom,rebin);
    h_turnon->SetLineWidth(2);
    h_turnon->SetLineColor(colors.at(ihist));
    h_turnon->SetMarkerColor(colors.at(ihist));
    h_turnon->SetMarkerStyle(20);
    if (ihist == 0) h_turnon->Draw("pe");
    else h_turnon->Draw("pe same");
    leg->AddEntry(h_turnon,labels.at(ihist).c_str(),"pe");
    if (outfile) {
      outfile->cd();
      h_turnon->Write();
      fin->cd();
    }
  } // loop on hists

  leg->Draw("same");

  if (savePlots_) c->SaveAs(Form("%s/%s.eps",outdir_.Data(),c->GetName()));

  return c;
}
