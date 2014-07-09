#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include <vector>

TH1D* rate_plot(TH1* h_in, float norm = 1., int rebin = 1);
TCanvas* rate_vs_pu(TH2* h_in, float norm = 1., int rebin = 1);
float ScaleFactor(float nZeroBias, float nBunches);
TCanvas* compare_rates(std::vector<TH1D*> hists, std::vector<std::string> labels, std::string name, int nBunches, int rebin = 1);
void dump_table(std::vector<TH1D*> hists, std::vector<std::string> labels, int nBunches);

//______________________________________________________
void make_rates(TString input_file = "histos_minbias_rates_PU40bx25_met_all.root", int bx = 25, bool doRateVsPU = false) {

  // bunches for 50ns spacing
  int nBunches = 1368;
  // bunches for 25ns spacing
  if (bx == 25) {
    nBunches = 2508; //2508 is what agreed with TSG for # bunches
  }

  TFile* f = new TFile(input_file);

  // --- dijets ----
  int rebin_jets = 4;

  TCanvas* c_dijet = new TCanvas("c_dijet","c_dijet");
  c_dijet->SetLogy();
  c_dijet->SetGrid(1,1);
  TH1D* h_dijet_eta30 = f->Get("h_jet2_pt_eta30");
  float norm = ScaleFactor(h_dijet_eta30->GetEntries(),nBunches);
  TH1D* h_dijet_eta30_rate = rate_plot(h_dijet_eta30,norm,rebin_jets);
  h_dijet_eta30_rate->SetLineColor(kRed);
  h_dijet_eta30_rate->SetMarkerColor(kRed);
  h_dijet_eta30_rate->SetMarkerStyle(20);
  h_dijet_eta30_rate->GetXaxis()->SetTitle("Dijet p_{T} [GeV]");
  h_dijet_eta30_rate->Draw("pe");

  TH1D* h_dijet_eta22 = f->Get("h_jet2_pt_eta22");
  TH1D* h_dijet_eta22_rate = rate_plot(h_dijet_eta22,norm,rebin_jets);
  h_dijet_eta22_rate->SetLineColor(kBlue);
  h_dijet_eta22_rate->SetMarkerColor(kBlue);
  h_dijet_eta22_rate->SetMarkerStyle(20);
  h_dijet_eta22_rate->Draw("pe same");

  TH2D* h_dijet_eta30_vs_dphi = f->Get("h_jet12_dphi_vs_jet2_pt_eta30");
  TH1D* h_dijet_eta30_dphi180 = h_dijet_eta30_vs_dphi->ProjectionX("h_dijet_eta30_dphi180",7,94,"e");
  TH1D* h_dijet_eta30_dphi180_rate = rate_plot(h_dijet_eta30_dphi180,norm,rebin_jets);
  h_dijet_eta30_dphi180_rate->SetLineColor(kGreen+2);
  h_dijet_eta30_dphi180_rate->SetMarkerColor(kGreen+2);
  h_dijet_eta30_dphi180_rate->SetMarkerStyle(20);
  h_dijet_eta30_dphi180_rate->Draw("pe same");

  TH1D* h_dijet_eta30_dphi160 = h_dijet_eta30_vs_dphi->ProjectionX("h_dijet_eta30_dphi160",12,89,"e");
  TH1D* h_dijet_eta30_dphi160_rate = rate_plot(h_dijet_eta30_dphi160,norm,rebin_jets);
  h_dijet_eta30_dphi160_rate->SetLineColor(kMagenta);
  h_dijet_eta30_dphi160_rate->SetMarkerColor(kMagenta);
  h_dijet_eta30_dphi160_rate->SetMarkerStyle(20);
  h_dijet_eta30_dphi160_rate->Draw("pe same");

  TLegend* leg_dijet = new TLegend(0.48,0.63,0.89,0.83);
  leg_dijet->SetLineColor(0);
  leg_dijet->SetFillColor(0);
  leg_dijet->SetShadowColor(0);

  leg_dijet->AddEntry(h_dijet_eta30_rate, "|#eta| < 3.0" ,"pe");
  leg_dijet->AddEntry(h_dijet_eta22_rate, "|#eta| < 2.2" ,"pe");
  leg_dijet->AddEntry(h_dijet_eta30_dphi180_rate, "|#eta| < 3.0, |#Delta#phi| < 9 bins" ,"pe");
  leg_dijet->AddEntry(h_dijet_eta30_dphi160_rate, "|#eta| < 3.0, |#Delta#phi| < 8 bins" ,"pe");
  leg_dijet->Draw("same");

  // --- trijets ----

  TCanvas* c_trijet = new TCanvas("c_trijet","c_trijet");
  c_trijet->SetLogy();
  c_trijet->SetGrid(1,1);
  TH1D* h_trijet_etaall = f->Get("h_jet3_pt");
  TH1D* h_trijet_etaall_rate = rate_plot(h_trijet_etaall,norm,rebin_jets);
  h_trijet_etaall_rate->SetMarkerStyle(20);
  h_trijet_etaall_rate->GetXaxis()->SetTitle("Trijet p_{T} [GeV]");
  h_trijet_etaall_rate->Draw("pe");

  TH1D* h_trijet_eta30 = f->Get("h_jet3_pt_eta30");
  TH1D* h_trijet_eta30_rate = rate_plot(h_trijet_eta30,norm,rebin_jets);
  h_trijet_eta30_rate->SetLineColor(kRed);
  h_trijet_eta30_rate->SetMarkerColor(kRed);
  h_trijet_eta30_rate->SetMarkerStyle(20);
  h_trijet_eta30_rate->Draw("pe same");

  TH1D* h_trijet_eta22 = f->Get("h_jet3_pt_eta22");
  TH1D* h_trijet_eta22_rate = rate_plot(h_trijet_eta22,norm,rebin_jets);
  h_trijet_eta22_rate->SetLineColor(kBlue);
  h_trijet_eta22_rate->SetMarkerColor(kBlue);
  h_trijet_eta22_rate->SetMarkerStyle(20);
  h_trijet_eta22_rate->Draw("pe same");

  TLegend* leg_trijet = new TLegend(0.48,0.63,0.89,0.83);
  leg_trijet->SetLineColor(0);
  leg_trijet->SetFillColor(0);
  leg_trijet->SetShadowColor(0);

  leg_trijet->AddEntry(h_trijet_etaall_rate, "all #eta" ,"pe");
  leg_trijet->AddEntry(h_trijet_eta30_rate, "|#eta| < 3.0" ,"pe");
  leg_trijet->AddEntry(h_trijet_eta22_rate, "|#eta| < 2.2" ,"pe");
  leg_trijet->Draw("same");


  // --- quadjets ----

  TCanvas* c_quadjet = new TCanvas("c_quadjet","c_quadjet");
  c_quadjet->SetLogy();
  c_quadjet->SetGrid(1,1);
  TH1D* h_quadjet_etaall = f->Get("h_jet4_pt");
  TH1D* h_quadjet_etaall_rate = rate_plot(h_quadjet_etaall,norm,rebin_jets);
  h_quadjet_etaall_rate->SetMarkerStyle(20);
  h_quadjet_etaall_rate->GetXaxis()->SetTitle("Quadjet p_{T} [GeV]");
  h_quadjet_etaall_rate->Draw("pe");

  TH1D* h_quadjet_eta30 = f->Get("h_jet4_pt_eta30");
  TH1D* h_quadjet_eta30_rate = rate_plot(h_quadjet_eta30,norm,rebin_jets);
  h_quadjet_eta30_rate->SetLineColor(kRed);
  h_quadjet_eta30_rate->SetMarkerColor(kRed);
  h_quadjet_eta30_rate->SetMarkerStyle(20);
  h_quadjet_eta30_rate->Draw("pe same");

  TH1D* h_quadjet_eta22 = f->Get("h_jet4_pt_eta22");
  TH1D* h_quadjet_eta22_rate = rate_plot(h_quadjet_eta22,norm,rebin_jets);
  h_quadjet_eta22_rate->SetLineColor(kBlue);
  h_quadjet_eta22_rate->SetMarkerColor(kBlue);
  h_quadjet_eta22_rate->SetMarkerStyle(20);
  h_quadjet_eta22_rate->Draw("pe same");

  TLegend* leg_quadjet = new TLegend(0.48,0.63,0.89,0.83);
  leg_quadjet->SetLineColor(0);
  leg_quadjet->SetFillColor(0);
  leg_quadjet->SetShadowColor(0);

  leg_quadjet->AddEntry(h_quadjet_etaall_rate, "all #eta" ,"pe");
  leg_quadjet->AddEntry(h_quadjet_eta30_rate, "|#eta| < 3.0" ,"pe");
  leg_quadjet->AddEntry(h_quadjet_eta22_rate, "|#eta| < 2.2" ,"pe");
  leg_quadjet->Draw("same");

  // --- HT110 + jets

  TCanvas* c_ht110jets = new TCanvas("c_ht110jets","c_ht110jets");
  c_ht110jets->SetLogy();
  c_ht110jets->SetGrid(1,1);

  TH2D* h_HT_vs_jet1pt_eta30 = f->Get("h_HT_vs_jet1pt_eta30");
  TH1D* h_jet1_ht110_eta30 = h_HT_vs_jet1pt_eta30->ProjectionX("h_jet1_ht110_eta30",111,1000,"e");
  TH1D* h_jet1_ht110_eta30_rate = rate_plot(h_jet1_ht110_eta30,norm,rebin_jets);
  h_jet1_ht110_eta30_rate->SetLineColor(kBlack);
  h_jet1_ht110_eta30_rate->SetMarkerColor(kBlack);
  h_jet1_ht110_eta30_rate->SetMarkerStyle(20);
  h_jet1_ht110_eta30_rate->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h_jet1_ht110_eta30_rate->Draw("pe");

  TH2D* h_HT_vs_jet2pt_eta30 = f->Get("h_HT_vs_jet2pt_eta30");
  TH1D* h_jet2_ht110_eta30 = h_HT_vs_jet2pt_eta30->ProjectionX("h_jet2_ht110_eta30",111,1000,"e");
  TH1D* h_jet2_ht110_eta30_rate = rate_plot(h_jet2_ht110_eta30,norm,rebin_jets);
  h_jet2_ht110_eta30_rate->SetLineColor(kRed);
  h_jet2_ht110_eta30_rate->SetMarkerColor(kRed);
  h_jet2_ht110_eta30_rate->SetMarkerStyle(20);
  h_jet2_ht110_eta30_rate->Draw("pe same");

  TH2D* h_HT_vs_jet2pt_eta30_dphiT = f->Get("h_HT_vs_jet2pt_eta30_dphiT");
  TH1D* h_jet2_ht110_eta30_dphiT = h_HT_vs_jet2pt_eta30_dphiT->ProjectionX("h_jet2_ht110_eta30_dphiT",111,1000,"e");
  TH1D* h_jet2_ht110_eta30_dphiT_rate = rate_plot(h_jet2_ht110_eta30_dphiT,norm,rebin_jets);
  h_jet2_ht110_eta30_dphiT_rate->SetLineColor(kMagenta);
  h_jet2_ht110_eta30_dphiT_rate->SetMarkerColor(kMagenta);
  h_jet2_ht110_eta30_dphiT_rate->SetMarkerStyle(20);
  h_jet2_ht110_eta30_dphiT_rate->Draw("pe same");

  TH2D* h_HT_vs_jet3pt_eta30 = f->Get("h_HT_vs_jet3pt_eta30");
  TH1D* h_jet3_ht110_eta30 = h_HT_vs_jet3pt_eta30->ProjectionX("h_jet3_ht110_eta30",111,1000,"e");
  TH1D* h_jet3_ht110_eta30_rate = rate_plot(h_jet3_ht110_eta30,norm,rebin_jets);
  h_jet3_ht110_eta30_rate->SetLineColor(kBlue);
  h_jet3_ht110_eta30_rate->SetMarkerColor(kBlue);
  h_jet3_ht110_eta30_rate->SetMarkerStyle(20);
  h_jet3_ht110_eta30_rate->Draw("pe same");

  TH2D* h_HT_vs_jet4pt_eta30 = f->Get("h_HT_vs_jet4pt_eta30");
  TH1D* h_jet4_ht110_eta30 = h_HT_vs_jet4pt_eta30->ProjectionX("h_jet4_ht110_eta30",111,1000,"e");
  TH1D* h_jet4_ht110_eta30_rate = rate_plot(h_jet4_ht110_eta30,norm,rebin_jets);
  h_jet4_ht110_eta30_rate->SetLineColor(kGreen+2);
  h_jet4_ht110_eta30_rate->SetMarkerColor(kGreen+2);
  h_jet4_ht110_eta30_rate->SetMarkerStyle(20);
  h_jet4_ht110_eta30_rate->Draw("pe same");

  TLegend* leg_ht110jets = new TLegend(0.48,0.63,0.89,0.83);
  leg_ht110jets->SetLineColor(0);
  leg_ht110jets->SetFillColor(0);
  leg_ht110jets->SetShadowColor(0);

  leg_ht110jets->AddEntry(h_jet1_ht110_eta30_rate, "Leading jet cut" ,"pe");
  leg_ht110jets->AddEntry(h_jet2_ht110_eta30_rate, "2nd jet cut" ,"pe");
  leg_ht110jets->AddEntry(h_jet2_ht110_eta30_dphiT_rate, "2nd jet cut, |#Delta#phi| < 8 bins" ,"pe");
  leg_ht110jets->AddEntry(h_jet3_ht110_eta30_rate, "3rd jet cut" ,"pe");
  leg_ht110jets->AddEntry(h_jet4_ht110_eta30_rate, "4th jet cut" ,"pe");
  leg_ht110jets->Draw("same");


  // --- HT125 + jets

  TCanvas* c_ht125jets = new TCanvas("c_ht125jets","c_ht125jets");
  c_ht125jets->SetLogy();
  c_ht125jets->SetGrid(1,1);

  TH1D* h_jet1_ht125_eta30 = h_HT_vs_jet1pt_eta30->ProjectionX("h_jet1_ht125_eta30",126,1000,"e");
  TH1D* h_jet1_ht125_eta30_rate = rate_plot(h_jet1_ht125_eta30,norm,rebin_jets);
  h_jet1_ht125_eta30_rate->SetLineColor(kBlack);
  h_jet1_ht125_eta30_rate->SetMarkerColor(kBlack);
  h_jet1_ht125_eta30_rate->SetMarkerStyle(20);
  h_jet1_ht125_eta30_rate->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h_jet1_ht125_eta30_rate->Draw("pe");

  TH1D* h_jet2_ht125_eta30 = h_HT_vs_jet2pt_eta30->ProjectionX("h_jet2_ht125_eta30",126,1000,"e");
  TH1D* h_jet2_ht125_eta30_rate = rate_plot(h_jet2_ht125_eta30,norm,rebin_jets);
  h_jet2_ht125_eta30_rate->SetLineColor(kRed);
  h_jet2_ht125_eta30_rate->SetMarkerColor(kRed);
  h_jet2_ht125_eta30_rate->SetMarkerStyle(20);
  h_jet2_ht125_eta30_rate->Draw("pe same");

  TH1D* h_jet2_ht125_eta30_dphiT = h_HT_vs_jet2pt_eta30_dphiT->ProjectionX("h_jet2_ht125_eta30_dphiT",126,1000,"e");
  TH1D* h_jet2_ht125_eta30_dphiT_rate = rate_plot(h_jet2_ht125_eta30_dphiT,norm,rebin_jets);
  h_jet2_ht125_eta30_dphiT_rate->SetLineColor(kMagenta);
  h_jet2_ht125_eta30_dphiT_rate->SetMarkerColor(kMagenta);
  h_jet2_ht125_eta30_dphiT_rate->SetMarkerStyle(20);
  h_jet2_ht125_eta30_dphiT_rate->Draw("pe same");

  TH1D* h_jet3_ht125_eta30 = h_HT_vs_jet3pt_eta30->ProjectionX("h_jet3_ht125_eta30",126,1000,"e");
  TH1D* h_jet3_ht125_eta30_rate = rate_plot(h_jet3_ht125_eta30,norm,rebin_jets);
  h_jet3_ht125_eta30_rate->SetLineColor(kBlue);
  h_jet3_ht125_eta30_rate->SetMarkerColor(kBlue);
  h_jet3_ht125_eta30_rate->SetMarkerStyle(20);
  h_jet3_ht125_eta30_rate->Draw("pe same");

  TH1D* h_jet4_ht125_eta30 = h_HT_vs_jet4pt_eta30->ProjectionX("h_jet4_ht125_eta30",126,1000,"e");
  TH1D* h_jet4_ht125_eta30_rate = rate_plot(h_jet4_ht125_eta30,norm,rebin_jets);
  h_jet4_ht125_eta30_rate->SetLineColor(kGreen+2);
  h_jet4_ht125_eta30_rate->SetMarkerColor(kGreen+2);
  h_jet4_ht125_eta30_rate->SetMarkerStyle(20);
  h_jet4_ht125_eta30_rate->Draw("pe same");

  TLegend* leg_ht125jets = new TLegend(0.48,0.63,0.89,0.83);
  leg_ht125jets->SetLineColor(0);
  leg_ht125jets->SetFillColor(0);
  leg_ht125jets->SetShadowColor(0);

  leg_ht125jets->AddEntry(h_jet1_ht125_eta30_rate, "Leading jet cut" ,"pe");
  leg_ht125jets->AddEntry(h_jet2_ht125_eta30_rate, "2nd jet cut" ,"pe");
  leg_ht125jets->AddEntry(h_jet2_ht125_eta30_dphiT_rate, "2nd jet cut, |#Delta#phi| < 8 bins" ,"pe");
  leg_ht125jets->AddEntry(h_jet3_ht125_eta30_rate, "3rd jet cut" ,"pe");
  leg_ht125jets->AddEntry(h_jet4_ht125_eta30_rate, "4th jet cut" ,"pe");
  leg_ht125jets->Draw("same");


  // --- HT110, eta < 2.2, reget20 + jets

  TCanvas* c_ht110_eta22_reget20_jets = new TCanvas("c_ht110_eta22_reget20_jets","c_ht110_eta22_reget20_jets");
  c_ht110_eta22_reget20_jets->SetLogy();
  c_ht110_eta22_reget20_jets->SetGrid(1,1);

  TH2D* h_HT_eta22_reget20_vs_jet1pt_eta22 = f->Get("h_HT_eta22_reget20_vs_jet1pt_eta22");
  TH1D* h_jet1_ht110_eta22_reget20 = h_HT_eta22_reget20_vs_jet1pt_eta22->ProjectionX("h_jet1_ht110_eta22_reget20",111,1000,"e");
  TH1D* h_jet1_ht110_eta22_reget20_rate = rate_plot(h_jet1_ht110_eta22_reget20,norm,rebin_jets);
  h_jet1_ht110_eta22_reget20_rate->SetLineColor(kBlack);
  h_jet1_ht110_eta22_reget20_rate->SetMarkerColor(kBlack);
  h_jet1_ht110_eta22_reget20_rate->SetMarkerStyle(20);
  h_jet1_ht110_eta22_reget20_rate->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
  h_jet1_ht110_eta22_reget20_rate->Draw("pe");

  TH2D* h_HT_eta22_reget20_vs_jet2pt_eta22 = f->Get("h_HT_eta22_reget20_vs_jet2pt_eta22");
  TH1D* h_jet2_ht110_eta22_reget20 = h_HT_eta22_reget20_vs_jet2pt_eta22->ProjectionX("h_jet2_ht110_eta22_reget20",111,1000,"e");
  TH1D* h_jet2_ht110_eta22_reget20_rate = rate_plot(h_jet2_ht110_eta22_reget20,norm,rebin_jets);
  h_jet2_ht110_eta22_reget20_rate->SetLineColor(kRed);
  h_jet2_ht110_eta22_reget20_rate->SetMarkerColor(kRed);
  h_jet2_ht110_eta22_reget20_rate->SetMarkerStyle(20);
  h_jet2_ht110_eta22_reget20_rate->Draw("pe same");

  TH2D* h_HT_eta22_reget20_vs_jet2pt_eta22_dphiT = f->Get("h_HT_eta22_reget20_vs_jet2pt_eta22_dphiT");
  TH1D* h_jet2_ht110_eta22_reget20_dphiT = h_HT_eta22_reget20_vs_jet2pt_eta22_dphiT->ProjectionX("h_jet2_ht110_eta22_reget20_dphiT",111,1000,"e");
  TH1D* h_jet2_ht110_eta22_reget20_dphiT_rate = rate_plot(h_jet2_ht110_eta22_reget20_dphiT,norm,rebin_jets);
  h_jet2_ht110_eta22_reget20_dphiT_rate->SetLineColor(kMagenta);
  h_jet2_ht110_eta22_reget20_dphiT_rate->SetMarkerColor(kMagenta);
  h_jet2_ht110_eta22_reget20_dphiT_rate->SetMarkerStyle(20);
  h_jet2_ht110_eta22_reget20_dphiT_rate->Draw("pe same");

  TH2D* h_HT_eta22_reget20_vs_jet3pt_eta22 = f->Get("h_HT_eta22_reget20_vs_jet3pt_eta22");
  TH1D* h_jet3_ht110_eta22_reget20 = h_HT_eta22_reget20_vs_jet3pt_eta22->ProjectionX("h_jet3_ht110_eta22_reget20",111,1000,"e");
  TH1D* h_jet3_ht110_eta22_reget20_rate = rate_plot(h_jet3_ht110_eta22_reget20,norm,rebin_jets);
  h_jet3_ht110_eta22_reget20_rate->SetLineColor(kBlue);
  h_jet3_ht110_eta22_reget20_rate->SetMarkerColor(kBlue);
  h_jet3_ht110_eta22_reget20_rate->SetMarkerStyle(20);
  h_jet3_ht110_eta22_reget20_rate->Draw("pe same");

  TH2D* h_HT_eta22_reget20_vs_jet4pt_eta22 = f->Get("h_HT_eta22_reget20_vs_jet4pt_eta22");
  TH1D* h_jet4_ht110_eta22_reget20 = h_HT_eta22_reget20_vs_jet4pt_eta22->ProjectionX("h_jet4_ht110_eta22_reget20",111,1000,"e");
  TH1D* h_jet4_ht110_eta22_reget20_rate = rate_plot(h_jet4_ht110_eta22_reget20,norm,rebin_jets);
  h_jet4_ht110_eta22_reget20_rate->SetLineColor(kGreen+2);
  h_jet4_ht110_eta22_reget20_rate->SetMarkerColor(kGreen+2);
  h_jet4_ht110_eta22_reget20_rate->SetMarkerStyle(20);
  h_jet4_ht110_eta22_reget20_rate->Draw("pe same");



  TLegend* leg_ht110_eta22_reget20_jets = new TLegend(0.48,0.63,0.89,0.83);
  leg_ht110_eta22_reget20_jets->SetLineColor(0);
  leg_ht110_eta22_reget20_jets->SetFillColor(0);
  leg_ht110_eta22_reget20_jets->SetShadowColor(0);

  leg_ht110_eta22_reget20_jets->AddEntry(h_jet1_ht110_eta22_reget20_rate, "Leading jet cut" ,"pe");
  leg_ht110_eta22_reget20_jets->AddEntry(h_jet2_ht110_eta22_reget20_rate, "2nd jet cut" ,"pe");
  leg_ht110_eta22_reget20_jets->AddEntry(h_jet2_ht110_eta22_reget20_dphiT_rate, "2nd jet cut, |#Delta#phi| < 8 bins" ,"pe");
  leg_ht110_eta22_reget20_jets->AddEntry(h_jet3_ht110_eta22_reget20_rate, "3rd jet cut" ,"pe");
  leg_ht110_eta22_reget20_jets->AddEntry(h_jet4_ht110_eta22_reget20_rate, "4th jet cut" ,"pe");
  leg_ht110_eta22_reget20_jets->Draw("same");


  std::vector<std::string> labels;
  labels.push_back("Region ET > 7 GeV");
  labels.push_back("Region ET > 10 GeV");
  labels.push_back("Region ET > 15 GeV");
  labels.push_back("Region ET > 20 GeV");

  int rebin_ht = 5;
 
  // ---- HT, eta 3.0

  std::vector<TH1D*> hists_eta30;
  hists_eta30.push_back((TH1D*)f->Get("h_HT_reget7"));
  hists_eta30.push_back((TH1D*)f->Get("h_HT_reget10"));
  hists_eta30.push_back((TH1D*)f->Get("h_HT_reget15"));
  hists_eta30.push_back((TH1D*)f->Get("h_HT_reget20"));

  TCanvas* c_ht_eta30 = compare_rates(hists_eta30,labels,"ht_eta30",nBunches,rebin_ht);

  std::cout << "--- table for |eta| < 3.0:" << std::endl;
  dump_table(hists_eta30,labels,nBunches);
  std::cout << std::endl;

  // ---- HT, eta 2.2

  std::vector<TH1D*> hists_eta22;
  hists_eta22.push_back((TH1D*)f->Get("h_HT_eta22_reget7"));
  hists_eta22.push_back((TH1D*)f->Get("h_HT_eta22_reget10"));
  hists_eta22.push_back((TH1D*)f->Get("h_HT_eta22_reget15"));
  hists_eta22.push_back((TH1D*)f->Get("h_HT_eta22_reget20"));

  TCanvas* c_ht_eta22 = compare_rates(hists_eta22,labels,"ht_eta22",nBunches,rebin_ht);

  std::cout << "--- table for |eta| < 2.2:" << std::endl;
  dump_table(hists_eta22,labels,nBunches);
  std::cout << std::endl;

  if (doRateVsPU) {
    TH2D* h_ht_vs_nPU = f->Get("h_caloHT2015_vs_nPU");
    TCanvas* c_ht_vs_nPU = rate_vs_pu(h_ht_vs_nPU);
  } // if doRateVsPU

  gPad->Modified();

  return;
}

//______________________________________________________
TH1D* rate_plot(TH1* h_in, float norm, int rebin) {

  //  h_in->Sumw2();
  if (rebin > 1) h_in->Rebin(rebin);
  TH1D* h_rate = h_in->Clone(Form("%s_rate",h_in->GetName()));

  for (unsigned int ibin = 1; ibin < h_in->GetNbinsX(); ++ibin) {
    Double_t val,err;
    val = h_in->IntegralAndError(ibin,-1,err);
    //    std::cout << "bin, val, err: " << ibin << ", " << val << ", " << err << std::endl;
    h_rate->SetBinContent(ibin,val);
    h_rate->SetBinError(ibin,err);
  }

  h_rate->Scale(norm);
  h_rate->GetXaxis()->SetTitle("H_{T} [GeV]");
  h_rate->GetYaxis()->SetTitle("Rate [Hz]");
  h_rate->SetMarkerStyle(20);

  return h_rate;
}

//______________________________________________________
TCanvas* rate_vs_pu(TH2* h_in, float norm, int rebin) {
  TString cname = Form("c_%s",h_in->GetName());
  TCanvas* c = new TCanvas(cname,cname);
  c->SetLogy();
  c->SetGrid(1,1);

  const unsigned int minPU = 20;
  const unsigned int maxPU = 50;
  const unsigned int stepPU = 5;

  int color = 2;
  for (unsigned int ipu = minPU; ipu < maxPU; ipu += stepPU) {
    TString hname = Form("%s_%dPU",h_in->GetName(),ipu);
    TH1D* h_ht_proj = h_in->ProjectionY(hname,ipu+1,ipu+1+stepPU,"e");
    TH1D* h_ht_proj_rate = rate_plot(h_ht_proj,norm);
    h_ht_proj_rate->SetLineColor(color);
    h_ht_proj_rate->SetMarkerColor(color);
    h_ht_proj_rate->SetMarkerStyle(20);
    if (ipu == minPU) {
      h_ht_proj_rate->Draw("pe");
    } else {
      h_ht_proj_rate->Draw("pe same");
    }
    ++color;
  } // loop over PU bins

  return c;
}

//______________________________________________________
// scale factor computed w.r.t. ZeroBias rate fratcion and # bunches 
// -- taken from BasicRatePlots.C
float ScaleFactor(float nZeroBias, float nBunches) {

  float scal = 11246.; // ZB per bunch in Hz
  scal /= nZeroBias;
  scal *= nBunches;

  return scal;
}

//______________________________________________________
TCanvas* compare_rates(std::vector<TH1D*> hists, std::vector<std::string> labels, std::string name, int nBunches, int rebin) {

  std::vector<int> colors;
  colors.push_back(kRed);
  colors.push_back(kBlue);
  colors.push_back(kGreen+1);
  colors.push_back(kMagenta);

  TCanvas* c = new TCanvas(Form("c_%s",name.c_str()),Form("c_%s",name.c_str()));
  c->SetLogy();
  c->SetGrid(1,1);

  TLegend* leg = new TLegend(0.48,0.63,0.89,0.83);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetShadowColor(0);

  for (unsigned int ihist = 0; ihist < hists.size(); ++ihist) {
    float norm = ScaleFactor(hists.at(ihist)->GetEntries(),nBunches);
    TH1D* h_rate = rate_plot(hists.at(ihist),norm,rebin);
    h_rate->SetLineColor(colors.at(ihist));
    h_rate->SetMarkerColor(colors.at(ihist));
    if (ihist == 0) h_rate->Draw("pe");
    else h_rate->Draw("pe same");
    leg->AddEntry(h_rate,labels.at(ihist).c_str(),"pe");
  } // loop on hists

  leg->Draw("same");

  return c;
}


//______________________________________________________
void dump_table(std::vector<TH1D*> hists, std::vector<std::string> labels, int nBunches) {

  std::vector<int> cuts;
  cuts.push_back(110);
  cuts.push_back(125);
  cuts.push_back(150);
  cuts.push_back(175);
  cuts.push_back(200);

  std::cout << " Cuts ";
  for (unsigned int icut = 0; icut < cuts.size(); ++icut) {
    std::cout << "& HT" << cuts.at(icut) << " ";
  }
  std::cout << " \\\\" << std::endl;

  for (unsigned int ihist = 0; ihist < hists.size(); ++ihist) {
    std::cout << " " << labels.at(ihist) << " ";
    for (unsigned int icut = 0; icut < cuts.size(); ++icut) {
      float norm = ScaleFactor(hists.at(ihist)->GetEntries(),nBunches);
      int bin = hists.at(ihist)->FindBin(cuts.at(icut)+1);
      double rate = hists.at(ihist)->Integral(bin,-1) * norm / 1000.;
      std::cout << "& " << Form("%.1f",rate) << " ";
    } // loop on cuts
    std::cout << " \\\\" << std::endl;
  } // loop on hists

  return;
}

