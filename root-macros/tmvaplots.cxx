#include <iostream>
#include <vector>
#include <list>


///////////////////////
// CLASS DEFINITIONS //
///////////////////////

class inputVariable {

public:
  TString name;
  int nbins;
  double low, high;
  void set_values(TString, int, double, double);
};

void inputVariable::set_values(TString name_, int nbins_, double low_, double high_) {
  name = name_;
  nbins = nbins_;
  low = low_;
  high = high_;
}


///////////////
// FUNCTIONS //
///////////////

void plotInputVariables (TString filepath, TString nbtags) {
   
  TFile* rootfile = TFile::Open(filepath, "READ");
  TString histpath = "Method_BDT/BDT/";

  //vector<inputVariable> variables;
  //variables.push_back();

  vector<TString> variableNames = {"deltaphi_bottomlepton", "deltaphi_bottomtoptag", "pseudotop_mass", "ptbalance", "wass_pt", "lepton_pt", "lepton_eta", "deltaphi_leptonnextjet"};
  if(nbtags == "2b") variableNames.push_back("n_btags");

  for (int i = 0; i < variableNames.size(); i++) {

    TH1F* hist_sig = (TH1F*)rootfile->Get(histpath+"tmva_"+variableNames.at(i)+"__Signal");
    TH1F* hist_bkg = (TH1F*)rootfile->Get(histpath+"tmva_"+variableNames.at(i)+"__Background");

    double binwidth = hist_sig->GetBinWidth(1);

    for (auto hist : {hist_sig, hist_bkg}) {

      hist->SetTitle("");
      TString yAxisTitle = hist->GetYaxis()->GetTitle();
      yAxisTitle.ReplaceAll("dN","(1/N) dN");
      hist->GetYaxis()->SetTitle(yAxisTitle);
      hist->GetYaxis()->SetTitleSize(0.07);
      hist->GetXaxis()->SetTitleSize(0.07);
      hist->GetYaxis()->SetLabelSize(0.07);
      hist->GetXaxis()->SetLabelSize(0.07);
      if(variableNames.at(i) == "wass_pt") hist->GetXaxis()->SetNdivisions(805);
      if(variableNames.at(i) == "lepton_eta") hist->GetXaxis()->SetNdivisions(1005);
      if(variableNames.at(i) == "pseudotop_mass") hist->GetXaxis()->SetNdivisions(1004);
      
      double integral = hist->Integral();
      hist->Scale(1./integral/binwidth);

      hist->SetLineWidth(3);
    }

    double hist_sig_max = hist_sig->GetMaximum();
    double hist_bkg_max = hist_bkg->GetMaximum();

    int cx = 450;
    int cy = 400;

    TCanvas* c = new TCanvas("c","",cx,cy);

    TPad* p = new TPad("p","pad",.0,0,1,1);
    Double_t margin = 0.1;
    p->SetMargin((60+cx-400)/((double)cx),20/((double)cx),0.18,0.02);
    p->Draw();
    p->cd();

    hist_sig->Draw("hist");
    hist_bkg->Draw("hist same");

    hist_sig->SetMaximum(1.1*max(hist_sig_max, hist_bkg_max));
    hist_sig->GetXaxis()->SetTitleOffset(1.2);
    hist_sig->GetYaxis()->SetTitleOffset(1.8);

    hist_sig->SetLineColor(597);
    hist_sig->SetFillColor(kBlue-9);
    hist_sig->SetLineStyle(0);

    hist_bkg->SetLineColor(810);
    hist_bkg->SetFillColor(810);
    hist_bkg->SetFillStyle(3004);

    gStyle->SetOptStat(kFALSE);

    TText *simulation = new TText(.2,.91,"Simulation");
    simulation->SetNDC();
    //simulation->SetTextAlign(12);
    simulation->SetTextFont(73);
    simulation->SetTextSize(24);
    //simulation->Draw();

    TLatex *lumi = new TLatex(.8,.91,"tW-channel (incl. decays)");
    lumi->SetNDC();
    lumi->SetTextAlign(31);
    lumi->SetTextFont(43);
    lumi->SetTextSize(24);
    //lumi->Draw();

    if(i == 0) {

      auto legend = new TLegend(0.2*cy/cx+(cx-400.)/400.+(nbtags == "1b" ? 0 : 0.25*cy/cx),0.65,0.55+(nbtags == "1b" ? 0 : 0.25*cy/cx),0.9);
      //legend->SetHeader("Training samples:");
      legend->AddEntry(hist_sig, "signal (tW)");
      legend->AddEntry(hist_bkg, "background (t#bar{t})");
      legend->Draw();
      legend->SetBorderSize(0);
      legend->SetTextSize(0.07);
    }

    c->SaveAs("plots/tmva_"+nbtags+"__"+variableNames.at(i)+".eps");
  }

}


double plotSignificance(TString channel, TString nbtags) {

  TH1F* eff_s = new TH1F("eff_s", "", 41, -1.025, 1.025);
  TH1F* eff_b = new TH1F("eff_b", "", 41, -1.025, 1.025);
  
  vector<TString> bkgs = {"SingleTop_tWch", "SingleTop_stch", "TTbar", "WJets", "DYJets", "Diboson"}; // QCD excluded since low statistics and negligible
  vector<TH1F*> hists_bkg;

  for (int i = 0; i < bkgs.size(); i++) {
    TString rootfile_path = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/uhh2.AnalysisModuleRunner.MC."+bkgs.at(i)+".root";
    TFile* rootfile = TFile::Open(rootfile_path, "READ");
    TString hist_path = "MVA_1t"+nbtags+"/mvaD";
    hists_bkg.push_back((TH1F*)rootfile->Get(hist_path));
    //rootfile->Close();
  }

  for (int i = 1; i <= 41; i++) {
    eff_s->SetBinContent(i,hists_bkg.at(0)->GetBinContent(i));
    double eff_b_bincontent = 0;
    for (int j = 1; j < bkgs.size(); j++) {
      eff_b_bincontent += hists_bkg.at(j)->GetBinContent(i);
    }
    eff_b->SetBinContent(i,eff_b_bincontent);
    //cout << eff_b_bincontent << endl;
  }

  //eff_s->Rebin(5);
  //eff_b->Rebin(5);

  double integral_s = eff_s->Integral();
  cout << "S = " << integral_s << endl;
  double integral_b = eff_b->Integral();
  cout << "B = " << integral_b << endl;

  TH1F* eff_s_save = (TH1F*)eff_s->Clone("eff_s_save");
  TH1F* eff_b_save = (TH1F*)eff_b->Clone("eff_b_save");

  int s = 0;
  for (auto hist : {eff_s, eff_b}) {
    double integral = hist->Integral();
    for (int i = 1; i <= 41; i++) {
      double efficiency = 0;
      for (int j = 1; j < i; j++) {
	if (s == 0) efficiency += eff_s_save->GetBinContent(j);
	else        efficiency += eff_b_save->GetBinContent(j);
      }
      efficiency = (integral-efficiency)/integral;
      hist->SetBinContent(i,efficiency);
    }
    s++;
  }
  
  double offset = 0.05;

  TCanvas* c = new TCanvas("canvas", "", 500, 500);

  c->cd();

  gStyle->SetOptStat(0);

  TPad* p1 = new TPad("p1", "", 0, .5+offset, 1, 1);
  p1->SetMargin(0.1+offset,0.1-offset,0.,(0.1-offset)/(0.5-offset));
  p1->Draw();
  p1->cd();
  p1->SetTicks(1,1);

  eff_s->Draw("l hist");
  eff_b->Draw("l hist same");
  eff_s->SetMinimum(0.01);
  eff_s->SetMaximum(1);
  eff_s->GetYaxis()->SetNdivisions(406);

  eff_s->SetLineColor(597);
  eff_s->SetLineWidth(2);
  eff_s->SetMarkerSize(0);
  eff_b->SetLineColor(810);
  eff_b->SetLineWidth(2);
  eff_b->SetLineStyle(2);
  eff_b->SetMarkerSize(0);
 
  eff_s->GetXaxis()->SetTitle("Cut value applied on BDT response");
  eff_s->GetYaxis()->SetTitle("Efficiency");
  eff_s->GetYaxis()->SetTitleSize(.1*2*(.5+offset));
  eff_s->GetYaxis()->SetLabelSize(.1*2*(.5+offset));
  eff_s->GetYaxis()->SetTitleOffset(0.6*2*(.5-offset));
  eff_s->GetYaxis()->SetTickSize(.02*2*(.5-offset));

  c->cd();

  TPad* p2 = new TPad("p2", "", 0, 0, 1, .5+offset);
  p2->SetMargin(0.1+offset,0.1-offset,(0.1+offset)/(0.5+offset),0.);
  p2->Draw();
  p2->cd();
  p2->SetTicks(1,1);

  TH1F* ssb = new TH1F("ssb", "", 41, -1.025, 1.025);
  TH1F* purity_s = new TH1F("purity_s", "", 41, -1.025, 1.025);

  double ssb_max = 0;
  double ssb_min = 100000;
  double significance = 0;
  double purity = 0;
  int ssb_max_i = 0;
  for (int i = 1; i <= 41; i++) {
    if (i == 41) {
      ssb->SetBinContent(i,0);
      purity_s->SetBinContent(i,1);
      break;
    }
    significance = integral_s*eff_s->GetBinContent(i)/sqrt(integral_s*eff_s->GetBinContent(i)+integral_b*eff_b->GetBinContent(i));
    purity = integral_s*eff_s->GetBinContent(i)/(integral_b*eff_b->GetBinContent(i));
    ssb->SetBinContent(i,significance);
    purity_s->SetBinContent(i,purity);
    if (significance > ssb_max) {
      ssb_max = significance;
      ssb_max_i = i;
    }
    if (significance < ssb_min) ssb_min = significance;
  }
  //eff_s->SetMaximum(ssb_max*1.3);
  //eff_s->SetMaximum(1);
  
  /*purity_s->Draw("l same");
  purity_s->SetLineColor(597);
  purity_s->SetLineWidth(2);
  purity_s->SetLineStyle(9);*/

  ssb->Draw("p l hist"); // try c istead of l for smooth line
  //ssb->Scale(1/ssb_max);

  ssb->SetMarkerStyle(8);
  ssb->SetMarkerSize(.6);

  ssb->SetLineColor(597);
  ssb->SetLineWidth(2);

  ssb->SetMaximum(ssb_max*1.15);
  ssb->SetMinimum(0);

  cout << ssb_min << ssb_max << endl;

  ssb->GetXaxis()->SetTitle("Cut value applied on BDT response");
  ssb->GetYaxis()->SetTitle("Signal significance"); //S/#sqrt{S+B}
  ssb->GetXaxis()->SetTitleSize(.1*2*(.5-offset));
  ssb->GetXaxis()->SetLabelSize(.1*2*(.5-offset));
  ssb->GetXaxis()->SetTitleOffset(1.2);
  ssb->GetYaxis()->SetLabelSize(.1*2*(.5-offset));
  ssb->GetYaxis()->SetTitleSize(.1*2*(.5-offset));
  ssb->GetYaxis()->SetTitleOffset(0.6*2*(.5+offset));
  ssb->GetYaxis()->SetTickSize(.02*2*(.5+offset));

  TLine *lineX = new TLine(ssb->GetBinCenter(ssb_max_i),0,ssb->GetBinCenter(ssb_max_i),100);
  lineX->SetLineStyle(2);
  lineX->SetLineWidth(2);
  lineX->Draw();

  TLine *lineY = new TLine(-1.025,ssb->GetBinContent(ssb_max_i),1.025,ssb->GetBinContent(ssb_max_i));
  lineY->SetLineStyle(2);
  lineY->SetLineWidth(2);
  lineY->Draw();

  TMarker *marker = new TMarker(ssb->GetBinCenter(ssb_max_i),ssb->GetBinContent(ssb_max_i),8);
  marker->Draw();

  p1->cd();

  auto legend = new TLegend(0.63,0.55,0.9,0.85);
  legend->AddEntry(eff_s, "signal (tW)");
  legend->AddEntry(eff_b, "all bkgs.");
  legend->Draw();
  legend->SetBorderSize(0);
  legend->SetTextSize(.08*2*(.5+offset)); // .1 as first factor is standard for all other texts in this plot

  c->SaveAs("plots/tmva_"+channel+"_"+nbtags+"__OptimalCutValue.eps");


  return ssb->GetBinCenter(ssb_max_i);
}


void Pal1()
{
   static Int_t  colors[50];
   static Bool_t initialized = kFALSE;
   Double_t Red[3]    = { 0.00, 0.00, 1.00};
   Double_t Green[3]  = { 0.00, 0.00, 0.00};
   Double_t Blue[3]   = { 1.00, 0.00, 0.00};
   Double_t Length[3] = { 0.00, 0.50, 1.00 };
   if(!initialized){
      Int_t FI = TColor::CreateGradientColorTable(3,Length,Red,Green,Blue,50);
      for (int i=0; i<50; i++) colors[i] = FI+i;
      initialized = kTRUE;
      return;
   }
   gStyle->SetPalette(50,colors);
}


void linearCorrelations(TString filepath, TString nbtags) {
   
  TFile* rootfile = TFile::Open(filepath, "READ");
  TH1F* hist_sig = (TH1F*)rootfile->Get("CorrelationMatrixS");
  TH1F* hist_bkg = (TH1F*)rootfile->Get("CorrelationMatrixB");
  //rootfile->Close();

  TString SB = "S";

  for(auto hist : {hist_sig, hist_bkg}) {

    TCanvas* c = new TCanvas("c","",700,500);

    c->cd();

    TPad* p = new TPad("p","",0,0,1,1);
    p->SetMargin(.1*7/5,.1*7/5,.12,.02);
    p->Draw();
    p->SetTicks(1,1);
    p->cd();

    //Pal1();
    gStyle->SetPalette(kTemperatureMap);

    hist->SetTitle("");
    hist->Draw("colz text");

    hist->SetMarkerColor(kBlack);

    c->SaveAs("plots/tmva_"+nbtags+"__LinearCorrelationMatrix"+SB+".eps");

    SB = "B";
  }

  //rootfile->Close();

}


void tmvaplots() {

  TString tmvafolder = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb/";

  for (auto nbtags : {"1b", "2b"}) {
    TString filepath = tmvafolder+"TMVAOutput_"+nbtags+".root";
    //plotInputVariables(filepath,nbtags);
    /*for (auto channel : {"Electron", "Muon"}) {
      double optimalCut = plotSignificance(channel,nbtags);
      cout << "Optimal BDT cut: " << optimalCut << endl;
      }*/
    linearCorrelations(filepath,nbtags);
  }
}
