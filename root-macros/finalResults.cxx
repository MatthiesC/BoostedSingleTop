#include <iostream>
#include <vector>

vector<double> recoEfficiency(TString channel) {
   
  TString input_path_rec = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/";
  TString input_path_gen = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop_tWch/";
  cout << "Make sure to already have run SFramePlotter over the input files to have the hadded tW sample (top+antitop)!" << endl;
  TFile* input_file_rec = TFile::Open(input_path_rec+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TFile* input_file_gen = TFile::Open(input_path_gen+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TH1F* hist_rec = (TH1F*)input_file_rec->Get("Matching_tW/topjet_pt_rebin");
  TString genhist;
  if(channel == "Electron") genhist = "nocuts/EleChannel_TopPt_GenCut";
  else if(channel == "Muon") genhist = "nocuts/MuoChannel_TopPt_GenCut";
  TH1F* hist_gen = (TH1F*)input_file_gen->Get(genhist);

  TCanvas* c = new TCanvas("c","",600,800);

  double splitter = 0.375; // decide where to put the border between the two pads
  
  //=======================//
  // UPPER RECO/GEN HISTOS //
  //=======================//

  TPad* upper = new TPad("upper","",0,splitter,1,1);

  upper->SetMargin(.15,.05,.0,.1*6/5);
  upper->Draw();
  upper->cd();
  hist_rec->Draw("hist e1");
  hist_gen->Draw("hist e1 same");

  gStyle->SetOptStat(kFALSE);
  
  hist_rec->SetTitle("");
  hist_rec->SetLineWidth(2);
  hist_rec->SetLineColor(597); // nice blue
  hist_rec->SetLineStyle(1);
  hist_rec->SetMinimum(0.4); // avoid showing the zero at y=0

  hist_gen->SetTitle("");
  hist_gen->SetLineWidth(2);
  hist_gen->SetLineColor(kBlack);
  hist_gen->SetLineStyle(2);

  hist_rec->GetYaxis()->SetTitle("Events");
  hist_rec->GetYaxis()->SetTitleOffset(1.5);
  hist_rec->GetYaxis()->SetTitleSize(0.05);
  hist_rec->GetYaxis()->SetLabelSize(0.05);
  hist_rec->GetXaxis()->SetTitle("");
  hist_rec->GetXaxis()->SetTitleOffset(1.4);
  hist_rec->SetMaximum(10.*hist_gen->GetBinContent(2));
  upper->SetLogy();

  //  c->cd();

  TText *simulation = new TText(.15,1-0.09*6/5,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();

  TLatex *lumi = new TLatex(.95,1-0.09*6/5,"35.9 fb^{-1} (13 TeV)");
  lumi->SetNDC();
  lumi->SetTextAlign(31);
  lumi->SetTextFont(43);
  lumi->SetTextSize(24);
  lumi->Draw();

  auto legend = new TLegend(0.65,0.52,0.92,0.80);
  legend->AddEntry(hist_rec, "reco. level");
  legend->AddEntry(hist_gen, "gen. level");
  if(channel == "Electron") legend->SetHeader("tW #rightarrow e+jets");
  else if(channel == "Muon") legend->SetHeader("tW #rightarrow #mu+jets");
  legend->Draw();
  legend->SetBorderSize(0);

  //==================//
  // EFFICIENCY HISTO //
  //==================//

  c->cd();
  TPad* lower = new TPad("lower", "",0,0,1,splitter);

  lower->SetMargin(.15,.05,.3,.0);
  lower->Draw();
  lower->cd();

  TH1F* hist_eff = (TH1F*)hist_gen->Clone();

  vector<double> result;

  hist_eff->SetBinContent(1,0);
  double uncertainty, bincontent;
  for(int i = 2; i < hist_eff->GetNbinsX()+1; i++) {
    bincontent = hist_rec->GetBinContent(i)/hist_gen->GetBinContent(i);
    hist_eff->SetBinContent(i,bincontent);
    uncertainty = bincontent*sqrt(pow(hist_rec->GetBinError(i)/hist_rec->GetBinContent(i),2)+pow(hist_gen->GetBinError(i)/hist_gen->GetBinContent(i),2));
    hist_eff->SetBinError(i,uncertainty);
    result.push_back(bincontent);
  }
  //hist_eff->SetMinimum(1.15);
  //hist_eff->SetMaximum(2.05);
  hist_eff->SetMaximum(0.14);
  hist_eff->GetYaxis()->SetNdivisions(503);
  hist_eff->GetYaxis()->SetLabelSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetLabelSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetTitleOffset(1.6);
  hist_eff->SetLineStyle(1);

  hist_eff->SetLineColor(kBlack);
  hist_eff->GetXaxis()->SetTitle("top-quark or top-jet p_{T} [GeV]");
  hist_eff->GetYaxis()->SetTitle("Efficiency");
  hist_eff->GetYaxis()->CenterTitle();
  hist_eff->GetYaxis()->SetTitleSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetTitleSize(0.05/(1-splitter));
  hist_eff->GetYaxis()->SetTitleOffset(0.9);

  hist_eff->Draw("e1");
 
  //==================//
  // SAVE AS EPS PLOT //
  //==================//

  c->cd();
  c->SaveAs("recoEfficiency_"+channel+".eps");

  return result;
}


void eventCountPlots(TString channel, vector<double> efficiencies) {

  TString input_path_measurement = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/theta-workdir/rootfiles/theta-output_stat.root";
  cout << "Make sure to run THETA to get your input histograms!" << endl;
  TFile* input_file_measurement = TFile::Open(input_path_measurement, "READ");
  TString input_hist_measurement = "";
  TH1F* hist_measured = (TH1F*)input_file_measurement->Get(input_hist_measurement);
}


void finalResults() {

  for(TString channel : {"Electron", "Muon"}) {
    auto efficiencies = recoEfficiency(channel);
    cout << "Results for the " << channel << " channel:" << endl;
    for(int i = 0; i < efficiencies.size(); i++) {
      cout << efficiencies.at(i) << endl;
    }
    eventCountPlots(channel, efficiencies);
  }
}
