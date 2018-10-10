#include <iostream>
#include <vector>

vector<double> recoEfficiency(TString channel, TString n_btags) {
   
  TString input_path_rec = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/";
  TString input_path_gen = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop_tWch/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded tW sample (top+antitop)!" << std::endl;
  TFile* input_file_rec = TFile::Open(input_path_rec+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TFile* input_file_gen = TFile::Open(input_path_gen+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TString input_hist_rec = "MVA_1t"+n_btags+"/nocuts_top_pt_rebin";
  TH1F* hist_rec = (TH1F*)input_file_rec->Get(input_hist_rec);
  TString genhist;
  if(channel == "Electron") genhist = "nocuts/EleChannel_TopPt_GenCut";
  else if(channel == "Muon") genhist = "nocuts/MuoChannel_TopPt_GenCut";
  TH1F* hist_gen_raw = (TH1F*)input_file_gen->Get(genhist);

  //double hf = 1./0.75;

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
  TH1F* hist_gen = (TH1F*)hist_rec->Clone();
  for(int i = 1; i < 5; i++) {
    hist_gen->SetBinContent(i,hist_gen_raw->GetBinContent(i+1));
    hist_gen->SetBinError(i,hist_gen_raw->GetBinError(i+1));
    }
  hist_gen->Draw("hist e1 same");

  gStyle->SetOptStat(kFALSE);
  
  hist_rec->SetTitle("");
  hist_rec->SetLineWidth(2);
  hist_rec->SetLineColor(597); // nice blue
  hist_rec->SetLineStyle(1);
  hist_rec->SetMinimum(4.0); // avoid showing the zero at y=0

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
  hist_rec->GetXaxis()->SetNdivisions(507);
  hist_rec->SetMaximum(10.*hist_gen->GetBinContent(1));
  upper->SetLogy();

  //  c->cd();

  TText *simulation = new TText(.15,1-0.09*6/5,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();

  TString category_1 = "tW #rightarrow ";
  TString category_2 = (channel == "Electron" ? "e+jets " : "#mu+jets ");
  TString category_3 = (n_btags == "1b" ? "(1 b-tag)" : "(#geq 2 b-tags)");

  TLatex *category = new TLatex(.2,1-0.09*6/5-0.1,category_1+category_2+category_3);
  category->SetNDC();
  //category->SetTextAlign(12);
  category->SetTextFont(43);
  category->SetTextSize(24);
  category->Draw();

  TLatex *lumi = new TLatex(.95,1-0.09*6/5,"35.9 fb^{-1} (13 TeV)");
  lumi->SetNDC();
  lumi->SetTextAlign(31);
  lumi->SetTextFont(43);
  lumi->SetTextSize(24);
  lumi->Draw();

  auto legend = new TLegend(0.65,0.52,0.92,0.80);
  legend->AddEntry(hist_rec, "reco. level");
  legend->AddEntry(hist_gen, "gen. level");
  legend->SetHeader("Selection type:");
  /*if(channel == "Electron" && n_btags == "1b") legend->SetHeader("tW #rightarrow e+jets (1 b-tag)");
  else if(channel == "Muon" && n_btags == "1b") legend->SetHeader("tW #rightarrow #mu+jets (1-btag)");
  else if(channel == "Electron" && n_btags == "2b") legend->SetHeader("tW #rightarrow e+jets (#geq 2 b-tags)");
  else if(channel == "Muon" && n_btags == "2b") legend->SetHeader("tW #rightarrow #mu+jets (#geq 2 b-tags)");*/
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
  for(int i = 1; i < hist_eff->GetNbinsX()+1; i++) {
    bincontent = hist_rec->GetBinContent(i)/hist_gen->GetBinContent(i);
    hist_eff->SetBinContent(i,bincontent);
    uncertainty = bincontent*sqrt(pow(hist_rec->GetBinError(i)/hist_rec->GetBinContent(i),2)+pow(hist_gen->GetBinError(i)/hist_gen->GetBinContent(i),2));
    hist_eff->SetBinError(i,uncertainty);
    result.push_back(bincontent);
  }
  hist_eff->SetMaximum(1.3*hist_eff->GetMaximum());
  hist_eff->SetMinimum(0.3*hist_eff->GetMinimum());
  if(n_btags == "1b") {
    hist_eff->SetMaximum(0.15);
    hist_eff->SetMinimum(0.05);
    hist_eff->GetYaxis()->SetNdivisions(405);
  }
  else if(n_btags == "2b") {
    hist_eff->SetMaximum(0.075);
    hist_eff->SetMinimum(0.005);
    hist_eff->GetYaxis()->SetNdivisions(404);
  }
  //hist_eff->GetYaxis()->SetNdivisions(503);
  hist_eff->GetYaxis()->SetLabelSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetLabelSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetTitleOffset(1.6);
  hist_eff->GetXaxis()->SetNdivisions(507);
  hist_eff->SetLineStyle(1);

  hist_eff->SetLineColor(kBlack);
  hist_eff->GetXaxis()->SetTitle("top-quark or top-jet p_{T} [GeV]");
  hist_eff->GetYaxis()->SetTitle("Efficiency");
  hist_eff->GetYaxis()->CenterTitle();
  hist_eff->GetYaxis()->SetTitleSize(0.05/(1-splitter));
  hist_eff->GetXaxis()->SetTitleSize(0.05/(1-splitter));
  hist_eff->GetYaxis()->SetTitleOffset(0.9);

  hist_eff->Draw("e1");

  TGaxis *axis = new TGaxis(200,hist_eff->GetMaximum(),1200,hist_eff->GetMaximum(),200,1200,507,"-");
  axis->Draw();
 
  //==================//
  // SAVE AS EPS PLOT //
  //==================//

  c->cd();
  c->SaveAs("recoEfficiency_"+channel+"_"+n_btags+".eps");

  return result;
}


vector<double> getDataPoints(TString channel, TString n_btags) {

  TString input_path_measurement = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/theta-workdir/rootfiles/theta-output_stat.root";
  cout << "Make sure to run THETA to get your input histograms!" << endl;
  TFile* input_file_measurement = TFile::Open(input_path_measurement, "READ");

  vector<double> result;

  for(TString topptbin : {"200to300", "300to400", "400to600", "600to1200"}) {

    TString input_hist_measurement = "BDT_"+topptbin+"_"+n_btags+(channel == "Electron" ? "_ele__" : "_muo__")+"SingleTop_tWch";
    TH1F* hist_measured = (TH1F*)input_file_measurement->Get(input_hist_measurement);

    hist_measured->Rebin(hist_measured->GetSize()-2);
    result.push_back(hist_measured->GetBinContent(1));
    cout << hist_measured->GetBinContent(1) << endl;
    result.push_back(hist_measured->GetBinError(1));
    cout << hist_measured->GetBinError(1) << endl;
  }

  return result;
}


void makeFinalPlots(TString channel, TString n_btags, vector<double> datapoints) {

  TString input_path_gen = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop_tWch/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded tW sample (top+antitop)!" << std::endl;
  TFile* input_file_gen = TFile::Open(input_path_gen+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TString genhist;
  if(channel == "Electron") genhist = "nocuts/EleChannel_TopPt_GenCut";
  else if(channel == "Muon") genhist = "nocuts/MuoChannel_TopPt_GenCut";
  TH1F* hist_gen = (TH1F*)input_file_gen->Get(genhist);

  TCanvas* c = new TCanvas("c","",600,600);

  hist_gen->Draw("hist");
  hist_gen->SetTitle("");
  hist_gen->GetXaxis()->SetTitle("top-quark p_{T} [GeV]");
  hist_gen->GetYaxis()->SetTitle("d#sigma/dp_{T} [pb]");
  hist_gen->SetLineColor(597);
  hist_gen->SetLineStyle(2);
  hist_gen->SetLineWidth(2);
  //hist_gen->GetXaxis()->SetRange(3,5); // used to show less top-pt bins, useful to increase visibility of higher-pt bins

  TH1F* hist_data = (TH1F*)hist_gen->Clone();
  for(int i = 2; i < 6; i++) {
    hist_data->SetBinContent(i,datapoints.at((i-2)*2));
    hist_data->SetBinError(i,datapoints.at((i-2)*2+1));
  }
  
  hist_data->Draw("p same e1");
  hist_data->SetLineColor(kBlack);
  hist_data->SetLineStyle(1);
  hist_data->SetMarkerStyle(8);

  double lumi = 1/35867.;
  //hist_gen->Scale(lumi);
  //hist_data->Scale(lumi);

  //c->SetLogy();

  //==================//
  // SAVE AS EPS PLOT //
  //==================//

  c->cd();
  c->SaveAs("diffXsection_"+channel+"_"+n_btags+".eps");

}


void finalResults() {

  for(TString channel : {"Electron", "Muon"}) {
    for(TString n_btags : {"1b", "2b"}) {
      auto efficiencies = recoEfficiency(channel, n_btags);
      cout << "Results for the " << channel << " channel (" << n_btags << " category):" << endl;
      for(int i = 0; i < efficiencies.size(); i++) {
	cout << efficiencies.at(i) << endl;
      }
      auto datapoints = getDataPoints(channel, n_btags);
      for(int i = 0; i < datapoints.size(); i++) {
	datapoints.at(i) = datapoints.at(i)/efficiencies.at(i/2);
      }
      makeFinalPlots(channel, n_btags, datapoints);
    }
  }
}
