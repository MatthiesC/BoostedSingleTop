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

  gStyle->SetEndErrorSize(0);

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
  hist_eff->GetXaxis()->SetTitleOffset(1.5);
  hist_eff->GetXaxis()->SetNdivisions(507);
  hist_eff->SetLineStyle(1);

  hist_eff->SetLineColor(kBlack);
  hist_eff->GetXaxis()->SetTitle("top-quark or top-jet p_{T} [GeV]");
  hist_eff->GetYaxis()->SetTitle("Correction factor");
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
  c->SaveAs("plots/recoEfficiency_"+channel+"_"+n_btags+".eps");

  return result;
}


// calculates the total uncertainty of a specific bin
double calcUncertainty(double mean, double mean_data, vector<double> variations, bool bStatOnly) {

  double uncertainty = pow(0.025*mean,2); // lumi unc.
  uncertainty += mean_data; // stat. unc. of real data
  for (auto var : variations) {
    uncertainty += pow(mean-var,2);
  }
  uncertainty = sqrt(uncertainty);
  if(bStatOnly) uncertainty = sqrt(mean_data);

  return uncertainty;
}


vector<double> getDataPoints(TString channel, TString n_btags, bool bStatOnly) {

  //TString input_path_measurement = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/theta-workdir/rootfiles/theta-output"++".root";
  //cout << "Make sure to run THETA to get your input histograms!" << endl;
  //TFile* input_file_measurement = TFile::Open(input_path_measurement, "READ");

  vector<double> result;

  for(TString topptbin : {"200to300", "300to400", "400to600", "600to1200"}) {

    TString input_path = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/theta-workdir/rootfiles/theta-output_"+channel+"_"+n_btags+"_"+topptbin+".root";
    TFile* input_file = TFile::Open(input_path, "READ");

    TString input_hist_mean = "BDT_"+topptbin+"_"+n_btags+(channel == "Electron" ? "_ele__" : "_muo__")+"SingleTop_tWch";
    TH1F* hist_mean = (TH1F*)input_file->Get(input_hist_mean);

    TString input_hist_data = "BDT_"+topptbin+"_"+n_btags+(channel == "Electron" ? "_ele__" : "_muo__")+"DATA";
    TH1F* hist_data = (TH1F*)input_file->Get(input_hist_data);

    hist_mean->Rebin(hist_mean->GetSize()-2);
    double mean = hist_mean->GetBinContent(1);
    result.push_back(mean);
    //result.push_back(hist_mean->GetBinError(1)); // MC stat. unc. --- Not what we want!!!

    hist_data->Rebin(hist_data->GetSize()-2);
    double mean_data = hist_data->GetBinContent(1);

    for (auto ud : {"plus", "minus"}) {
      vector<double> variations;
      for (auto unc : {"SingleTop_tWch_rate", "btagSFudsg", "btagSFbc"}) {
	TString input_hist_var = input_hist_mean+"__"+unc+"__"+ud;
	TH1F* hist_var = (TH1F*)input_file->Get(input_hist_var);
	hist_var->Rebin(hist_var->GetSize()-2);
	double var = hist_var->GetBinContent(1);
	variations.push_back(var);
      }
      result.push_back(calcUncertainty(mean, mean_data, variations, bStatOnly)); // syst. unc. if bStatOnly = false, else: tot. unc.
      result.push_back(calcUncertainty(mean, mean_data, variations, true)); // stat. unc.
    }
  }

  return result;
}


void printToScreen(float y[], float eyh[], float eyl[], float eyh_stat[], float eyl_stat[]) {

  cout << "All uncertainties are given relatively:" << endl;
  float syst_up, syst_down;
  for (int i = 0; i < 4; i++) {
    syst_up   = sqrt(pow(eyh[i],2)-pow(eyh_stat[i],2));
    syst_down = sqrt(pow(eyl[i],2)-pow(eyl_stat[i],2));
    cout << y[i] << endl;
    cout << " + " << eyh_stat[i]/y[i] << " - " << eyl_stat[i]/y[i] << " (stat.)" << endl;
    cout << " + " << syst_up/y[i] << " - " << syst_down/y[i] << " (syst.)" << endl;
    cout << " + " << eyh[i]/y[i] << " - " << eyl[i]/y[i] << " (tot.)" << endl;
  }
}


void makeFinalPlots(TString channel, TString n_btags, vector<double> datapoints, const int number_of_ptbins) {

  TString input_path_gen = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop_tWch/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded tW sample (top+antitop)!" << std::endl;
  TFile* input_file_gen = TFile::Open(input_path_gen+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TString genhist;
  if(channel == "Electron") genhist = "nocuts/EleChannel_TopPt_GenCut";
  else if(channel == "Muon") genhist = "nocuts/MuoChannel_TopPt_GenCut";
  TH1F* hist_gen = (TH1F*)input_file_gen->Get(genhist);

  TCanvas* c2 = new TCanvas("c2","",600,800);

  double splitter = 0.375; // decide where to put the border between the two pads
  
  //==============//
  // UPPER HISTOS //
  //==============//

  TPad* upper2 = new TPad("upper2","",0,splitter,1,1);

  upper2->SetMargin(.15,.05,.0,.1*6/5);
  upper2->Draw();
  upper2->cd();

  hist_gen->Draw("hist");
  hist_gen->SetTitle("");
  hist_gen->GetXaxis()->SetTitle("top-quark p_{T} [GeV]");
  hist_gen->GetYaxis()->SetTitle("d#sigma/dp_{T} #scale[0.7]{#left[#frac{pb}{100 GeV}#right]}");
  hist_gen->SetLineColor(597);
  hist_gen->SetLineStyle(2);
  hist_gen->SetLineWidth(2);
  hist_gen->GetXaxis()->SetRange(2,5); // used to show less top-pt bins, useful to increase visibility of higher-pt bins

  hist_gen->GetYaxis()->SetTitleOffset(1.3); // 1.5
  hist_gen->GetYaxis()->SetTitleSize(0.05);
  hist_gen->GetYaxis()->SetLabelSize(0.05);
  hist_gen->GetXaxis()->SetTitleOffset(1.4);
  hist_gen->GetXaxis()->SetNdivisions(507);

  // Divide by 100 GeV:
  hist_gen->SetBinContent(4,hist_gen->GetBinContent(4)/2.);
  hist_gen->SetBinError(4,hist_gen->GetBinError(4)/2.);
  hist_gen->SetBinContent(5,hist_gen->GetBinContent(5)/6.);
  hist_gen->SetBinError(5,hist_gen->GetBinError(5)/6.);

  TH1F* hist_data = (TH1F*)hist_gen->Clone();

  //https://root.cern.ch/root/roottalk/roottalk02/0249.html
  float x[number_of_ptbins];
  float y[number_of_ptbins];
  Float_t eyl[] = {0,0,0,0};
  Float_t eyh[] = {0,0,0,0};
  float exl[] = {50,50,100,300};
  float exh[] = {50,50,100,300};
  float eyl_stat[] = {0,0,0,0};
  float eyh_stat[] = {0,0,0,0};

  double luminosity = 35867.;

  for(int i = 2; i < 6; i++) {
    hist_data->SetBinContent(i,datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)));
    hist_data->SetBinError(i,datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)+1)); // only upper variations is taken as error for now
    y[i-2] = hist_data->GetBinContent(i)/luminosity*(1./(i == 4 ? 2 : (i == 5 ? 6 : 1)));
    x[i-2] = hist_data->GetXaxis()->GetBinCenter(i);
    eyh[i-2] = (datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)+1))/luminosity*(1./(i == 4 ? 2 : (i == 5 ? 6 : 1)));
    eyl[i-2] = (datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)+3))/luminosity*(1./(i == 4 ? 2 : (i == 5 ? 6 : 1)));
    eyh_stat[i-2] = (datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)+2))/luminosity*(1./(i == 4 ? 2 : (i == 5 ? 6 : 1)));
    eyl_stat[i-2] = (datapoints.at((i-2)*(datapoints.size()/number_of_ptbins)+4))/luminosity*(1./(i == 4 ? 2 : (i == 5 ? 6 : 1)));

    //printToScreen(y[i-2],eyh[i-2],eyl[i-2],eyh_stat[i-2],eyl_stat[i-2]);
  }

  printToScreen(y,eyh,eyl,eyh_stat,eyl_stat);

  graph = new TGraphAsymmErrors(number_of_ptbins,x,y,exl,exh,eyl,eyh);
  graph_stat = new TGraphAsymmErrors(number_of_ptbins,x,y,exl,exh,eyl_stat,eyh_stat);

  //hist_data->Draw("a same");
  hist_data->SetLineColor(kBlack);
  hist_data->SetLineStyle(1);
  hist_data->SetMarkerStyle(8);

  //double luminosity = 1/35867.;
  hist_gen->Scale(1./luminosity);
  //hist_data->Scale(luminosity);

  graph->SetMarkerStyle(8); // 21
  graph->SetLineWidth(2);
  graph->Draw("p");

  graph_stat->SetMarkerStyle(8); // 21
  graph_stat->SetLineWidth(2);
  graph_stat->Draw("p");

  upper2->SetLogy();

  gStyle->SetEndErrorSize(8);

  hist_gen->GetXaxis()->SetNdivisions(406);
  //hist_gen->SetMaximum(5.*hist_gen->GetMaximum());
  //hist_gen->SetMinimum(.1*hist_data->GetBinContent(5));
  hist_gen->GetYaxis()->SetRangeUser(0.0002,1.);

  /*TText *simulation = new TText(.15,1-0.09*6/5,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();*/

  TString category_1 = "tW #rightarrow ";
  TString category_2 = (channel == "Electron" ? "e+jets " : "#mu+jets ");
  TString category_3 = (n_btags == "1b" ? "(1 b-tag)" : "(#geq 2 b-tags)");

  TLatex *category = new TLatex(.15,1-0.09*6/5,category_1+category_2+category_3); //.2,1-0.09*6/5-0.1
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
  legend->AddEntry(hist_data, "Data");
  legend->AddEntry(hist_gen, "POWHEG");
  legend->Draw();
  legend->SetBorderSize(0);



  //=============//
  // RATIO HISTO //
  //=============//

  c2->cd();
  TPad* lower2 = new TPad("lower2", "",0,0,1,splitter);

  lower2->SetMargin(.15,.05,.3,.0);
  lower2->Draw();
  lower2->cd();

  TH1F* hist_ratio = (TH1F*)hist_gen->Clone();

  for(int i = 2; i < 6; i++) {
    float scalefactor = hist_ratio->GetBinContent(i);
    y[i-2] = y[i-2]/scalefactor;
    eyh[i-2] = eyh[i-2]/scalefactor;
    eyl[i-2] = eyl[i-2]/scalefactor;
    eyh_stat[i-2] = eyh_stat[i-2]/scalefactor;
    eyl_stat[i-2] = eyl_stat[i-2]/scalefactor;
    hist_ratio->SetBinContent(i,1); // to get a nice blue line in the center at 1
  }

  hist_ratio->Draw();

  hist_ratio->SetMaximum(2.7);

  hist_ratio->GetYaxis()->SetLabelSize(0.05/(1-splitter));
  hist_ratio->GetXaxis()->SetLabelSize(0.05/(1-splitter));
  hist_ratio->GetXaxis()->SetTitleOffset(1.5);
  hist_ratio->GetXaxis()->SetNdivisions(507);
  //hist_ratio->SetLineStyle(1);

  //hist_ratio->SetLineColor(kBlack);
  hist_ratio->GetXaxis()->SetTitle("top-quark p_{T} [GeV]");
  hist_ratio->GetYaxis()->SetTitle("Data / MC");
  hist_ratio->GetYaxis()->CenterTitle();
  hist_ratio->GetYaxis()->SetTitleSize(0.05/(1-splitter));
  hist_ratio->GetXaxis()->SetTitleSize(0.05/(1-splitter));
  hist_ratio->GetYaxis()->SetTitleOffset(0.9);

  ratio_graph = new TGraphAsymmErrors(number_of_ptbins,x,y,exl,exh,eyl,eyh);
  ratio_graph_stat = new TGraphAsymmErrors(number_of_ptbins,x,y,exl,exh,eyl_stat,eyh_stat);

  ratio_graph->SetMarkerStyle(8);
  ratio_graph->SetLineWidth(2);
  ratio_graph->Draw("p");

  ratio_graph_stat->SetMarkerStyle(8);
  ratio_graph_stat->SetLineWidth(2);
  ratio_graph_stat->Draw("p");

  TGaxis *axis = new TGaxis(200,hist_ratio->GetMaximum(),1200,hist_ratio->GetMaximum(),200,1200,507,"-");
  axis->Draw();

  //==================//
  // SAVE AS EPS PLOT //
  //==================//

  c2->cd();
  c2->SaveAs("plots/diffXsection_"+channel+"_"+n_btags+".eps");

}


void finalResults() {

  const int number_of_ptbins = 4;
  const bool only_stat_errors = false;

  for(TString channel : {"Electron", "Muon"}) {
    for(TString n_btags : {"1b", "2b"}) {
      auto efficiencies = recoEfficiency(channel, n_btags);
      cout << "Results for the " << channel << " channel (" << n_btags << " category):" << endl;
      for(int i = 0; i < efficiencies.size(); i++) {
	cout << efficiencies.at(i) << endl;
      }
      auto datapoints = getDataPoints(channel, n_btags, only_stat_errors);
      for(int i = 0; i < datapoints.size(); i++) {
	datapoints.at(i) = datapoints.at(i)/efficiencies.at(i/(datapoints.size()/number_of_ptbins));
      }
      makeFinalPlots(channel, n_btags, datapoints, number_of_ptbins);
    }
  }
}
