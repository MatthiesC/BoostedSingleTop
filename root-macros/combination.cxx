vector< vector<double> > calcCorrFactors(TString channel, TString nbtags) {

  cout << "Calculate correction factors for: " << channel << " " << nbtags << endl;

  // Get the gen and reco hists:

  TString input_path_rec = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/"+channel+"/";
  TString input_path_gen = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop_tWch/";
  TFile* input_file_rec = TFile::Open(input_path_rec+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TFile* input_file_gen = TFile::Open(input_path_gen+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch.root", "READ");
  TString input_hist_rec = "MVA_1t"+nbtags+"/nocuts_top_pt_rebin";
  TH1F* hist_rec = (TH1F*)input_file_rec->Get(input_hist_rec);
  TString input_hist_gen;
  if (channel == "Electron") input_hist_gen = "nocuts/EleChannel_TopPt_GenCut";
  else if (channel == "Muon") input_hist_gen = "nocuts/MuoChannel_TopPt_GenCut";
  TH1F* hist_gen = (TH1F*)input_file_gen->Get(input_hist_gen);

  // Calculate the correction factors:

  vector< vector<double> > result;
  vector<double> result_meanv, result_error;

  double rec, gen;
  double rec_err, gen_err;

  for (int i = 1; i <= 4; i++) {
    rec = hist_rec->GetBinContent(i);
    rec_err = hist_rec->GetBinError(i);
    gen = hist_gen->GetBinContent(i+1);
    gen_err = hist_gen->GetBinError(i+1);
    cout << rec/gen << endl;
    result_meanv.push_back(rec/gen);
    result_error.push_back(rec/gen*sqrt(rec_err*rec_err/(rec*rec)+gen_err*gen_err/(gen*gen)));
  }

  result.push_back(result_meanv); // corrFactors.at(0).at(i) --> the corr factors!
  result.push_back(result_error); // corrFactors.at(1).at(i) --> the bin error from MC stat uncertainty!

  // Close input files and return result:

  input_file_rec->Close();
  input_file_gen->Close();

  return result;
}


void plotCorrFactors(vector< vector< vector<double> > > vecCorrFactors) {

  double ptbins[] = {200,300,400,600,1200};
  TH1F* h_e1 = new TH1F("h_e1","e+jets (1b)",4,ptbins);
  TH1F* h_e2 = new TH1F("h_e2","e+jets (2b+)",4,ptbins);
  TH1F* h_m1 = new TH1F("h_m1","#mu+jets (1b)",4,ptbins);
  TH1F* h_m2 = new TH1F("h_m2","#mu+jets (2b+)",4,ptbins);

  vector<TH1F*> hists = {h_e1, h_e2, h_m1, h_m2};

  for (int j = 0; j < 4; j++) {
    for (int i = 1; i <=4; i++) {
      hists.at(j)->SetBinContent(i,vecCorrFactors.at(j).at(0).at(i-1));
      hists.at(j)->SetBinError(i,vecCorrFactors.at(j).at(1).at(i-1));
    }
  }

  double textsize = 0.05;

  TCanvas* c = new TCanvas("c","",600,600);
  c->cd();
  TPad* p = new TPad("p","",0,0,1,1);
  p->SetMargin(.15,.05,.15,.05);
  p->Draw();
  p->cd();

  TLegend* l = new TLegend(0.55,0.65,0.9,0.9);

  int x = 1;
  for (auto h : hists) {

    h->Draw("hist e1 same");

    h->GetXaxis()->SetTitle("top quark p_{T} [GeV]");
    h->GetYaxis()->SetTitle("Correction factor");
    h->GetXaxis()->SetLabelSize(textsize);
    h->GetYaxis()->SetLabelSize(textsize);
    h->GetXaxis()->SetTitleSize(textsize);
    h->GetYaxis()->SetTitleSize(textsize);
    h->GetXaxis()->SetTitleOffset(1.4);
    h->GetYaxis()->SetTitleOffset(1.6);
  
    h->GetXaxis()->SetNdivisions(406);

    h->SetMinimum(0.0001);
    h->SetMaximum(0.25);

    h->SetLineStyle(x);
    h->SetLineColor(x);
    h->SetLineWidth(3);

    l->AddEntry(h,h->GetTitle());

    x++;
  }

  l->SetBorderSize(0);
  l->Draw();
  h_e1->SetTitle("");

  gStyle->SetOptStat(0);

  c->cd();
  c->SaveAs("plots/corrFactors.eps");

}


vector< vector<double> > readDataPoints(TString channel, TString nbtags, vector<TString> uncertainties) {

  vector< vector<double> > resultAllFourBins;

  for(TString topptbin : {"200to300", "300to400", "400to600", "600to1200"}) {

    vector<double> result;
    double mean = 0;

    TString input_path = "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/theta-workdir/rootfiles/theta-output_"+channel+"_"+nbtags+"_"+topptbin+".root";
    TFile* input_file = TFile::Open(input_path, "READ");

    // get mean value:
    TString input_hist_mean = "BDT_"+topptbin+"_"+nbtags+(channel == "Electron" ? "_ele__" : "_muo__")+"SingleTop_tWch";
    TH1F* hist_mean = (TH1F*)input_file->Get(input_hist_mean);
    mean = hist_mean->Integral(0,hist_mean->GetSize()-1);
    result.push_back(mean);

    // get up/down values for each variation:
    for (auto unc : uncertainties) {
      vector<double> values;
      for (auto pm : {"plus", "minus"}) {
	TString input_hist_var = input_hist_mean+"__"+unc+"__"+pm;
	TH1F* hist_var = (TH1F*)input_file->Get(input_hist_var);
	values.push_back(hist_var->Integral(0,hist_var->GetSize()-1));
      }
      if(values.at(0) < values.at(1)) swap(values.at(0), values.at(1)); // We want the upper value first, the lower value second
      result.push_back(abs(mean-values.at(0)));
      //cout << abs(mean-values.at(0)) << endl;
      result.push_back(abs(mean-values.at(1)));
      //cout << abs(mean-values.at(1)) << endl;
    }

    resultAllFourBins.push_back(result);
  }

  return resultAllFourBins;
}


vector< vector<double> > calcCrossSections(vector< vector<double> > dataPoints, vector< vector<double> > corrFactors, vector<TString> uncertaintiesText) {
  
  double luminosity = 35867;

  auto result = dataPoints;

  cout << "Xsection values: " << endl;
  for (int i = 0; i < dataPoints.size(); i++) {
    for (int j = 0; j < dataPoints.at(i).size(); j++) {
      result.at(i).at(j) = result.at(i).at(j)/corrFactors.at(0).at(i);
      result.at(i).at(j) = result.at(i).at(j)/luminosity;
      if (i == 2) result.at(i).at(j) = result.at(i).at(j)/2.; // 400to600 bin
      if (i == 3) result.at(i).at(j) = result.at(i).at(j)/6.; // 600to1200 bin
    }
    double mean = result.at(i).at(0);
    cout << "Mean: " << mean << endl;
    cout << "Uncertainties: " << endl;
    for (int j = 1; j < dataPoints.at(i).size(); j=j+2) {
      cout << uncertaintiesText.at(j/2);
      cout << "+ " << result.at(i).at(j);
      cout << " ( + " << result.at(i).at(j)/mean*100 << " )  ";
      cout << " - " << result.at(i).at(j+1);
      cout << " ( - " << result.at(i).at(j+1)/mean*100 << " )  " << endl;;
    }
  }

  return result;
}


void combination() {

  cout << "Everything is happy." << endl;

  vector<TString> uncertainties = {"SingleTop_tWch_rate", "btagSFudsg", "btagSFbc", "jec", "jer", "pileup"};
  vector<TString> uncertaintiesText = {"tW rate:   ",
				       "b-mistag:  ",
				       "b-tag SF:  ",
				       "JEC:       ",
				       "JER:       ",
				       "Pile-up:   "};

  vector< vector< vector<double> > > vecCorrFactors;
  for (TString channel : {"Electron", "Muon"}) {
    for (TString nbtags : {"1b", "2b"}) {
      auto corrFactors = calcCorrFactors(channel, nbtags);
      vecCorrFactors.push_back(corrFactors);
      auto dataPoints = readDataPoints(channel, nbtags, uncertainties);
      auto crossSections = calcCrossSections(dataPoints, corrFactors, uncertaintiesText);
    }
  }
  plotCorrFactors(vecCorrFactors);

}
