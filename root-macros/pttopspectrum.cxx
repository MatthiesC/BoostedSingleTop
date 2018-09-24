#include <iostream>
#include <vector>

void pttopspectrum() {
   
  TString input_path = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded ones." << std::endl;
  std::vector<TFile*> input_files;
  std::vector<TH1F*> hists, hists2;
  for (TString channel : {"tWch_incl", "tch", "sch_incl"}) {
    TFile* input_file = TFile::Open(input_path+"uhh2.AnalysisModuleRunner.MC.SingleTop_"+channel+".root", "READ");
    input_files.push_back(input_file);
    TH1F* hist = (TH1F*)input_file->Get("singletopgen/Pt_top");
    TH1F* hist2 = (TH1F*)input_file->Get("singletopgen/Pt_top");
    hists.push_back(hist);
    hists2.push_back(hist2);
  }

  TCanvas* c = new TCanvas("c","",600,600);

  gStyle->SetOptStat(kFALSE);
  c->SetLogy();
  
  THStack *hs = new THStack("hs","");
  THStack *hs2 = new THStack("hs2","");

  for(Int_t i = 0; i < 3; i++) {
    hists.at(i)->SetLineWidth(3);
    hists2.at(i)->SetLineWidth(3);
    hists.at(i)->Rebin(50);
    hists2.at(i)->Rebin(50);
    hs->Add(hists.at(i));
    hs2->Add(hists2.at(i));
  }

  //hists.at(0)->SetFillColor(597);
  hists.at(0)->SetLineColor(597);
  //hists.at(1)->SetFillColor(kRed+1);
  hists.at(1)->SetLineColor(kRed+1);
  //hists.at(2)->SetFillColor(kOrange+7);
  hists.at(2)->SetLineColor(kOrange+7);

  //hists.at(0)->SetFillStyle(3354);
  //hists.at(1)->SetFillStyle(3345);
  //hists.at(2)->SetFillStyle(3354);

  hs->Draw("nostack hist");
  //hs2->Draw("nostack hist same");

  //gStyle->SetPalette(kOcean);
  hs->GetYaxis()->SetTitle("Events");
  hs->GetYaxis()->SetTitleOffset(1.4);
  hs->GetXaxis()->SetTitle("top-quark p_{T}");
  hs->GetXaxis()->SetTitleOffset(1.4);
  hs->SetMinimum(0.4);

  c->cd();

  TText *simulation = new TText(.1,.91,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();

  auto legend = new TLegend(0.7,0.7,0.9,0.9);
  legend->AddEntry(hists.at(0), "tW");
  legend->AddEntry(hists.at(1), "t-ch.");
  legend->AddEntry(hists.at(2), "s-ch.");
  legend->Draw();

  c->SaveAs("pttopspectrum.eps");

}
