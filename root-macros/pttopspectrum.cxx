#include <iostream>
#include <vector>

void pttopspectrum() {
   
  bool bNormalized = false;

  TString input_path = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded ones." << std::endl;
  std::vector<TFile*> input_files;
  std::vector<TH1F*> hists, hists2;
  for (TString channel : {"tWch_incl", "tch", "sch_incl"}) {
    TFile* input_file = TFile::Open(input_path+"uhh2.AnalysisModuleRunner.MC.SingleTop_"+channel+".root", "READ");
    input_files.push_back(input_file);
    TH1F* hist = (TH1F*)input_file->Get("singletopgen/Pt_top");
    TH1F* hist2 = (TH1F*)hist->Clone("hist2");
    //TF1 *f1 = new TF1("f1","1.2",0,2000);
    //hist2->Multiply(f1,1);
    //Double_t shift = 0.5*hist->GetBinWidth(1);
    //hist2->GetXaxis()->SetLimits(0+shift,2000+shift);
    hists.push_back(hist);
    hists2.push_back(hist2);
  }

  TCanvas* c = new TCanvas("c","",600,600);

  gStyle->SetOptStat(kFALSE);
  c->SetLogy();
  
  THStack *hs = new THStack("hs","");
  THStack *hs2 = new THStack("hs2","");

  for(Int_t i = 2; i > -1; i--) {
    hists.at(i)->SetLineWidth(2);
    //hists2.at(i)->SetLineWidth(3);
    hists.at(i)->Rebin(50);
    if(bNormalized){
      Double_t scale = 1./(hists.at(i)->Integral());
      std::cout << scale << std::endl;
      hists.at(i)->Scale(scale);
    }
    //hists2.at(i)->Rebin(50);
    hs->Add(hists.at(i));
    //hs2->Add(hists2.at(i));
  }

  hists.at(0)->SetLineColor(597);
  //hists.at(0)->SetLineStyle(1);
  hists.at(1)->SetLineColor(kOrange+7);
  hists.at(1)->SetLineStyle(2);
  hists.at(2)->SetLineColor(kGreen+2);
  hists.at(2)->SetLineStyle(5);

  //hists.at(0)->SetFillStyle(3354);
  //hists.at(1)->SetFillStyle(3345);
  //hists.at(2)->SetFillStyle(3354);

  hs->Draw("nostack hist e1");
  //hs2->Draw("nostack hist same");

  //gStyle->SetPalette(kOcean);
  hs->GetYaxis()->SetTitle((!bNormalized ? "Events" : "#DeltaN/N"));
  hs->GetYaxis()->SetTitleOffset(1.4);
  hs->GetXaxis()->SetTitle("top-quark p_{T} [GeV]");
  hs->GetXaxis()->SetTitleOffset(1.4);
  hs->SetMinimum(0.4);
  if(bNormalized) hs->SetMinimum(0.4e-7);


  c->cd();

  TText *simulation = new TText(.1,.91,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();

  TLatex *lumi = new TLatex(.9,.91,"35.9 fb^{-1} (13 TeV)");
  lumi->SetNDC();
  lumi->SetTextAlign(31);
  lumi->SetTextFont(43);
  lumi->SetTextSize(24);
  lumi->Draw();

  auto legend = new TLegend(0.65,0.65,0.85,0.85);
  legend->AddEntry(hists.at(0), "tW");
  legend->AddEntry(hists.at(1), "t-ch.");
  legend->AddEntry(hists.at(2), "s-ch.");
  legend->Draw();
  legend->SetBorderSize(0);

  if(bNormalized) c->SaveAs("plots/pttopspectrum_norm.eps");
  else c->SaveAs("plots/pttopspectrum.eps");

}
