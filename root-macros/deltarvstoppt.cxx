#include <iostream>
#include <vector>
#include <list>

void deltarvstoppt() {
   
  TString input_path = "/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Generator/SingleTop/";
  std::cout << "Make sure to already have run SFramePlotter over the input files to have the hadded ones." << std::endl;
  TFile* input_file = TFile::Open(input_path+"uhh2.AnalysisModuleRunner.MC.SingleTop_tWch_incl.root", "READ");
  TH2F* hist = (TH2F*)input_file->Get("singletopgen/Pt_top_vs_deltaR_top");
  hist->SetTitle("");
  TCanvas* c = new TCanvas("c","",800,600);

  TPad* p = new TPad("p","pad",.0,0,1,1);
  Double_t margin = 0.1;
  p->SetMargin(0.2,0.2,0.1,0.1);
  p->Draw();
  p->cd();
  //gPad->SetFrameFillColor(0);

  gStyle->SetOptStat(kFALSE);

  p->SetLogz();

  //gStyle->SetPalette(kCherry);






   const Int_t Number = 9;
   //Double_t Red[Number]    = { 37./255., 102./255., 157./255., 188./255., 196./255., 214./255., 223./255., 235./255., 251./255.};
   Double_t Red[Number] = {251./255., 235./255., 223./255., 214./255., 196./255., 188./255., 157./255., 102./255., 37./255.};
   //Double_t Green[Number]  = { 37./255.,  29./255.,  25./255.,  37./255.,  67./255.,  91./255., 132./255., 185./255., 251./255.};
   Double_t Green[Number]  = { 251./255., 185./255., 132./255., 91./255., 67./255., 37./255., 25./255., 29./255., 37./255.};
   //Double_t Blue[Number]   = { 37./255.,  32./255.,  33./255.,  45./255.,  66./255.,  98./255., 137./255., 187./255., 251./255.};
   Double_t Blue[Number]   = { 251./255., 187./255., 137./255., 98./255., 66./255., 45./255., 33./255., 32./255., 37./255. };
   Double_t Length[Number] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};

   Int_t nb=50;
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);



   //gStyle->SetPalette(kAquamarine);
  

  //TColor::InvertPalette();
  hist->GetYaxis()->SetTitle("max. #DeltaR_{ij} of top-quark decay products");
  hist->GetYaxis()->SetTitleOffset(1.4);
  hist->GetXaxis()->SetTitle("top-quark p_{T} [GeV]");
  hist->GetXaxis()->SetTitleOffset(1.4);

  hist->Rebin2D(2,2);

  TH2F* hist2 = (TH2F*)hist->Clone("hist2");

  Double_t rebinX[] = {0,100,200,300,400,600,800,1200,1600,2000};
  vector< vector<Int_t> > rebinX_ = {{1},{2},{3},{4},{5,6},{7,8},{9,10,11,12},{13,14,15,16},{17,18,19,20}}; // these are the bin numbers

  /*vector<Double_t> BinContentsX;
  for(Int_t i = 0; i < rebinX_.size(); i++)
    {
      Double_t BinContentX = 0;
      for(Int_t xbin : rebinX_.at(i))
	{
	  for(Int_t ybin = 0; ybin < 22; ybin++)
	    {
	      BinContentX += hist->GetBinContent(xbin,ybin);
	    }
	}
      BinContentsX.push_back(BinContentX);
      }*/

  /*for(int i = 0; i < 52; i++)
    {
      double scale = 0;
      for(int j = 0; j < 52; j++)
	{
	  scale += hist->GetBinContent(i,j);
	}
      for(int j = 0; j < 52; j++)
	{
	  hist2->SetBinContent(i,j,(hist->GetBinContent(i,j))/scale);
	}
	}*/

  hist2->Rebin2D(1,2);
  hist2->GetXaxis()->SetRangeUser(0,1000);
  hist2->Draw("colz");

  /*for(int i = 0; i < hist2->GetXaxis()->GetNbins()+1; i++)
    {
      double scale = 0;
      for(int j = 0; j < hist2->GetYaxis()->GetNbins()+1; j++)
	{
	  scale += hist->GetBinContent(i,j);
	}
      for(int j = 0; j < hist2->GetYaxis()->GetNbins()+1; j++)
	{
	  hist2->SetBinContent(i,j,(hist2->GetBinContent(i,j))/scale);
	}
    }
  */
  p->cd();

  TText *simulation = new TText(.2,.91,"Simulation");
  simulation->SetNDC();
  //simulation->SetTextAlign(12);
  simulation->SetTextFont(73);
  simulation->SetTextSize(24);
  simulation->Draw();

  TLatex *lumi = new TLatex(.8,.91,"tW-channel (incl. decays)");
  lumi->SetNDC();
  lumi->SetTextAlign(31);
  lumi->SetTextFont(43);
  lumi->SetTextSize(24);
  lumi->Draw();

  c->SaveAs("deltarvstoppt.eps");

}
