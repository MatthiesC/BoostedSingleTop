#include "UHH2/BoostedSingleTop/include/MVAHists.h"
#include "UHH2/core/include/Event.h"

//#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"

using namespace std;
using namespace uhh2;

MVAHists::MVAHists(Context & ctx, const string & dirname):
  Hists(ctx,dirname)
{
  h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");

  // Book all histograms here ...

  // BDT response
  double rebin_mvaD[] = {-1,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,1};
  unsigned int rebin_mvaD_bincount = sizeof(rebin_mvaD)/sizeof(*rebin_mvaD)-1;
  H_mvaD = book<TH1F>("mvaD", "BDT response", 40, -1, 1);
  H_mvaD_20bins = book<TH1F>("mvaD_20bins", "BDT response", 20, -1, 1);
  H_mvaD_rebin = book<TH1F>("mvaD_rebin", "BDT response", rebin_mvaD_bincount, rebin_mvaD);
  H_mvaD_1b = book<TH1F>("mvaD_1btags", "BDT response (1t1b)", 40, -1, 1);
  H_mvaD_2b = book<TH1F>("mvaD_2btags", "BDT response (1t2b)", 40, -1, 1);

  // Top-jet distributions after BDT cut
  double rebin_pt[5] = {200,300,400,600,1200};
  unsigned int rebin_pt_bincount = sizeof(rebin_pt)/sizeof(*rebin_pt)-1;
  double rebin_eta[11] = {-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5};
  unsigned int rebin_eta_bincount = sizeof(rebin_eta)/sizeof(*rebin_eta)-1;
  double rebin_eta_lept[9] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
  unsigned int rebin_eta_lept_bincount = sizeof(rebin_eta_lept)/sizeof(*rebin_eta_lept)-1;

  H_nocuts_top_pt_1000 = book<TH1F>("nocuts_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_nocuts_top_pt_rebin = book<TH1F>("nocuts_top_pt_rebin", "top-jet p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_nocuts_top_eta_1000 = book<TH1F>("nocuts_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_nocuts_top_eta_rebin = book<TH1F>("nocuts_top_eta_rebin", "top-jet #eta", rebin_eta_bincount, rebin_eta);
  H_nocuts_top_pt_vs_eta_1000 = book<TH2F>("nocuts_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_nocuts_top_pt_vs_eta_rebin = book<TH2F>("nocuts_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", rebin_eta_bincount, rebin_eta, rebin_pt_bincount, rebin_pt);

  H_nocuts_lept_pt = book<TH1F>("nocuts_lept_pt", "lepton p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_nocuts_lept_eta = book<TH1F>("nocuts_lept_eta", "lepton #eta", rebin_eta_lept_bincount, rebin_eta_lept);

  H_BDTcut_top_pt_1000 = book<TH1F>("BDTcut_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_BDTcut_top_pt_rebin = book<TH1F>("BDTcut_top_pt_rebin", "top-jet p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_BDTcut_top_eta_1000 = book<TH1F>("BDTcut_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_BDTcut_top_eta_rebin = book<TH1F>("BDTcut_top_eta_rebin", "top-jet #eta", rebin_eta_bincount, rebin_eta);
  H_BDTcut_top_pt_vs_eta_1000 = book<TH2F>("BDTcut_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_BDTcut_top_pt_vs_eta_rebin = book<TH2F>("BDTcut_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", rebin_eta_bincount, rebin_eta, rebin_pt_bincount, rebin_pt);

  H_BDTcut_lept_pt = book<TH1F>("BDTcut_lept_pt", "lepton p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_BDTcut_lept_eta = book<TH1F>("BDTcut_lept_eta", "lepton #eta", rebin_eta_lept_bincount, rebin_eta_lept);

  // now go through top_pt bins and write out BDT response for each
  H_200to300_BDT = book<TH1F>("200to300_mvaD", "BDT response (top-jet p_{T} #in (200,300] GeV)", 40, -1, 1);
  H_200to300_BDT_rebin = book<TH1F>("200to300_mvaD_rebin", "BDT response (top-jet p_{T} #in (200,300] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_200to300_top_eta = book<TH1F>("200to300_top_eta", "top-jet #eta (top-jet p_{T} #in (200,300] GeV)", rebin_eta_bincount, rebin_eta);
  H_200to300_lept_pt = book<TH1F>("200to300_lept_pt", "lepton p_{T} (top-jet p_{T} #in (200,300] GeV)", rebin_pt_bincount, rebin_pt);
  H_200to300_lept_eta = book<TH1F>("200to300_lept_eta", "lepton #eta (top-jet p_{T} #in (200,300] GeV)", rebin_eta_lept_bincount, rebin_eta_lept);

  H_300to400_BDT = book<TH1F>("300to400_mvaD", "BDT response (top-jet p_{T} #in (300,400] GeV)", 40, -1, 1);
  H_300to400_BDT_rebin = book<TH1F>("300to400_mvaD_rebin", "BDT response (top-jet p_{T} #in (300,400] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_300to400_top_eta = book<TH1F>("300to400_top_eta", "top-jet #eta (top-jet p_{T} #in (300,400] GeV)", rebin_eta_bincount, rebin_eta);
  H_300to400_lept_pt = book<TH1F>("300to400_lept_pt", "lepton p_{T} (top-jet p_{T} #in (300,400] GeV)", rebin_pt_bincount, rebin_pt);
  H_300to400_lept_eta = book<TH1F>("300to400_lept_eta", "lepton #eta (top-jet p_{T} #in (300,400] GeV)", rebin_eta_lept_bincount, rebin_eta_lept);

  H_400to600_BDT = book<TH1F>("400to600_mvaD", "BDT response (top-jet p_{T} #in (400,600] GeV)", 40, -1, 1);
  H_400to600_BDT_rebin = book<TH1F>("400to600_mvaD_rebin", "BDT response (top-jet p_{T} #in (400,600] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_400to600_top_eta = book<TH1F>("400to600_top_eta", "top-jet #eta (top-jet p_{T} #in (400,600] GeV)", rebin_eta_bincount, rebin_eta);
  H_400to600_lept_pt = book<TH1F>("400to600_lept_pt", "lepton p_{T} (top-jet p_{T} #in (400,600] GeV)", rebin_pt_bincount, rebin_pt);
  H_400to600_lept_eta = book<TH1F>("400to600_lept_eta", "lepton #eta (top-jet p_{T} #in (400,600] GeV)", rebin_eta_lept_bincount, rebin_eta_lept);

  H_600to1200_BDT = book<TH1F>("600to1200_mvaD", "BDT response (top-jet p_{T} #in (600,1200] GeV)", 40, -1, 1);
  H_600to1200_BDT_rebin = book<TH1F>("600to1200_mvaD_rebin", "BDT response (top-jet p_{T} #in (600,1200] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_600to1200_top_eta = book<TH1F>("600to1200_top_eta", "top-jet #eta (top-jet p_{T} #in (600,1200] GeV)", rebin_eta_bincount, rebin_eta);
  H_600to1200_lept_pt = book<TH1F>("600to1200_lept_pt", "lepton p_{T} (top-jet p_{T} #in (600,1200] GeV)", rebin_pt_bincount, rebin_pt);
  H_600to1200_lept_eta = book<TH1F>("600to1200_lept_eta", "lepton #eta (top-jet p_{T} #in (600,1200] GeV)", rebin_eta_lept_bincount, rebin_eta_lept);

}

void MVAHists::fill(const Event & e)
{
  double w = e.weight;
}

void MVAHists::fill_(const Event & e, const double & mvaD, const int & n_btags)
{
  double w = e.weight;
  TopJet topjet = (*e.topjets).at(0);
  const Particle& lepton = e.get(h_primlep);

  double top_pt = topjet.v4().pt();
  double top_eta = topjet.v4().eta();

  double lept_pt = lepton.v4().pt();
  double lept_eta = lepton.v4().eta();

  // Now fill the histograms ...
  H_mvaD->Fill(mvaD, w);
  H_mvaD_20bins->Fill(mvaD, w);
  H_mvaD_rebin->Fill(mvaD, w);
  if(n_btags == 1) H_mvaD_1b->Fill(mvaD, w);
  else if(n_btags == 2) H_mvaD_2b->Fill(mvaD, w);
  
  // before any BDT cuts ...
  H_nocuts_top_pt_1000->Fill(top_pt, w);
  H_nocuts_top_pt_rebin->Fill(top_pt, w);
  H_nocuts_top_eta_1000->Fill(top_eta, w);
  H_nocuts_top_eta_rebin->Fill(top_eta, w);
  H_nocuts_top_pt_vs_eta_1000->Fill(top_eta, top_pt, w);
  H_nocuts_top_pt_vs_eta_rebin->Fill(top_eta, top_pt, w);

  H_nocuts_lept_pt->Fill(lept_pt, w);
  H_nocuts_lept_eta->Fill(lept_eta, w);

  // now apply a BDT cut ...
  if(mvaD > 0.2)
    {
      H_BDTcut_top_pt_1000->Fill(top_pt, w);
      H_BDTcut_top_pt_rebin->Fill(top_pt, w);
      H_BDTcut_top_eta_1000->Fill(top_eta, w);
      H_BDTcut_top_eta_rebin->Fill(top_eta, w);
      H_BDTcut_top_pt_vs_eta_1000->Fill(top_eta, top_pt, w);
      H_BDTcut_top_pt_vs_eta_rebin->Fill(top_eta, top_pt, w);

      H_BDTcut_lept_pt->Fill(lept_pt, w);
      H_BDTcut_lept_eta->Fill(lept_eta, w);
    }

  // now go through top_pt bins and write out BDT response for each

  if(top_pt > 200 && top_pt <= 300)
    {
      H_200to300_BDT->Fill(mvaD, w);
      H_200to300_BDT_rebin->Fill(mvaD, w);
      H_200to300_top_eta->Fill(top_eta, w);
      H_200to300_lept_pt->Fill(lept_pt, w);
      H_200to300_lept_eta->Fill(lept_eta, w);
    }
  else if(top_pt > 300 && top_pt <= 400)
    {
      H_300to400_BDT->Fill(mvaD, w);
      H_300to400_BDT_rebin->Fill(mvaD, w);
      H_300to400_top_eta->Fill(top_eta, w);
      H_300to400_lept_pt->Fill(lept_pt, w);
      H_300to400_lept_eta->Fill(lept_eta, w);
    }
  else if(top_pt > 400 && top_pt <= 600)
    {
      H_400to600_BDT->Fill(mvaD, w);
      H_400to600_BDT_rebin->Fill(mvaD, w);
      H_400to600_top_eta->Fill(top_eta, w);
      H_400to600_lept_pt->Fill(lept_pt, w);
      H_400to600_lept_eta->Fill(lept_eta, w);
    }
  else if(top_pt > 600 && top_pt <= 1200)
    {
      H_600to1200_BDT->Fill(mvaD, w);
      H_600to1200_BDT_rebin->Fill(mvaD, w);
      H_600to1200_top_eta->Fill(top_eta, w);
      H_600to1200_lept_pt->Fill(lept_pt, w);
      H_600to1200_lept_eta->Fill(lept_eta, w);
    }
}

MVAHists::~MVAHists(){}
