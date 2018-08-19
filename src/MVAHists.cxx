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
  H_mvaD = book<TH1F>("mvaD", "BDT discriminator", 40, -1, 1);
  H_mvaD_20bins = book<TH1F>("mvaD_20bins", "BDT discriminator", 20, -1, 1);
  H_mvaD_1b = book<TH1F>("mvaD_1btags", "BDT discriminator (1 b-tag)", 40, -1, 1);
  H_mvaD_2b = book<TH1F>("mvaD_2btags", "BDT discriminator (2 b-tags)", 40, -1, 1);

  // Top-jet distributions after BDT cut
  double rebin_pt[7] = {200,300,400,500,600,900,1200};
  double rebin_eta[6] = {0,0.5,1,1.5,2,2.5};
  double rebin_eta_lept[5] = {0,0.6,1.2,1.8,2.4};

  H_nocuts_top_pt_1000 = book<TH1F>("nocuts_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_nocuts_top_pt_rebin = book<TH1F>("nocuts_top_pt_rebin", "top-jet p_{T} [GeV]", 6, rebin_pt);
  H_nocuts_top_eta_1000 = book<TH1F>("nocuts_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_nocuts_top_eta_rebin = book<TH1F>("nocuts_top_eta_rebin", "top-jet #eta", 5, rebin_eta);
  H_nocuts_top_pt_vs_eta_1000 = book<TH2F>("nocuts_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_nocuts_top_pt_vs_eta_rebin = book<TH2F>("nocuts_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", 5, rebin_eta, 6, rebin_pt);

  H_nocuts_lept_pt = book<TH1F>("nocuts_lept_pt", "lepton p_{T} [GeV]", 10, 200, 1200);
  H_nocuts_lept_eta = book<TH1F>("nocuts_lept_eta", "lepton #eta", 4, rebin_eta_lept);

  H_BDTcut_top_pt_1000 = book<TH1F>("BDTcut_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_BDTcut_top_pt_rebin = book<TH1F>("BDTcut_top_pt_rebin", "top-jet p_{T} [GeV]", 6, rebin_pt);
  H_BDTcut_top_eta_1000 = book<TH1F>("BDTcut_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_BDTcut_top_eta_rebin = book<TH1F>("BDTcut_top_eta_rebin", "top-jet #eta", 5, rebin_eta);
  H_BDTcut_top_pt_vs_eta_1000 = book<TH2F>("BDTcut_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_BDTcut_top_pt_vs_eta_rebin = book<TH2F>("BDTcut_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", 5, rebin_eta, 6, rebin_pt);

  H_BDTcut_lept_pt = book<TH1F>("BDTcut_lept_pt", "lepton p_{T} [GeV]", 10, 200, 1200);
  H_BDTcut_lept_eta = book<TH1F>("BDTcut_lept_eta", "lepton #eta", 4, rebin_eta_lept);
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

}

MVAHists::~MVAHists(){}
