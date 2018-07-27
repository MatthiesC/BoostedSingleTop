#include "UHH2/BoostedSingleTop/include/HOTVRHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace uhh2;

/*
 * WARNING: Fill this Hists only after cuts are applied to ensure
 * there is >=1 HOTVR TopJet.
 * 
 * This Hists Class implements Histograms with informations about the
 * physics performance of the HOTVR algorithm.
 * 
 */
HOTVRHists::HOTVRHists(Context & ctx, const string & dirname, const boost::optional<TopJetId> &topjetid): 
  Hists(ctx, dirname), m_topjetid(topjetid)
{
  std::string dataset_name = ctx.get("dataset_version");

  // book all histograms here

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",           "number of top-tagged jets", 10, -0.5, 9.5);

  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",          "top-tags jet p_{T} [GeV]", 32, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",         "top-tags #eta", 30, -3, 3);
  Pt_vs_Eta_HotvrTopjets  = book<TH2F>("Pt_vs_Eta_HOTVR",   "top-tags jet p_{T} [GeV] vs. #eta", 32, 0, 1600, 30, -3, 3);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",           "top-tags jet mass [GeV]", 24, 120, 240);
  A_HotvrTopjets          = book<TH1F>("A_HOTVR",           "top-tags jet area", 20, 0, 10);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",        "number of sub-jets", 10, -0.5, 9.5);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",         "leading sub-jet f_{pT}", 20, 0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",       "sub-jets m_{ij}^{min} [GeV]", 40, 0, 200);
  Tau32_HotvrTopjets      = book<TH1F>("tau32_HOTVR",       "#tau_{3}/#tau_{2}", 50, 0, 1);
  Pt_HotvrTopjets_Sub1    = book<TH1F>("Pt_HOTVR_Subjet1",  "p_{T}^{subjet1} [GeV]", 40, 0, 1000);
  A_HotvrTopjets_Sub1     = book<TH1F>("A_HOTVR_Subjet1",   "A_{subjet1}", 50, 0, 5);
  Pt_HotvrTopjets_Sub2    = book<TH1F>("Pt_HOTVR_Subjet2",  "p_{T}^{subjet2} [GeV]", 40, 0, 1000);
  A_HotvrTopjets_Sub2     = book<TH1F>("A_HOTVR_Subjet2",   "A_{subjet2}", 50, 0, 5);
  Pt_HotvrTopjets_Sub3    = book<TH1F>("Pt_HOTVR_Subjet3",  "p_{T}^{subjet3} [GeV]", 40, 0, 1000);
  A_HotvrTopjets_Sub3     = book<TH1F>("A_HOTVR_Subjet3",   "A_{subjet3}", 50, 0, 5);

  DeltaR_Lepton_HotvrTopjets   = book<TH1F>("DeltaR_Lepton_HOTVR",   "#DeltaR(lepton, top-tag)", 25, 0, 5);
  DeltaPhi_Lepton_HotvrTopjets = book<TH1F>("DeltaPhi_Lepton_HOTVR", "#Delta#phi(lepton, top-tag)", 32, 0, M_PI);

  Pt_HotvrTopjet1         = book<TH1F>("Pt_HOTVR1",         "top-tag jet p_{T} [GeV]", 32, 0, 1600);
  Eta_HotvrTopjet1        = book<TH1F>("Eta_HOTVR1",        "top-tag #eta", 30, -3, 3);
  Pt_vs_Eta_HotvrTopjet1  = book<TH2F>("Pt_vs_Eta_HOTVR1",  "top-tag jet p_{T} [GeV] vs. #eta", 32, 0, 1600, 30, -3, 3);
  M_HotvrTopjet1          = book<TH1F>("M_HOTVR1",          "top-tag jet mass [GeV]", 24, 120, 240);
  A_HotvrTopjet1          = book<TH1F>("A_HOTVR1",          "top-tag jet area", 20, 0, 10);
  NSub_HotvrTopjet1       = book<TH1F>("NSub_HOTVR1",       "number of sub-jets", 10, -0.5, 9.5);
  Fpt_HotvrTopjet1        = book<TH1F>("Fpt_HOTVR1",        "leading sub-jet f_{pT}", 20, 0, 1);
  Mpair_HotvrTopjet1      = book<TH1F>("Mpair_HOTVR1",      "sub-jets m_{ij}^{min} [GeV]", 40,  0, 200);
  Tau32_HotvrTopjet1      = book<TH1F>("tau32_HOTVR1",      "#tau_{3}/#tau_{2}", 50,  0, 1);
  Pt_HotvrTopjet1_Sub1    = book<TH1F>("Pt_HOTVR1_Subjet1", "p_{T}^{subjet1} [GeV]", 40, 0, 1000);
  A_HotvrTopjet1_Sub1     = book<TH1F>("A_HOTVR1_Subjet1",   "A_{subjet1}", 50, 0, 5);
  Pt_HotvrTopjet1_Sub2    = book<TH1F>("Pt_HOTVR1_Subjet2", "p_{T}^{subjet2} [GeV]", 40, 0, 1000);
  A_HotvrTopjet1_Sub2     = book<TH1F>("A_HOTVR1_Subjet2",   "A_{subjet2}", 50, 0, 5);
  Pt_HotvrTopjet1_Sub3    = book<TH1F>("Pt_HOTVR1_Subjet3", "p_{T}^{subjet3} [GeV]", 40, 0, 1000);
  A_HotvrTopjet1_Sub3     = book<TH1F>("A_HOTVR1_Subjet3",   "A_{subjet3}", 50, 0, 5);

  DeltaR_Lepton_HotvrTopjet1   = book<TH1F>("DeltaR_Lepton_HOTVR1",   "#DeltaR(lepton, top-tag)", 25, 0, 5);
  DeltaPhi_Lepton_HotvrTopjet1 = book<TH1F>("DeltaPhi_Lepton_HOTVR1", "#Delta#phi(lepton, top-tag)", 32, 0, M_PI);


  double pt_xbins[4] =   {200, 300, 400, 1600};
  Pt_rebin_HotvrTopjets   = book<TH1F>("Pt_rebin_HOTVR",   "p_{T}^{HOTVR} [GeV]", 3, pt_xbins); 
  double eta_xbins[5] =  {-2.5, -1.479, 0, 1.479, 2.5};
  EtaAbs_HotvrTopjets     = book<TH1F>("EtaAbs_HOTVR",     "|#eta|^{HOTVR}", 4, eta_xbins);
  Pt_vs_Eta_HotvrRebin    = book<TH2F>("Pt_vs_Eta_HOTVR_Rebin",   "p_{T}^{HOTVR} vs. #eta^{HOTVR}", 3, pt_xbins, 4, eta_xbins);
  if (m_topjetid)
    {
      NLeadingTopjet      = book<TH1F>("NLeadingTopjet",    "", 2, -0.5, 1.5);
    }

  // b-jets
  A_ak4                   = book<TH1F>("A_ak4",            "AK4 jet area", 50, 0, 5);
  Pt_ak4                  = book<TH1F>("pT_ak4",           "AK4 jet p_{T} [GeV]", 50, 0, 1000);
  N_bjets_loose           = book<TH1F>("N_bjets_loose",    "number of loose CSV b-tags", 10, -0.5, 9.5);
  N_bjets_medium          = book<TH1F>("N_bjets_medium",   "number of medium CSV b-tags", 10, -0.5, 9.5);
  N_bjets_tight           = book<TH1F>("N_bjets_tight",    "number of tight CSV b-tags", 10, -0.5, 9.5);

  // all tight b-jets
  bjets_pt          = book<TH1F>("BJets_Pt", "b-tags jet p_{T} [GeV]", 25, 0, 500);
  bjets_eta         = book<TH1F>("BJets_Eta", "b-tags #eta", 30, -3, 3);
  bjets_pt_vs_eta   = book<TH2F>("BJets_Pt_vs_Eta", "b-tags jet p_{T} [GeV] vs. #eta", 25, 0, 500, 30, -3, 3);
  bjets_m           = book<TH1F>("BJets_M", "b-tags jet mass [GeV]", 20, 0, 40);
  bjets_a           = book<TH1F>("BJets_A", "b-tags jet area", 50, 0, 5);

  // if N_bjets == 1
  bjet1_pt          = book<TH1F>("BJet1_Pt", "p_{T}^{b-tag} [GeV] for N_{b-tags} = 1", 25, 0, 500);
  bjet1_eta         = book<TH1F>("BJet1_Eta", "#eta_{b-tag} for N_{b-tags} = 1", 30, -3, 3);
  bjet1_pt_vs_eta   = book<TH2F>("BJet1_Pt_vs_Eta", "p_{T}^{b-tag} vs. #eta^{b-tags} for N_{b-tags} = 1", 25, 0, 500, 30, -3, 3);
  bjet1_m           = book<TH1F>("BJet1_M", "M_{b-tag} [GeV] for N_{b-tags} = 1", 20, 0, 40);
  bjet1_a           = book<TH1F>("BJet1_A", "A_{b-tag} for N_{b-tags} = 1", 50, 0, 5);
  bjet1_deltaR_top   = book<TH1F>("BJet1_deltaR_top", "#DeltaR(top-tag, b-tag) for N_{b-tags} = 1", 25, 0, 5);
  bjet1_deltaPhi_top = book<TH1F>("BJet1_deltaPhi_top", "#Delta#phi(top-tag, b-tag) for N_{b-tags} = 1", 32, 0, M_PI);
  bjet1_deltaR_lep   = book<TH1F>("BJet1_deltaR_lep", "#DeltaR(lepton, b-tag) for N_{b-tags} = 1", 25, 0, 5);
  bjet1_deltaPhi_lep = book<TH1F>("BJet1_deltaPhi_lep", "#Delta#phi(lepton, b-tag) for N_{b-tags} = 1", 32, 0, M_PI);

  // if N_bjets == 2
  bjet2_1_pt        = book<TH1F>("BJet2_1_Pt", "p_{T}^{b-tag1} [GeV] for N_{b-tags} = 2", 50, 0, 500);
  bjet2_2_pt        = book<TH1F>("BJet2_2_Pt", "p_{T}^{b-tag2} [GeV] for N_{b-tags} = 2", 50, 0, 500);
  bjet2_1_eta       = book<TH1F>("BJet2_1_Eta", "#eta_{b-tag1} for N_{b-tags} = 2", 30, -3, 3);
  bjet2_2_eta       = book<TH1F>("BJet2_2_Eta", "#eta_{b-tag2} for N_{b-tags} = 2", 30, -3, 3);
  bjet2_1_pt_vs_eta = book<TH2F>("BJet2_1_Pt_vs_Eta", "p_{T}^{b-tag1} vs. #eta^{b-tag1} for N_{b-tags} = 2", 20, 0, 500, 30, -3, 3);
  bjet2_2_pt_vs_eta = book<TH2F>("BJet2_2_Pt_vs_Eta", "p_{T}^{b-tag2} vs. #eta^{b-tag2} for N_{b-tags} = 2", 20, 0, 500, 30, -3, 3);
  bjet2_1_m         = book<TH1F>("BJet2_1_M", "M_{b-tag1} [GeV] for N_{b-tag1} = 2", 20, 0, 40);
  bjet2_2_m         = book<TH1F>("BJet2_2_M", "M_{b-tag2} [GeV] for N_{b-tag2} = 2", 20, 0, 40);
  bjet2_1_a         = book<TH1F>("BJet2_1_A", "A_{b-tag1} for N_{b-tag1} = 2", 50, 0, 5);
  bjet2_2_a         = book<TH1F>("BJet2_2_A", "A_{b-tag2} for N_{b-tag2} = 2", 50, 0, 5);

  // if N_bjets == 3
  bjet3_1_pt        = book<TH1F>("BJet3_1_Pt", "p_{T}^{b-tag1} [GeV] for N_{b-tags} = 3", 50, 0, 500);
  bjet3_2_pt        = book<TH1F>("BJet3_2_Pt", "p_{T}^{b-tag2} [GeV] for N_{b-tags} = 3", 50, 0, 500);
  bjet3_3_pt        = book<TH1F>("BJet3_3_Pt", "p_{T}^{b-tag3} [GeV] for N_{b-tags} = 3", 50, 0, 500);
  bjet3_1_eta       = book<TH1F>("BJet3_1_Eta", "#eta_{b-tag1} for N_{b-tags} = 3", 30, -3, 3);
  bjet3_2_eta       = book<TH1F>("BJet3_2_Eta", "#eta_{b-tag2} for N_{b-tags} = 3", 30, -3, 3);
  bjet3_3_eta       = book<TH1F>("BJet3_3_Eta", "#eta_{b-tag3} for N_{b-tags} = 3", 30, -3, 3);
  bjet3_1_pt_vs_eta = book<TH2F>("BJet3_1_Pt_vs_Eta", "p_{T}^{b-tag1} vs. #eta^{b-tag1} for N_{b-tags} = 3", 20, 0, 500, 30, -3, 3);
  bjet3_2_pt_vs_eta = book<TH2F>("BJet3_2_Pt_vs_Eta", "p_{T}^{b-tag2} vs. #eta^{b-tag2} for N_{b-tags} = 3", 20, 0, 500, 30, -3, 3);
  bjet3_3_pt_vs_eta = book<TH2F>("BJet3_3_Pt_vs_Eta", "p_{T}^{b-tag3} vs. #eta^{b-tag3} for N_{b-tags} = 3", 20, 0, 500, 30, -3, 3);
  bjet3_1_m         = book<TH1F>("BJet3_1_M", "M_{b-tag1} [GeV] for N_{b-tag1} = 3", 20, 0, 40);
  bjet3_2_m         = book<TH1F>("BJet3_2_M", "M_{b-tag2} [GeV] for N_{b-tag2} = 3", 20, 0, 40);
  bjet3_3_m         = book<TH1F>("BJet3_3_M", "M_{b-tag3} [GeV] for N_{b-tag3} = 3", 20, 0, 40);
  bjet3_1_a         = book<TH1F>("BJet3_1_A", "A_{b-tag1} for N_{b-tag1} = 3", 50, 0, 5);
  bjet3_2_a         = book<TH1F>("BJet3_2_A", "A_{b-tag2} for N_{b-tag2} = 3", 50, 0, 5);
  bjet3_3_a         = book<TH1F>("BJet3_3_A", "A_{b-tag3} for N_{b-tag3} = 3", 50, 0, 5);

  // TO DO: add hists for Delta Phi and Delta R between b-jet and top-jet, as well as b-jet and lepton, as well as min(Delta(lepton,b),Delta(top,b))
  // need to update HOTVRHists.h

  // if (N_bjets >= 1 && N_topjets == 1)
  //   {
  // top-jet
  h_top_pt = book<TH1F>("top_pt", "top-tag jet p_{T} [GeV]", 32, 0, 1600);
  h_top_eta = book<TH1F>("top_eta", "top-tag #eta", 30, -3, 3);
  h_top_pt_vs_eta = book<TH2F>("top_pt_vs_eta", "top-tag jet p_{T} vs. #eta", 32, 0, 1600, 30, -3, 3);
  h_top_m = book<TH1F>("top_m", "top-tag jet mass [GeV]", 24,  120, 240);
  h_top_a = book<TH1F>("top_a", "top-tag jet area", 20, 0, 10);
  h_top_nsub = book<TH1F>("top_nsub", "number of sub-jets", 10, -0.5, 9.5);
  h_top_fpt = book<TH1F>("top_fpt", "leading sub-jet f_{pT}", 20, 0, 1);
  h_top_mpair = book<TH1F>("top_mpair", "sub-jets m_{ij}^{min} [GeV]", 40, 0, 200);
  h_top_tau32 = book<TH1F>("top_tau32", "#tau_{3}/#tau_{2}", 50, 0, 1);

  // subjets (just filled if n_subjets >= 3)
  h_sub1_pt = book<TH1F>("sub1_pt", "p_{T}^{subjet1} [GeV]", 40, 0, 1000);
  h_sub1_a = book<TH1F>("sub1_a", "A_{subjet1}", 50, 0, 5);
  h_sub2_pt = book<TH1F>("sub2_pt", "p_{T}^{subjet2} [GeV]", 40, 0, 1000);
  h_sub2_a = book<TH1F>("sub2_a", "A_{subjet2}", 50, 0, 5);
  h_sub3_pt = book<TH1F>("sub3_pt", "p_{T}^{subjet3} [GeV]", 40, 0, 1000);
  h_sub3_a = book<TH1F>("sub3_a", "A_{subjet3}", 50, 0, 5);

  // top-jet/b-jet
  h_deltaR_top_bjet = book<TH1F>("deltaR_top_bjet", "#DeltaR(top-tag, b-tag)", 25, 0, 5);
  h_deltaPhi_top_bjet = book<TH1F>("deltaPhi_top_bjet", "#Delta#phi(top-tag, b-tag)", 32, 0, M_PI);
  // Wass
  h_wass_m = book<TH1F>("wass_m", "m^{reco}_{W} [GeV]", 50, 60, 160);
  h_wass_pt = book<TH1F>("wass_pt", "p_{T}^{W,reco} [GeV]", 32, 0, 1600);
  h_wass_eta = book<TH1F>("wass_eta", "#eta^{reco}_{W}", 30, -3, 3);
  // --- LEPTON
  h_lep_pt = book<TH1F>("lep_pt", "lepton p_{T} [GeV]", 32, 0, 1600);
  h_lep_eta = book<TH1F>("lep_eta", "lepton #eta", 30, -3, 3);
  h_lep_pt_vs_eta = book<TH2F>("lep_pt_vs_eta", "lepton p_{T} [GeV] vs. #eta", 32, 0, 1600, 30, -3, 3);
  // lepton/top-jet
  double dR_rebin[18] = {0,0.5,0.8,1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3,3.2,3.5,4,5};
  unsigned int dR_size = sizeof(dR_rebin)/sizeof(*dR_rebin)-1;
  h_deltaR_lep_top = book<TH1F>("deltaR_lep_top", "#DeltaR(lepton, top-tag)", 25, 0, 5);
  h_deltaPhi_lep_top = book<TH1F>("deltaPhi_lep_top", "#Delta#phi(lepton, top-tag)", 32, 0, M_PI);
  // b-jet to others
  h_deltaR_lep_bjet = book<TH1F>("deltaR_lep_bjet", "#DeltaR(lepton, b_{#DeltaR}^{tag})", 25, 0, 5);
  h_deltaPhi_lep_bjet = book<TH1F>("deltaPhi_lep_bjet", "#Delta#phi(lepton, b_{#DeltaR}^{tag})", 32, 0, M_PI);
  h_deltaPhis_topB_vs_lepB = book<TH2F>("deltaPhis_topB_vs_lepB", "#Delta#phi(top-tag, b_{#DeltaR}^{tag}) vs. #Delta#phi(lepton, b_{#DeltaR}^{tag})", 32, 0, M_PI, 32, 0, M_PI);
  // Lepton+Nu+BJet Reconstruction
  double m_rebin[14] = {0,80,120,160,200,240,280,320,360,400,450,500,600,1000};
  unsigned int m_size = sizeof(m_rebin)/sizeof(*m_rebin)-1;
  h_LepNuBSystem_M_combined = book<TH1F>("LepNuBSystem_M_combined", "reco. inv. lepton-neutrino-b_{#DeltaR}^{tag} mass [GeV]", 40, 0, 800);
  h_LepNuBSystem_Pt_combined = book<TH1F>("LepNuBSystem_Pt_combined", "reco. inv. lepton-neutrino-b_{#DeltaR}^{tag} p_{T} [GeV]", 32, 0, 1600);
  // Pt Balance
  double bal_rebin[10] = {-1.5,-1.0,-0.8,-0.6,-0.4,-0.2,0,0.25,0.5,1.5};
  unsigned int bal_bins = sizeof(bal_rebin)/sizeof(*bal_rebin)-1;
  h_ptBalance_combined = book<TH1F>("ptBalance_combined", "[ p_{T}^{W,reco} - p_{T}^{top-tag} ] / p_{T}^{top-tag}", 30, -1.5, 1.5);
  // Cuts
  double pt_rebin[6] = {200,300,400,600,1000,1600};
  unsigned int Nbins = sizeof(pt_rebin)/sizeof(*pt_rebin)-1; // so I don't have to compute that always by myself
  // cut combinatorics
  /*survivesMcut_top_pt = book<TH1F>("survivesMcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: M_{mean}^{reco}(lep, #nu, b-jet)", 32, 0, 1600);
    survivesDeltaRCut_top_pt = book<TH1F>("DeltaRCut_top_pt", "p_{T}^{top-jet} [GeV] after cut: #DeltaR(lep, b-jet)", 32, 0, 1600);*/
  /*
  survivesMcut_top_pt = book<TH1F>("survivesMcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: Mass", 32, 0, 1600);
  survivesMcut_top_pt_rebin = book<TH1F>("survivesMcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: Mass", Nbins, pt_rebin);
  survivesMcut_deltaR = book<TH1F>("survivesMcut_deltaR", "#DeltaR(lep, b-jet) after cut: Mass", dR_size, dR_rebin);
  survivesMcut_PtB = book<TH1F>("survivesMcut_PtB", "[ p_{T}^{reco}(W) - p_{T}^{reco}(t) ] / p_{T}^{reco}(t) after cut: Mass", bal_bins, bal_rebin);

  survivesMcut_survivesDeltaRcut_top_pt = book<TH1F>("survivesMcut_survivesDeltaRcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: Mass, DeltaR", 32, 0, 1600);
  survivesMcut_survivesDeltaRcut_top_pt_rebin = book<TH1F>("survivesMcut_survivesDeltaRcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: Mass, DeltaR", Nbins, pt_rebin);
  survivesMcut_survivesDeltaRcut_PtB = book<TH1F>("survivesMcut_survivesDeltaRcut_PtB", "[ p_{T}^{reco}(W) - p_{T}^{reco}(t) ] / p_{T}^{reco}(t) after cut: Mass, DeltaR", bal_bins, bal_rebin);

  survivesMcut_survivesPtBcut_top_pt = book<TH1F>("survivesMcut_survivesPtBcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: Mass, PtB", 32, 0, 1600);
  survivesMcut_survivesPtBcut_top_pt_rebin = book<TH1F>("survivesMcut_survivesPtBcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: Mass, PtB", Nbins, pt_rebin);
  survivesMcut_survivesPtBcut_deltaR = book<TH1F>("survivesMcut_survivesPtBcut_deltaR", "#DeltaR(lep, b-jet) after cut: Mass, PtB", dR_size, dR_rebin);

  survivesDeltaRcut_top_pt = book<TH1F>("survivesDeltaRcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: DeltaR", 32, 0, 1600);
  survivesDeltaRcut_top_pt_rebin = book<TH1F>("survivesDeltaRcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: DeltaR", Nbins, pt_rebin);
  survivesDeltaRcut_Mass = book<TH1F>("survivesDeltaRcut_Mass", "M_{mean}^{reco}(lep, #nu, b-jet) [GeV], #nu-sol^{s} combined, after cut: DeltaR", m_size, m_rebin);
  survivesDeltaRcut_PtB = book<TH1F>("survivesDeltaRcut_PtB", "[ p_{T}^{reco}(W) - p_{T}^{reco}(t) ] / p_{T}^{reco}(t) after cut: DeltaR", bal_bins, bal_rebin);

  survivesDeltaRcut_survivesMcut_top_pt = book<TH1F>("survivesDeltaRcut_survivesMcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: DeltaR, Mass", 32, 0, 1600);
  survivesDeltaRcut_survivesMcut_top_pt_rebin = book<TH1F>("survivesDeltaRcut_survivesMcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: DeltaR, Mass", Nbins, pt_rebin);
  survivesDeltaRcut_survivesMcut_PtB = book<TH1F>("survivesDeltaRcut_survivesMcut_PtB", "[ p_{T}^{reco}(W) - p_{T}^{reco}(t) ] / p_{T}^{reco}(t) after cut: DeltaR, Mass", bal_bins, bal_rebin);

  survivesDeltaRcut_survivesPtBcut_top_pt = book<TH1F>("survivesDeltaRcut_survivesPtBcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: DeltaR, PtB", 32, 0, 1600);
  survivesDeltaRcut_survivesPtBcut_top_pt_rebin = book<TH1F>("survivesDeltaRcut_survivesPtBcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: DeltaR, PtB", Nbins, pt_rebin);
  survivesDeltaRcut_survivesPtBcut_Mass = book<TH1F>("survivesDeltaRcut_survivesPtBcut_Mass", "M_{mean}^{reco}(lep, #nu, b-jet) [GeV], #nu-sol^{s} combined, after cut: DeltaR, PtB", m_size, m_rebin);

  survivesPtBcut_top_pt = book<TH1F>("survivesPtBcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: PtB", 32, 0, 1600);
  survivesPtBcut_top_pt_rebin = book<TH1F>("survivesPtBcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: PtB", Nbins, pt_rebin);
  survivesPtBcut_Mass = book<TH1F>("survivesPtBcut_Mass", "M_{mean}^{reco}(lep, #nu, b-jet) [GeV], #nu-sol^{s} combined, after cut: PtB", m_size, m_rebin);
  survivesPtBcut_deltaR = book<TH1F>("survivesPtBcut_deltaR", "#DeltaR(lep, b-jet) after cut: PtB", dR_size, dR_rebin);

  survivesPtBcut_survivesMcut_top_pt = book<TH1F>("survivesPtBcut_survivesMcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: PtB, Mass", 32, 0, 1600);
  survivesPtBcut_survivesMcut_top_pt_rebin = book<TH1F>("survivesPtBcut_survivesMcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: PtB, Mass", Nbins, pt_rebin);
  survivesPtBcut_survivesMcut_deltaR = book<TH1F>("survivesPtBcut_survivesMcut_deltaR", "#DeltaR(lep, b-jet) after cut: PtB, Mass", dR_size, dR_rebin);

  survivesPtBcut_survivesDeltaRcut_top_pt = book<TH1F>("survivesPtBcut_survivesDeltaRcut_top_pt", "p_{T}^{top-jet} [GeV] after cut: PtB, DeltaR", 32, 0, 1600);
  survivesPtBcut_survivesDeltaRcut_top_pt_rebin = book<TH1F>("survivesPtBcut_survivesDeltaRcut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: PtB, DeltaR", Nbins, pt_rebin);
  survivesPtBcut_survivesDeltaRcut_Mass = book<TH1F>("survivesPtBcut_survivesDeltaRcut_Mass", "M_{mean}^{reco}(lep, #nu, b-jet) [GeV], #nu-sol^{s} combined, after cut: PtB, DeltaR", m_size, m_rebin);
  */
  // all CUTS
  AllCuts_top_pt = book<TH1F>("AllCuts_top_pt", "all cuts: top-tag jet p_{T} [GeV]", 32, 0, 1600);
  AllCuts_top_pt_rebin = book<TH1F>("AllCuts_top_pt_rebin", "all cuts: top-tag jet p_{T} [GeV]", Nbins, pt_rebin);

  //AllCuts_2Btags = book<TH1F>("AllCuts_2Btags", "all cuts: top-tag jet p_{T} [GeV]", Nbins, pt_rebin);
  //AllCuts_1Btag = book<TH1F>("AllCuts_1Btag", "all cuts: top-tag jet p_{T} [GeV]", Nbins, pt_rebin);
  //AllCuts_1Btag_v2 = book<TH1F>("AllCuts_1Btag_v2", "all cuts: top-tag jet p_{T} [GeV]", Nbins, pt_rebin);

  /*h_BothCuts_top_pt = book<TH1F>("BothCuts_top_pt", "p_{T}^{top-jet} [GeV] after both cuts combined", 32, 0, 1600);
    h_MCut_top_pt_rebin = book<TH1F>("MCut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: M_{mean}^{reco}(lep, #nu, b-jet)", Nbins, pt_rebin);
    h_DeltaRCut_top_pt_rebin = book<TH1F>("DeltaRCut_top_pt_rebin", "p_{T}^{top-jet} [GeV] after cut: #DeltaR(lep, b-jet)", Nbins, pt_rebin);
    h_BothCuts_top_pt_rebin = book<TH1F>("BothCuts_top_pt_rebin", "p_{T}^{top-jet} [GeV] after both cuts combined", Nbins, pt_rebin);

    h_survivesMcut_deltaR_lep_bjet = book<TH1F>("SurvivesMCut_deltaR_lep_bjet", "#DeltaR(lep, top-jet) after cut: M_{mean}^{reco}(lep, #nu, b-jet)", dR_size, dR_rebin);
    h_survivesMcut_deltaPhi_lep_bjet = book<TH1F>("SurvivesMCut_deltaPhi_lep_bjet", "#Delta#phi(lep, b-jet) after cut: M_{mean}^{reco}(lep, #nu, b-jet)", 32, 0, M_PI);
    h_survivesDeltaRcut_M_combined = book<TH1F>("SurvivesDeltaRCut_M_combined", "M_{mean}^{reco}(lep, #nu, b-jet) [GeV], #nu-sol^{s} combined, after cut: #DeltaR(lep, b-jet)", m_size, m_rebin);*/
  //   }
}

void HOTVRHists::fill(const Event & event) {  

  double weight = event.weight; // event weight

  vector<TopJet> hotvrJets = *event.topjets; // HOTVR jets
  vector<Muon> muons = *event.muons;
  sort_by_pt<Muon>(muons);
  vector<Electron> electrons = *event.electrons;
  sort_by_pt<Electron>(electrons);
  vector<Jet> jets = *event.jets; // all jets
  vector<Jet> bjets; // b-tagged jets (tight)
  vector<TopJet> toptaggedjets; // HOTVR top-tagged jets
  for (TopJet hotvrjet : hotvrJets)
    {
      if (m_topjetid && (*m_topjetid)(hotvrjet, event))
	{
	  toptaggedjets.push_back(hotvrjet); // fill vector with all top-tagged HOTVR jets
	}
    }
  sort_by_pt<TopJet>(toptaggedjets);
  LorentzVector met = (*event.met).v4(); // used for the neutrino reconstruction and then mass reconstruction of the neutrino/lepton/b-jet system

  //b-tags
  const CSVBTag btag_loose(CSVBTag::WP_LOOSE);
  const CSVBTag btag_medium(CSVBTag::WP_MEDIUM);
  const CSVBTag btag_tight(CSVBTag::WP_TIGHT);

  // fill bjet hists

  int n_loose  = 0;
  int n_medium = 0;
  int n_tight  = 0;

  for (Jet jet : jets)
    {
      A_ak4->Fill(jet.jetArea(), weight);
      Pt_ak4->Fill(jet.pt(), weight);
      if (btag_tight(jet, event))
	{
	  ++n_tight;
	  bjets.push_back(jet); // fill b-jet vector (all bjets with at least looseID count)
	}
      if (btag_medium(jet, event))
	{
	  ++n_medium;
	}
      if (btag_loose(jet, event))
	{
	  ++n_loose;
	}
    }

  sort_by_pt<Jet>(bjets); // sorting AFTER filling the vector :-)

  // TO DO: sort bjets by R-distance to top-jet in the N_top == 1 case

  if(toptaggedjets.size() == 1 && electrons.size() + muons.size() == 1)
    {
      // Fill number of b-tagged ak4 jets
      N_bjets_loose->Fill(n_loose, weight);
      N_bjets_medium->Fill(n_medium, weight);
      N_bjets_tight->Fill(n_tight, weight);

      LorentzVector lepton;
      if (electrons.size() == 1) lepton = electrons.at(0).v4();
      else if (muons.size() == 1) lepton = muons.at(0).v4();
      LorentzVector topjet = toptaggedjets.at(0).v4();
      for (Jet bjet : bjets)
	{
	  bjets_pt->Fill(bjet.v4().pt(), weight);
	  bjets_eta->Fill(bjet.v4().eta(), weight);
	  bjets_pt_vs_eta->Fill(bjet.v4().pt(), bjet.v4().eta(), weight);
	  bjets_m->Fill(bjet.v4().M(), weight);
	  bjets_a->Fill(bjet.jetArea(), weight);
	}
      if (bjets.size() == 1)
	{
	  bjet1_pt->Fill(bjets.at(0).v4().pt(), weight);
	  bjet1_eta->Fill(bjets.at(0).v4().eta(), weight);
	  bjet1_pt_vs_eta->Fill(bjets.at(0).v4().pt(), bjets.at(0).v4().eta(), weight);
	  bjet1_m->Fill(bjets.at(0).v4().M(), weight);
	  bjet1_a->Fill(bjets.at(0).jetArea(), weight);
	  bjet1_deltaR_top->Fill(deltaR(bjets.at(0).v4(), topjet), weight);
	  bjet1_deltaPhi_top->Fill(deltaPhi(bjets.at(0).v4(), topjet), weight);
	  bjet1_deltaR_lep->Fill(deltaR(bjets.at(0).v4(), lepton), weight);
	  bjet1_deltaPhi_lep->Fill(deltaPhi(bjets.at(0).v4(), lepton), weight);
	}
      else if (bjets.size() == 2)
	{
	  // sorting b-jets by R-to-top: (those who are the farest away come first)
	  if (deltaR(bjets.at(0).v4(), topjet) < deltaR(bjets.at(1).v4(), topjet)) swap(bjets.at(0), bjets.at(1));

	  bjet2_1_pt->Fill(bjets.at(0).v4().pt(), weight);
	  bjet2_2_pt->Fill(bjets.at(1).v4().pt(), weight);
	  bjet2_1_eta->Fill(bjets.at(0).v4().eta(), weight);
	  bjet2_2_eta->Fill(bjets.at(1).v4().eta(), weight);
	  bjet2_1_pt_vs_eta->Fill(bjets.at(0).v4().pt(), bjets.at(0).v4().eta(), weight);
	  bjet2_2_pt_vs_eta->Fill(bjets.at(1).v4().pt(), bjets.at(1).v4().eta(), weight);
	  bjet2_1_m->Fill(bjets.at(0).v4().M(), weight);
	  bjet2_2_m->Fill(bjets.at(1).v4().M(), weight);
	  bjet2_1_a->Fill(bjets.at(0).jetArea(), weight);
	  bjet2_2_a->Fill(bjets.at(1).jetArea(), weight);
	}
      else if (bjets.size() == 3)
	{
	  // sorting b-jets by R-to-top:
	  if (deltaR(bjets.at(2).v4(), topjet) > deltaR(bjets.at(1), topjet)) swap(bjets.at(2), bjets.at(1));
	  if (deltaR(bjets.at(1).v4(), topjet) > deltaR(bjets.at(0), topjet)) swap(bjets.at(1), bjets.at(0));
	  if (deltaR(bjets.at(2).v4(), topjet) > deltaR(bjets.at(1), topjet)) swap(bjets.at(2), bjets.at(1)); // required 2 times!

	  bjet3_1_pt->Fill(bjets.at(0).v4().pt(), weight);
	  bjet3_2_pt->Fill(bjets.at(1).v4().pt(), weight);
	  bjet3_3_pt->Fill(bjets.at(2).v4().pt(), weight);
	  bjet3_1_eta->Fill(bjets.at(0).v4().eta(), weight);
	  bjet3_2_eta->Fill(bjets.at(1).v4().eta(), weight);
	  bjet3_3_eta->Fill(bjets.at(2).v4().eta(), weight);
	  bjet3_1_pt_vs_eta->Fill(bjets.at(0).v4().pt(), bjets.at(0).v4().eta(), weight);
	  bjet3_2_pt_vs_eta->Fill(bjets.at(1).v4().pt(), bjets.at(1).v4().eta(), weight);
	  bjet3_3_pt_vs_eta->Fill(bjets.at(2).v4().pt(), bjets.at(2).v4().eta(), weight);
	  bjet3_1_m->Fill(bjets.at(0).v4().M(), weight);
	  bjet3_2_m->Fill(bjets.at(1).v4().M(), weight);
	  bjet3_3_m->Fill(bjets.at(2).v4().M(), weight);
	  bjet3_1_a->Fill(bjets.at(0).jetArea(), weight);
	  bjet3_2_a->Fill(bjets.at(1).jetArea(), weight);
	  bjet3_3_a->Fill(bjets.at(2).jetArea(), weight);
	}
    }
  // fill Hists for 1Lep/1bjet/1top events
  if (toptaggedjets.size() == 1 && bjets.size() >= 1)
    {

      TopJet topjet = toptaggedjets.at(0);
      vector<Jet> subjets = topjet.subjets();

      double fpt = -1;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = (subjets.at(0).v4().pt()) / topjet.v4().pt();
	}
      double mpair = -1;
      double pt_sub1 = -1;
      double a_sub1 = -1;
      double pt_sub2 = -1;
      double a_sub2 = -1;
      double pt_sub3 = -1;
      double a_sub3 = -1;

      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	  pt_sub1 = subjets.at(0).v4().pt();
	  a_sub1  = subjets.at(0).jetArea();
	  pt_sub2 = subjets.at(1).v4().pt();
	  a_sub2  = subjets.at(1).jetArea();
	  pt_sub3 = subjets.at(2).v4().pt();
	  a_sub3  = subjets.at(2).jetArea();
	}

      // top-jet
      h_top_pt->Fill(topjet.v4().pt(), weight);
      h_top_eta->Fill(topjet.v4().eta(), weight);
      h_top_pt_vs_eta->Fill(topjet.v4().pt(), topjet.v4().eta(), weight);
      h_top_m->Fill(topjet.v4().M(), weight);
      h_top_a->Fill(topjet.jetArea(), weight);
      h_top_nsub->Fill(subjets.size(), weight);
      h_top_fpt->Fill(fpt, weight);
      h_top_mpair->Fill(mpair, weight);
      h_top_tau32->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);

      // subjets (just filled if n_subjets >= 3)
      h_sub1_pt->Fill(pt_sub1, weight);
      h_sub1_a->Fill(a_sub1, weight);
      h_sub2_pt->Fill(pt_sub2, weight);
      h_sub2_a->Fill(a_sub2, weight);
      h_sub3_pt->Fill(pt_sub3, weight);
      h_sub3_a->Fill(a_sub3, weight);

      // top-jet/b-jet
      double deltaR_top_bjet = deltaR(topjet.v4(), bjets.at(0).v4());
      double deltaPhi_top_bjet = deltaPhi(topjet.v4(), bjets.at(0).v4());
      h_deltaR_top_bjet->Fill(deltaR_top_bjet, weight);
      h_deltaPhi_top_bjet->Fill(deltaPhi_top_bjet, weight);

      if (electrons.size() + muons.size() == 1)
	{
	  LorentzVector lepton;
	  if (electrons.size() == 1) lepton = electrons.at(0).v4();
	  else if (muons.size() == 1) lepton = muons.at(0).v4();
	  vector<LorentzVector> nus = NeutrinoReconstruction(lepton, met);
	  vector<double> M_LepNuB, Pt_LepNuB, Eta_LepNuB;
	  vector<double> M_Wass, Pt_Wass, Eta_Wass;
	  for (LorentzVector neutrino : nus)
	    {
	      // LepNuB reconstruction
	      M_LepNuB.push_back((lepton + bjets.at(0).v4() + neutrino).M()); 
	      Pt_LepNuB.push_back((lepton + bjets.at(0).v4() + neutrino).Pt());
	      Eta_LepNuB.push_back((lepton + bjets.at(0).v4() + neutrino).Eta());
	      // Wass reconstruction
	      M_Wass.push_back((lepton + neutrino).M());
	      Pt_Wass.push_back((lepton + neutrino).Pt());
	      Eta_Wass.push_back((lepton + neutrino).Eta());
	    }
	  if(nus.size() == 2)
	    {
	      bool swap_solutions = M_Wass.at(0) > M_Wass.at(1);
	      if (swap_solutions) // choose the solution with the minimal reconstructed W mass               <---------------- very important note
		{
		  swap(M_Wass.at(0), M_Wass.at(1));
		  swap(Pt_Wass.at(0), Pt_Wass.at(1));
		  swap(Eta_Wass.at(0), Eta_Wass.at(1));
		  swap(M_LepNuB.at(0), M_LepNuB.at(1));
		  swap(Pt_LepNuB.at(0), Pt_LepNuB.at(1));
		  swap(Eta_LepNuB.at(0), Eta_LepNuB.at(1));
		}
	    }
	  // Wass
	  h_wass_m->Fill(M_Wass.at(0), weight);
	  h_wass_pt->Fill(Pt_Wass.at(0), weight);
	  h_wass_eta->Fill(Eta_Wass.at(0), weight);
	  // lepton
	  h_lep_pt->Fill(lepton.pt(), weight);
	  h_lep_eta->Fill(lepton.eta(), weight);
	  h_lep_pt_vs_eta->Fill(lepton.pt(), lepton.eta(), weight);
	  // lepton/top-jet
	  h_deltaR_lep_top->Fill(deltaR(lepton, topjet.v4()), weight);
	  h_deltaPhi_lep_top->Fill(deltaPhi(lepton, topjet.v4()), weight);
	  // lepton/b-jet (b-jet with largest R-distance from top-jet)                                                        <----------- very important note
	  double deltaR_lep_bjet = deltaR(lepton, bjets.at(0).v4());
	  double deltaPhi_lep_bjet = deltaPhi(lepton, bjets.at(0).v4());
	  //double deltaR_top_bjet = deltaR(topjet.v4(), bjets.at(0).v4());
	  double deltaPhi_top_bjet = deltaPhi(topjet.v4(), bjets.at(0).v4());
	  h_deltaR_lep_bjet->Fill(deltaR_lep_bjet, weight);
	  h_deltaPhi_lep_bjet->Fill(deltaPhi_lep_bjet, weight);
	  h_deltaPhis_topB_vs_lepB->Fill(deltaPhi_lep_bjet, deltaPhi_top_bjet, weight);
	  //h_deltaRmin_bjet_leptop->Fill(min(deltaR_lep_bjet, deltaR_top_bjet), weight);
	  //h_deltaPhimin_bjet_leptop->Fill(min(deltaPhi_lep_bjet, deltaPhi_top_bjet), weight);
	  //h_deltaRmax_bjet_leptop->Fill(max(deltaR_lep_bjet, deltaR_top_bjet), weight);
	  //h_deltaPhimax_bjet_leptop->Fill(max(deltaPhi_lep_bjet, deltaPhi_top_bjet), weight);
	  // LepNuB reconstruction & Pt Balance
	  double M_combined = M_LepNuB.at(0);
	  double Pt_combined = Pt_LepNuB.at(0);
	  double PtBal_combined = (Pt_Wass.at(0) - topjet.v4().Pt())/(topjet.v4().Pt());
	  h_LepNuBSystem_M_combined->Fill(M_combined, weight);
	  h_LepNuBSystem_Pt_combined->Fill(Pt_combined, weight);
	  h_ptBalance_combined->Fill(PtBal_combined, weight);

	  // CUTS
	  bool survives_Mcut = M_combined > 280; // just sense in case of N_bjets == 2 // not really.
	  bool survives_DeltaRcut = deltaR_lep_bjet > 1.6; // this cut makes physically just sense in the case of N_bjets == 2
	  bool survives_PtBcut = PtBal_combined > -0.4;
	  // cut combinatorics up to second order
	  /*if (survives_Mcut)
	    {
	      survivesMcut_top_pt->Fill(topjet.v4().pt(), weight);
	      survivesMcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
	      survivesMcut_deltaR->Fill(deltaR_lep_bjet, weight);
	      survivesMcut_PtB->Fill(PtBal_combined, weight);
	      if (survives_DeltaRcut)
		{
		  survivesMcut_survivesDeltaRcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesMcut_survivesDeltaRcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesMcut_survivesDeltaRcut_PtB->Fill(PtBal_combined, weight);
		}
	      if (survives_PtBcut)
		{
		  survivesMcut_survivesPtBcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesMcut_survivesPtBcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesMcut_survivesPtBcut_deltaR->Fill(deltaR_lep_bjet, weight);*/
		  /*if (deltaR_lep_bjet > 1.4)
		    {
		      AllCuts->Fill(topjet.v4().pt(), weight); // just a quick histo for monday slides
		      }*/
	  /*	}
	    }
	  if (survives_DeltaRcut)
	    {
	      survivesDeltaRcut_top_pt->Fill(topjet.v4().pt(), weight);
	      survivesDeltaRcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
	      survivesDeltaRcut_Mass->Fill(M_combined, weight);
	      survivesDeltaRcut_PtB->Fill(PtBal_combined, weight);
	      if (survives_Mcut)
		{
		  survivesDeltaRcut_survivesMcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesDeltaRcut_survivesMcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesDeltaRcut_survivesMcut_PtB->Fill(PtBal_combined, weight);
		}
	      if (survives_PtBcut)
		{
		  survivesDeltaRcut_survivesPtBcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesDeltaRcut_survivesPtBcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesDeltaRcut_survivesPtBcut_Mass->Fill(M_combined, weight);
		}
	    }
	  if (survives_PtBcut)
	    {
	      survivesPtBcut_top_pt->Fill(topjet.v4().pt(), weight);
	      survivesPtBcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
	      survivesPtBcut_Mass->Fill(M_combined, weight);
	      survivesPtBcut_deltaR->Fill(deltaR_lep_bjet, weight);
	      if (survives_Mcut)
		{
		  survivesPtBcut_survivesMcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesPtBcut_survivesMcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesPtBcut_survivesMcut_deltaR->Fill(deltaR_lep_bjet, weight);
		}
	      if (survives_DeltaRcut)
		{
		  survivesPtBcut_survivesDeltaRcut_top_pt->Fill(topjet.v4().pt(), weight);
		  survivesPtBcut_survivesDeltaRcut_top_pt_rebin->Fill(topjet.v4().pt(), weight);
		  survivesPtBcut_survivesDeltaRcut_Mass->Fill(M_combined, weight);
		  }
		  }*/
	  // all CUTS
	  if (survives_Mcut && survives_DeltaRcut && survives_PtBcut)
	    {
	      AllCuts_top_pt->Fill(topjet.v4().pt(), weight);
	      AllCuts_top_pt_rebin->Fill(topjet.v4().pt(), weight);
	    }
	  /*if (M_combined > 280 && deltaR_lep_bjet > 1.6 && PtBal_combined > -0.4)
	    {
	      AllCuts_2Btags->Fill(topjet.v4().pt(), weight);
	    }
	  if (M_combined > 320 && deltaR_lep_bjet > 1.8 && PtBal_combined > -0.2)
	    {
	      AllCuts_1Btag->Fill(topjet.v4().pt(), weight);
	    }
	  if (M_combined > 320 && deltaR_lep_bjet > 1.4 && deltaR_lep_bjet < 3.2 && PtBal_combined > -0.2)
	    {
	      AllCuts_1Btag_v2->Fill(topjet.v4().pt(), weight);
	      }*/
	}
    }
  
  int n_jets = 0;
  int jet_ind = 0;
  // fill HOTVR hists
  for (TopJet topjet : hotvrJets)
    {
      ++jet_ind;
      if (m_topjetid)
	{
	  if (!(*m_topjetid)(topjet, event)) continue;
	  if (jet_ind == 1 && n_jets == 0) NLeadingTopjet->Fill(0., weight);
	  else if(n_jets == 0) NLeadingTopjet->Fill(1., weight);
	}
      ++n_jets;
      vector<Jet> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();
      double a_topjet = topjet.jetArea();
      double fpt = -1;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = subjets.at(0).v4().pt() / pt_topjet;
	}
      double mpair = -1;
      double pt_sub1 = -1;
      double a_sub1 = -1;
      double pt_sub2 = -1;
      double a_sub2 = -1;
      double pt_sub3 = -1;
      double a_sub3 = -1;

      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	  pt_sub1 = subjets.at(0).v4().pt();
	  a_sub1  = subjets.at(0).jetArea();
	  pt_sub2 = subjets.at(1).v4().pt();
	  a_sub2  = subjets.at(1).jetArea();
	  for (unsigned int i = 0; i < subjets.size(); ++i)
	    {
	      pt_sub3 = subjets.at(i).v4().pt();
	      a_sub3  = subjets.at(i).jetArea();

	    }
	  pt_sub3 = subjets.at(2).v4().pt();
	  a_sub3  = subjets.at(2).jetArea();
	}

      // TH1Fs
      Pt_HotvrTopjets->Fill(pt_topjet, weight);
      Eta_HotvrTopjets->Fill(topjet.v4().eta(), weight);
      M_HotvrTopjets->Fill(topjet.v4().M(), weight);
      A_HotvrTopjets->Fill(a_topjet, weight);
      NSub_HotvrTopjets->Fill(subjets.size(), weight);
      Fpt_HotvrTopjets->Fill(fpt, weight);
      Mpair_HotvrTopjets->Fill(mpair, weight);
      Tau32_HotvrTopjets->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);
      Pt_HotvrTopjets_Sub1->Fill(pt_sub1, weight);
      A_HotvrTopjets_Sub1->Fill(a_sub1, weight);
      Pt_HotvrTopjets_Sub2->Fill(pt_sub2, weight);
      A_HotvrTopjets_Sub2->Fill(a_sub2, weight);
      Pt_HotvrTopjets_Sub3->Fill(pt_sub3, weight);
      A_HotvrTopjets_Sub3->Fill(a_sub3, weight);

      if (muons.size() > 0)
	{
	  DeltaR_Lepton_HotvrTopjets->Fill(deltaR(topjet.v4(), muons.at(0).v4()), weight);
	  DeltaPhi_Lepton_HotvrTopjets->Fill(deltaPhi(topjet.v4(), muons.at(0).v4()), weight);
	}
      else if (electrons.size() > 0)
	{
	  DeltaR_Lepton_HotvrTopjets->Fill(deltaR(topjet.v4(), electrons.at(0).v4()), weight);
	  DeltaPhi_Lepton_HotvrTopjets->Fill(deltaPhi(topjet.v4(), electrons.at(0).v4()), weight);
	}

      if (pt_topjet > 1550.) 
	{
	  Pt_rebin_HotvrTopjets->Fill(1550., weight);
	  Pt_vs_Eta_HotvrTopjets->Fill(1550., topjet.v4().eta(), weight);
	  Pt_vs_Eta_HotvrRebin->Fill(1550., topjet.v4().eta(), weight);
	}
      else
	{
	  Pt_rebin_HotvrTopjets->Fill(pt_topjet, weight);
	  Pt_vs_Eta_HotvrTopjets->Fill(pt_topjet, topjet.v4().eta(), weight);
	  Pt_vs_Eta_HotvrRebin->Fill(pt_topjet, topjet.v4().eta(), weight);
	}
      EtaAbs_HotvrTopjets->Fill(topjet.v4().eta(), weight);

      if (n_jets == 1)
	{
	  if (pt_topjet > 1550.) 
	    {
	      Pt_HotvrTopjet1->Fill(1550., weight);
	    }
	  else  
	    {
	      Pt_HotvrTopjet1->Fill(pt_topjet, weight);
	    }
	  Eta_HotvrTopjet1->Fill(topjet.v4().eta(), weight);
	  M_HotvrTopjet1->Fill(topjet.v4().M(), weight);
	  A_HotvrTopjet1->Fill(a_topjet, weight);
	  NSub_HotvrTopjet1->Fill(subjets.size(), weight);
	  Fpt_HotvrTopjet1->Fill(fpt, weight);
	  Mpair_HotvrTopjet1->Fill(mpair, weight);
	  Tau32_HotvrTopjet1->Fill(topjet.tau3_groomed()/topjet.tau2_groomed(), weight);
	  Pt_HotvrTopjet1_Sub1->Fill(pt_sub1, weight);
	  A_HotvrTopjet1_Sub1->Fill(a_sub1, weight);
	  Pt_HotvrTopjet1_Sub2->Fill(pt_sub2, weight);
	  A_HotvrTopjet1_Sub2->Fill(a_sub2, weight);
	  Pt_HotvrTopjet1_Sub3->Fill(pt_sub3, weight);
	  A_HotvrTopjet1_Sub3->Fill(a_sub3, weight);
	  if (muons.size() > 0)
	    {
	      DeltaR_Lepton_HotvrTopjet1->Fill(deltaR(topjet.v4(), muons.at(0).v4()), weight);
	      DeltaPhi_Lepton_HotvrTopjet1->Fill(deltaPhi(topjet.v4(), muons.at(0).v4()), weight);
	    }
	  else if (electrons.size() > 0)
	    {
	      DeltaR_Lepton_HotvrTopjet1->Fill(deltaR(topjet.v4(), electrons.at(0).v4()), weight);
	      DeltaPhi_Lepton_HotvrTopjet1->Fill(deltaPhi(topjet.v4(), electrons.at(0).v4()), weight);
	    }
	}
      
    }
  N_HotvrTopjets->Fill(n_jets, weight);

}

HOTVRHists::~HOTVRHists(){}

/*
 * WARNING: Fill this Hists only after cuts are applied to ensure
 * there is >=1 HOTVR TopJet.
 * 
 * This Hists Class implements Histograms with informations about the
 * physics performance of the HOTVR algorithm.
 * 
 */
HOTVRGenHists::HOTVRGenHists(Context & ctx, const string & dirname, const boost::optional<GenTopJetId> &gentopjetid): 
  Hists(ctx, dirname), 
  h_boostedsingletopgen(ctx.get_handle<BoostedSingleTopGen>("BoostedSingleTopgen")),
  m_gentopjetid(gentopjetid)
{
  // book all histograms here

  // HOTVR hists
  N_HotvrTopjets          = book<TH1F>("N_HOTVR",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjets         = book<TH1F>("Pt_HOTVR",         "p_{t}^{topjet} [GeV]", 40, 0, 1600);
  Eta_HotvrTopjets        = book<TH1F>("Eta_HOTVR",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjets          = book<TH1F>("M_HOTVR",          "M^{topjet} [GeV]", 40,  0, 400);
  R_HotvrTopjets          = book<TH1F>("R_HOTVR",          "R_{topjet}", 15,  0.1, 1.6);
  NSub_HotvrTopjets       = book<TH1F>("NSub_HOTVR",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjets        = book<TH1F>("Fpt_HOTVR",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjets      = book<TH1F>("Mpair_HOTVR",      "M_pair [GeV]", 40,  0, 200);
  Pt_HotvrTopjets_Sub1    = book<TH1F>("Pt_HOTVR_Subjet1",         "p_{t}^{subjet1} [GeV]", 100, 0, 1000);
  Pt_HotvrTopjets_Sub2    = book<TH1F>("Pt_HOTVR_Subjet2",         "p_{t}^{subjet2} [GeV]", 100, 0, 1000);
  Pt_HotvrTopjets_Sub3    = book<TH1F>("Pt_HOTVR_Subjet3",         "p_{t}^{subjet3} [GeV]", 100, 0, 1000);

  DeltaR_L_HotvrTopjets   = book<TH1F>("DeltaR_L_HOTVR",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjets = book<TH1F>("DeltaPhi_L_HOTVR", "#Delta #phi_{l,t}", 20, 0, 4);

  N_HotvrTopjet1          = book<TH1F>("N_HOTVR1",          "N_{topjets}", 10,  0, 10);
  Pt_HotvrTopjet1         = book<TH1F>("Pt_HOTVR1",         "p_{t}^{topjet} [GeV]", 40, 0, 1600);
  Eta_HotvrTopjet1        = book<TH1F>("Eta_HOTVR1",        "#eta^{topjet}", 30, -6, 6);
  M_HotvrTopjet1          = book<TH1F>("M_HOTVR1",          "M^{topjet} [GeV]", 40,  0, 400);
  R_HotvrTopjet1          = book<TH1F>("R_HOTVR1",           "R_{topjet}", 15,  0.1, 1.6);
  NSub_HotvrTopjet1       = book<TH1F>("NSub_HOTVR1",       "N_{subjets}", 10,  0, 10);
  Fpt_HotvrTopjet1        = book<TH1F>("Fpt_HOTVR1",        "f_{pt, 1}", 20,  0, 1);
  Mpair_HotvrTopjet1      = book<TH1F>("Mpair_HOTVR1",      "M_pair [GeV]", 40,  0, 200);
  Pt_HotvrTopjet1_Sub1    = book<TH1F>("Pt_HOTVR1_Subjet1",         "p_{t}^{subjet1} [GeV]", 100, 0, 1000);
  Pt_HotvrTopjet1_Sub2    = book<TH1F>("Pt_HOTVR1_Subjet2",         "p_{t}^{subjet2} [GeV]", 100, 0, 1000);
  Pt_HotvrTopjet1_Sub3    = book<TH1F>("Pt_HOTVR1_Subjet3",         "p_{t}^{subjet3} [GeV]", 100, 0, 1000);

  DeltaR_L_HotvrTopjet1   = book<TH1F>("DeltaR_L_HOTVR1",   "#Delta R_{l,t}", 20, 0, 4);
  DeltaPhi_L_HotvrTopjet1 = book<TH1F>("DeltaPhi_L_HOTVR1", "#Delta #phi_{l,t}", 20, 0, 4);

}

void HOTVRGenHists::fill(const Event & event) {  

  double weight = event.weight; // event weight
  vector<GenTopJet> hotvrJets = *event.gentopjets;
  BoostedSingleTopGen gen = event.get(h_boostedsingletopgen);
  LorentzVector muon = gen.ChargedLepton();
  
  
  int n_jets = 0;

  // fill HOTVR hists
  for (GenTopJet topjet : hotvrJets)
    {
      if (m_gentopjetid)
	{
	  if (!(*m_gentopjetid)(topjet, event)) continue;
	}
      ++n_jets;
      vector<Particle> subjets = topjet.subjets();

      double pt_topjet = topjet.v4().pt();
      double R_topjet  = 600/pt_topjet;
	if (R_topjet < 0.1) R_topjet = 0.1;
	else if (R_topjet > 1.5) R_topjet = 1.5;
      double fpt = -1;
      // fpt can only be calculated if there are subjets
      if (subjets.size() >= 1)
	{
	  fpt = subjets.at(0).v4().pt() / pt_topjet;
	}
      double mpair = -1;
      double pt_sub1 = -1;
      double pt_sub2 = -1;
      double pt_sub3 = -1;

      // mpair can only be calculated if there are at least 3 subjets
      if (subjets.size() >= 3)
	{
	  double m12 = (subjets.at(0).v4() + subjets.at(1).v4()).M();
	  double m13 = (subjets.at(0).v4() + subjets.at(2).v4()).M();
	  double m23 = (subjets.at(1).v4() + subjets.at(2).v4()).M();
	  mpair = min(min(m12, m13), m23);
	  pt_sub1 = subjets.at(0).v4().pt();
	  pt_sub2 = subjets.at(1).v4().pt();
	  pt_sub3 = subjets.at(2).v4().pt();
	}

      // TH1Fs
      Pt_HotvrTopjets->Fill(pt_topjet, weight);
      Eta_HotvrTopjets->Fill(topjet.v4().eta(), weight);
      M_HotvrTopjets->Fill(topjet.v4().M(), weight);
      R_HotvrTopjets->Fill(R_topjet, weight);
      NSub_HotvrTopjets->Fill(subjets.size(), weight);
      Fpt_HotvrTopjets->Fill(fpt, weight);
      Mpair_HotvrTopjets->Fill(mpair, weight);
      Pt_HotvrTopjets_Sub1->Fill(pt_sub1, weight);
      Pt_HotvrTopjets_Sub2->Fill(pt_sub2, weight);
      Pt_HotvrTopjets_Sub3->Fill(pt_sub3, weight);

      DeltaR_L_HotvrTopjets->Fill(deltaR(topjet.v4(), muon), weight);
      DeltaPhi_L_HotvrTopjets->Fill(deltaPhi(topjet.v4(), muon), weight);

      if (n_jets == 1)
	{
	  Pt_HotvrTopjet1->Fill(pt_topjet, weight);
	  Eta_HotvrTopjet1->Fill(topjet.v4().eta(), weight);
	  M_HotvrTopjet1->Fill(topjet.v4().M(), weight);
	  R_HotvrTopjet1->Fill(R_topjet, weight);
	  NSub_HotvrTopjet1->Fill(subjets.size(), weight);
	  Fpt_HotvrTopjet1->Fill(fpt, weight);
	  Mpair_HotvrTopjet1->Fill(mpair, weight);
	  Pt_HotvrTopjet1_Sub1->Fill(pt_sub1, weight);
	  Pt_HotvrTopjet1_Sub2->Fill(pt_sub2, weight);
	  Pt_HotvrTopjet1_Sub3->Fill(pt_sub3, weight);

	  DeltaR_L_HotvrTopjet1->Fill(deltaR(topjet.v4(), muon), weight);
	  DeltaPhi_L_HotvrTopjet1->Fill(deltaPhi(topjet.v4(), muon), weight);
	  //   }
	}
      
    }
  N_HotvrTopjets->Fill(n_jets, weight);


}

HOTVRGenHists::~HOTVRGenHists(){}


HOTVRPileUpHists::HOTVRPileUpHists(Context & ctx, const string & dirname):
  Hists(ctx, dirname){

  n = book<TH1F>("number", "N_{PV}", 10, 0, 50);
  u = book<TH1F>("uncorrected", "N_{PV}", 10, 0, 50);
  c = book<TH1F>("corrected", "N_{PV}", 10, 0, 50);

}
void HOTVRPileUpHists::fill(const Event &event) {

  double weight = event.weight;
  double rho = event.rho;
  bool u_flag = false;
  bool c_flag = false;
  vector<TopJet> topjets = *event.topjets;
  n->Fill(event.pvs->size(), weight);
  for (TopJet topjet : topjets)
    {
      // int i = 1;
      // cout << i << ": " << topjet.pt() << endl;
      LorentzVector temp_jet;
      
      if (topjet.pt() > 200. && abs(topjet.eta()) < 2.5)
	{
	  u_flag = true;
	}
      for (Jet subjet : topjet.subjets())
      	{
      	  double a = 1 - ((rho * subjet.jetArea()) / subjet.pt());
      	  temp_jet += subjet.v4() * a;
      	}
      if (temp_jet.Pt() > 200. && abs(temp_jet.Eta()) < 2.5)
      	{
      	  c_flag = true;
      	}
    }

  if (u_flag) u->Fill(event.pvs->size(), weight);
  if (c_flag) c->Fill(event.pvs->size(), weight);
}

HOTVRPileUpHists::~HOTVRPileUpHists(){}
