#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/core/include/Muon.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/common/include/JetIds.h"

#include "UHH2/BoostedSingleTop/include/HOTVRIds.h"
#include "UHH2/BoostedSingleTop/include/HOTVRGenIds.h"

#include <vector>

namespace uhh2 {

  /**  \brief Example class for booking and filling histograms
   * 
   * NOTE: This class uses the 'hist' method to retrieve histograms.
   * This requires a string lookup and is therefore slow if you have
   * many histograms. Therefore, it is recommended to use histogram
   * pointers as member data instead, like in 'common/include/ElectronHists.h'.
   */
  class HOTVRHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<TopJetId> &topjetid = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRHists();
  protected:
    // HOTVR
    TH1F *N_HotvrTopjets, *Pt_HotvrTopjets, *Eta_HotvrTopjets, *M_HotvrTopjets, *A_HotvrTopjets;
    TH1F *NSub_HotvrTopjets, *Fpt_HotvrTopjets, *Mpair_HotvrTopjets, *Tau32_HotvrTopjets;
    TH1F *DeltaR_Lepton_HotvrTopjets, *DeltaPhi_Lepton_HotvrTopjets;
    TH1F *Pt_HotvrTopjets_Sub1, *A_HotvrTopjets_Sub1, *Pt_HotvrTopjets_Sub2, *A_HotvrTopjets_Sub2, *Pt_HotvrTopjets_Sub3, *A_HotvrTopjets_Sub3;
    TH1F *N_HotvrTopjet1, *Pt_HotvrTopjet1,  *Eta_HotvrTopjet1, *M_HotvrTopjet1, *A_HotvrTopjet1;
    TH1F *NSub_HotvrTopjet1, *Fpt_HotvrTopjet1, *Mpair_HotvrTopjet1, *Tau32_HotvrTopjet1;
    TH1F *DeltaR_Lepton_HotvrTopjet1, *DeltaPhi_Lepton_HotvrTopjet1;
    TH1F *Pt_HotvrTopjet1_Sub1, *A_HotvrTopjet1_Sub1, *Pt_HotvrTopjet1_Sub2, *A_HotvrTopjet1_Sub2, *Pt_HotvrTopjet1_Sub3, *A_HotvrTopjet1_Sub3;
    TH2F *Pt_vs_Eta_HotvrTopjets, *Pt_vs_Eta_HotvrTopjet1, *Pt_vs_Eta_HotvrRebin;
    TH1F *Pt_rebin_HotvrTopjets, *EtaAbs_HotvrTopjets;
    TH1F *NLeadingTopjet;
    // b-jets
    TH1F *A_ak4, *Pt_ak4;
    TH1F *N_bjets_loose, *N_bjets_medium, *N_bjets_tight;
    TH1F *loose_bjets_pt, *loose_bjets_eta, *loose_bjets_m;
    TH2F *loose_bjets_pt_vs_eta;


    // all loose b-jets
    TH1F *bjets_pt;
    TH1F *bjets_eta;
    TH2F *bjets_pt_vs_eta;
    TH1F *bjets_m;
    TH1F *bjets_a;

    // if N_bjets == 1
    TH1F *bjet1_pt;
    TH1F *bjet1_eta;
    TH2F *bjet1_pt_vs_eta;
    TH1F *bjet1_m;
    TH1F *bjet1_a;
    TH1F *bjet1_deltaR_top;
    TH1F *bjet1_deltaPhi_top;
    TH1F *bjet1_deltaR_lep;
    TH1F *bjet1_deltaPhi_lep;

    // if N_bjets == 2
    TH1F *bjet2_1_pt;
    TH1F *bjet2_2_pt;
    TH1F *bjet2_1_eta;
    TH1F *bjet2_2_eta;
    TH2F *bjet2_1_pt_vs_eta;
    TH2F *bjet2_2_pt_vs_eta;
    TH1F *bjet2_1_m;
    TH1F *bjet2_2_m;
    TH1F *bjet2_1_a;
    TH1F *bjet2_2_a;

    // if N_bjets == 3
    TH1F *bjet3_1_pt;
    TH1F *bjet3_2_pt;
    TH1F *bjet3_3_pt;
    TH1F *bjet3_1_eta;
    TH1F *bjet3_2_eta;
    TH1F *bjet3_3_eta;
    TH2F *bjet3_1_pt_vs_eta;
    TH2F *bjet3_2_pt_vs_eta;
    TH2F *bjet3_3_pt_vs_eta;
    TH1F *bjet3_1_m;
    TH1F *bjet3_2_m;
    TH1F *bjet3_3_m;
    TH1F *bjet3_1_a;
    TH1F *bjet3_2_a;
    TH1F *bjet3_3_a;

    // if (N_bjets == 1 && N_topjets == 1)
    //   {
    // top-jet
    TH1F *h_top_pt;
    TH1F *h_top_eta;
    TH2F *h_top_pt_vs_eta;
    TH1F *h_top_m;
    TH1F *h_top_a;
    TH1F *h_top_nsub;
    TH1F *h_top_fpt;
    TH1F *h_top_mpair;
    TH1F *h_top_tau32;

    // subjets (just filled if n_subjets >= 3)
    TH1F *h_sub1_pt;
    TH1F *h_sub1_a;
    TH1F *h_sub2_pt;
    TH1F *h_sub2_a;
    TH1F *h_sub3_pt;
    TH1F *h_sub3_a;

    // top-jet/b-jet
    TH1F *h_deltaR_top_bjet;
    TH1F *h_deltaPhi_top_bjet;

    // Wass
    TH1F *h_wass_m;
    TH1F *h_wass_pt;
    TH1F *h_wass_eta;
    // lepton
    TH1F *h_lep_pt;
    TH1F *h_lep_eta;
    TH2F *h_lep_pt_vs_eta;
    // lepton/top-jet
    TH1F *h_deltaR_lep_top;
    TH1F *h_deltaPhi_lep_top;
    // lepton/b-jet
    TH1F *h_deltaR_lep_bjet;
    TH1F *h_deltaPhi_lep_bjet;
    TH2F *h_deltaPhis_topB_vs_lepB;
    // Lepton+Nu+BJet Reconstruction
    TH1F *h_LepNuBSystem_M_combined;
    TH1F *h_LepNuBSystem_Pt_combined;
    // Pt Balance
    TH1F *h_ptBalance_combined;
    // Cuts
    /* TH1F *h_MCut_top_pt, *h_MCut_top_pt_rebin;
    TH1F *h_DeltaRCut_top_pt, *h_DeltaRCut_top_pt_rebin;
    TH1F *h_BothCuts_top_pt, *h_BothCuts_top_pt_rebin;
    TH1F *h_survivesMcut_deltaR_lep_bjet, *h_survivesMcut_deltaPhi_lep_bjet, *h_survivesDeltaRcut_M_combined;*/
    /*
    TH1F *survivesMcut_top_pt;
    TH1F *survivesMcut_top_pt_rebin;
    TH1F *survivesMcut_deltaR;
    TH1F *survivesMcut_PtB;

    TH1F *survivesMcut_survivesDeltaRcut_top_pt;
    TH1F *survivesMcut_survivesDeltaRcut_top_pt_rebin;
    TH1F *survivesMcut_survivesDeltaRcut_PtB;

    TH1F *survivesMcut_survivesPtBcut_top_pt;
    TH1F *survivesMcut_survivesPtBcut_top_pt_rebin;
    TH1F *survivesMcut_survivesPtBcut_deltaR;

    TH1F *survivesDeltaRcut_top_pt;
    TH1F *survivesDeltaRcut_top_pt_rebin;
    TH1F *survivesDeltaRcut_Mass;
    TH1F *survivesDeltaRcut_PtB;

    TH1F *survivesDeltaRcut_survivesMcut_top_pt;
    TH1F *survivesDeltaRcut_survivesMcut_top_pt_rebin;
    TH1F *survivesDeltaRcut_survivesMcut_PtB;

    TH1F *survivesDeltaRcut_survivesPtBcut_top_pt;
    TH1F *survivesDeltaRcut_survivesPtBcut_top_pt_rebin;
    TH1F *survivesDeltaRcut_survivesPtBcut_Mass;

    TH1F *survivesPtBcut_top_pt;
    TH1F *survivesPtBcut_top_pt_rebin;
    TH1F *survivesPtBcut_Mass;
    TH1F *survivesPtBcut_deltaR;

    TH1F *survivesPtBcut_survivesMcut_top_pt;
    TH1F *survivesPtBcut_survivesMcut_top_pt_rebin;
    TH1F *survivesPtBcut_survivesMcut_deltaR;

    TH1F *survivesPtBcut_survivesDeltaRcut_top_pt;
    TH1F *survivesPtBcut_survivesDeltaRcut_top_pt_rebin;
    TH1F *survivesPtBcut_survivesDeltaRcut_Mass;
    */
    // all CUTS
    TH1F *AllCuts_top_pt;
    TH1F *AllCuts_top_pt_rebin;
    //TH1F *AllCuts_2Btags, *AllCuts_1Btag, *AllCuts_1Btag_v2;

    //   }

    boost::optional<TopJetId> m_topjetid;
  };

  class HOTVRGenHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRGenHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<GenTopJetId> &gentopjetid = boost::none);

    virtual void fill(const uhh2::Event & event) override;
    virtual ~HOTVRGenHists();
  protected:
    // HOTVR
    TH1F *N_HotvrTopjets, *Pt_HotvrTopjets, *Eta_HotvrTopjets, *M_HotvrTopjets, *R_HotvrTopjets;
    TH1F *NSub_HotvrTopjets, *Fpt_HotvrTopjets, *Mpair_HotvrTopjets, *Tau32_HotvrTopjets;
    TH1F *DeltaR_L_HotvrTopjets, *DeltaPhi_L_HotvrTopjets;
    TH1F *Pt_HotvrTopjets_Sub1, *Pt_HotvrTopjets_Sub2, *Pt_HotvrTopjets_Sub3;
    TH1F *N_HotvrTopjet1, *Pt_HotvrTopjet1, *Eta_HotvrTopjet1, *M_HotvrTopjet1, *R_HotvrTopjet1;
    TH1F *NSub_HotvrTopjet1, *Fpt_HotvrTopjet1, *Mpair_HotvrTopjet1, *Tau32_HotvrTopjet1;
    TH1F *DeltaR_L_HotvrTopjet1, *DeltaPhi_L_HotvrTopjet1;
    TH1F *Pt_HotvrTopjet1_Sub1, *Pt_HotvrTopjet1_Sub2, *Pt_HotvrTopjet1_Sub3;
    // b-jets
    TH1F *N_bjets_loose, *N_bjets_medium, *N_bjets_tight;
    uhh2::Event::Handle<BoostedSingleTopGen> h_boostedsingletopgen;
    boost::optional<GenTopJetId> m_gentopjetid;
  };


  class HOTVRPileUpHists: public uhh2::Hists {
  public:
    // use the same constructor arguments as Hists for forwarding:
    HOTVRPileUpHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;
    virtual ~HOTVRPileUpHists();
  protected:

    TH1F *n,*u,*c;

  };
   
}
