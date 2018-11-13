#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"

#include <vector>

#include "TH1F.h"
#include "TH2F.h"

namespace uhh2 {

  class MVAHists: public uhh2::Hists {
  public:
    MVAHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & e) override;
    void fill_(const uhh2::Event & e, const double & mvaD, const std::vector<double> & mva_inputvars, const double & BDTcut);
    virtual ~MVAHists();

  protected:

    uhh2::Event::Handle<FlavorParticle> h_primlep;
    uhh2::Event::Handle<TopJet> h_primtopjet;

    // Declare all histograms here ...
    TH1F *H_mvaD;
    TH1F *H_mvaD_20bins;
    TH1F *H_mvaD_rebin;
    TH1F *H_mvaD_rebin01;
    TH1F *H_mvaD_1b;
    TH1F *H_mvaD_2b;


    TH1F* H_nocuts_top_pt_1000;
    TH1F* H_nocuts_top_pt_rebin;
    TH1F* H_nocuts_top_eta_1000;
    TH1F* H_nocuts_top_eta_rebin;
    TH2F* H_nocuts_top_pt_vs_eta_1000;
    TH2F* H_nocuts_top_pt_vs_eta_rebin;
    TH1F* H_nocuts_top_mass;

    TH1F* H_nocuts_tmva_n_btags;
    TH1F* H_nocuts_tmva_deltaPhi_leptonbottom;
    TH1F* H_nocuts_tmva_deltaPhi_topbottom;
    TH1F* H_nocuts_tmva_mass_LepNuB;
    TH1F* H_nocuts_tmva_pt_balance;
    TH1F* H_nocuts_tmva_Wass_pt;
    TH1F* H_nocuts_tmva_lepton_pt;
    TH1F* H_nocuts_tmva_lepton_eta;
    TH1F* H_nocuts_tmva_deltaPhi_leptonnextjet;
    TH1F* H_nocuts_tmva_met_pt;


    TH1F* H_BDTcut_top_pt_1000;
    TH1F* H_BDTcut_top_pt_rebin;
    TH1F* H_BDTcut_top_eta_1000;
    TH1F* H_BDTcut_top_eta_rebin;
    TH2F* H_BDTcut_top_pt_vs_eta_1000;
    TH2F* H_BDTcut_top_pt_vs_eta_rebin;
    TH1F* H_BDTcut_top_mass;

    TH1F* H_BDTcut_tmva_n_btags;
    TH1F* H_BDTcut_tmva_deltaPhi_leptonbottom;
    TH1F* H_BDTcut_tmva_deltaPhi_topbottom;
    TH1F* H_BDTcut_tmva_mass_LepNuB;
    TH1F* H_BDTcut_tmva_pt_balance;
    TH1F* H_BDTcut_tmva_Wass_pt;
    TH1F* H_BDTcut_tmva_lepton_pt;
    TH1F* H_BDTcut_tmva_lepton_eta;
    TH1F* H_BDTcut_tmva_deltaPhi_leptonnextjet;
    TH1F* H_BDTcut_tmva_met_pt;


    TH1F* H_200to300_BDT;
    TH1F* H_200to300_BDT_20bins;
    TH1F* H_200to300_BDT_rebin;
    TH1F* H_200to300_BDT_rebin01;
    TH1F* H_200to300_top_eta;

    TH1F* H_200to300_tmva_n_btags;
    TH1F* H_200to300_tmva_deltaPhi_leptonbottom;
    TH1F* H_200to300_tmva_deltaPhi_topbottom;
    TH1F* H_200to300_tmva_mass_LepNuB;
    TH1F* H_200to300_tmva_pt_balance;
    TH1F* H_200to300_tmva_Wass_pt;
    TH1F* H_200to300_tmva_lepton_pt;
    TH1F* H_200to300_tmva_lepton_eta;
    TH1F* H_200to300_tmva_deltaPhi_leptonnextjet;
    TH1F* H_200to300_tmva_met_pt;


    TH1F* H_300to400_BDT;
    TH1F* H_300to400_BDT_20bins;
    TH1F* H_300to400_BDT_rebin;
    TH1F* H_300to400_BDT_rebin01;
    TH1F* H_300to400_top_eta;

    TH1F* H_300to400_tmva_n_btags;
    TH1F* H_300to400_tmva_deltaPhi_leptonbottom;
    TH1F* H_300to400_tmva_deltaPhi_topbottom;
    TH1F* H_300to400_tmva_mass_LepNuB;
    TH1F* H_300to400_tmva_pt_balance;
    TH1F* H_300to400_tmva_Wass_pt;
    TH1F* H_300to400_tmva_lepton_pt;
    TH1F* H_300to400_tmva_lepton_eta;
    TH1F* H_300to400_tmva_deltaPhi_leptonnextjet;
    TH1F* H_300to400_tmva_met_pt;


    TH1F* H_400to600_BDT;
    TH1F* H_400to600_BDT_20bins;
    TH1F* H_400to600_BDT_rebin;
    TH1F* H_400to600_BDT_rebin01;
    TH1F* H_400to600_top_eta;

    TH1F* H_400to600_tmva_n_btags;
    TH1F* H_400to600_tmva_deltaPhi_leptonbottom;
    TH1F* H_400to600_tmva_deltaPhi_topbottom;
    TH1F* H_400to600_tmva_mass_LepNuB;
    TH1F* H_400to600_tmva_pt_balance;
    TH1F* H_400to600_tmva_Wass_pt;
    TH1F* H_400to600_tmva_lepton_pt;
    TH1F* H_400to600_tmva_lepton_eta;
    TH1F* H_400to600_tmva_deltaPhi_leptonnextjet;
    TH1F* H_400to600_tmva_met_pt;


    TH1F* H_600to1200_BDT;
    TH1F* H_600to1200_BDT_20bins;
    TH1F* H_600to1200_BDT_rebin;
    TH1F* H_600to1200_BDT_rebin01;
    TH1F* H_600to1200_top_eta;

    TH1F* H_600to1200_tmva_n_btags;
    TH1F* H_600to1200_tmva_deltaPhi_leptonbottom;
    TH1F* H_600to1200_tmva_deltaPhi_topbottom;
    TH1F* H_600to1200_tmva_mass_LepNuB;
    TH1F* H_600to1200_tmva_pt_balance;
    TH1F* H_600to1200_tmva_Wass_pt;
    TH1F* H_600to1200_tmva_lepton_pt;
    TH1F* H_600to1200_tmva_lepton_eta;
    TH1F* H_600to1200_tmva_deltaPhi_leptonnextjet;
    TH1F* H_600to1200_tmva_met_pt;

  };
   
}
