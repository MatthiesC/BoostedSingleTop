#pragma once

#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Hists.h"

#include "TH1F.h"
#include "TH2F.h"

namespace uhh2 {

  class MVAHists: public uhh2::Hists {
  public:
    MVAHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & e) override;
    void fill_(const uhh2::Event & e, const double & mvaD, const int & n_btags);
    virtual ~MVAHists();

  protected:

    uhh2::Event::Handle<FlavorParticle> h_primlep;

    // Declare all histograms here ...
    TH1F *H_mvaD;
    TH1F *H_mvaD_20bins;
    TH1F *H_mvaD_rebin;
    TH1F *H_mvaD_1b;
    TH1F *H_mvaD_2b;

    TH1F* H_nocuts_top_pt_1000;
    TH1F* H_nocuts_top_pt_rebin;
    TH1F* H_nocuts_top_eta_1000;
    TH1F* H_nocuts_top_eta_rebin;
    TH2F* H_nocuts_top_pt_vs_eta_1000;
    TH2F* H_nocuts_top_pt_vs_eta_rebin;

    TH1F* H_nocuts_lept_pt;
    TH1F* H_nocuts_lept_eta;

    TH1F* H_BDTcut_top_pt_1000;
    TH1F* H_BDTcut_top_pt_rebin;
    TH1F* H_BDTcut_top_eta_1000;
    TH1F* H_BDTcut_top_eta_rebin;
    TH2F* H_BDTcut_top_pt_vs_eta_1000;
    TH2F* H_BDTcut_top_pt_vs_eta_rebin;

    TH1F* H_BDTcut_lept_pt;
    TH1F* H_BDTcut_lept_eta;

    TH1F* H_200to300_BDT;
    TH1F* H_200to300_BDT_rebin;
    TH1F* H_200to300_top_eta;
    TH1F* H_200to300_lept_pt;
    TH1F* H_200to300_lept_eta;

    TH1F* H_300to400_BDT;
    TH1F* H_300to400_BDT_rebin;
    TH1F* H_300to400_top_eta;
    TH1F* H_300to400_lept_pt;
    TH1F* H_300to400_lept_eta;

    TH1F* H_400to600_BDT;
    TH1F* H_400to600_BDT_rebin;
    TH1F* H_400to600_top_eta;
    TH1F* H_400to600_lept_pt;
    TH1F* H_400to600_lept_eta;

    TH1F* H_600to1200_BDT;
    TH1F* H_600to1200_BDT_rebin;
    TH1F* H_600to1200_top_eta;
    TH1F* H_600to1200_lept_pt;
    TH1F* H_600to1200_lept_eta;

  };
   
}
