#pragma once

#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/TopJet.h"
#include "UHH2/core/include/Electron.h"
#include "UHH2/core/include/Muon.h"
#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/common/include/JetIds.h"

#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWch.h"
#include "UHH2/BoostedSingleTop/include/TTbarGen_.h"

#include "TH1F.h"
#include "TH2F.h"

#include <vector>

namespace uhh2 {

  class MatchHists: public uhh2::Hists {
  public:
    MatchHists(uhh2::Context & ctx, const std::string & dirname, const boost::optional<TopJetId> & topjetid = boost::none);

    virtual void fill(const uhh2::Event & e) override;
    virtual ~MatchHists();

    bool is_tw, is_tt;

  protected:
    uhh2::Event::Handle<SingleTopGen_tWch> h_singletopgen_twch;
    uhh2::Event::Handle<TTbarGen_> h_ttbargen;
    boost::optional<TopJetId> m_topjetid;

    // Declare all histograms here ...
    TH1F *H_shat, *H_DecayChannel;

    TH1F *H_topgen_pt, *H_lepgen_pt, *H_topgen_eta, *H_lepgen_eta, *H_topjet_pt, *H_leprec_pt, *H_topjet_eta, *H_leprec_eta;
    TH1F *H_match_deltaR_lep, *H_match_deltaR_top;
    TH1F *H_match_deltaPhi_lep, *H_match_deltaPhi_top, *H_match_deltaPhi_neu;

  };
   
}
