#pragma once

#include "UHH2/core/include/fwd.h"
#include "UHH2/core/include/Selection.h"

#include "UHH2/common/include/TopJetIds.h"

#include "UHH2/BoostedSingleTop/include/BoostedSingleTopGen.h"
#include "UHH2/BoostedSingleTop/include/BoostedSingleTopHypothesis.h"

namespace uhh2 {

  /** 
   * This class selects event by requiring there to be at least n_min
   * and maximal n_max topjets, with pt > pt_min and eta < eta_max
   */
  class NHotvrSelection: public uhh2::Selection {
  public:
    NHotvrSelection(unsigned int n_min_, unsigned int n_max_, double pt_min_, double eta_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    unsigned int n_min, n_max;
    double pt_min, eta_max;
  };

  /**
   * This class selects events by requiring MET > met_min
   */
  class METSelection: public uhh2::Selection {
  public:
    METSelection(double met_min_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double met_min;
  };

  /**
   * This class selects events by requiring ST > ST_min
   */
  class STSelection: public uhh2::Selection {
  public:
    STSelection(double st_min);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_st_min;
  };

  class Chi2Selection: public uhh2::Selection {
  public:
    Chi2Selection(uhh2::Context &ctx, std::string label, double chi2_max);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    double m_chi2_max;
    uhh2::Event::Handle<std::vector<BoostedSingleTopHypothesis>> h_hyp;    

  };

  class LeadingTopJetSelection: public uhh2::Selection {
  public:
    LeadingTopJetSelection(TopJetId id_topjet);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    TopJetId m_id_topjet;
  };

  class NGenJetSelection: public uhh2::Selection {
  public:
    NGenJetSelection(unsigned int n_min_, unsigned int n_max_);
    virtual bool passes(const uhh2::Event &event) override;
  private:
    unsigned int n_min, n_max;
  };

  class HOTVRLeptonCleaner: public uhh2::AnalysisModule {
  public:
    HOTVRLeptonCleaner(double deltaRmin);
    virtual bool process(uhh2::Event &event) override;
  private:
    double m_deltaRmin;
  };

}
