#pragma once

#include "UHH2/core/include/Hists.h"

#include "TH1F.h"

namespace uhh2 {

  class MVAHists: public uhh2::Hists {
  public:
    MVAHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & e) override;
    void fill_(const uhh2::Event & e, const double & mvaD, const int & n_btags);
    virtual ~MVAHists();

  protected:

    // Declare all histograms here ...
    TH1F *H_mvaD, *H_mvaD_1b, *H_mvaD_2b;

  };
   
}
