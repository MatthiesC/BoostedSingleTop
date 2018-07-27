#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BoostedSingleTop/include/SingleTopGen.h"

/** \brief Histograms for singletop quantities on generator (parton) level
 * 
 * SingleTopGen container has to be filled before calling this histogram class
 */
class SingleTopGenHists: public uhh2::Hists {
public:
    SingleTopGenHists(uhh2::Context & ctx, const std::string & dirname);

    virtual void fill(const uhh2::Event & ev) override;

protected:
    TH1F* shat, *M_top, *Pt_top, *Pt_top_over_shat, *eta_top, *y_top, *deltaR_top_decays, *deltaR_W_decays, *Pt_b, *eta_b, *y_b, *M_b, *Pt_W, *eta_W, *y_W, *M_W, *DecayChannel;

    TH2F* shat_vs_deltaR_top, *Pt_top_vs_deltaR_top, *shat_vs_Pt_top, *Pt_top_vs_Pt_W;
    uhh2::Event::Handle<SingleTopGen> h_singletopgen;
};
