#include "UHH2/core/include/Hists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWch.h"

/** \brief Histograms for singletop quantities on generator (parton) level
 * 
 * SingleTopGen container has to be filled before calling this histogram class
 */
class SingleTopGen_tWchHists: public uhh2::Hists {
 public:
  SingleTopGen_tWchHists(uhh2::Context & ctx, const std::string & dirname);

  virtual void fill(const uhh2::Event & ev) override;

 protected:
  TH1F *shat, *DecayChannel;

  TH1F *Pt_top, *eta_top, *y_top, *M_top;
  TH1F *Pt_bTop, *eta_bTop, *y_bTop, *M_bTop;
  TH1F *Pt_WTop, *eta_WTop, *y_WTop, *M_WTop;
  TH1F *Pt_WAss, *eta_WAss, *y_WAss, *M_WAss;
  TH1F *Pt_bAss, *eta_bAss, *y_bAss, *M_bAss;
  TH1F *Pt_gluonAss, *eta_gluonAss, *y_gluonAss, *M_gluonAss;
  TH1F *Pt_WTopDecay1, *eta_WTopDecay1, *y_WTopDecay1, *M_WTopDecay1;
  TH1F *Pt_WTopDecay2, *eta_WTopDecay2, *y_WTopDecay2, *M_WTopDecay2;
  TH1F *Pt_WAssDecay1, *eta_WAssDecay1, *y_WAssDecay1, *M_WAssDecay1;
  TH1F *Pt_WAssDecay2, *eta_WAssDecay2, *y_WAssDecay2, *M_WAssDecay2;
  
  TH1F *Pt_top_over_shat, *Pt_bTop_over_shat, *Pt_WTop_over_shat, *Pt_WAss_over_shat, *Pt_WTopDecay1_over_shat, *Pt_WTopDecay2_over_shat,  *Pt_WAssDecay1_over_shat, *Pt_WAssDecay2_over_shat;

  TH1F *deltaR_top_decays;
  TH2F *shat_vs_deltaR_top, *DeltaR_top_vs_Pt_top;

  TH1F *deltaR_WTop_decays;
  TH2F *shat_vs_deltaR_WTop, *Pt_WTop_vs_deltaR_WTop;

  TH1F *deltaR_WAss_decays;
  TH2F *shat_vs_deltaR_WAss, *Pt_WAss_vs_deltaR_WAss;

  TH2F *shat_vs_Pt_top, *shat_vs_Pt_WAss, *Pt_top_vs_Pt_WAss, *Pt_top_vs_Pt_LeptAss;

  TH1F *Pt_WAss_over_Pt_top;

  TH1F *deltaR_bAss_top, *deltaR_bAss_top_decays_min, *deltaR_gluonAss_top, *deltaR_gluonAss_top_decays_min, *deltaR_bAss_WAss_decays_min, *deltaR_gluonAss_WAss_decays_min;
  TH2F *deltaR_bAss_top_vs_deltaR_top, *deltaR_bAss_top_decays_min_vs_deltaR_top, *deltaR_gluonAss_top_vs_deltaR_top, *deltaR_gluonAss_top_decays_min_vs_deltaR_top;
  
  // Electron Channel

  TH1F *EleChannel_lep_pt, *EleChannel_lep_eta, *EleChannel_top_pt, *EleChannel_top_eta;
  TH2F *EleChannel_pt_topVSlep;
  
  TH1F *EleChannel_lep_pt_topPt0to100, *EleChannel_lep_pt_topPt100to200, *EleChannel_lep_pt_topPt200to300, *EleChannel_lep_pt_topPt300to400, *EleChannel_lep_pt_topPt400to500, *EleChannel_lep_pt_topPt500to600, *EleChannel_lep_pt_topPt600to700, *EleChannel_lep_pt_topPt700toInf;

  TH1F *EleChannel_nu_pt;
  TH2F *EleChannel_pt_topVSnu;

  // Muon Channel

  TH1F *MuoChannel_lep_pt, *MuoChannel_lep_eta, *MuoChannel_top_pt, *MuoChannel_top_eta;
  TH2F *MuoChannel_pt_topVSlep;

  TH1F *MuoChannel_lep_pt_topPt0to100, *MuoChannel_lep_pt_topPt100to200, *MuoChannel_lep_pt_topPt200to300, *MuoChannel_lep_pt_topPt300to400, *MuoChannel_lep_pt_topPt400to500, *MuoChannel_lep_pt_topPt500to600, *MuoChannel_lep_pt_topPt600to700, *MuoChannel_lep_pt_topPt700toInf;

  TH1F *MuoChannel_nu_pt;
  TH2F *MuoChannel_pt_topVSnu;

  uhh2::Event::Handle<SingleTopGen_tWch> h_singletopgen_twch;
};
