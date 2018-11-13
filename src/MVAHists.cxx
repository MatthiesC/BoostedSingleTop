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
  h_primtopjet = ctx.get_handle<TopJet>("PrimaryTopJet");

  // Book all histograms here ...

  // BDT response
  //double rebin_mvaD[] = {-1,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,1};
  double rebin_mvaD[] = {-1,-0.6,-0.45,-0.35,-0.25,-0.15,0,0.2,0.4,0.6,1};
  unsigned int rebin_mvaD_bincount = sizeof(rebin_mvaD)/sizeof(*rebin_mvaD)-1;
  double rebin01_mvaD[] = {-1,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1};
  unsigned int rebin01_mvaD_bincount = sizeof(rebin01_mvaD)/sizeof(*rebin01_mvaD)-1;
  H_mvaD = book<TH1F>("mvaD", "BDT response", 40, -1, 1);
  H_mvaD_20bins = book<TH1F>("mvaD_20bins", "BDT response", 20, -1, 1);
  H_mvaD_rebin = book<TH1F>("mvaD_rebin", "BDT response", rebin_mvaD_bincount, rebin_mvaD);
  H_mvaD_rebin01 = book<TH1F>("mvaD_rebin01", "BDT response", rebin01_mvaD_bincount, rebin01_mvaD);
  H_mvaD_1b = book<TH1F>("mvaD_1btags", "BDT response (1t1b)", 40, -1, 1);
  H_mvaD_2b = book<TH1F>("mvaD_2btags", "BDT response (1t2b)", 40, -1, 1);

  double rebin_pt[] = {200,300,400,600,1200};
  unsigned int rebin_pt_bincount = sizeof(rebin_pt)/sizeof(*rebin_pt)-1;
  double rebin_eta[] = {-2.5,-2,-1.5,-1,-0.5,0,0.5,1,1.5,2,2.5};
  unsigned int rebin_eta_bincount = sizeof(rebin_eta)/sizeof(*rebin_eta)-1;
  double rebin_eta_lept[] = {-2.4,-1.8,-1.2,-0.6,0,0.6,1.2,1.8,2.4};
  unsigned int rebin_eta_lept_bincount = sizeof(rebin_eta_lept)/sizeof(*rebin_eta_lept)-1;
  //char[] text_extension ("");
  
  // Top-jet distributions after BDT cut
  H_nocuts_top_pt_1000 = book<TH1F>("nocuts_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_nocuts_top_pt_rebin = book<TH1F>("nocuts_top_pt_rebin", "top-jet p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_nocuts_top_eta_1000 = book<TH1F>("nocuts_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_nocuts_top_eta_rebin = book<TH1F>("nocuts_top_eta_rebin", "top-jet #eta", rebin_eta_bincount, rebin_eta);
  H_nocuts_top_pt_vs_eta_1000 = book<TH2F>("nocuts_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_nocuts_top_pt_vs_eta_rebin = book<TH2F>("nocuts_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", rebin_eta_bincount, rebin_eta, rebin_pt_bincount, rebin_pt);
  H_nocuts_top_mass = book<TH1F>("nocuts_top_mass", "top-jet mass [GeV]", 16, 140, 220);

  H_nocuts_tmva_n_btags = book<TH1F>("nocuts_tmva_n_btags", "number of medium CSV b-tags", 6, 0.5, 6.5);
  H_nocuts_tmva_deltaPhi_leptonbottom = book<TH1F>("nocuts_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad]", 40, 0, M_PI);
  H_nocuts_tmva_deltaPhi_topbottom = book<TH1F>("nocuts_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad]", 40, 0, M_PI);
  H_nocuts_tmva_mass_LepNuB = book<TH1F>("nocuts_tmva_mass_LepNuB", "l#nub invariant mass [GeV]", 40, 0, 1000);
  H_nocuts_tmva_pt_balance = book<TH1F>("nocuts_tmva_pt_balance", "p_{T} balance", 40, -2, 2);
  H_nocuts_tmva_Wass_pt = book<TH1F>("nocuts_tmva_Wass_pt", "associated W p_{T} [GeV]", 40, 0, 1000);
  H_nocuts_tmva_lepton_pt = book<TH1F>("nocuts_tmva_lepton_pt", "lepton p_{T} [GeV]", 48, 0, 600);
  H_nocuts_tmva_lepton_eta = book<TH1F>("nocuts_tmva_lepton_eta", "lepton #eta", 25, -2.5, 2.5);
  H_nocuts_tmva_deltaPhi_leptonnextjet = book<TH1F>("nocuts_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad]", 40, 0, M_PI);
  H_nocuts_tmva_met_pt = book<TH1F>("nocuts_tmva_met_pt", "missing E_{T} [GeV]", 48, 0, 600);

  H_BDTcut_top_pt_1000 = book<TH1F>("BDTcut_top_pt_1000", "top-jet p_{T} [GeV]", 1000, 200, 1200);
  H_BDTcut_top_pt_rebin = book<TH1F>("BDTcut_top_pt_rebin", "top-jet p_{T} [GeV]", rebin_pt_bincount, rebin_pt);
  H_BDTcut_top_eta_1000 = book<TH1F>("BDTcut_top_eta_1000", "top-jet #eta", 1000, -2.5, 2.5);
  H_BDTcut_top_eta_rebin = book<TH1F>("BDTcut_top_eta_rebin", "top-jet #eta", rebin_eta_bincount, rebin_eta);
  H_BDTcut_top_pt_vs_eta_1000 = book<TH2F>("BDTcut_top_pt_vs_eta_1000", "top-jet p_{T} [GeV] vs. #eta", 1000, -2.5, 2.5, 1000, 200, 1200);
  H_BDTcut_top_pt_vs_eta_rebin = book<TH2F>("BDTcut_top_pt_vs_eta_rebin", "top-jet p_{T} [GeV] vs. #eta", rebin_eta_bincount, rebin_eta, rebin_pt_bincount, rebin_pt);
  H_BDTcut_top_mass = book<TH1F>("BDTcut_top_mass", "top-jet mass [GeV]", 16, 140, 220);

  H_BDTcut_tmva_n_btags = book<TH1F>("BDTcut_tmva_n_btags", "number of medium CSV b-tags", 6, 0.5, 6.5);
  H_BDTcut_tmva_deltaPhi_leptonbottom = book<TH1F>("BDTcut_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad]", 40, 0, M_PI);
  H_BDTcut_tmva_deltaPhi_topbottom = book<TH1F>("BDTcut_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad]", 40, 0, M_PI);
  H_BDTcut_tmva_mass_LepNuB = book<TH1F>("BDTcut_tmva_mass_LepNuB", "l#nub invariant mass [GeV]", 40, 0, 1000);
  H_BDTcut_tmva_pt_balance = book<TH1F>("BDTcut_tmva_pt_balance", "p_{T} balance", 40, -2, 2);
  H_BDTcut_tmva_Wass_pt = book<TH1F>("BDTcut_tmva_Wass_pt", "associated W p_{T} [GeV]", 40, 0, 1000);
  H_BDTcut_tmva_lepton_pt = book<TH1F>("BDTcut_tmva_lepton_pt", "lepton p_{T} [GeV]", 48, 0, 600);
  H_BDTcut_tmva_lepton_eta = book<TH1F>("BDTcut_tmva_lepton_eta", "lepton #eta", 25, -2.5, 2.5);
  H_BDTcut_tmva_deltaPhi_leptonnextjet = book<TH1F>("BDTcut_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad]", 40, 0, M_PI);
  H_BDTcut_tmva_met_pt = book<TH1F>("BDTcut_tmva_met_pt", "missing E_{T} [GeV]", 48, 0, 600);

  // now go through top_pt bins and write out BDT response for each
  //text_extension = " (top-jet p_{T} #in (200,300] GeV)";
  H_200to300_BDT = book<TH1F>("200to300_mvaD", "BDT response (top-jet p_{T} #in (200,300] GeV)", 40, -1, 1);
  H_200to300_BDT_20bins = book<TH1F>("200to300_mvaD_20bins", "BDT response (top-jet p_{T} #in (200,300] GeV)", 20, -1, 1);
  H_200to300_BDT_rebin = book<TH1F>("200to300_mvaD_rebin", "BDT response (top-jet p_{T} #in (200,300] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_200to300_BDT_rebin01 = book<TH1F>("200to300_mvaD_rebin01", "BDT response (top-jet p_{T} #in (200,300] GeV)", rebin01_mvaD_bincount, rebin01_mvaD);
  H_200to300_top_eta = book<TH1F>("200to300_top_eta", "top-jet #eta (top-jet p_{T} #in (200,300] GeV)", rebin_eta_bincount, rebin_eta);
  H_200to300_tmva_n_btags = book<TH1F>("200to300_tmva_n_btags", "number of medium CSV b-tags (top-jet p_{T} #in (200,300] GeV)", 6, 0.5, 6.5);
  H_200to300_tmva_deltaPhi_leptonbottom = book<TH1F>("200to300_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad] (top-jet p_{T} #in (200,300] GeV)", 40, 0, M_PI);
  H_200to300_tmva_deltaPhi_topbottom = book<TH1F>("200to300_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad] (top-jet p_{T} #in (200,300] GeV)", 40, 0, M_PI);
  H_200to300_tmva_mass_LepNuB = book<TH1F>("200to300_tmva_mass_LepNuB", "l#nub invariant mass [GeV] (top-jet p_{T} #in (200,300] GeV)", 40, 0, 1000);
  H_200to300_tmva_pt_balance = book<TH1F>("200to300_tmva_pt_balance", "p_{T} balance (top-jet p_{T} #in (200,300] GeV)", 40, -2, 2);
  H_200to300_tmva_Wass_pt = book<TH1F>("200to300_tmva_Wass_pt", "associated W p_{T} [GeV] (top-jet p_{T} #in (200,300] GeV)", 40, 0, 1000);
  H_200to300_tmva_lepton_pt = book<TH1F>("200to300_tmva_lepton_pt", "lepton p_{T} [GeV] (top-jet p_{T} #in (200,300] GeV)", 48, 0, 600);
  H_200to300_tmva_lepton_eta = book<TH1F>("200to300_tmva_lepton_eta", "lepton #eta (top-jet p_{T} #in (200,300] GeV)", 25, -2.5, 2.5);
  H_200to300_tmva_deltaPhi_leptonnextjet = book<TH1F>("200to300_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad] (top-jet p_{T} #in (200,300] GeV)", 40, 0, M_PI);
  H_200to300_tmva_met_pt = book<TH1F>("200to300_tmva_met_pt", "missing E_{T} [GeV] (top-jet p_{T} #in (200,300] GeV)", 48, 0, 600);

  //text_extension = " (top-jet p_{T} #in (300,400] GeV)";
  H_300to400_BDT = book<TH1F>("300to400_mvaD", "BDT response (top-jet p_{T} #in (300,400] GeV)", 40, -1, 1);
  H_300to400_BDT_20bins = book<TH1F>("300to400_mvaD_20bins", "BDT response (top-jet p_{T} #in (300,400] GeV)", 20, -1, 1);
  H_300to400_BDT_rebin = book<TH1F>("300to400_mvaD_rebin", "BDT response (top-jet p_{T} #in (300,400] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_300to400_BDT_rebin01 = book<TH1F>("300to400_mvaD_rebin01", "BDT response (top-jet p_{T} #in (300,400] GeV)", rebin01_mvaD_bincount, rebin01_mvaD);
  H_300to400_top_eta = book<TH1F>("300to400_top_eta", "top-jet #eta (top-jet p_{T} #in (300,400] GeV)", rebin_eta_bincount, rebin_eta);
  H_300to400_tmva_n_btags = book<TH1F>("300to400_tmva_n_btags", "number of medium CSV b-tags (top-jet p_{T} #in (300,400] GeV)", 6, 0.5, 6.5);
  H_300to400_tmva_deltaPhi_leptonbottom = book<TH1F>("300to400_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad] (top-jet p_{T} #in (300,400] GeV)", 40, 0, M_PI);
  H_300to400_tmva_deltaPhi_topbottom = book<TH1F>("300to400_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad] (top-jet p_{T} #in (300,400] GeV)", 40, 0, M_PI);
  H_300to400_tmva_mass_LepNuB = book<TH1F>("300to400_tmva_mass_LepNuB", "l#nub invariant mass [GeV] (top-jet p_{T} #in (300,400] GeV)", 40, 0, 1000);
  H_300to400_tmva_pt_balance = book<TH1F>("300to400_tmva_pt_balance", "p_{T} balance (top-jet p_{T} #in (300,400] GeV)", 40, -2, 2);
  H_300to400_tmva_Wass_pt = book<TH1F>("300to400_tmva_Wass_pt", "associated W p_{T} [GeV] (top-jet p_{T} #in (300,400] GeV)", 40, 0, 1000);
  H_300to400_tmva_lepton_pt = book<TH1F>("300to400_tmva_lepton_pt", "lepton p_{T} [GeV] (top-jet p_{T} #in (300,400] GeV)", 48, 0, 600);
  H_300to400_tmva_lepton_eta = book<TH1F>("300to400_tmva_lepton_eta", "lepton #eta (top-jet p_{T} #in (300,400] GeV)", 25, -2.5, 2.5);
  H_300to400_tmva_deltaPhi_leptonnextjet = book<TH1F>("300to400_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad] (top-jet p_{T} #in (300,400] GeV)", 40, 0, M_PI);
  H_300to400_tmva_met_pt = book<TH1F>("300to400_tmva_met_pt", "missing E_{T} [GeV] (top-jet p_{T} #in (300,400] GeV)", 48, 0, 600);

  //text_extension = " (top-jet p_{T} #in (400,600] GeV)";
  H_400to600_BDT = book<TH1F>("400to600_mvaD", "BDT response (top-jet p_{T} #in (400,600] GeV)", 40, -1, 1);
  H_400to600_BDT_20bins = book<TH1F>("400to600_mvaD_20bins", "BDT response (top-jet p_{T} #in (400,600] GeV)", 20, -1, 1);
  H_400to600_BDT_rebin = book<TH1F>("400to600_mvaD_rebin", "BDT response (top-jet p_{T} #in (400,600] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_400to600_BDT_rebin01 = book<TH1F>("400to600_mvaD_rebin01", "BDT response (top-jet p_{T} #in (400,600] GeV)", rebin01_mvaD_bincount, rebin01_mvaD);
  H_400to600_top_eta = book<TH1F>("400to600_top_eta", "top-jet #eta (top-jet p_{T} #in (400,600] GeV)", rebin_eta_bincount, rebin_eta);
  H_400to600_tmva_n_btags = book<TH1F>("400to600_tmva_n_btags", "number of medium CSV b-tags (top-jet p_{T} #in (400,600] GeV)", 6, 0.5, 6.5);
  H_400to600_tmva_deltaPhi_leptonbottom = book<TH1F>("400to600_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad] (top-jet p_{T} #in (400,600] GeV)", 40, 0, M_PI);
  H_400to600_tmva_deltaPhi_topbottom = book<TH1F>("400to600_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad] (top-jet p_{T} #in (400,600] GeV)", 40, 0, M_PI);
  H_400to600_tmva_mass_LepNuB = book<TH1F>("400to600_tmva_mass_LepNuB", "l#nub invariant mass [GeV] (top-jet p_{T} #in (400,600] GeV)", 40, 0, 1000);
  H_400to600_tmva_pt_balance = book<TH1F>("400to600_tmva_pt_balance", "p_{T} balance (top-jet p_{T} #in (400,600] GeV)", 40, -2, 2);
  H_400to600_tmva_Wass_pt = book<TH1F>("400to600_tmva_Wass_pt", "associated W p_{T} [GeV] (top-jet p_{T} #in (400,600] GeV)", 40, 0, 1000);
  H_400to600_tmva_lepton_pt = book<TH1F>("400to600_tmva_lepton_pt", "lepton p_{T} [GeV] (top-jet p_{T} #in (400,600] GeV)", 48, 0, 600);
  H_400to600_tmva_lepton_eta = book<TH1F>("400to600_tmva_lepton_eta", "lepton #eta (top-jet p_{T} #in (400,600] GeV)", 25, -2.5, 2.5);
  H_400to600_tmva_deltaPhi_leptonnextjet = book<TH1F>("400to600_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad] (top-jet p_{T} #in (400,600] GeV)", 40, 0, M_PI);
  H_400to600_tmva_met_pt = book<TH1F>("400to600_tmva_met_pt", "missing E_{T} [GeV] (top-jet p_{T} #in (400,600] GeV)", 48, 0, 600);

  //text_extension = " (top-jet p_{T} #in (600,1200] GeV)";
  H_600to1200_BDT = book<TH1F>("600to1200_mvaD", "BDT response (top-jet p_{T} #in (600,1200] GeV)", 40, -1, 1);
  H_600to1200_BDT_20bins = book<TH1F>("600to1200_mvaD_20bins", "BDT response (top-jet p_{T} #in (600,1200] GeV)", 20, -1, 1);
  H_600to1200_BDT_rebin = book<TH1F>("600to1200_mvaD_rebin", "BDT response (top-jet p_{T} #in (600,1200] GeV)", rebin_mvaD_bincount, rebin_mvaD);
  H_600to1200_BDT_rebin01 = book<TH1F>("600to1200_mvaD_rebin01", "BDT response (top-jet p_{T} #in (600,1200] GeV)", rebin01_mvaD_bincount, rebin01_mvaD);
  H_600to1200_top_eta = book<TH1F>("600to1200_top_eta", "top-jet #eta (top-jet p_{T} #in (600,1200] GeV)", rebin_eta_bincount, rebin_eta);
  H_600to1200_tmva_n_btags = book<TH1F>("600to1200_tmva_n_btags", "number of medium CSV b-tags (top-jet p_{T} #in (600,1200] GeV)", 6, 0.5, 6.5);
  H_600to1200_tmva_deltaPhi_leptonbottom = book<TH1F>("600to1200_tmva_deltaPhi_leptonbottom", "#Delta#phi(l,b) [rad] (top-jet p_{T} #in (600,1200] GeV)", 40, 0, M_PI);
  H_600to1200_tmva_deltaPhi_topbottom = book<TH1F>("600to1200_tmva_deltaPhi_topbottom", "#Delta#phi(t,b) [rad] (top-jet p_{T} #in (600,1200] GeV)", 40, 0, M_PI);
  H_600to1200_tmva_mass_LepNuB = book<TH1F>("600to1200_tmva_mass_LepNuB", "l#nub invariant mass [GeV] (top-jet p_{T} #in (600,1200] GeV)", 40, 0, 1000);
  H_600to1200_tmva_pt_balance = book<TH1F>("600to1200_tmva_pt_balance", "p_{T} balance (top-jet p_{T} #in (600,1200] GeV)", 40, -2, 2);
  H_600to1200_tmva_Wass_pt = book<TH1F>("600to1200_tmva_Wass_pt", "associated W p_{T} [GeV] (top-jet p_{T} #in (600,1200] GeV)", 40, 0, 1000);
  H_600to1200_tmva_lepton_pt = book<TH1F>("600to1200_tmva_lepton_pt", "lepton p_{T} [GeV] (top-jet p_{T} #in (600,1200] GeV)", 48, 0, 600);
  H_600to1200_tmva_lepton_eta = book<TH1F>("600to1200_tmva_lepton_eta", "lepton #eta (top-jet p_{T} #in (600,1200] GeV)", 25, -2.5, 2.5);
  H_600to1200_tmva_deltaPhi_leptonnextjet = book<TH1F>("600to1200_tmva_deltaPhi_leptonnextjet", "#Delta#phi(l,jet) [rad] (top-jet p_{T} #in (600,1200] GeV)", 40, 0, M_PI);
  H_600to1200_tmva_met_pt = book<TH1F>("600to1200_tmva_met_pt", "missing E_{T} [GeV] (top-jet p_{T} #in (600,1200] GeV)", 48, 0, 600);

}

void MVAHists::fill(const Event & e)
{
  double w = e.weight;
}

void MVAHists::fill_(const Event & e, const double & mvaD, const vector<double> & mva_inputvars, const double & BDTcut)
{
  double w = e.weight;
  //TopJet topjet = (*e.topjets).at(0); // that's wrong!
  const TopJet& topjet = e.get(h_primtopjet);
  const Particle& lepton = e.get(h_primlep);

  double top_pt = topjet.v4().pt();
  double top_eta = topjet.v4().eta();
  double top_mass = topjet.v4().M();

  // import tmva's input variables
  /*double tmva_n_btags, tmva_deltaPhi_leptonbottom, tmva_deltaPhi_topbottom, tmva_mass_LepNuB, tmva_pt_balance, tmva_Wass_pt, tmva_lepton_pt, tmva_lepton_eta, tmva_deltaPhi_leptonnextjet, tmva_met_pt;
  tmva_n_btags                = mva_inputvars.at(0);
  tmva_deltaPhi_leptonbottom  = mva_inputvars.at(1);
  tmva_deltaPhi_topbottom     = mva_inputvars.at(2);
  tmva_mass_LepNuB            = mva_inputvars.at(3);
  tmva_pt_balance             = mva_inputvars.at(4);
  tmva_Wass_pt                = mva_inputvars.at(5);
  tmva_lepton_pt              = mva_inputvars.at(6);
  tmva_lepton_eta             = mva_inputvars.at(7);
  tmva_deltaPhi_leptonnextjet = mva_inputvars.at(8);
  tmva_met_pt                 = mva_inputvars.at(9);*/
  

  // Now fill the histograms ...
  H_mvaD->Fill(mvaD, w);
  H_mvaD_20bins->Fill(mvaD, w);
  H_mvaD_rebin->Fill(mvaD, w);
  H_mvaD_rebin01->Fill(mvaD, w);
  if(abs(mva_inputvars.at(0)-1) < 0.01) H_mvaD_1b->Fill(mvaD, w);
  else if(abs(mva_inputvars.at(0)-2) < 0.01) H_mvaD_2b->Fill(mvaD, w);
  
  // before any BDT cuts ...
  H_nocuts_top_pt_1000->Fill(top_pt, w);
  H_nocuts_top_pt_rebin->Fill(top_pt, w);
  H_nocuts_top_eta_1000->Fill(top_eta, w);
  H_nocuts_top_eta_rebin->Fill(top_eta, w);
  H_nocuts_top_pt_vs_eta_1000->Fill(top_eta, top_pt, w);
  H_nocuts_top_pt_vs_eta_rebin->Fill(top_eta, top_pt, w);
  H_nocuts_top_mass->Fill(top_mass, w);

  H_nocuts_tmva_n_btags->Fill(mva_inputvars.at(0), w);
  H_nocuts_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
  H_nocuts_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
  H_nocuts_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
  H_nocuts_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
  H_nocuts_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
  H_nocuts_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
  H_nocuts_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
  H_nocuts_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
  H_nocuts_tmva_met_pt->Fill(mva_inputvars.at(9), w);

  // now apply a BDT cut ...
  if(mvaD > BDTcut)
    {
      H_BDTcut_top_pt_1000->Fill(top_pt, w);
      H_BDTcut_top_pt_rebin->Fill(top_pt, w);
      H_BDTcut_top_eta_1000->Fill(top_eta, w);
      H_BDTcut_top_eta_rebin->Fill(top_eta, w);
      H_BDTcut_top_pt_vs_eta_1000->Fill(top_eta, top_pt, w);
      H_BDTcut_top_pt_vs_eta_rebin->Fill(top_eta, top_pt, w);
      H_BDTcut_top_mass->Fill(top_mass, w);

      H_BDTcut_tmva_n_btags->Fill(mva_inputvars.at(0), w);
      H_BDTcut_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
      H_BDTcut_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
      H_BDTcut_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
      H_BDTcut_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
      H_BDTcut_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
      H_BDTcut_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
      H_BDTcut_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
      H_BDTcut_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
      H_BDTcut_tmva_met_pt->Fill(mva_inputvars.at(9), w);
    }

  // now go through top_pt bins and write out BDT response for each

  if(top_pt > 200 && top_pt <= 300)
    {
      H_200to300_BDT->Fill(mvaD, w);
      H_200to300_BDT_20bins->Fill(mvaD, w);
      H_200to300_BDT_rebin->Fill(mvaD, w);
      H_200to300_BDT_rebin01->Fill(mvaD, w);
      H_200to300_top_eta->Fill(top_eta, w);
      H_200to300_tmva_n_btags->Fill(mva_inputvars.at(0), w);
      H_200to300_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
      H_200to300_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
      H_200to300_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
      H_200to300_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
      H_200to300_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
      H_200to300_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
      H_200to300_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
      H_200to300_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
      H_200to300_tmva_met_pt->Fill(mva_inputvars.at(9), w);
    }
  else if(top_pt > 300 && top_pt <= 400)
    {
      H_300to400_BDT->Fill(mvaD, w);
      H_300to400_BDT_20bins->Fill(mvaD, w);
      H_300to400_BDT_rebin->Fill(mvaD, w);
      H_300to400_BDT_rebin01->Fill(mvaD, w);
      H_300to400_top_eta->Fill(top_eta, w);
      H_300to400_tmva_n_btags->Fill(mva_inputvars.at(0), w);
      H_300to400_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
      H_300to400_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
      H_300to400_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
      H_300to400_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
      H_300to400_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
      H_300to400_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
      H_300to400_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
      H_300to400_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
      H_300to400_tmva_met_pt->Fill(mva_inputvars.at(9), w);
    }
  else if(top_pt > 400 && top_pt <= 600)
    {
      H_400to600_BDT->Fill(mvaD, w);
      H_400to600_BDT_20bins->Fill(mvaD, w);
      H_400to600_BDT_rebin->Fill(mvaD, w);
      H_400to600_BDT_rebin01->Fill(mvaD, w);
      H_400to600_top_eta->Fill(top_eta, w);
      H_400to600_tmva_n_btags->Fill(mva_inputvars.at(0), w);
      H_400to600_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
      H_400to600_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
      H_400to600_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
      H_400to600_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
      H_400to600_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
      H_400to600_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
      H_400to600_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
      H_400to600_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
      H_400to600_tmva_met_pt->Fill(mva_inputvars.at(9), w);
    }
  else if(top_pt > 600 && top_pt <= 1200)
    {
      H_600to1200_BDT->Fill(mvaD, w);
      H_600to1200_BDT_20bins->Fill(mvaD, w);
      H_600to1200_BDT_rebin->Fill(mvaD, w);
      H_600to1200_BDT_rebin01->Fill(mvaD, w);
      H_600to1200_top_eta->Fill(top_eta, w);
      H_600to1200_tmva_n_btags->Fill(mva_inputvars.at(0), w);
      H_600to1200_tmva_deltaPhi_leptonbottom->Fill(mva_inputvars.at(1), w);
      H_600to1200_tmva_deltaPhi_topbottom->Fill(mva_inputvars.at(2), w);
      H_600to1200_tmva_mass_LepNuB->Fill(mva_inputvars.at(3), w);
      H_600to1200_tmva_pt_balance->Fill(mva_inputvars.at(4), w);
      H_600to1200_tmva_Wass_pt->Fill(mva_inputvars.at(5), w);
      H_600to1200_tmva_lepton_pt->Fill(mva_inputvars.at(6), w);
      H_600to1200_tmva_lepton_eta->Fill(mva_inputvars.at(7), w);
      H_600to1200_tmva_deltaPhi_leptonnextjet->Fill(mva_inputvars.at(8), w);
      H_600to1200_tmva_met_pt->Fill(mva_inputvars.at(9), w);
    }
}

MVAHists::~MVAHists(){}
