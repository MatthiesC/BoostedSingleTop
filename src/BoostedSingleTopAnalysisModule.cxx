#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/common/include/CommonModules.h"
#include "UHH2/common/include/CleaningModules.h"
#include "UHH2/common/include/MCWeight.h"
#include "UHH2/common/include/TopJetIds.h"
#include "UHH2/common/include/JetIds.h"
#include "UHH2/common/include/NSelections.h"
#include "UHH2/common/include/TTbarReconstruction.h"
#include "UHH2/common/include/JetHists.h" // contains BTagMCEfficiencyHists
#include "UHH2/common/include/TopPtReweight.h"
#include "UHH2/common/include/MuonIds.h"
#include "UHH2/common/include/ElectronIds.h"

#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWch.h"
#include "UHH2/BoostedSingleTop/include/TTbarGen_.h"

#include "UHH2/BoostedSingleTop/include/HOTVRHists.h"
#include "UHH2/BoostedSingleTop/include/HOTVRIds.h"
#include "UHH2/BoostedSingleTop/include/HOTVRScaleFactor.h"
#include "UHH2/BoostedSingleTop/include/HOTVRPileUpUncertainty.h"

#include "UHH2/BoostedSingleTop/include/BoostedSingleTopGen.h"
#include "UHH2/BoostedSingleTop/include/BoostedSingleTopSelections.h"
#include "UHH2/BoostedSingleTop/include/BoostedSingleTopHypothesisDiscriminators.h"
#include "UHH2/BoostedSingleTop/include/AndHists.h"
#include "UHH2/BoostedSingleTop/include/MatchHists.h"
#include "UHH2/BoostedSingleTop/include/MVAHists.h"

#include "UHH2/BoostedSingleTop/include/MVADiscriminator.h"

using namespace std;
using namespace uhh2;

namespace uhh2 {

  class BoostedSingleTopAnalysisModule: public AnalysisModule {
  public:
    
    explicit BoostedSingleTopAnalysisModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:  
    std::string dataset_name;

    // --- Additional Modules --- //
    // Common Modules
    std::unique_ptr<CommonModules> common;
    // Reconsturction
    std::unique_ptr<AnalysisModule> primary_lep;
    std::unique_ptr<AnalysisModule> SingleTopGen_tWchProd, TTbarGenProd;
    Event::Handle<FlavorParticle> h_primlep;
    Event::Handle<TopJet> h_primtopjet;

    // Scale factors & Uncertainties
    std::unique_ptr<AnalysisModule> sf_muon_id, sf_muon_iso, sf_muon_trigger, sf_muon_trk;
    std::unique_ptr<AnalysisModule> sf_ele_id, sf_ele_rec;
    std::unique_ptr<AnalysisModule> sf_btag_medium;

    std::string sys_muon_id, sys_muon_iso, sys_muon_trigger, sys_muon_trk;
    std::string sys_ele_id, sys_ele_rec;
    std::string sys_btag_medium;

    std::unique_ptr<AnalysisModule> sf_top_pt_reweight;
    bool do_top_pt_reweight;
    std::unique_ptr<AnalysisModule> sf_toptag;
    std::string sys_toptag;

    std::unique_ptr<AnalysisModule> scale_variation;
    std::unique_ptr<AnalysisModule> L1_variation;
    string sys_L1;
    bool do_scale_variation;
    string sys_pu;

    std::unique_ptr<Selection> sel_nbjetcut_loose, sel_nbjetcut_medium, sel_1bjetcut_medium, sel_2bjetcut_medium, sel_toptags_1, sel_toptags_0; //sel_2bjetcut_tight, sel_1bjetcut_tight;

    bool is_ele, is_muo;

    TopJetId id_toptag;

    // Cleaner chris
    // std::unique_ptr<AnalysisModule> cl_muon, cl_ele;
  

    // --- Selections and Histogramms --- //
    std::unique_ptr<AndHists> hist_presel, hist_presel_btagSF, hist_nbjetcut_loose, hist_nbjetcut_medium; //hist_2bjetcut_tight, hist_1bjetcut_tight; 
    std::unique_ptr<MatchHists> hist_match_tw, hist_match_tw_1b, hist_match_tw_2b, hist_match_tt;
    std::unique_ptr<MVAHists> hist_mva_0t1b, hist_mva_0t2b, hist_mva_1t, hist_mva_1t1b, hist_mva_1t2b;
    std::unique_ptr<Hists> hist_BTagMCEfficiency;

    // --- Declare new Output for TMVA --- //
    Event::Handle<double> h_tmva_weight;
    Event::Handle<double> h_tmva_n_btags;
    Event::Handle<double> h_tmva_pseudotop_mass;
    Event::Handle<double> h_tmva_deltaphi_bottomlepton;
    Event::Handle<double> h_tmva_deltaphi_bottomtoptag;
    Event::Handle<double> h_tmva_wass_pt;
    Event::Handle<double> h_tmva_ptbalance;
    Event::Handle<double> h_tmva_lepton_pt;
    Event::Handle<double> h_tmva_lepton_eta;
    Event::Handle<double> h_tmva_deltaphi_leptonnextjet;
    Event::Handle<double> h_tmva_top_pt;
    Event::Handle<double> h_tmva_top_eta;
    Event::Handle<double> h_tmva_ht_jets;
    Event::Handle<double> h_tmva_met_pt;

    bool do_mva;
    //MVADiscriminator *discr_ele_, *discr_muo_;
    MVADiscriminator *discr_BDT_, *discr_BDT_2b;
    MVADiscriminator_1b *discr_BDT_1b;
    Event::Handle<double> h_mvadiscr, h_mvadiscr_1b, h_mvadiscr_2b;

  };

  BoostedSingleTopAnalysisModule::BoostedSingleTopAnalysisModule(Context & ctx) {

    is_ele = ctx.get("analysis_channel") == "ELECTRON";
    is_muo = ctx.get("analysis_channel") == "MUON";

    dataset_name = ctx.get("dataset_version");

    
    do_scale_variation = (ctx.get("ScaleVariationMuR") == "up" || ctx.get("ScaleVariationMuR") == "down") || (ctx.get("ScaleVariationMuF") == "up" || ctx.get("ScaleVariationMuF") == "down");
    do_top_pt_reweight = ctx.get("b_TopPtReweight") == "true";
    do_mva = ctx.get("do_mva") == "true";

    sys_muon_id      = ctx.get("Systematic_MuonID");
    sys_muon_trigger = ctx.get("Systematic_MuonTrigger");
    sys_muon_iso     = ctx.get("Systematic_MuonIso");
    sys_muon_trk     = ctx.get("Systematic_MuonTrk");

    sys_ele_id       = ctx.get("Systematic_EleID");
    sys_ele_rec      = ctx.get("Systematic_EleRec");
 
    sys_btag_medium  = ctx.get("Systematic_BTag");
    sys_toptag       = ctx.get("Systematic_TopTag");
    sys_L1           = ctx.get("Systematic_L1");
    sys_pu           = ctx.get("Systematic_PU");

    // --- Selection Variables --- //
    double deltaPhi_min = M_PI/2;  // minimum delta phi between lepton and top

    double top_fpt_max   = 0.8;    // maximum pt fraction of leading subjet
    double top_m_min     = 140.;   // minimum topjet mass
    double top_m_max     = 220.;   // maximum topjet mass
    double top_mpair_min = 50.;    // minimum pairwise mass of first three subjets
    double top_tau32_max = 0.56;   // maximum nsubjetiness tau_3/2
    //double top_pt_min    = 200.0;
    //double top_eta_max   = 2.5;

    //double chi2_max      = 20;     // maximum chi2 of the reconstructed bstar hypothesis

    /*double lep_eta_max = 2.4;
    double lepveto_pt_min = 30.;
    double lep_pt_min  = 150.; // chris
    double muo_iso_max = 0.15;*/



    // - Object IDs - //
    // TopJetId id_toptag = HOTVRTopTag(top_fpt_max, 0, FLT_MAX, top_mpair_min); // hotvr top tag without mass cut
    TopJetId id_topjet =  DeltaPhiCut(deltaPhi_min);
    //cl_topjet.reset(new TopJetCleaner(ctx, id_topjet));
    //sel_ntop.reset(new NTopJetSelection(1, -1));

    /*// chris                                                                                                                                                  
    MuonId id_muo_veto = AndId<Muon>(MuonIDLoose(), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    MuonId id_muo_tight = AndId<Muon>(MuonIDTight(), PtEtaCut(lep_pt_min, lep_eta_max), MuonIso(muo_iso_max));
    ElectronId id_ele_veto = AndId<Electron>(ElectronID_Spring16_veto, PtEtaCut(lepveto_pt_min, lep_eta_max));
    ElectronId id_ele_tight = AndId<Electron>(ElectronID_Spring16_tight, PtEtaCut(lep_pt_min, lep_eta_max));
    cl_muon.reset(new MuonCleaner(id_muo_tight));
    cl_ele.reset(new ElectronCleaner(id_ele_tight));*/
 




    TopJetId id_toptag_without_tau32 = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), DeltaPhiCut(deltaPhi_min));
    TopJetId id_toptag_without_deltaphi = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max));
    TopJetId id_toptag_only_HOTVR = HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min);
    id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max), DeltaPhiCut(deltaPhi_min)); // hotvr top tag with tau_3/2 and delta phi
    // TopJetId id_toptag = AndId<TopJet>(HOTVRTopTag(top_fpt_max, top_m_min, top_m_max, top_mpair_min), Tau32Groomed(top_tau32_max)); // hotvr top tag with tau_3/2 and without delta phi

    CSVBTag::wp btag_wp_loose = CSVBTag::WP_LOOSE;
    JetId id_btag_loose = CSVBTag(btag_wp_loose);
    CSVBTag::wp btag_wp_medium = CSVBTag::WP_MEDIUM;
    JetId id_btag_medium = CSVBTag(btag_wp_medium);
    CSVBTag::wp btag_wp_tight = CSVBTag::WP_TIGHT;
    JetId id_btag_tight = CSVBTag(btag_wp_tight);

    sel_toptags_1.reset(new NTopJetSelection(1, 1, id_toptag));
    sel_toptags_0.reset(new NTopJetSelection(0, 0, id_toptag));

    // --- Setup Additional Modules --- //
    // - Common Modules - //
    common.reset(new CommonModules());
    common->disable_jec();
    common->disable_jersmear();
    common->init(ctx, sys_pu);

    // - Reconstruction - //
    primary_lep.reset(new PrimaryLepton(ctx));
    SingleTopGen_tWchProd.reset(new SingleTopGen_tWchProducer(ctx, "h_singletopgen_twch"));
    TTbarGenProd.reset(new TTbarGen_Producer(ctx, "h_ttbargen"));
    h_primlep = ctx.get_handle<FlavorParticle>("PrimaryLepton");
    h_primtopjet = ctx.get_handle<TopJet>("PrimaryTopJet");

    // - Scale Factors & Uncertainties - //
    if(is_muo) {
      sf_muon_id.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonID_EfficienciesAndSF_average_RunBtoH.root", "MC_NUM_TightID_DEN_genTracks_PAR_pt_eta", 1, "tightID", true, sys_muon_id));
      sf_muon_iso.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonIso_EfficienciesAndSF_average_RunBtoH.root", "TightISO_TightID_pt_eta", 1, "iso", true, sys_muon_iso));
      sf_muon_trigger.reset(new MCMuonScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/MuonTrigger_EfficienciesAndSF_average_RunBtoH.root", "IsoMu24_OR_IsoTkMu24_PtEtaBins", 0.5, "trigger", true, sys_muon_trigger));
      sf_muon_trk.reset(new MCMuonTrkScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/Tracking_EfficienciesAndSF_BCDEFGH.root", 1, "track", sys_muon_trk, "muons"));
    } else if(is_ele) {
      sf_ele_id.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_CutBased_Tight_ID.root", 1, "tightID", sys_ele_id));
      sf_ele_rec.reset(new MCElecScaleFactor(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/common/data/egammaEffi.txt_EGM2D_RecEff_Moriond17.root", 1, "recEff", sys_ele_rec));
    }

    sf_top_pt_reweight.reset(new TopPtReweight(ctx, 0.0615, -0.0005, "", "", do_top_pt_reweight, 1.0)); // https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
   
    sf_btag_medium.reset(new MCBTagScaleFactor(ctx, btag_wp_medium, "jets", sys_btag_medium, "mujets", "incl", "MCBtagEfficienciesMedium"));


    scale_variation.reset(new MCScaleVariation(ctx));
    L1_variation.reset(new HOTVRPileUpUncertainty(ctx, "/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/config/HOTVR_L1Uncertainty.root", sys_L1 ));

    // --- Selections and Histogramms --- //
    hist_presel.reset(new AndHists(ctx, "PreSel"));
    hist_presel->add_hist(new HOTVRHists(ctx, "PreSel_HOTVRtagged", id_toptag));

    hist_presel_btagSF.reset(new AndHists(ctx, "PreSel_btagSF"));
    hist_presel_btagSF->add_hist(new HOTVRHists(ctx, "PreSel_btagSF_HOTVRtagged", id_toptag));

    hist_BTagMCEfficiency.reset(new BTagMCEfficiencyHists(ctx, "BTagMCEfficiency", btag_wp_medium)); // btag_wp_tight, btag_wp_medium

    sel_nbjetcut_medium.reset(new NJetSelection(1, -1, id_btag_medium)); // id_btag_tight
    hist_nbjetcut_medium.reset(new AndHists(ctx, "NBJetCutMedium"));
    hist_nbjetcut_medium->add_hist(new HOTVRHists(ctx, "NBJetCutMedium_HOTVRtagged", id_toptag));

    hist_match_tw.reset(new MatchHists(ctx, "Matching_tW", id_toptag));
    hist_match_tw_1b.reset(new MatchHists(ctx, "Matching_tW_1b", id_toptag));
    hist_match_tw_2b.reset(new MatchHists(ctx, "Matching_tW_2b", id_toptag));
    hist_match_tt.reset(new MatchHists(ctx, "Matching_tt", id_toptag));

    // --- Declare new Output for TMVA --- //
    // ctx.undeclare_all_event_output();

    h_tmva_weight = ctx.declare_event_output<double>("tmva_weight");
    h_tmva_n_btags = ctx.declare_event_output<double>("tmva_n_btags");
    h_tmva_pseudotop_mass = ctx.declare_event_output<double>("tmva_pseudotop_mass");
    h_tmva_deltaphi_bottomlepton = ctx.declare_event_output<double>("tmva_deltaphi_bottomlepton");
    h_tmva_deltaphi_bottomtoptag = ctx.declare_event_output<double>("tmva_deltaphi_bottomtoptag");
    h_tmva_wass_pt = ctx.declare_event_output<double>("tmva_wass_pt");
    h_tmva_ptbalance = ctx.declare_event_output<double>("tmva_ptbalance");
    h_tmva_lepton_pt = ctx.declare_event_output<double>("tmva_lepton_pt");
    h_tmva_lepton_eta = ctx.declare_event_output<double>("tmva_lepton_eta");
    h_tmva_deltaphi_leptonnextjet = ctx.declare_event_output<double>("tmva_deltaphi_leptonnextjet");
    h_tmva_top_pt = ctx.declare_event_output<double>("tmva_top_pt");
    h_tmva_top_eta = ctx.declare_event_output<double>("tmva_top_eta");
    h_tmva_ht_jets = ctx.declare_event_output<double>("tmva_ht_jets");
    h_tmva_met_pt = ctx.declare_event_output<double>("tmva_met_pt");

    //if(do_mva && is_ele)      discr_ele_ = new MVADiscriminator("");
    //else if(do_mva && is_muo) discr_muo_ = new MVADiscriminator("/nfs/dust/cms/user/matthies/BoostedSingleTop/RunII_80X_v5/Analysis/tmva/"); // needs to be updated once the weight file is produced

    if(do_mva) {
      discr_BDT_ =   new MVADiscriminator("/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb/weights_nbtags/TMVAClassification_BDT.weights.xml");
      discr_BDT_1b = new MVADiscriminator_1b("/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb/weights_1btag/TMVAClassification_BDT_100_25_01_10.weights.xml");
      discr_BDT_2b = new MVADiscriminator("/nfs/dust/cms/user/matthies/Analysis_80x_v5/CMSSW_8_0_24_patch1/src/UHH2/BoostedSingleTop/tmva_comb/weights_2btags/TMVAClassification_BDT_50_25_03_10.weights.xml");
      h_mvadiscr    = ctx.declare_event_output<double>("mvadiscr");
      h_mvadiscr_1b = ctx.declare_event_output<double>("mvadiscr_1b");
      h_mvadiscr_2b = ctx.declare_event_output<double>("mvadiscr_2b");
    }

    hist_mva_0t1b.reset(new MVAHists(ctx, "MVA_0t1b"));
    hist_mva_0t2b.reset(new MVAHists(ctx, "MVA_0t2b"));

    hist_mva_1t.reset(new MVAHists(ctx, "MVA_1t"));

    sel_1bjetcut_medium.reset(new NJetSelection(1, 1, id_btag_medium));
    sel_2bjetcut_medium.reset(new NJetSelection(2, -1, id_btag_medium));

    hist_mva_1t1b.reset(new MVAHists(ctx, "MVA_1t1b"));
    hist_mva_1t2b.reset(new MVAHists(ctx, "MVA_1t2b"));

  }

  bool BoostedSingleTopAnalysisModule::process(Event & event) {

     //==============//
    // After PreSel //
   //==============//

    // Apply scale factors for leptons and trigger
    // Muon channel
    if(is_muo) {
      //      cl_muon->process(event); //chris
      sf_muon_id->process(event);
      sf_muon_iso->process(event);
      sf_muon_trigger->process(event);
      sf_muon_trk->process(event);
    }
    // Electron channel
    else if(is_ele) {
      //      cl_ele->process(event); //chris
      sf_ele_id->process(event);
      sf_ele_rec->process(event);
    }

    //L1_variation->process(event);
    if(do_scale_variation) scale_variation->process(event);
    if(do_top_pt_reweight && dataset_name.find("TTbar") == 0) sf_top_pt_reweight->process(event); // just tested for TTbar!!! Do not apply to other MC samples!!!

    common->process(event);
    primary_lep->process(event); // for the identification of the leading lepton -> see TTbarReconstruction.h

    hist_presel->fill(event);


     //===============================================//
    // Cuts on events regarding the number of b-jets //
   //===============================================//

    hist_BTagMCEfficiency->fill(event);

    sf_btag_medium->process(event);

    hist_presel_btagSF->fill(event);

    if(!sel_nbjetcut_medium->passes(event)) return false;
    //sf_btag_medium->process(event);
    hist_nbjetcut_medium->fill(event);


     //=====================================================//
    // Reconstruction: Matching between MC-Gen and MC-Reco //
   //=====================================================//

    /*if(dataset_name.find("SingleTop_T_tWch") == 0 || dataset_name.find("SingleTop_Tbar_tWch") == 0)
      {
	SingleTopGen_tWchProd->process(event);
	hist_match_tw->fill(event);
      }
    else if(dataset_name.find("TTbar") == 0) // find "xyz" at position 0 of name
      {
	TTbarGenProd->process(event);
	hist_match_tt->fill(event);
      }
    */

     //=================================================//
    // Cuts on events regarding the number of top-tags //
   //=================================================//

    if(!(sel_toptags_1->passes(event) || sel_toptags_0->passes(event))) return false;


     //===========================================//
    // Set Output for TMVA (BDT Input Variables) //
   //===========================================//

    vector<TopJet> hotvrJets = *event.topjets;
    sort_by_pt<TopJet>(hotvrJets);
    vector<Jet> jets = *event.jets;
    vector<Jet> bjets;
    vector<TopJet> toptaggedjets;
    for (TopJet hotvrjet : hotvrJets)
      {
	if (id_toptag(hotvrjet, event))
	  {
	    toptaggedjets.push_back(hotvrjet);
	  }
      }
    if(toptaggedjets.size() == 0) toptaggedjets.push_back(hotvrJets.at(0)); // if no top-tag then use leading HOTVR jet!
    if(toptaggedjets.size() != 1) throw logic_error("BoostedSingleTopAnalysisModule: Event has no top-tagged HOTVR jet despite prior selections (should be == 1)");
    TopJet topjet = toptaggedjets.at(0);
    event.set(h_primtopjet, toptaggedjets.at(0)); // VERY IMPORTANT
    LorentzVector met = (*event.met).v4();
    double met_pt = met.pt();
    const Particle& lepton = event.get(h_primlep);
    Jet nextjet;
    nextjet = *nextJet(lepton, *event.jets);

    //b-tags
    const CSVBTag btag_medium(CSVBTag::WP_MEDIUM); // WP_TIGHT
    int n_btags_medium = 0;
    for(Jet jet : jets)
      {
	if(btag_medium(jet, event))
	  {
	    ++n_btags_medium;
	    bjets.push_back(jet);
	  }
      }
    sort_by_pt<Jet>(bjets);
    Jet bjet0; // the b-tag with the largest deltaR to the top-tag
    double deltaR_top_bjet = 0;
    for (Jet jet : bjets)
      {
	if(deltaR_top_bjet < deltaR(jet, topjet))
	  {
	    deltaR_top_bjet = deltaR(jet, topjet);
	    bjet0 = jet;
	  }
      }

    vector<LorentzVector> nus = NeutrinoReconstruction(lepton.v4(), met);
    vector<double> M_LepNuB;
    vector<double> M_Wass, Pt_Wass;
    LorentzVector neutrinoRec;
    for (LorentzVector neutrino : nus)
      {
	// LepNuB reconstruction
	M_LepNuB.push_back((lepton.v4() + bjet0.v4() + neutrino).M()); 
	// Wass reconstruction
	M_Wass.push_back((lepton.v4() + neutrino).M());
	Pt_Wass.push_back((lepton.v4() + neutrino).Pt());
	neutrinoRec = neutrino;
      }
    if(nus.size() == 2)
      {
	bool swap_solutions = M_Wass.at(0) > M_Wass.at(1);
	if (swap_solutions) // choose the solution with the minimal reconstructed W mass               <---------------- very important note
	  {
	    swap(Pt_Wass.at(0), Pt_Wass.at(1));
	    swap(M_LepNuB.at(0), M_LepNuB.at(1));
	    neutrinoRec = nus.at(1);
	  }
      }

    event.set(h_tmva_weight, event.weight);
    event.set(h_tmva_n_btags, n_btags_medium);
    event.set(h_tmva_pseudotop_mass, M_LepNuB.at(0));
    event.set(h_tmva_deltaphi_bottomlepton, deltaPhi(lepton.v4(), bjet0.v4()));
    event.set(h_tmva_deltaphi_bottomtoptag, deltaPhi(topjet, bjet0.v4()));
    event.set(h_tmva_wass_pt, Pt_Wass.at(0));
    double pt_balance = (Pt_Wass.at(0) - topjet.v4().pt())/(topjet.v4().pt());
    event.set(h_tmva_ptbalance, pt_balance);
    event.set(h_tmva_lepton_pt, lepton.v4().pt());
    event.set(h_tmva_lepton_eta, lepton.v4().eta());
    event.set(h_tmva_deltaphi_leptonnextjet, deltaPhi(lepton, nextjet));
    event.set(h_tmva_top_pt, topjet.v4().pt());
    event.set(h_tmva_top_eta, topjet.v4().eta());
    double ht_jets = 0;
    for(Jet jet : *event.jets) ht_jets += jet.v4().pt();
    event.set(h_tmva_ht_jets, ht_jets);
    event.set(h_tmva_met_pt, met_pt);

    std::vector<double> mva_inputvars {double(n_btags_medium), deltaPhi(lepton.v4(), bjet0.v4()), deltaPhi(topjet, bjet0.v4()), M_LepNuB.at(0), pt_balance, Pt_Wass.at(0), lepton.v4().pt(), lepton.v4().eta(), deltaPhi(lepton, nextjet), met_pt};

     //==================//
    // TMVA Application //
   //==================//

    double mvaD = -100;
    double mvaD_1b = -100;
    double mvaD_2b = -100;
    if(do_mva) {
      mvaD = discr_BDT_->eval(float(n_btags_medium), float(deltaPhi(lepton.v4(), bjet0.v4())), float(deltaPhi(topjet, bjet0.v4())), float(M_LepNuB.at(0)), float(pt_balance), float(Pt_Wass.at(0)), float(lepton.v4().pt()), float(lepton.v4().eta()), float(deltaPhi(lepton, nextjet)));//, float(met_pt));
      event.set(h_mvadiscr, mvaD);
      if(mvaD > 0.999) mvaD = 0.999; // if this is not set, purely signal-like events are not shown in the last bin but migrate into the overflow bin...
      mvaD_1b = discr_BDT_1b->eval(float(deltaPhi(lepton.v4(), bjet0.v4())), float(deltaPhi(topjet, bjet0.v4())), float(M_LepNuB.at(0)), float(pt_balance), float(Pt_Wass.at(0)), float(lepton.v4().pt()), float(lepton.v4().eta()), float(deltaPhi(lepton, nextjet)));
      if(mvaD_1b > 0.999) mvaD_1b = 0.999;
      event.set(h_mvadiscr_1b, mvaD_1b);
      mvaD_2b = discr_BDT_2b->eval(float(n_btags_medium), float(deltaPhi(lepton.v4(), bjet0.v4())), float(deltaPhi(topjet, bjet0.v4())), float(M_LepNuB.at(0)), float(pt_balance), float(Pt_Wass.at(0)), float(lepton.v4().pt()), float(lepton.v4().eta()), float(deltaPhi(lepton, nextjet)));
      if(mvaD_2b > 0.999) mvaD_2b = 0.999;
      event.set(h_mvadiscr_2b, mvaD_2b);
    }

    // Get migration plots
    if(sel_toptags_1->passes(event) && (sel_1bjetcut_medium->passes(event) || sel_2bjetcut_medium->passes(event))) {
      if(dataset_name.find("SingleTop_T_tWch") == 0 || dataset_name.find("SingleTop_Tbar_tWch") == 0)
	{
	  SingleTopGen_tWchProd->process(event); // needs to be called before and match hists can be filled!!!
	  hist_match_tw->fill_(event,neutrinoRec);
	  if(sel_1bjetcut_medium->passes(event)) hist_match_tw_1b->fill_(event,neutrinoRec);
	  else if(sel_2bjetcut_medium->passes(event)) hist_match_tw_2b->fill_(event,neutrinoRec);
	}
      else if(dataset_name.find("TTbar") == 0) // find "xyz" at position 0 of name
	{
	  TTbarGenProd->process(event);
	  hist_match_tt->fill_(event,neutrinoRec);
	}
    }

    // Get MVA plots
    if(sel_1bjetcut_medium->passes(event)) {
      if(sel_toptags_0->passes(event)) {
	hist_mva_0t1b->fill_(event, mvaD_1b, mva_inputvars, 0);
	return false;
      } else if(sel_toptags_1->passes(event)) {
	hist_mva_1t1b->fill_(event, mvaD_1b, mva_inputvars, (is_ele ? 0.15 : 0.10));
	//return false;
      }
    }
    else if(sel_2bjetcut_medium->passes(event)) {
      if(sel_toptags_0->passes(event)) {
	hist_mva_0t2b->fill_(event, mvaD_2b, mva_inputvars, 0);
	return false;
      } else if(sel_toptags_1->passes(event)) {
	hist_mva_1t2b->fill_(event, mvaD_2b, mva_inputvars, (is_ele ? 0.60 : 0.45));
	return false;
      }
    }




     //======//
    // Done //
   //======//

    return true;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BoostedSingleTopAnalysisModule)

}
