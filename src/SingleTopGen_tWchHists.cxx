#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWchHists.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;

SingleTopGen_tWchHists::SingleTopGen_tWchHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

  // General Event Properties

  shat = book< TH1F>( "shat", "#hat{s} [GeV]", 1000, 0, 5000 ) ;
  DecayChannel = book< TH1F>( "DecayChannel", "decay channel", 17, -0.5, 16.5 ) ;
  
  // GenParticles

  Pt_top = book<TH1F>( "Pt_top", "p_{T}^{top} [GeV]", 1000, 0, 2000 );
  eta_top = book<TH1F>( "eta_top", "#eta_{top}", 1000, -5, 5 ) ;
  y_top = book<TH1F>( "y_top", "y_{top}", 1000, -5, 5 ) ;
  M_top = book<TH1F>("M_top", "M_{top} [GeV]", 1000, 150, 200) ;
       
  Pt_bTop = book<TH1F>( "Pt_bTop","p_{T}^{bTop} [GeV]",1000,0,1000);
  eta_bTop = book<TH1F>( "eta_bTop","#eta_{bTop}",1000,-5,5);
  y_bTop = book<TH1F>( "y_bTop","y_{bTop}",1000,-5,5);
  M_bTop = book<TH1F>( "M_bTop","M_{bTop} [GeV]",1000,0,10);

  Pt_WTop = book<TH1F>( "Pt_WTop","p_{T}^{WTop} [GeV]",1000,0,1000);
  eta_WTop = book<TH1F>( "eta_WTop","#eta_{WTop}",1000,-5,5);
  y_WTop = book<TH1F>( "y_WTop","y_{WTop}",1000,-5,5);
  M_WTop = book<TH1F>( "M_WTop","M_{WTop} [GeV]",1000,50,100);

  Pt_WAss = book<TH1F>( "Pt_WAss","p_{T}^{WAss} [GeV]",1000,0,1000);
  eta_WAss = book<TH1F>( "eta_WAss","#eta_{WAss}",1000,-5,5);
  y_WAss = book<TH1F>( "y_WAss","y_{WAss}",1000,-5,5);
  M_WAss = book<TH1F>( "M_WAss","M_{WAss} [GeV]",1000,50,100);

  Pt_bAss = book<TH1F>( "Pt_bAss","p_{T}^{bAss} [GeV]",1000,0,1000);
  eta_bAss = book<TH1F>( "eta_bAss","#eta_{bAss}",1000,-5,5);
  y_bAss = book<TH1F>( "y_bAss","y_{bAss}",1000,-5,5);
  M_bAss = book<TH1F>( "M_bAss","M_{bAss} [GeV]",1000,0,10);
    
  Pt_gluonAss = book<TH1F>( "Pt_gluonAss","p_{T}^{gluonAss} [GeV]",1000,0,1000);
  eta_gluonAss = book<TH1F>( "eta_gluonAss","#eta_{gluonAss}",1000,-5,5);
  y_gluonAss = book<TH1F>( "y_gluonAss","y_{gluonAss}",1000,-5,5);
  M_gluonAss = book<TH1F>( "M_gluonAss","M_{gluonAss} [GeV]",1000,0,10);

  Pt_WTopDecay1 = book<TH1F>( "Pt_WTopDecay1","p_{T}^{WTopDecay1} [GeV]",1000,0,1000);
  eta_WTopDecay1 = book<TH1F>( "eta_WTopDecay1","#eta_{WTopDecay1}",1000,-5,5);
  y_WTopDecay1 = book<TH1F>( "y_WTopDecay1","y_{WTopDecay1}",1000,-5,5);
  M_WTopDecay1 = book<TH1F>( "M_WTopDecay1","M_{WTopDecay1} [GeV]",1000,0,10);

  Pt_WTopDecay2 = book<TH1F>( "Pt_WTopDecay2","p_{T}^{WTopDecay2} [GeV]",1000,0,1000);
  eta_WTopDecay2 = book<TH1F>( "eta_WTopDecay2","#eta_{WTopDecay2}",1000,-5,5);
  y_WTopDecay2 = book<TH1F>( "y_WTopDecay2","y_{WTopDecay2}",1000,-5,5);
  M_WTopDecay2 = book<TH1F>( "M_WTopDecay2","M_{WTopDecay2} [GeV]",1000,0,10);

  Pt_WAssDecay1 = book<TH1F>( "Pt_WAssDecay1","p_{T}^{WAssDecay1} [GeV]",1000,0,1000);
  eta_WAssDecay1 = book<TH1F>( "eta_WAssDecay1","#eta_{WAssDecay1}",1000,-5,5);
  y_WAssDecay1 = book<TH1F>( "y_WAssDecay1","y_{WAssDecay1}",1000,-5,5);
  M_WAssDecay1 = book<TH1F>( "M_WAssDecay1","M_{WAssDecay1} [GeV]",1000,0,10);

  Pt_WAssDecay2 = book<TH1F>( "Pt_WAssDecay2","p_{T}^{WAssDecay2} [GeV]",1000,0,1000);
  eta_WAssDecay2 = book<TH1F>( "eta_WAssDecay2","#eta_{WAssDecay2}",1000,-5,5);
  y_WAssDecay2 = book<TH1F>( "y_WAssDecay2","y_{WAssDecay2}",1000,-5,5);
  M_WAssDecay2 = book<TH1F>( "M_WAssDecay2","M_{WAssDecay2} [GeV]",1000,0,10);

  // more

  Pt_top_over_shat = book<TH1F>( "Pt_top_over_shat", "p_{T}^{top}/#hat{s}", 1000, 0, 1 );
  Pt_bTop_over_shat = book<TH1F>( "Pt_bTop_over_shat", "p_{T}^{bTop}/#hat{s}", 1000, 0, 1 );
  Pt_WTop_over_shat = book<TH1F>( "Pt_WTop_over_shat", "p_{T}^{WTop}/#hat{s}", 1000, 0, 1 );
  Pt_WAss_over_shat = book<TH1F>( "Pt_WAss_over_shat", "p_{T}^{WAss}/#hat{s}", 1000, 0, 1 );
  Pt_WTopDecay1_over_shat = book<TH1F>( "Pt_WTopDecay1_over_shat", "p_{T}^{WTopDecay1}/#hat{s}", 1000, 0, 1 );
  Pt_WTopDecay2_over_shat = book<TH1F>( "Pt_WTopDecay2_over_shat", "p_{T}^{WTopDecay2}/#hat{s}", 1000, 0, 1 );
  Pt_WAssDecay1_over_shat = book<TH1F>( "Pt_WAssDecay1_over_shat", "p_{T}^{WAssDecay1}/#hat{s}", 1000, 0, 1 );
  Pt_WAssDecay2_over_shat = book<TH1F>( "Pt_WAssDecay2_over_shat", "p_{T}^{WAssDecay2}/#hat{s}", 1000, 0, 1 );

  double ptRebin[16] = {0,50,100,150,200,250,300,350,400,450,500,600,800,1000,1500,2000};
  double deltaRRebin[26] = {0.0,0.2,0.4,0.6,0.8,
			    1.0,1.2,1.4,1.6,1.8,
			    2.0,2.2,2.4,2.6,2.8,
			    3.0,3.2,3.4,3.6,3.8,
			    4.0,4.2,4.4,4.6,4.8,5.0};
  unsigned int ptRebin_size = sizeof(ptRebin)/sizeof(*ptRebin)-1;
  unsigned int deltaRRebin_size = sizeof(deltaRRebin)/sizeof(*deltaRRebin)-1;

  deltaR_top_decays = book<TH1F>( "deltaR_top_decays", "#DeltaR_{max}(top decay prod.)",1000,0,5);   
  shat_vs_deltaR_top = book<TH2F>( "shat_vs_deltaR_top", "#hat{s} [GeV] vs. #DeltaR_{max}(top decay prod.)",500,0,5000,500,0,5);
  DeltaR_top_vs_Pt_top = book<TH2F>( "DeltaR_top_vs_Pt_top", "#DeltaR_{max}(top decay prod.) vs. p_{T}^{top} [GeV]",ptRebin_size,ptRebin,deltaRRebin_size,deltaRRebin);
 
  deltaR_WTop_decays = book<TH1F>( "deltaR_WTop_decays", "#DeltaR(WTop decay prod.)",1000,0,5);
  shat_vs_deltaR_WTop = book<TH2F>( "shat_vs_deltaR_WTop", "#hat{s} [GeV] vs. #DeltaR(WTop decay prod.)",500,0,5000,500,0,5);
  Pt_WTop_vs_deltaR_WTop = book<TH2F>( "Pt_WTop_vs_deltaR_WTop", "p_{T}^{WTop} [GeV] vs. #DeltaR(WTop decay prod.)",500,0,2000,500,0,5);

  deltaR_WAss_decays = book<TH1F>( "deltaR_WAss_decays", "#DeltaR(WAss decay prod.)",1000,0,5);
  shat_vs_deltaR_WAss = book<TH2F>( "shat_vs_deltaR_WAss", "#hat{s} [GeV] vs. #DeltaR(WAss decay prod.)",500,0,5000,500,0,5);
  Pt_WAss_vs_deltaR_WAss = book<TH2F>( "Pt_WAss_vs_deltaR_WAss", "p_{T}^{WAss} [GeV] vs. #DeltaR(WAss decay prod.)",500,0,2000,500,0,5);

  shat_vs_Pt_top = book<TH2F>( "shat_vs_Pt_top", "#hat{s} [GeV] vs. p_{T}^{top} [GeV]",500,0,5000,500,0,2000);
  shat_vs_Pt_WAss = book<TH2F>( "shat_vs_Pt_WAss", "#hat{s} [GeV] vs. p_{T}^{WAss} [GeV]",500,0,5000,500,0,2000);
  Pt_top_vs_Pt_WAss = book<TH2F>( "Pt_top_vs_Pt_WAss", "p_{T}^{top} [GeV] vs. p_{T}^{WAss} [GeV] (all decay channels)",2000,0,2000,2000,0,2000);
  Pt_top_vs_Pt_LeptAss = book<TH2F>( "Pt_top_vs_Pt_LeptAss", "p_{T}^{top} [GeV] vs. p_{T}^{lep} [GeV] (all semi-leptonic channels with hadronic top)",2000,0,2000,2000,0,2000);

  Pt_WAss_over_Pt_top = book<TH1F>( "Pt_WAss_over_Pt_top", "p_{T}^{WAss}/p_{T}^{top}", 1000, 0, 10 ) ;

  deltaR_bAss_top = book<TH1F>( "deltaR_bAss_top", "#DeltaR(bAss, top)",1000,0,5);  
  deltaR_bAss_top_decays_min = book<TH1F>( "deltaR_bAss_top_decays_min", "#DeltaR_{min}(bAss, top decay prod.)",1000,0,5);  
  deltaR_bAss_top_vs_deltaR_top = book<TH2F>( "deltaR_bAss_top_vs_deltaR_top", "#DeltaR(bAss, top) vs #DeltaR_{max}(top decay prod.)",500,0,5,500,0,5);
  deltaR_bAss_top_decays_min_vs_deltaR_top = book<TH2F>( "deltaR_bAss_top_decays_min_vs_deltaR_top", "#DeltaR_{min}(bAss, top decay prod.) vs #DeltaR_{max}(top decay prod.)",500,0,5,500,0,5);
  deltaR_bAss_WAss_decays_min = book<TH1F>( "deltaR_bAss_WAss_decays_min", "#DeltaR_{min}(bAss, WAss decay prod.)",1000,0,5);  

  deltaR_gluonAss_top = book<TH1F>( "deltaR_gluonAss_top", "#DeltaR(gluonAss, top)",1000,0,5);  
  deltaR_gluonAss_top_decays_min = book<TH1F>( "deltaR_gluonAss_top_decays_min", "#DeltaR_{min}(gluonAss, top decay prod.)",1000,0,5);  
  deltaR_gluonAss_top_vs_deltaR_top = book<TH2F>( "deltaR_gluonAss_top_vs_deltaR_top", "#DeltaR(gluonAss, top) vs #DeltaR_{max}(top decay prod.)",500,0,5,500,0,5);
  deltaR_gluonAss_top_decays_min_vs_deltaR_top = book<TH2F>( "deltaR_gluonAss_top_decays_min_vs_deltaR_top", "#DeltaR_{min}(gluonAss, top decay prod.) vs #DeltaR_{max}(top decay prod.)",500,0,5,500,0,5);
  deltaR_gluonAss_WAss_decays_min = book<TH1F>( "deltaR_gluonAss_WAss_decays_min", "#DeltaR_{min}(gluonAss, WAss decay prod.)",1000,0,5);  

  // Electron Channel

  EleChannel_lep_pt = book<TH1F>( "EleChannel_lep_pt", "p_{T}^{lep} [GeV] (e-channel)", 2000, 0, 2000 );
  EleChannel_lep_eta = book<TH1F>( "EleChannel_lep_eta", "#eta_{lep} (e-channel)", 2000, -5, 5 ) ;
  EleChannel_top_pt = book<TH1F>( "EleChannel_top_pt", "p_{T}^{top} [GeV] (e-channel)", 2000, 0, 2000 );
  EleChannel_top_eta = book<TH1F>( "EleChannel_top_eta", "#eta_{top} (e-channel)", 2000, -5, 5 ) ;
  EleChannel_pt_topVSlep = book<TH2F>( "EleChannel_pt_topVSlep", "p_{T}^{top} [GeV] vs. p_{T}^{lep} [GeV] (e-channel)",2000,0,2000,2000,0,2000);

  EleChannel_lep_pt_topPt0to100 = book<TH1F>( "EleChannel_lep_pt_topPt0to100", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [0, 100) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt100to200 = book<TH1F>( "EleChannel_lep_pt_topPt100to200", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [100, 200) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt200to300 = book<TH1F>( "EleChannel_lep_pt_topPt200to300", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [200, 300) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt300to400 = book<TH1F>( "EleChannel_lep_pt_topPt300to400", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [300, 400) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt400to500 = book<TH1F>( "EleChannel_lep_pt_topPt400to500", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [400, 500) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt500to600 = book<TH1F>( "EleChannel_lep_pt_topPt500to600", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [500, 600) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt600to700 = book<TH1F>( "EleChannel_lep_pt_topPt600to700", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [600, 700) GeV)", 2000, 0, 2000 );
  EleChannel_lep_pt_topPt700toInf = book<TH1F>( "EleChannel_lep_pt_topPt700toInf", "p_{T}^{lep} [GeV] (e-channel, p_{T}^{top} #in [700, #infty) GeV)", 2000, 0, 2000 );

  EleChannel_nu_pt = book<TH1F>( "EleChannel_nu_pt", "p_{T}^{#nu} [GeV] (e-channel)", 2000, 0, 2000 );
  EleChannel_pt_topVSnu = book<TH2F>( "EleChannel_pt_topVSnu", "p_{T}^{top} [GeV] vs. p_{T}^{#nu} [GeV] (e-channel)",2000,0,2000,2000,0,2000);

  // Muon Channel

  MuoChannel_lep_pt = book<TH1F>( "MuoChannel_lep_pt", "p_{T}^{lep} [GeV] (#mu-channel)", 2000, 0, 2000 );
  MuoChannel_lep_eta = book<TH1F>( "MuoChannel_lep_eta", "#eta_{lep} (#mu-channel)", 2000, -5, 5 ) ;
  MuoChannel_top_pt = book<TH1F>( "MuoChannel_top_pt", "p_{T}^{top} [GeV] (#mu-channel)", 2000, 0, 2000 );
  MuoChannel_top_eta = book<TH1F>( "MuoChannel_top_eta", "#eta_{top} (#mu-channel)", 2000, -5, 5 ) ;
  MuoChannel_pt_topVSlep = book<TH2F>( "MuoChannel_pt_topVSlep", "p_{T}^{top} [GeV] vs. p_{T}^{lep} [GeV] (#mu-channel)",2000,0,2000,2000,0,2000);

  MuoChannel_lep_pt_topPt0to100 = book<TH1F>( "MuoChannel_lep_pt_topPt0to100", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [0, 100) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt100to200 = book<TH1F>( "MuoChannel_lep_pt_topPt100to200", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [100, 200) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt200to300 = book<TH1F>( "MuoChannel_lep_pt_topPt200to300", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [200, 300) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt300to400 = book<TH1F>( "MuoChannel_lep_pt_topPt300to400", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [300, 400) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt400to500 = book<TH1F>( "MuoChannel_lep_pt_topPt400to500", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [400, 500) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt500to600 = book<TH1F>( "MuoChannel_lep_pt_topPt500to600", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [500, 600) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt600to700 = book<TH1F>( "MuoChannel_lep_pt_topPt600to700", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [600, 700) GeV)", 2000, 0, 2000 );
  MuoChannel_lep_pt_topPt700toInf = book<TH1F>( "MuoChannel_lep_pt_topPt700toInf", "p_{T}^{lep} [GeV] (#mu-channel, p_{T}^{top} #in [700, #infty) GeV)", 2000, 0, 2000 );

  MuoChannel_nu_pt = book<TH1F>( "MuoChannel_nu_pt", "p_{T}^{#nu} [GeV] (#mu-channel)", 2000, 0, 2000 );
  MuoChannel_pt_topVSnu = book<TH2F>( "MuoChannel_pt_topVSnu", "p_{T}^{top} [GeV] vs. p_{T}^{#nu} [GeV] (#mu-channel)",2000,0,2000,2000,0,2000);

  h_singletopgen_twch = ctx.get_handle<SingleTopGen_tWch>("singletopgen_twch");
    
}


void SingleTopGen_tWchHists::fill(const uhh2::Event & e){
  //do not fill histograms if singletopgen information has not been filled
  if(!e.is_valid(h_singletopgen_twch)){
    return;
  }
  const auto & singletopgen_twch = e.get(h_singletopgen_twch);
    
  LorentzVector top = singletopgen_twch.Top().v4();
  LorentzVector bTop = singletopgen_twch.bTop().v4();
  LorentzVector WTop = singletopgen_twch.WTop().v4();
  LorentzVector WTopDecay1 = singletopgen_twch.WTopDecay1().v4();
  LorentzVector WTopDecay2 = singletopgen_twch.WTopDecay2().v4();
  LorentzVector WAss = singletopgen_twch.WAss().v4();
  LorentzVector WAssDecay1 = singletopgen_twch.WAssDecay1().v4();
  LorentzVector WAssDecay2 = singletopgen_twch.WAssDecay2().v4();
  LorentzVector LeptAss, NuAss;
  if(singletopgen_twch.IsAssLeptonicDecay())
    {
      LeptAss = singletopgen_twch.LeptAss().v4();
      NuAss = singletopgen_twch.NuAss().v4();
    }
  LorentzVector bAss = singletopgen_twch.bAss().v4();
  LorentzVector gluonAss = singletopgen_twch.gluonAss().v4();

  double sh = (e.genparticles->at(0).v4()+ e.genparticles->at(1).v4()).M();


  // General Event Properties

  shat->Fill(sh, e.weight);
  DecayChannel->Fill(singletopgen_twch.DecayChannel(), e.weight);
 
  // GenParticles

  Pt_top->Fill( top.Pt(), e.weight); 
  eta_top->Fill( top.eta(), e.weight);
  y_top->Fill( top.Rapidity(), e.weight);
  M_top->Fill( top.M(), e.weight);

  Pt_bTop->Fill( bTop.pt(), e.weight);
  eta_bTop->Fill( bTop.eta(), e.weight);
  y_bTop->Fill( bTop.Rapidity(), e.weight);
  M_bTop->Fill( bTop.M(), e.weight);

  Pt_WTop->Fill( WTop.pt(), e.weight);
  eta_WTop->Fill( WTop.eta(), e.weight);
  y_WTop->Fill( WTop.Rapidity(), e.weight);
  if(singletopgen_twch.WTop().v4().isTimelike())
    M_WTop->Fill( WTop.M(), e.weight);
    
  Pt_WAss->Fill( WAss.pt(), e.weight);
  eta_WAss->Fill( WAss.eta(), e.weight);
  y_WAss->Fill( WAss.Rapidity(), e.weight);
  if(singletopgen_twch.WAss().v4().isTimelike())
    M_WAss->Fill( WAss.M(), e.weight);

  Pt_WTopDecay1->Fill( WTopDecay1.pt(), e.weight);
  eta_WTopDecay1->Fill( WTopDecay1.eta(), e.weight);
  y_WTopDecay1->Fill( WTopDecay1.Rapidity(), e.weight);
  M_WTopDecay1->Fill( WTopDecay1.M(), e.weight);

  Pt_WTopDecay2->Fill( WTopDecay2.pt(), e.weight);
  eta_WTopDecay2->Fill( WTopDecay2.eta(), e.weight);
  y_WTopDecay2->Fill( WTopDecay2.Rapidity(), e.weight);
  M_WTopDecay2->Fill( WTopDecay2.M(), e.weight);

  Pt_WAssDecay1->Fill( WAssDecay1.pt(), e.weight);
  eta_WAssDecay1->Fill( WAssDecay1.eta(), e.weight);
  y_WAssDecay1->Fill( WAssDecay1.Rapidity(), e.weight);
  M_WAssDecay1->Fill( WAssDecay1.M(), e.weight);

  Pt_WAssDecay2->Fill( WAssDecay2.pt(), e.weight);
  eta_WAssDecay2->Fill( WAssDecay2.eta(), e.weight);
  y_WAssDecay2->Fill( WAssDecay2.Rapidity(), e.weight);
  M_WAssDecay2->Fill( WAssDecay2.M(), e.weight);

  if(singletopgen_twch.HasAssociatedBottom()){
    Pt_bAss->Fill( bAss.pt(), e.weight);
    eta_bAss->Fill( bAss.eta(), e.weight);
    y_bAss->Fill( bAss.Rapidity(), e.weight);
    M_bAss->Fill( bAss.M(), e.weight);
  }
    
  if(singletopgen_twch.HasAssociatedGluon()){
    Pt_gluonAss->Fill( gluonAss.pt(), e.weight);
    eta_gluonAss->Fill( gluonAss.eta(), e.weight);
    y_gluonAss->Fill( gluonAss.Rapidity(), e.weight);
    M_gluonAss->Fill( gluonAss.M(), e.weight);
  }

  // more

  Pt_top_over_shat->Fill( top.Pt()/sh, e.weight);
  Pt_bTop_over_shat->Fill( bTop.Pt()/sh, e.weight);
  Pt_WTop_over_shat->Fill( WTop.Pt()/sh, e.weight);
  Pt_WAss_over_shat->Fill( WAss.Pt()/sh, e.weight);
  Pt_WTopDecay1_over_shat->Fill( WTopDecay1.Pt()/sh, e.weight);
  Pt_WTopDecay2_over_shat->Fill( WTopDecay2.Pt()/sh, e.weight);
  Pt_WAssDecay1_over_shat->Fill( WAssDecay1.Pt()/sh, e.weight);
  Pt_WAssDecay2_over_shat->Fill( WAssDecay2.Pt()/sh, e.weight);


  double deltaR_top = std::max(
			       std::max(
					uhh2::deltaR(singletopgen_twch.bTop(), singletopgen_twch.WTopDecay1()), uhh2::deltaR(singletopgen_twch.bTop(), singletopgen_twch.WTopDecay2())
					),
			       uhh2::deltaR(singletopgen_twch.WTopDecay1(), singletopgen_twch.WTopDecay2())
			       );
  deltaR_top_decays->Fill(deltaR_top, e.weight);
  shat_vs_deltaR_top->Fill(sh, deltaR_top, e.weight);
  DeltaR_top_vs_Pt_top->Fill(top.Pt(), deltaR_top, e.weight);

  double deltaR_WTop = uhh2::deltaR(singletopgen_twch.WTopDecay1(), singletopgen_twch.WTopDecay2());
  deltaR_WTop_decays->Fill(deltaR_WTop, e.weight);
  shat_vs_deltaR_WTop->Fill(sh, deltaR_WTop, e.weight);
  Pt_WTop_vs_deltaR_WTop->Fill(top.Pt(), deltaR_WTop, e.weight);

  double deltaR_WAss = uhh2::deltaR(singletopgen_twch.WAssDecay1(), singletopgen_twch.WAssDecay2());
  deltaR_WAss_decays->Fill(deltaR_WAss, e.weight);
  shat_vs_deltaR_WAss->Fill(sh, deltaR_WAss, e.weight);
  Pt_WAss_vs_deltaR_WAss->Fill(top.Pt(), deltaR_WAss, e.weight);


  shat_vs_Pt_top->Fill(sh, top.Pt(), e.weight);
  shat_vs_Pt_WAss->Fill(sh, WAss.Pt(), e.weight);
  Pt_top_vs_Pt_WAss->Fill( WAss.Pt(), top.Pt(), e.weight);
  if(singletopgen_twch.IsAssLeptonicDecay()){
    Pt_top_vs_Pt_LeptAss->Fill( LeptAss.Pt(), top.Pt(), e.weight);
  }

  Pt_WAss_over_Pt_top->Fill( WAss.Pt()/top.Pt(), e.weight);

  if(singletopgen_twch.HasAssociatedBottom()){
    double deltaR_bAss_top_variable = uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.Top());
    double deltaR_bAss_top_decays = std::min(
					     std::min(
						      uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.bTop()), uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.WTopDecay1())
						      ),
					     uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.WTopDecay2())
					     );
    double deltaR_bAss_WAss_decays = std::min(
					      uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.WAssDecay1()), uhh2::deltaR(singletopgen_twch.bAss(), singletopgen_twch.WAssDecay2())
					      );
    deltaR_bAss_top->Fill(deltaR_bAss_top_variable, e.weight);
    deltaR_bAss_top_decays_min->Fill(deltaR_bAss_top_decays, e.weight);
    deltaR_bAss_top_vs_deltaR_top->Fill(deltaR_bAss_top_variable, deltaR_top, e.weight);
    deltaR_bAss_top_decays_min_vs_deltaR_top->Fill(deltaR_bAss_top_decays, deltaR_top, e.weight);
    deltaR_bAss_WAss_decays_min->Fill(deltaR_bAss_WAss_decays, e.weight);
    
  }
  if(singletopgen_twch.HasAssociatedGluon()){
    double deltaR_gluonAss_top_variable = uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.Top());
    double deltaR_gluonAss_top_decays = std::min(
						 std::min(
							  uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.bTop()), uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.WTopDecay1())
							  ),
						 uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.WTopDecay2())
						 );
    double deltaR_gluonAss_WAss_decays = std::min(
					      uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.WAssDecay1()), uhh2::deltaR(singletopgen_twch.gluonAss(), singletopgen_twch.WAssDecay2())
					      );
    deltaR_gluonAss_top->Fill(deltaR_gluonAss_top_variable, e.weight);
    deltaR_gluonAss_top_decays_min->Fill(deltaR_gluonAss_top_decays, e.weight);
    deltaR_gluonAss_top_vs_deltaR_top->Fill(deltaR_gluonAss_top_variable, deltaR_top, e.weight);
    deltaR_gluonAss_top_decays_min_vs_deltaR_top->Fill(deltaR_gluonAss_top_decays, deltaR_top, e.weight);
    deltaR_gluonAss_WAss_decays_min->Fill(deltaR_gluonAss_WAss_decays, e.weight);
  }
  
  if(singletopgen_twch.IsTopHadronicDecay() && singletopgen_twch.IsAssToElectronDecay()){
    EleChannel_lep_pt->Fill(LeptAss.Pt(), e.weight);
    EleChannel_lep_eta->Fill(LeptAss.Eta(), e.weight);
    EleChannel_top_pt->Fill(top.Pt(), e.weight);
    EleChannel_top_eta->Fill(top.Eta(), e.weight);
    EleChannel_pt_topVSlep->Fill(LeptAss.Pt(), top.Pt(), e.weight);
    if(top.Pt() >= 0 && top.Pt() < 100) EleChannel_lep_pt_topPt0to100->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 100 && top.Pt() < 200) EleChannel_lep_pt_topPt100to200->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 200 && top.Pt() < 300) EleChannel_lep_pt_topPt200to300->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 300 && top.Pt() < 400) EleChannel_lep_pt_topPt300to400->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 400 && top.Pt() < 500) EleChannel_lep_pt_topPt400to500->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 500 && top.Pt() < 600) EleChannel_lep_pt_topPt500to600->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 600 && top.Pt() < 700) EleChannel_lep_pt_topPt600to700->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() > 700) EleChannel_lep_pt_topPt700toInf->Fill(LeptAss.Pt(), e.weight);

    EleChannel_nu_pt->Fill(NuAss.Pt(), e.weight);
    EleChannel_pt_topVSnu->Fill(NuAss.Pt(), top.Pt(), e.weight);
  }

  if(singletopgen_twch.IsTopHadronicDecay() && singletopgen_twch.IsAssToMuonDecay()){
    MuoChannel_lep_pt->Fill(LeptAss.Pt(), e.weight);
    MuoChannel_lep_eta->Fill(LeptAss.Eta(), e.weight);
    MuoChannel_top_pt->Fill(top.Pt(), e.weight);
    MuoChannel_top_eta->Fill(top.Eta(), e.weight);
    MuoChannel_pt_topVSlep->Fill(LeptAss.Pt(), top.Pt(), e.weight);
    if(top.Pt() >= 0 && top.Pt() < 100) MuoChannel_lep_pt_topPt0to100->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 100 && top.Pt() < 200) MuoChannel_lep_pt_topPt100to200->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 200 && top.Pt() < 300) MuoChannel_lep_pt_topPt200to300->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 300 && top.Pt() < 400) MuoChannel_lep_pt_topPt300to400->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 400 && top.Pt() < 500) MuoChannel_lep_pt_topPt400to500->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 500 && top.Pt() < 600) MuoChannel_lep_pt_topPt500to600->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 600 && top.Pt() < 700) MuoChannel_lep_pt_topPt600to700->Fill(LeptAss.Pt(), e.weight);
    else if(top.Pt() >= 700) MuoChannel_lep_pt_topPt700toInf->Fill(LeptAss.Pt(), e.weight);

    MuoChannel_nu_pt->Fill(NuAss.Pt(), e.weight);
    MuoChannel_pt_topVSnu->Fill(NuAss.Pt(), top.Pt(), e.weight);
  }

}
