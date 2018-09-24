#include "UHH2/BoostedSingleTop/include/SingleTopGenHists.h"
#include "TH1F.h"
#include "TH2F.h"

using namespace uhh2;

SingleTopGenHists::SingleTopGenHists(uhh2::Context & ctx, const std::string & dirname): Hists(ctx, dirname){

    shat = book< TH1F>( "shat", "#hat{s} [GeV]", 500, 0, 5000 ) ;
    DecayChannel = book< TH1F>( "DecayChannel", "decay channel", 11, -0.5, 10.5 ) ;

    M_top = book< TH1F>("M_top", "top-quark mass [GeV]", 600, 50, 350) ;

    Pt_top = book< TH1F>( "Pt_top", "top-quark p_{T} [GeV]", 1000, 0, 2000 ) ;
    
    Pt_top_over_shat = book< TH1F>( "Pt_top_over_shat", "top-quark p_{T}/#hat{s}", 1000, 0, 1 ) ;
    
    eta_top = book< TH1F>( "eta_top", "top-quark #eta", 1000, -5, 5 ) ;
    y_top = book< TH1F>( "y_top", "top-quark rapidity", 1000, -5, 5 ) ;
    phi_top = book<TH1F> ( "phi_top", "top-quark #phi [rad]", 500, -M_PI, M_PI);
    top_eta_vs_phi = book<TH2F> ("", "top-quark #eta-#phi-plane", 500, -M_PI, M_PI, 1000, -5, 5 );
     
    deltaR_top_decays = book<TH1F>( "deltaR_top_decays", "max. #DeltaR between top-quark decay products",1000,0,5);
     
    shat_vs_deltaR_top = book<TH2F>( "shat_vs_deltaR_top", "#hat{s} [GeV] vs. #DeltaR(t decay prod.)",500,0,5000,500,0,5);
 
    Pt_top_vs_deltaR_top = book<TH2F>( "Pt_top_vs_deltaR_top", "P_{T,t} [GeV/c] vs. #DeltaR(t decay prod.)",500,0,2000,500,0,5);
    Pt_top_vs_deltaR_top_hadronic = book<TH2F>( "Pt_top_vs_deltaR_top_hadronic", "P_{T,t} [GeV/c] vs. #DeltaR(t decay prod.)",500,0,2000,500,0,5);
 
    deltaR_W_decays = book<TH1F>( "deltaR_W_decays", "max. #DeltaR between W-boson decay products",1000,0,5);

    shat_vs_Pt_top = book<TH2F>( "shat_vs_Pt_top", "#hat{s} [GeV] vs P_{T,t} [GeV/c]",500,0,5000,500,0,2000); 
    
    Pt_b = book<TH1F>( "Pt_b","P_{T,b} [GeV/c]",1000,0,1000);
    eta_b = book<TH1F>( "eta_b","#eta_{b}",1000,-5,5);
    y_b = book<TH1F>( "y_b","y_{b}",1000,-5,5);
    M_b = book<TH1F>( "M_b","M_{b} [GeV/c^{2}]",1000,0,10);

    Pt_W = book<TH1F>( "Pt_W","P_{T,W} [GeV/c]",1000,0,1000);
    eta_W = book<TH1F>( "eta_W","#eta_{W}",1000,-5,5);
    y_W = book<TH1F>( "y_W","y_{W}",1000,-5,5);
    M_W = book<TH1F>( "M_W","M_{W} [GeV/c^{2}]",1000,50,100);

    Pt_top_vs_Pt_W = book<TH2F>( "Pt_top_vs_Pt_W", "P_{T,t} [GeV] vs P_{T,W}",500,0,2000,500,0,2000);

    h_singletopgen = ctx.get_handle<SingleTopGen>("singletopgen");
    
}


void SingleTopGenHists::fill(const uhh2::Event & e){
    //do not fill histograms if singletopgen information has not been filled
    if(!e.is_valid(h_singletopgen)){
      return;
    }
    const auto & singletopgen = e.get(h_singletopgen);
    
    LorentzVector top = singletopgen.Top().v4();
 
    double sh = (e.genparticles->at(0).v4()+ e.genparticles->at(1).v4()).M();

    DecayChannel->Fill(singletopgen.DecayChannel(), e.weight);

    shat->Fill(sh, e.weight);
    
    M_top->Fill( top.M(), e.weight);
    Pt_top->Fill( top.Pt(), e.weight);

    Pt_top_over_shat->Fill( top.Pt()/sh, e.weight);

    eta_top->Fill( top.eta(), e.weight);
    phi_top->Fill( top.phi(), e.weight);
    top_eta_vs_phi->Fill(top.phi(), top.eta(), e.weight);
    y_top->Fill( top.Rapidity(), e.weight);

    double deltaR_topdecays = std::max (std::max( uhh2::deltaR(singletopgen.bTop(), singletopgen.Wdecay1() ), 
				  uhh2::deltaR(singletopgen.bTop(), singletopgen.Wdecay2() ) )
			     , uhh2::deltaR(singletopgen.Wdecay1(), singletopgen.Wdecay2() ) );


    deltaR_top_decays->Fill(deltaR_topdecays,e.weight);
    
    double deltaR_Wdecays = uhh2::deltaR(singletopgen.Wdecay1(), singletopgen.Wdecay2());

    deltaR_W_decays->Fill(deltaR_Wdecays,e.weight);

    Pt_top_vs_deltaR_top->Fill(top.Pt(), deltaR_topdecays, e.weight);
    if(singletopgen.IsTopHadronicDecay()) Pt_top_vs_deltaR_top_hadronic->Fill(top.Pt(), deltaR_topdecays, e.weight);

    shat_vs_deltaR_top->Fill(sh, deltaR_topdecays, e.weight);
    shat_vs_Pt_top->Fill(sh, top.Pt(),  e.weight);

    Pt_b->Fill( singletopgen.bTop().pt(), e.weight);
    eta_b->Fill( singletopgen.bTop().eta(), e.weight);
    y_b->Fill( singletopgen.bTop().v4().Rapidity(), e.weight);
    M_b->Fill( singletopgen.bTop().v4().M(), e.weight);

    Pt_W->Fill( singletopgen.WTop().pt(), e.weight);
    eta_W->Fill( singletopgen.WTop().eta(), e.weight);
    y_W->Fill( singletopgen.WTop().v4().Rapidity(), e.weight);
    if(singletopgen.WTop().v4().isTimelike())
      M_W->Fill( singletopgen.WTop().v4().M(), e.weight);

    Pt_top_vs_Pt_W->Fill( top.Pt(), singletopgen.WTop().pt(), e.weight);

}
