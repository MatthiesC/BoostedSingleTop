#include "UHH2/BoostedSingleTop/include/MatchHists.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWch.h"
#include "UHH2/BoostedSingleTop/include/TTbarGen_.h"

#include "UHH2/common/include/Utils.h"

#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
#include <vector>

using namespace std;
using namespace uhh2;

// Use this Histogram Class only for single-top and ttbar MC samples after all selections are applied (most important: N_topjet >= 1)!
// Also make sure that also the respective GenProducers are applied beforehand!

MatchHists::MatchHists(Context & ctx, const string & dirname, const boost::optional<TopJetId> & topjetid):
  Hists(ctx,dirname), m_topjetid(topjetid)
{
  std::string dataset_name = ctx.get("dataset_version");
  is_tw = (dataset_name.find("SingleTop_T_tWch") == 0 || dataset_name.find("SingleTop_Tbar_tWch") == 0);
  is_tt = (dataset_name.find("TTbar") == 0);
  //is_tt = true;
  //is_tw = true;

  if(is_tw)
    {
      h_singletopgen_twch = ctx.get_handle<SingleTopGen_tWch>("h_singletopgen_twch");
    }
  else if(is_tt)
    {
      h_ttbargen = ctx.get_handle<TTbarGen_>("h_ttbargen");
    }
  
  // Book all histograms here ...
  H_shat = book<TH1F>("shat", "#hat{s} [GeV]", 50, 0, 5000);
  H_DecayChannel = book<TH1F>("DecayChannel", "decay channel", 17, -0.5, 16.5);

  H_topgen_pt = book<TH1F>("topgen_pt", "top quark p_{T}^{gen} [GeV]", 16, 0, 1600);
  H_lepgen_pt = book<TH1F>("lepgen_pt", "lepton p_{T}^{gen} [GeV]", 16, 0, 1600);

  H_topjet_pt = book<TH1F>("topjet_pt", "top quark p_{T}^{reco} [GeV]", 16, 0, 1600);
  H_leprec_pt = book<TH1F>("leprec_pt", "lepton p_{T}^{reco} [GeV]", 16, 0, 1600);

  H_topgen_eta = book<TH1F>("topgen_eta", "top quark #eta^{gen} [GeV]", 20, -4, 4);
  H_lepgen_eta = book<TH1F>("lepgen_eta", "lepton #eta^{gen} [GeV]", 20, -4, 4);

  H_topjet_eta = book<TH1F>("topjet_eta", "top quark #eta^{reco} [GeV]", 20, -4, 4);
  H_leprec_eta = book<TH1F>("leprec_eta", "lepton #eta^{reco} [GeV]", 20, -4, 4);

  H_match_deltaR_lep = book<TH1F>("match_deltaR_lep", "#DeltaR(l^{reco}, l^{gen})", 10, 0, 5);
  H_match_deltaR_top = book<TH1F>("match_deltaR_top", "#DeltaR(t^{reco}, t^{gen})", 10, 0, 5);
  H_match_deltaPhi_lep = book<TH1F>("match_Phi_lep", "#Delta#phi(l^{reco}, l^{gen})", 10, -M_PI, M_PI);
  H_match_deltaPhi_top = book<TH1F>("match_Phi_top", "#Delta#phi(t^{reco}, t^{gen})", 10, -M_PI, M_PI);
  H_match_deltaPhi_neu = book<TH1F>("match_Phi_neu", "#Delta#phi(#nu^{reco}, #nu^{gen})", 10, -M_PI, M_PI);

}

void MatchHists::fill(const Event & e)
{
  double w = e.weight;

  // Get all necessary event information on reco level ...
  vector<TopJet> topjets = *e.topjets;
  vector<Muon> muons = *e.muons;
  sort_by_pt<Muon>(muons);
  vector<Electron> electrons = *e.electrons;
  sort_by_pt<Electron>(electrons);
  vector<Jet> jets = *e.jets;
  sort_by_pt<Jet>(jets);
  vector<Jet> bjets;
  vector<TopJet> toptaggedjets;
  for (TopJet topjet : topjets)
    {
      if (m_topjetid && (*m_topjetid)(topjet, e))
	{
	  toptaggedjets.push_back(topjet);
	}
    }
  sort_by_pt<TopJet>(toptaggedjets);
  MET met = *e.met;
  const CSVBTag btag_tight(CSVBTag::WP_TIGHT);
  int n_btags = 0;
  for (Jet jet : jets)
    {
      if (btag_tight(jet, e))
	{
	  bjets.push_back(jet);
	  ++n_btags;
	}
    }
  sort_by_pt<Jet>(bjets);

  // Check whether we actually have the signal signature ...
  if(!(toptaggedjets.size() == 1 && (muons.size() + electrons.size() == 1) && n_btags > 0)) return;

  // Identify some important particles on reco level ...
  TopJet hotvrjet = toptaggedjets.at(0);
  FlavorParticle primlept;
  if(electrons.size() > 0) primlept = electrons.at(0);
  else if(muons.size() > 0) primlept = muons.at(0);

  // Sort b-jets by deltaR to top-tagged jet ...
  /*vector<double> deltaR_bjet_topjet;
  for(Jet bjet : bjets)
    {
      deltaR_bjet_topjet.push_back(deltaR(bjet,hotvrjet));
    }
  vector<int> index(bjets.size(), 0);
  for(unsigned int i = 0; i < bjets.size(); ++i)
    {
      index.at(i) = i;
    }
  sort(index.begin(), index.end(),
       [&](const int & a, const int & b) {
	 return(deltaR_bjet_topjet.at(a) > deltaR_bjet_topjet.at(b));
       }
       );
  vector<Jet> bjets_sorted; // new vector containing all b-jets sorted by deltaR (descending) to hotvrjet
  bjets_sorted.reserve(bjets.size());
  if(n_btags > 0)
    {
      for(unsigned int i = 0; i < bjets.size(); ++i)
	{
	  cout << bjets_sorted[i].v4().pt() << endl;
	  cout << index[i] << endl;
	  bjets_sorted[i] = bjets[index[i]];
	}
	}*/

  // Get all necessary event information on gen level ...
  LorentzVector top, lepton, neutrino;
  vector<LorentzVector> bottoms;
  vector<LorentzVector> topdaughters;
  double shat = (e.genparticles->at(0).v4() + e.genparticles->at(1).v4()).M();
  int DecayChannel = 0;

  bool is_signal;

  if(is_tw)
    {
      if(!e.is_valid(h_singletopgen_twch))
	{
	  //return; // Do not fill histograms if gen_tw information has not been filled
	}
      const auto & gen_tw = e.get(h_singletopgen_twch);
      DecayChannel = gen_tw.DecayChannel();
      is_signal = (gen_tw.IsTopHadronicDecay() && (gen_tw.IsAssToElectronDecay() || gen_tw.IsAssToMuonDecay()));
      if(is_signal)
	{
	  top = gen_tw.Top().v4();
	  lepton = gen_tw.LeptAss().v4();
	  neutrino = gen_tw.NuAss().v4();
	  topdaughters.push_back(gen_tw.bTop().v4());
	  topdaughters.push_back(gen_tw.WTopDecay1().v4());
	  topdaughters.push_back(gen_tw.WTopDecay2().v4());
	  for(GenParticle gp : *e.genparticles)
	    {
	      if(abs(gp.pdgId()) == 5)
		{
		  bottoms.push_back(gp.v4());
		}
	    }
	}
    }
  else if(is_tt)
    {
      if(!e.is_valid(h_ttbargen))
	{
	  //return; // Do not fill histograms if gen_tt information has not been filled
	}
      const auto & gen_tt = e.get(h_ttbargen);
      DecayChannel = gen_tt.DecayChannel();
      is_signal = (gen_tt.IsSemiLeptonicDecay() && (gen_tt.ChargedLepton().pdgId() == 11 || gen_tt.ChargedLepton().pdgId() == 13));
      if(is_signal)
	{
	  if(gen_tt.IsTopHadronicDecay())
	    {
	      top = gen_tt.Top().v4();
	      topdaughters.push_back(gen_tt.bTop().v4());
	      topdaughters.push_back(gen_tt.Wdecay1().v4());
	      topdaughters.push_back(gen_tt.Wdecay2().v4());
	    }
	  else if(gen_tt.IsAntiTopHadronicDecay())
	    {
	      top = gen_tt.Antitop().v4();
	      topdaughters.push_back(gen_tt.bAntitop().v4());
	      topdaughters.push_back(gen_tt.WMinusdecay1().v4());
	      topdaughters.push_back(gen_tt.WMinusdecay2().v4());
	    }
	  lepton = gen_tt.ChargedLepton().v4();
	  neutrino = gen_tt.Neutrino().v4();
	  for(GenParticle gp : *e.genparticles)
	    {
	      if(abs(gp.pdgId()) == 5)
		{
		  bottoms.push_back(gp.v4());
		}
	    }
	}
    }

  // Now fill the histograms ...
  H_shat->Fill(shat, w);
  H_DecayChannel->Fill(DecayChannel, w);

  H_topgen_pt->Fill(top.Pt(), w);
  H_lepgen_pt->Fill(lepton.Pt(), w);

  H_topjet_pt->Fill(hotvrjet.v4().Pt(), w);
  H_leprec_pt->Fill(primlept.v4().Pt(), w);

  H_topgen_eta->Fill(top.Eta(), w);
  H_lepgen_eta->Fill(lepton.Eta(), w);

  H_topjet_eta->Fill(hotvrjet.v4().Eta(), w);
  H_leprec_eta->Fill(primlept.v4().Eta(), w);

  H_match_deltaR_lep->Fill(deltaR(lepton, primlept), w);
  H_match_deltaR_top->Fill(deltaR(top, hotvrjet), w);
  H_match_deltaPhi_lep->Fill(deltaPhi(lepton, primlept), w);
  H_match_deltaPhi_top->Fill(deltaPhi(top, hotvrjet), w);
  H_match_deltaPhi_neu->Fill(deltaPhi(met, neutrino), w);
}

MatchHists::~MatchHists(){}
