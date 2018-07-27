#include "UHH2/BoostedSingleTop/include/TTbarGen_.h"

using namespace std;
using namespace uhh2;

TTbarGen_::TTbarGen_(const vector<GenParticle> & genparticles, bool throw_on_failure): m_type(e_notfound) {    
    int n_top = 0, n_antitop = 0;
    GenParticle top;
    GenParticle antitop;
    for(unsigned int i=0; i<genparticles.size(); ++i) {
        const GenParticle & genp = genparticles[i];
        if (abs(genp.pdgId()) == 6){
	    if(top.pdgId() !=6 && genp.pdgId()==6) top = genp;
	    if(antitop.pdgId() !=-6 && genp.pdgId()==-6) antitop = genp;
            auto w = genp.daughter(&genparticles, 1);
            auto b = genp.daughter(&genparticles, 2);
            if(!w || !b){
	      if(throw_on_failure) throw runtime_error("TTbarGen: top has not ==2 daughters");
              return;
            }
            if(abs(w->pdgId()) != 24){
	      std::swap(w, b);
            }
	    /* It rarely happens that the list of genparticles contains 4 or more (instead of 2) particles which reckon the same top
	       as their mother although each particle including the tops can just have two daughters. E.g. if the top emits a photon
	       before decaying into b and W, this photon may split up into two leptons which reckon the top as their mother, too.
	       Therefore, it may happen that those leptons are considered as the top daughters whereas b and W are "ignored" and cannot
	       be found. This workaround fixes that issue: */
	    if(abs(w->pdgId()) != 24) {
	      for(unsigned int j = 0; j < genparticles.size(); ++j) {
		const GenParticle & gp = genparticles[j];
		auto m1 = gp.mother(&genparticles, 1);
		auto m2 = gp.mother(&genparticles, 2);
		bool has_top_mother = ((m1 && m1->index() == genp.index()) || (m2 && m2->index() == genp.index()));
		if(has_top_mother && (abs(gp.pdgId()) == 24)) {
		  w = &gp;
		  break;
		}
	      }
	    }
            if(abs(w->pdgId()) != 24){
	        if(throw_on_failure) throw runtime_error("TTbarGen: top has no W daughter");
                return;
            }
            
            // NOTE: here, we could skip over intermediate W bosons. However,
            // this Pythia8-related problem is now fixed when creating ntuples already,
            // so this should not be necessary.
            
	    /* Do a similar workaround as above if the expected b daughter has not been found yet */
	    if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3 && abs(b->pdgId()) != 1) {
	      for(unsigned int j = 0; j < genparticles.size(); ++j) {
		const GenParticle & gp = genparticles[j];
		auto m1 = gp.mother(&genparticles, 1);
		auto m2 = gp.mother(&genparticles, 2);
		bool has_top_mother = ((m1 && m1->index() == genp.index()) || (m2 && m2->index() == genp.index()));
		if(has_top_mother && (abs(gp.pdgId()) == 5 || abs(gp.pdgId()) == 3 || abs(gp.pdgId()) == 1)) {
		  b = &gp;
		  break;
		}
	      }
	    }
            if(abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3   && abs(b->pdgId()) != 1){
	      if(throw_on_failure) throw runtime_error("TTbarGen: top has no b daughter");
	      return;
            }
            // now get W daughters:

	    int n_wdau=0;

	    auto wd1 = w->daughter(&genparticles, 1);
	    auto wd2 = w->daughter(&genparticles, 2);

	    while (n_wdau!=2){

	      if(wd1 && !wd2){
		w = wd1;
		wd1 = w->daughter(&genparticles, 1);
		wd2 = w->daughter(&genparticles, 2);
	      }
	      else if(wd1 && wd2){
		n_wdau=2;
	      }

	      else{
		if(throw_on_failure) throw runtime_error("TTbarGen: W has no daughters");
		return;
	      }

	    }
	    if(!wd1 || !wd2){
	      if(throw_on_failure) throw runtime_error("TTbarGen: W has not ==2 daughters");
	      return;
            }
            
            // now that we collected everything, fill the member variables. 
            // Use different member variables according to top quark charge.
            if(genp.pdgId() == 6){
                m_Top = top;
                m_WTop = *w;
                m_bTop = *b;
                m_Wdecay1 = *wd1;
                m_Wdecay2 = *wd2;
                ++n_top;
            }
            else{
                m_Antitop = antitop;
                m_WAntitop = *w;
                m_bAntitop = *b;
                m_WMinusdecay1 = *wd1;
                m_WMinusdecay2 = *wd2;
                ++n_antitop;
            }
        }
    }

    if(n_top != 1 || n_antitop != 1){
        if(throw_on_failure)  throw runtime_error("TTbarGen: did not find exactly one top and one antitop in the event");
        return;
    }
    
    // calculate decay channel by counting the number of charged leptons
    // in the W daughters:
    int n_e = 0, n_m = 0, n_t = 0;
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        int id = abs(wd.pdgId());
        if(id == 11) ++n_e;
        else if(id == 13) ++n_m;
        else if(id == 15) ++n_t;
    }

    // dilepton channels:
    if(n_e == 2){
        m_type = e_ee;
    }
    else if(n_e == 1 && n_m == 1){
        m_type = e_emu;
    }
    else if(n_e == 1 && n_t == 1){
        m_type = e_etau;
    }
    else if(n_m == 2){
        m_type = e_mumu;
    }
    else if(n_m == 1 && n_t == 1){
        m_type = e_mutau;
    }
    else if(n_t == 2){
        m_type = e_tautau;
    }
    // lepton+jet channels:
    else if(n_e == 1){
        m_type = e_ehad;
    }
    else if(n_m == 1){
        m_type = e_muhad;
    }
    else if(n_t == 1){
        m_type = e_tauhad;
    }
    // hadronic:
    else{
        m_type = e_had;
    }
}   
 

GenParticle TTbarGen_::Top() const{
    return m_Top;
}

GenParticle TTbarGen_::Antitop() const{
    return m_Antitop;
} 

GenParticle TTbarGen_::WTop() const{
    return m_WTop;
}

GenParticle TTbarGen_::WAntitop() const{
    return m_WAntitop;
}

GenParticle TTbarGen_::bTop() const{
    return m_bTop;
}

GenParticle TTbarGen_::bAntitop() const{
    return m_bAntitop;
} 

GenParticle TTbarGen_::Wdecay1() const{
    return m_Wdecay1;
} 

GenParticle TTbarGen_::Wdecay2() const{
    return m_Wdecay2;
} 

GenParticle TTbarGen_::WMinusdecay1() const{
    return m_WMinusdecay1;
} 

GenParticle TTbarGen_::WMinusdecay2() const{
    return m_WMinusdecay2;
} 

TTbarGen_::E_DecayChannel TTbarGen_::DecayChannel()  const{  
    return m_type;
}

bool TTbarGen_::IsTopHadronicDecay()  const{
    return abs(m_Wdecay1.pdgId()) <= 5;
}

bool TTbarGen_::IsAntiTopHadronicDecay()  const{
    return abs(m_WMinusdecay1.pdgId()) <= 5;
}

bool TTbarGen_::IsSemiLeptonicDecay() const{
   return m_type == e_ehad ||  m_type == e_muhad ||m_type == e_tauhad;  
}


namespace {
    
bool is_charged_lepton(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 11 || id == 13 || id == 15;
}

bool is_neutrino(const GenParticle & gp){
    int id = abs(gp.pdgId());
    return id == 12 || id == 14 || id == 16;
}

}

GenParticle TTbarGen_::ChargedLepton() const{
   if (m_type != e_ehad &&  m_type != e_muhad && m_type!= e_tauhad){
        throw runtime_error("TTbarGen_::ChargedLepton called, but this is no l+jets ttbar event!");
    }
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        if(is_charged_lepton(wd)) return wd;
    }
    throw logic_error("logic error in TTbarGen_::ChargedLepton");
}

GenParticle TTbarGen_::Neutrino() const{
   if (m_type != e_ehad &&  m_type != e_muhad && m_type!= e_tauhad){
        throw runtime_error("TTbarGen_::ChargedLepton called, but this is no l+jets ttbar event!");
    }
    for(const auto & wd : {m_Wdecay1, m_Wdecay2, m_WMinusdecay1, m_WMinusdecay2}){
        if(is_neutrino(wd)) return wd;
    }
    throw logic_error("logic error in TTbarGen_::Neutrino");
}

GenParticle TTbarGen_::TopLep() const{
    if(ChargedLepton().charge()>0) return Top();
    else return Antitop();
}

GenParticle TTbarGen_::TopHad() const{
    if(ChargedLepton().charge()<0) return Top();
    else return Antitop();
}

GenParticle TTbarGen_::BLep() const{
    if(ChargedLepton().charge()>0) return bTop();
    else return bAntitop();
}

GenParticle TTbarGen_::BHad() const{
    if(ChargedLepton().charge()<0) return bTop();
    else return bAntitop();
}

GenParticle TTbarGen_::WLep() const{
    if(ChargedLepton().charge()>0) return WTop();
    else return WAntitop();
}

GenParticle TTbarGen_::WHad() const{
    if(ChargedLepton().charge()<0) return WTop();
    else return WAntitop();
}

GenParticle TTbarGen_::Q1() const{
    if(ChargedLepton().charge()>0) return WMinusdecay1();
    else return Wdecay1();
}

GenParticle TTbarGen_::Q2() const{
    if(ChargedLepton().charge()>0) return WMinusdecay2();
    else return Wdecay2();
}


TTbarGen_Producer::TTbarGen_Producer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
  h_ttbargen = ctx.get_handle<TTbarGen_>(name);
}

bool TTbarGen_Producer::process(Event & event){
    event.set(h_ttbargen, TTbarGen_(*event.genparticles, throw_on_failure));
    return true;
}
