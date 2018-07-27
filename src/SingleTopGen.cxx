// based on TTbarGen.cxx
// edited by Christopher Matthies, 02/2018

#include "UHH2/BoostedSingleTop/include/SingleTopGen.h"
#include "UHH2/core/include/Utils.h"

using namespace std;
using namespace uhh2;

SingleTopGen::SingleTopGen(const vector<GenParticle> & genparticles, bool throw_on_failure): m_type(e_notfound) {    

  int n_top = 0;
  GenParticle top;

  for(unsigned int i = 0; i < genparticles.size(); ++i) {
    const GenParticle & genp = genparticles[i];
    if (abs(genp.pdgId()) == 6){ // We don't distinguish between top and antitop. May be changed later but for now it's okay
      top = genp;
      auto w = genp.daughter(&genparticles, 1);
      auto b = genp.daughter(&genparticles, 2);
      if(!w || !b) {
        if(throw_on_failure) throw runtime_error("SingleTopGen: top has not ==2 daughters");
        continue;
      }
      if (abs(w->pdgId()) != 24) {
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
      if (abs(w->pdgId()) != 24) {
        if(throw_on_failure) throw runtime_error("SingleTopGen: top has no W daughter");
        continue;
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
      if (abs(b->pdgId()) != 5 && abs(b->pdgId()) != 3   && abs(b->pdgId()) != 1) {
        if(throw_on_failure) throw runtime_error("SingleTopGen: top has no b daughter");
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
          if(throw_on_failure) throw runtime_error("SingleTopGen: W has no daughters");
          return;
        }

      }
      if(!wd1 || !wd2){
        if(throw_on_failure) throw runtime_error("SingleTopGen: W has not ==2 daughters");
        return;
      }
            
      // now that we collected everything, fill the member variables. 
      // Unlike in TTbarGen.cxx, where we use different member variables according to top quark charge, we do not distinguish them here because there is just one top
      m_Top = top;
      m_WTop = *w;
      m_bTop = *b;
      m_Wdecay1 = *wd1;
      m_Wdecay2 = *wd2;
      ++n_top;

    }
  }

  if(n_top != 1){
    if(throw_on_failure)  throw runtime_error("SingleTopGen: did not find exactly one (anti)top in the event");
    return;
  }
    
  // calculate decay channel by counting the number of charged leptons
  // in the W daughters:
  int n_e = 0, n_m = 0, n_t = 0;
  for(const auto & wd : {m_Wdecay1, m_Wdecay2}){
    int id = abs(wd.pdgId());
    if(id == 11) ++n_e;
    else if(id == 13) ++n_m;
    else if(id == 15) ++n_t;
  }

  // leptonic:
  if(n_e == 1){
    m_type = e_ele;
  }
  else if(n_m == 1){
    m_type = e_muo;
  }
  else if(n_t == 1){
    m_type = e_tau;
  }
  // hadronic:
  else{
    m_type = e_had;
  }
}   
 

GenParticle SingleTopGen::Top() const{
  return m_Top;
}

GenParticle SingleTopGen::WTop() const{
  return m_WTop;
}

GenParticle SingleTopGen::bTop() const{
  return m_bTop;
}

GenParticle SingleTopGen::Wdecay1() const{
  return m_Wdecay1;
} 

GenParticle SingleTopGen::Wdecay2() const{
  return m_Wdecay2;
} 

SingleTopGen::E_DecayChannel SingleTopGen::DecayChannel()  const{  
  return m_type;
}

bool SingleTopGen::IsTopHadronicDecay() const{
  return abs(m_Wdecay1.pdgId()) <= 5;
}

bool SingleTopGen::IsTopLeptonicDecay() const{
  return m_type == e_ele || m_type == e_muo || m_type == e_tau;
}

bool SingleTopGen::IsTopToElectronDecay() const{
  return m_type == e_ele;
}

bool SingleTopGen::IsTopToMuonDecay() const{
  return m_type == e_muo;
}

bool SingleTopGen::IsTopToTauonDecay() const{
  return m_type == e_tau;
}

SingleTopGenProducer::SingleTopGenProducer(uhh2::Context & ctx, const std::string & name, bool throw_on_failure_): throw_on_failure(throw_on_failure_){
  h_singletopgen = ctx.get_handle<SingleTopGen>(name);
}

bool SingleTopGenProducer::process(Event & event){
  event.set(h_singletopgen, SingleTopGen(*event.genparticles, throw_on_failure));
  return true;
}
