#include "UHH2/BoostedSingleTop/include/myutils.h"
#include "UHH2/core/include/LorentzVector.h"


// credits to Torben
bool TopJetGroomer::process(uhh2::Event & event) { 

  assert(event.topjets);

  for(auto & topjet : *event.topjets){
    LorentzVector subjet_sum(0,0,0,0);

    for(const auto & subjet : topjet.subjets()){
      subjet_sum += subjet.v4();
          }
    topjet.set_v4(subjet_sum);
  }

  return true;
}
