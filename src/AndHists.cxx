#include "UHH2/BoostedSingleTop/include/AndHists.h"
#include "UHH2/BoostedSingleTop/include/BoostedSingleTopHists.h"
#include "UHH2/BoostedSingleTop/include/HOTVRHists.h"
#include "UHH2/common/include/LuminosityHists.h"
#include "UHH2/common/include/EventHists.h"
#include "UHH2/common/include/JetHists.h"

#include "UHH2/BoostedSingleTop/include/MyElectronHists.h"
#include "UHH2/BoostedSingleTop/include/MyMuonHists.h"


using namespace std;
using namespace uhh2;
AndHists::AndHists(Context &ctx, const string & dirname, const boost::optional<TopJetId> &id_topjet):
  Hists(ctx, dirname+"_Counter")
{
  // Counting Hist
  nevt = book<TH1F>("NEvt", "", 1,0,1);

  // Add common hists to vector
  hists_vector.push_back(new BoostedSingleTopHists(ctx, dirname + "_BoostedSingleTop"));
  if (id_topjet) hists_vector.push_back(new HOTVRHists(ctx, dirname + "_HOTVR", *id_topjet));
  else hists_vector.push_back(new HOTVRHists(ctx, dirname + "_HOTVR"));
  hists_vector.push_back(new LuminosityHists(ctx, dirname + "_Lumi"));
  hists_vector.push_back(new EventHists(ctx, dirname + "_Event"));
  hists_vector.push_back(new MyMuonHists(ctx, dirname + "_Muon"));
  hists_vector.push_back(new MyElectronHists(ctx, dirname + "_Electron"));
  hists_vector.push_back(new JetHists(ctx, dirname + "_Jet"));
}

void AndHists::fill(const Event & event) {  
  nevt->Fill(0., event.weight);
  for (Hists *hist : hists_vector)
    {
      hist->fill(event);
    }
}

void AndHists::add_hist(Hists *hist) { 
  hists_vector.push_back(hist);
}

AndHists::~AndHists() {
  for (unsigned int i = 0; i < hists_vector.size(); ++i)
    {
      delete hists_vector[i];
    }
  hists_vector.clear();
}
