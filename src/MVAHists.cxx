#include "UHH2/BoostedSingleTop/include/MVAHists.h"
#include "UHH2/core/include/Event.h"

//#include "UHH2/common/include/Utils.h"

#include "TH1F.h"

using namespace std;
using namespace uhh2;

MVAHists::MVAHists(Context & ctx, const string & dirname):
  Hists(ctx,dirname)
{
  // Book all histograms here ...
  H_mvaD = book<TH1F>("mvaD", "BDT discriminator", 40, -1, 1);
  H_mvaD_1b = book<TH1F>("mvaD_1btags", "BDT discriminator (1 b-tag)", 40, -1, 1);
  H_mvaD_2b = book<TH1F>("mvaD_2btags", "BDT discriminator (2 b-tags)", 40, -1, 1);

}

void MVAHists::fill(const Event & e)
{
  double w = e.weight;
}

void MVAHists::fill_(const Event & e, const double & mvaD, const int & n_btags)
{
  double w = e.weight;

  // Now fill the histograms ...
  H_mvaD->Fill(mvaD, w);
  if(n_btags == 1) H_mvaD_1b->Fill(mvaD, w);
  else if(n_btags == 2) H_mvaD_2b->Fill(mvaD, w);
  
}

MVAHists::~MVAHists(){}
