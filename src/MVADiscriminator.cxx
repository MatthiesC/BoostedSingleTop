#include <cmath>
#include "UHH2/BoostedSingleTop/include/MVADiscriminator.h"
#include "TMVA/Reader.h"

using namespace std;

MVADiscriminator::MVADiscriminator(const string & path_weightfile)
{
  reader_ = new TMVA::Reader("Silent"); // "Silent" to suppress all screen output

  reader_->AddVariable("tmva_n_btags",&floatvars_[0]);
  reader_->AddVariable("tmva_deltaphi_bottomlepton",&floatvars_[1]);
  reader_->AddVariable("tmva_deltaphi_bottomtoptag",&floatvars_[2]);
  reader_->AddVariable("tmva_pseudotop_mass",&floatvars_[3]);
  reader_->AddVariable("tmva_ptbalance",&floatvars_[4]);
  reader_->AddVariable("tmva_wass_pt",&floatvars_[5]);
  reader_->AddVariable("tmva_lepton_pt",&floatvars_[6]);
  reader_->AddVariable("tmva_lepton_eta",&floatvars_[7]);
  reader_->AddVariable("tmva_deltaphi_leptonnextjet",&floatvars_[8]);
  //reader_->AddVariable("tmva_met_pt",&floatvars_[9]);
  /*reader_->AddVariable("tmva_top_pt",&floatvars_[9]);
  reader_->AddVariable("tmva_top_eta",&floatvars_[10]);
  reader_->AddVariable("tmva_ht_jets",&floatvars_[11]);*/
  
  reader_->BookMVA("BDT",path_weightfile);

}

float MVADiscriminator::eval(float var1, float var2, float var3, float var4, float var5, float var6, float var7, float var8, float var9)
{

  floatvars_[0] = var1;
  floatvars_[1] = var2;
  floatvars_[2] = var3;
  floatvars_[3] = var4;
  floatvars_[4] = var5;
  floatvars_[5] = var6;
  floatvars_[6] = var7;
  floatvars_[7] = var8;
  floatvars_[8] = var9;
  //floatvars_[9] = var10;
  /*floatvars_[10] = var11;
  floatvars_[11] = var12;*/
      
  return reader_->EvaluateMVA("BDT");
}

MVADiscriminator::~MVADiscriminator()
{
  delete reader_;
}












// for 1btag region:


MVADiscriminator_1b::MVADiscriminator_1b(const string & path_weightfile)
{
  reader_ = new TMVA::Reader("Silent"); // "Silent" to suppress all screen output

  //reader_->AddVariable("tmva_n_btags",&floatvars_[0]);
  reader_->AddVariable("tmva_deltaphi_bottomlepton",&floatvars_[0]);
  reader_->AddVariable("tmva_deltaphi_bottomtoptag",&floatvars_[1]);
  reader_->AddVariable("tmva_pseudotop_mass",&floatvars_[2]);
  reader_->AddVariable("tmva_ptbalance",&floatvars_[3]);
  reader_->AddVariable("tmva_wass_pt",&floatvars_[4]);
  reader_->AddVariable("tmva_lepton_pt",&floatvars_[5]);
  reader_->AddVariable("tmva_lepton_eta",&floatvars_[6]);
  reader_->AddVariable("tmva_deltaphi_leptonnextjet",&floatvars_[7]);
  //reader_->AddVariable("tmva_met_pt",&floatvars_[9]);
  /*reader_->AddVariable("tmva_top_pt",&floatvars_[9]);
  reader_->AddVariable("tmva_top_eta",&floatvars_[10]);
  reader_->AddVariable("tmva_ht_jets",&floatvars_[11]);*/
  
  reader_->BookMVA("BDT",path_weightfile);

}

float MVADiscriminator_1b::eval(float var1, float var2, float var3, float var4, float var5, float var6, float var7, float var8)
{

  floatvars_[0] = var1;
  floatvars_[1] = var2;
  floatvars_[2] = var3;
  floatvars_[3] = var4;
  floatvars_[4] = var5;
  floatvars_[5] = var6;
  floatvars_[6] = var7;
  floatvars_[7] = var8;
  //floatvars_[8] = var9;
  //floatvars_[9] = var10;
  /*floatvars_[10] = var11;
  floatvars_[11] = var12;*/
      
  return reader_->EvaluateMVA("BDT");
}

MVADiscriminator_1b::~MVADiscriminator_1b()
{
  delete reader_;
}
