#include <cmath>
#include "BoostedDiscriminatorMVADesy.h"

BoostedDiscriminatorMVADesy::BoostedDiscriminatorMVADesy()
{
  reader_  = new TMVA::Reader("!Color:!Silent");

  reader_->AddVariable("jettau1Leading",&var_[0]);
  reader_->AddVariable("jettau2Leading",&var_[1]);
  reader_->AddVariable("jettau3Leading",&var_[2]);
  //reader_->AddVariable("tau31Leading",&var_[0]);
  //reader_->AddVariable("tau32Leading",&var_[1]);
  //reader_->AddVariable("ptLeading",&var_[5]);
  //reader_->AddVariable("SDmassLeading",&var_[6]);
  //reader_->AddVariable("MassLeading",&var_[2]);
  //reader_->AddVariable("MassSubJet0Leading",&var_[3]);
  //reader_->AddVariable("MassSubJet1Leading",&var_[9]);
  //reader_->AddVariable("jetSubJetpt0Leading",&var_[10]);
  //reader_->AddVariable("jetSubJetpt1Leading",&var_[11]);
  //reader_->AddVariable("DeltaRSubJetsLeading",&var_[12]);

  reader_->AddVariable("jettau1SubLeading",&var_[3]);
  reader_->AddVariable("jettau2SubLeading",&var_[4]);
  reader_->AddVariable("jettau3SubLeading",&var_[5]);
  //reader_->AddVariable("tau31SubLeading",&var_[4]);
  //reader_->AddVariable("tau32SubLeading",&var_[5]);
  //reader_->AddVariable("ptSubLeading",&var_[18]);
  //reader_->AddVariable("SDmassSubLeading",&var_[19]);
  //reader_->AddVariable("MassSubLeading",&var_[6]);
  //reader_->AddVariable("MassSubJet0SubLeading",&var_[7]);
  //reader_->AddVariable("MassSubJet1SubLeading",&var_[22]);
  //reader_->AddVariable("jetSubJetpt0SubLeading",&var_[23]);
  //reader_->AddVariable("jetSubJetpt1SubLeading",&var_[24]);
  //reader_->AddVariable("DeltaRSubJetsSubLeading",&var_[25]);
  
  //reader_->BookMVA("LD","../test/weights/TMVAClassification_LD.weights.xml");  //4 variables
  reader_->BookMVA("LD","../test/weights/TMVAClassification_MLP.weights.xml");  //8 variables
  //reader_->BookMVA("LD","../test/weights/Daniela/TMVAClassification_BDT.weights.xml");  //Daniela

}
//-------------------------------------------------------------
//float BoostedDiscriminatorMVADesy::eval(float var0,float var1,float var2,float var3,float var4,float var5,float var6,float var7)
float BoostedDiscriminatorMVADesy::eval(float var0,float var1, float var2,float var3, float var4,float var5)
{

  var_[0] = var0;
  var_[1] = var1;
  var_[2] = var2;
  var_[3] = var3;
  var_[4] = var4;
  var_[5] = var5;
  //var_[6] = var6;
  //var_[7] = var7;
      
  return reader_->EvaluateMVA("LD");
}
//-------------------------------------------------------------
BoostedDiscriminatorMVADesy::~BoostedDiscriminatorMVADesy()
{
  delete reader_;
}
