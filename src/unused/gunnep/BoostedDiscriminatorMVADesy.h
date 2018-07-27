#ifndef BoostedDiscriminatorMVADesy_h
#define BoostedDiscriminatorMVADesy_h

#include "TMVA/Reader.h"

class BoostedDiscriminatorMVADesy
{
  public:
    BoostedDiscriminatorMVADesy();
    ~BoostedDiscriminatorMVADesy();
    //float eval(float var0,float var1,float var2,float var3,float var4,float var5,float var6,float var7);
    float eval(float var0,float var1, float var2,float var3,float var4, float var5);
        
  private:
    TMVA::Reader *reader_;
    float var_[6];
};
#endif
