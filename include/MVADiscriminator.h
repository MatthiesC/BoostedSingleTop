#pragma once

#include "TMVA/Reader.h"

class MVADiscriminator
{
 public:
  MVADiscriminator(const std::string&);
  ~MVADiscriminator();
  float eval(float var1, float var2, float var3, float var4, float var5, float var6, float var7, float var8, float var9);
  
 private:
  TMVA::Reader *reader_;
  float floatvars_[9];
};
