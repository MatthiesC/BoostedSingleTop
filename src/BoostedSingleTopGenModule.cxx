#include <iostream>
#include <memory>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
//#include "UHH2/common/include/CommonModules.h"

#include "UHH2/BoostedSingleTop/include/SingleTopGen.h"
#include "UHH2/BoostedSingleTop/include/SingleTopGenHists.h"
//#include "UHH2/BoostedSingleTop/include/SingleTopGen_tWch.h"
//#include "UHH2/BoostedSingleTop/include/SingleTopGenHists_tWch.h"

using namespace std;
using namespace uhh2;

namespace uhh2examples {

  class BoostedSingleTopGenModule: public AnalysisModule {
  public:
    
    explicit BoostedSingleTopGenModule(Context & ctx);
    virtual bool process(Event & event) override;

  private:
    
    std::unique_ptr<AnalysisModule> SingleTopGenProd;
    
    std::unique_ptr<Hists> h_singletopgen;
  };


  BoostedSingleTopGenModule::BoostedSingleTopGenModule(Context & ctx){
    
    cout << "Hello World from BoostedSingleTopGenModule!" << endl;
    
    string testvalue = ctx.get("TestKey", "<not set>");
    cout << "TestKey in the configuration was: " << testvalue << endl;
    
    for(auto & kv : ctx.get_all()){
      cout << " " << kv.first << " = " << kv.second << endl;
    }
    
    //common.reset(new CommonModules());
    // TODO: configure common here, e.g. by 
    // calling common->set_*_id or common->disable_*
    //common->disable_jec();
    //common->init(ctx);

    SingleTopGenProd.reset(new SingleTopGenProducer(ctx, "singletopgen"));
    h_singletopgen.reset(new SingleTopGenHists(ctx, "singletopgen"));
  }


  bool BoostedSingleTopGenModule::process(Event & event) {

    //common->process(event);
    SingleTopGenProd->process(event);
    h_singletopgen->fill(event);
    
    return false;
  }

  UHH2_REGISTER_ANALYSIS_MODULE(BoostedSingleTopGenModule)

}
