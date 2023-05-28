#include "detVar.h"
#include "GENIE_total.h"


#include "shared0.h"
using namespace shared;

void all(){


  std::cout<<"Running Detector Variations"<<std::endl;
  detVar detvar;
  detvar.main();

  std::cout<<"Running GENIE Total"<<std::endl;
  GENIE_total genie;
  genie.main();

}
