#define systematics_cxx
#include "systematics.h"

void systematics::main(){


  //Grabbing the number of iterations
  //////////////////////////////////

  
  //Run the Detector Variations
  /////////////////////////////
  std::cout<<"[RUNNING THE DETECTOR SYSTEMATICS]"<<std::endl;
  //detvar.main();
  std::cout<<"[FINISHED RUNNING THE DETECTOR SYSTEMATICS]"<<std::endl;
  
  //Run the Multisims and Unisims
  //////////////////////////////
  std::cout<<"[RUNNING THE ALL THE MULTISIMS]"<<std::endl;
  //multisims.main(); 
  std::cout<<"[FINISHED RUNNING THE ALL THE MULTISIMS]"<<std::endl;
  
  //Run the Dirt Systematics
  //////////////////////////
  std::cout<<"[RUNNING THE DIRT SYSTEMATICS]"<<std::endl;
  dirt.main();
  std::cout<<"[FINISHED RUNNING THE DIRT SYSTEMATICS]"<<std::endl;
  

} //end of code
