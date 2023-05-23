#ifndef systematics_h
#define systematics_h

//C++ includes
#include <iostream>
#include <string>
#include <fstream>

//#include "detvar.h"
//#include "multisims.h"
#include "dirt.h"

class systematics{

 public:
  virtual void main();

  //Classes
  //////////////
  //DetVar detvar;
  //Process_Multisims multisims; 
  Dirt dirt;
  
};

#endif
#ifdef systematics_cxx

/*
void systematics::total_uncertainty(){

  double flux_uncertainty = std::pow(0.02,2); //Uncertainty on the flux number in the denominator of the cross-section 
  double N_targets_uncertainty = std::pow(0.01,2); //Uncertainty on the number of targets in the denominator of the cross-section

  double n_bins = ;

  for(int i=1; i < n_bins+1; i++){

    double detector_systematic = std::pow(,2);
    double dirt_systematic =  std::pow(,2);
    double flux_systematic = std::pow(,2);
    double reint_systematic = std::pow(,2);
    double GENIE_multisims = std::pow(,2);
    double AxFFCCQEshape = std::pow(,2);
    double DecayAngMEC = std::pow(,2);
    double NormCCCOH = std::pow(,2);
    double NormNCCOH = std::pow(,2);
    double ThetaDelta2NRad = std::pow(,2);
    double Theta_Delta2Npi =std::pow(,2);
    double VecFFCCQEshape = std::pow(,2);
    double XSecShape_CCMEC = std::pow(,2);
    double xsr_scc_Fa3_SCC = std::pow(,2);
    double xsr_scc_Fv3_SCC = std::pow(,2);
    double RPA_CCQE = std::pow(,2);
    double RootinoFix = std::pow(,2);
    

    double value = std::sqrt(flux_uncertainty + N_targets_uncertainty +
			     detector_systematic + );
    


  }
}
*/
#endif
