#ifndef CONSTANTS_H
#define CONSTANTS_H
#include "RooUnfoldResponse.h"
#include "RooUnfoldBayes.h"
#include "TColor.h"
#include "TStyle.h"
#include "neutrino_flux.h"

namespace Constants{

  //Pot num and sample name
  ///////////////////////
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";
  
  // //XSec Plot Colors                                                                                                                                                       
  ////////////////////////                                                                                                                                                 
  TColor tcolor;
  Int_t color1 = tcolor.GetColor("#FFCE54"); //yellow                                                                                                                      
  Int_t color2 = tcolor.GetColor("#A0D568"); //green                                                                                                                       
  Int_t color3 = tcolor.GetColor("#4FC1E8"); //blue                                                                                                                        
  Int_t color4 = tcolor.GetColor("#ED5564"); //red                                                                                                                         
  Int_t color5 = tcolor.GetColor("#AC92EB"); //purple  TColor* color0 = kBlack;                                                                                            
  Int_t xsec_colors[] = {color1,color2,color3,color4,color5};

  ////////////////////
  //Neutrino Flux gunk
  ////////////////////
  Flux flux;
  std::pair<double,double> p = flux.calculate_flux();
  double N_targets = p.first;
  double flux_value = p.second;

}

#endif
