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
  
  ////////////////////
  //Neutrino Flux gunk
  ////////////////////
  Flux flux;
  std::pair<double,double> p = flux.calculate_flux();
  double N_targets = p.first;
  double flux_value = p.second;

  ////////////////////
  //For the MC Plots
  ////////////////////
  
  //double sigma_empirical = 30.2249 * 1E-38; //steven's value
  //double sigma_empirical = 2.73478 * 1E-37; //New value as of 12/7/2021
  double sigma_empirical = 37.2171E-38; //value from Afro
  
  double sigma_nieves = 27.0266 * 1E-38; //steven's value of nieves
  //double sigma_nieves = 27.1682 * 1E-38;  //my value of nieves

  double sigma_susa = 38.1417 * 1E-38; //my value

  //N is the total number of events generated in each MEC type
  ///////////////////////////////////////////////////////////
  int N_empirical = 4000000;
  //int N_empirical = 1000000; //new value as of 12/7/2021
  int N_nieves = 2400000;
  int N_susa = 3800000;
  int N_GCF = 500000;

  //XSec Plot Colors
  ////////////////////////
  TColor tcolor;
  Int_t color1 = tcolor.GetColor("#FFCE54"); //yellow
  Int_t color2 = tcolor.GetColor("#A0D568"); //green
  Int_t color3 = tcolor.GetColor("#4FC1E8"); //blue
  Int_t color4 = tcolor.GetColor("#ED5564"); //red
  Int_t color5 = tcolor.GetColor("#AC92EB"); //purple  TColor* color0 = kBlack;
  Int_t xsec_colors[] = {color1,color2,color3,color4,color5};
  
  /////////////////////////////////
  //Stuff Needed for the Overlay Plots
  ////////////////////////////////
  static const int num_channels = 2;
  const char* channel[num_channels] = {"_total", "_cc2p0pi"};
  
  ///////////////////////////
  //Stuff for the particle plots
  //////////////////////////////
  static const int num_var = 3;
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta) (Rad.)","#phi (Rad.)"};

  static const int num_particles = 3;
  const char* particles[num_particles] = {"_muon","_leading","_recoil"};
  const char* particles_titles[num_particles] = {"Muon","Leading Proton","Recoiling Proton"};

  //Muon
  ///////////////
  TH1D* h_bnb_muon[num_var];
  TH1D* h_ext_muon[num_var];
  TH1D* h_dirt_muon[num_var];
  TH1D* h_overlay_muon[num_var][num_channels];

  TH1D* h_muon_eff[num_var]; //efficiency 
  TH2D* h_muon_matrices[num_var]; //smearing matrix
  TH1D* h_muon_num[num_var]; //numberator (for double checking)
  TH1D* h_muon_denom[num_var]; //denominator (for double checking)

  TH1D* h_nuwro_muon[num_var]; //nuwro events
  TH1D* h_empirical_muon[num_var]; //empirical events
  TH1D* h_nieves_muon[num_var]; //nieves events
  TH1D* h_susa_muon[num_var]; //susa v2 events
  TH1D* h_GCF_muon[num_var]; //GCF events

  TCanvas* canv_muon_eff[num_var];
  TCanvas* canv_muon_2d[num_var];

  double muon_max[num_var] = {15000,15000,4000};
  double muon_max_theory[num_var] = {0.5,0.7,0.15};
  //double muon_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  double muon_xsec_max[num_var] = {17,12,2};
  bool flip_legend_muon[num_var] = {false,true,false};
  int muon_num_iterations[num_var] = {3,2,2};//{3,2,2};
  TH1D* h_muon_systematic[num_var];

  //Leading Proton
  /////////////////
  TH1D* h_bnb_leading[num_var];
  TH1D* h_ext_leading[num_var];
  TH1D* h_dirt_leading[num_var];
  TH1D* h_overlay_leading[num_var][num_channels];

  TH1D* h_leading_eff[num_var]; //efficiency 
  TH2D* h_leading_matrices[num_var]; //smearing matrix
  TH1D* h_leading_num[num_var]; //numberator (for double checking)
  TH1D* h_leading_denom[num_var]; //denominator (for double checking)

  TH1D* h_nuwro_leading[num_var]; //nuwro events
  TH1D* h_empirical_leading[num_var]; //empirical events
  TH1D* h_nieves_leading[num_var]; //nieves events
  TH1D* h_susa_leading[num_var]; //susa v2 events
  TH1D* h_GCF_leading[num_var]; //GCF events

  TCanvas* canv_leading_eff[num_var];
  TCanvas* canv_leading_2d[num_var];

  double leading_max[num_var] = {15000,15000,4000};
  double leading_max_theory[num_var] = {0.5,0.7,0.15};
  //double leading_xsec_max[num_var] = {25E-38,12E-38,3E-38};
  double leading_xsec_max[num_var] = {25,12,3};
  bool flip_legend_leading[num_var] = {false,true,false};
  int leading_num_iterations[num_var] = {2,3,2};//{2,2,2}
  TH1D* h_leading_systematic[num_var];
  
  //Recoil Proton
  ////////////////////
  TH1D* h_bnb_recoil[num_var];
  TH1D* h_ext_recoil[num_var];
  TH1D* h_dirt_recoil[num_var];
  TH1D* h_overlay_recoil[num_var][num_channels];

  TH1D* h_recoil_eff[num_var]; //efficiency 
  TH2D* h_recoil_matrices[num_var]; //smearing matrix
  TH1D* h_recoil_num[num_var]; //numberator (for double checking)
  TH1D* h_recoil_denom[num_var]; //denominator (for double checking)

  TH1D* h_nuwro_recoil[num_var]; //nuwro events
  TH1D* h_empirical_recoil[num_var]; //empirical events
  TH1D* h_nieves_recoil[num_var]; //nieves events
  TH1D* h_susa_recoil[num_var]; //susa v2 events
  TH1D* h_GCF_recoil[num_var]; //GCF events

  TCanvas* canv_recoil_eff[num_var];
  TCanvas* canv_recoil_2d[num_var];

  double recoil_max[num_var] = {15000,15000,4000};
  double recoil_max_theory[num_var] = {0.8,0.5,0.15};
  //double recoil_xsec_max[num_var] = {35E-38,8E-38,2E-38};
  double recoil_xsec_max[num_var] = {35,8,2};
  bool flip_legend_recoil[num_var] = {false,true,false};
  int recoil_num_iterations[num_var] = {2,2,2};//{2,2,2}
  TH1D* h_recoil_systematic[num_var];
  
  ////////////////////////////
  //Stuff for the other plots
  //////////////////////////

  
  static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
  					  "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
  						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
  						 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};

  //const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
  //					  "_delta_PT","_delta_alphaT","_delta_phiT","_nu_E","_pn"};
  //const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
  //						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
  //						 "Estimated Neutrino Energy (GeV/c)","Estimated p_{n} (GeV/c)"};
  
  /*static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                               "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  */

  TH1D* h_bnb_other[num_other_var];
  TH1D* h_ext_other[num_other_var];
  TH1D* h_dirt_other[num_other_var];
  TH1D* h_overlay_other[num_other_var][num_channels];

  TH1D* h_other_eff[num_other_var]; //efficiency 
  TH2D* h_other_matrices[num_other_var]; //smearing matrix
  TH1D* h_other_num[num_other_var]; //numberator (for double checking)
  TH1D* h_other_denom[num_other_var]; //denominator (for double checking)

  TH1D* h_nuwro_other[num_other_var]; //nuwro events
  TH1D* h_empirical_other[num_other_var]; //empirical events
  TH1D* h_nieves_other[num_other_var]; //nieves events
  TH1D* h_susa_other[num_other_var]; //susa v2 events
  TH1D* h_GCF_other[num_other_var]; //GCF events

  TH1D* h_nuwro_other_true[num_other_var]; //nuwro events
  TH1D* h_empirical_other_true[num_other_var]; //empirical events
  TH1D* h_nieves_other_true[num_other_var]; //nieves events
  TH1D* h_susa_other_true[num_other_var]; //susa v2 events
  TH1D* h_GCF_other_true[num_other_var]; //GCF events  
  
  TCanvas* canv_other_eff[num_other_var];
  TCanvas* canv_other_2d[num_other_var];

  double other_max[num_other_var] = {15000,15000,6000,15000,15000,15000,15000};
  double other_max_theory[num_other_var] = {0.3,0.4,0.3,0.4,0.7,0.7,1.0};
  //double other_xsec_max[num_other_var] = {8E-38,8E-38,25E-38,0.1E-38,0.1E-38,35E-38,15E-38};
  double other_xsec_max[num_other_var] = {8,8,20,0.1,0.08,35,15};
  bool flip_legend_other[num_other_var] = {false,true,false,false,false,false,false};
  int other_num_iterations[num_other_var] = {2,2,3,3,2,3,4}; //{2,2,3,3,2,3,5}
  TH1D* h_other_systematic[num_other_var];
  
}

#endif
