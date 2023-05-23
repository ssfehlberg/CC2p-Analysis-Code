#ifndef HISTOGRAMS_H
#define HISTOGRAMS_H

namespace Histograms{

  //CV Files we will grab CV plots from
  /////////////////////////////////
  TFile* file_overlay = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_eff = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_mc_eff.root"));
  TFile* file_bnb = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));

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
  TH1D* h_bnb_muon[num_var]; //CV BNB
  TH1D* h_ext_muon[num_var]; //CV EXT
  TH1D* h_dirt_muon[num_var]; //CV Dirt
  TH1D* h_overlay_muon[num_var][num_channels]; //CV Overlay for total and cc2p
  
  //cross-section drawing stuffs
  TH1D* h_muon_num[num_var]; //numberator (for double checking)
  TH1D* h_muon_denom[num_var]; //denominator (for double checking)
  TH1D* h_muon_eff[num_var]; //efficiency 
  TCanvas* canv_muon_eff[num_var];

  TH2D* h_muon_matrices[num_var]; //smearing matrix
  TCanvas* canv_muon_2d[num_var];

  TH1D* h_muon_xsec[1000][num_var]; //cross-section in each universe
  double muon_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_muon[num_var] = {false,true,false};
  int muon_iter[num_var]; //number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_muon[num_var]; 
  TLegend* legend_muon[num_var];
  double ymax_muon[num_var] = {2000,2000,2000};
  double ymin_muon[num_var] = {0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_muon[num_var];
  TCanvas* canv_2D_muon[num_var];

  //Error drawing stuff
  TH1D* h_muon_error[num_var];
  TCanvas* canv_muon_error[num_var];

  //fractional error stuff
  TH1D* h_muon_fractional_error[num_var];
  TCanvas* canv_muon_fractional_error[num_var];

  //Leading Proton
  /////////////////
  TH1D* h_bnb_leading[num_var]; //CV BNB
  TH1D* h_ext_leading[num_var]; //CV EXT
  TH1D* h_dirt_leading[num_var]; //CV Dirt
  TH1D* h_overlay_leading[num_var][num_channels]; //CV Overlay for total and cc2p
  
  //cross-section drawing stuff
  TH1D* h_leading_num[num_var]; //numberator (for double checking)
  TH1D* h_leading_denom[num_var]; //denominator (for double checking)
  TH1D* h_leading_eff[num_var]; //efficiency 
  TCanvas* canv_leading_eff[num_var];

  TH2D* h_leading_matrices[num_var]; //smearing matrix
  TCanvas* canv_leading_2d[num_var];

  TH1D* h_leading_xsec[1000][num_var]; //cross-section in each universe
  double leading_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_leading[num_var] = {false,true,false};
  int leading_iter[num_var];//number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_leading[num_var]; 
  TLegend* legend_leading[num_var];
  double ymax_leading[num_var] = {2000,2000,2000};
  double ymin_leading[num_var] = {0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_leading[num_var];
  TCanvas* canv_2D_leading[num_var];

  //Error drawing stuff
  TH1D* h_leading_error[num_var];
  TCanvas* canv_leading_error[num_var];

  //fractional error stuff
  TH1D* h_leading_fractional_error[num_var];
  TCanvas* canv_leading_fractional_error[num_var];

  //Recoil Proton
  /////////////////
  TH1D* h_bnb_recoil[num_var]; //CV BNB
  TH1D* h_ext_recoil[num_var]; //CV EXT
  TH1D* h_dirt_recoil[num_var]; //CV Dirt
  TH1D* h_overlay_recoil[num_var][num_channels]; //CV Overlay for total and cc2p
  
  //cross-section drawing stuff
  TH1D* h_recoil_num[num_var]; //numberator (for double checking)
  TH1D* h_recoil_denom[num_var]; //denominator (for double checking)
  TH1D* h_recoil_eff[num_var]; //efficiency 
  TCanvas* canv_recoil_eff[num_var];

  TH2D* h_recoil_matrices[num_var]; //smearing matrix
  TCanvas* canv_recoil_2d[num_var];

  TH1D* h_recoil_xsec[1000][num_var]; //cross-section in each universe
  double recoil_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_recoil[num_var] = {false,true,false};
  int recoil_iter[num_var]; //number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_recoil[num_var]; 
  TLegend* legend_recoil[num_var];
  double ymax_recoil[num_var] = {2000,2000,2000};
  double ymin_recoil[num_var] = {0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_recoil[num_var];
  TCanvas* canv_2D_recoil[num_var];

  //Error drawing stuff
  TH1D* h_recoil_error[num_var];
  TCanvas* canv_recoil_error[num_var];

  //fractional error stuff
  TH1D* h_recoil_fractional_error[num_var];
  TCanvas* canv_recoil_fractional_error[num_var];


  ///////////////////////
  //Stuff for other plots
  //////////////////////
  
  static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
					  "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
						 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};
  /*
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  */

  TH1D* h_bnb_other[num_other_var]; //CV BNB
  TH1D* h_ext_other[num_other_var]; //CV EXT
  TH1D* h_dirt_other[num_other_var]; //CV Dirt
  TH1D* h_overlay_other[num_other_var][num_channels]; //CV Overlay for total and cc2p
  
  //cross-section drawing stuff
  TH1D* h_other_num[num_other_var]; //numberator (for double checking)
  TH1D* h_other_denom[num_other_var]; //denominator (for double checking)
  TH1D* h_other_eff[num_other_var]; //efficiency 
  TCanvas* canv_other_eff[num_other_var];

  TH2D* h_other_matrices[num_other_var]; //smearing matrix
  TCanvas* canv_other_2d[num_other_var];

  TH1D* h_other_xsec[1000][num_other_var]; //cross-section in each universe
  double other_xsec_max[num_other_var] = {8E-38,8E-38,25E-38,0.18E-38,0.1E-38,2E-38,2E-38};    
  bool flip_legend_other[num_other_var] = {false,true,false,false,false,false,false};
  int other_iter[num_other_var] = {2,2,1,1,2,1,1}; //number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_other[num_other_var]; 
  TLegend* legend_other[num_other_var];
  double ymax_other[num_other_var] = {2000,2000,2000,2000,2000,2000,2000};
  double ymin_other[num_other_var] = {0,0,0,0,0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_other[num_other_var];
  TCanvas* canv_2D_other[num_other_var];

  //Error drawing stuff
  TH1D* h_other_error[num_other_var];
  TCanvas* canv_other_error[num_other_var];

  //fractional error stuff
  TH1D* h_other_fractional_error[num_other_var];
  TCanvas* canv_other_fractional_error[num_other_var];
  
}

#endif
  
