#ifndef CONSTANTS_H
#define CONSTANTS_H
namespace Constants{

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};
  
  //Muon
  TH1D* h_muon[1000][num_var];
  TH1D* h_muon_CV[num_var];
  TH1D* h_muon_BNB[num_var];
  TH1D* h_muon_EXT[num_var];
  TH1D* h_muon_Dirt[num_var];
  double ymax_muon[num_var] = {2000,2000,2000};
  double ymin_muon[num_var] = {0,0,0};
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
  TH2D* h_2D_muon[num_var];
  TCanvas* canv_2D_muon[num_var];
  TH1D* h_muon_error[num_var];
  TCanvas* canv_muon_error[num_var];
  TH1D* h_muon_fractional_error[num_var];
  TCanvas* canv_muon_fractional_error[num_var];
  
  //Lead
  TH1D* h_leading[1000][num_var];
  TH1D* h_leading_BNB[num_var];
  TH1D* h_leading_CV[num_var];
  TH1D* h_leading_EXT[num_var];
  TH1D* h_leading_Dirt[num_var];
  double ymax_leading[num_var] = {2000,2000,2000};
  double ymin_leading[num_var] = {0,0,0};
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  TH2D* h_2D_leading[num_var];
  TCanvas* canv_2D_leading[num_var];
  TH1D* h_leading_error[num_var];
  TCanvas* canv_leading_error[num_var];
   TH1D* h_leading_fractional_error[num_var];
  TCanvas* canv_leading_fractional_error[num_var];
  
  //Recoil
  TH1D* h_recoil[1000][num_var];
  TH1D* h_recoil_BNB[num_var];
  TH1D* h_recoil_CV[num_var];
  TH1D* h_recoil_EXT[num_var];
  TH1D* h_recoil_Dirt[num_var];
  double ymax_recoil[num_var] = {2000,2000,2000};
  double ymin_recoil[num_var] = {0,0,0};
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  TH2D* h_2D_recoil[num_var];
  TCanvas* canv_2D_recoil[num_var];
  TH1D* h_recoil_error[num_var];
  TCanvas* canv_recoil_error[num_var];
  TH1D* h_recoil_fractional_error[num_var];
  TCanvas* canv_recoil_fractional_error[num_var];
  
  //Other Variables
  /////////////////////////
  /*  static const int num_other_var = 9;
      const	char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading",
      "_opening_angle_mu_both","_delta_PT","_delta_alphaT",
      "_delta_phiT","_pn","_nu_E"};
      const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})","cos(#gamma_{#mu,p_{L}})",
      "cos(#gamma_{#mu,p_{L} + p_{R}})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)",
      "#delta #phi_{T} (Deg.)","Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"}; 
  */
  static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
					  "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
						 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};
  
  TH1D* h_other[1000][num_other_var];
  TH1D* h_other_CV[num_other_var];
  TH1D* h_other_BNB[num_other_var];
  TH1D* h_other_EXT[num_other_var];
  TH1D* h_other_Dirt[num_other_var];
  double ymax_other[num_other_var] = {2000,2000,2000,2000,2000,2000,2000};
  double ymin_other[num_other_var] = {0,0,0,0,0,0,0};
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  TH2D* h_2D_other[num_other_var];
  TCanvas* canv_2D_other[num_other_var];
  TH1D* h_other_error[num_other_var];
  TCanvas* canv_other_error[num_other_var];
   TH1D* h_other_fractional_error[num_other_var];
  TCanvas* canv_other_fractional_error[num_other_var];
}

#endif
