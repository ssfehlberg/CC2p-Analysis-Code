#ifndef other_detvar_h
#define other_detvar_h
#include "plotting.h"
#include "tools/xsec.h"
#include "tools/Covariance_Matrix.h"

class other_detvar{

public:
  virtual void main(); //main function
  
 private:
  virtual void  Grab_Histograms(const char* input_file, int num_universes);
  virtual TH2D* make_covariance_matrix(TH1D* h_CV, std::vector<TH1D*> h_universe, int num_universes,const char* name);
  virtual void  plot_covariance(const char* variable, const char* title, TCanvas* canv, TH2D* hist, TH1D* hist_CV, const char* directory);
  virtual TH1D* make_error_histogram(const char* name, TH1D* hist_CV,TH2D* h_covariance);
  virtual TH1D* make_fractional_error_histogram(const char* name, TH1D* hist_CV,TH1D* hist_error);


  //Plotting class definition
  /////////////////
  plotting_tools plot;

  //XSec Class Definition
  ///////////////////////
  xsec Xsec;

  //number of universes
  //all of these samples are unisims i.e. # universes = 1
  static const int universes = 1; 

  //CV stuff
  /////////////
  static const int num_CV = 2;
  const char* CV_files[num_CV] = {"detVar_CV","detVar_Overlay"};
  const char* CV_titles[num_CV] = {"Central Value", "CV Overlay"};
  int which_CV; //which CV sample are we using
  TFile* file_CV; //TFIle for Central Value stuff
  TFile* file_CV_eff;
  
  //File stuff
  ////////////
  static const int num_files = 9;
  const char* directory_name;
  const char* directory_name_list[num_files] = {"detVar_LY_Attenuation","detVar_LY_Down","detVar_LY_Rayleigh",
						"detVar_ThetaXZ", "detVar_ThetaYZ",
						"detVar_X","detVar_YZ",
						"detVar_Recombination","detVar_SCE"};

  const char* sample_titles[num_files] = {"LY Attenuation","LY Down","LY Rayleigh",
					  "ThetaXZ","ThetaYZ",
					  "X","YZ",
					  "Recombination","SCE"};



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
  
   //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  ///////////////
  TH1D* h_bnb_muon[num_var]; //CV BNB
  TH1D* h_ext_muon[num_var]; //CV EXT
  TH1D* h_dirt_muon[num_var]; //CV Dirt
  
  TH1D* h_overlay_muon[num_var][num_channels][num_CV]; //CV Overlay for total and cc2p
  TH1D* h_muon_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_muon_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_muon_num_CV[num_var][num_CV]; //CV muon numerator of efficiency
  TH1D* h_muon_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_muon_denom_CV[num_var][num_CV]; //CV muon denominator of the efficiency
  TH1D* h_muon_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_muon_matrices_CV[num_var][num_CV]; //the CV muon matrices
  TH2D* h_muon_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_muon_2d[num_var];
  
  TH1D* h_muon_xsec[1000][num_var]; //cross-section in each universe
  double muon_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_muon[num_var] = {false,true,false};
  int muon_iter[num_var] = {2,2,2}; //number of iterations used in the unfolding. Really dumbly defined

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
  
  /* //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon0[num_var];
  TH1D* h_muon_total_up[num_var];
  TH1D* h_muon_total_down[num_var];
  TH1D* h_muon_CV[num_CV][num_var];
  TH2D* h_2D_muon[num_var];
  TH1D* h_muon_error[num_var];
  TH1D* h_muon_fractional_error[num_var];
  TCanvas* canv_muon[num_var];
  TCanvas* canv_muon_error[num_var];
  TCanvas* canv_muon_fractional_error[num_var];
  std::vector<TH2D*> muon_2D;*/

  //Lead
  /////////////
  TH1D* h_bnb_leading[num_var]; //CV BNB
  TH1D* h_ext_leading[num_var]; //CV EXT
  TH1D* h_dirt_leading[num_var]; //CV Dirt
  
  TH1D* h_overlay_leading[num_var][num_channels][num_CV]; //CV Overlay for total and cc2p
  TH1D* h_leading_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_leading_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_leading_num_CV[num_var][num_CV]; //CV leading numerator of efficiency
  TH1D* h_leading_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_leading_denom_CV[num_var][num_CV]; //CV leading denominator of the efficiency
  TH1D* h_leading_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_leading_matrices_CV[num_var][num_CV]; //the CV leading matrices
  TH2D* h_leading_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_leading_2d[num_var];
  
  TH1D* h_leading_xsec[1000][num_var]; //cross-section in each universe
  double leading_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_leading[num_var] = {false,true,false};
  int leading_iter[num_var] = {2,2,2}; //number of iterations used in the unfolding. Really dumbly defined

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
  
  
  //Recoil
  ///////////////////////////////
  TH1D* h_bnb_recoil[num_var]; //CV BNB
  TH1D* h_ext_recoil[num_var]; //CV EXT
  TH1D* h_dirt_recoil[num_var]; //CV Dirt
  
  TH1D* h_overlay_recoil[num_var][num_channels][num_CV]; //CV Overlay for total and cc2p
  TH1D* h_recoil_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_recoil_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_recoil_num_CV[num_var][num_CV]; //CV recoil numerator of efficiency
  TH1D* h_recoil_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_recoil_denom_CV[num_var][num_CV]; //CV recoil denominator of the efficiency
  TH1D* h_recoil_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_recoil_matrices_CV[num_var][num_CV]; //the CV recoil matrices
  TH2D* h_recoil_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_recoil_2d[num_var];
  
  TH1D* h_recoil_xsec[1000][num_var]; //cross-section in each universe
  double recoil_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_recoil[num_var] = {false,true,false};
  int recoil_iter[num_var] = {2,2,2}; //number of iterations used in the unfolding. Really dumbly defined

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
  
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
					  "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)",
                                                 "#delta #phi_{T} (Deg.)"};//,"Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};

  TH1D* h_bnb_other[num_other_var]; //CV BNB
  TH1D* h_ext_other[num_other_var]; //CV EXT
  TH1D* h_dirt_other[num_other_var]; //CV Dirt
  
  TH1D* h_overlay_other[num_other_var][num_channels][num_CV]; //CV Overlay for total and cc2p
  TH1D* h_other_total[1000][num_other_var]; //total MC contribution in a universe
  TH1D* h_other_cc2p[1000][num_other_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_other_num_CV[num_other_var][num_CV]; //CV other numerator of efficiency
  TH1D* h_other_num[1000][num_other_var]; //variations of numberator of efficiency
  TH1D* h_other_denom_CV[num_other_var][num_CV]; //CV other denominator of the efficiency
  TH1D* h_other_denom[1000][num_other_var]; //variations of denominator of efficiency
  
  TH2D* h_other_matrices_CV[num_other_var][num_CV]; //the CV other matrices
  TH2D* h_other_matrices[1000][num_other_var]; //the variation matrices
  TCanvas* canv_other_2d[num_other_var];
  
  TH1D* h_other_xsec[1000][num_other_var]; //cross-section in each universe
  double other_xsec_max[num_other_var] = {8E-38,8E-38,25E-38,0.18E-38,0.1E-38};
  bool flip_legend_other[num_other_var] = {false,true,false,false,false};
  int other_iter[num_other_var] = {2,4,2,2,2}; //number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_other[num_other_var]; 
  TLegend* legend_other[num_other_var];
  double ymax_other[num_other_var] = {2000,2000,2000,2000,2000};
  double ymin_other[num_other_var] = {0,0,0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_other[num_other_var];
  TCanvas* canv_2D_other[num_other_var];

  //Error drawing stuff
  TH1D* h_other_error[num_other_var];
  TCanvas* canv_other_error[num_other_var];

  //fractional error stuff
  TH1D* h_other_fractional_error[num_other_var];
  TCanvas* canv_other_fractional_error[num_other_var];
  


}; //end of class definition

#endif
#ifdef other_detvar_cxx

void other_detvar::Grab_Histograms(const char* input_file, int num_universes){
  
  std::cout<<"[GRAABBING ALL THE FILES]"<<std::endl;

  ////////////////////////////////////////////////////////////////////////
  //Grab the overlay smearing matrix and efficiency for all the variables
  /////////////////////////////////////////////////////////////////////

  for(int j=0; j < num_CV; j++){
    file_CV = new TFile(Form("root_files/detVar/histograms_pelee_xsec_%s.root",CV_files[j]));	
    file_CV_eff = new TFile(Form("root_files/detVar/histograms_%s_mc_eff.root",CV_files[j]));	
     
    for(int i=0; i < num_var; i++){
      h_bnb_muon[i] = (TH1D*)file_bnb->Get(Form("h_muon%s_bnb",var0[i])); //bnb
      h_ext_muon[i] = (TH1D*)file_ext->Get(Form("h_muon%s_ext",var0[i])); //ext
      h_dirt_muon[i] = (TH1D*)file_dirt->Get(Form("h_muon%s_dirt_wgt",var0[i])); //dirt
      
      h_bnb_leading[i] = (TH1D*)file_bnb->Get(Form("h_leading%s_bnb",var0[i])); //bnb
      h_ext_leading[i] = (TH1D*)file_ext->Get(Form("h_leading%s_ext",var0[i])); //ext
      h_dirt_leading[i] = (TH1D*)file_dirt->Get(Form("h_leading%s_dirt_wgt",var0[i])); //dirt
      
      h_bnb_recoil[i] = (TH1D*)file_bnb->Get(Form("h_recoil%s_bnb",var0[i])); //bnb
      h_ext_recoil[i] = (TH1D*)file_ext->Get(Form("h_recoil%s_ext",var0[i])); //ext
      h_dirt_recoil[i] = (TH1D*)file_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[i])); //dirt
    
      h_overlay_muon[i][0][j] = (TH1D*)file_CV->Get(Form("h_muon%s_total_0",var0[i])); //total CV overlay
      h_overlay_muon[i][1][j] = (TH1D*)file_CV->Get(Form("h_muon%s_cc2p0pi_0",var0[i])); //cc2p0pi CV overlay
      h_muon_matrices_CV[i][j] =  (TH2D*)file_CV_eff->Get(Form("h_particle_matrices_muon_all%s_0",var[i])); //CV matrix
      h_muon_num_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_num_muon_all%s_0",var[i])); //CV numerator 
      h_muon_denom_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_denom_muon_all%s_0",var[i])); //CV denominator

      h_overlay_leading[i][0][j] = (TH1D*)file_CV->Get(Form("h_leading%s_total_0",var0[i])); //total CV overlay
      h_overlay_leading[i][1][j] = (TH1D*)file_CV->Get(Form("h_leading%s_cc2p0pi_0",var0[i])); //cc2p0pi CV overlay
      h_leading_matrices_CV[i][j] =  (TH2D*)file_CV_eff->Get(Form("h_particle_matrices_lead_proton%s_0",var[i])); //CV matrix
      h_leading_num_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_num_lead_proton%s_0",var[i])); //CV numerator 
      h_leading_denom_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_denom_lead_proton%s_0",var[i])); //CV denominator

      h_overlay_recoil[i][0][j] = (TH1D*)file_CV->Get(Form("h_recoil%s_total_0",var0[i])); //total CV overlay
      h_overlay_recoil[i][1][j] = (TH1D*)file_CV->Get(Form("h_recoil%s_cc2p0pi_0",var0[i])); //cc2p0pi CV overlay
      h_recoil_matrices_CV[i][j] =  (TH2D*)file_CV_eff->Get(Form("h_particle_matrices_recoil_proton%s_0",var[i])); //CV matrix
      h_recoil_num_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_num_recoil_proton%s_0",var[i])); //CV numerator 
      h_recoil_denom_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_particle_denom_recoil_proton%s_0",var[i])); //CV denominator
      
    }
  

    for(int i=0; i < num_other_var; i++){
      h_bnb_other[i] = (TH1D*)file_bnb->Get(Form("h%s_bnb",other_var[i])); //bnb                                                                                       
      h_ext_other[i] = (TH1D*)file_ext->Get(Form("h%s_ext",other_var[i])); //ext                                                                                       
      h_dirt_other[i] = (TH1D*)file_dirt->Get(Form("h%s_dirt_wgt",other_var[i])); //dirt  
      
      h_overlay_other[i][0][j] = (TH1D*)file_CV->Get(Form("h%s_total_0",other_var[i])); //total CV overlay
      h_overlay_other[i][1][j] = (TH1D*)file_CV->Get(Form("h%s_cc2p0pi_0",other_var[i])); //cc2p0pi CV overlay
      h_other_matrices_CV[i][j] =  (TH2D*)file_CV_eff->Get(Form("h_other_matrices%s_0",other_var[i])); //CV matrix
      h_other_num_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_other_eff_num%s_0",other_var[i])); //CV numerator 
      h_other_denom_CV[i][j] = (TH1D*)file_CV_eff->Get(Form("h_other_eff_denom%s_0",other_var[i])); //CV denominator
      
    }
  }


  /////////////////////////////////////////
  //Now we have to grab the universe values
  /////////////////////////////////////////
  TFile* file = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/detVar/histograms_pelee_xsec_%s.root",input_file));
  TFile* file_efficiency = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/detVar/histograms_%s_mc_eff.root",input_file));

  //Values are stored in _0      
  for(int j=0; j < num_var; j++){
    h_muon_total[0][j] = (TH1D*)file->Get(Form("h_muon%s_total_0",var0[j]));
    h_muon_cc2p[0][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_0",var0[j]));
    h_muon_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_0",var[j])); //CV matrix
    h_muon_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_0",var[j])); //CV numerator 
    h_muon_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_0",var[j])); //CV denominator

    h_leading_total[0][j] = (TH1D*)file->Get(Form("h_leading%s_total_0",var0[j]));
    h_leading_cc2p[0][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_0",var0[j]));
    h_leading_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_0",var[j])); //CV matrix
    h_leading_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_0",var[j])); //CV numerator 
    h_leading_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_0",var[j])); //CV denominator

    h_recoil_total[0][j] = (TH1D*)file->Get(Form("h_recoil%s_total_0",var0[j]));
    h_recoil_cc2p[0][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_0",var0[j]));
    h_recoil_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_0",var[j])); //CV matrix
    h_recoil_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_0",var[j])); //CV numerator 
    h_recoil_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_0",var[j])); //CV denominator
    
  }
  
  for(int j=0; j < num_other_var; j++){
    h_other_total[0][j] = (TH1D*)file->Get(Form("h%s_total_0",other_var[j]));
    h_other_cc2p[0][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_0",other_var[j]));
    h_other_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_0",other_var[j])); //CV matrix
    h_other_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_0",other_var[j])); //CV numerator 
    h_other_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_0",other_var[j])); //CV denominator
  }
  
 std::cout<<"[FINISHED GRAABBING ALL THE FILES]"<<std::endl;
  
}//end of grab_histograms


TH2D* other_detvar::make_covariance_matrix(TH1D* h_CV, std::vector<TH1D*> h_universe, int num_universes,const char* name){
  
  ////////////////////////////////////////////////////
  //First let us graab nu_CV and the universe value
  ////////////////////////////////////////////////////
  int nbins = h_CV->GetNbinsX();
  double nu_CV[nbins]; //array of length nbins that will contain the central values
  double nu_universe[nbins][num_universes]; //array of length nbins that will contain all of the universe values

  for(int reco_bin = 1; reco_bin <nbins+1; reco_bin++){
    double CV_value = (h_CV->GetBinContent(reco_bin)); //Total Number of MC events in Reco Bin
    nu_CV[reco_bin-1] =  CV_value;
    
    for(int u = 0; u < num_universes; u++){
      double universe_value = h_universe[u]->GetBinContent(reco_bin); //Total Number of MC events in Reco Bin
      nu_universe[reco_bin-1][u] = universe_value;
    }
  }
  
  ////////////////////////////////////
  //Now to Make our Covariance Matrix
  ///////////////////////////////////
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = h_CV->GetBinLowEdge(i+1);
  }

  TH2D* covariance_matrix = new TH2D(Form("h_2D_covariance%s",name),Form("h_2D_covariance%s",name),nbins,edges,nbins,edges);
  for(int a = 1; a < nbins+1; a++){
    for(int b=1; b < nbins+1; b++){
	double sum = 0;

	for(int u = 0; u < num_universes; u++){
	  double x = nu_CV[a-1] - nu_universe[a-1][u];
	  double y = nu_CV[b-1] - nu_universe[b-1][u];
	  double product = x * y;
	  sum += product;
	}
	double value = std::pow(num_universes,-1) * sum;
	covariance_matrix->SetBinContent(a,b,value);

    } //end of loop over b                                                                                                                                                                
  } //end of loop over a                                                                                                                                                                  

  return covariance_matrix;

} //end of make covariance matrix

void other_detvar::plot_covariance(const char* variable, const char* title, TCanvas* canv, TH2D* hist, TH1D* hist_CV, const char* directory){

  canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  hist->Draw("colz text");
  
  hist->SetTitle(Form("Covariance Matrix: %s",title)); //title
  hist->SetXTitle("Reco. Bin a");
  hist->GetXaxis()->SetTitleSize(40);
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(30);
  hist->GetXaxis()->SetTickSize(0);
  hist->SetYTitle("Reco. Bin b");
  hist->GetYaxis()->SetTitleSize(40);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(30);
  hist->GetYaxis()->SetTickSize(0);

  int nbins = hist->GetNbinsX();   
  for(int i=1; i < nbins+1; i++){
    hist->GetXaxis()->SetBinLabel(i,Form("%d",i));
    hist->GetYaxis()->SetBinLabel(i,Form("%d",i));
  }
  
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist_CV->GetBinLowEdge(i+1);
  }

  TLine* lx;
  TLine* ly;

  for (int bin = 0; bin < nbins+1; bin++) {
    lx = new TLine(edges[bin],edges[0],edges[bin],edges[nbins]);
    lx->Draw("SAME");
    lx->SetLineColor(kBlack);
    lx->SetLineWidth(2);
    ly = new TLine(edges[0],edges[bin],edges[nbins],edges[bin]);
    ly->Draw("SAME");
    ly->SetLineColor(kBlack);
    ly->SetLineWidth(2);  
  }
  gStyle->SetOptStat(0);
  canv->Print(Form("images/%s/Covariance_Matrices/_2D%s.png",directory,variable));
}

TH1D* other_detvar::make_error_histogram(const char* name, TH1D* hist_CV,TH2D* h_covariance){

  int nbins = hist_CV->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist_CV->GetBinLowEdge(i+1);
  }

  TH1D* hist_errors = new TH1D(Form("hist_errors%s",name),Form("hist_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    for(int j=1; j < nbins+1; j++){
      if(i == j){
	double bin_content = h_covariance->GetBinContent(i,j);
	double value = std::pow(bin_content,0.5);
	hist_errors->SetBinContent(i,value);
      }
    }
  }

  return hist_errors;
  
} //end of make_error_histograms

TH1D* other_detvar::make_fractional_error_histogram(const char* name, TH1D* hist_CV,TH1D* hist_error){

  int nbins = hist_CV->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist_CV->GetBinLowEdge(i+1);
  }

  TH1D* hist_fractional_errors = new TH1D(Form("hist_fractional_errors%s",name),Form("hist_fractional_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    double bin_content_error = hist_error->GetBinContent(i);
    double bin_content_CV = hist_CV->GetBinContent(i);
    double value = bin_content_error / bin_content_CV;
    hist_fractional_errors->SetBinContent(i,value);
  }
  
  return hist_fractional_errors;
}//end of make_fractional_error_histograms

#endif
