#define all_plot_cxx
#include "tools/plotting.h"

void plot_GENIE_total(){

  //Stuff for Drawing
  //////////////////////
  Color_t colors[] = {kRed, kRed, kOrange+6, kOrange+6, kYellow-3, kYellow-3, kGreen+2 ,kGreen+2, kBlue, kBlue, kViolet+1, kViolet+1,kMagenta};
  int line_style[] = {1,9,1,9,1,9,1,9,1,9,1,9,1};

  
  static const int num_files = 13;
  TFile* file[num_files];
  const char* file_names[num_files] = {"All_UBGenie","AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
				       "NormNCCOH_UBGenie","RPA_CCQE_UBGenie",
				       "ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
				       "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie",
				       "xsr_scc_Fa3_SCC","xsr_scc_Fv3_SCC","RootinoFix_UBGenie"};

  const char* file_titles[num_files] = {"GENIE Multisims", "AxFFCCQE Shape","DecayAngleMEC","NormCCCOH",
					"NormNCCOH","RPA_CCQE",
					"ThetaDelta2NRad","Theta_Delta2Npi",
					"VecFFCCQEshape","XSecShape_CCMEC",
					"SCC: F^{3}_{V}","SCC: F^{3}_{A}","Rootino Fix"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon_total[num_var];
  TH2D* h_muon_covariance[num_files][num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
  bool flip_legend_muon[num_var] = {true,false,false};
  
  //Leading
  TH1D* h_leading[num_files][num_var];
  TH1D* h_leading_total[num_var];
  TH2D* h_leading_covariance[num_files][num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH1D* h_recoil[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  TH2D* h_recoil_covariance[num_files][num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};
  
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TH2D* h_other_covariance[num_files][num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};
  
  //Grab Files
  ////////////////
  for(int f=0; f < num_files; f++){
    file[f] =  new TFile(Form("root_files/%s/systematics.root",file_names[f]));

    for(int j=0; j < num_var; j++){
      h_muon[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_muon%s",var0[j]));
      h_muon_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancemuon_%s",var0[j]));
      
      h_leading[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_leading%s",var0[j]));
      h_leading_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covarianceleading_%s",var0[j]));
       
      h_recoil[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j]));
      h_recoil_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancerecoil_%s",var0[j]));
      
    }

    for(int j=0; j< num_other_var; j++){
      h_other[f][j] =  (TH1D*)file[f]->Get(Form("hist_fractional_errors%s",other_var[j]));
      h_other_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariance%s",other_var[j]));
    }
  }


  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("root_files/Genie_Total/systematics.root"),"RECREATE"); //output root file 

  
  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){

    std::vector<TH1D*> muon;
    std::vector<TH2D*> muon_2D;
    std::vector<TH1D*> leading;
    std::vector<TH2D*> leading_2D;
    std::vector<TH1D*> recoil;
    std::vector<TH2D*> recoil_2D;

    for(int f=0; f< num_files; f++){
      muon.push_back(h_muon[f][j]);
      muon_2D.push_back(h_muon_covariance[f][j]);
      plot_covariance(Form("_muon%s",var0[j]),Form("%s: Muon %s",file_titles[f],var_titles[j]),h_muon_covariance[f][j],file_names[f]);

      leading.push_back(h_leading[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      plot_covariance(Form("_leading%s",var0[j]),Form("%s: Leading Proton %s",file_titles[f],var_titles[j]),h_leading_covariance[f][j],file_names[f]);

      recoil.push_back(h_recoil[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]);
      plot_covariance(Form("_recoil%s",var0[j]),Form("%s: Recoil Proton %s",file_titles[f],var_titles[j]),h_recoil_covariance[f][j],file_names[f]);

    }

    plot_GENIE_total_error(Form("_muon%s",var0[j]),canv_muon[j],flip_legend_muon[j],legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    make_total_covariance(Form("_muon%s",var0[j]),Form("Total GENIE Error: Muon %s",var_titles[j]),muon_2D,"GENIE_Total");

    plot_GENIE_total_error(Form("_leading%s",var0[j]),canv_leading[j],flip_legend_leading[j],legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    make_total_covariance(Form("_leading%s",var0[j]),Form("Total GENIE Error: Leading Proton %s",var_titles[j]),leading_2D,"GENIE_Total");

    plot_GENIE_total_error(Form("_recoil%s",var0[j]),canv_recoil[j],flip_legend_recoil[j],legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    make_total_covariance(Form("_recoil%s",var0[j]),Form("Total GENIE Error: Recoil Proton %s",var_titles[j]),recoil_2D,"GENIE_Total");
    
  }

  for(int j=0; j < num_other_var; j++){

    std::vector<TH1D*> other;
    std::vector<TH2D*> other_2D;

    for(int f=0; f < num_files; f++){
      other.push_back(h_other[f][j]);
      other_2D.push_back(h_other_covariance[f][j]);
      plot_covariance(Form("%s",other_var[j]),Form("%s: %s",file_titles[f],other_var_titles[j]),h_other_covariance[f][j],file_names[f]);
    }

    plot_GENIE_total_error(Form("%s",other_var[j]),canv_other[j],flip_legend_other[j],legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
    make_total_covariance(Form("%s",other_var[j]),Form("Total GENIE Error: %s",other_var_titles[j]),other_2D,"GENIE_Total"); 
  }

  tfile_mine->Close();
  
}

void all_plot(){

  //plot_GENIE_total();
  std::cout<<"Finished Running Plot Genie Total."<<std::endl;
  
  //Files we are going to need for the Dirt and EXT Error
  ///////////////////////////////////////////////////////
  TFile* file_cv = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));
   
  static const int num_files = 7;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt","Iteration","Statistical"};
  const char* file_titles[num_files] = {"GENIE Total","Detector Variations","Flux Multisims","G4 Multisims","Dirt (Sys.)","Iteration","Statistical Error"};
  
  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon_total[num_var];
  TH2D* h_muon_covariance[num_files][num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
  bool flip_legend_muon[num_var] = {true,false,false};

  //Leading
  TH1D* h_leading[num_files][num_var];
  TH1D* h_leading_total[num_var];
  TH2D* h_leading_covariance[num_files][num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH1D* h_recoil[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  TH2D* h_recoil_covariance[num_files][num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};
  
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{p_{L},p_{R}})","cos(#gamma_{#mu,p_{sum}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TH2D* h_other_covariance[num_files][num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};

  //Grab Hstograms
  ////////////////
  for(int f=0; f < num_files; f++){
    file[f] =  new TFile(Form("root_files/%s/systematics.root",file_names[f]));

    for(int j=0; j < num_var; j++){
      h_muon[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_muon%s",var0[j]));
      h_muon_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancemuon_%s",var0[j]));
      
      h_leading[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_leading%s",var0[j]));
      h_leading_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covarianceleading_%s",var0[j]));
       
      h_recoil[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j]));
      h_recoil_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancerecoil_%s",var0[j]));
    }

    for(int j=0; j< num_other_var; j++){
      h_other[f][j] =  (TH1D*)file[f]->Get(Form("hist_fractional_errors%s",other_var[j]));
      h_other_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariance%s",other_var[j]));
    }
  }

  
  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("root_files/total_error.root"),"RECREATE"); //output root file      

  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){
    std::cout<<"Value of j: "<<j<<std::endl;
    
    std::vector<TH1D*> muon;
    std::vector<TH2D*> muon_2D;
    std::vector<TH1D*> leading;
    std::vector<TH2D*> leading_2D;
    std::vector<TH1D*> recoil;
    std::vector<TH2D*> recoil_2D;

    for(int f=0; f< num_files; f++){
      std::cout<<"Value of f: "<<f<<std::endl;
      muon.push_back(h_muon[f][j]);
      leading.push_back(h_leading[f][j]);
      recoil.push_back(h_recoil[f][j]);

      muon_2D.push_back(h_muon_covariance[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]);
      plot_covariance(Form("_muon%s",var0[j]),Form("%s: Muon %s",file_titles[f],var_titles[j]),h_muon_covariance[f][j],file_names[f]);
      plot_covariance(Form("_leading%s",var0[j]),Form("%s: Leading Proton %s",file_titles[f],var_titles[j]),h_leading_covariance[f][j],file_names[f]);
      plot_covariance(Form("_recoil%s",var0[j]),Form("%s: Recoil Proton %s",file_titles[f],var_titles[j]),h_recoil_covariance[f][j],file_names[f]);
    }

    plot_total_error(Form("_muon%s",var0[j]),canv_muon[j],flip_legend_muon[j],legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    make_total_covariance(Form("_muon%s",var0[j]),Form("Total Error: Muon %s",var_titles[j]),muon_2D,"Total_Error");
    
    plot_total_error(Form("_leading%s",var0[j]),canv_leading[j],flip_legend_leading[j],legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    make_total_covariance(Form("_leading%s",var0[j]),Form("Total Error: Leading Proton %s",var_titles[j]),leading_2D,"Total_Error");
    
    plot_total_error(Form("_recoil%s",var0[j]),canv_recoil[j],flip_legend_recoil[j],legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    make_total_covariance(Form("_recoil%s",var0[j]),Form("Total Error: Recoil Proton %s",var_titles[j]),recoil_2D,"Total_Error");
  }
  
  for(int j=0; j < num_other_var; j++){
    std::vector<TH1D*> other;
    std::vector<TH2D*> other_2D;

    for(int f=0; f < num_files; f++){
      other.push_back(h_other[f][j]);
      other_2D.push_back(h_other_covariance[f][j]);
      plot_covariance(Form("%s",other_var[j]),Form("%s: %s",file_titles[f],other_var_titles[j]),h_other_covariance[f][j],file_names[f]);      
    }
    
    plot_total_error(Form("%s",other_var[j]),canv_other[j],flip_legend_other[j],legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
    make_total_covariance(Form("%s",other_var[j]),Form("Total Error: %s",other_var_titles[j]),other_2D,"Total_Error");
  }

  tfile_mine->Close();
  
}
