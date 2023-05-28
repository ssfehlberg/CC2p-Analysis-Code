#define all_plot_cxx
#include "plotting.h"

void all_plot(){

  //Plotting class
  ////////////////
  plotting_tools plot;
  
  //Files we are going to need for the Dirt and EXT Error
  ///////////////////////////////////////////////////////
  TFile* file_cv = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));
   
  static const int num_files = 5;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt"};//,"Statistical"};
  
  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon_total[num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
  bool flip_legend_muon[num_var] = {true,false,false};

  //Leading
  TH1D* h_leading[num_files][num_var];
  TH1D* h_leading_total[num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH1D* h_recoil[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};
  
  /*static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
                                                 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};
  bool flip_legend_other[num_other_var] = {true,false,false,false,false,false,false};*/

    static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};
  
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];

  //Grab Hstograms
  ////////////////
  for(int f=0; f < num_files; f++){

    file[f] =  new TFile(Form("../root_files/%s/systematics.root",file_names[f]));
    
    for(int j=0; j < num_var; j++){
      h_muon[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_muon%s",var0[j]));
      h_leading[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_leading%s",var0[j]));
      h_recoil[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j]));
    }

    for(int j=0; j< num_other_var; j++){
      h_other[f][j] =  (TH1D*)file[f]->Get(Form("hist_fractional_errors%s",other_var[j]));
    }
  }

  h_muon[0][0]->Draw("hist");
  
  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/total_error.root"),"RECREATE"); //output root file      

  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){
    
    std::vector<TH1D*> muon;
    std::vector<TH1D*> leading;
    std::vector<TH1D*> recoil;

    for(int f=0; f < num_files; f++){
      muon.push_back(h_muon[f][j]);
      leading.push_back(h_leading[f][j]);
      recoil.push_back(h_recoil[f][j]);
    }

    plot.plot_total_error(Form("_muon%s",var0[j]),canv_muon[j],flip_legend_muon[j],legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    plot.plot_total_error(Form("_leading%s",var0[j]),canv_leading[j],flip_legend_leading[j],legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    plot.plot_total_error(Form("_recoil%s",var0[j]),canv_recoil[j],flip_legend_recoil[j],legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    
  }

  for(int j=0; j < num_other_var; j++){

    std::vector<TH1D*> other;
    for(int f=0; f < num_files; f++){
      other.push_back(h_other[f][j]);
    }

    plot.plot_total_error(Form("%s",other_var[j]),canv_other[j],flip_legend_other[j],legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
  }

  tfile_mine->Close();
  
}
