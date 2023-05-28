#include "plotting.h"

class detvar_plot{

 public:
  virtual void main();
 private:

  virtual void Grab_Histograms();
  
  static const int num_files = 9;
  TFile* file[num_files];
  const char* file_names[num_files] = {"detVar_LY_Attenuation","detVar_LY_Down","detVar_LY_Rayleigh",
				       "detVar_ThetaXZ", "detVar_ThetaYZ",
				       "detVar_X","detVar_YZ",
				       "detVar_Recombination","detVar_SCE"};

  const char* sample_titles[num_files] = {"LY Attenuation","LY Down","LY Rayleigh",
                                          "ThetaXZ","ThetaYZ",
                                          "X","YZ",
                                          "Recombination","SCE"};

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
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TH2D* h_other_covariance[num_files][num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};

}; //end of class definition

void detvar_plot::Grab_Histograms(){

  //Grab Files
  ////////////////
  for(int f=0; f < num_files; f++){
    file[f] =  new TFile(Form("../root_files/detVar/systematics_%s.root",file_names[f]));
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

} //end of grab histograms

void detvar_plot::main(){


  plotting_tools plot;

  //Grab Files
  ////////////////
  Grab_Histograms();

  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/detVar/systematics.root"),"RECREATE"); //output root file 
  
  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){

    std::vector<TH1D*> muon;
    std::vector<TH2D*> muon_2D;
    std::vector<TH1D*> leading;
    std::vector<TH2D*> leading_2D;
    std::vector<TH1D*> recoil;
    std::vector<TH2D*> recoil_2D;

    for(int f=0; f < num_files; f++){
      muon.push_back(h_muon[f][j]);
      muon_2D.push_back(h_muon_covariance[f][j]);
      leading.push_back(h_leading[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      recoil.push_back(h_recoil[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]);
    }

    plot.plot_detvar_total_error(Form("_muon%s",var0[j]),canv_muon[j],flip_legend_muon[j],legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    plot.make_total_covariance(Form("_muon%s",var0[j]),Form("Total DetVar Error: Muon %s",var_titles[j]),muon_2D,"detVar/Total");

    plot.plot_detvar_total_error(Form("_leading%s",var0[j]),canv_leading[j],flip_legend_leading[j],legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    plot.make_total_covariance(Form("_leading%s",var0[j]),Form("Total DetVar Error: Leading Proton %s",var_titles[j]),leading_2D,"detVar/Total");
    
    plot.plot_detvar_total_error(Form("_recoil%s",var0[j]),canv_recoil[j],flip_legend_recoil[j],legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    plot.make_total_covariance(Form("_recoil%s",var0[j]),Form("Total DetVar Error: Recoil Proton %s",var_titles[j]),recoil_2D,"detVar/Total");
    
  }
  
  for(int j=0; j < num_other_var; j++){

    std::vector<TH1D*> other;
    std::vector<TH2D*> other_2D;
    
    for(int f=0; f<num_files; f++){
      other.push_back(h_other[f][j]);
      other_2D.push_back(h_other_covariance[f][j]);
    }

    plot.plot_detvar_total_error(Form("%s",other_var[j]),canv_other[j],flip_legend_other[j],legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
    plot.make_total_covariance(Form("%s",other_var[j]),Form("Total DetVar Error: %s",other_var_titles[j]),other_2D,"detVar/Total");
  }
  
  tfile_mine->Close();

} //end of program
