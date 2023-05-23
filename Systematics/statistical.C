#define statistical_cxx
#include "tools/helpers.h"
#include "tools/plotting.h"

TH2D* make_stat_covariance(const char* name, TH1D* h_mc_tot, TH1D* h_bnb,TH1D* h_ext, TH1D* h_dirt){
  
  int nbins = h_mc_tot->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = h_mc_tot->GetBinLowEdge(i+1);
  }

  TH2D* covariance_matrix = new TH2D(Form("h_2D_covariance%s",name),Form("h_2D_covariance%s",name),nbins,edges,nbins,edges);

  for(int i=1; i < nbins+1; i++){
    for(int j=1; j < nbins+1; j++){
      if(i==j){
	double mc = h_mc_tot->GetBinContent(i);
	double bnb = 0; h_bnb->GetBinContent(i);
	double ext = h_ext->GetBinContent(i);
	double dirt = h_dirt->GetBinContent(i);
	double value = mc + bnb + ext + dirt;
	covariance_matrix->SetBinContent(i,j,value);
	
      } else {
	covariance_matrix->SetBinContent(i,j,0);
      }
    }
  }

  return covariance_matrix;
  
}


void statistical(){


  const char* directory_name = "Statistical";
  
  //Have to define the histograms cause I did this dumbly
  //////////////////////////////////////////////////////
  TFile* file_cv = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));
  TFile* file_bnb = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root "));


  //DEFINE BNB SHIT IN THE CONSTANTS FILE
  
  for(int j=0; j < num_var; j++){

    h_muon_CV[j] = (TH1D*)file_cv->Get(Form("h_muon%s_total",var0[j]));
    h_muon_BNB[j] = (TH1D*)file_bnb->Get(Form("h_muon%s_bnb",var0[j]));
    h_muon_EXT[j] = (TH1D*)file_ext->Get(Form("h_muon%s_ext",var0[j]));
    h_muon_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_muon%s_dirt_wgt",var0[j]));

    h_leading_CV[j] = (TH1D*)file_cv->Get(Form("h_leading%s_total",var0[j]));
    h_leading_BNB[j] = (TH1D*)file_bnb->Get(Form("h_leading%s_bnb",var0[j]));
    h_leading_EXT[j] = (TH1D*)file_ext->Get(Form("h_leading%s_ext",var0[j]));
    h_leading_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_leading%s_dirt_wgt",var0[j]));

    h_recoil_CV[j] = (TH1D*)file_cv->Get(Form("h_recoil%s_total",var0[j]));
    h_recoil_BNB[j] = (TH1D*)file_bnb->Get(Form("h_recoil%s_bnb",var0[j]));
    h_recoil_EXT[j] = (TH1D*)file_ext->Get(Form("h_recoil%s_ext",var0[j]));
    h_recoil_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[j]));

  }

  for(int j=0; j < num_other_var; j++){
    h_other_CV[j] = (TH1D*)file_cv->Get(Form("h%s_total",other_var[j]));
    h_other_BNB[j] = (TH1D*)file_bnb->Get(Form("h%s_bnb",other_var[j]));
    h_other_EXT[j] = (TH1D*)file_ext->Get(Form("h%s_ext",other_var[j]));
    h_other_Dirt[j] = (TH1D*)file_dirt->Get(Form("h%s_dirt_wgt",other_var[j]));
  }

  //File to save shit to
  TFile *tfile_mine = new TFile(Form("root_files/%s/systematics.root",directory_name),"RECREATE"); //output root file   
  
  for(int j=0; j < num_var; j++){

    //Muon
    //////////////
    
    //Make covariance matrix
    h_2D_muon[j] = make_stat_covariance(Form("_muon_%s",var0[j]), h_muon_CV[j],h_muon_BNB[j],h_muon_EXT[j],h_muon_Dirt[j]);
    plot_covariance(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_2D_muon[j],h_2D_muon[j],h_muon_CV[j],directory_name);
    h_2D_muon[j]->Write();
    
    //Getting the Error
    h_muon_error[j] = make_error_histogram(Form("_muon%s",var0[j]),h_muon_CV[j],h_2D_muon[j]);
    plot_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_error[j],h_muon_error[j],directory_name);
    h_muon_error[j]->Write();
    
    //Getting the fractional error
    h_muon_fractional_error[j] = make_fractional_error_histogram(Form("_muon%s",var0[j]),h_muon_CV[j],h_muon_error[j]);
    plot_fractional_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_fractional_error[j],h_muon_fractional_error[j],directory_name);
    h_muon_fractional_error[j]->Write();

    //Leading
    //////////////                                                                                                                                                                                                             
    //Make covariance matrix
    h_2D_leading[j] = make_stat_covariance(Form("_leading_%s",var0[j]), h_leading_CV[j],h_leading_BNB[j],h_leading_EXT[j],h_leading_Dirt[j]);
    plot_covariance(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_2D_leading[j],h_2D_leading[j],h_leading_CV[j], directory_name);
    h_2D_leading[j]->Write();
    
    //Getting the Error
    h_leading_error[j] = make_error_histogram(Form("_leading%s",var0[j]),h_leading_CV[j],h_2D_leading[j]);
    plot_error(Form("_leading%s",var0[j]), Form("Leading Proton %s",var_titles[j]),canv_leading_error[j],h_leading_error[j],directory_name);
    h_leading_error[j]->Write();

    //Getting the fractional error
    h_leading_fractional_error[j] = make_fractional_error_histogram(Form("_leading%s",var0[j]),h_leading_CV[j],h_leading_error[j]);
    plot_fractional_error(Form("_leading%s",var0[j]), Form("Leading Proton %s",var_titles[j]),canv_leading_fractional_error[j],h_leading_fractional_error[j],directory_name);
    h_leading_fractional_error[j]->Write();

    //Recoil
    //////////////                                                                                                                                                                                                             
    //Make covariance matrix
    h_2D_recoil[j] = make_stat_covariance(Form("_recoil_%s",var0[j]), h_recoil_CV[j],h_recoil_BNB[j],h_recoil_EXT[j],h_recoil_Dirt[j]);
    plot_covariance(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_2D_recoil[j],h_2D_recoil[j],h_recoil_CV[j],directory_name);
    h_2D_recoil[j]->Write();
    
    //Getting the Error
    h_recoil_error[j] = make_error_histogram(Form("_recoil%s",var0[j]),h_recoil_CV[j],h_2D_recoil[j]);
    plot_error(Form("_recoil%s",var0[j]), Form("Recoil Proton %s",var_titles[j]),canv_recoil_error[j],h_recoil_error[j],directory_name);
    h_recoil_error[j]->Write();

    //Getting the fractional error
    h_recoil_fractional_error[j] = make_fractional_error_histogram(Form("_recoil%s",var0[j]),h_recoil_CV[j],h_recoil_error[j]);
    plot_fractional_error(Form("_recoil%s",var0[j]), Form("Recoil Proton %s",var_titles[j]),canv_recoil_fractional_error[j],h_recoil_fractional_error[j],directory_name);
    h_recoil_fractional_error[j]->Write();

  }

  //Other Variables
  ////////////////
  
  for(int j=0; j < num_other_var; j++){

    //Make covariance matrix
    h_2D_other[j] = make_stat_covariance(Form("%s",other_var[j]), h_other_CV[j],h_other_BNB[j],h_other_EXT[j],h_other_Dirt[j]);
    plot_covariance(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_2D_other[j],h_2D_other[j],h_other_CV[j],directory_name);
    h_2D_other[j]->Write();
      
    //Getting the Error
    h_other_error[j] = make_error_histogram(Form("%s",other_var[j]),h_other_CV[j],h_2D_other[j]);
    plot_error(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_other_error[j],h_other_error[j],directory_name);
    h_other_error[j]->Write();

    //Getting the fractional error
    h_other_fractional_error[j] = make_fractional_error_histogram(Form("%s",other_var[j]),h_other_CV[j],h_other_error[j]);
    plot_fractional_error(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_other_fractional_error[j],h_other_fractional_error[j],directory_name);
    h_other_fractional_error[j]->Write();

  }

  tfile_mine->Close();

  
}


