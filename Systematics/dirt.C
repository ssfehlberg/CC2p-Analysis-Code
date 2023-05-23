#define dirt_cxx
#include "tools/helpers.h"
#include "tools/plotting.h"

void dirt(){

  helper help; //helper class
  int universes =2;
  const char* directory_name = {"Dirt"};

  //Have to define the histograms cause I did this dumbly
  //////////////////////////////////////////////////////
  TFile* file_cv = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));


  for(int j=0; j < num_var; j++){

    h_muon_CV[j] = (TH1D*)file_cv->Get(Form("h_muon%s_total",var0[j]));
    h_muon_EXT[j] = (TH1D*)file_ext->Get(Form("h_muon%s_ext",var0[j]));
    h_muon_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_muon%s_dirt_wgt",var0[j]));

    h_leading_CV[j] = (TH1D*)file_cv->Get(Form("h_leading%s_total",var0[j]));
    h_leading_EXT[j] = (TH1D*)file_ext->Get(Form("h_leading%s_ext",var0[j]));
    h_leading_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_leading%s_dirt_wgt",var0[j]));

    h_recoil_CV[j] = (TH1D*)file_cv->Get(Form("h_recoil%s_total",var0[j]));
    h_recoil_EXT[j] = (TH1D*)file_ext->Get(Form("h_recoil%s_ext",var0[j]));
    h_recoil_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[j]));

  }

  for(int j=0; j < num_other_var; j++){
    h_other_CV[j] = (TH1D*)file_cv->Get(Form("h%s_total",other_var[j]));
    h_other_EXT[j] = (TH1D*)file_ext->Get(Form("h%s_ext",other_var[j]));
    h_other_Dirt[j] = (TH1D*)file_dirt->Get(Form("h%s_dirt_wgt",other_var[j]));
  }


  //Here is where we make the dirt universes:
  //[0] -> Means CV (i.e. 100%)
  //[1] -> Means Universe with 130%
  //////////////////////////////////////////
  for(int j=0; j < num_var; j++){
    h_muon[0][j] = (TH1D*)h_muon_Dirt[j]->Clone();
    h_muon[1][j] = (TH1D*)h_muon_Dirt[j]->Clone();
    h_muon[1][j]->Scale(1.3);

    h_leading[0][j] = (TH1D*)h_leading_Dirt[j]->Clone();
    h_leading[1][j] = (TH1D*)h_leading_Dirt[j]->Clone();
    h_leading[1][j]->Scale(1.3);

    h_recoil[0][j] = (TH1D*)h_recoil_Dirt[j]->Clone();
    h_recoil[1][j] = (TH1D*)h_recoil_Dirt[j]->Clone();
    h_recoil[1][j]->Scale(1.3);
  }

  for(int j=0; j < num_other_var; j++){
    h_other[0][j] = (TH1D*)h_other_Dirt[j]->Clone();
    h_other[1][j] = (TH1D*)h_other_Dirt[j]->Clone();
    h_other[1][j]->Scale(1.3);
  }
    

  ///////////////////////////
  //Now to fucking plot shit
  //////////////////////////

  //File to save shit to
  TFile *tfile_mine = new TFile(Form("root_files/%s/systematics.root",directory_name),"RECREATE"); //output root file   
  
    
  for(int j=0; j < num_var; j++){
    
    //Muon
    ///////////////////////////////////
    std::vector<TH1D*> h_muon_universes;
    for(int k=0; k < universes; k++){
      h_muon_universes.push_back(h_muon[k][j]);
    }

    //multisim
    plot_multisims(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]), canv_muon[j], legend_muon[j], h_muon_universes, h_muon_CV[j], ymin_muon[j], ymax_muon[j], universes, directory_name);
    
    //Covariance shit
    h_2D_muon[j] = make_covariance_matrix(h_muon_CV[j], h_muon_EXT[j], h_muon_Dirt[j],h_muon_universes, universes, Form("muon_%s",var0[j]),true);
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
    
    //Leading Proton
    ////////////////
    std::vector<TH1D*> h_leading_universes;
    for(int k=0; k < universes; k++){
      h_leading_universes.push_back(h_leading[k][j]);
    }
      
    //multisim
    plot_multisims(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]), canv_leading[j], legend_leading[j], h_leading_universes, h_leading_CV[j], ymin_leading[j], ymax_leading[j], universes,directory_name);
    
    //Covariance shit
    h_2D_leading[j] = make_covariance_matrix(h_leading_CV[j], h_leading_EXT[j], h_leading_Dirt[j],h_leading_universes, universes, Form("leading_%s",var0[j]),true);
    plot_covariance(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_2D_leading[j],h_2D_leading[j],h_leading_CV[j], directory_name);
    h_2D_leading[j]->Write();
    
    //Getting the Error
    h_leading_error[j] = make_error_histogram(Form("_leading%s",var0[j]),h_leading_CV[j],h_2D_leading[j]);
    plot_error(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_leading_error[j],h_leading_error[j],directory_name);
    h_leading_error[j]->Write();
    
    //Getting the fractional error
    h_leading_fractional_error[j] = make_fractional_error_histogram(Form("_leading%s",var0[j]),h_leading_CV[j],h_leading_error[j]);
    plot_fractional_error(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_leading_fractional_error[j],h_leading_fractional_error[j],directory_name);
    h_leading_fractional_error[j]->Write();
      
    //Recoil Proton
    ///////////////
    std::vector<TH1D*> h_recoil_universes;
    for(int k=0; k < universes; k++){
      h_recoil_universes.push_back(h_recoil[k][j]);
    }
      
    //multisim
    plot_multisims(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]), canv_recoil[j], legend_recoil[j], h_recoil_universes, h_recoil_CV[j], ymin_recoil[j], ymax_recoil[j], universes,directory_name);
    
    //Covariance shit
    h_2D_recoil[j] = make_covariance_matrix(h_recoil_CV[j], h_recoil_EXT[j], h_recoil_Dirt[j],h_recoil_universes, universes, Form("recoil_%s",var0[j]),true);
    plot_covariance(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_2D_recoil[j],h_2D_recoil[j],h_recoil_CV[j],directory_name);
    h_2D_recoil[j]->Write();
    
    //Getting the Error
    h_recoil_error[j] = make_error_histogram(Form("_recoil%s",var0[j]),h_recoil_CV[j],h_2D_recoil[j]);
    plot_error(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_recoil_error[j],h_recoil_error[j],directory_name);
    h_recoil_error[j]->Write();
    
    //Getting the fractional error
    h_recoil_fractional_error[j] = make_fractional_error_histogram(Form("_recoil%s",var0[j]),h_recoil_CV[j],h_recoil_error[j]);
    plot_fractional_error(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_recoil_fractional_error[j],h_recoil_fractional_error[j],directory_name);
    h_recoil_fractional_error[j]->Write();
      
  } //end of loop over particles
    
  for(int j=0; j <num_other_var; j++){
    std::vector<TH1D*> h_other_universes;
    for(int k=0; k < universes; k++){
      h_other_universes.push_back(h_other[k][j]);
    }
    
    //multisim
    plot_multisims(Form("%s",other_var[j]), Form("%s",other_var_titles[j]), canv_other[j], legend_other[j], h_other_universes, h_other_CV[j], ymin_other[j], ymax_other[j], universes,directory_name);
    
    //Covariance shit
    h_2D_other[j] = make_covariance_matrix(h_other_CV[j], h_other_EXT[j], h_other_Dirt[j],h_other_universes, universes, Form("%s",other_var[j]),true);
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
    
      
  } //end of loop over other var
    
  tfile_mine->Close();
  
}//end of program
