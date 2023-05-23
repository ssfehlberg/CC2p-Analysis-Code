#define make_plots_cxx
#include "tools/helpers.h"
#include "tools/plotting.h"

void make_plots(){

  helper help; //helper class
  
  static const int num_files = 15;
  const char* directory_name;
  const char* directory_name_list[num_files] = {"flux_all","reint_all","All_UBGenie",
						"AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
						"NormNCCOH_UBGenie","ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
						"VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","xsr_scc_Fa3_SCC",
						"xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie","RootinoFix_UBGenie"};
  int universes;
  std::vector<int> uni = {1000,1000,500,
			  1,1,1,
			  1,1,1,
			  1,1,1,
			  1,2,1};
  bool unisim[num_files] = {false,false,false,
			    true,true,true,
			    true,true,true,
			    true,true,true,
			    true,true,true};

  for(int s = 0; s < num_files; s++){

    directory_name = directory_name_list[s];
    universes = uni[s];
    help.Grab_Histograms(directory_name,universes,unisim[s],directory_name);
  
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
      h_2D_muon[j] = make_covariance_matrix(h_muon_CV[j], h_muon_EXT[j], h_muon_Dirt[j],h_muon_universes, universes, Form("muon_%s",var0[j]));
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
      h_2D_leading[j] = make_covariance_matrix(h_leading_CV[j], h_leading_EXT[j], h_leading_Dirt[j],h_leading_universes, universes, Form("leading_%s",var0[j]));
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
      h_2D_recoil[j] = make_covariance_matrix(h_recoil_CV[j], h_recoil_EXT[j], h_recoil_Dirt[j],h_recoil_universes, universes, Form("recoil_%s",var0[j]));
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
      h_2D_other[j] = make_covariance_matrix(h_other_CV[j], h_other_EXT[j], h_other_Dirt[j],h_other_universes, universes, Form("%s",other_var[j]));
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

  } //end of loop over directories
}//end of program
