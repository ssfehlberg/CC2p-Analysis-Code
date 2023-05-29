//////////////////////////////////////////////////////////////////////////////////////
// This .C file is meant to process the detector variation samples in a unisim fashion.
// This differs from the treatment used in detVar.C which ignores bin to bin correlations
///////////////////////////////////////////////////////////////////////////////////////
#define other_detvar_cxx
#include "other_detvar.h"

void other_detvar::main(){

  //Grab the CV histograms. Note there are 2 CV samples:
  //One for most of the samples and a smaller one for recombination + SCE
  ////////////////////////////////
  //Grab_CV();

  //Now to loop over the files
  ////////////////////////////
  for(int i=0; i < num_files; i++){
    directory_name = directory_name_list[i];
    Grab_Histograms(directory_name, universes);
    TFile *tfile_mine = new TFile(Form("root_files/detVar/systematics_%s.root",directory_name),"RECREATE"); //output root file 

    if(i == 7|| i == 8){
      which_CV = 1;
    } else {
      which_CV = 0;
    }
    
    ///////////
    //Particles
    ///////////
    for(int j=0; j< num_var; j++){

      //Muon                                                                                                                                                         
      ////////////////////////////////////                                                                                                                           
      std::vector<TH1D*> h_muon_universes;
      h_muon_xsec[0][j] =  Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_overlay_muon[j][0][which_CV],h_overlay_muon[j][1][which_CV],h_bnb_muon[j],
                                              h_muon_matrices_CV[j][which_CV],muon_iter[j],h_muon_num_CV[j][which_CV], h_muon_denom_CV[j][which_CV],
                                              muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_CV",var[j]),Form("detVar/%s",directory_name),false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
        h_muon_xsec[k+1][j] = Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_muon_total[k][j],h_muon_cc2p[k][j],h_bnb_muon[j],
						 h_muon_matrices[k][j],muon_iter[j],h_muon_num[k][j], h_muon_denom[k][j],
						 muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_%d",var[j],k),Form("detVar/%s",directory_name),false,true);
	h_muon_universes.push_back(h_muon_xsec[k+1][j]); //h_muon_xsec[0][j] should be same as h_muon_xsec_CV[j]                                                     
      }

      //multisim                                                                                                                                                     
      //plot_multisims(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]), h_muon_universes, h_muon_xsec[0][j], ymin_muon[j], muon_xsec_max[j], universes, Form("detVar/%s",directory_name));

      //Covariance shit                                                                                                                                              
      h_2D_muon[j] = make_covariance_matrix(h_muon_xsec[0][j],h_muon_universes, universes, Form("muon_%s",var0[j]));
      plot.plot_covariance(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_2D_muon[j],h_2D_muon[j],h_muon_xsec[0][j],Form("detVar/%s",directory_name));
      h_2D_muon[j]->Write();

      //Getting the Error                                                                                                                                            
      h_muon_error[j] = make_error_histogram(Form("_muon%s",var0[j]),h_muon_xsec[0][j],h_2D_muon[j]);
      plot.plot_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_error[j],h_muon_error[j],Form("detVar/%s",directory_name));
      h_muon_error[j]->Write();

      //Getting the fractional error                                                                                                                                 
      h_muon_fractional_error[j] = make_fractional_error_histogram(Form("_muon%s",var0[j]),h_muon_xsec[0][j],h_muon_error[j]);
      plot.plot_fractional_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_fractional_error[j],h_muon_fractional_error[j],Form("detVar/%s",directory_name));
      h_muon_fractional_error[j]->Write();

      //Leading                                                                                                                                                                                                                         
      ///////////////////////////////////                                                                                                                                                                                               
      std::vector<TH1D*> h_leading_universes;
      h_leading_xsec[0][j] =  Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_overlay_leading[j][0][which_CV],h_overlay_leading[j][1][which_CV],h_bnb_leading[j],
                                              h_leading_matrices_CV[j][which_CV],leading_iter[j],h_leading_num_CV[j][which_CV], h_leading_denom_CV[j][which_CV],
                                              leading_xsec_max[j],Form("True Leading Proton %s",var_titles[j]),Form("_leading%s_CV",var[j]),Form("detVar/%s",directory_name),false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
        h_leading_xsec[k+1][j] = Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_leading_total[k][j],h_leading_cc2p[k][j],h_bnb_leading[j],
						 h_leading_matrices[k][j],leading_iter[j],h_leading_num[k][j], h_leading_denom[k][j],
						 leading_xsec_max[j],Form("True Leading Proton %s",var_titles[j]),Form("_leading%s_%d",var[j],k),Form("detVar/%s",directory_name),false,true);
	h_leading_universes.push_back(h_leading_xsec[k+1][j]); //h_leading_xsec[0][j] should be same as h_leading_xsec_CV[j]                                                     
      }

      //multisim                                                                                                                                                     
      //plot_multisims(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]), h_leading_universes, h_leading_xsec[0][j], ymin_leading[j], leading_xsec_max[j], universes, Form("detVar/%s",directory_name));

      //Covariance shit                                                                                                                                              
      h_2D_leading[j] = make_covariance_matrix(h_leading_xsec[0][j],h_leading_universes, universes, Form("leading_%s",var0[j]));
      plot.plot_covariance(Form("_leading%s",var0[j]), Form("Leading Proton %s",var_titles[j]),canv_2D_leading[j],h_2D_leading[j],h_leading_xsec[0][j],Form("detVar/%s",directory_name));
      h_2D_leading[j]->Write();

      //Getting the Error                                                                                                                                            
      h_leading_error[j] = make_error_histogram(Form("_leading%s",var0[j]),h_leading_xsec[0][j],h_2D_leading[j]);
      plot.plot_error(Form("_leading%s",var0[j]), Form("Leading Proton %s",var_titles[j]),canv_leading_error[j],h_leading_error[j],Form("detVar/%s",directory_name));
      h_leading_error[j]->Write();

      //Getting the fractional error                                                                                                                                 
      h_leading_fractional_error[j] = make_fractional_error_histogram(Form("_leading%s",var0[j]),h_leading_xsec[0][j],h_leading_error[j]);
      plot.plot_fractional_error(Form("_leading%s",var0[j]), Form("Leading Proton %s",var_titles[j]),canv_leading_fractional_error[j],h_leading_fractional_error[j],Form("detVar/%s",directory_name));
      h_leading_fractional_error[j]->Write();
      
      //Recoil
      ///////////////////////////////////                                                                                                                                                                                               
      std::vector<TH1D*> h_recoil_universes;
      h_recoil_xsec[0][j] =  Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_overlay_recoil[j][0][which_CV],h_overlay_recoil[j][1][which_CV],h_bnb_recoil[j],
                                              h_recoil_matrices_CV[j][which_CV],recoil_iter[j],h_recoil_num_CV[j][which_CV], h_recoil_denom_CV[j][which_CV],
                                              recoil_xsec_max[j],Form("True Recoil Proton %s",var_titles[j]),Form("_recoil%s_CV",var[j]),Form("detVar/%s",directory_name),false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
        h_recoil_xsec[k+1][j] = Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_recoil_total[k][j],h_recoil_cc2p[k][j],h_bnb_recoil[j],
						 h_recoil_matrices[k][j],recoil_iter[j],h_recoil_num[k][j], h_recoil_denom[k][j],
						 recoil_xsec_max[j],Form("True Recoil Proton %s",var_titles[j]),Form("_recoil%s_%d",var[j],k),Form("detVar/%s",directory_name),false,true);
	h_recoil_universes.push_back(h_recoil_xsec[k+1][j]); //h_recoil_xsec[0][j] should be same as h_recoil_xsec_CV[j]                                                     
      }

      //multisim                                                                                                                                                     
      //plot_multisims(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]), h_recoil_universes, h_recoil_xsec[0][j], ymin_recoil[j], recoil_xsec_max[j], universes, Form("detVar/%s",directory_name));

      //Covariance shit                                                                                                                                              
      h_2D_recoil[j] = make_covariance_matrix(h_recoil_xsec[0][j],h_recoil_universes, universes, Form("recoil_%s",var0[j]));
      plot.plot_covariance(Form("_recoil%s",var0[j]), Form("Recoil Proton %s",var_titles[j]),canv_2D_recoil[j],h_2D_recoil[j],h_recoil_xsec[0][j],Form("detVar/%s",directory_name));
      h_2D_recoil[j]->Write();

      //Getting the Error                                                                                                                                            
      h_recoil_error[j] = make_error_histogram(Form("_recoil%s",var0[j]),h_recoil_xsec[0][j],h_2D_recoil[j]);
      plot.plot_error(Form("_recoil%s",var0[j]), Form("Recoil Proton %s",var_titles[j]),canv_recoil_error[j],h_recoil_error[j],Form("detVar/%s",directory_name));
      h_recoil_error[j]->Write();

      //Getting the fractional error                                                                                                                                 
      h_recoil_fractional_error[j] = make_fractional_error_histogram(Form("_recoil%s",var0[j]),h_recoil_xsec[0][j],h_recoil_error[j]);
      plot.plot_fractional_error(Form("_recoil%s",var0[j]), Form("Recoil Proton %s",var_titles[j]),canv_recoil_fractional_error[j],h_recoil_fractional_error[j],Form("detVar/%s",directory_name));
      h_recoil_fractional_error[j]->Write();
      
    } //end of loop over particle variables

    ////////////////
    //Other variables
    /////////////////

    std::cout<<"Right before other variables"<<std::endl;
    for(int j=0; j< num_other_var; j++){

      std::vector<TH1D*> h_other_universes;
      h_other_xsec[0][j] =  Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_overlay_other[j][0][which_CV],h_overlay_other[j][1][which_CV],h_bnb_other[j],
						 h_other_matrices_CV[j][which_CV],other_iter[j],h_other_num_CV[j][which_CV], h_other_denom_CV[j][which_CV],
						 other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("%s_CV",other_var[j]),Form("detVar/%s",directory_name),false,true); //CV Montecarlo values
      
      for(int k = 0; k < universes; k++){
        h_other_xsec[k+1][j] = Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_other_total[k][j],h_other_cc2p[k][j],h_bnb_other[j],
						    h_other_matrices[k][j],other_iter[j],h_other_num[k][j], h_other_denom[k][j],
						    other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("%s_%d",other_var[j],k),Form("detVar/%s",directory_name),false,true);
	h_other_universes.push_back(h_other_xsec[k+1][j]); //h_other_xsec[0][j] should be same as h_other_xsec_CV[j]                                                     
      }

      //multisim                                                                                                                                                     
      //plot_multisims(Form("_other%s",var0[j]), Form("Other %s",var_titles[j]), h_other_universes, h_other_xsec[0][j], ymin_other[j], other_xsec_max[j], universes, Form("detVar/%s",directory_name));
      
      //Covariance shit                                                                                                                                              
      h_2D_other[j] = make_covariance_matrix(h_other_xsec[0][j],h_other_universes, universes, Form("%s",other_var[j]));
      plot.plot_covariance(Form("%s",other_var[j]), Form(" %s",other_var_titles[j]),canv_2D_other[j],h_2D_other[j],h_other_xsec[0][j],Form("detVar/%s",directory_name));
      h_2D_other[j]->Write();

      //Getting the Error                                                                                                                                            
      h_other_error[j] = make_error_histogram(Form("_other%s",other_var[j]),h_other_xsec[0][j],h_2D_other[j]);
      plot.plot_error(Form("%s",other_var[j]), Form(" %s",other_var_titles[j]),canv_other_error[j],h_other_error[j],Form("detVar/%s",directory_name));
      h_other_error[j]->Write();

      //Getting the fractional error                                                                                                                                 
      h_other_fractional_error[j] = make_fractional_error_histogram(Form("%s",other_var[j]),h_other_xsec[0][j],h_other_error[j]);
      plot.plot_fractional_error(Form("%s",other_var[j]), Form(" %s",other_var_titles[j]),canv_other_fractional_error[j],h_other_fractional_error[j],Form("detVar/%s",directory_name));
      h_other_fractional_error[j]->Write();

    } //end loop over other
      
    tfile_mine->Close();
    
  } //end of loop over number of files
} //end of main
