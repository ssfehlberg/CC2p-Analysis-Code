//////////////////
// S.Sword-Fehlberg, November 2021
// Extracts the Cross-Section for a variety of variables
// How to Run:
// root -b
// gSystem->Load("RooUnfold/libRooUnfold.so"); (change if your .so file is located elsewhere)
// .L xsec.C
// xsec t
// t.main()
///////////////////////////

#define xsec_cxx
#include "xsec.h"

void xsec::main(){

  //Stupid stuff for drawing
  /////////////////////////
  gStyle->SetEndErrorSize(4);
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3g");
  gStyle->SetEndErrorSize(10);
  TLatex* t;
  
  //define the classes
  ///////////////////
  closure closure_test;
  iterations iteration_test;
  MC_Comparison mc;

  //File to save the number of iterations to
  ///////////////////////////////////////////
  ofstream myfile;
  myfile.open("num_iterations.csv");
  
  //Grab the histograms
  //////////////////////
  Grab_Histograms();
  
  /////////////////////
  //Particles First
  //////////////////////
   for(int i=0; i < num_var; i++){

    
    //Muon
    /////////////

    //Make the Efficiency Histogram
    h_muon_eff[i] = make_efficiency_plot(h_muon_num[i],h_muon_denom[i],Form("True Muon %s",var_titles[i]),Form("_muon%s",var[i]));

    //Plot the matrix and make the smearing matrix
    plot_matrices(h_muon_matrices[i],Form("True Muon %s",var_titles[i]),Form("_muon%s",var[i]));
   
    //Closure test
    closure_test.Closure_Test(canv_muon_eff[i],h_muon_eff[i],h_muon_num[i],h_muon_denom[i],Form("Muon %s",var_titles[i]),Form("_muon%s",var[i]),
			      canv_muon_2d[i],h_muon_matrices[i],
			      h_overlay_muon[i][1],muon_max[i],false);  

    //Iterations Test using NuWro Data as input
    iteration_test.Iteration_Test(11,h_bnb_muon[i],h_ext_muon[i],h_overlay_muon[i][0], h_overlay_muon[i][1],h_dirt_muon[i],
      				  h_nuwro_muon[i],h_muon_eff[i],h_muon_denom[i],h_muon_matrices[i],Form("_muon%s",var[i]), Form("Muon %s",var_titles[i]),true);

    //Make the MC Comparison Plots
    mc.model_comparison(h_empirical_muon[i], h_nieves_muon[i], h_susa_muon[i], h_GCF_muon[i],muon_max_theory[i],true,Form("True Muon %s",var_titles[i]),Form("_muon%s",var[i]));

    //Now to extract the cross-section
    cross_section(h_empirical_muon[i], h_nieves_muon[i], h_susa_muon[i], h_GCF_muon[i],
		  h_ext_muon[i], h_dirt_muon[i],h_overlay_muon[i][0], h_overlay_muon[i][1], h_bnb_muon[i],
		  h_muon_matrices[i], muon_num_iterations[i], h_muon_eff[i], h_muon_denom[i],h_muon_nuwro_denom[i],
		  h_muon_systematic[i], muon_xsec_max[i],Form("True Muon %s",var_titles[i]),Form("_muon%s",var[i]),flip_legend_muon[i],true);
   

    
    //Leading
    /////////////

    //Make the Efficiency Histogram
    h_leading_eff[i] = make_efficiency_plot(h_leading_num[i],h_leading_denom[i],Form("True Leading Proton %s",var_titles[i]),Form("_leading%s",var[i]));

    //Plot the matrix and make the smearing matrix
    plot_matrices(h_leading_matrices[i],Form("True Leading Proton %s",var_titles[i]),Form("_leading%s",var[i]));
   
    //Closure test
    closure_test.Closure_Test(canv_leading_eff[i],h_leading_eff[i],h_leading_num[i],h_leading_denom[i],Form("Leading Proton %s",var_titles[i]),Form("_leading%s",var[i]),
			      canv_leading_2d[i],h_leading_matrices[i],
			      h_overlay_leading[i][1],leading_max[i],false);  

    //Iterations Test using NuWro Data as input
    iteration_test.Iteration_Test(11,h_bnb_leading[i],h_ext_leading[i],h_overlay_leading[i][0],h_overlay_leading[i][1],h_dirt_leading[i],
    				  h_nuwro_leading[i],h_leading_eff[i],h_leading_denom[i],h_leading_matrices[i], Form("_leading%s",var[i]), Form("True Leading Proton %s",var_titles[i]),true);

    //Make the MC Comparison Plots
    mc.model_comparison(h_empirical_leading[i], h_nieves_leading[i], h_susa_leading[i], h_GCF_leading[i],leading_max_theory[i],true,Form("True Leading Proton %s",var_titles[i]),Form("_leading%s",var[i]));

    //Now to extract the cross-section
    cross_section(h_empirical_leading[i], h_nieves_leading[i], h_susa_leading[i], h_GCF_leading[i],
		  h_ext_leading[i], h_dirt_leading[i],h_overlay_leading[i][0], h_overlay_leading[i][1], h_bnb_leading[i],
		  h_leading_matrices[i], leading_num_iterations[i], h_leading_eff[i], h_leading_denom[i],h_leading_nuwro_denom[i],
		  h_leading_systematic[i],leading_xsec_max[i], Form("True Leading Proton %s",var_titles[i]),Form("_leading%s",var[i]),flip_legend_leading[i],false);
    

    //Recoil
    /////////////
 
    //Make the Efficiency Histogram
    h_recoil_eff[i] = make_efficiency_plot(h_recoil_num[i],h_recoil_denom[i],Form("True Recoil Proton %s",var_titles[i]),Form("_recoil%s",var[i]));

    //Plot the matrix and make the smearing matrix
    plot_matrices(h_recoil_matrices[i],Form("True Recoil Proton %s",var_titles[i]),Form("_recoil%s",var[i]));
   
    //Closure test
    closure_test.Closure_Test(canv_recoil_eff[i],h_recoil_eff[i],h_recoil_num[i],h_recoil_denom[i],Form("Recoil Proton %s",var_titles[i]),Form("_recoil%s",var[i]),
			      canv_recoil_2d[i],h_recoil_matrices[i],
			      h_overlay_recoil[i][1],recoil_max[i],false);  

    //Iterations Test using NuWro Data as input
    iteration_test.Iteration_Test(11,h_bnb_recoil[i],h_ext_recoil[i],h_overlay_recoil[i][0],h_overlay_recoil[i][1],h_dirt_recoil[i],
                                h_nuwro_recoil[i],h_recoil_eff[i],h_recoil_denom[i],h_recoil_matrices[i], Form("_recoil%s",var[i]), Form("True Recoil Proton %s",var_titles[i]),true);

    //Make the MC Comparison Plots
    mc.model_comparison(h_empirical_recoil[i], h_nieves_recoil[i], h_susa_recoil[i], h_GCF_recoil[i],recoil_max_theory[i],true,Form("True Recoil Proton %s",var_titles[i]),Form("_recoil%s",var[i]));

    //Now to extract the cross-section
    cross_section(h_empirical_recoil[i], h_nieves_recoil[i], h_susa_recoil[i], h_GCF_recoil[i],
    h_ext_recoil[i], h_dirt_recoil[i],h_overlay_recoil[i][0], h_overlay_recoil[i][1], h_bnb_recoil[i],
		  h_recoil_matrices[i], recoil_num_iterations[i], h_recoil_eff[i], h_recoil_denom[i],h_recoil_nuwro_denom[i],
		  h_recoil_systematic[i],recoil_xsec_max[i], Form("True Recoil Proton %s",var_titles[i]), Form("_recoil%s",var[i]),flip_legend_recoil[i],false);
		  
   }
  
  
  ////////////////////////////
  //Now for the Physics Plots
  ////////////////////////////////
  
  for(int i=0; i < num_other_var; i++){
  
    //Make the Efficiency Histogram
    h_other_eff[i] = make_efficiency_plot(h_other_num[i],h_other_denom[i],Form("True %s",other_var_titles[i]), Form("%s",other_var[i]));
  
    //Plot the matrix and make the smearing matrix
    plot_matrices(h_other_matrices[i],Form("True %s",other_var_titles[i]), Form("%s",other_var[i]));
  
    //Closure test
    closure_test.Closure_Test(canv_other_eff[i],h_other_eff[i],h_other_num[i],h_other_denom[i],Form(" %s",other_var_titles[i]), Form("%s",other_var[i]),
			      canv_other_2d[i],h_other_matrices[i],
			      h_overlay_other[i][1],other_max[i],false);  

    //Iterations Test using NuWro Data as input
    iteration_test.Iteration_Test(11,h_bnb_other[i],h_ext_other[i],h_overlay_other[i][0],h_overlay_other[i][1],h_dirt_other[i],
    				  h_nuwro_other[i],h_other_eff[i],h_other_denom[i],h_other_matrices[i],  Form("%s",other_var[i]), Form("True %s",other_var_titles[i]),true);

    if(i == 5 || i == 6){
      mc.make_plot_weird(h_empirical_other[i], h_empirical_other_true[i],
			 h_nieves_other[i], h_nieves_other_true[i],
			 h_susa_other[i], h_susa_other_true[i],
			 h_GCF_other[i],h_GCF_other_true[i],
			 1.0,Form("%s",other_var_titles[i]), Form("%s",other_var[i]));
    }
    
    //Make the MC Comparison Plots
    mc.model_comparison(h_empirical_other[i], h_nieves_other[i], h_susa_other[i], h_GCF_other[i],other_max_theory[i],true,Form("True %s",other_var_titles[i]), Form("%s",other_var[i]));

    
    //Now to extract the cross-section
    cross_section(h_empirical_other[i], h_nieves_other[i], h_susa_other[i], h_GCF_other[i],
    		  h_ext_other[i], h_dirt_other[i],h_overlay_other[i][0], h_overlay_other[i][1], h_bnb_other[i],
    		  h_other_matrices[i], other_num_iterations[i], h_other_eff[i], h_other_denom[i], h_other_nuwro_denom[i],
		  h_other_systematic[i],other_xsec_max[i], Form("True %s",other_var_titles[i]), Form("%s",other_var[i]),flip_legend_other[i],false);

		  
    /*    if(i == 1){
      iteration_test.BNB_Testing(2,3,h_bnb_other[i],h_other_eff[i],h_other_matrices[i],"",other_var[i]);
    } else if(i == 2){
      iteration_test.BNB_Testing(3,7,h_bnb_other[i],h_other_eff[i],h_other_matrices[i],"",other_var[i]);
    } else if(i == 3){
      iteration_test.BNB_Testing(2,3,h_bnb_other[i],h_other_eff[i],h_other_matrices[i],"",other_var[i]);
    } else if(i == 4){
      iteration_test.BNB_Testing(2,3,h_bnb_other[i],h_other_eff[i],h_other_matrices[i],"",other_var[i]);
      }*/
   
  }

  myfile << muon_num_iterations[0] << ","<< muon_num_iterations[1] << "," << muon_num_iterations[2] <<"\n";
  myfile << leading_num_iterations[0] << ","<< leading_num_iterations[1] << ","<< leading_num_iterations[2] <<"\n";
  myfile << recoil_num_iterations[0] << ","<< recoil_num_iterations[1] << ","<< recoil_num_iterations[2] <<"\n";
  myfile << other_num_iterations[0] << "," << other_num_iterations[1] << ","<< other_num_iterations[2] << "," << other_num_iterations[3] << ","<< other_num_iterations[4] << "," << other_num_iterations[5] << ","<< other_num_iterations[6] <<"\n";
  myfile.close(); //close my csv file
    
  
}//end of program

