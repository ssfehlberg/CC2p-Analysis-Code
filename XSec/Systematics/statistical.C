#define statistical_cxx
#include "tools/xsec.h"
#include "tools/histograms.h"
using namespace Histograms;

class statistical{

 public:
  virtual void main();

 private:
  virtual void Grab_Histograms();


  //Xsec Class definition
  xsec Xsec;

  //Extra Histograms that we need
  //////////////////////////////
  TH1D* h_muon_xsec[num_var];
  TH1D* h_muon_stat_error[num_var];
  TH1D* h_muon_stat_fractional_error[num_var];
  TH2D* h_muon_stat_matrix[num_var];
  
  TH1D* h_leading_xsec[num_var];
  TH1D* h_leading_stat_error[num_var];
  TH1D* h_leading_stat_fractional_error[num_var];
  TH2D* h_leading_stat_matrix[num_var];
    
  TH1D* h_recoil_xsec[num_var];
  TH1D* h_recoil_stat_error[num_var];
  TH1D* h_recoil_stat_fractional_error[num_var];
  TH2D* h_recoil_stat_matrix[num_var];

  TH1D* h_other_xsec[num_other_var];
  TH1D* h_other_stat_error[num_other_var];
  TH1D* h_other_stat_fractional_error[num_other_var];
  TH2D* h_other_stat_matrix[num_other_var];
  
};//end of class

void statistical::main(){

  Grab_Histograms();
  TFile *tfile_mine = new TFile("root_files/Statistical/systematics.root","RECREATE"); //output root file 

  for(int j=0; j < num_var; j++){
    
    //Muon
    /////////////////////
    h_muon_xsec[j] =  Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_overlay_muon[j][0],h_overlay_muon[j][1],h_bnb_muon[j],
					 h_muon_matrices[j],muon_iter[j],h_muon_num[j], h_muon_denom[j],
					 muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_CV",var0[j]),Form("Statistical"),false,true); //CV Montecarlo values
    h_muon_stat_error[j] = (TH1D*)h_muon_xsec[j]->Clone(Form("hist_errors_muon%s",var0[j]));
    h_muon_stat_fractional_error[j] = (TH1D*)h_overlay_muon[j][0]->Clone();//Form("hist_fractional_errors_muon%s",var0[j]));
    h_muon_stat_matrix[j] = (TH2D*)h_muon_matrices[j]->Clone(Form("h_2D_covariancemuon_%s",var0[j]));
    
    for(int x = 1; x < h_muon_xsec[j]->GetNbinsX()+1;  x++){
      double bin_content = h_muon_xsec[j]->GetBinContent(x);
      double bin_stat_error = h_muon_xsec[j]->GetBinError(x);

      std::cout<<"Value of bin Content: "<<bin_content<<std::endl;
      std::cout<<"Value of stat error: "<<bin_stat_error<<std::endl;
      std::cout<<"Value of fractional error: "<<bin_stat_error/bin_content<<std::endl;
      
      h_muon_stat_error[j]->SetBinContent(x,bin_stat_error);
      h_muon_stat_fractional_error[j]->SetBinContent(x,(bin_stat_error/bin_content));
      
      for(int y=1; y < h_muon_xsec[j]->GetNbinsX()+1; y++){ //symmetric matrix
	if(x == y){
	  h_muon_stat_matrix[j]->SetBinContent(x,y,std::pow(bin_stat_error,2));
	} else {
	  h_muon_stat_matrix[j]->SetBinContent(x,y,0.0);
	} //end of else
      }	//end loop over y bins
    }//end loop over x bins

    std::cout<<"Value of bin 1 for fractional error: "<<h_muon_stat_fractional_error[j]->GetBinContent(1)<<std::endl;
    h_muon_stat_error[j]->Write();
    //h_muon_stat_matrix->SetTitle(Form("Muon %s",var_titles[j]));
    //h_muon_stat_matrix->GetXaxis()->SetTitle("Reco. Bin a");
    //h_muon_stat_matrix->GetYaxis()->SetTitle("Reco. Bin b");
    h_muon_stat_matrix[j]->Write();
    h_muon_stat_fractional_error[j]->Write(Form("hist_fractional_errors_muon%s",var0[j]));

    //Leading
    /////////////////////
    h_leading_xsec[j] =  Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_overlay_leading[j][0],h_overlay_leading[j][1],h_bnb_leading[j],
					 h_leading_matrices[j],leading_iter[j],h_leading_num[j], h_leading_denom[j],
					 leading_xsec_max[j],Form("True Leading Proton %s",var_titles[j]),Form("_leading%s_CV",var[j]),Form("Statistical"),false,true); //CV Montecarlo values
    h_leading_stat_error[j] = (TH1D*)h_leading_xsec[j]->Clone(Form("hist_errors_leading%s",var0[j]));
    h_leading_stat_fractional_error[j] = (TH1D*)h_overlay_leading[j][0]->Clone();
    h_leading_stat_matrix[j] = (TH2D*)h_leading_matrices[j]->Clone(Form("h_2D_covarianceleading_%s",var0[j]));
    
    for(int x = 1; x < h_leading_xsec[j]->GetNbinsX()+1;  x++){
      double bin_content = h_leading_xsec[j]->GetBinContent(x);
      double bin_stat_error = h_leading_xsec[j]->GetBinError(x);
      h_leading_stat_error[j]->SetBinContent(x,bin_stat_error);
      h_leading_stat_fractional_error[j]->SetBinContent(x,(bin_stat_error/bin_content));

      for(int y=1; y < h_leading_xsec[j]->GetNbinsX()+1; y++){ //symmetric matrix
	if(x == y){
	  h_leading_stat_matrix[j]->SetBinContent(x,y,std::pow(bin_stat_error,2));
	} else {
	  h_leading_stat_matrix[j]->SetBinContent(x,y,0.0);
	} //end of else
      }	//end loop over y bins
    }//end loop over x bins

    h_leading_stat_error[j]->Write();
    h_leading_stat_fractional_error[j]->Write(Form("hist_fractional_errors_leading%s",var0[j]));
    h_leading_stat_matrix[j]->Write();

    //Recoil
    /////////////////////
    h_recoil_xsec[j] =  Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_overlay_recoil[j][0],h_overlay_recoil[j][1],h_bnb_recoil[j],
					 h_recoil_matrices[j],recoil_iter[j],h_recoil_num[j], h_recoil_denom[j],
					 recoil_xsec_max[j],Form("True Recoil Proton %s",var_titles[j]),Form("_recoil%s_CV",var[j]),Form("Statistical"),false,true); //CV Montecarlo values
    h_recoil_stat_error[j] = (TH1D*)h_recoil_xsec[j]->Clone(Form("hist_errors_recoil%s",var0[j]));
    h_recoil_stat_fractional_error[j] = (TH1D*)h_overlay_recoil[j][0]->Clone();
    h_recoil_stat_matrix[j] = (TH2D*)h_recoil_matrices[j]->Clone(Form("h_2D_covariancerecoil_%s",var0[j]));
    
    for(int x = 1; x < h_recoil_xsec[j]->GetNbinsX()+1;  x++){
      double bin_content = h_recoil_xsec[j]->GetBinContent(x);
      double bin_stat_error = h_recoil_xsec[j]->GetBinError(x);
      h_recoil_stat_error[j]->SetBinContent(x,bin_stat_error);
      h_recoil_stat_fractional_error[j]->SetBinContent(x,(bin_stat_error/bin_content));
      
      for(int y=1; y < h_recoil_xsec[j]->GetNbinsX()+1; y++){ //symmetric matrix
	if(x == y){
	  h_recoil_stat_matrix[j]->SetBinContent(x,y,std::pow(bin_stat_error,2));
	} else {
	  h_recoil_stat_matrix[j]->SetBinContent(x,y,0.0);
	} //end of else
      }	//end loop over y bins
    }//end loop over x bins

    h_recoil_stat_error[j]->Write();
    h_recoil_stat_fractional_error[j]->Write(Form("hist_fractional_errors_recoil%s",var0[j]));
    h_recoil_stat_matrix[j]->Write();

  } //end loop over num var


  for(int j=0; j < num_other_var; j++){

    h_other_xsec[j] =  Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_overlay_other[j][0],h_overlay_other[j][1],h_bnb_other[j],
					 h_other_matrices[j],other_iter[j],h_other_num[j], h_other_denom[j],
					 other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("_%s_CV",other_var[j]),Form("Statistical"),false,true); //CV Montecarlo values
    h_other_stat_error[j] = (TH1D*)h_other_xsec[j]->Clone(Form("hist_errors%s",other_var[j]));
    h_other_stat_fractional_error[j] = (TH1D*)h_overlay_other[j][0]->Clone();
    h_other_stat_matrix[j] = (TH2D*)h_other_matrices[j]->Clone(Form("h_2D_covariance%s",other_var[j]));
    
    for(int x = 1; x < h_other_xsec[j]->GetNbinsX()+1;  x++){
      double bin_content = h_other_xsec[j]->GetBinContent(x);
      double bin_stat_error = h_other_xsec[j]->GetBinError(x);
      h_other_stat_error[j]->SetBinContent(x,bin_stat_error);
      h_other_stat_fractional_error[j]->SetBinContent(x,(bin_stat_error/bin_content));
      
      for(int y=1; y < h_other_xsec[j]->GetNbinsX()+1; y++){ //symmetric matrix
	if(x == y){
	  h_other_stat_matrix[j]->SetBinContent(x,y,std::pow(bin_stat_error,2));
	} else {
	  h_other_stat_matrix[j]->SetBinContent(x,y,0.0);
	} //end of else
      }	//end loop over y bins
    }//end loop over x bins

    h_other_stat_error[j]->Write();
    h_other_stat_fractional_error[j]->Write(Form("hist_fractional_errors%s",other_var[j]));
    h_other_stat_matrix[j]->Write();
    

  } //end of loop over other variables
  
  tfile_mine->Close();

}//end of main program

void statistical::Grab_Histograms(){

  //BNB,EXT,Dirt, and Overlay Products
  ////////////////////////////////
  TFile* f_bnb = new TFile(Form("../../root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root"));
  TFile* f_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* f_dirt = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root");
  TFile* f_overlay = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root");
  TFile* f_matrices = new TFile("../../root_files/pelee/Run_all/histograms_mc_eff.root");
  TFile* f_eff = new TFile("../../root_files/pelee/Run_all/histograms_mc_eff.root");
 
  for(int i=0; i < num_var; i++){

    //BNB, EXT, Dirt Data Products
    h_bnb_muon[i] = (TH1D*)f_bnb->Get(Form("h_muon%s_bnb",var0[i])); //bnb
    h_ext_muon[i] = (TH1D*)f_ext->Get(Form("h_muon%s_ext",var0[i])); //ext
    h_dirt_muon[i] = (TH1D*)f_dirt->Get(Form("h_muon%s_dirt_wgt",var0[i])); //dirt

    h_bnb_leading[i] = (TH1D*)f_bnb->Get(Form("h_leading%s_bnb",var0[i])); //bnb
    h_ext_leading[i] = (TH1D*)f_ext->Get(Form("h_leading%s_ext",var0[i])); //ext
    h_dirt_leading[i] = (TH1D*)f_dirt->Get(Form("h_leading%s_dirt_wgt",var0[i])); //dirt  

    h_bnb_recoil[i] = (TH1D*)f_bnb->Get(Form("h_recoil%s_bnb",var0[i])); //bnb
    h_ext_recoil[i] = (TH1D*)f_ext->Get(Form("h_recoil%s_ext",var0[i])); //ext
    h_dirt_recoil[i] = (TH1D*)f_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[i])); //dirt   

    //Overlay Products
    for(int j=0; j < num_channels; j++){
      h_overlay_muon[i][j] = (TH1D*)f_overlay->Get(Form("h_muon%s%s",var0[i],channel[j])); //overlay products
      h_overlay_leading[i][j] = (TH1D*)f_overlay->Get(Form("h_leading%s%s",var0[i],channel[j])); //overlay products
      h_overlay_recoil[i][j] = (TH1D*)f_overlay->Get(Form("h_recoil%s%s",var0[i],channel[j])); //overlay products 
    }

    //Matrices, Numerator, and Denominator
    h_muon_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_muon_all%s",(var[i]))); //smearing matrix
    h_muon_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_muon_all%s",var[i])); //numerator of efficiency
    h_muon_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_muon_all%s",var[i])); //denom of efficiency

    h_leading_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_lead_proton%s",(var[i]))); //smearing matrix
    h_leading_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_lead_proton%s",var[i])); //numerator of efficiency
    h_leading_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_lead_proton%s",var[i])); //denom of efficiency

    h_recoil_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_recoil_proton%s",(var[i]))); //smearing matrix
    h_recoil_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_recoil_proton%s",var[i])); //numerator of efficiency
    h_recoil_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_recoil_proton%s",var[i])); //denom of efficiency  

  } //end of loop over num_var

  for(int i=0; i < num_other_var; i++){

    //BNB, EXT, Dirt Data Products
    h_bnb_other[i] = (TH1D*)f_bnb->Get(Form("h%s_bnb",other_var[i])); //bnb
    h_ext_other[i] = (TH1D*)f_ext->Get(Form("h%s_ext",other_var[i])); //ext
    h_dirt_other[i] = (TH1D*)f_dirt->Get(Form("h%s_dirt_wgt",other_var[i])); //dirt

    //Overlay Products
    for(int j=0; j < num_channels; j++){
      h_overlay_other[i][j] = (TH1D*)f_overlay->Get(Form("h%s%s",other_var[i],channel[j])); //overlay products
    }
    //Matrices, Numerator, and Denominator
    h_other_matrices[i] = (TH2D*)f_matrices->Get(Form("h_other_matrices%s",other_var[i])); //smearing matrix
    h_other_num[i] = (TH1D*)f_eff->Get(Form("h_other_eff_num%s",other_var[i])); //numerator of efficiency
    h_other_denom[i] = (TH1D*)f_eff->Get(Form("h_other_eff_denom%s",other_var[i])); //denom of efficiency
    
  } //end of other loop
  
} //end of grab_histograms







