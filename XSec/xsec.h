#ifndef xsec_h
#define xsec_h

//C++ includes
#include <iostream>
#include <string>
#include <fstream>

//My includes
#include "closure_test.h"
#include "iterations.h"
#include "mc_model_comparison.h"
//#include "systematics.h"
#include "constants.h"
using namespace Constants;

class xsec{
 public:
  virtual void main();
  virtual void Grab_Histograms();
  virtual TH1D* make_efficiency_plot(TH1D* h_num, TH1D* h_denom,const char* title, const char* name);
  virtual void plot_matrices(TH2D* h_matrix,const char* title, const char* name);
  virtual TH1D* make_bnb_plot(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input, const char* name);
  virtual std::vector<TH1D*> mc_plots(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_nuisance, TH1D* h_GCF, bool print_contents = true);
  virtual void Fix_Systematic(TH1D* h,TH1D* h_MC_CV,TH1D* h_systematic, bool print_contents);
  virtual void Get_Systematic(TH1D* h,TH1D* h_stat,TH2D* h_covariance, bool print_contents);
  virtual double Chi2(TH1D* h_bnb,TH1D* h_sys,TH1D* hist);
  virtual double Chi2_Cov(TH1D* h_data, TH1D* h_mc, TH2D* h_cov);
  virtual double CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi);
  virtual void Ratio_Plot(TH1D* h_xsec, TH1D* h_uboone, TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_nuwro, TH1D* h_systematic,int iter, const char* title, const char* name);
  virtual void cross_section(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_nuisance, TH1D* h_GCF,
			     TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			     TH2D* h_smearing, int iter, TH1D* h_eff_input, TH1D* h_denom_input,TH1D* h_denom_nuwro_input,
			     TH1D* h_systematic,TH2D* h_covariance,double maximum, const char* title, const char* name,bool flip_legend,bool print_contents = true);


private:

  const char* unit_title; //units for the cross-section plots 
  
  TH1D* h_muon_nuwro_denom[num_var];
  TH1D* h_leading_nuwro_denom[num_var];
  TH1D* h_recoil_nuwro_denom[num_var];
  TH1D* h_other_nuwro_denom[num_other_var];
  
};

#endif
#ifdef xsec_cxx

void xsec::Grab_Histograms(){

  //BNB,EXT,Dirt, and Overlay Products
  ////////////////////////////////
  TFile* f_bnb = new TFile(Form("../root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root"));
  TFile* f_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* f_dirt = new TFile("../root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root");
  TFile* f_overlay = new TFile("../root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root");
  TFile* f_matrices = new TFile("../root_files/pelee/Run_all/histograms_mc_eff.root");
  TFile* f_eff = new TFile("../root_files/pelee/Run_all/histograms_mc_eff.root");

  TFile* f_nuwro = new TFile("../root_files/nuwro/Run_all/histograms_nuwro_xsec_overlay_wgt.root");
  TFile* f_empirical = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/MEC/hists_empirical_lwellyn_fsi.root");
  TFile* f_nieves = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/MEC/hists_nieves_fsi.root");
  TFile* f_susa = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/MEC/hists_susav2_fsi.root");
  TFile* f_nuisance = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/MEC/hists_Nusiance_fsi_xsec.root");
  TFile* f_GCF = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/GCF/hists_GCF_CCQE_fsi.root");

  TFile* f_covar = new TFile("Systematics/root_files/Total_Error/total_covariance_matrices.root"); //total covariance matrices
  
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

    //MC Models
    h_nuwro_muon[i] = (TH1D*)f_nuwro->Get(Form("h_muon%s_total",var0[i])); //nuwro
    h_empirical_muon[i] = (TH1D*)f_empirical->Get(Form("h_muon%s_lead_cut",var0[i])); //empirical
    h_nieves_muon[i] = (TH1D*)f_nieves->Get(Form("h_muon%s_lead_cut",var0[i])); //nieves
    h_susa_muon[i] = (TH1D*)f_susa->Get(Form("h_muon%s_lead_cut",var0[i])); //susa
    h_nuisance_muon[i] = (TH1D*)f_nuisance->Get(Form("h_muon%s_lead_cut",var0[i]));
    h_GCF_muon[i] = (TH1D*)f_GCF->Get(Form("h_muon%s_lead_cut",var0[i])); //gcf

    h_nuwro_leading[i] = (TH1D*)f_nuwro->Get(Form("h_leading%s_total",var0[i])); //nuwro
    h_empirical_leading[i] = (TH1D*)f_empirical->Get(Form("h_leading%s_lead_cut",var0[i])); //empirical
    h_nieves_leading[i] = (TH1D*)f_nieves->Get(Form("h_leading%s_lead_cut",var0[i])); //nieves
    h_susa_leading[i] = (TH1D*)f_susa->Get(Form("h_leading%s_lead_cut",var0[i])); //susa
    h_nuisance_leading[i] = (TH1D*)f_nuisance->Get(Form("h_leading%s_lead_cut",var0[i]));
    h_GCF_leading[i] = (TH1D*)f_GCF->Get(Form("h_leading%s_lead_cut",var0[i])); //gcf

    h_nuwro_recoil[i] = (TH1D*)f_nuwro->Get(Form("h_recoil%s_total",var0[i])); //nuwro
    h_empirical_recoil[i] = (TH1D*)f_empirical->Get(Form("h_recoil%s_lead_cut",var0[i])); //empirical
    h_nieves_recoil[i] = (TH1D*)f_nieves->Get(Form("h_recoil%s_lead_cut",var0[i])); //nieves
    h_susa_recoil[i] = (TH1D*)f_susa->Get(Form("h_recoil%s_lead_cut",var0[i])); //susa
    h_nuisance_recoil[i] = (TH1D*)f_nuisance->Get(Form("h_recoil%s_lead_cut",var0[i]));
    h_GCF_recoil[i] = (TH1D*)f_GCF->Get(Form("h_recoil%s_lead_cut",var0[i])); //gcf

    //Total Covariance Matrices
    h_covar_muon[i] = (TH2D*)f_covar->Get(Form("h_2D_covariancemuon_%s",var0[i]));
    h_covar_leading[i] = (TH2D*)f_covar->Get(Form("h_2D_covarianceleading_%s",var0[i]));
    h_covar_recoil[i] = (TH2D*)f_covar->Get(Form("h_2D_covariancerecoil_%s",var0[i]));
    
  }

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

    //MC Models
    h_nuwro_other[i] = (TH1D*)f_nuwro->Get(Form("h%s_total",other_var[i])); //nuwro
    h_empirical_other[i] = (TH1D*)f_empirical->Get(Form("h%s_lead_cut",other_var[i])); //empirical
    h_nieves_other[i] = (TH1D*)f_nieves->Get(Form("h%s_lead_cut",other_var[i])); //nieves
    h_susa_other[i] = (TH1D*)f_susa->Get(Form("h%s_lead_cut",other_var[i])); //susa
    h_nuisance_other[i] = (TH1D*)f_nuisance->Get(Form("h%s_lead_cut",other_var[i]));
    h_GCF_other[i] = (TH1D*)f_GCF->Get(Form("h%s_lead_cut",other_var[i])); //gcf

    //total covariance matrices
    h_covar_other[i] = (TH2D*)f_covar->Get(Form("h_2D_covariance%s",other_var[i]));
    
    //MC Models: Truth level for neutrino energy and pn
    if(i == 5 || i == 6){

      h_empirical_other_true[i] = (TH1D*)f_empirical->Get(Form("h%s_lead_cut",other_var[i])); //empirical
      h_nieves_other_true[i] = (TH1D*)f_nieves->Get(Form("h%s_lead_cut",other_var[i])); //nieves
      h_susa_other_true[i] = (TH1D*)f_susa->Get(Form("h%s_lead_cut",other_var[i])); //susa
      h_nuisance_other[i] = (TH1D*)f_nuisance->Get(Form("h%s_lead_cut",other_var[i]));
      h_GCF_other_true[i] = (TH1D*)f_GCF->Get(Form("h%s_lead_cut",other_var[i])); //gcf  

      h_empirical_other[i] = (TH1D*)f_empirical->Get(Form("h%s_true_lead_cut",other_var[i])); //empirical
      h_nieves_other[i] = (TH1D*)f_nieves->Get(Form("h%s_true_lead_cut",other_var[i])); //nieves
      h_susa_other[i] = (TH1D*)f_susa->Get(Form("h%s_true_lead_cut",other_var[i])); //susa
      h_nuisance_other[i] = (TH1D*)f_nuisance->Get(Form("h%s_lead_cut",other_var[i]));
      h_GCF_other[i] = (TH1D*)f_GCF->Get(Form("h%s_true_lead_cut",other_var[i])); //gcf  
    }
  }

  //Systematic Uncertainty
  ////////////////////////
  TFile* f_systematic = new TFile("Systematics/root_files/total_error.root"); //this total error includes flux, reint, GENIE total, detvar, dirt, and iteration uncertainty

   for(int i = 0; i < num_var; i++){
     h_muon_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_muon%s",var0[i]));
     h_leading_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_leading%s",var0[i]));
     h_recoil_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_recoil%s",var0[i]));
   }
   
   for(int i=0; i < num_other_var; i++){
     h_other_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors%s",other_var[i]));
   }


   //Grabing the NuWro Shit
   //////////////////////////
   TFile* f_nuwro_matrices = new TFile("../root_files/nuwro/Run_all/histograms_mc_eff.root");
   
   for(int i = 0; i < num_var; i++){
     h_muon_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_muon_all%s",var[i])); //true nuwro cc2p0pi
     h_leading_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_lead_proton%s",var[i])); //true nuwro cc2p0pi
     h_recoil_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_recoil_proton%s",var[i])); //true nuwro cc2p0pi
   }
   for(int i=0; i < num_other_var; i++){
     h_other_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_other_eff_denom%s",other_var[i])); //true nuwro cc2p0pi
   }
 
   
} //end of grab histograms

//Makes the efficiency plot from the numerator and denominator
TH1D* xsec::make_efficiency_plot(TH1D* h_num_input, TH1D* h_denom_input,const char* title, const char* name){

  TCanvas* canv_eff = new TCanvas("canv_eff","canv_eff",2000,1500);
  TH1D* h_num = (TH1D*)h_num_input->Clone();
  TH1D*	h_denom = (TH1D*)h_denom_input->Clone();
  h_num->Divide(h_num,h_denom,1.0,1.0, "cp");
  h_num->Draw("1e1p");
  h_num->SetTitle(Form("%s ; %s ; Efficiency",title,title));
  h_num->SetLineColor(kViolet);
  h_num->SetMaximum(0.4);
  h_num->SetMinimum(0);

  //TLatex* t;
  //t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",title));
  //t->DrawLatex(0.3,0.92,Form("%s",pot_num));
  //t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");

  //a[i] = new TLine(xlim_eff[i],0,xlim_eff[i],1);
  //a[i]->Draw("same");
  //a[i]->SetLineColor(kBlack);
  //a[i]->SetLineWidth(4);
  //a1[i] = new TLine(xlim_high_eff[i],0,xlim_high_eff[i],1);
  //a1[i]->Draw("same");
  //a1[i]->SetLineColor(kBlack);
  //a1[i]->SetLineWidth(4);

  canv_eff->Print(Form("images/Prep/%s_eff.png",name)); 

  return h_num;
  
}

//Plots the raw 2D matrix of reco versus truth
//Then normalizes the columns to 1. This is the same thing as creating the smearing matrix
//Then plots the smearing matrix
//a) h_matrix = raw Th2D
//b) title = histogram titles
//c) name = name to save plot as
///////////////////////////////////////////////////////////
void xsec::plot_matrices(TH2D* h_matrix,const char* title, const char* name){

  //Drawing
  TCanvas* canv_matrix = new TCanvas("canv_matrix","canv_matrix",2000,1500);
  h_matrix->Draw("colz text");
  h_matrix->SetTitle(Form("%s ; True %s; Reco. %s",title,title,title));

  /* TText* t;
  t->DrawText(0.3,0.92,Form("%s",pot_num0));
  t->DrawText(0.82,0.92,Form("%s",sample_name0));*/
  canv_matrix->Update();
  canv_matrix->Print(Form("images/Prep/%s_matrix.png",name));
  
  //checking normalization
  TCanvas* canv1_matrix = new TCanvas("canv1_matrix","canv1_matrix",2000,1500);
  TH2D* h_matrix_normalized = (TH2D*)h_matrix->Clone();
  
  int NBinsX = h_matrix_normalized->GetXaxis()->GetNbins();
  int NBinsY = h_matrix_normalized->GetYaxis()->GetNbins();
  
  for(int bin_true = 0; bin_true < NBinsX; bin_true++){
    double NEventsInColumn = 0;
    
    for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
      double bin_content = h_matrix_normalized->GetBinContent(bin_true+1,bin_reco+1);
      NEventsInColumn += bin_content;
    }
    
    for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
      if(NEventsInColumn > 0){
	double FracErr =  h_matrix_normalized->GetBinError(bin_true+1,bin_reco+1) / h_matrix_normalized->GetBinContent(bin_true+1,bin_reco+1);
	double CV = double(  h_matrix_normalized->GetBinContent(bin_true+1,bin_reco+1))/ double(NEventsInColumn);
	h_matrix_normalized->SetBinContent(bin_true+1,bin_reco+1,CV);
	
	double error = CV * TMath::Sqrt( TMath::Power(FracErr,2.) + 1./double(NEventsInColumn) );
	h_matrix_normalized->SetBinError(bin_true+1,bin_reco+1,error) ; 
	
      } else { 
	h_matrix_normalized->SetBinContent(bin_true+1,bin_reco+1,0); 
	h_matrix_normalized->SetBinError(bin_true+1,bin_reco+1,0); 
      } //ends else  
    }
  }
      
  h_matrix_normalized->Draw("colz text");
  h_matrix_normalized->SetTitle(Form("%s; True %s; Reco. %s",title,title,title));

  /* t->DrawText(0.3,0.92,Form("%s",pot_num));
     t->DrawText(0.82,0.92,Form("%s",sample_name));*/
  canv1_matrix->Update();
  canv1_matrix->Print(Form("images/Prep/%s_smearing_matrix.png",name));
      
}



TH1D* xsec::make_bnb_plot(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input, const char* name){

  //Creating the Background: add the backgrounds together: h_ext is what we will subtract
  ////////////////////////////////////////////////////////////////////////////////////////
  TCanvas* c_background = new TCanvas("c6","Total Background",800,600); //you have to define a canvas to draw.

  TH1D* h_ext = (TH1D*)h_ext_input->Clone(); //ext background
  std::cout<<"Num entries in h_ext: "<<h_ext->Integral()<<std::endl;

  TH1D* h_dirt = (TH1D*)h_dirt_input->Clone(); //dirt background
  std::cout<<"Num entries in h_dirt: "<<h_dirt->Integral()<<std::endl;

  h_ext->Add(h_dirt);
  std::cout<<"Num entries in h_ext After Adding Dirt: "<<h_ext->Integral()<<std::endl;

  TH1D* h_overlay_total = (TH1D*)h_overlay_total_input->Clone(); //overlay background
  std::cout<<"Number of Entries in h_overlay: "<<h_overlay_total->Integral()<<std::endl;
  TH1D* h_overlay_cc2p = (TH1D*)h_overlay_cc2p_input->Clone();
  std::cout<<"Number of Entries in h_overlay_cc2p: "<<h_overlay_cc2p->Integral()<<std::endl;

  h_overlay_total->Add(h_overlay_cc2p, -1); // make sure to taake the cc2p out of the total
  std::cout<<"Num entries in h_overlay after Subtraction: "<<h_overlay_total->Integral()<<std::endl;

  h_ext->Add(h_overlay_total);
  std::cout<<"Num entries in h_ext after everything: "<<h_ext->Integral()<<std::endl;
  h_ext->Draw("HIST");
  c_background->Print(Form("images/BNB_Checks/ext_after_subtraction%s.png",name));

  //Draw the subtracted backgrounds as a check
  TCanvas* canv_bnb = new TCanvas("canv_bnb","BNB with Background Subtracted",800,600); //you have to define a canvas to draw.
  TH1D* h_bnb = (TH1D*)h_bnb_input->Clone();
  std::cout<<"Num entries in h_bnb: "<<h_bnb->Integral()<<std::endl;
  h_bnb->Add(h_ext,-1);
  std::cout<<"Num entries in h_bnb After Subtraction: "<<h_bnb->Integral()<<std::endl;
  h_bnb->Draw("HIST");
  canv_bnb->Print(Form("images/BNB_Checks/bnb_after_everything_subtracted%s.png",name));

  return h_bnb;
      
}

std::vector<TH1D*> xsec::mc_plots(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D*  h_nuisance, TH1D* h_GCF, bool print_contents = true){

  //Before we can plot the MEC models on the same plot as the BNB, we have to "normalize" using the procedure outlined
  // in section 7.4 of this document: https://arxiv.org/pdf/2101.11867.pdf
  // We are specifically refering to equations 98: MC Estimate of the XSec and 99: The SD of this calculation
  // (98) dsigma/dx = <sigma> * n/ (N* deltax) 
  //where deltax is bin width, n is # of events in a bin, N is total # of events, and sigma is the flux weighted xsec
  //(99) SD(dsigma/dx) = <sigma>/(deltax * N) * (((N-n)n)/N)^(1/2)

  std::vector<TH1D*> mc_vector;
  
  //Here is were we calculate the dsigma/dx and the SD of dsigma/dx
  ////////////////////////////////////////////////////////////////
  int n_bins = h_empirical->GetNbinsX();
  for(int i=1; i < n_bins+1; i++){
      
    double delta_x = h_empirical->GetBinWidth(i);
    std::cout<<"Width of Bin: "<<delta_x<<std::endl;
    
    //empirical
    double n_empirical = h_empirical->GetBinContent(i);
    double value_empirical = (sigma_empirical * n_empirical)/( N_empirical * delta_x);
    double SD_empirical = (sigma_empirical)/(delta_x*N_empirical) * std::sqrt(((N_empirical - n_empirical)*n_empirical)/(N_empirical));
    h_empirical->SetBinContent(i,value_empirical);
    h_empirical->SetBinError(i,SD_empirical);
    
    //nieves
    double n_nieves = h_nieves->GetBinContent(i);
    double value_nieves = (sigma_nieves * n_nieves)/(N_nieves * delta_x);
    double SD_nieves = (sigma_nieves)/(delta_x*N_nieves) * std::sqrt(((N_nieves - n_nieves)*n_nieves)/(N_nieves));
    h_nieves->SetBinContent(i,value_nieves);
    h_nieves->SetBinError(i,SD_nieves);
    
    //susa
    double n_susa = h_susa->GetBinContent(i);
    double value_susa = (sigma_susa * n_susa)/(N_susa * delta_x);
    double SD_susa = (sigma_susa)/(delta_x*N_susa) * std::sqrt(((N_susa - n_susa)*n_susa)/(N_susa));
    h_susa->SetBinContent(i,value_susa);
    h_susa->SetBinError(i,SD_susa);

    //nuisance
    double n_nuisance = h_nuisance->GetBinContent(i);
    double value_nuisance = (sigma_nuisance * n_nuisance)/(N_nuisance * delta_x);
    double SD_nuisance = (sigma_nuisance)/(delta_x*N_nuisance) * std::sqrt(((N_nuisance - n_nuisance)*n_nuisance)/(N_nuisance));
    h_nuisance->SetBinContent(i,value_nuisance);
    h_nuisance->SetBinError(i,SD_nuisance);

    //Dealing with the GCF is a bit trickier
    //Each event in the GCF is assigned a weight which modifies the CCQE differential cross section
    double n_GCF = h_GCF->GetBinContent(i); 
    double value_GCF = (n_GCF)/(N_GCF*delta_x);
    double SD_GCF = (1)/(delta_x*N_GCF) * std::sqrt(((N_GCF - n_GCF)*n_GCF)/(N_GCF));
    h_GCF->SetBinContent(i,value_GCF);
    h_GCF->SetBinError(i,SD_GCF);

    if(print_contents){

      std::cout<<"Width of Bin: "<<delta_x<<std::endl;

      std::cout<<"n_empirical: "<<n_empirical<<std::endl;
      std::cout<<sigma_empirical*n_empirical<<std::endl;
      std::cout<<N_empirical * delta_x<<std::endl;
      std::cout<<"Value_empirical: "<<value_empirical<<std::endl;

      std::cout<<"n_nieves: "<<n_nieves<<std::endl;
      std::cout<<sigma_nieves*n_nieves<<std::endl;
      std::cout<<N_nieves * delta_x<<std::endl;
      std::cout<<"Value_nieves: "<<value_nieves<<std::endl;

      std::cout<<"n_susa: "<<n_susa<<std::endl;
      std::cout<<sigma_susa*n_susa<<std::endl;
      std::cout<<N_susa * delta_x<<std::endl;
      std::cout<<"Value_susa: "<<value_susa<<std::endl;

      std::cout<<"n_nuisance: "<<n_nuisance<<std::endl;
      std::cout<<sigma_nuisance*n_nuisance<<std::endl;
      std::cout<<N_nuisance * delta_x<<std::endl;
      std::cout<<"Value_nuisance: "<<value_nuisance<<std::endl;
      
      std::cout<<"n_GCF: "<<n_GCF<<std::endl;
      std::cout<<"Value_GCF: "<<value_GCF<<std::endl;

    } 
  }

  mc_vector.push_back(h_empirical);
  mc_vector.push_back(h_nieves);
  mc_vector.push_back(h_susa);
  mc_vector.push_back(h_nuisance);
  mc_vector.push_back(h_GCF);

  return mc_vector;

  
} //end of mc_plots

//Changes the systematics from fraactional uncertainty to cross-section uncertainty
void xsec::Fix_Systematic(TH1D* h,TH1D* h_MC_CV,TH1D* h_systematic, bool print_contents){

  double nbins = h_systematic->GetXaxis()->GetNbins();
  for(int i = 1; i < nbins+1; i++){
    double MC_CV = h_MC_CV->GetBinContent(i);
    double systematic = h_systematic->GetBinContent(i);
    double value = MC_CV * systematic;
    if(print_contents){
      std::cout<<"Value of MC_CV: "<<MC_CV<<std::endl;
      std::cout<<"Value of systematic: "<<systematic<<std::endl;
      std::cout<<"Value of Value: "<<value<<std::endl;
    }
    h->SetBinError(i,value);
  }
}

//Changes the systematics from fraactional uncertainty to cross-section uncertainty
void xsec::Get_Systematic(TH1D* h,TH1D* h_stat, TH2D* h_covariance, bool print_contents){

  double nbins = h_stat->GetXaxis()->GetNbins();
  TH2D* h_covar_clone = (TH2D*)h_covariance->Clone();
  //h_covar_clone->Scale(1E76);
  
  for(int i = 1; i < nbins+1; i++){
    double stat_unc = h_stat->GetBinError(i);
    double covar_element = h_covar_clone->GetBinContent(i,i);
    double value = std::sqrt(covar_element - std::pow(stat_unc,2));
    h->SetBinError(i,value);

    if(print_contents){
      std::cout<<"Value of i: "<<i<<std::endl;
      std::cout<<"Value of statistical uncertainty: "<<stat_unc<<std::endl;
      std::cout<<"Value of Covar_element for i i: "<<covar_element<<std::endl;
      std::cout<<"Value of Value: "<<value<<std::endl;
    }
  }
  delete h_covar_clone;
  
}

//Chi2: Calculates the chi2 using chi2 = sum_{i = 1,# of x bins} (Data_i - MC_i)^2/(Total_Error_i)
//NOTE: This particuolar version of chi2 assumes no correlations between bins in the covariance matrix
// i.e. the covariance matrix is diagonal. In general, this is not true, so we must use the weighted least squares calculation
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double xsec::Chi2(TH1D* h_bnb,TH1D* h_sys,TH1D* hist){

  int nbins = h_bnb->GetNbinsX();
  double chi2 = 0;
  
  for(int i = 1; i < nbins+1; i++){

    double observed_i = h_bnb->GetBinContent(i);
    double expected_i = hist->GetBinContent(i);

    double stat_error = h_bnb->GetBinError(i);
    double sys_error = h_sys->GetBinError(i);
    double bin_error = std::pow(stat_error,2) + std::pow(sys_error,2);

    if(bin_error == 0.0){
      chi2 += 0;
    }
    else{
      chi2 += std::pow((observed_i - expected_i),2)/(bin_error);
      //std::cout<<"Value of obsesrved in bin: "<<observed_i<<std::endl;
      //std::cout<<"Value of expected in bin: "<<expected_i<<std::endl;
      //std::cout<<"Value of stat error: "<<stat_error<<std::endl;
      //std::cout<<"Value of sys error: "<<sys_error<<std::endl;
      //std::cout<<"Value of bin error: "<<bin_error<<std::endl;
      //std::cout<<"Value of chi2 in the loop: "<<chi2<<std::endl;
    } 
  }

  //std::cout<<"value of chi2: "<<chi2<<std::endl;
  
  return chi2;
  
}


//Chi2_Cov: Calculates the chi2 using weighted least squares i.e. chi2 = r^T * W * r
// where r is the residuals and W is the inverse of the covariance matrix.
// When talking about cross sections this equation becomes the follwoing:
// chi2 = sum_i sum_j (Data_i - MC_i)*(Covariance_ij)^(-1)*(Data_j - MC_j)
// where i and j run from 1 to # of bins. Note that the Covariance Matrix is the
// total covariance matrix (i.e. systematic + statistical). 
////////////////////////////////////////////////////////////////////////
double xsec::Chi2_Cov(TH1D* h_data, TH1D* h_mc, TH2D* h_cov){

  TH1D* h_model_clone = (TH1D*)h_mc->Clone();
  TH1D* h_data_clone  = (TH1D*)h_data->Clone();
  TH2D* h_cov_clone   = (TH2D*)h_cov->Clone();
  
  int nbins_x = h_cov_clone->GetNbinsX();
  int nbins_y = h_cov_clone->GetNbinsY();
  double chi2 = 0;

  //Creating the covariance matrix from the TH2D object
  ////////////////////////////////////////////////////
  TMatrixD covariance_matrix;
  covariance_matrix.Clear();
  covariance_matrix.ResizeTo(nbins_x,nbins_y);
  for(int i=0; i < nbins_x; i++){
    for(int j=0; j < nbins_y; j++){
      covariance_matrix[i][j] = h_cov_clone->GetBinContent(i+1,j+1);
    }
  }

  TMatrixD copy_covar_matrix = covariance_matrix; //Copy of the covariance matrix
  TMatrixD inver_covar_matrix = covariance_matrix.Invert(); //inverse of the covariance matrix

  //Checking to make sure we get the identiy matrix
  /*if(_debug){
    TMatrixD Iden;
    Iden.Mult(covariance_matrix,inver_covar_matrix); //should be the identity matrix
    for(int i=0; i < nbins_x+1; i++){
      for(int j=0; j < nbins_y+1; j++){
	std::cout<<Form("Value of Identity Matrix at i=%d and j=%d: %f",i,j,Iden[i][j])<<std::endl;
      }
    }    
    } //end of debug loop*/

  //Now to calculate the chi2
  ///////////////////////////
  for(int i = 0; i < nbins_x; i++){
    for(int j = 0; j < nbins_x; j++){
      double diff_i = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
      double diff_j = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);

      if(i==j){
	double LocalChi = diff_i * inver_covar_matrix[i][j] * diff_j; 
	chi2 += LocalChi;
      }
    } 
  }
  
  delete h_model_clone;
  delete h_data_clone;
  delete h_cov_clone;

  
  if(_debug) std::cout<<Form("value of chi2: %f",chi2)<<std::endl;
  return chi2;
  
}

//Afro's Chi2 calculator
double xsec::CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi){//, int &ndof, double &pval) {

  std::cout<<"[CalcChiSquared]: Beginning Calculation"<<std::endl;
  
  // Clone them so we can scale them 
    TH1D* h_model_clone = (TH1D*)h_model->Clone();
    TH1D* h_data_clone  = (TH1D*)h_data->Clone();
    TH2D* h_cov_clone   = (TH2D*)cov->Clone();
    //h_cov_clone->Scale(1E76);
    
    int NBins = h_cov_clone->GetNbinsX();
    // Getting covariance matrix in TMatrix form
    TMatrixD cov_m;
    cov_m.Clear();
    cov_m.ResizeTo(NBins,NBins);
    // loop over rows
    for (int i = 0; i < NBins; i++) {           
        // loop over columns
        for (int j = 0; j < NBins; j++) {
            cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
        }
    }
    TMatrixD copy_cov_m = cov_m;
    // Inverting the covariance matrix
    TMatrixD inverse_cov_m = cov_m.Invert();
    // Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
    // x = data, mu = model, E^(-1) = inverted covariance matrix 
    chi = 0.;
    
    for (int i = 0; i < NBins; i++) {
        //double XWidth = h_data_clone->GetBinWidth(i+1);
        for (int j = 0; j < NBins; j++) {
            //double YWidth = h_data_clone->GetBinWidth(i+1);
            double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
            double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
	    double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
	    chi += LocalChi;
	    std::cout<<"[CalcChiSquared]:Value of diffi for i="<<i<<": "<<diffi<<std::endl;
            std::cout<<"[CalcChiSquared]:Value of diffj for j="<<j<<": "<<diffj<<std::endl;
	    std::cout<<"[CalcChiSquared]:Value of inverse_cov_m for i="<<i<<" j="<<j<<": "<<inverse_cov_m[i][j]<<std::endl;
	    std::cout<<"[CalcChiSquared]:Value of LocalChi for i="<<i<<" j="<<j<<": "<<LocalChi<<std::endl;
	    std::cout<<"[CalcChiSquared]:Value of SumChi for i="<<i<<" j="<<j<<": "<<chi<<std::endl;
	    
        }
    }
    delete h_model_clone;
    delete h_data_clone;
    delete h_cov_clone;

    std::cout<<"[CalcChiSquared]:Value of Total Chi: "<<chi<<std::endl;
    return chi;
}

void xsec::Ratio_Plot(TH1D* h_xsec, TH1D* h_uboone, TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_nuwro,TH1D* h_systematic,int iter, const char* title, const char* name){

  gStyle->SetPaintTextFormat("4.2f");
  
  TH1D* h_BNB = (TH1D*)h_xsec->Clone();
  TH1D* h_bnb = (TH1D*)h_xsec->Clone();
  TH1D* h_UB = (TH1D*)h_uboone->Clone();
  TH1D* h_Empirical = (TH1D*)h_empirical->Clone();
  TH1D* h_Nieves = (TH1D*)h_nieves->Clone();
  TH1D* h_Susa = (TH1D*)h_susa->Clone();
  TH1D* h_Nuwro = (TH1D*)h_nuwro->Clone();
 
  //Statistical Error
  TH1D* h_stat_up = (TH1D*)h_BNB->Clone();
  int nbins = h_BNB->GetNbinsX();
  for(int i =1; i < nbins+1; i++){
    double error = h_BNB->GetBinError(i);
    double content = h_BNB->GetBinContent(i);
    double frac_error = error/content;
    h_stat_up->SetBinContent(i, frac_error);
  }
  TH1D* h_stat_down = (TH1D*)h_stat_up->Clone();
  h_stat_down->Scale(-1.0);

  //Systematic Error
  TH1D* h_sys_up = (TH1D*)h_systematic->Clone();
  TH1D* h_sys_down = (TH1D*)h_systematic->Clone();
  h_sys_down->Scale(-1.0);
  
  h_BNB->Divide(h_BNB,h_bnb,1,1,"b");
  h_UB->Divide(h_bnb,h_UB,1,1,"b");
  h_Empirical->Divide(h_bnb,h_Empirical,1,1,"b");
  h_Nieves->Divide(h_bnb,h_Nieves,1,1,"b");
  h_Susa->Divide(h_bnb,h_Susa,1,1,"b");
  h_Nuwro->Divide(h_bnb,h_Nuwro,1,1,"b");

  for(int i=1; i < nbins+1; i++){
    h_BNB->SetBinContent(i, (h_BNB->GetBinContent(i) - 1.0));
    h_UB->SetBinContent(i, (h_UB->GetBinContent(i) - 1.0));
    h_Empirical->SetBinContent(i, (h_Empirical->GetBinContent(i) - 1.0));
    h_Nieves->SetBinContent(i, (h_Nieves->GetBinContent(i) - 1.0));
    h_Susa->SetBinContent(i, (h_Susa->GetBinContent(i) - 1.0));
    h_Nuwro->SetBinContent(i, (h_Nuwro->GetBinContent(i) - 1.0));  
  }

  TCanvas* c_ratio = new TCanvas("c_ratio","c_ratio",2000,1500);
  //c_ratio->SetGridx();
  h_BNB->Draw("hist");
  h_BNB->SetLineColor(kBlack);
  h_BNB->SetLineWidth(4);
  h_BNB->SetXTitle(Form("%s",title));
  h_BNB->SetMaximum(3.5);
  h_BNB->SetMinimum(-1);
  h_BNB->GetXaxis()->SetTitleSize(50); //35
  h_BNB->GetXaxis()->SetTitleFont(43);
  h_BNB->GetXaxis()->SetTitleOffset(1.3);
  h_BNB->GetXaxis()->SetLabelFont(43);
  h_BNB->GetXaxis()->SetLabelSize(50);
  h_BNB->SetYTitle("(Measured-Theoretical)/Theoretical");
  h_BNB->GetYaxis()->SetTitleSize(50);
  h_BNB->GetYaxis()->SetTitleFont(43); //4 = hevelatica normal 3 = precision
  h_BNB->GetYaxis()->SetTitleOffset(1.3);
  h_BNB->GetYaxis()->SetLabelFont(43);
  h_BNB->GetYaxis()->SetLabelSize(50);

  h_sys_up->Draw("hist same");
  h_sys_up->SetFillColorAlpha(kGray,0.4);
  h_sys_up->SetLineWidth(0);
  h_sys_down->Draw("hist same");
  h_sys_down->SetFillColorAlpha(kGray,0.4);
  h_sys_down->SetLineWidth(0);

  /*h_stat_up->Draw("hist same");
  h_stat_up->SetFillColorAlpha(kGray+2,0.4);
  h_stat_up->SetLineWidth(0);
  h_stat_down->Draw("hist same");
  h_stat_down->SetFillColorAlpha(kGray+2,0.4);
  h_stat_down->SetLineWidth(0);*/
  
  h_UB->Draw("E same");
  h_UB->SetLineColor(color1);
  h_UB->SetLineWidth(4);
  h_UB->SetMarkerColor(color1);
  //h_UB->SetMarkerSize(0);

  h_Empirical->Draw("E same");
  h_Empirical->SetLineColor(color2);
  h_Empirical->SetLineWidth(4);
  h_Empirical->SetMarkerColor(color2);

  h_Nieves->Draw("E same");
  h_Nieves->SetLineColor(color3);
  h_Nieves->SetLineWidth(4);
  h_Nieves->SetMarkerColor(color3);
  
  h_Susa->Draw("E same");
  h_Susa->SetLineColor(color4);
  h_Susa->SetLineWidth(4);
  h_Susa->SetMarkerColor(color4);
  
  h_Nuwro->Draw("E same");
  h_Nuwro->SetLineColor(color5);
  h_Nuwro->SetLineWidth(4);
  h_Nuwro->SetMarkerColor(color5);

  h_BNB->Draw("hist same");
  
  //legend stuff
  TLegend* legend_ratio = new TLegend(0.109, 0.5, 0.859, 0.89);
  legend_ratio->SetNColumns(2);
  legend_ratio->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C"); // option "C" allows to center the header 
  legend_ratio->AddEntry(h_BNB,Form("%d Iterations Unfolded Data",iter),"L");
  //legend_ratio->AddEntry(h_stat_up,"Stat. Unc.","F");
    /*legend_ratio->AddEntry(h_sys_up,"Systematic Unc.","F");
  legend_ratio->AddEntry(h_UB,"LFG + RPA QE + 2p2h + Tuned CCQE + Tunned CCMEC (Stat Unc.)","EL"); //microboone tune
  legend_ratio->AddEntry(h_Nieves,"LFG + RPA QE + 2p2h","EL"); //nieves
  legend_ratio->AddEntry(h_Empirical,"BRFG + Standard Model QE + Empirical MEC (Stat Unc.)","EL"); //empirical
  legend_ratio->AddEntry(h_Susa,"RMFA + Superscaling QE + Superscaling MEC (Stat Unc.)","EL"); //susa
  legend_ratio->AddEntry(h_Nuwro,"LFG + Standard Model QE + 2p2h (Stat Unc.)","EL"); //nuwro */

  legend_ratio->AddEntry(h_UB,"MicroBooNE Tune","EL"); //microboone tune
  legend_ratio->AddEntry(h_Nieves,"Nieves QE + MEC (Stat. Unc)","EL"); //nieves
  legend_ratio->AddEntry(h_Empirical,"Llewellyn QE + Empirical MEC (Stat Unc.)","EL"); //empirical
  legend_ratio->AddEntry(h_Susa,"SuSAv2 QE + MEC (Stat Unc.)","EL"); //susa
  legend_ratio->AddEntry(h_Nuwro,"NuWro QE + MEC (Stat Unc.)","EL"); //nuwro

  legend_ratio->SetLineWidth(0);
  legend_ratio->SetFillColor(kWhite);
  //legend_ratio->SetTextSize(0.027);
  legend_ratio->Draw("SAME");

  c_ratio->Print(Form("images/XSec/ratio%s.png",name));
  //c_ratio->Print(Form("images/XSec/ratio%s.pdf",name));

}


void xsec::cross_section(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D*  h_nuisance, TH1D* h_GCF,
			 TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			 TH2D* h_smearing, int iter, TH1D* h_eff_input, TH1D* h_denom_input,TH1D* h_denom_nuwro_input,
			 TH1D* h_systematic,TH2D* h_covariance,double maximum, const char* title, const char* name,bool flip_legend,bool print_contents = true){


  //gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetEndErrorSize(10);
  //gStyle->SetOptStat(0);
  
  //First we have to get the MC Vector
  /////////////////////////////////////
  std::vector<TH1D*> mc_plot_vector = mc_plots(h_empirical,h_nieves,h_susa, h_nuisance, h_GCF, print_contents);

  //Here is where the magic happens!
  //////////////////////////////////
  TH1D* h_bnb = make_bnb_plot(h_ext_input,h_dirt_input,h_overlay_total_input,h_overlay_cc2p_input,h_bnb_input,name);
  RooUnfoldResponse Response_Matrix(0,0,h_smearing,"Response_Matrix");
  RooUnfoldBayes RooUnfoldBayes_Data(&Response_Matrix, h_bnb, iter);
  TH1D* h_unfolded_signal =(TH1D* )RooUnfoldBayes_Data.Hreco();

  //One without the constants applied
  /////////////////////////////////////////
  TCanvas* c_xsec_no_const = new TCanvas("c_xsec_no_const","c_xsec_no_const",800,600);
  TH1D* h_xsec_no_const = (TH1D*)h_unfolded_signal->Clone();
  h_xsec_no_const->SetTitle(Form("%s: Unfolded Only",title));
  h_xsec_no_const->SetXTitle(Form("%s: Unfolded Only",title));
  h_xsec_no_const->SetYTitle("Arb. Units");
  h_xsec_no_const->Draw("1e1p");
  c_xsec_no_const->Print(Form("images/XSec/xsec_without_constants%s.png",name));
  
  //One with the efficiency applied
  /////////////////////////////////
  TCanvas* c_xsec_eff = new TCanvas("c_xsec_eff","c_xsec_eff",800,600);
  TH1D* h_xsec_eff = (TH1D*)h_unfolded_signal->Clone();
  TH1D* h_eff = (TH1D*) h_eff_input->Clone(); 
  h_xsec_eff->Divide(h_unfolded_signal,h_eff,1,1,"cp");
  h_xsec_eff->SetTitle(Form("%s: Efficiency Applied",title));
  h_xsec_eff->SetXTitle(Form("%s: Efficiency Applied",title));
  h_xsec_eff->SetYTitle("Arb. Units");
  h_xsec_eff->Draw("1e1p");
  c_xsec_eff->Print(Form("images/XSec/xsec_with_eff%s.png",name));

  //One with everything applied
  /////////////////////////////
  TCanvas* c_xsec = new TCanvas("c_xsec","c_xsec",2000,1500);
  //c_xsec->SetGridx();
  c_xsec->SetRightMargin(0.03);
  c_xsec->SetLeftMargin(0.13);
  c_xsec->SetBottomMargin(0.13);
  c_xsec->SetTopMargin(0.02); 
  TH1D* h_xsec = (TH1D*)h_unfolded_signal->Clone(); //extracted cross-section
  TH1D* h_denom = (TH1D*) h_denom_input->Clone(); //denominator of our efficiency aka the uboone tune prediction

  //Define a function with the bin widths
  TH1F* h_width=(TH1F*)h_xsec->Clone();
  for(int i=1; i <= h_xsec->GetNbinsX() ;i++){
      h_width->SetBinContent(i,h_xsec->GetBinWidth(i));
      h_width->SetBinError(i,0);
    }
  
  //divide BNB by efficiency
  h_xsec->Divide(h_unfolded_signal,h_eff,1,1,"cp");

  //Divide BNB by Bin Width
  h_xsec->Divide(h_width);
  h_denom->Divide(h_width);

  //scale Both MC and BNB by flux and number of targets
  double scale_value = (1)/(1E-38*N_targets*flux_value);
  h_xsec->Scale(scale_value);
  h_denom->Scale(scale_value);

  //Drawing all the MC curves
  //////////////////////////////////////////////////////////////////
  h_denom->Draw("HIST"); //Microboone Tune output aka denominator of our efficiency
  h_denom->SetLineColor(color1);
  h_denom->SetLineWidth(4);
  
  mc_plot_vector[0]->Scale(1/1E-38); //empirical
  mc_plot_vector[0]->Draw("HIST SAME");
  mc_plot_vector[0]->SetLineColor(color2);
  mc_plot_vector[0]->SetLineWidth(4);

  mc_plot_vector[1]->Scale(1/1E-38); //nieves
  mc_plot_vector[1]->Draw("hist same");
  mc_plot_vector[1]->SetLineColor(color3);
  mc_plot_vector[1]->SetLineWidth(4);

  /*  mc_plot_vector[2]->Scale(1/1E-38); //susav2: Test tag from Afro and Steven
  mc_plot_vector[2]->Draw("hist same");
  mc_plot_vector[2]->SetLineColor(color4);
  mc_plot_vector[2]->SetLineWidth(4);*/

  mc_plot_vector[3]->Scale(1/1E-38); //SuSAv2: GENIE v3.2.0 tag
  mc_plot_vector[3]->Draw("hist same");
  mc_plot_vector[3]->SetLineColor(kRed);
  mc_plot_vector[3]->SetLineWidth(4);

  //mc_plot_vector[4]->Scale(1/1E-38); //GCF
  //mc_plot_vector[4]->Draw("hist same");
  //mc_plot_vector[4]->SetLineColor(color5);
  //mc_plot_vector[4]->SetLineWidth(4);

  TH1D* h_denom_nuwro=(TH1D*)h_denom_nuwro_input->Clone(); //nuwro prediction
  h_denom_nuwro->Divide(h_width);
  h_denom_nuwro->Scale(scale_value);
  h_denom_nuwro->Draw("HIST SAME");
  h_denom_nuwro->SetLineColor(color5);
  h_denom_nuwro->SetLineWidth(8);
  h_denom_nuwro->SetLineStyle(10);
 
  //Now we can draw the data and the systematic
  ///////////////////////////////////////
  h_xsec->Draw("1e1p SAME");
  h_xsec->SetLineColor(kBlack);
  h_xsec->SetLineWidth(4);
  h_xsec->SetMarkerStyle(8);
  h_xsec->SetMarkerSize(2);

  TH1D* h_xsec_sys = (TH1D*)h_xsec->Clone();
  TH1D* h_xsec_CV = (TH1D*)h_xsec->Clone();
  Get_Systematic(h_xsec_sys,h_xsec_CV,h_covariance,true); //gets the systematic from the covariance matrix
  h_xsec_sys->Draw("1e1p SAME");
  h_xsec_sys->SetLineColor(kBlack);
  h_xsec_sys->SetLineWidth(4);
  h_xsec_sys->SetMarkerStyle(8);
  h_xsec_sys->SetMarkerSize(2);

  //Set titles and stuff cause root is dumb
  ////////////////////////////////////////
  //h_denom->SetTitle(Form("%s", title));
  h_denom->SetTitle("");
  h_denom->SetMaximum(maximum);
  h_denom->SetMinimum(0);
  h_denom->SetXTitle(Form("%s",title));
  h_denom->GetXaxis()->SetTitleSize(60); //35
  h_denom->GetXaxis()->SetTitleFont(43);
  h_denom->GetXaxis()->SetTitleOffset(1.3);
  h_denom->GetXaxis()->SetLabelFont(43); //4 = hevelatica normal 3 = precision
  h_denom->GetXaxis()->SetLabelSize(60);
  if(std::strcmp(name,"_muon_mom") ==  0 || std::strcmp(name,"_leading_mom") ==  0 || std::strcmp(name,"_recoil_mom") ==  0 || std::strcmp(name,"_delta_PT") == 0){ unit_title = "GeV/c";}
  else{ unit_title ="";}
  h_denom->SetYTitle(Form("Differential Cross-Section [#frac{10^{-38} cm^{2}}{%s Ar}]",unit_title));
  h_denom->GetYaxis()->SetTitleSize(60);
  h_denom->GetYaxis()->SetTitleFont(43); //4 = hevelatica normal 3 = precision
  h_denom->GetYaxis()->SetTitleOffset(1.4);
  h_denom->GetYaxis()->SetLabelFont(43);
  h_denom->GetYaxis()->SetLabelSize(60);
  
  //legend stuff
  /////////////////////////
  TLegend* legend_xsec;
  int nbins = h_xsec->GetNbinsX();
  if(flip_legend){
    legend_xsec = new TLegend(0.16, 0.56, 0.92, 0.96);
  }else{
    legend_xsec = new TLegend(0.16, 0.56, 0.92, 0.96);
  }
  legend_xsec->SetHeader("#bf{MicroBooNE 6.79 x 10^{20} POT}","C"); // option "C" allows to center the header 
  legend_xsec->AddEntry(h_xsec,Form("%d Iterations Unfolded Data (Stat. + Sys.)",iter),"p");
  double chi_uboone = 0;
  double chi_emp = 0;
  double chi_nieves = 0;
  double chi_susa = 0;
  double chi_nuwro = 0;
  legend_xsec->AddEntry(h_denom,Form("GENIE MicroBooNE Tune: Tuned Nieves QE + MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(h_denom,h_xsec,h_covariance,chi_uboone),nbins),"L");
  legend_xsec->AddEntry(mc_plot_vector[0],Form("GENIE Empirical: Llewellyn Smith QE + Empirical MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(mc_plot_vector[0],h_xsec,h_covariance,chi_emp),nbins),"L");
  legend_xsec->AddEntry(mc_plot_vector[1],Form("GENIE Nieves: Nieves QE + MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(mc_plot_vector[1],h_xsec,h_covariance,chi_nieves),nbins),"L");
  //legend_xsec->AddEntry(mc_plot_vector[2],Form("GENIE SuSAv2: SuSAv2 QE +MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(mc_plot_vector[2], h_xsec, h_covariance,chi_susa),nbins),"L"); //test version from Steven and Afro. Depricated 9/21/22
  legend_xsec->AddEntry(mc_plot_vector[3],Form("GENIE SuSAv2: SuSAv2 QE +MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(mc_plot_vector[3], h_xsec, h_covariance,chi_susa),nbins),"L"); //what I call "nuisance", but is really SuSAv2 from GENIE v3.2.0
  legend_xsec->AddEntry(h_denom_nuwro,Form("NuWro: Llewellyn Smith QE + Nieves MEC (#chi^{2}/DoF: %4.1f/%d)",CalcChiSquared(h_denom_nuwro, h_xsec, h_covariance,chi_nuwro),nbins),"L");
  legend_xsec->SetLineWidth(0);
  legend_xsec->SetFillColor(kWhite);
  legend_xsec->SetMargin(0.1); //sets margin for legend
  legend_xsec->Draw("SAME");

  TText *t1 = new TText();
  t1->SetNDC();
  t1->SetTextFont(43);
  t1->SetTextSize(60);
  t1->SetTextAlign(22);
  t1->SetTextAngle(0);
  
  if(std::strcmp(title,"Unfolded cos(#gamma_{#vec{p}_{L}, #vec{p}_{R}})") == 0){
    t1->DrawText(0.95,0.95,"(a)");
  }else if(std::strcmp(title,"Unfolded cos(#gamma_{#vec{p}_{#mu}, #vec{p}_{sum}})") == 0){
    t1->DrawText(0.95,0.95,"(b)");
  }
  
  c_xsec->Print(Form("images/XSec/xsec%s.png",name));
  delete t1;

  Ratio_Plot(h_xsec,h_denom,  mc_plot_vector[0],  mc_plot_vector[1], mc_plot_vector[2], h_denom_nuwro, h_systematic, iter, title,name);

  std::cout<<"-------------------------------------------"<<std::endl;
  std::cout<<"Variable Name: "<< title << " Number of Bins: "<< nbins <<std::endl; 
  for(int i=1; i < nbins+1; i++){
    double xlow = h_xsec->GetBinLowEdge(i);
    double xhigh = h_xsec->GetBinLowEdge(i+1);
    double bin_xsec = h_xsec->GetBinContent(i);
    double bin_statistical = h_xsec->GetBinError(i);
    double bin_sys = h_xsec_sys->GetBinError(i);
    double sys_total = std::pow((std::pow(h_xsec->GetBinError(i),2)+std::pow(h_xsec_sys->GetBinError(i),2)),0.5);
    std::cout<<"Bin Number: "<<i<<Form(" Bin Range: [%f-%f]",xlow,xhigh)<<" Value of XSec in this Bin: "<<bin_xsec<<" Value of Stat. "<<bin_statistical<<" Value of Sys: "<<bin_sys<<" Value of Total Uncertainty: "<<sys_total<<std::endl;
  }
  std::cout<<"-------------------------------------------"<<std::endl;


  h_xsec->Write(Form("h_data_xsec%s",name));
  h_denom->Write(Form("h_uboone_tune_xsec%s",name));
  mc_plot_vector[0]->Write(Form("h_empirical_xsec%s",name));
  mc_plot_vector[1]->Write(Form("h_nieves_xsec%s",name));
  mc_plot_vector[2]->Write(Form("h_susav2_xsec%s",name));
  mc_plot_vector[3]->Write(Form("h_nuisance_xsec%s",name));
  h_xsec_sys->Write(Form("h_total_systematic_uncertainty%s",name));
  TH2D* h_covar = (TH2D*)h_covariance->Clone();
  h_covar->Write(Form("h_2D_total_covariance_matrix%s",name));

  mc_plot_vector.clear();
	     
}

#endif

