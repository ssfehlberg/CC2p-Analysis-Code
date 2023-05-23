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
  virtual std::vector<TH1D*> mc_plots(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_GCF, bool print_contents = true);
  virtual void Fix_Systematic(TH1D* h,TH1D* h_MC_CV,TH1D* h_systematic, bool print_contents);
  virtual void cross_section(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_GCF,
			     TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			     TH2D* h_smearing, int iter, TH1D* h_eff_input, TH1D* h_denom_input,TH1D* h_denom_nuwro_input,
			     TH1D* h_systematic,double maximum, const char* title, const char* name,bool flip_legend,bool print_contents = true);


private:

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
  TFile* f_GCF = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/GCF/hists_GCF_CCQE_fsi.root");
  
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
    h_GCF_muon[i] = (TH1D*)f_GCF->Get(Form("h_muon%s_lead_cut",var0[i])); //gcf

    h_nuwro_leading[i] = (TH1D*)f_nuwro->Get(Form("h_leading%s_total",var0[i])); //nuwro
    h_empirical_leading[i] = (TH1D*)f_empirical->Get(Form("h_leading%s_lead_cut",var0[i])); //empirical
    h_nieves_leading[i] = (TH1D*)f_nieves->Get(Form("h_leading%s_lead_cut",var0[i])); //nieves
    h_susa_leading[i] = (TH1D*)f_susa->Get(Form("h_leading%s_lead_cut",var0[i])); //susa
    h_GCF_leading[i] = (TH1D*)f_GCF->Get(Form("h_leading%s_lead_cut",var0[i])); //gcf

    h_nuwro_recoil[i] = (TH1D*)f_nuwro->Get(Form("h_recoil%s_total",var0[i])); //nuwro
    h_empirical_recoil[i] = (TH1D*)f_empirical->Get(Form("h_recoil%s_lead_cut",var0[i])); //empirical
    h_nieves_recoil[i] = (TH1D*)f_nieves->Get(Form("h_recoil%s_lead_cut",var0[i])); //nieves
    h_susa_recoil[i] = (TH1D*)f_susa->Get(Form("h_recoil%s_lead_cut",var0[i])); //susa
    h_GCF_recoil[i] = (TH1D*)f_GCF->Get(Form("h_recoil%s_lead_cut",var0[i])); //gcf

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
    h_GCF_other[i] = (TH1D*)f_GCF->Get(Form("h%s_lead_cut",other_var[i])); //gcf

    //MC Models: Truth level for neutrino energy and pn
    if(i == 5 || i == 6){

      h_empirical_other_true[i] = (TH1D*)f_empirical->Get(Form("h%s_lead_cut",other_var[i])); //empirical
      h_nieves_other_true[i] = (TH1D*)f_nieves->Get(Form("h%s_lead_cut",other_var[i])); //nieves
      h_susa_other_true[i] = (TH1D*)f_susa->Get(Form("h%s_lead_cut",other_var[i])); //susa
      h_GCF_other_true[i] = (TH1D*)f_GCF->Get(Form("h%s_lead_cut",other_var[i])); //gcf  

      h_empirical_other[i] = (TH1D*)f_empirical->Get(Form("h%s_true_lead_cut",other_var[i])); //empirical
      h_nieves_other[i] = (TH1D*)f_nieves->Get(Form("h%s_true_lead_cut",other_var[i])); //nieves
      h_susa_other[i] = (TH1D*)f_susa->Get(Form("h%s_true_lead_cut",other_var[i])); //susa
      h_GCF_other[i] = (TH1D*)f_GCF->Get(Form("h%s_true_lead_cut",other_var[i])); //gcf  
    }
    
  }

  //Systematic Uncertainty
  ////////////////////////
  TFile* f_systematic = new TFile("Systematics/root_files/total_error.root");

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

  TLatex* t;
  t->DrawLatex(0.3,0.92,Form("%s",pot_num));
  t->DrawLatex(0.82,0.92,Form("%s",sample_name));
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

  t->DrawLatex(0.3,0.92,Form("%s",pot_num));
  t->DrawLatex(0.82,0.92,Form("%s",sample_name));
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

std::vector<TH1D*> xsec::mc_plots(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_GCF, bool print_contents = true){

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
      
      std::cout<<"n_GCF: "<<n_GCF<<std::endl;
      std::cout<<"Value_GCF: "<<value_GCF<<std::endl;

    } 
  }

  mc_vector.push_back(h_empirical);
  mc_vector.push_back(h_nieves);
  mc_vector.push_back(h_susa);
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

void xsec::cross_section(TH1D* h_empirical, TH1D* h_nieves, TH1D* h_susa, TH1D* h_GCF,
			 TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			 TH2D* h_smearing, int iter, TH1D* h_eff_input, TH1D* h_denom_input,TH1D* h_denom_nuwro_input,
			 TH1D* h_systematic,double maximum, const char* title, const char* name,bool flip_legend,bool print_contents = true){


  std::cout<<"Name at beginning: "<<name<<std::endl;
  
  //First we have to get the MC Vector
  /////////////////////////////////////
  std::vector<TH1D*> mc_plot_vector = mc_plots(h_empirical,h_nieves,h_susa,h_GCF,print_contents);

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

  std::cout<<"Name after no Constants: "<<name<<std::endl;
  
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
  c_xsec->SetGridx();
  TH1D* h_xsec = (TH1D*)h_unfolded_signal->Clone();
  h_xsec->Draw("1e1p");

  //drawing denom on the same plot
  TH1D* h_denom = (TH1D*) h_denom_input->Clone();
  h_denom->Draw("HIST SAME");
  h_denom->SetLineColor(color1);
  h_denom->SetLineWidth(4);

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

  //Setting all the values now that we are done with all the manipulation
  //////////////////////////////////////////////////////////////////
  h_xsec->SetLineColor(kBlack);
  h_xsec->SetLineWidth(4);
  h_xsec->SetMarkerSize(1);
  h_xsec->SetTitle(Form("%s", title));
  h_xsec->SetMaximum(maximum);
  h_xsec->SetMinimum(0);

  h_xsec->SetXTitle(Form("%s",title));
  h_xsec->GetXaxis()->SetTitleSize(50); //35                                                                                                                                                                                                                              
  h_xsec->GetXaxis()->SetTitleFont(43);
  h_xsec->GetXaxis()->SetTitleOffset(1.3);
  h_xsec->GetXaxis()->SetLabelFont(43);
  h_xsec->GetXaxis()->SetLabelSize(50);

  h_xsec->SetYTitle("Differential Cross-Section [10^{-38} cm^{2} / Argon]");
  h_xsec->GetYaxis()->SetTitleSize(50);
  h_xsec->GetYaxis()->SetTitleFont(43); //4 = hevelatica normal 3 = precision                                                                                                                                                                                             
  h_xsec->GetYaxis()->SetTitleOffset(1.3);
  h_xsec->GetYaxis()->SetLabelFont(43);
  h_xsec->GetYaxis()->SetLabelSize(50);

  //Systematics
  TH1D* h_xsec_sys = (TH1D*)h_xsec->Clone();
  TH1D* h_xsec_CV = (TH1D*)h_xsec->Clone();

  Fix_Systematic(h_xsec_sys,h_xsec_CV,h_systematic,true);
  h_xsec_sys->Draw("1e1p SAME");
  h_xsec_sys->SetLineColor(kBlack);
  h_xsec_sys->SetLineWidth(4);
  h_xsec_sys->SetMarkerSize(1);

  //drawing the different models on the same plot
  //////////////////////////////////////////////
  mc_plot_vector[0]->Draw("hist same"); //empirical
  mc_plot_vector[0]->SetLineColor(color2);
  mc_plot_vector[0]->SetLineWidth(4);
  mc_plot_vector[0]->Scale(1/1E-38);

  mc_plot_vector[1]->Draw("hist same"); //nieves
  mc_plot_vector[1]->SetLineColor(color3);
  mc_plot_vector[1]->SetLineWidth(4);
  mc_plot_vector[1]->Scale(1/1E-38);

  mc_plot_vector[2]->Draw("hist same"); //susa
  mc_plot_vector[2]->SetLineColor(color4);
  mc_plot_vector[2]->SetLineWidth(4);
  mc_plot_vector[2]->Scale(1/1E-38);

  //mc_plot_vector[3]->Draw("hist same"); GCF
  //mc_plot_vector[3]->SetLineColor(color5);
  //mc_plot_vector[3]->SetLineWidth(4);
  //mc_plot_vector[3]->Scale(1/1E-38);

  TH1D* h_denom_nuwro=(TH1D*)h_denom_nuwro_input->Clone();
  h_denom_nuwro->Divide(h_width);
  h_denom_nuwro->Scale(scale_value);
  h_denom_nuwro->Draw("hist same");
  h_denom_nuwro->SetLineColor(color5);
  h_denom_nuwro->SetLineWidth(4);
  
  
  //legend stuff
  TLegend* legend_xsec;
  if(flip_legend){
    legend_xsec = new TLegend(0.105, 0.5, 0.584, 0.89);
  }else{
    legend_xsec = new TLegend(0.37, 0.5, 0.859, 0.89);
  }
  legend_xsec->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C"); // option "C" allows to center the header 
  legend_xsec->AddEntry(h_xsec,Form("%d Iterations Unfolded Data (Stat. + Sys.)",iter),"lepf");
  legend_xsec->AddEntry(h_denom,"MicroBooNE Tune","L");
  legend_xsec->AddEntry(mc_plot_vector[0],"Empirical MEC + Lwellyn Smith QE + FSI","L");
  legend_xsec->AddEntry(mc_plot_vector[1],"Nieves MEC & QE + FSI","L");
  legend_xsec->AddEntry(mc_plot_vector[2],"SuSav2 MEC + QE + FSI","L");
  //legend_xsec->AddEntry(mc_plot_vector[3],"GCF + FSI","L");
  legend_xsec->AddEntry(h_denom_nuwro,"NuWro","L");
  legend_xsec->SetLineWidth(0);
  legend_xsec->SetFillColor(kWhite);
  legend_xsec->SetTextSize(0.03);
  legend_xsec->Draw("SAME");

  std::cout<<"Name: "<<name<<std::endl;
  c_xsec->Print(Form("images/XSec/xsec%s.png",name));
  mc_plot_vector.clear();
	     
}

#endif

