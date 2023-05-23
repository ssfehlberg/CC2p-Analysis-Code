#ifndef analysis_h
#define analysis_h

#include "tools/helper_funcs.h"
#include "tools/paul_tol_colors.hpp"
#include "tools/constants.h"
using namespace Constants;

#include <iostream>
#include <ctime>
#include <string>

//////////////////////////////
// Functions to perform plotting are defined here.
// Some helper functions can be found in helper_funcs.h
// Histograms are also defined here to create global variables
//////////////////////////////

class nuwro_analysis{

public:
  virtual void main();


private:
  virtual const char* which_sample();
  virtual std::pair<const char*, const char*> which_run();
  virtual void Define_Parameters(const char* run,const char* sample);
  virtual void Grab_Histograms(TFile* f_overlay,  TFile* f_bnb,  TFile* f_ext, TFile* f_dirt, TFile* f_eff, TFile* f_mom_thresholds);
  virtual void Plot_Histograms(const char* run, char const* pot_num, const char* sample_name,Color_t colors[], std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2,  TH1D* h_statistical, TH1D* h_systematic, TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1, TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend, int ymax, int ymin, int num_channels, const char* titles, string path,  TLine* a2,const char* titles2 = "", const char* plots="", const char* cut = "", bool plot_total = false, bool plot_ccnue = false, bool flip_legend = false,double pad_lim = 0.0, double pad0_lim = 0.19,double bnb_min = 0.0, double bnb_max = 2.7, Double_t xlim = -9999.0, Double_t xlim_up = +9999.0);
  
}; //end of class

#endif
#ifdef nuwro_analysis_cxx

//////////////////////////
//FUNCTIONS
////////////////////////

const char* nuwro_analysis::which_sample(){

  //Get the correct Sample                                                                                                                                                                                                                                                     
  char response;
  const char* Sample;
  std::cout<<"Which Sample is This?"<<std::endl;
  std::cout<<" 0 = nuwro \n 1 = nuwro with xsec binning"<<std::endl;
  std::cin>>response;

  if(response == '0'){
    Sample = "";
  }else if (response == '1'){
    Sample = "_xsec";
  }else{
    std::cout<<"Invalid Response. Please Type 0, 1, or 2 for pelee, filtered, and unfiltered samples respectively."<<std::endl;
  }

  return Sample;

}

std::pair<const char*, const char*> nuwro_analysis::which_run(){

  //Get the Correct Run                                                                                                                                                                                                                                                        
  char response1;
  const char* Run;
  const char* Run_Title;
  std::cout<<"Which Run is This?"<<std::endl;
  std::cout<<" 0 = Jan \n 1 = Run 1 \n 2 = Run 2 \n 3 = Run 3 \n 4 = Runs 1+2+3"<<std::endl;
  std::cin>>response1;

  if(response1 == '0'){
    Run = "Jan";
    Run_Title = "January Sample";
  } else if(response1 == '1'){
    Run = "Run1";
    Run_Title = "Run 1";
  } else if(response1 == '2') {
    Run = "Run2";
    Run_Title = "Run 2";
  } else if (response1 == '3'){
    Run = "Run3";
    Run_Title = "Run 3";
  }else if(response1 == '4'){
    Run = "Run_all";
    Run_Title = "Runs 1+2+3";
  }else{
      std::cout<<"Invalid Response. Please Type 0, 1, 2 ,3, 4, or 5 for January, Run 1, Run2 , Run 3, Runs 1+2+3, respectively."<<std::endl;
  }

  return std::make_pair(Run,Run_Title);

}

void nuwro_analysis::Define_Parameters(const char* run, const char* sample){

  //latex and style stuff
  ////////////////////////
  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  gStyle->SetTitleSize(0.1);
  gStyle->SetTitleAlign(23);
  gStyle->SetEndErrorSize(10);
  t->SetNDC();
  t->SetTextAlign(22);

  //POT number and In-Progress
  ////////////////////////////
  if(strcmp(run,"Jan") == 0){
    pot_num="#scale[0.6]{January Sample Accumulated POT: 4.54e+19}";//pot number printed on the plots: January sample
    run_num = 0;
  }else if(strcmp(run,"Run1") == 0){
    pot_num="#scale[0.6]{Run 1 Accumulated POT: 1.62e+20}";//pot number printed on the plots: Run 1
    run_num = 1;
  }else if(strcmp(run,"Run2") == 0){
    pot_num="#scale[0.6]{Run 2 Accumulated POT: 2.62e+20}";//pot number printed on the plots: Run 2
    run_num = 2;
  }else if(strcmp(run,"Run3") == 0){
    pot_num="#scale[0.6]{Run 3 Accumulated POT: 2.55e+20}";//pot number printed on the plots: Run 3
    run_num = 3;
  }else if(strcmp(run,"Run_all") == 0 &&  strcmp(sample,"") == 0){
    pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";//pot number printed on the plots: Run 3
    run_num = 4;
  } else if(strcmp(run,"Run_all") == 0 &&  strcmp(sample,"_xsec") == 0){
    pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";//pot number printed on the plots: Run 3
    run_num = 5;
    }

  //sample_name="#scale[0.6]{MicroBooNE In-Progress}";//sample name printed on the plots
  sample_name="#scale[0.6]{MicroBooNE Preliminary}";//sample name printed on the plots 
  
} //end of define parameters


void nuwro_analysis::Grab_Histograms( TFile* f_overlay,  TFile* f_bnb,  TFile* f_ext, TFile* f_dirt, TFile* f_eff, TFile* f_mom_thresholds){

  for(int i = 0; i < num_cuts; i++){
    for(int j = 0; j < num_variables; j++){
      h_bnb[i][j][0] = (TH1D*)f_bnb->Get(Form("h%s%s_bnb",plots[j],cut[i]));
      h_ext[i][j][0] = (TH1D*)f_ext->Get(Form("h%s%s_ext",plots[j],cut[i]));
      h_dirt[i][j][0] = (TH1D*)f_dirt->Get(Form("h%s%s_dirt_wgt",plots[j],cut[i]));
      for(int k = 0; k < num_channels; k++){
	h_overlay[i][j][k] = (TH1D*)f_overlay->Get(Form("h%s%s%s",plots[j],cut[i],channel[k]));
      }
      for(int k=0; k < num_channels_raquel; k++){
	h_overlay_raquel[i][j][k] = (TH1D*)f_overlay->Get(Form("h%s_raquel%s%s",plots[j],cut[i],channel_raquel[k]));
      }
    }
    
    //grabbing truth stuff
    for(int j=0; j < num_truth; j++){
      for(int k=0; k < num_channels; k++){
	h_mc[i][j][k] = (TH1D*)f_overlay->Get(Form("h%s%s%s",truth[j],cut[i],channel[k]));
      }
    } 
  }

  //grabbing things related the PFP's
  for(int i=0; i < num_group; i++){
    h_overlay_pfp[i] = (TH1D*)f_overlay->Get(Form("h_%s_overlay",group[i]));
    h_bnb_pfp[i] = (TH1D*)f_bnb->Get(Form("h_%s_bnb",group[i]));
    h_ext_pfp[i] =(TH1D*)f_ext->Get(Form("h_%s_ext",group[i]));
  }
  
  //grabbing the 2D histograms
  for(int i=0; i < num_group2d; i++){
    h_overlay2D[i] = (TH2D*)f_overlay->Get(Form("h_correlation_overlay_%s",group2d[i])); 
  }

  //grabbing efficiency plots:
  for(int i=0; i < num_eff; i++){
    h_num[i] =  (TH1D*)f_mom_thresholds->Get(Form("h_mom_threshold_num_%s",eff[i]));
    h_denom[i] =  (TH1D*)f_mom_thresholds->Get(Form("h_mom_threshold_denom_%s",eff[i]));
  }

  eff_graph = (TGraph*)f_overlay->Get("eff_graph");
  pur_graph = (TGraph*)f_overlay->Get("pur_graph");

   for(int i=0; i < num_particles_eff; i++){
      for(int j=0; j < num_particles_eff_plots; j++){
        h_particle_num[i][j] = (TH1D*)f_eff->Get(Form("h_particle_num%s%s",particles_eff[i],particles_eff_var[j]));
	h_particle_denom[i][j] = (TH1D*)f_eff->Get(Form("h_particle_denom%s%s",particles_eff[i],particles_eff_var[j]));

      }
    }

  for(int i = 0; i < num_other_eff; i++){
    h_other_eff_num[i] = (TH1D*)f_eff->Get(Form("h_other_eff_num%s",other_eff[i]));
    h_other_eff_denom[i] = (TH1D*)f_eff->Get(Form("h_other_eff_denom%s",other_eff[i]));
  }
     
  //random track variables
  for(int i =0; i < num_track; i++){
    for(int k=0; k < track_cut; k++){
      h_track_bnb[i][k][0] = (TH1D*)f_bnb->Get(Form("h_track%s%s",variable[i],which_track_cut[k]));
      h_track_ext[i][k][0] = (TH1D*)f_ext->Get(Form("h_track%s%s",variable[i],which_track_cut[k]));
      h_track_dirt[i][k][0] = (TH1D*)f_dirt->Get(Form("h_track%s%s",variable[i],which_track_cut[k]));
      for(int j = 0; j < num_particles; j++){
	h_track_overlay[i][k][j] = (TH1D*)f_overlay->Get(Form("h_track%s%s%s",variable[i],which_track_cut[k],particles[j]));
      }
    }
  }
  
  //grabbing particle plots
  for(int i = 0; i < num_var; i++){
    h_muon_ext[i][0] = (TH1D*)f_ext->Get(Form("h_muon%s_ext",var[i]));
    h_recoil_ext[i][0] = (TH1D*)f_ext->Get(Form("h_recoil%s_ext",var[i]));
    h_leading_ext[i][0] = (TH1D*)f_ext->Get(Form("h_leading%s_ext",var[i]));
    h_muon_bnb[i][0] = (TH1D*)f_bnb->Get(Form("h_muon%s_bnb",var[i]));
    h_recoil_bnb[i][0] = (TH1D*)f_bnb->Get(Form("h_recoil%s_bnb",var[i]));
    h_leading_bnb[i][0] = (TH1D*)f_bnb->Get(Form("h_leading%s_bnb",var[i]));
    h_muon_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_muon%s_dirt_wgt",var[i]));
    h_recoil_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_recoil%s_dirt_wgt",var[i]));
    h_leading_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h_leading%s_dirt_wgt",var[i]));

    for(int j =0; j < num_channels; j++){
      h_muon_overlay[i][j] = (TH1D*)f_overlay->Get(Form("h_muon%s%s",var[i],channel[j]));
      h_recoil_overlay[i][j] = (TH1D*)f_overlay->Get(Form("h_recoil%s%s",var[i],channel[j]));
      h_leading_overlay[i][j] = (TH1D*)f_overlay->Get(Form("h_leading%s%s",var[i],channel[j]));
    }

    for(int j=0; j < num_channels_raquel; j++){
      h_muon_overlay_raquel[i][j] = (TH1D*)f_overlay->Get(Form("h_muon_raquel%s%s",var[i],channel_raquel[j]));
      h_recoil_overlay_raquel[i][j] = (TH1D*)f_overlay->Get(Form("h_recoil_raquel%s%s",var[i],channel_raquel[j]));
      h_leading_overlay_raquel[i][j] = (TH1D*)f_overlay->Get(Form("h_leading_raquel%s%s",var[i],channel_raquel[j]));

    }
    
  }

  //grabbing physics plots
  for(int i=0; i < num_phys; i++){
    //std::cout<<"Physics[i]: "<<physics[i]<<std::endl;
    h_phys_ext[i][0] = (TH1D*)f_ext->Get(Form("h%s_ext",physics[i]));
    h_phys_bnb[i][0] = (TH1D*)f_bnb->Get(Form("h%s_bnb",physics[i]));
    h_phys_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h%s_dirt_wgt",physics[i]));
    for(int j=0; j < num_channels; j++){
      //std::cout<<"Channel[j]: "<<channel[j]<<std::endl;
      h_phys_overlay[i][j] = (TH1D*)f_overlay->Get(Form("h%s%s",physics[i],channel[j])); 
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_phys_overlay_raquel[i][j] = (TH1D*)f_overlay->Get(Form("h%s_raquel%s",physics[i],channel_raquel[j]));
    }  
  }

  //grabbing stvs plots
  for(int i=0; i < num_stv; i++){
    h_stv_ext[i][0] = (TH1D*)f_ext->Get(Form("h%s_ext",stv[i]));
    h_stv_bnb[i][0] = (TH1D*)f_bnb->Get(Form("h%s_bnb",stv[i]));
    h_stv_dirt[i][0] = (TH1D*)f_dirt->Get(Form("h%s_dirt_wgt",stv[i]));
    for(int j=0; j < num_channels; j++){
      h_stv_overlay[i][j] = (TH1D*)f_overlay->Get(Form("h%s%s",stv[i],channel[j])); 
    }
   for(int j=0; j < num_channels_raquel; j++){
      h_stv_overlay_raquel[i][j] = (TH1D*)f_overlay->Get(Form("h%s_raquel%s",stv[i],channel_raquel[j]));
    }    
  }

  //Grabbing the systematic and statistical uncertainty plots
  ////////////////////////////////////////////////////////////
  TFile* f_stat = new TFile("../Systematics/root_files/Statistical/systematics.root");
  TFile* f_systematic = new TFile("../Systematics/root_files/total_error.root");
  
 for(int i = 0; i < num_var; i++){
    h_muon_stat[i] = (TH1D*)f_stat->Get(Form("hist_fractional_errors_muon%s",var[i]));
    h_muon_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_muon%s",var[i]));
    h_leading_stat[i] = (TH1D*)f_stat->Get(Form("hist_fractional_errors_leading%s",var[i]));
    h_leading_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_leading%s",var[i]));
    h_recoil_stat[i] = (TH1D*)f_stat->Get(Form("hist_fractional_errors_recoil%s",var[i]));
    h_recoil_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors_recoil%s",var[i]));
 }

 std::cout<<"Grabbed the particles"<<std::endl;
 
 for(int i=0; i < num_phys; i++){
    h_phys_stat[i] = (TH1D*)f_stat->Get(Form("hist_fractional_errors%s",physics[i]));
    h_phys_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors%s",physics[i]));
 }
 
  std::cout<<"Grabbed the physics"<<std::endl;

 for(int i=0; i < num_stv; i++){
   h_stv_stat[i] = (TH1D*)f_stat->Get(Form("hist_fractional_errors%s",stv[i]));
   h_stv_systematic[i] = (TH1D*)f_systematic->Get(Form("hist_fractional_errors%s",stv[i]));
 }

  std::cout<<"Grabbed the stvs"<<std::endl;
 
}

void nuwro_analysis::Plot_Histograms(char const* run, char const* pot_num, const char* sample_name,Color_t colors[],std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2, TH1D* h_statistical, TH1D* h_systematic, TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1,TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend,int ymax, int ymin, int num_channels, const char* titles, string path,  TLine* a2,const char* titles2 ="", const char* plots="", const char* cut = "",bool plot_total = false, bool plot_ccnue = false, bool flip_legend = false, double pad_lim = 0.0, double pad0_lim = 0.19, double bnb_min = 0.0, double bnb_max = 2.7, Double_t xlim = -9999.,Double_t xlim_up = +9999.0)
{      

  h_ext1 = (TH1D*)h_ext->Clone();
  h_ext2 = (TH1D*)h_ext->Clone();
  h_dirt1 = (TH1D*)h_dirt->Clone();
  h_dirt2 = (TH1D*)h_dirt->Clone();
  h_bnb1 = (TH1D*)h_bnb->Clone();
  h_overlay00 = (TH1D*)h_overlay[0]->Clone();
  h_overlay1 = (TH1D*)h_overlay[0]->Clone();
  h_overlay2 = (TH1D*)h_overlay[0]->Clone();
  
  canv = new TCanvas(Form("C%s%s",plots,cut),Form("C%s%s",plots,cut),2000,1500);
  canv->SetGridx();
  h = new THStack(Form("h%s%s",plots,cut),Form("h%s%s",plots,cut));
   
  //Stacked Histrogram parameters
  h->Draw("HIST");
  h->SetTitle("");
  h->SetMaximum(ymax); 
  h->SetMinimum(ymin);
  
  //MC
  int z;
  int f;
  if(plot_ccnue){
    z = num_channels;
    f = num_channels;
  }else{
    z = num_channels-1;
    f = num_channels-1;
  }
          
  for(int k=1; k < z ; k++){
    if(k%2 != 0){
      h_overlay[k]->SetLineColor(colors[k]);
      h_overlay[k]->SetFillColor(colors[k]);
      h_overlay[k]->SetLineWidth(1);
      h->Add(h_overlay[k]);
      }
  }
  
   for(int k=1; k < z ; k++){
    if(k%2 == 0){
      h_overlay[k]->SetLineColor(colors[k]);
      h_overlay[k]->SetFillColor(colors[k]);
      h_overlay[k]->SetLineWidth(1);
      h->Add(h_overlay[k]);
    }
  }
  
  //Dirt
  h->Add(h_dirt);
  h_dirt->SetFillColor(kOrange-8);
  h_dirt->SetLineColor(kOrange-8);
  h_dirt->SetLineWidth(1);
  
  //EXT
  h->Add(h_ext);
  h_ext->SetFillColor(kViolet-7);
  h_ext->SetFillStyle(3005);
  h_ext->SetLineWidth(1);
  
  //BNB
  /* std::cout<<"Right before systematic"<<std::endl;
  Fix_Systematic(h_bnb1,h_overlay00,h_systematic);
  h_bnb1->Draw("1e1p SAME");
  h_bnb1->SetLineColor(kBlack);
  h_bnb1->SetLineWidth(4);
  h_bnb1->SetMarkerSize(1);
  */
  
  h_bnb->Draw("1e1p SAME");
  h_bnb->SetLineColor(kBlack);
  h_bnb->SetLineWidth(4);
  h_bnb->SetMarkerSize(1);

  h->GetXaxis()->SetTitle(Form("%s",titles));
  h->GetXaxis()->SetTitleSize(50); //35
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(50);

  h->GetYaxis()->SetTitle("No. Events");
  h->GetYaxis()->SetTitleSize(50);
  h->GetYaxis()->SetTitleFont(43); //4 = hevelatica normal 3 = precision
  h->GetYaxis()->SetTitleOffset(1.3);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(50);
  
  //Drawing cut lines if needed
  a2 = new TLine(xlim,ymin,xlim,ymax);  
  a2->Draw("SAME");
  a2->SetLineColor(kBlack);
  a2->SetLineWidth(3);
  
  a3 = new TLine(xlim_up,ymin,xlim_up,ymax);  
  a3->Draw("SAME");
  a3->SetLineColor(kBlack);
  a3->SetLineWidth(3);
  
  //if you want to plot the total for sanity sake:
  if(plot_total){
    h_overlay[0]->Draw("SAME");
  }

  //Make sure to plot the statistical uncertainty
  Fix_Statistical(h_overlay1,h_overlay00,h_ext1,h_dirt1);
  h_overlay1->Draw("e2SAME");
  h_overlay1->SetLineColor(kBlack);
  h_overlay1->SetFillColor(kBlack);
  h_overlay1->SetFillStyle(3004);
  h_overlay1->SetMarkerSize(0);
  h_overlay1->SetLineWidth(1);

  if(flip_legend){
    legend = new TLegend(0.075, 0.56, 0.584, 0.87);
  }else{
    legend = new TLegend(0.37, 0.56, 0.889, 0.87);
  }

  legend->SetNColumns(2);
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C"); // option "C" allows to center the header
  legend->AddEntry(h_bnb,"BNB Data (Stat)","lepf");
  legend->AddEntry(h_ext,"EXT Data","f");
  legend->AddEntry(h_overlay1,"MC Stat. Unc.","f");
  legend->AddEntry(h_dirt,"Dirt","f");
  for(int k =1; k < z; k++){	
    legend->AddEntry(h_overlay[k],Form("%s",channel_legend[k]),"f");
  }
  legend->SetLineWidth(2);
  legend->SetLineColor(kGray);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.03);
  legend->Draw("same");
  t->SetNDC();
  t->SetTextAlign(22);
  t->DrawLatex(0.515,0.95,Form("#scale[1.0]{%s %s}",titles,titles2));

  canv->Print(Form("%s%s%s.png",path.c_str(),plots,cut));
  canv->Print(Form("%s%s%s.pdf",path.c_str(),plots,cut));
  
} //end of plot histograms

#endif

