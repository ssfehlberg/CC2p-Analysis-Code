#include "constants.h"
using namespace Constants;

class xsec{

public:
  virtual TH1D* make_efficiency_plot(TH1D* h_num_input, TH1D* h_denom_input,const char* title, const char* name,const char* directory);
  //virtual void plot_matrices(TH2D* h_matrix,const char* title, const char* name);
  virtual TH1D* make_bnb_plot(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input, const char* name,const char* directory);
  virtual TH1D* cross_section(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			     TH2D* h_smearing, int iter, TH1D* h_eff_input, TH1D* h_denom_input,
			      double maximum, const char* title, const char* name,const char* directory,bool flip_legend,bool print_contents = true);
};

//Makes the efficiency plot from the numerator and denominator
TH1D* xsec::make_efficiency_plot(TH1D* h_num_input, TH1D* h_denom_input,const char* title, const char* name,const char* directory){

  TCanvas* canv_eff = new TCanvas("canv_eff","canv_eff",2000,1500);
  TH1D* h_num = (TH1D*)h_num_input->Clone();
  TH1D*	h_denom = (TH1D*)h_denom_input->Clone();
  h_num->Divide(h_num,h_denom,1.0,1.0, "cp");
  h_num->Draw("1e1p");
  h_num->SetTitle(Form("%s ; %s ; Efficiency",title,title));
  h_num->SetLineColor(kViolet);
  h_num->SetMaximum(0.4);
  h_num->SetMinimum(0);

  /*TLatex* t;
  t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",title));
  t->DrawLatex(0.3,0.92,Form("%s",pot_num));
  t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
  */
  
  //a[i] = new TLine(xlim_eff[i],0,xlim_eff[i],1);
  //a[i]->Draw("same");
  //a[i]->SetLineColor(kBlack);
  //a[i]->SetLineWidth(4);
  //a1[i] = new TLine(xlim_high_eff[i],0,xlim_high_eff[i],1);
  //a1[i]->Draw("same");
  //a1[i]->SetLineColor(kBlack);
  //a1[i]->SetLineWidth(4);

  canv_eff->Print(Form("images/%s/Prep/%s_eff.png",directory,name)); 

  return h_num;
  
}

/*
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
*/

TH1D* xsec::make_bnb_plot(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input, const char* name,const char* directory){

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
  c_background->Print(Form("images/%s/BNB_Checks/ext_after_subtraction%s.png",directory,name));

  //Draw the subtracted backgrounds as a check
  TCanvas* canv_bnb = new TCanvas("canv_bnb","BNB with Background Subtracted",800,600); //you have to define a canvas to draw.
  TH1D* h_bnb = (TH1D*)h_bnb_input->Clone();
  std::cout<<"Num entries in h_bnb: "<<h_bnb->Integral()<<std::endl;
  h_bnb->Add(h_ext,-1);
  std::cout<<"Num entries in h_bnb After Subtraction: "<<h_bnb->Integral()<<std::endl;
  h_bnb->Draw("HIST");
  canv_bnb->Print(Form("images/%s/BNB_Checks/bnb_after_everything_subtracted%s.png",directory,name));

  return h_bnb;
      
}

TH1D* xsec::cross_section(TH1D* h_ext_input, TH1D* h_dirt_input,TH1D* h_overlay_total_input, TH1D* h_overlay_cc2p_input, TH1D* h_bnb_input,
			 TH2D* h_smearing, int iter, TH1D* h_num_input, TH1D* h_denom_input,
			  double maximum, const char* title, const char* name,const char* directory,bool flip_legend,bool print_contents = true){

  //Make the Efficiency Plot
  /////////////////////////
  TH1D* h_eff_input = make_efficiency_plot(h_num_input, h_denom_input,title,name,directory);
  
  //Make the 2D Plot
  ////////////////
  //plot_matrices(TH2D* h_smearing,title,name)
    
  //Here is where the magic happens!
  //////////////////////////////////
  TH1D* h_bnb = make_bnb_plot(h_ext_input,h_dirt_input,h_overlay_total_input,h_overlay_cc2p_input,h_bnb_input,name,directory);
  RooUnfoldResponse Response_Matrix(0,0,h_smearing,"Response_Matrix");
  RooUnfoldBayes RooUnfoldBayes_Data(&Response_Matrix, h_bnb, iter);
  TH1D* h_unfolded_signal =(TH1D* )RooUnfoldBayes_Data.Hreco();

  //One without the constants applied
  /////////////////////////////////////////
  TCanvas* c_xsec_no_const = new TCanvas("c_xsec_no_const","Unfolded Data ONLY",800,600);
  TH1D* h_xsec_no_const = (TH1D*)h_unfolded_signal->Clone();
  h_xsec_no_const->SetTitle(Form("%s: Unfolded Only",title));
  h_xsec_no_const->SetXTitle(Form("%s: Unfolded Only",title));
  h_xsec_no_const->SetYTitle("Arb. Units");
  h_xsec_no_const->Draw("1e1p");
  c_xsec_no_const->Print(Form("images/%s/XSec/xsec_w_o_constants%s.png",directory,name));

  //One with the efficiency applied
  /////////////////////////////////
  TCanvas* c_xsec_eff = new TCanvas("c_xsec_eff","Unfolded Data w/ Efficiency",800,600);
  TH1D* h_xsec_eff = (TH1D*)h_unfolded_signal->Clone();
  TH1D* h_eff = (TH1D*) h_eff_input->Clone(); 
  h_xsec_eff->Divide(h_unfolded_signal,h_eff,1,1,"cp");
  h_xsec_eff->SetTitle(Form("%s: Efficiency Applied",title));
  h_xsec_eff->SetXTitle(Form("%s: Efficiency Applied",title));
  h_xsec_eff->SetYTitle("Arb. Units");
  h_xsec_eff->Draw("1e1p");
  c_xsec_eff->Print(Form("images/%s/XSec/xsec_w_eff%s.png",directory,name));

  //One with everything applied
  /////////////////////////////
  TCanvas* c_xsec = new TCanvas("c_xsec","Unfolded Data w Everything",2000,1500);
  c_xsec->SetGridx();
  TH1D* h_xsec = (TH1D*)h_unfolded_signal->Clone();
  h_xsec->Draw("1e1p");

  //drawing denom on the same plot
  TH1D* h_denom = (TH1D*) h_denom_input->Clone();
  h_denom->Draw("HIST SAME") ;
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
  double scale_value = 1/(N_targets*flux_value);
  h_xsec->Scale(scale_value);
  h_denom->Scale(scale_value);

  //Setting all the values now that we are done with all the manipulation
  //////////////////////////////////////////////////////////////////
  h_xsec->SetLineColor(kBlack);
  h_xsec->SetLineWidth(4);
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
    
  //legend stuff
  TLegend* legend_xsec;
  if(flip_legend){
    legend_xsec = new TLegend(0.105, 0.5, 0.584, 0.89);
  }else{
    legend_xsec = new TLegend(0.37, 0.5, 0.859, 0.89);
  }
  legend_xsec->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C"); // option "C" allows to center the header 
  legend_xsec->AddEntry(h_xsec,Form("%d Iterations Unfolded Data (Stat. Error)",iter),"lepf");
  legend_xsec->AddEntry(h_denom,"uB Tune v2","L");
  legend_xsec->SetLineWidth(0);
  legend_xsec->SetFillColor(kWhite);
  legend_xsec->SetTextSize(0.03);
  legend_xsec->Draw("SAME");
  c_xsec->Print(Form("images/%s/XSec/xsec%s.png",directory,name));

  return h_xsec;
 
}
