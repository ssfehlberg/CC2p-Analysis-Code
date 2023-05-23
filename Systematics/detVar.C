#define detVar_cxx
#include "../paul_tol_colors.hpp"

void area_normalize(TH1D* h_CV, TH1D* hist){
  hist->Add(hist,h_CV,-1);
  hist->Divide(hist,h_CV,1,1,"B");
}

void set_Total(std::vector<TH1D*> h, TH1D* h_total_up){

  double n_bins = h[0]->GetNbinsX();

  for(int i = 1; i < n_bins+1; i++){
    double value_LY_Atten = std::pow(h[1]->GetBinContent(i),2);
    double value_LY_Down = std::pow(h[2]->GetBinContent(i),2);
    double value_LY_Raliegh = std::pow(h[3]->GetBinContent(i),2);
    double value_Recombination = std::pow(h[4]->GetBinContent(i),2);
    double value_SCE = std::pow(h[5]->GetBinContent(i),2);
    double value_ThetaXZ = std::pow(h[6]->GetBinContent(i),2);
    double value_ThetaYZ = std::pow(h[7]->GetBinContent(i),2);
    double value_X = std::pow(h[9]->GetBinContent(i),2);
    double value_YZ = std::pow(h[9]->GetBinContent(i),2);

    double value = std::sqrt(value_LY_Atten      + value_LY_Down + value_LY_Raliegh + 
			     value_Recombination + value_SCE     + value_ThetaXZ +
			     value_ThetaYZ       + value_X       + value_YZ);

    h_total_up->SetBinContent(i,value);

  }  
} //end of set_total


void draw_plot(std::vector<TH1D*> h_vector ,const char* title, const char* name){

  Color_t colors[] = {kBlack,kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet-8, kMagenta, kOrange+7,kGray+2};
  static const int number_of_files = 11;
  const char* sample_titles[number_of_files] = {"Central Value","LY Attenuation",
						"LY Down","LY Rayleigh",
						"Recombination","SCE","ThetaXZ",
						"ThetaYZ","X","YZ","Overlay"};

  //Draw plot with just the cross-sections
  ////////////////////////////////////////
  TCanvas* canv_xsec = new TCanvas("canv_xsec","canv_xsec",2000,1500);
  canv_xsec->SetGridx(); 
  h_vector[0]->Draw("hist");
  h_vector[0]->SetLineColor(kBlack);
  
  h_vector[0]->SetLineWidth(3);
  h_vector[0]->SetTitle(title);
  h_vector[0]->SetXTitle(title);
  h_vector[0]->GetXaxis()->SetTitleSize(40);
  h_vector[0]->GetXaxis()->SetTitleFont(43);
  h_vector[0]->GetXaxis()->SetTitleOffset(1.5);
  h_vector[0]->GetXaxis()->SetLabelFont(43);
  h_vector[0]->GetXaxis()->SetLabelSize(30);
      
  h_vector[0]->SetYTitle("Number of Events");
  h_vector[0]->GetYaxis()->SetTitleSize(40);
  h_vector[0]->GetYaxis()->SetTitleFont(43);
  h_vector[0]->GetYaxis()->SetTitleOffset(1.5);
  h_vector[0]->GetYaxis()->SetLabelFont(43);
  h_vector[0]->GetYaxis()->SetLabelSize(30);
  
  h_vector[0]->SetTitle(title);
  
  TLegend* legend_xsec = new TLegend(0.48, 0.65, 0.89, 0.89);
  legend_xsec->SetNColumns(3);
  legend_xsec->AddEntry(h_vector[0],Form("%s",sample_titles[0]),"L");
  
  for(int i = 1; i < h_vector.size()-1; i++){
    h_vector[i]->Draw("hist same");
    h_vector[i]->SetLineColor(colors[i]);
    h_vector[i]->SetLineWidth(3);
    legend_xsec->AddEntry(h_vector[i],Form("%s",sample_titles[i]),"L");
  }
  
  legend_xsec->Draw("same");
  canv_xsec->Print(Form("images/detVar/num_events%s.png",name));

  //Draw plot with the fractional uncertainty
  ///////////////////////////////////////////
  TH1D* h_CV = (TH1D*)h_vector[0]->Clone();

  for(int i=0; i < h_vector.size(); i++){ //should be h_vector.size()-1. don't want to consider the overlay sample
    //When we look at the recombination or the space chaarge, we have to compare to the Overlay sample
    if(i == 8 || i == 9 || i == 10){
     area_normalize(h_vector[10],h_vector[i]);  
    }else{
      area_normalize(h_CV,h_vector[i]);
    }
  }
  
  TCanvas* canv = new TCanvas("canv","canv",2000,1500);
  canv->SetGridx(); 
  h_vector[0]->Draw("hist");
  h_vector[0]->SetLineColor(kBlack);
  
  h_vector[0]->SetLineWidth(3);
  h_vector[0]->SetTitle(title);
  h_vector[0]->SetXTitle(title);
  h_vector[0]->GetXaxis()->SetTitleSize(40);
  h_vector[0]->GetXaxis()->SetTitleFont(43);
  h_vector[0]->GetXaxis()->SetTitleOffset(1.5);
  h_vector[0]->GetXaxis()->SetLabelFont(43);
  h_vector[0]->GetXaxis()->SetLabelSize(30);
      
  h_vector[0]->SetYTitle("Fractional Uncertainty");
  h_vector[0]->GetYaxis()->SetTitleSize(40);
  h_vector[0]->GetYaxis()->SetTitleFont(43);
  h_vector[0]->GetYaxis()->SetTitleOffset(1.5);
  h_vector[0]->GetYaxis()->SetLabelFont(43);
  h_vector[0]->GetYaxis()->SetLabelSize(30);
  
  h_vector[0]->SetTitle(title);
  h_vector[0]->SetMaximum(1.0);
  h_vector[0]->SetMinimum(-0.6);
  
  TLegend* legend = new TLegend(0.45, 0.62, 0.89, 0.89);
  legend->SetNColumns(3);
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C");
  legend->AddEntry(h_vector[0],Form("%s",sample_titles[0]),"L");
  
  for(int i = 1; i < h_vector.size() -1 ; i++){
    h_vector[i]->Draw("hist same");
    h_vector[i]->SetLineColor(colors[i]);
    h_vector[i]->SetLineWidth(3);
    legend->AddEntry(h_vector[i],Form("%s",sample_titles[i]),"L");
  }
  
  TH1D* h_total_up = (TH1D*)h_vector[0]->Clone();
  set_Total(h_vector, h_total_up);
  h_total_up->Draw("hist same text");
  h_total_up->SetLineColor(kBlack);
  h_total_up->SetLineStyle(10);
  h_total_up->Write(Form("hist_fractional_errors%s",name));
  
  TH1D* h_total_down = (TH1D*)h_total_up->Clone();
  h_total_down->Scale(-1.0);
  h_total_down->Draw("hist same");
  h_total_down->SetLineColor(kBlack);
  h_total_down->SetLineStyle(4);
  
  legend->AddEntry(h_total_up,"Total Uncertainty (+)","lf");
  legend->AddEntry(h_total_down,"Total Uncertainty (-)", "lf");
  legend->Draw("same");
  
  canv->Print(Form("images/detVar/%s.png",name));

}



void detVar(){

  //Detector Variation Samples
  /////////////////////////////
  static const int number_of_files = 11;
  const char* samples[number_of_files] = {"detVar_CV",
					  "detVar_LY_Attenuation","detVar_LY_Down","detVar_LY_Rayleigh",
					  "detVar_Recombination","detVar_SCE","detVar_ThetaXZ", "detVar_ThetaYZ","detVar_X","detVar_YZ","detVar_Overlay"};
  const char* sample_titles[number_of_files] = {"Central Value",
						"LY Attenuation","LY Down","LY Rayleigh",
						"Recombination","SCE","ThetaXZ","ThetaYZ","X","YZ", "CV Overlay"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};
  
  //Muon
  TH1D* h_muon[number_of_files][num_var];
  TH1D* h_muon0[num_var];
  TH1D* h_muon_total_up[num_var];
  TH1D* h_muon_total_down[num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];
 
  //Lead
  TH1D* h_leading[number_of_files][num_var];
  TH1D* h_leading0[num_var];
  TH1D* h_leading_total_up[num_var];
  TH1D* h_leading_total_down[num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  
  //Recoil
  TH1D* h_recoil[number_of_files][num_var];
  TH1D* h_recoil0[num_var];
  TH1D* h_recoil_total_up[num_var];
  TH1D* h_recoil_total_down[num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  
  //Other Variables
  /////////////////////////
  /*  static const int num_other_var = 9;
  const	char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading",
                                          "_opening_angle_mu_both","_delta_PT","_delta_alphaT",
					  "_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})","cos(#gamma_{#mu,p_{L}})",
						 "cos(#gamma_{#mu,p_{L} + p_{R}})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)",
						 "#delta #phi_{T} (Deg.)","Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"}; 
  */
  
  static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
					"_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})""cos(#gamma_{#mu,p_{L} + p_{R}})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)",
                                                 "#delta #phi_{T} (Deg.)","Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};


  TH1D* h_other[number_of_files][num_other_var];
  TH1D* h_other0[num_other_var];
  TH1D* h_other_total_up[num_other_var];
  TH1D* h_other_total_down[num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];

  
  //Empty file to grab everything
  TFile* file;
  
  //Grab the Files
  ///////////////////////////
  for(int i=0; i < number_of_files; i++){
    file = new TFile(Form("root_files/detVar/histograms_pelee_xsec_%s.root",samples[i]));

    for(int j=0; j < num_var; j++){
      h_muon[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_0",var0[j]));
      h_leading[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_0",var0[j]));
      h_recoil[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_0",var0[j]));
    } //end of loop over particle variables

    for(int j=0; j < num_other_var; j++){
      h_other[i][j] = (TH1D*)file->Get(Form("h%s_total_0",other_var[j]));
    } //end of loop over other variables

  }//end of loop over number of files

  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("root_files/detVar/systematics.root"),"RECREATE"); //output root file     
  

  ///////////////////////////////
  //Now to draw all of this garbage
  ////////////////////////////////
  

  for(int j=0; j < num_var; j++){

    //Muon
    std::vector<TH1D*> muon_vec;
    for(int k=0; k < number_of_files; k++){
      muon_vec.push_back(h_muon[k][j]);
   }
    draw_plot(muon_vec , Form("Muon %s",var_titles[j]), Form("_muon%s",var0[j]));

    //Leading
    std::vector<TH1D*> leading_vec;
    for(int k=0; k < number_of_files; k++){
      leading_vec.push_back(h_leading[k][j]);
    }
    draw_plot(leading_vec , Form("Leading Proton %s",var_titles[j]), Form("_leading%s",var0[j]));

    //Recoil
    std::vector<TH1D*> recoil_vec;
    for(int k=0; k < number_of_files; k++){
      recoil_vec.push_back(h_recoil[k][j]);
   }
    draw_plot(recoil_vec , Form("Recoil Proton %s",var_titles[j]), Form("_recoil%s",var0[j]));
  }
    
  
  //Now for the Other Variables
  ///////////////////////////
  for(int j=0; j < num_other_var; j++){
    std::vector<TH1D*> other_vec;
    for(int k=0; k < number_of_files; k++){
      other_vec.push_back(h_other[k][j]);
   }
    draw_plot(other_vec, Form("%s",other_var_titles[j]), Form("%s",other_var[j]));
   
  }
  
  tfile_mine->Close();

} //end of program
