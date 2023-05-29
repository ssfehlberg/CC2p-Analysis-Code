#include "shared.h"

class detvar_plot{

 public:
  virtual void main();

 private:
  virtual void Grab_Histograms();
  virtual void plot_total_error(const char* variable, bool flip_legend, std::vector<TH1D*> hist,const char* title, TH1D* h_total);


  //File Specifications
  /////////////////////////
  static const int num_files = 9;
  TFile* file[num_files];
  const char* file_names[num_files] = {"detvar_plot_LY_Attenuation","detvar_plot_LY_Down","detvar_plot_LY_Rayleigh",
				       "detvar_plot_ThetaXZ", "detvar_plot_ThetaYZ",
				       "detvar_plot_X","detvar_plot_YZ",
				       "detvar_plot_Recombination","detvar_plot_SCE"};

  const char* sample_titles[num_files] = {"LY Attenuation","LY Down","LY Rayleigh",
                                          "ThetaXZ","ThetaYZ",
                                          "X","YZ",
                                          "Recombination","SCE"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon_total[num_var];
  TH2D* h_muon_covariance[num_files][num_var];
  bool flip_legend_muon[num_var] = {true,false,false};
  
  //Leading
  TH1D* h_leading[num_files][num_var];
  TH1D* h_leading_total[num_var];
  TH2D* h_leading_covariance[num_files][num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH1D* h_recoil[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  TH2D* h_recoil_covariance[num_files][num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};
  
  /*static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
                                                 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};
  bool flip_legend_other[num_other_var] = {true,false,false,false,false,false,false};*/

  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TH2D* h_other_covariance[num_files][num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};
   
}; //end of class definition

//////////
//[main]
////////
void detvar_plot::main(){
  
  //Grab the histograms
  /////////////////////
  Grab_Histograms();

  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/detvar_plot/systematics.root"),"RECREATE"); //output root file 

  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){

    std::vector<TH1D*> muon;
    std::vector<TH2D*> muon_2D;
    std::vector<TH1D*> leading;
    std::vector<TH2D*> leading_2D;
    std::vector<TH1D*> recoil;
    std::vector<TH2D*> recoil_2D;

    for(int f=0; f < num_files; f++){
      muon.push_back(h_muon[f][j]);
      muon_2D.push_back(h_muon_covariance[f][j]);
      leading.push_back(h_leading[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      recoil.push_back(h_recoil[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]);
    }

    plot_total_error(Form("_muon%s",var0[j]),flip_legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    make_total_covariance_matrix(Form("_muon%s",var0[j]),Form("Total DetVar Error: Muon %s",var_titles[j]),muon_2D,"detvar_plot/Total");

    plot_total_error(Form("_leading%s",var0[j]),flip_legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    make_total_covariance_matrix(Form("_leading%s",var0[j]),Form("Total DetVar Error: Leading Proton %s",var_titles[j]),leading_2D,"detvar_plot/Total");
    
    plot_total_error(Form("_recoil%s",var0[j]),flip_legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    make_total_covariance_matrix(Form("_recoil%s",var0[j]),Form("Total DetVar Error: Recoil Proton %s",var_titles[j]),recoil_2D,"detvar_plot/Total");
    
  }
  
  for(int j=0; j < num_other_var; j++){

    std::vector<TH1D*> other;
    std::vector<TH2D*> other_2D;
    
    for(int f=0; f<num_files; f++){
      other.push_back(h_other[f][j]);
      other_2D.push_back(h_other_covariance[f][j]);
    }

    plot_total_error(Form("%s",other_var[j]),flip_legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
    make_total_covariance_matrix(Form("%s",other_var[j]),Form("Total DetVar Error: %s",other_var_titles[j]),other_2D,"detvar_plot/Total");
  }
  
  tfile_mine->Close();
    
}


////////////////////////
//[Grab_Histograms]
//Grabs all the histograms from the appropriate files
////////////////////////
void detvar_plot::Grab_Histograms(){
                                                                                                                                                                                                   
  for(int f=0; f < num_files; f++){
    file[f] =  new TFile(Form("../root_files/detvar_plot/systematics_%s.root",file_names[f]));
    for(int j=0; j < num_var; j++){
      h_muon[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_muon%s",var0[j]));
      h_muon_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancemuon_%s",var0[j]));

      h_leading[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_leading%s",var0[j]));
      h_leading_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covarianceleading_%s",var0[j]));

      h_recoil[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j]));
      h_recoil_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancerecoil_%s",var0[j]));
    }

    for(int j=0; j< num_other_var; j++){
      h_other[f][j] =  (TH1D*)file[f]->Get(Form("hist_fractional_errors%s",other_var[j]));
      h_other_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariance%s",other_var[j]));
    }
  }

} //end of grab histograms

///////////////
//[plot_total_error]
//Plots the total error for the detector variation samples ONLY
// variables: which variable we are using
//flip_legend indicates wether legend is on left or right
//hist: vector of TH1D that are in the order of asdlkfj. NOTE: Assumes they are fractional uncertainty
//title: title of the plot
//h_total: will be the sum (in quadrature) of all the contributions in hist
//////////////
void detvar_plot::plot_total_error(const char* variable, bool flip_legend, std::vector<TH1D*> hist,const char* title, TH1D* h_total){

  //Stuff for Drawing
  //////////////////////
  Color_t colors[] = {kRed, kOrange+6, kYellow-3, kGreen+2, kCyan, kBlue, kViolet+1, kMagenta, kGray};
  int line_style[] = {1,9,1,9,1,9,1,9,1,9,1,9,1};
  static const int num_files = 9;
  const char* file_titles[num_files] = {"LY Attenuation","LY Down","LY Rayleigh",
					"ThetaXZ","ThetaYZ",
					"X","YZ",
					"Recombination","SCE"};
  //Now to plot
  //////////////
  TCanvas* canv = new TCanvas(Form("canv%s",variable),Form("canv%s",variable),2000,1500);
  TLegend* legend;
  if(flip_legend == true){
    legend = new TLegend(0.115, 0.65, 0.525, 0.89);
  } else {
    legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  }
  legend->SetNColumns(2);
  
  hist[0]->Draw("HIST TEXT00");
  hist[0]->SetLineColor(colors[0]);
  hist[0]->SetLineWidth(3);
  hist[0]->SetLineStyle(line_style[0]);
  hist[0]->SetMarkerColor(colors[0]);  
  hist[0]->SetTitle(Form("%s",title)); //title
  hist[0]->SetXTitle(Form("%s",title)); //xtitle
  hist[0]->GetXaxis()->SetTitleSize(40);
  hist[0]->GetXaxis()->SetTitleFont(43);
  hist[0]->GetXaxis()->SetTitleOffset(1.5);
  hist[0]->GetXaxis()->SetLabelFont(43);
  hist[0]->GetXaxis()->SetLabelSize(30);
  hist[0]->SetYTitle("Fractional Uncertainty"); //Y title
  hist[0]->GetYaxis()->SetTitleSize(40);
  hist[0]->GetYaxis()->SetTitleFont(43);
  hist[0]->GetYaxis()->SetTitleOffset(1.5);
  hist[0]->GetYaxis()->SetLabelFont(43);
  hist[0]->GetYaxis()->SetLabelSize(30);
  hist[0]->SetMaximum(1.0); //max
  hist[0]->SetMinimum(0); //min  
  
  legend->AddEntry(hist[0],Form("%s",file_titles[0]),"lpf");
  for(int f=1; f < hist.size(); f++){
    hist[f]->Draw("HIST TEXT00 SAME");
    hist[f]->SetLineColor(colors[f]);
    hist[f]->SetLineWidth(3);
    hist[f]->SetLineStyle(line_style[f]);
    hist[f]->SetMarkerColor(colors[f]);
    legend->AddEntry(hist[f],Form("%s",file_titles[f]),"lpf");
  }
    
  h_total = total_error(Form("%s",variable),hist);
  h_total->Draw("HIST TEXT00 SAME");
  h_total->SetLineColor(kBlack);
  h_total->SetLineWidth(5);
  h_total->SetLineStyle(1);
  h_total->SetMarkerColor(kBlack);
  h_total->Write();
  
  legend->AddEntry(h_total,"Total detvar_plot Uncertainty","lpf");
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, In-Progress","C");
  legend->Draw("same");
  gStyle->SetOptStat(0);
  canv->Print(Form("../images/detVar/Total/_total_error%s.png",variable));
 
} //end of plot_total_error
  


  

