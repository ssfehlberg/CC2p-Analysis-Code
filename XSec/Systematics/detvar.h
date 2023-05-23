#include "tools/paul_tol_colors.hpp"
#include "tools/xsec.h"
#include "tools/histograms.h"
using namespace Histograms;

class DetVar{

 public:

  virtual void Grab_Files();
  virtual void Grab_Iterations();
  virtual void area_normalize(TH1D* h_CV, TH1D* hist);
  virtual void set_Total(std::vector<TH1D*> h, TH1D* h_total_up);
  virtual void draw_plot(std::vector<TH1D*> h_vector ,const char* title, const char* name);
  virtual void main();
  
  //XSEc class
  //////////////
  xsec Xsec;
  
  //Detector Variation Samples
  /////////////////////////////
  static const int number_of_files = 11;
  const char* samples[number_of_files] = {"detVar_CV","detVar_LY_Attenuation","detVar_LY_Down","detVar_LY_Rayleigh",
					  "detVar_ThetaXZ", "detVar_ThetaYZ","detVar_X","detVar_YZ",
					  "detVar_Recombination","detVar_SCE","detVar_Overlay"};

  const char* sample_titles[number_of_files] = {"Central Value","LY Attenuation","LY Down","LY Rayleigh",
					        "ThetaXZ","ThetaYZ","X","YZ",
						"Recombination","SCE","Overlay"};

  //Individual Particles Plots
  /////////////////////////////

  //Muon
  TH2D* h_muon_matrices[number_of_files][num_var];
  TH1D* h_muon_num[number_of_files][num_var];
  TH1D* h_muon_denom[number_of_files][num_var];
  TH1D* h_overlay_muon_total[number_of_files][num_var];
  TH1D* h_overlay_muon_cc2p[number_of_files][num_var];
  TH1D* h_muon_xsec[number_of_files][num_var];

  //Lead
  TH2D* h_leading_matrices[number_of_files][num_var];
  TH1D* h_leading_num[number_of_files][num_var];
  TH1D* h_leading_denom[number_of_files][num_var];
  TH1D* h_overlay_leading_total[number_of_files][num_var];
  TH1D* h_overlay_leading_cc2p[number_of_files][num_var];
  TH1D* h_leading_xsec[number_of_files][num_var];

  //Recoil
  TH2D* h_recoil_matrices[number_of_files][num_var];
  TH1D* h_recoil_num[number_of_files][num_var];
  TH1D* h_recoil_denom[number_of_files][num_var];
  TH1D* h_overlay_recoil_total[number_of_files][num_var];
  TH1D* h_overlay_recoil_cc2p[number_of_files][num_var];
  TH1D* h_recoil_xsec[number_of_files][num_var];

  //Other Variables
  /////////////////
  TH2D* h_other_matrices[number_of_files][num_other_var];
  TH1D* h_other_num[number_of_files][num_other_var];
  TH1D* h_other_denom[number_of_files][num_other_var];
  TH1D* h_overlay_other_total[number_of_files][num_other_var];
  TH1D* h_overlay_other_cc2p[number_of_files][num_other_var];
  TH1D* h_other_xsec[number_of_files][num_other_var];

  //Empty files to grab everything
  ///////////////////////////////
  TFile* file;
  TFile* file_mc;
  
};

//Graab the files from the relevant samples
///////////////////////////////////////////
void DetVar::Grab_Files(){

  //Grab the Files
  ///////////////////////////

  //Grab the BNB, EXT, and,DIRT
  TFile* f_bnb = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run3/histograms_pelee_xsec_bnb.root"));
  TFile* f_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run3/histograms_pelee_xsec_ext.root "));
  TFile* f_dirt = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run3/histograms_pelee_xsec_dirt_wgt.root");
  
  for(int i=0; i < num_var; i++){
    h_bnb_muon[i] = (TH1D*)f_bnb->Get(Form("h_muon%s_bnb",var0[i])); //bnb
    h_ext_muon[i] = (TH1D*)f_ext->Get(Form("h_muon%s_ext",var0[i])); //ext
    h_dirt_muon[i] = (TH1D*)f_dirt->Get(Form("h_muon%s_dirt_wgt",var0[i])); //dirt

    h_bnb_leading[i] = (TH1D*)f_bnb->Get(Form("h_leading%s_bnb",var0[i])); //bnb
    h_ext_leading[i] = (TH1D*)f_ext->Get(Form("h_leading%s_ext",var0[i])); //ext
    h_dirt_leading[i] = (TH1D*)f_dirt->Get(Form("h_leading%s_dirt_wgt",var0[i])); //dirt

    h_bnb_recoil[i] = (TH1D*)f_bnb->Get(Form("h_recoil%s_bnb",var0[i])); //bnb
    h_ext_recoil[i] = (TH1D*)f_ext->Get(Form("h_recoil%s_ext",var0[i])); //ext
    h_dirt_recoil[i] = (TH1D*)f_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[i])); //dirt
  }

  for(int i=0; i < num_other_var; i++){
    h_bnb_other[i] = (TH1D*)f_bnb->Get(Form("h%s_bnb",other_var[i])); //bnb
    h_ext_other[i] = (TH1D*)f_ext->Get(Form("h%s_ext",other_var[i])); //ext
    h_dirt_other[i] = (TH1D*)f_dirt->Get(Form("h%s_dirt_wgt",other_var[i])); //dirt
  }

  for(int i=0; i < number_of_files; i++){
    file = new TFile(Form("root_files/detVar/histograms_pelee_xsec_%s.root",samples[i]));
    file_mc =  new TFile(Form("root_files/detVar/histograms_%s_mc_eff.root",samples[i]));
    
    for(int j=0; j < num_var; j++){
      h_overlay_muon_total[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_0",var0[j]));
      h_overlay_muon_cc2p[i][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_0",var0[j]));
      h_muon_matrices[i][j] = (TH2D*)file_mc->Get(Form("h_particle_matrices_muon_all%s_0",var[j]));
      h_muon_num[i][j] = (TH1D*)file_mc->Get(Form("h_particle_num_muon_all%s_0",var[j]));
      h_muon_denom[i][j] = (TH1D*)file_mc->Get(Form("h_particle_denom_muon_all%s_0",var[j]));

      h_overlay_leading_total[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_0",var0[j]));
      h_overlay_leading_cc2p[i][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_0",var0[j]));
      h_leading_matrices[i][j] = (TH2D*)file_mc->Get(Form("h_particle_matrices_lead_proton%s_0",var[j]));
      h_leading_num[i][j] = (TH1D*)file_mc->Get(Form("h_particle_num_lead_proton%s_0",var[j]));
      h_leading_denom[i][j] = (TH1D*)file_mc->Get(Form("h_particle_denom_lead_proton%s_0",var[j]));

      h_overlay_recoil_total[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_0",var0[j]));
      h_overlay_recoil_cc2p[i][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_0",var0[j]));
      h_recoil_matrices[i][j] = (TH2D*)file_mc->Get(Form("h_particle_matrices_recoil_proton%s_0",var[j]));
      h_recoil_num[i][j] = (TH1D*)file_mc->Get(Form("h_particle_num_recoil_proton%s_0",var[j]));
      h_recoil_denom[i][j] = (TH1D*)file_mc->Get(Form("h_particle_denom_recoil_proton%s_0",var[j]));
      
    } //end of loop over particle variables

    for(int j=0; j < num_other_var; j++){

      h_overlay_other_total[i][j] = (TH1D*)file->Get(Form("h%s_total_0",other_var[j]));
      h_overlay_other_cc2p[i][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_0",other_var[j]));
      h_other_matrices[i][j] = (TH2D*)file_mc->Get(Form("h_other_matrices%s_0",other_var[j]));
      h_other_num[i][j] = (TH1D*)file_mc->Get(Form("h_other_eff_num%s_0",other_var[j]));
      h_other_denom[i][j] = (TH1D*)file_mc->Get(Form("h_other_eff_denom%s_0",other_var[j]));

    } //end of loop over other variables

  }//end of loop over number of files


}

void DetVar::Grab_Iterations(){

  ifstream file;
  file.open("../num_iterations.csv");

  vector<vector<string>> content;
  vector<string> row;
  string line, word;

  while(getline(file, line)){
    row.clear();
    stringstream str(line);
    
    while(getline(str, word, ',')){
      row.push_back(word);
    }
    content.push_back(row);
  }

  for(int i=0; i < content.size();i++){
    for(int j=0;j<content[i].size();j++){
      int value = std::stoi(content[i][j]);
      
      if(i == 0){
	muon_iter[j] = value;
      }else if(i == 1){
	leading_iter[j] = value;
      } else if(i == 2){
	recoil_iter[j] = value;
      } else if(i == 3){
	other_iter[j] = value;
      }
      
    }
  }

}//end of grab iterations

//Area normalize the histograms
///////////////////////////////
void DetVar::area_normalize(TH1D* h_CV, TH1D* hist){
  hist->Add(hist,h_CV,-1);
  hist->Divide(hist,h_CV,1,1,"cp");
}

//Gets the total uncertainty
//I know it is not really well written
///////////////////////////////////////////
void DetVar::set_Total(std::vector<TH1D*> h, TH1D* h_total_up){

  double n_bins = h[0]->GetNbinsX();

  for(int i = 1; i < n_bins+1; i++){

    double value_LY_Atten = std::pow(h[1]->GetBinContent(i),2);
    double value_LY_Down = std::pow(h[2]->GetBinContent(i),2);
    double value_LY_Raliegh = std::pow(h[3]->GetBinContent(i),2);
    double value_Recombination = std::pow(h[4]->GetBinContent(i),2);
    double value_SCE = std::pow(h[5]->GetBinContent(i),2);
    double value_ThetaXZ = std::pow(h[6]->GetBinContent(i),2);
    double value_ThetaYZ = std::pow(h[7]->GetBinContent(i),2);
    double value_X = std::pow(h[8]->GetBinContent(i),2);
    double value_YZ = std::pow(h[9]->GetBinContent(i),2);

    Double_t value = std::sqrt(value_LY_Atten      + value_LY_Down + value_LY_Raliegh + 
			     value_Recombination + value_SCE     + value_ThetaXZ +
			     value_ThetaYZ       + value_X       + value_YZ);


    h_total_up->SetBinContent(i,value);

  }  
} //end of set_total

void DetVar::draw_plot(std::vector<TH1D*> h_vector ,const char* title, const char* name){

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
      
  h_vector[0]->SetYTitle("Differential Cross-Section [10^{-38} cm^{2} / Argon]");
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
  canv_xsec->Print(Form("images/detVar/xsec%s.png",name));

  //Draw plot with the fractional uncertainty
  ///////////////////////////////////////////
  TH1D* h_xsec_CV = (TH1D*)h_vector[0]->Clone();

  for(int i=0; i < h_vector.size(); i++){ //should be h_vector.size()-1. don't want to consider the overlay sample
    //When we look at the recombination or the space chaarge, we have to compare to the Overlay sample
    if(i == 8 || i == 9 || i == 10){
     area_normalize(h_vector[10],h_vector[i]);  
    }else{
      area_normalize(h_xsec_CV,h_vector[i]);
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

void DetVar::main(){

  //Stuff for Drawing
  //////////////////
  tolcols::init();
  Color_t colors[] = {0,kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet-8, kMagenta, kGray+2, kOrange+7};
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3g");     

  //Grab number of iterations
  ////////////////////////////
  Grab_Iterations();
  
  //Grab the files
  ///////////////
  Grab_Files();
  
  //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("root_files/detVar/systematics.root"),"RECREATE"); //output root file     
  
  ///////////////////////////////
  //Now to draw all of this garbage
  ////////////////////////////////

  //Particles
  ///////////////
  for(int j=0; j < num_var; j++){

    //Muon
    std::vector<TH1D*> muon_vec;
    for(int i = 0; i < number_of_files; i++){
      h_muon_xsec[i][j] = Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_overlay_muon_total[i][j],h_overlay_muon_cc2p[i][j],h_bnb_muon[j],
					     h_muon_matrices[i][j],muon_iter[j],h_muon_num[i][j], h_muon_denom[i][j],
					     muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_%s",var[j],samples[i]),"detVar",false,true);
      muon_vec.push_back(h_muon_xsec[i][j]);
    }

    std::cout<<"Right before draw plot for muon"<<std::endl;
    draw_plot(muon_vec , Form("True Muon %s",var_titles[j]), Form("_muon%s",var0[j]));
    std::cout<<"After draw plot for muon"<<std::endl;
    muon_vec.clear();

    //Leading Proton
    std::vector<TH1D*> leading_vec;
    for(int i = 0; i < number_of_files; i++){
      h_leading_xsec[i][j] = Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_overlay_leading_total[i][j],h_overlay_leading_cc2p[i][j],h_bnb_leading[j],
						h_leading_matrices[i][j],leading_iter[j],h_leading_num[i][j], h_leading_denom[i][j],
						leading_xsec_max[j],Form("True Leading Proton %s",var_titles[j]),Form("_leading%s_%s",var[j],samples[i]),"detVar",false,true);
      leading_vec.push_back(h_leading_xsec[i][j]);
    }

    draw_plot(leading_vec , Form("True Leading %s",var_titles[j]), Form("_leading%s",var0[j]));
    leading_vec.clear();

    //Recoil Proton
    std::vector<TH1D*> recoil_vec;
    for(int i = 0; i < number_of_files; i++){
      h_recoil_xsec[i][j] = Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_overlay_recoil_total[i][j],h_overlay_recoil_cc2p[i][j],h_bnb_recoil[j],
					     h_recoil_matrices[i][j],recoil_iter[j],h_recoil_num[i][j], h_recoil_denom[i][j],
					     recoil_xsec_max[j],Form("True Recoil Proton %s",var_titles[j]),Form("_recoil%s_%s",var[j],samples[i]),"detVar",false,true);
      recoil_vec.push_back(h_recoil_xsec[i][j]);
    }

    draw_plot(recoil_vec , Form("True Recoil %s",var_titles[j]), Form("_recoil%s",var0[j]));
    recoil_vec.clear();
      
  } //end of loop over number of variaables


  //Other variables
  /////////////////
  for(int j=0; j < num_other_var; j++){
    
    std::vector<TH1D*> other_vec;
    for(int i = 0; i < number_of_files; i++){
      h_other_xsec[i][j] = Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_overlay_other_total[i][j],h_overlay_other_cc2p[i][j],h_bnb_other[j],
					      h_other_matrices[i][j],other_iter[j],h_other_num[i][j], h_other_denom[i][j],
					      other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("%s_%s",other_var[j],samples[i]),"detVar",false,true);
      other_vec.push_back(h_other_xsec[i][j]);
    }
    
    draw_plot(other_vec , Form("True %s",other_var_titles[j]), Form("%s",other_var[j]));
    other_vec.clear();
    
  } //end of loop over other variables

  tfile_mine->Close();

} //end of program
