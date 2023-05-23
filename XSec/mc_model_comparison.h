#include "constants.h"
using namespace Constants;

class MC_Comparison{

public:
  virtual void area_normalize(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF);
  virtual void cross_section(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF);
  virtual void make_plot(TH1D* h_empirical0,TH1D* h_nieves0,TH1D* h_susa0,TH1D* h_GCF0, double max_value, bool yes_fsi, const char* title,const char* name);
  virtual void make_plot_weird(TH1D* h_empirical0,TH1D* h_empirical0_true, TH1D* h_nieves0,TH1D* h_nieves0_true,TH1D* h_susa0,TH1D* h_susa0_true, TH1D* h_GCF0,TH1D* h_GCF0_true, double max_value, const char* title,const char* name);
  virtual void model_comparison(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF, double max_value, bool yes_fsi, const char* title,const char* name);
 
};

void MC_Comparison::area_normalize(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF){

  double n_bins = h_empirical->GetNbinsX();
  double empirical_integral = h_empirical->Integral();
  double nieves_integral = h_nieves->Integral();
  double susa_integral = h_susa->Integral();
  double GCF_integral = h_GCF->Integral();

  for(int i=1; i < n_bins+1; i++){
    double empirical_bin_content = h_empirical->GetBinContent(i);
    double empirical_value = empirical_bin_content/empirical_integral;
    h_empirical->SetBinContent(i,empirical_value);
      
    double nieves_bin_content = h_nieves->GetBinContent(i);
    double nieves_value = nieves_bin_content/nieves_integral;
    h_nieves->SetBinContent(i,nieves_value);
    
    double susa_bin_content = h_susa->GetBinContent(i);
    double susa_value = susa_bin_content/susa_integral;
    h_susa->SetBinContent(i,susa_value);
    
    double GCF_bin_content = h_GCF->GetBinContent(i);
    double GCF_value = GCF_bin_content/GCF_integral;
    h_GCF->SetBinContent(i,GCF_value);
  }


}//end of area_normalize


void MC_Comparison::cross_section(TH1D* h_empirical,TH1D* h_nieves,TH1D* h_susa,TH1D* h_GCF){

  //here are the sigmas. Taken from the GENIE splines
  double sigma_empirical = 3.02249 * 1E-37;
  double sigma_nieves = 27.1682 * 1E-38;
  double sigma_susa = 38.1417 * 1E-38;

  //N is the total number of events generated in each MEC type
  double N_empirical = 4000000;
  double N_nieves = 2400000;
  double N_susa = 3800000;
  double N_GCF = 500000;

  //Here is were we calculate the dsigma/dx and the SD of dsigma/dx
  double n_bins = h_empirical->GetNbinsX();

  for(int i=1; i < n_bins+1; i++){
    
    double delta_x = h_empirical->GetBinWidth(i);
    
    //empirical
    double n_empirical = h_empirical->GetBinContent(i) * N_empirical;
    double value_empirical = (sigma_empirical * n_empirical)/( N_empirical * delta_x);
    double SD_empirical = (sigma_empirical)/(delta_x*N_empirical) * std::sqrt(((N_empirical - n_empirical)*n_empirical)/(N_empirical));
    h_empirical->SetBinContent(i,value_empirical);
    h_empirical->SetBinError(i,SD_empirical);
    std::cout<<"EMPIRICAL value: "<<value_empirical<<std::endl;
	  
    //nieves
    double n_nieves = h_nieves->GetBinContent(i);
    double value_nieves = (sigma_nieves * n_nieves)/(N_nieves * delta_x);
    double SD_nieves = (sigma_nieves)/(delta_x*N_nieves) * std::sqrt(((N_nieves - n_nieves)*n_nieves)/(N_nieves));
    h_nieves->SetBinContent(i,value_nieves);
    h_nieves->SetBinError(i,SD_nieves);
    std::cout<<"NIEVES value: "<<value_nieves<<std::endl;
    
    //susa
    double n_susa = h_susa->GetBinContent(i);
    double value_susa = (sigma_susa * n_susa)/(N_susa * delta_x);
    double SD_susa = (sigma_susa)/(delta_x*N_susa) * std::sqrt(((N_susa - n_susa)*n_susa)/(N_susa));
    h_susa->SetBinContent(i,value_susa);
    h_susa->SetBinError(i,SD_susa);
    std::cout<<"SUSA value: "<<value_susa<<std::endl;
                                                                           
    //Dealing with the GCF is a bit trickier
    //Each event in the GCF is assigned a weight which modifies the CCQE differential cross section
    //While we hav
    double n_GCF = h_GCF->GetBinContent(i);
    double value_GCF = (n_GCF)/(N_GCF*delta_x);
    double SD_GCF = (1)/(delta_x*N_GCF) * std::sqrt(((N_GCF - n_GCF)*n_GCF)/(N_GCF));
    h_GCF->SetBinContent(i,value_GCF);
    h_GCF->SetBinError(i,SD_GCF);
    std::cout<<"GCF value: "<<value_GCF<<std::endl;

  } //end of for loop

} //end of cross-section

void MC_Comparison::make_plot(TH1D* h_empirical0,TH1D* h_nieves0,TH1D* h_susa0,TH1D* h_GCF0, double max_value, bool yes_fsi, const char* title,const char* name){

  TH1D* h_empirical = (TH1D*)h_empirical0->Clone();
  TH1D* h_nieves = (TH1D*)h_nieves0->Clone();
  TH1D* h_susa = (TH1D*)h_susa0->Clone();
  TH1D* h_GCF = (TH1D*)h_GCF0->Clone();

  //the canvas
  TCanvas* canv_theory = new TCanvas("canv_theory","canv_theory",1000,700);

  h_empirical->Draw("hist");
  h_empirical->SetLineColor(color2);
  h_empirical->SetTitle(Form("%s",title));
  h_empirical->SetXTitle(Form("%s",title));
  h_empirical->GetXaxis()->SetTitleSize(28);
  h_empirical->GetXaxis()->SetTitleFont(43);
  h_empirical->SetYTitle("Fractional Number of Events");
  h_empirical->GetYaxis()->SetTitleSize(28);
  h_empirical->GetYaxis()->SetTitleFont(43);
  h_empirical->SetMaximum(max_value);

  h_nieves->Draw("hist same");
  h_nieves->SetLineColor(color3);

  h_susa->Draw("hist same");
  h_susa->SetLineColor(color4);

  //h_GCF->Draw("hist same");
  //h_GCF->SetLineColor(color5);

  TLegend* legend_theory = new TLegend(0.14, 0.72, 0.87, 0.87);

  if(yes_fsi == false){
    legend_theory->AddEntry(h_empirical,"Empirical MEC + Lwellyn Smith QE","L");
    legend_theory->AddEntry(h_nieves,"Nieves MEC + Nieves QE","L");
    legend_theory->AddEntry(h_susa,"SuSAv2 MEC + SuSAv2 QE","L");
    //legend_theory->AddEntry(h_GCF,"GCF","L");
  }else{
    legend_theory->AddEntry(h_empirical,"Empirical MEC + Lwellyn Smith QE + GENIE hA2018 FSI","L");
    legend_theory->AddEntry(h_nieves,"Nieves MEC + Nieves QE + GENIE hA2018 FSI","L");
    legend_theory->AddEntry(h_susa,"SuSAv2 MEC + SuSAv2 QE + GENIE hN2018 FSI","L");
    //legend_theory->AddEntry(h_GCF,"GCF + GENIE hA2018 FSI","L") ;
  }
  
  legend_theory->SetLineWidth(0);
  legend_theory->SetFillColor(kWhite);
  legend_theory->SetTextSize(0.03);
  legend_theory->Draw("SAME");
  canv_theory->Print(Form("images/Model_Comparisons/%s.png",name));

  
} //end of make plot


void MC_Comparison::make_plot_weird(TH1D* h_empirical0,TH1D* h_empirical0_true, TH1D* h_nieves0,TH1D* h_nieves0_true,TH1D* h_susa0,TH1D* h_susa0_true, TH1D* h_GCF0,TH1D* h_GCF0_true, double max_value, const char* title,const char* name){

  TH1D* h_empirical = (TH1D*)h_empirical0->Clone();
  TH1D* h_nieves = (TH1D*)h_nieves0->Clone();
  TH1D* h_susa = (TH1D*)h_susa0->Clone();
  TH1D* h_GCF = (TH1D*)h_GCF0->Clone();

  area_normalize(h_empirical,h_nieves,h_susa,h_GCF);
  
  TH1D* h_empirical_true = (TH1D*)h_empirical0_true->Clone();
  TH1D* h_nieves_true = (TH1D*)h_nieves0_true->Clone();
  TH1D* h_susa_true = (TH1D*)h_susa0_true->Clone();
  TH1D* h_GCF_true = (TH1D*)h_GCF0_true->Clone();

  area_normalize(h_empirical_true,h_nieves_true,h_susa_true,h_GCF_true);

  //the canvas
  TCanvas* canv_theory = new TCanvas("canv_theory","canv_theory",1000,700);

  h_empirical->Draw("hist");
  h_empirical->SetLineColor(color2);
  h_empirical->SetMaximum(1.0);
  h_empirical->SetTitle(Form("%s",title));
  h_empirical->SetXTitle(Form("%s",title));
  h_empirical->GetXaxis()->SetTitleSize(28);
  h_empirical->GetXaxis()->SetTitleFont(43);
  h_empirical->SetYTitle("Fractional Number of Events");
  h_empirical->GetYaxis()->SetTitleSize(28);
  h_empirical->GetYaxis()->SetTitleFont(43);
  h_empirical->SetMaximum(1.0);

  h_empirical_true->Draw("hist");
  h_empirical_true->SetLineColor(color2);
  h_empirical_true->SetLineStyle(2);

  h_nieves->Draw("hist same");
  h_nieves->SetLineColor(color3);
  h_nieves_true->Draw("hist same");
  h_nieves_true->SetLineColor(color3);
  h_nieves_true->SetLineStyle(2);

  h_susa->Draw("hist same");
  h_susa->SetLineColor(color4);
  h_susa_true->Draw("hist same");
  h_susa_true->SetLineColor(color4);
  h_susa_true->SetLineStyle(2);

  //h_GCF->Draw("hist same");
  //h_GCF->SetLineColor(color5);
  //h_GCF_true->Draw("hist same");
  //h_GCF_true->SetLineColor(color5);
  //h_GCF_true->SetLineStyle(2);

  TLegend* legend_theory = new TLegend(0.14, 0.72, 0.87, 0.87);
  legend_theory->AddEntry(h_empirical,"Empirical: Estimate","L");
  legend_theory->AddEntry(h_empirical_true,"Empirical: Truth Value","L");
  legend_theory->AddEntry(h_nieves,"Nieves: Estimate","L");
  legend_theory->AddEntry(h_nieves_true,"Nieves: Truth Value","L");
  legend_theory->AddEntry(h_susa,"SuSAv2: Estimate","L");
  legend_theory->AddEntry(h_susa_true,"SuSAv2: Truth Value","L");
  //legend_theory->AddEntry(h_GCF,"GCF: Estimate","L") ;
  //legend_theory->AddEntry(h_GCF_true,"GCF: Truth Value","L") ;
  legend_theory->SetLineWidth(0);
  legend_theory->SetFillColor(kWhite);
  legend_theory->SetTextSize(0.03);
  legend_theory->Draw("SAME");
  canv_theory->Print(Form("images/Model_Comparisons/%s_both.png",name));

  
} //end of make plot


void MC_Comparison::model_comparison(TH1D* h_empirical0,TH1D* h_nieves0,TH1D* h_susa0,TH1D* h_GCF0, double max_value, bool yes_fsi, const char* title,const char* name){

  TH1D* h_empirical = (TH1D*)h_empirical0->Clone();
  TH1D* h_nieves = (TH1D*)h_nieves0->Clone();
  TH1D* h_susa = (TH1D*)h_susa0->Clone();
  TH1D* h_GCF = (TH1D*)h_GCF0->Clone();
  
  area_normalize(h_empirical,h_nieves,h_susa,h_GCF);
  make_plot(h_empirical,h_nieves,h_susa,h_GCF,max_value,yes_fsi,title,name);
  cross_section(h_empirical,h_nieves,h_susa,h_GCF);

}//end of make model comparison

