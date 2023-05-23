#include "constants.h"
using namespace Constants;

class NuWro{

 public:
  virtual void main();

 private:

  virtual void Grab_Histograms();
  virtual void Normalize_Matrix(TH2D* h_matrix);
  virtual void Plot_Matrices(TH2D* h_genie_matrix, TH2D* h_nuwro_matrix, const char* title, const char* name);
  virtual void Iteration_Testing(int iter,TH2D* h_matrix, TH1D* h_num, TH1D* h_denom, TH1D* h_nuwro, int max, int max_chi2, const char* title, const char* name, bool flip_legend = false);

  TH2D* h_muon_genie_matrices[num_var];
  TH2D* h_muon_nuwro_matrices[num_var];
  TH1D* h_muon_nuwro_num[num_var];
  TH1D* h_muon_nuwro_denom[num_var];
  int muon_nuwro_max[num_var] = {8000,6000,2500};
  int muon_nuwro_chi2_max[num_var] = {150,6,2};

  TH2D* h_leading_genie_matrices[num_var];
  TH2D* h_leading_nuwro_matrices[num_var];
  TH1D* h_leading_nuwro_num[num_var];
  TH1D* h_leading_nuwro_denom[num_var];
  int leading_nuwro_max[num_var] = {6000,6000,2500};
  int leading_nuwro_chi2_max[num_var] = {150,150,10};
  
  TH2D* h_recoil_genie_matrices[num_var];
  TH2D* h_recoil_nuwro_matrices[num_var];
  TH1D* h_recoil_nuwro_num[num_var];
  TH1D* h_recoil_nuwro_denom[num_var];
  int recoil_nuwro_max[num_var] = {6000,6000,2500};
  int recoil_nuwro_chi2_max[num_var] = {25,20,2};

  TH2D* h_other_genie_matrices[num_other_var];
  TH2D* h_other_nuwro_matrices[num_other_var];
  TH1D* h_other_nuwro_num[num_other_var];
  TH1D* h_other_nuwro_denom[num_other_var];
  int other_nuwro_max[num_other_var] = {6000,6000,5000,12000,10000,10000,20000};
  int other_nuwro_chi2_max[num_other_var] = {20,10,150,150,150,150,1000};
  

};

void NuWro::main(){

  Grab_Histograms();

  for(int i =0; i < num_var; i++){
    Plot_Matrices(h_muon_genie_matrices[i],h_muon_nuwro_matrices[i], Form("Muon %s",var_titles[i]),Form("_muon%s",var[i]));
    Plot_Matrices(h_leading_genie_matrices[i],h_leading_nuwro_matrices[i], Form("Leading %s",var_titles[i]),Form("_leading%s",var[i]));
    Plot_Matrices(h_recoil_genie_matrices[i],h_recoil_nuwro_matrices[i], Form("Recoil %s",var_titles[i]),Form("_recoil%s",var[i]));
    
    Iteration_Testing(9,h_muon_genie_matrices[i],h_muon_nuwro_num[i],h_muon_nuwro_denom[i],h_nuwro_muon[i], muon_nuwro_max[i],muon_nuwro_chi2_max[i],Form("Muon %s",var_titles[i]),Form("_muon%s",var[i]), flip_legend_muon[i]);
    //Iteration_Testing(2,h_muon_nuwro_matrices[i],h_muon_nuwro_num[i],h_muon_nuwro_denom[i],h_nuwro_muon[i], muon_nuwro_max[i],muon_nuwro_chi2_max[i],Form("Muon %s",var_titles[i]),Form("_muon%s",var[i]), flip_legend_muon[i]);


    Iteration_Testing(9,h_leading_genie_matrices[i],h_leading_nuwro_num[i],h_leading_nuwro_denom[i],h_nuwro_leading[i], leading_nuwro_max[i],leading_nuwro_chi2_max[i],Form("Leading %s",var_titles[i]),Form("_leading%s",var[i]), flip_legend_leading[i]);
    //Iteration_Testing(2,h_leading_nuwro_matrices[i],h_leading_nuwro_num[i],h_leading_nuwro_denom[i],h_nuwro_leading[i], leading_nuwro_max[i],leading_nuwro_chi2_max[i],Form("Leading %s",var_titles[i]),Form("_leading%s",var[i]), flip_legend_leading[i]);

    Iteration_Testing(9,h_recoil_genie_matrices[i],h_recoil_nuwro_num[i],h_recoil_nuwro_denom[i],h_nuwro_recoil[i], recoil_nuwro_max[i],recoil_nuwro_chi2_max[i],Form("Recoil %s",var_titles[i]),Form("_recoil%s",var[i]), flip_legend_recoil[i]);
    //Iteration_Testing(2,h_recoil_nuwro_matrices[i],h_recoil_nuwro_num[i],h_recoil_nuwro_denom[i],h_nuwro_recoil[i], recoil_nuwro_max[i],recoil_nuwro_chi2_max[i],Form("Recoil %s",var_titles[i]),Form("_recoil%s",var[i]), flip_legend_recoil[i]);
  }

  for(int i=0; i < num_other_var; i++){
    Plot_Matrices(h_other_genie_matrices[i],h_other_nuwro_matrices[i], Form("%s",other_var_titles[i]),Form("%s",other_var[i]));
    Iteration_Testing(9,h_other_genie_matrices[i],h_other_nuwro_num[i],h_other_nuwro_denom[i],h_nuwro_other[i], other_nuwro_max[i],other_nuwro_chi2_max[i],Form("Other %s",other_var_titles[i]),Form("_other%s",other_var[i]), flip_legend_other[i]);
    //Iteration_Testing(2,h_other_nuwro_matrices[i],h_other_nuwro_num[i],h_other_nuwro_denom[i],h_nuwro_other[i], other_nuwro_max[i],other_nuwro_chi2_max[i],Form("Other %s",other_var_titles[i]),Form("_other%s",other_var[i]), flip_legend_other[i]);
  }

  
}//end of main


void NuWro::Grab_Histograms(){

  TFile* f_genie_matrices = new TFile("../root_files/pelee/Run_all/histograms_mc_eff.root");
  TFile* f_nuwro_matrices = new TFile("../root_files/nuwro/Run_all/histograms_mc_eff.root");
  TFile* f_nuwro = new TFile("../root_files/nuwro/Run_all/histograms_nuwro_xsec_overlay_wgt.root");

  for(int i=0; i < num_var; i++){

    h_muon_genie_matrices[i] = (TH2D*)f_genie_matrices->Get(Form("h_particle_matrices_muon_all%s",(var[i]))); //genie migration matrix
    h_muon_nuwro_matrices[i] = (TH2D*)f_nuwro_matrices->Get(Form("h_particle_matrices_muon_all%s",(var[i]))); //nuwro migration matrix  
    h_muon_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_muon_all%s",var[i])); //true nuwro cc2p0pi
    h_muon_nuwro_num[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_num_muon_all%s",var[i])); //true nuwro cc2p0pi
    h_nuwro_muon[i] = (TH1D*)f_nuwro->Get(Form("h_muon%s_cc2p0pi",var0[i])); //reconstructed nuwro cc2p0pi

    h_leading_genie_matrices[i] = (TH2D*)f_genie_matrices->Get(Form("h_particle_matrices_lead_proton%s",(var[i]))); //genie migration matrix
    h_leading_nuwro_matrices[i] = (TH2D*)f_nuwro_matrices->Get(Form("h_particle_matrices_lead_proton%s",(var[i]))); //nuwro migration matrix  
    h_leading_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_lead_proton%s",var[i])); //true nuwro cc2p0pi
    h_leading_nuwro_num[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_num_lead_proton%s",var[i])); //true nuwro cc2p0pi
    h_nuwro_leading[i] = (TH1D*)f_nuwro->Get(Form("h_leading%s_cc2p0pi",var0[i])); //reconstructed nuwro cc2p0pi

    h_recoil_genie_matrices[i] = (TH2D*)f_genie_matrices->Get(Form("h_particle_matrices_recoil_proton%s",(var[i]))); //genie migration matrix
    h_recoil_nuwro_matrices[i] = (TH2D*)f_nuwro_matrices->Get(Form("h_particle_matrices_recoil_proton%s",(var[i]))); //nuwro migration matrix  
    h_recoil_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_denom_recoil_proton%s",var[i])); //true nuwro cc2p0pi
    h_recoil_nuwro_num[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_particle_num_recoil_proton%s",var[i])); //true nuwro cc2p0pi
    h_nuwro_recoil[i] = (TH1D*)f_nuwro->Get(Form("h_recoil%s_cc2p0pi",var0[i])); //reconstructed nuwro cc2p0pi
    
  }

  for(int i=0; i < num_other_var; i++){


    std::cout<<Form("h_other_matrices%s",(other_var[i]))<<std::endl;

    h_other_genie_matrices[i] = (TH2D*)f_genie_matrices->Get(Form("h_other_matrices%s",(other_var[i]))); //genie migration matrix
    h_other_nuwro_matrices[i] = (TH2D*)f_nuwro_matrices->Get(Form("h_other_matrices%s",(other_var[i]))); //nuwro migration matrix  
    h_other_nuwro_denom[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_other_eff_denom%s",other_var[i])); //true nuwro cc2p0pi
    h_other_nuwro_num[i] = (TH1D*)f_nuwro_matrices->Get(Form("h_other_eff_num%s",other_var[i])); //true nuwro cc2p0pi
    h_nuwro_other[i] = (TH1D*)f_nuwro->Get(Form("h%s_cc2p0pi",other_var[i])); //reconstructed nuwro cc2p0pi
    
  }


  
} //end of grab histograms


void NuWro::Normalize_Matrix(TH2D* h_matrix){


  int NBinsX = h_matrix->GetXaxis()->GetNbins();
  int NBinsY = h_matrix->GetYaxis()->GetNbins();
  
  for(int bin_true = 1; bin_true < NBinsX+1; bin_true++){
    double NEventsInColumn = 0;
    
    for(int bin_reco = 1; bin_reco < NBinsY+1; bin_reco++){
      double bin_content = h_matrix->GetBinContent(bin_true,bin_reco);
      NEventsInColumn += bin_content;
    }
    
    for(int bin_reco = 1; bin_reco < NBinsY+1; bin_reco++){
      if(NEventsInColumn > 0){
	double FracErr =  h_matrix->GetBinError(bin_true,bin_reco) / h_matrix->GetBinContent(bin_true,bin_reco);
	double CV = double( h_matrix->GetBinContent(bin_true,bin_reco))/ double(NEventsInColumn);
	h_matrix->SetBinContent(bin_true,bin_reco,CV);
	
	double error = CV * TMath::Sqrt( TMath::Power(FracErr,2.) + 1./double(NEventsInColumn) );
	h_matrix->SetBinError(bin_true,bin_reco,error) ;
	
      } else {
	h_matrix->SetBinContent(bin_true,bin_reco);
	h_matrix->SetBinError(bin_true,bin_reco);
      } //ends else                                                                                                                                                                                            
    }
  }
  

} //end of normalize_matrix

void NuWro::Plot_Matrices(TH2D* h_genie_matrix, TH2D* h_nuwro_matrix, const char* title, const char* name){
  
  Normalize_Matrix(h_genie_matrix);
  Normalize_Matrix(h_nuwro_matrix);
  
  TCanvas* c_2d =  new TCanvas("c_2d","Smearing Matrix",800,1350);
  c_2d->SetRightMargin(0.09);
  c_2d->SetLeftMargin(0.15);
  c_2d->SetBottomMargin(0.15);
  c_2d->Divide(1,2);

  c_2d->cd(1);
  TH2D* h_smearing_genie = (TH2D*)h_genie_matrix->Clone();
  h_smearing_genie->Draw("colz text");
  h_smearing_genie->SetTitle(Form("GENIE %s",title));
  h_smearing_genie->GetXaxis()->SetTitle(Form("True %s",title));
  h_smearing_genie->GetYaxis()->SetTitle(Form("Reco. %s",title));
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);


  c_2d->cd(2);
  TH2D* h_smearing_nuwro = (TH2D*)h_nuwro_matrix->Clone();
  h_smearing_nuwro->Draw("colz text");
  h_smearing_nuwro->SetTitle(Form("NuWro %s",title));
  h_smearing_nuwro->GetXaxis()->SetTitle(Form("True %s",title));
  h_smearing_nuwro->GetYaxis()->SetTitle(Form("Reco. %s",title));
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);
  
  c_2d->Print(Form("images/NuWro/smearing_matrices_both%s.png",name));
} //end of plot matrices


void NuWro::Iteration_Testing(int iter,TH2D* h_matrix, TH1D* h_num, TH1D* h_denom, TH1D* h_nuwro, int max,int max_chi2, const char* title, const char* name, bool flip_legend = false){


  Color_t colors[] = {kBlack, kRed, kOrange+1, kYellow-3, kGreen, kCyan, kBlue, kViolet-9, kMagenta, kPink+-9, kGray,};
  std::vector<TH1D*> hist_vec;

  TH1D* h_unfolded_signal_test[iter+1];
  TH1D* h_eff = (TH1D*)h_num->Clone();
  h_eff->Divide(h_eff,h_denom,1,1);
  
  TCanvas* canv_nuwro =  new TCanvas(Form("canv_nuwro_0"),Form("canv_nuwro_0"),800,600);
  h_denom->Draw("hist");
  h_denom->SetLineColor(colors[0]);
  h_denom->SetTitle(title);
  h_denom->GetXaxis()->SetTitle(title);
  h_denom->GetYaxis()->SetTitle("No. Events");
  h_denom->SetMaximum(max);
  
  h_nuwro->Draw("hist same");
  h_nuwro->SetLineColor(colors[1]);
  hist_vec.push_back(h_nuwro);

  TLegend* legend_nuwro;
  if(flip_legend){
    legend_nuwro = new TLegend(0.175, 0.56, 0.684, 0.87);
  }else{
    legend_nuwro = new TLegend(0.37, 0.56, 0.889, 0.87);
  }
  legend_nuwro->SetNColumns(2);
  legend_nuwro->AddEntry(h_denom,"True 1#mu2p Events (NuWro)","lf");
  legend_nuwro->AddEntry(h_nuwro, "Reco. 1#mu2p Events (NuWro)","lf");

  RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  Response_Matrix.UseOverflow(kTRUE); //specify to include contributions from overflow bins. 
  
  for(int i=1; i < iter; i++){
    RooUnfoldBayes RooUnfoldBayes_Test(&Response_Matrix, h_nuwro, i); //unfold the nuwro data set using our response matrix and number of iterations                                                               
    h_unfolded_signal_test[i] = (TH1D*) RooUnfoldBayes_Test.Hreco();
    h_unfolded_signal_test[i]->Divide(h_unfolded_signal_test[i],h_eff,1,1,"cp");
    hist_vec.push_back(h_unfolded_signal_test[i]);
    
    h_unfolded_signal_test[i]->Draw("SAME 1e1p");
    h_unfolded_signal_test[i]->SetLineColor(colors[i+1]);
    legend_nuwro->AddEntry(h_unfolded_signal_test[i], Form("%d Iteration Unfolded NuWro",i),"lf");

  }

  legend_nuwro->Draw("same");
  canv_nuwro->Print(Form("images/NuWro/%s_nuwro_all.png",name));

  //Now here is where we calculate the chi2 plot
  /////////////////////////////////////////////
  const int nbins = h_denom->GetNbinsX();
  Int_t n = iter; //number of iterations we are testing
  Double_t iterations[iter]; //vector of the number of iterations so we can make the graph
  Double_t chi2_value[iter]; //chi2 value as a function of the number of iterations
  Double_t chi2_dof_value[iter];

  for(int j=0; j < iter; j++){
    double chi2 = 0;

    for(int i=1; i < nbins+1; i++){
      
      double unfolded = hist_vec[j]->GetBinContent(i);
      double truth = h_denom->GetBinContent(i);

      if(truth == 0.0){
	chi2 += 0;
      } else{
	chi2+= std::pow((unfolded - truth),2) / truth;
      }
      
    } //end of loop over bins

    iterations[j] = j;
    chi2_value[j] = chi2;
    chi2_dof_value[j] = chi2/nbins;
    
  } //end of loop over iterations

  
  TCanvas* canv_chi2 = new TCanvas("canv_chi2","canv_chi2",800,600);
  canv_chi2->SetGrid();
  canv_chi2->SetRightMargin(0.09);
  canv_chi2->SetLeftMargin(0.15);
  canv_chi2->SetBottomMargin(0.15);

  TGraph* gra_chi2 =  new TGraph(n,iterations,chi2_value);
  gra_chi2->Draw("AL*");
  gra_chi2->SetMarkerColor(kBlue);
  gra_chi2->SetMarkerSize(2);

  TLine* a = new TLine(0,1,iter-1,1);
  a->Draw("same");
  a->SetLineColor(kRed);

  gra_chi2->GetXaxis()->SetTitle("Number of Iterations");
  gra_chi2->GetYaxis()->SetTitle("#chi^{2} = #sum_{i=1}^{NBins+1} #frac{(Unfolded_{i} - Truth_{i})^{2}}{Truth_{i}}");
  gra_chi2->SetTitle(Form("Unfolded CC2p NuWro: %s",title));
  gra_chi2->SetMinimum(-10);
  //gra_chi2->SetMaximum(100);
  canv_chi2->Print(Form("images/NuWro/%s_chi2.png",name));

  TCanvas* canv_chi2_dof = new TCanvas("canv_chi2_dof","canv_chi2_dof",800,600);
  canv_chi2_dof->SetGrid();
  canv_chi2_dof->SetRightMargin(0.09);
  canv_chi2_dof->SetLeftMargin(0.15);
  canv_chi2_dof->SetBottomMargin(0.15);

  TGraph* gra_chi2_dof =  new TGraph(n,iterations,chi2_dof_value);
  gra_chi2_dof->Draw("AL*");
  gra_chi2_dof->SetMarkerColor(kBlue);
  gra_chi2_dof->SetMarkerSize(2);

  TLine* a1 = new TLine(0,1,iter-1,1);
  a1->Draw("same");
  a1->SetLineColor(kRed);

  gra_chi2_dof->GetXaxis()->SetTitle("Number of Iterations");
  gra_chi2_dof->GetYaxis()->SetTitle("#chi^{2} / No. Bins");
  gra_chi2_dof->SetTitle(Form("Unfolded CC2p NuWro: %s",title));
  gra_chi2_dof->SetMinimum(-10);
  //gra_chi2_dof->SetMaximum(20);
  canv_chi2_dof->Print(Form("images/NuWro/%s_chi2_dof.png",name));

  Double_t x[iter-1];
  Double_t y[iter-1];
	    
  for(int i = 0; i < iter-1; i++){
    x[i] = i;
    y[i] = chi2_dof_value[i+1] - chi2_dof_value[i];
    std::cout<<"Value of y[i] for "<<i<<" iterations: "<<y[i]<<std::endl;
  }

  TCanvas* canv_chi2_diff = new TCanvas("canv_chi2_diff","canv_chi2_diff",800,600);
  canv_chi2_diff->SetGrid();
  canv_chi2_diff->SetRightMargin(0.09);
  canv_chi2_diff->SetLeftMargin(0.15);
  canv_chi2_diff->SetBottomMargin(0.15);
  
  TGraph* gra_chi2_diff =  new TGraph(n-1,x,y);
  gra_chi2_diff->Draw("AL*");
  gra_chi2_diff->SetMarkerColor(kBlue);
  gra_chi2_diff->SetMarkerSize(2);

  TLine* a2 = new TLine(0,0,iter-2,0);
  a2->Draw("same");
  a2->SetLineColor(kRed);

  gra_chi2_diff->GetXaxis()->SetTitle("Number of Iterations");
  gra_chi2_diff->GetYaxis()->SetTitle("#chi^{2}_{Iter+1} - #chi^{2}_{Iter}");
  gra_chi2_diff->SetTitle(Form("Unfolded CC2p NuWro: %s",title));
  gra_chi2_diff->SetMinimum(-2);
  gra_chi2_diff->SetMaximum(max_chi2);
  canv_chi2_diff->Print(Form("images/NuWro/%s_chi2_diff.png",name));
 
} //end of iteration testing
