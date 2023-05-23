#include "constants.h"
using namespace Constants;

class iterations{

 public:
  virtual void Iteration_Test(int iterations, TH1D* h_bnb1, TH1D* h_ext1, TH1D* h_overlay_total1, TH1D* h_overlay_cc2p1, TH1D* h_dirt1,TH1D* h_nuwro, TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, const char* title,bool print_contents);

private:

  //virtual void Grab_Histograms(const char* particle, const char* variable, const char* var);
  //virtual std::vector<TH1D*> Mean_Generation(const char* particle, const char* variable,const char* var,int iterations, TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents);
  virtual std::vector<TH1D*> NuWro_Generation(int iterations, TH1D* h_nuwro,TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents);
  virtual std::vector<TH1D*> GENIE_CV_Generation(int iterations,TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents);
  virtual std::vector<TH1D*> BNB_Generation(int iterations,TH1D* h_bnb1, TH1D* h_ext1, TH1D* h_overlay_total1, TH1D* h_overlay_cc2p1, TH1D* h_dirt1,TH1D* h_eff, TH2D* h_matrix, const char* name);
  //virtual void BNB_Testing(int iterations,TH1D* h_bnb1, TH1D* h_ext1, TH1D* h_overlay_total1, TH1D* h_overlay_cc2p1, TH1D* h_dirt1,TH1D* h_eff, TH2D* h_matrix, const char* particle, const char* variable);
  
  static const int num_universes = 500;
  TH2D* h_genie_matrices[500];
  TH1D* h_genie_multisims[num_universes];
  TH1D* h_genie_num[num_universes];
  TH1D* h_genie_denom[num_universes];
  TH1D* h_genie_eff[num_universes];
 
};

/*void iterations::Grab_Histograms(const char* particle, const char* variable, const char* var){

  //Grab all of our multisims
  /////////////////////////////////
  TFile* file = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/All_UBGenie/histograms_pelee_xsec.root");
  TFile* file_efficiency = new TFile("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/All_UBGenie/histograms_mc_eff.root");
  
  for(Int_t i=0; i < num_universes; i++){    

    if(std::strcmp(particle,"muon") ==0){
      h_genie_matrices[i] = (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_%d",var,i));
      h_genie_multisims[i] = (TH1D*)file->Get(Form("h_%s%s_cc2p0pi_%d",particle,variable,i));
      h_genie_num[i] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_%d",var,i));
      h_genie_denom[i] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_%d",var,i));

    } else if (std::strcmp(particle,"leading") == 0){
      h_genie_matrices[i] = (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_%d",var,i));
      h_genie_multisims[i] = (TH1D*)file->Get(Form("h_%s%s_cc2p0pi_%d",particle,variable,i));
      h_genie_num[i] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_%d",var,i));
      h_genie_denom[i] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_%d",var,i));

    } else if (std::strcmp(particle,"recoil") == 0){
       h_genie_matrices[i] = (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_%d",var,i));
      h_genie_multisims[i] = (TH1D*)file->Get(Form("h_%s%s_cc2p0pi_%d",particle,variable,i));
      h_genie_num[i] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_%d",var,i));
      h_genie_denom[i] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_%d",var,i));
      
    } else {
      h_genie_matrices[i] = (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_%d",var,i));
      h_genie_multisims[i] = (TH1D*)file->Get(Form("h%s_cc2p0pi_%d",variable,i));
      h_genie_num[i] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_%d",var,i));
      h_genie_denom[i] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_%d",var,i));
    }

    h_genie_eff[i] = (TH1D*)h_genie_num[i]->Clone();
    h_genie_eff[i]->Divide(h_genie_eff[i],h_genie_denom[i],1,1,"cp");
  }
  
}//end of grab histograms
*/

/*std::vector<TH1D*> iterations::Mean_Generation(const char* particle, const char* variable,const char* var,int iterations, TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents){

  static const int iter = iterations;

  //Grab all of our multisims
  /////////////////////////////////
  Grab_Histograms(particle,variable,var);
  
  //Defining all the gunk
  //////////////////////////////////////////////
  TH1D* h_unfolded_signal_test[iter+1][num_universes];
  TH1D* h_mean[iter+1]; //mean histogram that will be a function of the number of iterations
  int n_bins = h_denom->GetNbinsX();
  std::vector<TH1D*> hist_mean;

  //Make the 0th order plot
  /////////////////////////////
  h_mean[0] = (TH1D*)h_denom->Clone();

  for(int k=1; k < n_bins+1; k++){
    double sum = 0;
    
    for(int j=0; j < num_universes; j++){ 
      double value = h_genie_multisims[j]->GetBinContent(k);
      sum += value;
    }
    
    double mean = sum / num_universes;
    h_mean[0]->SetBinContent(k,mean);
  }
  hist_mean.push_back(h_mean[0]);
  
  TCanvas* canv_mean =  new TCanvas(Form("canv_mean_0"),Form("canv_mean_0"),800,600);
  h_denom->Draw("hist");
  h_denom->SetLineColor(kBlack);
  
  for(int j=0; j < num_universes; j ++){
    h_genie_multisims[j]->Draw("hist same");
    h_genie_multisims[j]->SetLineColor(kGreen);
  }
  
  h_denom->Draw("hist same");
  h_denom->SetLineColor(kBlack);
  
  h_mean[0]->Draw("hist p same");
  h_mean[0]->SetLineColor(kBlue);
  h_mean[0]->SetMarkerStyle(20);
  h_mean[0]->SetMarkerSize(1);
  
  TLegend* legend_mean = new TLegend(0.71, 0.54, 0.899, 0.89);
  legend_mean->AddEntry(h_denom,"All Generated 1#mu2p Events","lf");
  legend_mean->AddEntry(h_genie_multisims[0], "Multisims","lf");
  legend_mean->AddEntry(h_mean[0],"Unfolded Multisims Average","p");
  legend_mean->Draw("same");
  canv_mean->Print(Form("images/iteration_test/%s%s_multisims_0.png",particle,variable));  

  
  //Now to make the iteration plots
  //////////////////////////////////
  RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  
  for(int i=1; i < iter+1; i++){
    for(int j=0; j < num_universes; j++){
      RooUnfoldBayes RooUnfoldBayes_Test(&Response_Matrix, h_genie_multisims[j], i); //unfold the nuwro data set using our response matrix and number of iterations
      h_unfolded_signal_test[i][j] = (TH1D*) RooUnfoldBayes_Test.Hreco();
      h_unfolded_signal_test[i][j]->Divide(h_unfolded_signal_test[i][j],h_eff,1,1,"cp");
    }
    
    h_mean[i] = (TH1D*)h_unfolded_signal_test[i][0]->Clone();
    
    for(int k=1; k < n_bins+1; k++){
      double sum = 0;
    
      for(int j=0; j < num_universes; j++){ 
	double value = h_unfolded_signal_test[i][j]->GetBinContent(k);
	sum += value;
      }

      double mean = sum / num_universes;
      h_mean[i]->SetBinContent(k,mean);
    }

    hist_mean.push_back(h_mean[i]);
    
    TCanvas* canv_mean =  new TCanvas(Form("canv_mean_%d",i),Form("canv_mean_%d",i),800,600);
    h_denom->Draw("hist");
    h_denom->SetLineColor(kBlack);

    for(int j=0; j < num_universes; j ++){
      h_unfolded_signal_test[i][j]->Draw("hist same");
      h_unfolded_signal_test[i][j]->SetLineColor(kGreen);
    }
    
    h_denom->Draw("hist same");
    h_denom->SetLineColor(kBlack);

    h_mean[i]->Draw("hist p same");
    h_mean[i]->SetLineColor(kBlue);
    h_mean[i]->SetMarkerStyle(20);
    h_mean[i]->SetMarkerSize(1);
    
    TLegend* legend_mean = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_mean->AddEntry(h_denom,"All Generated 1#mu2p Events","lf");
    legend_mean->AddEntry(h_unfolded_signal_test[i][0], Form("%d Iteration Unfolded Multisims",i),"lf");
    legend_mean->AddEntry(h_mean[i],"Unfolded Multisims Average","p");
    legend_mean->Draw("same");
    canv_mean->Print(Form("images/iteration_test/%s%s_multisims_%d.png",particle,variable,i));  

  } //end of loop over iterations

  return hist_mean;

} //end of mean_generation 
*/

std::vector<TH1D*> iterations::NuWro_Generation(int iterations, TH1D* h_nuwro,TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents){
  
  //Defining all the gunk
  ///////////////////////
  static const int iter = iterations;
  TH1D* h_unfolded_signal_test[iter+1];
  std::vector<TH1D*> hist_nuwro;

  //Make the 0th order plot
  /////////////////////////////
  hist_nuwro.push_back(h_nuwro);
  
  TCanvas* canv_nuwro =  new TCanvas(Form("canv_nuwro_0"),Form("canv_nuwro_0"),800,600);
  h_denom->Draw("hist");
  h_denom->SetLineColor(kBlack);
  h_nuwro->Draw("hist same");
  h_nuwro->SetLineColor(kGreen);
  
  TLegend* legend_nuwro = new TLegend(0.71, 0.54, 0.899, 0.89);
  legend_nuwro->AddEntry(h_denom,"All Generated 1#mu2p Events in Overlay","lf");
  legend_nuwro->AddEntry(h_nuwro, "All Generated 1#mu2p Events in NuWro","lf");
  legend_nuwro->Draw("same");
  canv_nuwro->Print(Form("images/iteration_test/%s_nuwro_0.png",name));  

  //Now to do the unfolding
  //////////////////////////
  RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  Response_Matrix.UseOverflow(kTRUE); //specify to include contributions from overflow bins
  
  for(int i=1; i < iter+1; i++){
    RooUnfoldBayes RooUnfoldBayes_Test(&Response_Matrix, h_nuwro, i); //unfold the nuwro data set using our response matrix and number of iterations
    h_unfolded_signal_test[i] = (TH1D*) RooUnfoldBayes_Test.Hreco();
    h_unfolded_signal_test[i]->Divide(h_unfolded_signal_test[i],h_eff,1,1,"cp");
    hist_nuwro.push_back(h_unfolded_signal_test[i]);
  

    TCanvas* canv_nuwro =  new TCanvas(Form("canv_nuwro_%d",i),Form("canv_nuwro_%d",i),800,600);
    h_denom->Draw("hist");
    h_denom->SetLineColor(kBlack);
    h_unfolded_signal_test[i]->Draw("SAME HIST");
    h_unfolded_signal_test[i]->SetLineColor(kGreen);
  
    TLegend* legend_nuwro = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_nuwro->AddEntry(h_denom,"All Generated 1#mu2p Events","lf");
    legend_nuwro->AddEntry(h_unfolded_signal_test[i], Form("%d Iteration Unfolded NuWro",i),"lf");
    legend_nuwro->Draw("same");
    canv_nuwro->Print(Form("images/iteration_test/%s_nuwro_%d.png",name,i));  
  } //end of loop over iterations

  return hist_nuwro;
  hist_nuwro.clear();
  
} //end of NuWro_Genration

std::vector<TH1D*> iterations::GENIE_CV_Generation(int iterations,TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, bool print_contents){

  //Defining all the gunk
  ///////////////////////
   static const int iter = iterations;
  TH1D* h_unfolded_signal_test[iter+1];
  std::vector<TH1D*> hist_genie_CV;

  //Pushing back the 0th order histogram
  /////////////////////////////
  hist_genie_CV.push_back(h_denom);
  
  //Now to do the unfolding
  //////////////////////////
  RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  Response_Matrix.UseOverflow(kTRUE); //specify to include contributions from overflow bins
  
  for(int i=1; i < iter+1; i++){
    RooUnfoldBayes RooUnfoldBayes_Test(&Response_Matrix, h_denom, i); //unfold the nuwro data set using our response matrix and number of iterations
    h_unfolded_signal_test[i] = (TH1D*) RooUnfoldBayes_Test.Hreco();
    h_unfolded_signal_test[i]->Divide(h_unfolded_signal_test[i],h_eff,1,1,"cp");
    hist_genie_CV.push_back(h_unfolded_signal_test[i]);
  
    TCanvas* canv_CV =  new TCanvas(Form("canv_CV_%d",i),Form("canv_CV_%d",i),800,600);
    h_denom->Draw("hist");
    h_denom->SetLineColor(kBlack);
    h_unfolded_signal_test[i]->Draw("SAME HIST");
    h_unfolded_signal_test[i]->SetLineColor(kGreen);
  
    TLegend* legend_CV = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_CV->AddEntry(h_denom,"All Generated 1#mu2p Events","lf");
    legend_CV->AddEntry(h_unfolded_signal_test[i], Form("%d Iteration Unfolded GENIE CV",i),"lf");
    legend_CV->Draw("same");
    canv_CV->Print(Form("images/iteration_test/%s_CV_%d.png",name,i));  
  } //end of loop over iterations

  return hist_genie_CV;
  hist_genie_CV.clear();
}

/*void iterations::BNB_Testing(int num_iter_1, int num_iter_2, TH1D* h_bnb, TH1D* h_eff, TH2D* h_matrix, const char* particle, const char* variable){
  
  //Response Matrix
   RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  
  TH1D* h_num_1;
  RooUnfoldBayes RooUnfoldBayes_Num1(&Response_Matrix, h_bnb, num_iter_1); 
  h_num_1 = (TH1D*) RooUnfoldBayes_Num1.Hreco();
  h_num_1->Divide(h_num_1,h_eff,1,1,"cp");

  TH1D* h_num_2;
  RooUnfoldBayes RooUnfoldBayes_Num2(&Response_Matrix, h_bnb, num_iter_2); 
  h_num_2 = (TH1D*) RooUnfoldBayes_Num2.Hreco();
  h_num_2->Divide(h_num_2,h_eff,1,1,"cp");


  TCanvas* canv_compare = new TCanvas("canv_compare","canv_compare",800,600);
  h_num_1->Draw("hist");
  h_num_1->SetLineColor(kRed);
  h_num_2->Draw("same hist");
  h_num_2->SetLineColor(kBlue);
  TLegend* legend_compare = new TLegend(0.71, 0.54, 0.899, 0.89);
  legend_compare->AddEntry(h_num_1,Form("%d Iterations Unfolded",num_iter_1),"lpf");
  legend_compare->AddEntry(h_num_2,Form("%d Iterations Unfolded",num_iter_2),"lpf");
  legend_compare->Draw("same");
  canv_compare->Print(Form("images/iteration_test/%s%s_bnb.png",particle,variable));

  const Int_t n_bins = h_bnb->GetNbinsX();
  double chi2;

  TH1D* h_diff = (TH1D*)h_num_1->Clone();
  TH1D* h_BNB = (TH1D*)h_bnb->Clone();
  h_BNB->Divide(h_BNB,h_eff,1,1,"cp");


  for(int i=1; i < n_bins+1; i++){
    double E = h_num_2->GetBinContent(i);
    double O = h_num_1->GetBinContent(i);
    double BNB = h_BNB->GetBinContent(i);
    if(E == 0){
      chi2 = 0;
    }else{
      chi2 += std::pow(O - E,2) / O;
    }

    double value = (E-O)/BNB;
    h_diff->SetBinContent(i,value);
  } //loop over the bins
  
  TCanvas* canv_diff = new TCanvas("canv_diff","canv_diff",800,600);
  h_diff->Draw("hist text");
  h_diff->GetXaxis()->SetTitle(Form("h(%d) - h(%d) / h_BNB",num_iter_2,num_iter_1));
  canv_diff->Print(Form("images/iteration_test/%s%s_bnb_diff.png",particle,variable));
  
  std::cout<<"[Iterations.BNB_Testing] Value of the Chi2/N_bins Between "<<num_iter_1<<" iterations and "<<num_iter_2<<" iterations: "<<chi2/n_bins<<std::endl;
  
  }*/

std::vector<TH1D*> iterations::BNB_Generation(int iterations, TH1D* h_bnb1, TH1D* h_ext1, TH1D* h_overlay_total1, TH1D* h_overlay_cc2p1, TH1D* h_dirt1,TH1D* h_eff, TH2D* h_matrix, const char* name){

  //So the first thing we have to do is remove the background events from the bnb
  ///////////////////////////////////////////////////////////////////////////////
  TH1D* h_bnb = (TH1D*)h_bnb1->Clone();
  TH1D* h_ext = (TH1D*)h_ext1->Clone();
  TH1D* h_overlay_total = (TH1D*)h_overlay_total1->Clone();
  TH1D* h_overlay_cc2p = (TH1D*)h_overlay_cc2p1->Clone();
  TH1D* h_dirt = (TH1D*)h_dirt1->Clone();

  h_ext->Add(h_dirt);
  h_overlay_total->Add(h_overlay_cc2p, -1);
  h_ext->Add(h_overlay_total);
  h_bnb->Add(h_ext,-1); 
  
  //Defining all the gunk
  ///////////////////////
   static const int iter = iterations;
  TH1D* h_unfolded_signal_test[iter+1];
  std::vector<TH1D*> hist_BNB;

  //Pushing back the 0th order histogram
  /////////////////////////////
  hist_BNB.push_back(h_bnb);
  
  //Now to do the unfolding
  //////////////////////////
  RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  Response_Matrix.UseOverflow(kTRUE); //specify to include contributions from overflow bins
  
  for(int i=1; i < iter+1; i++){
    RooUnfoldBayes RooUnfoldBayes_Test(&Response_Matrix, h_bnb, i); //unfold the nuwro data set using our response matrix and number of iterations
    h_unfolded_signal_test[i] = (TH1D*) RooUnfoldBayes_Test.Hreco();
    h_unfolded_signal_test[i]->Divide(h_unfolded_signal_test[i],h_eff,1,1,"cp");
    hist_BNB.push_back(h_unfolded_signal_test[i]);
  
    TCanvas* canv_bnb =  new TCanvas(Form("canv_bnb_%d",i),Form("canv_bnb_%d",i),800,600);
    h_bnb->Draw("hist");
    h_bnb->SetLineColor(kBlack);
    h_unfolded_signal_test[i]->Draw("SAME HIST");
    h_unfolded_signal_test[i]->SetLineColor(kGreen);
  
    TLegend* legend_bnb = new TLegend(0.71, 0.54, 0.899, 0.89);
    legend_bnb->AddEntry(h_bnb,"N_{BNB} - N_{Background}","lf");
    legend_bnb->AddEntry(h_unfolded_signal_test[i], Form("%d Iteration Unfolded Data",i),"lf");
    legend_bnb->Draw("same");
    canv_bnb->Print(Form("images/iteration_test/%s_bnb_%d.png",name,i));  
  } //end of loop over iterations

  return hist_BNB;
  hist_BNB.clear();

  /*
  //Response Matrix
   RooUnfoldResponse Response_Matrix(0,0,h_matrix,"Response_Matrix"); //the response matrix
  
  TH1D* h_num_1;
  RooUnfoldBayes RooUnfoldBayes_Num1(&Response_Matrix, h_bnb, num_iter_1); 
  h_num_1 = (TH1D*) RooUnfoldBayes_Num1.Hreco();
  h_num_1->Divide(h_num_1,h_eff,1,1,"cp");

  TH1D* h_num_2;
  RooUnfoldBayes RooUnfoldBayes_Num2(&Response_Matrix, h_bnb, num_iter_2); 
  h_num_2 = (TH1D*) RooUnfoldBayes_Num2.Hreco();
  h_num_2->Divide(h_num_2,h_eff,1,1,"cp");


  TCanvas* canv_compare = new TCanvas("canv_compare","canv_compare",800,600);
  h_num_1->Draw("hist");
  h_num_1->SetLineColor(kRed);
  h_num_2->Draw("same hist");
  h_num_2->SetLineColor(kBlue);
  TLegend* legend_compare = new TLegend(0.71, 0.54, 0.899, 0.89);
  legend_compare->AddEntry(h_num_1,Form("%d Iterations Unfolded",num_iter_1),"lpf");
  legend_compare->AddEntry(h_num_2,Form("%d Iterations Unfolded",num_iter_2),"lpf");
  legend_compare->Draw("same");
  canv_compare->Print(Form("images/iteration_test/%s%s_bnb.png",particle,variable));

  const Int_t n_bins = h_bnb->GetNbinsX();
  double chi2;

  TH1D* h_diff = (TH1D*)h_num_1->Clone();
  TH1D* h_BNB = (TH1D*)h_bnb->Clone();
  h_BNB->Divide(h_BNB,h_eff,1,1,"cp");


  for(int i=1; i < n_bins+1; i++){
    double E = h_num_2->GetBinContent(i);
    double O = h_num_1->GetBinContent(i);
    double BNB = h_BNB->GetBinContent(i);
    if(E == 0){
      chi2 = 0;
    }else{
      chi2 += std::pow(O - E,2) / O;
    }

    double value = (E-O)/BNB;
    h_diff->SetBinContent(i,value);
  } //loop over the bins
  
  TCanvas* canv_diff = new TCanvas("canv_diff","canv_diff",800,600);
  h_diff->Draw("hist text");
  h_diff->GetXaxis()->SetTitle(Form("h(%d) - h(%d) / h_BNB",num_iter_2,num_iter_1));
  canv_diff->Print(Form("images/iteration_test/%s%s_bnb_diff.png",particle,variable));
  
  std::cout<<"[Iterations.BNB_Testing] Value of the Chi2/N_bins Between "<<num_iter_1<<" iterations and "<<num_iter_2<<" iterations: "<<chi2/n_bins<<std::endl;
  */
}

void iterations::Iteration_Test(int iterations, TH1D* h_bnb1, TH1D* h_ext1, TH1D* h_overlay_total1, TH1D* h_overlay_cc2p1, TH1D* h_dirt1,TH1D* h_nuwro, TH1D* h_eff,TH1D* h_denom, TH2D* h_matrix, const char* name, const char* title,bool print_contents){

  const Int_t n_bins = h_denom->GetNbinsX();
  const Int_t n = iterations+1;
  Double_t chi2_values[n];
  Double_t chi2_reg_values[n];
  std::vector<TH1D*> histogram_vec;
  static const int num_types = 3; 
  const char* type_name_vec[num_types] = {"nuwro","genie_CV","bnb"};
  const char* theory_vec[num_types] = {"NuWro","Overlay","BNB"};
  
  for(int type = 2 ; type < 3; type++){

    const char* type_name = type_name_vec[type];
    const char* theory = theory_vec[type];
    
    if(type == 0){
      histogram_vec = NuWro_Generation(iterations, h_nuwro, h_eff,h_denom, h_matrix, name, print_contents); //if you want to use Nuwro predictions
      
    } else if (type == 1){
      histogram_vec = GENIE_CV_Generation(iterations, h_eff, h_denom, h_matrix, name, print_contents); //if you want to use our CV predictions

    } else if (type == 2){
      histogram_vec = BNB_Generation(iterations,h_bnb1, h_ext1, h_overlay_total1, h_overlay_cc2p1, h_dirt1, h_eff, h_matrix, name); //if you want to use our CV predictions
    }
    
    /*else if(type == 3){
      histogram_vec = Mean_Generation(particle,  variable, var, iterations, h_eff,h_denom, h_matrix,  name, print_contents);
      type_name = "multisim";
    }
    */

    for(int x=1; x < iterations+1; x++){
      double chi2 = 0;
      TH1D* h_Observed = (TH1D*)histogram_vec[x-1]->Clone();
      TH1D* h_Observed_x = (TH1D*)histogram_vec[x]->Clone();
      TH1D* h_Observed_x_1 = (TH1D*)histogram_vec[x-1]->Clone();
      
      for(int i=1; i < n_bins+1; i++){
	double O_X = h_Observed_x->GetBinContent(i);
	double O_X_1 = h_Observed_x_1->GetBinContent(i);
	if(O_X == 0){
	  chi2 = 0;
	}else{
	  chi2 += std::pow(O_X - O_X_1,2) / O_X;
	}
      } //loop over the bins
      
      double chi2_reg = 0;
      for(int i=2; i < n_bins+1; i++){
	double x = h_Observed->GetBinContent(i);
	double y = h_Observed->GetBinContent(i-1);
	double value = std::pow(x-y,2);
	chi2_reg += value;
      }
      
      chi2_reg_values[x-1] = chi2_reg;
      chi2_values[x-1] = chi2/n_bins;
      if(print_contents) std::cout<<"[Iterations.chi_square_test] Value of chi2/n_bins at Iteration "<<x<<": "<<chi2/n_bins<<std::endl;
      if(print_contents) std::cout<<"[Iterations.chi_square_test] Value of chi2 at Iteration "<<x<<": "<<chi2 <<std::endl;
      if(print_contents) std::cout<<"[Iterations.chi_square_test] Value of chi2 reg values at Iteration "<<x<<": "<<chi2_reg<<std::endl;
      
    } //loop over iteratiions 
    
    
    //Set the x and y for the TGraph
    ////////////////////////////
    Double_t x[n-1];
    Double_t y[n-1];
    
    for(int i=1; i< iterations+1; i++){    
      x[i-1] = chi2_values[i];
      y[i-1] = log(chi2_reg_values[i]);
    }
    
    //Now to actually Draw
    ////////////////////////
    Color_t colors[10] = {kBlack,kRed, kOrange, kYellow, kGreen, kCyan, kBlue, kViolet-3, kMagenta, kGray};
    
    TCanvas* canv_L_curve = new TCanvas("canv_L_curve","canv_L_curve",800,600);
    canv_L_curve->SetGrid();
    TGraph* gra =  new TGraph(n-2,x,y);
    gra->Draw("AL*");
    gra->SetMarkerColor(kBlue);
    gra->SetMarkerSize(2);
    gra->GetXaxis()->SetTitle("#chi^{2}_{Fit}");
    gra->GetYaxis()->SetTitle("Spikey-ness of  Result");
    gra->SetTitle(Form("Unfolded %s CC2p: %s",theory,title));
    canv_L_curve->Print(Form("images/iteration_test/%s_%s_L_Curve.png",name,type_name));


    /*
    //Now to determine the point closest to the origin
    //////////////////////////////////////////////////
    std::vector<std::pair<double,int>> zipped1;

    for(int i=1; i< iterations+1; i++){
      double norm_x = log(chi2_values[i]);
      double norm_y = log(chi2_reg_values[i]);
      double distance = std::pow((std::pow(norm_x,2) + std::pow(norm_y,2)),2);
      std::cout<<"Value of x[i]: "<<chi2_values[i]<<" Value of x[n-2]: "<<chi2_values[iterations-1]<< "Value of x[i]/x[n-2]: "<<norm_x<<std::endl;
      std::cout<<"Value of y[i]: "<<chi2_reg_values[i]<<" Value of y[n-2]: "<<chi2_reg_values[iterations-1]<< "Value of y[i]/y[n-2]: "<<norm_y<<std::endl;
      std::cout<<"Value of distance at iteration "<<i<<" :"<<distance<<std::endl;
      zipped1.push_back(std::make_pair(distance,i));
    }

    std::sort(zipped1.begin(), zipped1.end());
    double distance_value;
    int iterations_value;
    
    if(zipped1[0].first != 0){
      distance_value = zipped1[0].first;
      iterations_value = zipped1[0].second;
    } else{
      distance_value = zipped1[1].first;
      iterations_value = zipped1[1].second;
    }
    std::cout<<"Value of FINAL distance value at iteration "<<iterations_value<<" : "<<distance_value<<std::endl; 
    zipped1.clear();
    */


    //Going to try and determmine the gradient
    //////////////////////////////////////////
    std::vector<std::pair<double,int>> zipped_gradient;
    
    for(int i=1; i <iterations-1; i++){
      double gradient_up = (x[i] - x[i+1]) / (y[i] - y[i+1]);
      double gradient_down = (x[i] - x[i-1]) / (y[i] - y[i-1]);
      double gradient_diff = gradient_up - gradient_down;
      //std::cout<<"Value of i: "<<i<<std::endl;
      //std::cout<<"Value of x[i]: "<<x[i]<<"Value of x[i-1]: "<<x[i-1]<<"Value of x[i+1]: "<<x[i+1]<<std::endl;
      //std::cout<<"Value of y[i]: "<<y[i]<<"Value of y[i-1]: "<<y[i-1]<<"Value of y[i+1]: "<<y[i+1]<<std::endl;
      //std::cout<<"Value of gradient_diff: "<<gradient_diff<<std::endl;
      zipped_gradient.push_back(std::make_pair(gradient_diff,i));
    }
    
    std::sort(zipped_gradient.begin(), zipped_gradient.end(),greater());
    double gradient_value = zipped_gradient[0].first;
    int iterations_value = zipped_gradient[0].second + 1;

    std::cout<<"Value of FINAL distance value at iteration "<<iterations_value<<" : "<<gradient_value<<std::endl; 
    
  } //loop over the type

  histogram_vec.clear();
  
} //end of iteration test
