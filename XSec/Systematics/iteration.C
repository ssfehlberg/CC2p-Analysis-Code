#define iteration_cxx
#include "tools/histograms.h"
using namespace Histograms;

class iteration{

 public:
  virtual void Grab_Iterations();
  virtual void Grab_Histograms();
  virtual TH1D* make_plot(TH1D* h_one_iter, TH1D* h_iter, int iter, const char* directory, const char* title,const char* variable);
  virtual TH1D* total_error(const char* name,std::vector<TH1D*> hist);
  virtual void main();
    
 private:

  static const int num_files = 17;
  const char* directory_name;
  const char* directory_name_list[num_files] = {"detVar","Dirt","flux_all","reint_all","All_UBGenie",
						"AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
						"NormNCCOH_UBGenie","ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
						"VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","xsr_scc_Fa3_SCC",
						"xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie","RootinoFix_UBGenie"};//,"detVar"};
  //Histograams
  //////////////////x
  static const int num_samples = 2;
 
  //Muon
  TH1D* h_muon[num_files][num_var][num_samples];
  TH1D* h_muon_diff[num_files][num_var];
  TH1D* h_muon_total[num_var];
  int muon_iter[num_var];

  //Leading
  TH1D* h_leading[num_files][num_var][num_samples];
  TH1D* h_leading_diff[num_files][num_var];
  TH1D* h_leading_total[num_var];
  int leading_iter[num_var];

  //recoil
  TH1D* h_recoil[num_files][num_var][num_samples];
  TH1D* h_recoil_diff[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  int recoil_iter[num_var];

  
  /* static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
  "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};*/
  
  //Other
  TH1D* h_other[num_files][num_other_var][num_samples];
  TH1D* h_other_diff[num_files][num_var];
  TH1D* h_other_total[num_var];
  int other_iter[num_var];


  TFile* file_one_iter[num_files];
  TFile* file_iter[num_files];

};


void iteration::Grab_Iterations(){

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

void iteration::Grab_Histograms(){

  for(int f=0; f < num_files; f++){
    file_one_iter[f] =  new TFile(Form("root_files/%s/systematics_one_iter.root",directory_name_list[f]));
    file_iter[f] =  new TFile(Form("root_files/%s/systematics.root",directory_name_list[f]));

    for(int j=0; j < num_var; j++){
      h_muon[f][j][0] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_muon%s",var0[j])); //one iteration unfolding
      h_muon[f][j][1] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_muon%s",var0[j])); //correct number of iterations

      h_leading[f][j][0] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_leading%s",var0[j])); //one iteration unfolding
      h_leading[f][j][1] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_leading%s",var0[j])); //correct number of iterations

      h_recoil[f][j][0] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j])); //one iteration unfolding
      h_recoil[f][j][1] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j])); //correct number of iterations

    } //end of num_var loop

    for( int j=0; j <num_other_var; j++){
      h_other[f][j][0] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors%s",other_var[j])); //one iteration unfolding
      h_other[f][j][1] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors%s",other_var[j])); //correct number of iterations
 } //end of num_other_var

  } //end of num_files loop
} //end of Grab_Histograms


TH1D* iteration::make_plot(TH1D* h_one_iter, TH1D* h_iter, int iter, const char* directory, const char* title,const char* variable){

  TH1D* h_diff = (TH1D*)h_iter->Clone();
  TH1D* h_one_iter0 = (TH1D*)h_one_iter->Clone();

  TCanvas* canv = new TCanvas("canv","canv",2000,1500);
  h_one_iter->Draw("hist");
  h_one_iter->SetLineColor(kRed);
  h_one_iter->SetLineWidth(3);
  h_one_iter->SetLineStyle(0);
  h_one_iter->SetTitle(Form("1 Iteration Unfolding vs. %d Iterations Unfolding: %s",iter, variable));
  h_one_iter->GetYaxis()->SetTitle("Fractional Uncertainty");
  
  h_iter->Draw("same hist");
  h_iter->SetLineColor(kBlue);
  h_iter->SetLineStyle(0);
  h_iter->SetLineWidth(3);

  int nbins = h_iter->GetNbinsX();
  double xlow = h_iter->GetBinLowEdge(1);
  double xhigh = h_iter->GetBinLowEdge(nbins+1);

  TLine* a0 = new TLine(xlow,0,xhigh,0);
  a0->Draw("same");
  a0->SetLineColor(kBlack);
  a0->SetLineWidth(2);
  
  TLegend* legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  legend->AddEntry(h_one_iter,"1 Iteration Unfolding","L");
  legend->AddEntry(h_iter,Form("%d Iterations Unfolding",iter),"L");
  legend->SetLineWidth(0);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.03);
  legend->Draw("same");

  canv->Update();
  canv->Print(Form("images/%s/Iteration/%s.png",directory,variable));
  
  TCanvas* canv_diff = new TCanvas("canv_diff","canv_diff",2000,1500);
  h_diff->Add(h_one_iter0, -1.0);
  h_diff->Draw("HIST");
  h_diff->SetLineColor(kBlack);
  h_diff->SetLineStyle(0);
  h_diff->SetLineWidth(3);
  h_diff->GetYaxis()->SetTitle("Fractional Uncertainty");

  TLine* a = new TLine(xlow,0,xhigh,0);
  a->Draw("same");
  a->SetLineColor(kRed);
  a->SetLineWidth(2);
  
  TLegend* legend_diff = new TLegend(0.48, 0.65, 0.89, 0.89);
  legend_diff->AddEntry(h_one_iter,Form("%d-1 Iteration Unfolding",iter),"L");
  legend_diff->SetLineWidth(0);
  legend_diff->SetFillColor(kWhite);
  legend_diff->SetTextSize(0.03);
  legend_diff->Draw("same");

  canv_diff->Update();
  canv_diff->Print(Form("images/%s/Iteration/%s_diff.png",directory,variable));

  return h_diff;
  
} //end of make plot

TH1D* iteration::total_error(const char* name,std::vector<TH1D*> hist){

  int nbins = hist[0]->GetNbinsX();
  double xlow = hist[0]->GetBinLowEdge(1);
  double xhigh = hist[0]->GetBinLowEdge(nbins+1);
  TH1D* h_total = (TH1D*)hist[0]->Clone();

  for(int i=1; i < nbins+1; i++){
    double sum = 0;

    for(int j=0; j < hist.size(); j++){
      double error = std::pow(hist[j]->GetBinContent(i),2);
      sum += error;
    }
    double total_error = std::sqrt(sum);
    h_total->SetBinContent(i,total_error);
  }

  TCanvas* canv_total = new TCanvas("canv_total","canv_total",2000,1500);
  h_total->Draw("HIST");
  h_total->SetLineColor(kRed);
  h_total->SetLineStyle(0);
  h_total->SetLineWidth(3);
  h_total->GetYaxis()->SetTitle("Fractional Uncertainty");
  h_total->SetMaximum(0.3);
  h_total->SetMinimum(-0.01);

  TLine* a = new TLine(xlow,0,xhigh,0);
  a->Draw("same");
  a->SetLineColor(kBlack);
  a->SetLineWidth(2);

  canv_total->Update();
  canv_total->Print(Form("images/Iteration/%s.png",name));

  return h_total;

} //end of total error

void iteration::main(){

  Grab_Iterations(); //grabs the number of iterations
  Grab_Histograms(); //grabs the histograms

  TFile* tfile_out = new TFile(Form("root_files/Iteration/total_iter_error.root"),"RECREATE"); //output root file
  
  for(int j=0; j < num_var; j++){
    std::vector<TH1D*> muon_vec;
    std::vector<TH1D*> leading_vec;
    std::vector<TH1D*> recoil_vec;

    for(int f=0; f< num_files; f++){
      directory_name = directory_name_list[f];
      h_muon_diff[f][j] = make_plot(h_muon[f][j][0],h_muon[f][j][1], muon_iter[j],directory_name,Form("True Muon %s",var_titles[j]),Form("_muon%s",var0[j])); //makes the two diff plots
      h_leading_diff[f][j] = make_plot(h_leading[f][j][0],h_leading[f][j][1], leading_iter[j],directory_name,Form("True Leading %s",var_titles[j]),Form("_leading%s",var0[j])); //makes the two diff plots
      h_recoil_diff[f][j] = make_plot(h_recoil[f][j][0],h_recoil[f][j][1], recoil_iter[j],directory_name,Form("True Recoil %s",var_titles[j]),Form("_recoil%s",var0[j])); //makes the two diff plots
      muon_vec.push_back(h_muon_diff[f][j]);
      leading_vec.push_back(h_leading_diff[f][j]);
      recoil_vec.push_back(h_recoil_diff[f][j]);
    }//end of loop over num_files

    h_muon_total[j] = total_error(Form("_muon%s",var0[j]),muon_vec);
    h_muon_total[j]->Write(Form("hist_fractional_errors_muon%s",var0[j]));
    muon_vec.clear();

    h_leading_total[j] = total_error(Form("_leading%s",var0[j]),leading_vec);
    h_leading_total[j]->Write(Form("hist_fractional_errors_leading%s",var0[j]));
    leading_vec.clear();

    h_recoil_total[j] = total_error(Form("_recoil%s",var0[j]),recoil_vec);
    h_recoil_total[j]->Write(Form("hist_fractional_errors_recoil%s",var0[j]));
    recoil_vec.clear();
    
  } //end of loop num_var

  for(int j=0; j < num_other_var; j++){
    std::vector<TH1D*> other_vec;
    for(int f=0; f< num_files; f++){
      directory_name = directory_name_list[f];
      h_other_diff[f][j] = make_plot(h_other[f][j][0],h_other[f][j][1], other_iter[j],directory_name,Form("True %s",other_var_titles[j]),Form("_muon%s",other_var[j])); //makes the two diff plots
      other_vec.push_back(h_other_diff[f][j]);
    }//end of loop over num_files
    
  h_other_total[j] = total_error(Form("_%s",other_var[j]),other_vec);
  h_other_total[j]->Write(Form("hist_fractional_errors%s",other_var[j]));
  other_vec.clear();
    
  } //end of loop over num_other_var

  tfile_out->Close();
  
}//end of main
