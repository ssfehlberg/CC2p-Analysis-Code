#define plot_GENIE_all_cxx
#include "../paul_tol_colors.hpp"
#include <iostream>
#include <ctime>
#include <string>

TH1D* total_error(const char* name,std::vector<TH1D*> hist){

  int nbins = hist[0]->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist[0]->GetBinLowEdge(i+1);
  }

  TH1D* h_total_error = new TH1D(Form("hist_fractional_errors%s",name),Form("hist_fractional_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    double sum = 0;

    for(int j=0; j < hist.size(); j++){
      double error = std::pow(hist[j]->GetBinContent(i),2);
      sum += error;
    }

    double total_error = std::sqrt(sum);
    h_total_error->SetBinContent(i,total_error);

  }

  return h_total_error;

}

void plot_GENIE_all(){

  //Stuff for Drawing
  //////////////////////
  TLatex* t = new TLatex(); //latex
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  t->SetNDC();
  t->SetTextAlign(22);
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";
  tolcols::init();
  Color_t colors[] = {kRed, kOrange, kOrange+7, kYellow, kGreen, kTeal,kCyan, kBlue-10,kBlue, kViolet-8, kMagenta, kGray+2};

  static const int num_files = 12;
  TFile* file[num_files];
  const char* file_names[num_files] = {"All_UBGenie","AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
				       "NormNCCOH_UBGenie","RPA_CCQE_UBGenie",
				       "ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
				       "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie",
				       "xsr_scc_Fa3_SCC","xsr_scc_Fv3_SCC"};



  const char* file_titles[num_files] = {"GENIE Multisims", "AxFFCCQE Shape","DecayAngleMEC","NormCCCOH",
					"NormNCCOH","RPA_CCQE",
					"ThetaDelta2NRad","Theta_Delta2Npi",
					"VecFFCCQEshape","XSecShape_CCMEC",
					"SCC: F^{3}_{V}","SCC: F^{3}_{A}"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH1D* h_muon[num_files][num_var];
  TH1D* h_muon_total[num_var];
  TCanvas* canv_muon[num_var];
  TLegend* legend_muon[num_var];

  //Leading
  TH1D* h_leading[num_files][num_var];
  TH1D* h_leading_total[num_var];
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];

  //Recoil
  TH1D* h_recoil[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];

  static const int num_other_var = 7;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)",
                                                 "Estimated p_{n} (GeV/c)","Estimated Neutrino Energy (GeV/c)"};
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  

  //Grab Files
  ////////////////
  for(int f=0; f < num_files; f++){
    file[f] =  new TFile(Form("root_files/%s/systematics.root",file_names[f]));
    for(int j=0; j < num_var; j++){
       h_muon[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_muon%s",var0[j]));
       h_leading[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_leading%s",var0[j]));
       h_recoil[f][j] = (TH1D*)file[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j]));
    }
  }


    //File to save stuff
  ////////////////////
  TFile *tfile_mine = new TFile(Form("root_files/Genie_Total/systematics.root"),"RECREATE"); //output root file 

  
  for(int j=0; j < num_var; j++){

    canv_muon[j] = new TCanvas(Form("canv_muon%s",var0[j]),Form("canv_muon%s",var0[j]),2000,1500);
    legend_muon[j] = new TLegend(0.48, 0.65, 0.89, 0.89);
    legend_muon[j]->SetNColumns(2);
    
    h_muon[0][j]->Draw("hist");
    h_muon[0][j]->SetLineColor(colors[0]);
    h_muon[0][j]->SetLineWidth(2);
    h_muon[0][j]->SetLineStyle(1);
    
    h_muon[0][j]->SetTitle(Form("Muon %s",var_titles[j])); //title
    h_muon[0][j]->SetXTitle(Form("Muon %s",var_titles[j])); //xtitle
    h_muon[0][j]->GetXaxis()->SetTitleSize(40);
    h_muon[0][j]->GetXaxis()->SetTitleFont(43);
    h_muon[0][j]->GetXaxis()->SetTitleOffset(1.5);
    h_muon[0][j]->GetXaxis()->SetLabelFont(43);
    h_muon[0][j]->GetXaxis()->SetLabelSize(30);
    h_muon[0][j]->SetYTitle("Fractional Error"); //Y title
    h_muon[0][j]->GetYaxis()->SetTitleSize(40);
    h_muon[0][j]->GetYaxis()->SetTitleFont(43);
    h_muon[0][j]->GetYaxis()->SetTitleOffset(1.5);
    h_muon[0][j]->GetYaxis()->SetLabelFont(43);
    h_muon[0][j]->GetYaxis()->SetLabelSize(30);
    h_muon[0][j]->SetMaximum(1.0); //max
    h_muon[0][j]->SetMinimum(0); //min  
 
    std::vector<TH1D*> muon;
    muon.push_back(h_muon[0][j]);

    legend_muon[j]->AddEntry(h_muon[0][j],Form("%s",file_titles[0]),"lpf");
    for(int f=1; f <num_files; f++){
      h_muon[f][j]->Draw("hist same");
      h_muon[f][j]->SetLineColor(colors[f]);
      h_muon[f][j]->SetLineWidth(2);
      h_muon[f][j]->SetLineStyle(1);
      legend_muon[j]->AddEntry(h_muon[f][j],Form("%s",file_titles[f]),"lpf");
      muon.push_back(h_muon[f][j]);
    }
    
    h_muon_total[j] = total_error(Form("_muon%s",var0[j]),muon);
    h_muon_total[j]->Draw("hist same text");
    h_muon_total[j]->SetLineColor(kBlack);
    h_muon_total[j]->SetLineWidth(2);
    h_muon_total[j]->SetLineStyle(1);
    h_muon_total[j]->Write();

    legend_muon[j]->AddEntry(h_muon_total[j],"Total","lpf");
    legend_muon[j]->Draw("same");
    gStyle->SetOptStat(0);
    canv_muon[j]->Print(Form("images/GENIE_Total/_total_error_muon%s.png",var0[j]));

  }

  tfile_mine->Close();
  
}
