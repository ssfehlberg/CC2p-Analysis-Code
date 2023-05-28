#define iteration_other_cxx
#include "plotting.h"
#include "tools/Covariance_Matrix.h"
#include "tools/histograms.h"
#include "tools/xsec.h"
using namespace Histograms;

class iteration_other{

public:
  virtual void main();

private:
  virtual void Grab_Histograms();
  virtual TH1D* make_plot(TH1D* h_one_iter, TH1D* h_iter, int iter, const char* directory, const char* title,const char* variable);
  virtual void Make_Non_Fractional(TH1D* h_Xsec,TH1D* h_systematic);
  virtual TH1D* total_error(const char* name,std::vector<TH1D*> hist);
  virtual void total_covariance(const char* variable, const char* title,std::vector<TH2D*> hist_2D, const char* directory);

  //Plotting class definition
  /////////////////
  plotting_tools plot;
  xsec Xsec;

  //Files we have to evaluate
  static const int num_files = 17;
  const char* directory_name;
  const char* directory_name_list[num_files] = {"detVar","Dirt","flux_all","reint_all","All_UBGenie",
						"AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
						"NormNCCOH_UBGenie","ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
						"VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","xsr_scc_Fa3_SCC",
						"xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie","RootinoFix_UBGenie"};

  //number of universes
  //all of these samples are unisims i.e. # universes = 1
  static const int universes = 1; 

  //Extra Histograms that we need
  //////////////////////////////
  static const int num_samples = 2;

  //Muon
  TH1D* h_muon[num_files][num_var][num_samples];
  TH2D* h_2D_muon[num_files][num_var];
  TH1D* h_muon_xsec[num_var];
  TH1D* h_muon_stat_error[num_var];
  TH1D* h_muon_stat_fractional_error[num_var];
  TH2D* h_muon_stat_matrix[num_var];
  TH1D* h_muon_diff[num_files][num_var];
  TH1D* h_muon_total[num_var];

  //Leading
  TH1D* h_leading[num_files][num_var][num_samples];
  TH2D* h_2D_leading[num_files][num_var];
  TH1D* h_leading_xsec[num_var];
  TH1D* h_leading_stat_error[num_var];
  TH1D* h_leading_stat_fractional_error[num_var];
  TH2D* h_leading_stat_matrix[num_var];
  TH1D* h_leading_diff[num_files][num_var];
  TH1D* h_leading_total[num_var];


  //recoil
  TH1D* h_recoil[num_files][num_var][num_samples];
  TH2D* h_2D_recoil[num_files][num_var];
  TH1D* h_recoil_xsec[num_var];
  TH1D* h_recoil_stat_error[num_var];
  TH1D* h_recoil_stat_fractional_error[num_var];
  TH2D* h_recoil_stat_matrix[num_var];
  TH1D* h_recoil_diff[num_files][num_var];
  TH1D* h_recoil_total[num_var];
  
  //Other
  TH1D* h_other[num_files][num_other_var][num_samples];
  TH2D* h_2D_other[num_files][num_other_var];
  TH1D* h_other_xsec[num_other_var];
  TH1D* h_other_stat_error[num_other_var];
  TH1D* h_other_stat_fractional_error[num_other_var];
  TH2D* h_other_stat_matrix[num_other_var];
  TH1D* h_other_diff[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];

  TFile* file_one_iter[num_files];
  TFile* file_iter[num_files];
  
}; //end of class definition

void iteration_other::main(){

  Grab_Histograms();
  TFile* tfile_mine = new TFile(Form("root_files/Iteration/systematics.root"),"RECREATE"); //output root file                                                                                        

  //particle first
  ////////////////////////////////
  for(int j=0; j < num_var; j++){

    //Muon
    std::vector<TH1D*> muon_vec;
    std::vector<TH2D*> muon_matrices;
    h_muon_xsec[j] =  Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_overlay_muon[j][0],h_overlay_muon[j][1],h_bnb_muon[j],
                                         h_muon_matrices[j],muon_iter[j],h_muon_num[j], h_muon_denom[j],
                                         muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_CV",var0[j]),Form("Iteration"),false,true); //CV Montecarlo values

    //Leading proton
    std::vector<TH1D*> leading_vec;
    std::vector<TH2D*> leading_matrices;
    h_leading_xsec[j] =  Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_overlay_leading[j][0],h_overlay_leading[j][1],h_bnb_leading[j],
                                         h_leading_matrices[j],leading_iter[j],h_leading_num[j], h_leading_denom[j],
                                         leading_xsec_max[j],Form("True Leading Proton %s",var_titles[j]),Form("_leading%s_CV",var0[j]),Form("Iteration"),false,true); //CV Montecarlo values

    //recoil proton
    std::vector<TH1D*> recoil_vec;
    std::vector<TH2D*> recoil_matrices;
    h_recoil_xsec[j] =  Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_overlay_recoil[j][0],h_overlay_recoil[j][1],h_bnb_recoil[j],
                                         h_recoil_matrices[j],recoil_iter[j],h_recoil_num[j], h_recoil_denom[j],
                                         recoil_xsec_max[j],Form("True Recoil Proton %s",var_titles[j]),Form("_recoil%s_CV",var0[j]),Form("Iteration"),false,true); //CV Montecarlo values


    //Loop over the number of Files
    //////////////////////////////
    for(int f=0; f< num_files; f++){
      directory_name = directory_name_list[f];

      //Muon
      /////////////////
      
      //difference plot
      h_muon_diff[f][j] = make_plot(h_muon[f][j][0],h_muon[f][j][1], muon_iter[j],directory_name,Form("True Muon %s",var_titles[j]),Form("_muon%s",var0[j])); //makes the two diff plots                 
      muon_vec.push_back(h_muon_diff[f][j]); //push back diff plot to make total fractional error plot

      //covariance matrrix
      std::vector<TH1D*> muon_universes;
      Make_Non_Fractional(h_muon_xsec[j],h_muon[f][j][0]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      Make_Non_Fractional(h_muon_xsec[j],h_muon[f][j][1]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      muon_universes.push_back(h_muon[f][j][1]); //This has one iteration
      h_2D_muon[f][j] = make_covariance_matrix(h_muon[f][j][0], muon_universes, universes, Form("muon_%s_%s",var0[j],directory_name)); //make the covariance matrix for this particular file
      muon_matrices.push_back(h_2D_muon[f][j]); //vector of all the covariance matrices
      plot.plot_covariance(Form("_muon%s_%s",var0[j],directory_name), Form("Muon %s",var_titles[j]),canv_2D_muon[j],h_2D_muon[f][j],h_muon[f][j][0],Form("Iteration/%s",directory_name)); //plot that matrix
      muon_universes.clear();

      //Leading
      /////////////////
      
      //difference plot
      h_leading_diff[f][j] = make_plot(h_leading[f][j][0],h_leading[f][j][1], leading_iter[j],directory_name,Form("True Leading Proton%s",var_titles[j]),Form("_leading%s",var0[j])); //makes the two diff plots                 
      leading_vec.push_back(h_leading_diff[f][j]); //push back diff plot to make total fractional error plot

      //covariance matrrix
      std::vector<TH1D*> leading_universes;
      Make_Non_Fractional(h_leading_xsec[j],h_leading[f][j][0]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      Make_Non_Fractional(h_leading_xsec[j],h_leading[f][j][1]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      leading_universes.push_back(h_leading[f][j][1]); //This has one iteration
      h_2D_leading[f][j] = make_covariance_matrix(h_leading[f][j][0], leading_universes, universes, Form("leading_%s_%s",var0[j],directory_name)); //make the covariance matrix for this particular file
      leading_matrices.push_back(h_2D_leading[f][j]); //vector of all the covariance matrices
      plot.plot_covariance(Form("_leading%s_%s",var0[j],directory_name), Form("Leading Proton %s",var_titles[j]),canv_2D_leading[j],h_2D_leading[f][j],h_leading[f][j][0],Form("Iteration/%s",directory_name)); //plot that matrix
      leading_universes.clear();

      //Recoil
      /////////////////
      
      //difference plot
      h_recoil_diff[f][j] = make_plot(h_recoil[f][j][0],h_recoil[f][j][1], recoil_iter[j],directory_name,Form("True Recoil Proton%s",var_titles[j]),Form("_recoil%s",var0[j])); //makes the two diff plots                 
      recoil_vec.push_back(h_recoil_diff[f][j]); //push back diff plot to make total fractional error plot

      //covariance matrrix
      std::vector<TH1D*> recoil_universes;
      Make_Non_Fractional(h_recoil_xsec[j],h_recoil[f][j][0]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      Make_Non_Fractional(h_recoil_xsec[j],h_recoil[f][j][1]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      recoil_universes.push_back(h_recoil[f][j][1]); //This has one iteration
      h_2D_recoil[f][j] = make_covariance_matrix(h_recoil[f][j][0], recoil_universes, universes, Form("recoil_%s_%s",var0[j],directory_name)); //make the covariance matrix for this particular file
      recoil_matrices.push_back(h_2D_recoil[f][j]); //vector of all the covariance matrices
      plot.plot_covariance(Form("_recoil%s_%s",var0[j],directory_name), Form("Recoil Proton %s",var_titles[j]),canv_2D_recoil[j],h_2D_recoil[f][j],h_recoil[f][j][0],Form("Iteration/%s",directory_name)); //plot that matrix
      recoil_universes.clear();
      
    }//end of loop over num_files                                                                                                                                                                        

    //Total histograms
    ///////////////////
    h_muon_total[j] = total_error(Form("_muon%s",var0[j]),muon_vec);
    h_muon_total[j]->Write(Form("hist_fractional_errors_muon%s",var0[j]));
    total_covariance(Form("muon_%s",var0[j]), Form("Muon %s",var_titles[j]),muon_matrices,"Iteration");
    muon_vec.clear();
    muon_matrices.clear();

    h_leading_total[j] = total_error(Form("_leading%s",var0[j]),leading_vec);
    h_leading_total[j]->Write(Form("hist_fractional_errors_leading%s",var0[j]));
    total_covariance(Form("leading_%s",var0[j]), Form("Leading %s",var_titles[j]),leading_matrices,"Iteration");
    leading_vec.clear();
    leading_matrices.clear();

    h_recoil_total[j] = total_error(Form("_recoil%s",var0[j]),recoil_vec);
    h_recoil_total[j]->Write(Form("hist_fractional_errors_recoil%s",var0[j]));
    total_covariance(Form("recoil_%s",var0[j]), Form("Recoil %s",var_titles[j]),recoil_matrices,"Iteration");
    recoil_vec.clear();
    recoil_matrices.clear();

    
  } //end loop over num of variables

  ////////////////
  //Now for the other variables
  ///////////////////
  for(int j=0; j < num_other_var; j++){

    std::vector<TH1D*> other_vec;
    std::vector<TH2D*> other_matrices;
    h_other_xsec[j] =  Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_overlay_other[j][0],h_overlay_other[j][1],h_bnb_other[j],
					  h_other_matrices[j],other_iter[j],h_other_num[j], h_other_denom[j],
					  other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("_%s_CV",other_var[j]),Form("Iteration"),false,true); //CV Montecarlo values

    //Loop over the number of Files
    //////////////////////////////
    for(int f=0; f< num_files; f++){
      directory_name = directory_name_list[f];

      //difference plot
      h_other_diff[f][j] = make_plot(h_other[f][j][0],h_other[f][j][1], other_iter[j],directory_name,Form("True %s",other_var_titles[j]),Form("_%s",other_var[j])); //makes the two diff plots                 
      other_vec.push_back(h_other_diff[f][j]); //push back diff plot to make total fractional error plot

      //covariance matrrix
      std::vector<TH1D*> other_universes;
      Make_Non_Fractional(h_other_xsec[j],h_other[f][j][0]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      Make_Non_Fractional(h_other_xsec[j],h_other[f][j][1]);//be sure to multiple the histogram by the CV to get it in terms of the cross-section values
      other_universes.push_back(h_other[f][j][1]); //This has one iteration
      h_2D_other[f][j] = make_covariance_matrix(h_other[f][j][0], other_universes, universes, Form("_%s_%s",other_var[j],directory_name)); //make the covariance matrix for this particular file
      other_matrices.push_back(h_2D_other[f][j]); //vector of all the covariance matrices
      plot.plot_covariance(Form("_other%s_%s",other_var[j],directory_name), Form("%s",var_titles[j]),canv_2D_other[j],h_2D_other[f][j],h_other[f][j][0],Form("Iteration/%s",directory_name)); //plot that matrix
      other_universes.clear();
      
    }//end of loop over number of files

    h_other_total[j] = total_error(Form("_%s",other_var[j]),other_vec);
    h_other_total[j]->Write(Form("hist_fractional_errors%s",other_var[j]));
    total_covariance(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),other_matrices,"Iteration");
    other_vec.clear();
    other_matrices.clear();

  } //end of loop over number of other var

    tfile_mine->Close();
   
} //end of main function

  
void iteration_other::Grab_Histograms(){

   for(int f=0; f < num_files; f++){
    file_one_iter[f] =  new TFile(Form("root_files/%s/systematics_one_iter.root",directory_name_list[f]));
    file_iter[f] =  new TFile(Form("root_files/%s/systematics.root",directory_name_list[f]));

    for(int j=0; j < num_var; j++){
      h_muon[f][j][1] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_muon%s",var0[j])); //one iteration unfolding
      h_muon[f][j][0] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_muon%s",var0[j])); //correct number of iterations

      h_leading[f][j][1] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_leading%s",var0[j])); //one iteration unfolding
      h_leading[f][j][0] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_leading%s",var0[j])); //correct number of iterations

      h_recoil[f][j][1] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j])); //one iteration unfolding
      h_recoil[f][j][0] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors_recoil%s",var0[j])); //correct number of iterations

    } //end of num_var loop

    for( int j=0; j <num_other_var; j++){
      h_other[f][j][1] = (TH1D*)file_one_iter[f]->Get(Form("hist_fractional_errors%s",other_var[j])); //one iteration unfolding
      h_other[f][j][0] = (TH1D*)file_iter[f]->Get(Form("hist_fractional_errors%s",other_var[j])); //correct number of iterations
    } //end of num_other_var

   } //end of num_files loop
   
  
  //BNB,EXT,Dirt, and Overlay Products
  ////////////////////////////////
  TFile* f_bnb = new TFile(Form("../../root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root"));
  TFile* f_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* f_dirt = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root");
  TFile* f_overlay = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root");
  TFile* f_matrices = new TFile("../../root_files/pelee/Run_all/histograms_mc_eff.root");
  TFile* f_eff = new TFile("../../root_files/pelee/Run_all/histograms_mc_eff.root");
  
  for(int i=0; i < num_var; i++){
    
    //BNB, EXT, Dirt Data Products
    h_bnb_muon[i] = (TH1D*)f_bnb->Get(Form("h_muon%s_bnb",var0[i])); //bnb
    h_ext_muon[i] = (TH1D*)f_ext->Get(Form("h_muon%s_ext",var0[i])); //ext
    h_dirt_muon[i] = (TH1D*)f_dirt->Get(Form("h_muon%s_dirt_wgt",var0[i])); //dirt

    h_bnb_leading[i] = (TH1D*)f_bnb->Get(Form("h_leading%s_bnb",var0[i])); //bnb
    h_ext_leading[i] = (TH1D*)f_ext->Get(Form("h_leading%s_ext",var0[i])); //ext
    h_dirt_leading[i] = (TH1D*)f_dirt->Get(Form("h_leading%s_dirt_wgt",var0[i])); //dirt  

    h_bnb_recoil[i] = (TH1D*)f_bnb->Get(Form("h_recoil%s_bnb",var0[i])); //bnb
    h_ext_recoil[i] = (TH1D*)f_ext->Get(Form("h_recoil%s_ext",var0[i])); //ext
    h_dirt_recoil[i] = (TH1D*)f_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[i])); //dirt   

    //Overlay Products
    for(int j=0; j < num_channels; j++){
      h_overlay_muon[i][j] = (TH1D*)f_overlay->Get(Form("h_muon%s%s",var0[i],channel[j])); //overlay products
      h_overlay_leading[i][j] = (TH1D*)f_overlay->Get(Form("h_leading%s%s",var0[i],channel[j])); //overlay products
      h_overlay_recoil[i][j] = (TH1D*)f_overlay->Get(Form("h_recoil%s%s",var0[i],channel[j])); //overlay products 
    }

    //Matrices, Numerator, and Denominator
    h_muon_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_muon_all%s",(var[i]))); //smearing matrix
    h_muon_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_muon_all%s",var[i])); //numerator of efficiency
    h_muon_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_muon_all%s",var[i])); //denom of efficiency

    h_leading_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_lead_proton%s",(var[i]))); //smearing matrix
    h_leading_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_lead_proton%s",var[i])); //numerator of efficiency
    h_leading_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_lead_proton%s",var[i])); //denom of efficiency

    h_recoil_matrices[i] = (TH2D*)f_matrices->Get(Form("h_particle_matrices_recoil_proton%s",(var[i]))); //smearing matrix
    h_recoil_num[i] = (TH1D*)f_eff->Get(Form("h_particle_num_recoil_proton%s",var[i])); //numerator of efficiency
    h_recoil_denom[i] = (TH1D*)f_eff->Get(Form("h_particle_denom_recoil_proton%s",var[i])); //denom of efficiency  

  } //end of loop over num_var

  for(int i=0; i < num_other_var; i++){

    //BNB, EXT, Dirt Data Products
    h_bnb_other[i] = (TH1D*)f_bnb->Get(Form("h%s_bnb",other_var[i])); //bnb
    h_ext_other[i] = (TH1D*)f_ext->Get(Form("h%s_ext",other_var[i])); //ext
    h_dirt_other[i] = (TH1D*)f_dirt->Get(Form("h%s_dirt_wgt",other_var[i])); //dirt

    //Overlay Products
    for(int j=0; j < num_channels; j++){
      h_overlay_other[i][j] = (TH1D*)f_overlay->Get(Form("h%s%s",other_var[i],channel[j])); //overlay products
    }
    //Matrices, Numerator, and Denominator
    h_other_matrices[i] = (TH2D*)f_matrices->Get(Form("h_other_matrices%s",other_var[i])); //smearing matrix
    h_other_num[i] = (TH1D*)f_eff->Get(Form("h_other_eff_num%s",other_var[i])); //numerator of efficiency
    h_other_denom[i] = (TH1D*)f_eff->Get(Form("h_other_eff_denom%s",other_var[i])); //denom of efficiency
    
  } //end of other loop
   
  
} //end of grab histograms


TH1D* iteration_other::make_plot(TH1D* h_one_iter, TH1D* h_iter, int iter, const char* directory, const char* title,const char* variable){

  TH1D* h_diff = (TH1D*)h_iter->Clone();
  TH1D* h_one_iter0 = (TH1D*)h_one_iter->Clone();

  TCanvas* canv = new TCanvas("canv","canv",2000,1500);
  h_one_iter->Draw("hist text");
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
  h_diff->Draw("HIST text");
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

void iteration_other::Make_Non_Fractional(TH1D* h_Xsec,TH1D* h_systematic){

  double nbins = h_Xsec->GetNbinsX();
  for(int x = 1; x < nbins+1; x++){
    double Xsec = h_Xsec->GetBinContent(x);
    double systematic = h_systematic->GetBinContent(x);
    double value = Xsec * systematic;
    h_systematic->SetBinContent(x,value);
  }//end loop over x bins
  
} //end of make non fractional

TH1D* iteration_other::total_error(const char* name,std::vector<TH1D*> hist){

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
  h_total->Draw("HIST TEXT");
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

void iteration_other::total_covariance(const char* variable, const char* title,std::vector<TH2D*> hist_2D, const char* directory){

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size(); i++){
    hist->Add(hist_2D[i]);
  }
  hist->Draw("colz text");
  hist->SetTitle(Form("%s",title)); //title                                                                                                                                                              
  hist->SetXTitle("Reco. Bin a");
  hist->GetXaxis()->SetTitleSize(40);
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleOffset(1.5);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(30);
  hist->GetXaxis()->SetTickSize(0);
  hist->SetYTitle("Reco. Bin b");
  hist->GetYaxis()->SetTitleSize(40);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleOffset(1.5);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(30);
  hist->GetYaxis()->SetTickSize(0);
  hist->Write(Form("h_2D_covariance%s",variable));
  gStyle->SetOptStat(0);
  canv->Print(Form("images/%s/Covariance_Matrices/_2D_%s.png",directory,variable));

}

