//#include "shared.h"


class plot_all{

 public:
  virtual void main();

  
 private:
  virtual void Grab_Histograms();
  virtual void plot_total_error(const char* variable, bool flip_legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total);
  
  //Names of all the systematic files
  ////////////////////////////////////
  static const int num_files = 5;//5,6,7;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt"};//,"Iteration","Statistical"};
  const char* file_titles[num_files] = {"GENIE Total","Detector Variations","Flux Multisims","G4 Multisims","Dirt (Sys.)"};//,"Iteration","Statistical Error"};
  
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
  bool flip_legend_recoil[num_var] = {true,false,false};


  //Other variables
  /////////////////////////////
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  //Other
  TH1D* h_other[num_files][num_other_var];
  TH1D* h_other_total[num_other_var];
  TH2D* h_other_covariance[num_files][num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};


  
}; //end of class definition


//////////////
//[main]
//////////////////
void plot_all::main(){

  /* //Defining all the classes
  //////////////////////////
  detVar detvar;
  GENIE_total genie;

  //plots the detector variations and creates the total histogram
  std::cout<<"[MAKING DETECTOR SYSTEMATICS PLOTS]"<<std::endl;
  detvar.main();
  std::cout<<"[FINISHEDMAKING DETECTOR SYSTEMATICS PLOTS]"<<std::endl;

  //makes the total GENIE plot
  std::cout<<"[MAKING MC SYSTEMATICS PLOTS]"<<std::endl;
  genie.main();
  std::cout<<"[MAKING MC SYSTEMATICS PLOTS]"<<std::endl;*/

  Grab_Histograms();

  TFile *tfile_mine = new TFile(Form("../root_files/Total_Error/total_error.root"),"RECREATE"); //output root file

  //Plotting
  //////////////////
  for(int j=0; j < num_var; j++){
    
    std::vector<TH1D*> muon;
    std::vector<TH2D*> muon_2D;
    std::vector<TH1D*> leading;
    std::vector<TH2D*> leading_2D;
    std::vector<TH1D*> recoil;
    std::vector<TH2D*> recoil_2D;

    for(int f=0; f< num_files; f++){
      muon.push_back(h_muon[f][j]);
      leading.push_back(h_leading[f][j]);
      recoil.push_back(h_recoil[f][j]);
      muon_2D.push_back(h_muon_covariance[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]);

      plot_covariance(Form("_muon%s",var0[j]),Form("%s: Muon %s",file_titles[f],var_titles[j]),h_muon_covariance[f][j],file_names[f]);
      plot_covariance(Form("_leading%s",var0[j]),Form("%s: Leading Proton %s",file_titles[f],var_titles[j]),h_leading_covariance[f][j],file_names[f]);
      plot_covariance(Form("_recoil%s",var0[j]),Form("%s: Recoil Proton %s",file_titles[f],var_titles[j]),h_recoil_covariance[f][j],file_names[f]);
    }

    plot_total_error(Form("_muon%s",var0[j]),flip_legend_muon[j],muon,Form("Muon %s",var_titles[j]),h_muon_total[j]);
    make_total_covariance_matrix(Form("_muon%s",var0[j]),Form("Total Error: Muon %s",var_titles[j]),muon_2D,"Total_Error");

    plot_total_error(Form("_leading%s",var0[j]),flip_legend_leading[j],leading,Form("Leading Proton %s",var_titles[j]),h_leading_total[j]);
    make_total_covariance_matrix(Form("_leading%s",var0[j]),Form("Total Error: Leading Proton %s",var_titles[j]),leading_2D,"Total_Error");
    
    plot_total_error(Form("_recoil%s",var0[j]),flip_legend_recoil[j],recoil,Form("Recoil Proton %s",var_titles[j]),h_recoil_total[j]);
    make_total_covariance_matrix(Form("_recoil%s",var0[j]),Form("Total Error: Recoil Proton %s",var_titles[j]),recoil_2D,"Total_Error");
  }
  
  for(int j=0; j < num_other_var; j++){
    std::vector<TH1D*> other;
    std::vector<TH2D*> other_2D;

    for(int f=0; f < num_files; f++){
      other.push_back(h_other[f][j]);
 
      other_2D.push_back(h_other_covariance[f][j]);
      plot_covariance(Form("%s",other_var[j]),Form("%s: %s",file_titles[f],other_var_titles[j]),h_other_covariance[f][j],file_names[f]);
    }
    
    plot_total_error(Form("%s",other_var[j]),flip_legend_other[j],other,Form("%s",other_var_titles[j]),h_other_total[j]);
    make_total_covariance_matrix(Form("%s",other_var[j]),Form("Total Error: %s",other_var_titles[j]),other_2D,"Total_Error");
  }

  tfile_mine->Close();

} //end of main function


///////////
//[Grab_Histograms]
//Grabs histograms from relevant files
///////////////////
void plot_all::Grab_Histograms(){

  for(int f=0; f < num_files; f++){
    
    if(f == 5){
      file[f] =  new TFile(Form("../root_files/Iteration/total_iter_error.root"));
    }else{
      file[f] =  new TFile(Form("../root_files/%s/systematics.root",file_names[f]));
    }
    
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
} //end of Grab_histograms

///////////
//[plot_total_error]
//Makes the total error plot for the multiple sources of error
//variable: name of variable we are plotting
//flip legend: do you want legend on left or right
//hist: vector of TH1D that represent the fractional error for all the files listed in file_names NOTE: ASSUMES FRACTIONAL ERROR
//title: title of plot and x axis
//h_total: addition in quadrature of all sources in hist.  
//////////////////////////
void plot_all::plot_total_error(const char* variable, bool flip_legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total){

  //Stuff for Drawing
  /////////////////////
  static const int num_files = 6;
  const char* file_titles[num_files] = {"Cross-Section","Detector Variations","Flux","Reinteractions","Dirt","Statistical Error"};
  Color_t colors[] = {kRed, kOrange+6,kYellow-3, kGreen, kBlue,kViolet+1};

  //Now to plot
  /////////////
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
   hist[0]->SetLineStyle(1);
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
     hist[f]->SetLineStyle(1);
     hist[f]->SetMarkerColor(colors[f]);
     legend->AddEntry(hist[f],Form("%s",file_titles[f]),"l");
   }
   
   h_total = total_error(Form("%s",variable),hist);
   h_total->Draw("HIST TEXT00 SAME");
   h_total->SetLineColor(kBlack);
   h_total->SetLineWidth(3);
   h_total->SetLineStyle(1);
   h_total->SetMarkerColor(kBlack);
   h_total->Write();
   
   legend->AddEntry(h_total,"Total Uncertainty","lpf");
   legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, In-Progress","C");
   legend->Draw("same");
   gStyle->SetOptStat(0);
   canv->Print(Form("../images/Total_Error/_total_error%s.png",variable));
  
}//end of plot_total_error
