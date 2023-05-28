///////////////////////////////////////////////////////////////////////////////////////////////////////////
//July 11th, 2022
//Author: Samantha Sword-Fehlberg
//Purpose of this code is to produce the total correlation matrices (statistical+systematic)
//for all of the variables of interest. Systematic sources include tota GENIE uncertainty (as calculated in plot_all class (plot_GENIE_total)
//detector variations (produced in detector.C), Flux Multisims (), G4 multisim, and dirt systematic
//Total correlation matrices are saved as TH2D inside of ../root_files/total_correlation_matrix.root
//Plot images are saved inside of ../images/Total_Error/Correlation_Matrices
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class total_correlation_matrix{

 public:
  virtual void main();
  
 private:
  
  virtual void Grab_Histograms();
  virtual TH2D* calc_correlation_mat(TH2D* h_covariance);
  virtual void plot_correlation_matrix(const char* variable, const char* title, TH2D* hist, const char* directory);
  virtual void plot_total_correlation_matrix(const char* variable, const char* title,std::vector<TH2D*> hist_2D);

  const char* ztitle; //to make the z axis title
  
  //File names
  ///////////////////////////////////////
  static const int num_files = 7;//,6,7;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt","Iteration","Statistical"};
  const char* file_titles[num_files] = {"GENIE Total","Detector Variations","Flux Multisims","G4 Multisims","Dirt (Sys.)","Iteration","Statistical Error"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH2D* h_muon_covariance[num_files][num_var];
  TH2D* h_muon_correlation[num_files][num_var];
  bool flip_legend_muon[num_var] = {true,false,false};

  //Leading
  TH2D* h_leading_covariance[num_files][num_var];
  TH2D* h_leading_correlation[num_files][num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH2D* h_recoil_covariance[num_files][num_var];
  TH2D* h_recoil_correlation[num_files][num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};


  //Other plots
  ////////////////////
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{#vec{p}_{L}, #vec{p}_{R}})","cos(#gamma_{#vec{p}_{#mu}, #vec{p}_{sum}})",
						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  
  TH2D* h_other_covariance[num_files][num_other_var];
  TH2D* h_other_correlation[num_files][num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};

}; //end of total_correlation_matrix class definition

void total_correlation_matrix::main(){

  Grab_Histograms();

  //File where I will save the TH2Ds
  ///////////////////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/Total_Error/total_correlation_matrices.root"),"RECREATE"); //output root file 
  
  for(int j=0; j < num_var; j++){
    
    std::vector<TH2D*> muon_2D;
    std::vector<TH2D*> leading_2D;
    std::vector<TH2D*> recoil_2D;

    
    for(int f=0; f < num_files; f++){
      
      h_muon_correlation[f][j] = calc_correlation_mat(h_muon_covariance[f][j]);
      h_muon_correlation[f][j]->Write(Form("h_2D_correlation_%s_muon%s",file_titles[f],var0[j]));
      plot_correlation_matrix(Form("_muon%s",var0[j]),Form("%s: Muon %s",file_titles[f],var_titles[j]),h_muon_covariance[f][j],file_names[f]);
      muon_2D.push_back(h_muon_correlation[f][j]);

      h_leading_correlation[f][j] = calc_correlation_mat(h_leading_covariance[f][j]);
      h_leading_correlation[f][j]->Write(Form("h_2D_correlation_%s_leading%s",file_titles[f],var0[j]));
      plot_correlation_matrix(Form("_leading%s",var0[j]),Form("%s: Leading Proton %s",file_titles[f],var_titles[j]),h_leading_covariance[f][j],file_names[f]);
      leading_2D.push_back(h_leading_correlation[f][j]);


      h_recoil_correlation[f][j] = calc_correlation_mat(h_recoil_covariance[f][j]);
      h_recoil_correlation[f][j]->Write(Form("h_2D_correlation_%s_recoil%s",file_titles[f],var0[j]));
      plot_correlation_matrix(Form("_recoil%s",var0[j]),Form("%s: Recoil Proton %s",file_titles[f],var_titles[j]),h_recoil_covariance[f][j],file_names[f]);
      recoil_2D.push_back(h_recoil_correlation[f][j]);

    }

    plot_total_correlation_matrix(Form("_muon%s",var0[j]),Form("Total Correlation Matrix: Muon %s",var_titles[j]),muon_2D);
    plot_total_correlation_matrix(Form("_leading%s",var0[j]),Form("Total Correlation Matrix: Leading Proton %s",var_titles[j]),leading_2D);
    plot_total_correlation_matrix(Form("_recoil%s",var0[j]),Form("Total Correlation Matrix: Recoil Proton %s",var_titles[j]),recoil_2D);
  }

  for(int j=0; j < num_other_var; j++){

    std::vector<TH2D*> other_2D;

    for(int f=0; f < num_files; f++){
     h_other_correlation[f][j] = calc_correlation_mat(h_other_covariance[f][j]);
     h_other_correlation[f][j]->Write(Form("h_2D_correlation_%s%s",file_titles[f],other_var[j]));
     plot_correlation_matrix(Form("%s",other_var[j]),Form("%s: %s",file_titles[f],other_var_titles[j]),h_other_covariance[f][j],file_names[f]);
     other_2D.push_back(h_other_correlation[f][j]);
    }
  
    plot_total_correlation_matrix(Form("%s",other_var[j]),Form("Total Correlation Matrix: %s",other_var_titles[j]),other_2D);
  
    }
  
  tfile_mine->Close();
  
}//end of main


/////////////////////////////////////////////////////////////////
//Grab_Histograms()
//Grabs correlation matrices for all variables from specific files
/////////////////////////////////////////////////////////////////
void total_correlation_matrix::Grab_Histograms(){

  for(int f=0; f < num_files; f++){ 
    file[f] =  new TFile(Form("../root_files/%s/systematics.root",file_names[f]));
    
    for(int j=0; j < num_var; j++){
      h_muon_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancemuon_%s",var0[j]));
      h_leading_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covarianceleading_%s",var0[j]));
      h_recoil_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariancerecoil_%s",var0[j]));
    }
    for(int j=0; j< num_other_var; j++){
      h_other_covariance[f][j] = (TH2D*)file[f]->Get(Form("h_2D_covariance%s",other_var[j]));
    }
  } //end of loop over numbers of files
  
}//end of grab histograms

TH2D* total_correlation_matrix::calc_correlation_mat(TH2D* h_covariance){

TH2D* h_correlation = (TH2D*)h_covariance->Clone();
  int nbinsX = h_correlation->GetNbinsX();
  int nbinsY = h_correlation->GetNbinsY();
  
  for(int i=1; i< nbinsX+1; i++){
    for(int j=1; j < nbinsY+1; j++){

      double binValue = h_covariance->GetBinContent(i,j);
      double iBinValue =  h_covariance->GetBinContent(i,i);
      double jBinValue =  h_covariance->GetBinContent(j,j);
      double value = binValue / (iBinValue * jBinValue);
      h_correlation->SetBinContent(i,j,value);
        
    }
  }

  return h_correlation;

} //end of calc_correlation_mat

void total_correlation_matrix::plot_correlation_matrix(const char* variable, const char* title, TH2D* hist, const char* directory){


  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  hist->Draw("colz text");
  
  hist->SetTitle(Form("Correlation Matrix: %s",title)); //title
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

  int nbins = hist->GetNbinsX();   
  for(int i=1; i < nbins+1; i++){
    hist->GetXaxis()->SetBinLabel(i,Form("%d",i));
    hist->GetYaxis()->SetBinLabel(i,Form("%d",i));
  }
  
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist->GetXaxis()->GetBinLowEdge(i+1);
  }

  TLine* lx;
  TLine* ly;

  for (int bin = 0; bin < nbins+1; bin++) {
    lx = new TLine(edges[bin],edges[0],edges[bin],edges[nbins]);
    lx->Draw("SAME");
    lx->SetLineColor(kBlack);
    lx->SetLineWidth(2);
    ly = new TLine(edges[0],edges[bin],edges[nbins],edges[bin]);
    ly->Draw("SAME");
    ly->SetLineColor(kBlack);
    ly->SetLineWidth(2);  
  }
  gStyle->SetOptStat(0);
  canv->Print(Form("../images/%s/Correlation_Matrices/_2D%s.png",directory,variable));
  
} //end of plot_correlation_matrix



/////////////////////////////////////////
//plot_total_correlation_matrix():
//makes the actual total correlation matrix
//variable: variable we are plotting
//title: title for the plot
//hist_2D: vector of all the correlation matrices listed in file_names
//directory: where to save the image
/////////////////////////////////////////////
void total_correlation_matrix::plot_total_correlation_matrix(const char* variable, const char* title, std::vector<TH2D*> hist_2D){

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  canv->SetRightMargin(0.15);
  //canv->SetLeftMargin(0.1);
  canv->SetBottomMargin(0.1);
  canv->SetTopMargin(0.1);
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size(); i++){
    hist->Add(hist_2D[i]);
  }

  hist->Scale(1E76);
  hist->SetTitle(Form("%s",title));//title
  int nbins = hist->GetNbinsX();
  
  for(int a = 1; a < nbins+1; a++){
    double variance = hist->GetBinContent(a,a);
    
    for(int b=1; b < nbins+1; b++){
      double value = hist->GetBinContent(a,b);
      hist->SetBinContent(a,b,value/variance);
    }
  }
    
  //X title
  hist->GetXaxis()->SetTitle("Bin Number");
  hist->GetXaxis()->SetTitleSize(50);
  hist->GetXaxis()->SetTitleFont(43);
  hist->GetXaxis()->SetTitleOffset(1.3);
  hist->GetXaxis()->SetLabelFont(43);
  hist->GetXaxis()->SetLabelSize(50);
  hist->GetXaxis()->SetTickSize(0);

  //Y Title
  hist->GetYaxis()->SetTitle("Bin Number");
  hist->GetYaxis()->SetTitleSize(50);
  hist->GetYaxis()->SetTitleFont(43);
  hist->GetYaxis()->SetTitleOffset(1.3);
  hist->GetYaxis()->SetLabelFont(43);
  hist->GetYaxis()->SetLabelSize(50);
  hist->GetYaxis()->SetTickSize(0);

  //Z title
  /*if(std::strcmp(variable,"_muon_mom") ==  0 || std::strcmp(variable,"_leading_mom") ==  0 || std::strcmp(variable,"_recoil_mom") ==  0 || std::strcmp(variable,"_delta_PT") == 0){ ztitle = "/(GeV/c)";}
    else{ ztitle ="";}*/
  ztitle = "";
  hist->GetZaxis()->SetTitle(Form("Correlation Element [10^{-76}%s]",ztitle));
  hist->GetZaxis()->CenterTitle(true);
  hist->GetZaxis()->SetTitleSize(50);
  hist->GetZaxis()->SetTitleFont(43);
  hist->GetZaxis()->SetTitleOffset(1.3);
  hist->GetZaxis()->SetLabelFont(43);
  hist->GetZaxis()->SetLabelSize(50);
  hist->GetZaxis()->SetMaxDigits(2);
  hist->GetZaxis()->SetTickSize(0);
  //hist->GetZaxis()->SetRangeUser(hist->GetMinimum(),hist->GetMaximum());
  hist->GetZaxis()->SetRangeUser(-1.0 * hist->GetMaximum(),hist->GetMaximum());
  //hist->GetZaxis()->SetNdivisions(-5);
  
  hist->Draw("COLZ TEXT L");
  gStyle->SetPaintTextFormat("4.3f");
  gStyle->SetPalette(kThermometer);
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");

  //the following lines moe the paletter. Choose the values you need for the position.
   palette->SetX1NDC(0.86);
   palette->SetX2NDC(0.89);
   palette->SetY1NDC(0.1);
   palette->SetY2NDC(0.9);
   gPad->Modified();
   gPad->Update();

  hist->Write();
  canv->Print(Form("../images/Total_Error/Correlation_Matrices/total_correlation_matrix%s.png",variable));
  
} //end of plot_total_correlation


