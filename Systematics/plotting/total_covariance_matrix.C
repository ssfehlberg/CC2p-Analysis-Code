///////////////////////////////////////////////////////////////////////////////////////////////////////////
//July 11th, 2022
//Author: Samantha Sword-Fehlberg
//Purpose of this code is to produce the total covariance matrices (statistical+systematic)
//for all of the variables of interest. Systematic sources include tota GENIE uncertainty (as calculated in plot_all class (plot_GENIE_total)
//detector variations (produced in detector.C), Flux Multisims (), G4 multisim, and dirt systematic
//Total covariance matrices are saved as TH2D inside of ../root_files/total_covariance_matrix.root
//Plot images are saved inside of ../images/Total_Error/Covariance_Matrices
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class total_covariance_matrix{

 public:
  virtual void main();
  
 private:
  
  virtual void Grab_Histograms();
  virtual void plot_total_covariance_matrix(const char* variable, const char* title,std::vector<TH2D*> hist_2D);

  const char* ztitle; //to make the z axis title
  
  //File names
  ///////////////////////////////////////
  static const int num_files = 6;//,6,7;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt","Statistical"};
  const char* file_titles[num_files] = {"GENIE Total","Detector Variations","Flux Multisims","G4 Multisims","Dirt (Sys.)","Statistical Error"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH2D* h_muon_covariance[num_files][num_var];
  bool flip_legend_muon[num_var] = {true,false,false};

  //Leading
  TH2D* h_leading_covariance[num_files][num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH2D* h_recoil_covariance[num_files][num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};


  //Other plots
  ////////////////////
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{#vec{p}_{L}, #vec{p}_{R}})","cos(#gamma_{#vec{p}_{#mu}, #vec{p}_{sum}})",
						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  
  TH2D* h_other_covariance[num_files][num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};

}; //end of total_covariance_matrix class definition

void total_covariance_matrix::main(){

  Grab_Histograms();

  //File where I will save the TH2Ds
  ///////////////////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/Total_Error/total_covariance_matrices.root"),"RECREATE"); //output root file 

  for(int j=0; j < num_var; j++){

    std::vector<TH2D*> muon_2D;
    std::vector<TH2D*> leading_2D;
    std::vector<TH2D*> recoil_2D;
    
    for(int f=0; f < num_files; f++){
      muon_2D.push_back(h_muon_covariance[f][j]);
      leading_2D.push_back(h_leading_covariance[f][j]);
      recoil_2D.push_back(h_recoil_covariance[f][j]); 
    }
    
    plot_total_covariance_matrix(Form("_muon%s",var0[j]),Form("Total Covariance Matrix: Muon %s",var_titles[j]),muon_2D);
    plot_total_covariance_matrix(Form("_leading%s",var0[j]),Form("Total Covariance Matrix: Leading Proton %s",var_titles[j]),leading_2D);
    plot_total_covariance_matrix(Form("_recoil%s",var0[j]),Form("Total Covariance Matrix: Recoil Proton %s",var_titles[j]),recoil_2D);
  }

  for(int j=0; j < num_other_var; j++){

    std::vector<TH2D*> other_2D;

    for(int f=0; f < num_files; f++){
     other_2D.push_back(h_other_covariance[f][j]);
    }
  
    plot_total_covariance_matrix(Form("%s",other_var[j]),Form("Total Covariance Matrix: %s",other_var_titles[j]),other_2D);
  
  }
  
  tfile_mine->Close();
  
}//end of main


/////////////////////////////////////////////////////////////////
//Grab_Histograms()
//Grabs covariance matrices for all variables from specific files
/////////////////////////////////////////////////////////////////
void total_covariance_matrix::Grab_Histograms(){

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

/////////////////////////////////////////
//plot_total_covariance_matrix():
//makes the actual total covariance matrix
//variable: variable we are plotting
//title: title for the plot
//hist_2D: vector of all the covariance matrices listed in file_names
//directory: where to save the image
/////////////////////////////////////////////
void total_covariance_matrix::plot_total_covariance_matrix(const char* variable, const char* title, std::vector<TH2D*> hist_2D){

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  canv->SetRightMargin(0.15);
  //canv->SetLeftMargin(0.1);
  canv->SetBottomMargin(0.1);
  canv->SetTopMargin(0.1);
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size(); i++){
    hist->Add(hist_2D[i]);
  }

  hist->SetTitle(Form("%s",title));//title

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
  hist->GetZaxis()->SetTitle(Form("Covariance Element [10^{-76}%s]",ztitle));
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
  canv->Print(Form("../images/Total_Error/Covariance_Matrices/total_covariance_matrix%s.png",variable));
  
} //end of plot_total_covariance


