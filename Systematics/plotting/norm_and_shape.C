///////////////////////////////////////////////////////////////////////////////////////////////////////////
//July 11th, 2022
//Author: Samantha Sword-Fehlberg
//Purpose of this code is to produce the total covariance matrices (statistical+systematic)
//for all of the variables of interest. Systematic sources include tota GENIE uncertainty (as calculated in plot_all class (plot_GENIE_total)
//detector variations (produced in detector.C), Flux Multisims (), G4 multisim, and dirt systematic
//Total covariance matrices are saved as TH2D inside of ../root_files/total_covariance_matrix.root
//Plot images are saved inside of ../images/Total_Error/Covariance_Matrices
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
class norm_and_shape{

 public:
  virtual void main();
  
 private:
  
  virtual void Grab_Histograms();
  virtual TH1D* calculate_prediction(TH1D* h_overlay, TH1D* h_ext, TH1D* h_dirt);
  virtual TMatrixD convert_to_Tmatrix(TH2D* h_cov);
  virtual TVectorD convert_to_Tvector(TH1D* h_input);
  virtual std::vector<TMatrixD> MatrixDecomp(int nbins,TVectorD matrix_pred,TMatrixD matrix_syst);
  //virtual void plot_norm_and_shape(const char* variable, const char* title,std::vector<TH2D*> hist_2D);

  const char* ztitle; //to make the z axis title
  
  //File names
  ///////////////////////////////////////
  static const int num_files = 5;//,6,7;
  TFile* file[num_files];
  const char* file_names[num_files] = {"GENIE_total","detVar","flux_all","reint_all","Dirt"};
  const char* file_titles[num_files] = {"GENIE Total","Detector Variations","Flux Multisimse","G4 Multisims","Dirt (Sys.)"};

  //Individual Particles Plots
  /////////////////////////////
  static const int num_var = 3;
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};

  //Muon
  TH2D* h_muon_covariance[num_var]; //total covariance matrix not includeing statistical
  TH1D* h_muon_bnb[num_var];
  TH1D* h_muon_ext[num_var];
  TH1D* h_muon_dirt[num_var];
  TH1D* h_muon_overlay[num_var];
  TH1D* h_muon_prediction[num_var];
  TH1D* h_muon_bnb_bkg_sub[num_var];
  TH2D h_muon_shape[num_var];
  TH2D h_muon_norm[num_var];
  bool flip_legend_muon[num_var] = {true,false,false};

  //Leading
  TH2D* h_leading_covariance[num_var];  //total covariance matrix not includeing statistical
  TH1D* h_leading_bnb[num_var];
  TH1D* h_leading_ext[num_var];
  TH1D* h_leading_dirt[num_var];
  TH1D* h_leading_overlay[num_var];
  TH1D* h_leading_prediction[num_var];
  TH1D* h_leading_bnb_bkg_sub[num_var];
  TH2D h_leading_shape[num_var];
  TH2D h_leading_norm[num_var];
  bool flip_legend_leading[num_var] = {true,false,false};
 
  //Recoil
  TH2D* h_recoil_covariance[num_var];  //total covariance matrix not includeing statistical
  TH1D* h_recoil_bnb[num_var];
  TH1D* h_recoil_ext[num_var];
  TH1D* h_recoil_dirt[num_var];
  TH1D* h_recoil_overlay[num_var];
  TH1D* h_recoil_prediction[num_var];
  TH1D* h_recoil_bnb_bkg_sub[num_var];
  TH2D h_recoil_shape[num_var];
  TH2D h_recoil_norm[num_var];
  bool flip_legend_recoil[num_var] = {true,false,false};


  //Other plots
  ////////////////////
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{#vec{p}_{L}, #vec{p}_{R}})","cos(#gamma_{#vec{p}_{#mu}, #vec{p}_{sum}})",
						 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};
  
  TH2D* h_other_covariance[num_other_var];
  TH1D* h_other_bnb[num_other_var];
  TH1D* h_other_ext[num_other_var];
  TH1D* h_other_dirt[num_other_var];
  TH1D* h_other_overlay[num_other_var];
  TH1D* h_other_prediction[num_other_var];
  TH1D* h_other_bnb_bkg_sub[num_other_var];
  TH2D h_other_shape[num_other_var];
  TH2D h_other_norm[num_other_var];
  bool flip_legend_other[num_other_var] = {true,false,false,false,false};

}; //end of norm_and_shape class definition

void norm_and_shape::main(){

  Grab_Histograms();

  //File where I will save the TH2Ds
  ///////////////////////////////////
  TFile *tfile_mine = new TFile(Form("../root_files/Total_Error/shape_and_normalization.root"),"RECREATE"); //output root file 

  for(int j=0; j < num_var; j++){

    h_muon_prediction[j] = calculate_prediction(h_muon_overlay[j],h_muon_ext[j],h_muon_dirt[j]);
    h_muon_prediction[j]->Write(Form("h_muon%s_prediction",var0[j]));
    h_leading_prediction[j] = calculate_prediction(h_leading_overlay[j],h_leading_ext[j],h_leading_dirt[j]);
    h_leading_prediction[j]->Write(Form("h_leading%s_prediction",var0[j]));
    h_recoil_prediction[j] = calculate_prediction(h_recoil_overlay[j],h_recoil_ext[j],h_recoil_dirt[j]);
    h_recoil_prediction[j]->Write(Form("h_recoil%s_prediction",var0[j]));
    
    h_muon_bnb_bkg_sub[j] = (TH1D*)h_muon_bnb[j]->Clone();
    h_muon_bnb_bkg_sub[j]->Add(h_muon_prediction[j],-1);
    h_muon_bnb_bkg_sub[j]->Write(Form("h_muon%s_bnb_bkg_sub",var0[j]));
    h_leading_bnb_bkg_sub[j] = (TH1D*)h_leading_bnb[j]->Clone();
    h_leading_bnb_bkg_sub[j]->Add(h_leading_prediction[j],-1);
    h_leading_bnb_bkg_sub[j]->Write(Form("h_leading%s_bnb_bkg_sub",var0[j]));
    h_recoil_bnb_bkg_sub[j] = (TH1D*)h_recoil_bnb[j]->Clone();
    h_recoil_bnb_bkg_sub[j]->Add(h_recoil_prediction[j],-1);
    h_recoil_bnb_bkg_sub[j]->Write(Form("h_recoil%s_bnb_bkg_sub",var0[j]));
    
    int nbins_muon = h_muon_bnb_bkg_sub[j]->GetNbinsX();
    int nbins_leading = h_leading_bnb_bkg_sub[j]->GetNbinsX();
    int nbins_recoil = h_recoil_bnb_bkg_sub[j]->GetNbinsX();
    
    TVectorD muon_bnb = convert_to_Tvector(h_muon_bnb_bkg_sub[j]);
    TVectorD leading_bnb = convert_to_Tvector(h_leading_bnb_bkg_sub[j]);
    TVectorD recoil_bnb = convert_to_Tvector(h_recoil_bnb_bkg_sub[j]);
    
    TMatrixD muon_mat = convert_to_Tmatrix(h_muon_covariance[j]);
    TMatrixD leading_mat = convert_to_Tmatrix(h_leading_covariance[j]);
    TMatrixD recoil_mat = convert_to_Tmatrix(h_recoil_covariance[j]);
 
    std::vector<TMatrixD> muon_result = MatrixDecomp(nbins_muon,muon_bnb,muon_mat);
    h_muon_shape[j] = TH2D(muon_result[0]);
    h_muon_shape[j].Write(Form("h_2D_muon%s_shape",var0[j]));
    h_muon_norm[j] = TH2D(muon_result[1]);
    h_muon_norm[j].Write(Form("h_2D_muon%s_norm",var0[j]));

    std::vector<TMatrixD> leading_result = MatrixDecomp(nbins_leading,leading_bnb,leading_mat);
    h_leading_shape[j] = TH2D(leading_result[0]);
    h_leading_shape[j].Write(Form("h_2D_leading%s_shape",var0[j]));
    h_leading_norm[j] = TH2D(leading_result[1]);
    h_leading_norm[j].Write(Form("h_2D_leading%s_norm",var0[j]));

    std::vector<TMatrixD> recoil_result = MatrixDecomp(nbins_recoil,recoil_bnb,recoil_mat);
    h_recoil_shape[j] = TH2D(recoil_result[0]);
    h_recoil_shape[j].Write(Form("h_2D_recoil%s_shape",var0[j]));
    h_recoil_norm[j] = TH2D(recoil_result[1]);
    h_recoil_norm[j].Write(Form("h_2D_recoil%s_norm",var0[j]));

    //plot_norm_and_shape(Form("_muon%s",var0[j]),Form("Total Covariance Matrix: Muon %s",var_titles[j]),muon_2D);
    //plot_norm_and_shape(Form("_leading%s",var0[j]),Form("Total Covariance Matrix: Leading Proton %s",var_titles[j]),leading_2D);
    //plot_norm_and_shape(Form("_recoil%s",var0[j]),Form("Total Covariance Matrix: Recoil Proton %s",var_titles[j]),recoil_2D);
  }

  for(int j=0; j < num_other_var; j++){

    h_other_prediction[j] = calculate_prediction(h_other_overlay[j],h_other_ext[j],h_other_dirt[j]);
    h_other_prediction[j]->Write(Form("h%s_prediction",other_var[j]));

    h_other_bnb_bkg_sub[j] = (TH1D*)h_other_bnb[j]->Clone();
    h_other_bnb_bkg_sub[j]->Add(h_other_prediction[j],-1);
    h_other_bnb_bkg_sub[j]->Write(Form("h%s_bnb_bkg_sub",other_var[j]));
    
    int nbins_other = h_other_bnb_bkg_sub[j]->GetNbinsX();
    TVectorD other_bnb = convert_to_Tvector(h_other_bnb_bkg_sub[j]);
    TMatrixD other_mat = convert_to_Tmatrix(h_other_covariance[j]);
    std::vector<TMatrixD> other_result = MatrixDecomp(nbins_other,other_bnb,other_mat);
    h_other_shape[j] = TH2D(other_result[0]);
    h_other_shape[j].Write(Form("h_2D%s_shape",other_var[j]));
    h_other_norm[j] = TH2D(other_result[1]);
    h_other_norm[j].Write(Form("h_2D%s_norm",other_var[j]));
    
    //plot_norm_and_shape(Form("%s",other_var[j]),Form("Total Covariance Matrix: %s",other_var_titles[j]),other_2D);
  
  }
  
  tfile_mine->Close();
  
}//end of main


/////////////////////////////////////////////////////////////////
//Grab_Histograms()
//Grabs covariance matrices for all variables from specific files
/////////////////////////////////////////////////////////////////
void norm_and_shape::Grab_Histograms(){

  
  //Grabbing the total covariance matrix (these ones do not include the statistical uncertainty
  ///////////////////////////////////////////////////////////////////////////////////////
  TFile* file = new TFile("../root_files/Total_Error/total_covariance_matrices_no_stat.root");
  TFile* fBNB = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root");
  TFile* fEXT = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_ext.root");
  TFile* fOverlay = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root");
  TFile* fDirt = new TFile("../../root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root");
  
  for(int i=0; i < num_var; i++){

    h_muon_covariance[i] = (TH2D*)file->Get(Form("h_2D_covariancemuon_%s",var0[i]));
    h_muon_bnb[i] = (TH1D*)fBNB->Get(Form("h_muon%s_bnb",var0[i]));
    h_muon_ext[i] = (TH1D*)fEXT->Get(Form("h_muon%s_ext",var0[i]));
    h_muon_overlay[i] = (TH1D*)fOverlay->Get(Form("h_muon%s_total",var0[i]));
    h_muon_dirt[i] = (TH1D*)fDirt->Get(Form("h_muon%s_dirt_wgt",var0[i]));

    h_leading_covariance[i] = (TH2D*)file->Get(Form("h_2D_covarianceleading_%s",var0[i]));
    h_leading_bnb[i] = (TH1D*)fBNB->Get(Form("h_leading%s_bnb",var0[i]));
    h_leading_ext[i] = (TH1D*)fEXT->Get(Form("h_leading%s_ext",var0[i]));
    h_leading_overlay[i] = (TH1D*)fOverlay->Get(Form("h_leading%s_total",var0[i]));
    h_leading_dirt[i] = (TH1D*)fDirt->Get(Form("h_leading%s_dirt_wgt",var0[i]));
    
    h_recoil_covariance[i] = (TH2D*)file->Get(Form("h_2D_covariancerecoil_%s",var0[i]));
    h_recoil_bnb[i] = (TH1D*)fBNB->Get(Form("h_recoil%s_bnb",var0[i]));
    h_recoil_ext[i] = (TH1D*)fEXT->Get(Form("h_recoil%s_ext",var0[i]));
    h_recoil_overlay[i] = (TH1D*)fOverlay->Get(Form("h_recoil%s_total",var0[i]));
    h_recoil_dirt[i] = (TH1D*)fDirt->Get(Form("h_recoil%s_dirt_wgt",var0[i]));
  }
    
  for(int j=0; j< num_other_var; j++){
    h_other_covariance[j] = (TH2D*)file->Get(Form("h_2D_covariance%s",other_var[j]));
    h_other_bnb[j] = (TH1D*)fBNB->Get(Form("h%s_bnb",other_var[j]));
    h_other_ext[j] = (TH1D*)fEXT->Get(Form("h%s_ext",other_var[j]));
    h_other_overlay[j] = (TH1D*)fOverlay->Get(Form("h%s_total",other_var[j]));
    h_other_dirt[j] = (TH1D*)fDirt->Get(Form("h%s_dirt_wgt",other_var[j]));
  }
  
}//end of grab histograms

////////////////////////
// Calculate Prediction
/////////////////////////
TH1D* norm_and_shape::calculate_prediction(TH1D* h_overlay, TH1D* h_ext, TH1D* h_dirt){
  
  TH1D* h_prediction = (TH1D*)h_overlay->Clone();
  h_prediction->Add(h_ext);
  h_prediction->Add(h_dirt);

  h_overlay->Write();
  h_ext->Write();
  h_dirt->Write();
  
  return h_prediction;
  
} //end of calculate_prediction

//Creating the covariance matrix from the TH2D object                                                                                                                                                                              /
////////////////////////////////////////////////// 
TMatrixD norm_and_shape::convert_to_Tmatrix(TH2D* h_cov){

  TH2D* h_cov_clone   = (TH2D*)h_cov->Clone();
  int nbins_x = h_cov_clone->GetNbinsX();
  int nbins_y = h_cov_clone->GetNbinsY();
  TMatrixD covariance_matrix;

  covariance_matrix.Clear();
  covariance_matrix.ResizeTo(nbins_x,nbins_y);
  for(int i=0; i < nbins_x; i++){
    for(int j=0; j < nbins_y; j++){
      covariance_matrix[i][j] = h_cov_clone->GetBinContent(i+1,j+1);
    }
  }

  return covariance_matrix;
}

//Creating the covariance matrix from the TH2D object                                                                                                                                                                              /
////////////////////////////////////////////////// 
TVectorD norm_and_shape::convert_to_Tvector(TH1D* h_input){

  TH1D* h_clone   = (TH1D*)h_input->Clone();
  int nbins = h_clone->GetNbinsX();
  TVectorD output;

  output.Clear();
  output.ResizeTo(nbins);
  for(int i=0; i < nbins; i++){
      output[i] = h_clone->GetBinContent(i);
  }

  return output;
}

////////////////////////////////
//MatrixDecomp: Calculates Shape and Normalization Matrices
//////////////////////////////////////////////////////////

std::vector<TMatrixD> norm_and_shape::MatrixDecomp(int nbins, TVectorD matrix_pred,TMatrixD matrix_syst){

  // MiniBooNE note from Mike Schaevitz
  // https://microboone-docdb.fnal.gov/cgi-bin/sso/RetrieveFile?docid=5926&filename=tn253.pdf&version=1
	
  TMatrixD matrix_shape(nbins, nbins);
  TMatrixD matrix_mixed(nbins, nbins);
  TMatrixD matrix_norm(nbins, nbins);
  ///
  double N_T = 0;
  for (int idx = 0; idx < nbins; idx++) { N_T += matrix_pred(idx); }
  ///
  double M_kl = 0;
  for (int i = 0; i < nbins; i++) {
		
    for (int j = 0; j < nbins; j++) {
      
			
      M_kl += matrix_syst(i,j);
	
    }
  }
  ///
  for (int i = 0; i < nbins; i++) {
    for (int j = 0; j < nbins; j++) {	
  
      double N_i = matrix_pred(i);
      double N_j = matrix_pred(j);
      double M_ij = matrix_syst(i,j);	  
      double M_ik = 0; for(int k=0; k<nbins; k++) M_ik += matrix_syst(i,k);
      double M_kj = 0; for(int k=0; k<nbins; k++) M_kj += matrix_syst(k,j);
      matrix_shape(i,j) = M_ij - N_j*M_ik/N_T - N_i*M_kj/N_T + N_i*N_j*M_kl/N_T/N_T;
      matrix_mixed(i,j) = N_j*M_ik/N_T + N_i*M_kj/N_T - 2*N_i*N_j*M_kl/N_T/N_T;	
      matrix_norm(i,j) = N_i*N_j*M_kl/N_T/N_T;
    }
  }
  std::vector<TMatrixD> NormShapeVector = {matrix_norm+matrix_mixed,matrix_shape};
  return NormShapeVector;

} //end of MatrixDecomp



/////////////////////////////////////////
//plot_norm_and_shape():
//makes the actual total covariance matrix
//variable: variable we are plotting
//title: title for the plot
//hist_2D: vector of all the covariance matrices listed in file_names
//directory: where to save the image
/////////////////////////////////////////////
/*void norm_and_shape::plot_norm_and_shape(const char* variable, const char* title, std::vector<TH2D*> hist_2D){


  //Calculate total systematic uncertainty matrix
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size()-1; i++){
    hist->Add(hist_2D[i]);
  }

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  canv->SetRightMargin(0.15);
  //canv->SetLeftMargin(0.1);
  canv->SetBottomMargin(0.1);
  canv->SetTopMargin(0.1);
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size()-1; i++){
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
  if(std::strcmp(variable,"_muon_mom") ==  0 || std::strcmp(variable,"_leading_mom") ==  0 || std::strcmp(variable,"_recoil_mom") ==  0 || std::strcmp(variable,"_delta_PT") == 0){ ztitle = "/(GeV/c)";}
    else{ ztitle ="";}
  ztitle = "";
  //hist->GetZaxis()->SetTitle(Form("Covariance Element [10^{-76}%s]",ztitle));
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
  canv->Print(Form("../images/Total_Error/Covariance_Matrices/totaal_covariance_matrix_no_stat%s.png",variable));
  
  } //end of plot_total_covariance*/


