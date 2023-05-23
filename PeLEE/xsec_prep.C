#define xsec_prep_cxx
#include "xsec_prep.h"
#include "tools/helper_funcs.h"

void xsec_prep::main(){

  TFile *f_overlay=new TFile("root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root");//overlay histograms.
  TFile *f_dirt=new TFile("root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root");//dirt histograms 
  TFile *f_bnb=new TFile("root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root");//bnb histograms  
  TFile *f_ext=new TFile("root_files/pelee/Run_all/histograms_pelee_xsec_ext.root");//extbnb histograms 
  TFile *f_eff=new TFile("root_files/pelee/Run_all/histograms_mc_eff.root");//extbnb histograms  
  
  TFile *tfile = new TFile("root_files/pelee/Run_all/xsec_extraction.root","RECREATE"); //output root file
   
  //Grabbing Correct Date and Time and Directory
  ///////////////////////////////////////////////////
  const char* pathname = Form("images/pelee_xsec/Run_all/");
  string path(pathname);
  int dir_exists = dirExists(pathname);
  if(dir_exists == 0){
    mkdir(pathname,0777);
    std::cout<<"New Directory Succesfully Created"<<std::endl;
  } else if (dir_exists < 0){
    std::cout<<"An Error has Occured. Please Check the Input Path Name."<<std::endl;
  }else if(dir_exists > 0){
    std::cout<<"Directory Already Exists. Continuing with Analysis"<<std::endl;
    }

  //Stuff for Drawing
  //////////////////////
  TLatex* t = new TLatex(); //latex
  gStyle->SetPaintTextFormat("4.2f");gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  t->SetNDC();
  t->SetTextAlign(22);
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";
  
  //Calculate the Constants for the XSec
  ///////////////////////////////////////////
  Define_Constants();

  //Get the Histograms and Define the XSec ones to Save
  ///////////////////////////////////////////////////
  Grab_Histograms(f_bnb,f_ext,f_overlay,f_dirt, f_eff);


  //Determine Bin Size Test
  //Determine_Bin_Size(h_other_matrices[4]);

  
  //Plotting the Migration Matrices and checking normalization in each column (i.e. each truth bin)
  //////////////////////////////////////////////////////////////////////////////////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){

      //Drawing
      canv_particle_matrices[i][j] = new TCanvas(Form("canv_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),Form("canv_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),2000,1500);
      h_particle_matrices[i][j]->Draw("colz text");
      h_particle_matrices[i][j]->SetTitle(Form("%s: %s ; True %s; Reco. %s",titles_particles_matrices[i],titles_particle_matrices[j],titles_particle_matrices[j],titles_particle_matrices[j]));
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      h_particle_matrices[i][j]->Write();
      canv_particle_matrices[i][j]->Update();
      canv_particle_matrices[i][j]->Print(Form("%s%s%s_matrix.pdf",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      canv_particle_matrices[i][j]->Print(Form("%s%s%s_matrix.png",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));

      //checking normalization
      canv1_particle_matrices[i][j] = new TCanvas(Form("canv1_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),Form("canv1_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]),2000,1500);
      h_particle_matrices_normalized[i][j] = (TH2D*)h_particle_matrices[i][j]->Clone(Form("h_particle_matrices%s%s_normalized",particles_matrices[i],particles_matrices_plots[j]));

      int NBinsX =  h_particle_matrices_normalized[i][j]->GetXaxis()->GetNbins();
      int NBinsY = h_particle_matrices_normalized[i][j]->GetYaxis()->GetNbins();
      
      for(int bin_true = 0; bin_true < NBinsX; bin_true++){
	double NEventsInColumn = 0;

	for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
	  double bin_content = h_particle_matrices_normalized[i][j]->GetBinContent(bin_true+1,bin_reco+1);
	  NEventsInColumn += bin_content;
	}

	for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
	  if(NEventsInColumn > 0){
	    double FracErr =  h_particle_matrices_normalized[i][j]->GetBinError(bin_true+1,bin_reco+1) / h_particle_matrices_normalized[i][j]->GetBinContent(bin_true+1,bin_reco+1);
	    double CV = double(  h_particle_matrices_normalized[i][j]->GetBinContent(bin_true+1,bin_reco+1))/ double(NEventsInColumn);
	    h_particle_matrices_normalized[i][j]->SetBinContent(bin_true+1,bin_reco+1,CV);
	    
	    double error = CV * TMath::Sqrt( TMath::Power(FracErr,2.) + 1./double(NEventsInColumn) );
	    h_particle_matrices_normalized[i][j]->SetBinError(bin_true+1,bin_reco+1,error) ; 
	    
	  } else { 
	    h_particle_matrices_normalized[i][j]->SetBinContent(bin_true+1,bin_reco+1,0); 
	    h_particle_matrices_normalized[i][j]->SetBinError(bin_true+1,bin_reco+1,0); 
	  } //ends else  
	}
      }
      
      h_particle_matrices_normalized[i][j]->Draw("colz text");
      h_particle_matrices_normalized[i][j]->SetTitle(Form("%s: %s ; True %s; Reco. %s",titles_particles_matrices[i],titles_particle_matrices[j],titles_particle_matrices[j],titles_particle_matrices[j]));
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      h_particle_matrices_normalized[i][j]->Write();
      canv1_particle_matrices[i][j]->Update();
      canv1_particle_matrices[i][j]->Print(Form("%s%s%s_normalization.pdf",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      canv1_particle_matrices[i][j]->Print(Form("%s%s%s_normalization.png",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      
    }
  }

  for(int i=0; i < num_other_matrices; i++){

    //Drawing
    canv_other_matrices[i] = new TCanvas(Form("canv_other_matrices%s",other_matrices[i]),Form("canv_other_matrices%s",other_matrices[i]),2000,1500);
    h_other_matrices[i]->Draw("colz text"); 
    h_other_matrices[i]->SetTitle(Form("%s ; True %s; Reco. %s",titles_other_matrices[i],titles_other_matrices[i],titles_other_matrices[i]));
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
    canv_other_matrices[i]->Update();
    h_other_matrices[i]->Write();
    canv_other_matrices[i]->Print(Form("%s%s_matrix.pdf",path.c_str(),other_matrices[i]));
    canv_other_matrices[i]->Print(Form("%s%s_matrix.png",path.c_str(),other_matrices[i]));

    //checking normalization
    canv1_other_matrices[i] = new TCanvas(Form("canv1_other_matrices%s",other_matrices[i]),Form("canv1_particle_matrices%s",other_matrices[i]),2000,1500);
    h_other_matrices_normalized[i] = (TH2D*)h_other_matrices[i]->Clone(Form("h_other_matrices%s_normalized",other_matrices[i]));

    int NBinsX =  h_other_matrices_normalized[i]->GetXaxis()->GetNbins();
    int NBinsY = h_other_matrices_normalized[i]->GetYaxis()->GetNbins();
    
    for(int bin_true = 0; bin_true < NBinsX; bin_true++){
      double NEventsInColumn = 0;
      
      for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
	double bin_content = h_other_matrices_normalized[i]->GetBinContent(bin_true+1,bin_reco+1);
	NEventsInColumn += bin_content;
      }

      for(int bin_reco = 0; bin_reco < NBinsY; bin_reco++){
	if(NEventsInColumn > 0){
	  double FracErr =  h_other_matrices_normalized[i]->GetBinError(bin_true+1,bin_reco+1) / h_other_matrices_normalized[i]->GetBinContent(bin_true+1,bin_reco+1);
	  double CV = double(  h_other_matrices_normalized[i]->GetBinContent(bin_true+1,bin_reco+1))/ double(NEventsInColumn);
	  h_other_matrices_normalized[i]->SetBinContent(bin_true+1,bin_reco+1,CV);
	  
	  double error = CV * TMath::Sqrt( TMath::Power(FracErr,2.) + 1./double(NEventsInColumn) );
	  h_other_matrices_normalized[i]->SetBinError(bin_true+1,bin_reco+1,error) ; 
	  
	} else { 
	  h_other_matrices_normalized[i]->SetBinContent(bin_true+1,bin_reco+1,0); 
	  h_other_matrices_normalized[i]->SetBinError(bin_true+1,bin_reco+1,0); 
	} //ends else  
      } //end reco bin
    } //end true bin
      
    h_other_matrices_normalized[i]->Draw("colz text");
    h_other_matrices_normalized[i]->SetTitle(Form("%s ; True %s; Reco. %s",titles_other_matrices[i],titles_other_matrices[i],titles_other_matrices[i]));
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
    h_other_matrices_normalized[i]->Write();
    canv1_other_matrices[i]->Update();
    canv1_other_matrices[i]->Print(Form("%s%s_normalization.pdf",path.c_str(),other_matrices[i]));
    canv1_other_matrices[i]->Print(Form("%s%s_normalization.png",path.c_str(),other_matrices[i]));
    
  } 
  

  //Have to make the efficiency histograms and save them to the file for smearing matrix
  ////////////////////////////////////////////////////////////////////////////////////

  //Momentum theta and phi
  for(int i =0; i < num_particles_eff; i++){
    for (int j=0; j <num_particles_eff_plots; j++){
      canv_particle_eff[i][j] = new TCanvas(Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),2000,1500);
      h_particle_num[i][j]->Divide(h_particle_num[i][j],h_particle_denom[i][j],1.0,1.0, "B");
      h_particle_num[i][j]->Draw("1e1p");
      h_particle_num[i][j]->SetTitle(Form(" ; %s ; Efficiency",particles_eff_var_titles[j]));
      h_particle_num[i][j]->SetLineColor(kViolet);
      h_particle_num[i][j]->SetMaximum(0.4);
      h_particle_num[i][j]->SetMinimum(0);
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s of %s}",particles_eff_var_titles[j],particles_eff_titles[i]));
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
      //a[i] = new TLine(xlim_eff[i],0,xlim_eff[i],1);
      //a[i]->Draw("same");
      //a[i]->SetLineColor(kBlack);
      //a[i]->SetLineWidth(4);
      //a1[i] = new TLine(xlim_high_eff[i],0,xlim_high_eff[i],1);
      //a1[i]->Draw("same");
      //a1[i]->SetLineColor(kBlack);
      //a1[i]->SetLineWidth(4);
      
      h_particle_num[i][j]->Write();
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff_4_smearing.png",path.c_str(),particles_eff[i],particles_eff_var[j]));
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff_4_smearing.pdf",path.c_str(),particles_eff[i],particles_eff_var[j])); 
    }
  }

  //all other variables
  for(int i =0; i < num_other_eff; i++){
    canv_other_eff[i] = new TCanvas(Form("C_other_eff%s",other_eff[i]),Form("C_other_eff%s",other_eff[i]),2000,1500);
    h_other_eff_num[i]->Divide(h_other_eff_num[i],h_other_eff_denom[i],1.0,1.0, "B");
    h_other_eff_num[i]->Draw("1e1p");
    h_other_eff_num[i]->SetTitle(Form(" ; %s ; Efficiency",other_eff_titles[i]));
    h_other_eff_num[i]->SetLineColor(kViolet);
    h_other_eff_num[i]->SetMaximum(0.4);
    h_other_eff_num[i]->SetMinimum(0);
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",other_eff_titles[i]));
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
    h_other_eff_num[i]->Write();
    canv_other_eff[i]->Print(Form("%s%s_eff_4_smearing.png",path.c_str(),other_eff[i]));
    canv_other_eff[i]->Print(Form("%s%s_eff_4_smearing.pdf",path.c_str(),other_eff[i]));

  }


  //Creating the Smearing Matrices
  //////////////////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){

      canv_particle_smearing[i][j] =  new TCanvas(Form("canv_particle_smearing%s%s",particles_matrices[i],particles_matrices_plots[j]),Form("canv_particle_smearing%s%s",particles_matrices[i],particles_matrices_plots[j]),2000,1500);
      h_particle_smearing[i][j] = (TH2D*)h_particle_matrices[i][j]->Clone(Form("h_particle_matrices%s%s_smearing",particles_matrices[i],particles_matrices_plots[j]));
      double num_bins_true = double(h_particle_smearing[i][j]->GetXaxis()->GetNbins());
      double num_bins_reco = double(h_particle_smearing[i][j]->GetYaxis()->GetNbins()); 
      
      for(int jj=1; jj < num_bins_reco+1; jj++){
	double sum = 0;
	for(int k=1; k < num_bins_true+1; k++){
	  double bin_content = h_particle_smearing[i][j]->GetBinContent(k,jj);
	  sum += 0 + bin_content;
	}
	for(int k=1; k < num_bins_true+1; k++){
	  double bin_content = h_particle_smearing[i][j]->GetBinContent(k,jj);
	  h_particle_smearing[i][j]->SetBinContent(k,jj,(bin_content/sum));
	}
      } //end of loop over true bins

      h_particle_smearing[i][j]->Draw("colz text"); //will add text later
      h_particle_smearing[i][j]->SetTitle(Form("%s: %s ; True %s; Reco. %s",titles_particles_matrices[i],titles_particle_matrices[j],titles_particle_matrices[j],titles_particle_matrices[j]));
      h_particle_smearing[i][j]->Write();
      t->DrawLatex(0.3,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      canv_particle_smearing[i][j]->Update();
      canv_particle_smearing[i][j]->Print(Form("%s%s%s_smearing.pdf",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));
      canv_particle_smearing[i][j]->Print(Form("%s%s%s_smearing.png",path.c_str(),particles_matrices[i],particles_matrices_plots[j]));        
    }
  }

  for(int i=0; i <num_other_matrices; i++){

    canv_other_smearing[i] = new TCanvas(Form("canv_other_smearing%s",other_matrices[i]),Form("canv_other_smearing%s",other_matrices[i]),2000,1500);
    h_other_smearing[i] = (TH2D*)h_other_matrices[i]->Clone(Form("h_other_matrices%s_smearing",other_matrices[i]));
    double num_bins_true = double(h_other_smearing[i]->GetXaxis()->GetNbins());
    double num_bins_reco = double(h_other_smearing[i]->GetYaxis()->GetNbins());
  
    for(int jj=1; jj < num_bins_reco+1; jj++){
      double sum = 0;
      for(int k=1; k < num_bins_true+1; k++){
	double bin_content = h_other_smearing[i]->GetBinContent(k,jj);
	sum += 0 + bin_content;
      }
      for(int k=1; k < num_bins_true+1; k++){
	double bin_content = h_other_smearing[i]->GetBinContent(k,jj);
	double eff = h_other_eff_num[i]->GetBinContent(k);
	double eff_correction;
	if(eff == 0.0 | eff < 0.025){
	  eff_correction = 0.0;
	}else{
	  eff_correction = 1.0/eff;
	}
	h_other_smearing[i]->SetBinContent(k,jj,(bin_content/sum)); //remember to apply eff_correction
      }
    } //end of loop over reco bins

    h_other_smearing[i]->Draw("colz text"); //add text later
    h_other_smearing[i]->SetTitle(Form("%s ; True %s; Reco. %s",titles_other_matrices[i],titles_other_matrices[i],titles_other_matrices[i]));
    h_other_smearing[i]->Write();
    t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.82,0.92,Form("%s",sample_name));
    canv_other_smearing[i]->Update();
    canv_other_smearing[i]->Print(Form("%s%s_smearing.pdf",path.c_str(),other_matrices[i]));
    canv_other_smearing[i]->Print(Form("%s%s_smearing.png",path.c_str(),other_matrices[i]));
  }
  

  //Close File Before Finishing:
  /////////////////////////////////////
  tfile->Close();
  
} //end of program
