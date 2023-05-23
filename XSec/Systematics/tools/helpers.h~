#include "Covariance_Matrix.h"
#include "constants.h"
using namespace Constants;

class helper{

 public:
  
  void Grab_Histograms(const char* input_file, int universes, bool genie_unisim, const char* file_type = "");


}; //end of class helper

////////////////////////////////////////////////////
//Function has 2 parts
//1) grabs the centerval values. Same for all the files
//2) grabs the specific histograms from a specific file
///////////////////////////////////////////////////
void helper::Grab_Histograms(const char* input_file, int num_universes, bool genie_unisim, const char* file_type = ""){

  //CV values
  ///////////
  TFile* file_cv = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));
  
  for(int j=0; j < num_var; j++){
    
    h_muon_CV[j] = (TH1D*)file_cv->Get(Form("h_muon%s_total",var0[j]));
    h_muon_EXT[j] = (TH1D*)file_ext->Get(Form("h_muon%s_ext",var0[j]));
    h_muon_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_muon%s_dirt_wgt",var0[j]));
    
    h_leading_CV[j] = (TH1D*)file_cv->Get(Form("h_leading%s_total",var0[j]));
    h_leading_EXT[j] = (TH1D*)file_ext->Get(Form("h_leading%s_ext",var0[j]));
    h_leading_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_leading%s_dirt_wgt",var0[j]));

    h_recoil_CV[j] = (TH1D*)file_cv->Get(Form("h_recoil%s_total",var0[j]));
    h_recoil_EXT[j] = (TH1D*)file_ext->Get(Form("h_recoil%s_ext",var0[j]));
    h_recoil_Dirt[j] = (TH1D*)file_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[j]));

  } 

  for(int j=0; j < num_other_var; j++){
    h_other_CV[j] = (TH1D*)file_cv->Get(Form("h%s_total",other_var[j]));
    h_other_EXT[j] = (TH1D*)file_ext->Get(Form("h%s_ext",other_var[j]));
    h_other_Dirt[j] = (TH1D*)file_dirt->Get(Form("h%s_dirt_wgt",other_var[j]));
  }

  //Universe values
  ///////////////
  TFile* file = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/Systematics/root_files/%s/histograms_pelee_xsec.root",input_file));

  //Unisims
  ///////////
  if(genie_unisim == true){

    if(strncmp(file_type,"RPA_CCQE_UBGenie",16) == 0){

      //values stored in _0 and _1
      for(int i = 0; i < num_universes; i++){
	for(int j=0; j < num_var; j++){
	  h_muon[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_%d",var0[j],i));
	  h_leading[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_%d",var0[j],i));
	  h_recoil[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_%d",var0[j],i));
	}
	for(int j=0; j < num_other_var; j++){
	  h_other[i][j] = (TH1D*)file->Get(Form("h%s_total_%d",other_var[j],i));
	}
      }

     
    } else if (strncmp(file_type,"XSecShape_CCMEC_UBGenie",23) == 0){

      //Values are stored in _1
      for(int j=0; j < num_var; j++){
	h_muon[0][j] = (TH1D*)file->Get(Form("h_muon%s_total_1",var0[j]));
	h_leading[0][j] = (TH1D*)file->Get(Form("h_leading%s_total_1",var0[j]));
	h_recoil[0][j] = (TH1D*)file->Get(Form("h_recoil%s_total_1",var0[j]));
      }
      for(int j=0; j < num_other_var; j++){
	h_other[0][j] = (TH1D*)file->Get(Form("h%s_total_1",other_var[j]));
      }
      
    } else {
    
      //Values are stored in _0      
      for(int j=0; j < num_var; j++){
	h_muon[0][j] = (TH1D*)file->Get(Form("h_muon%s_total_0",var0[j]));
	h_leading[0][j] = (TH1D*)file->Get(Form("h_leading%s_total_0",var0[j]));
	h_recoil[0][j] = (TH1D*)file->Get(Form("h_recoil%s_total_0",var0[j]));
      }
      for(int j=0; j < num_other_var; j++){
	h_other[0][j] = (TH1D*)file->Get(Form("h%s_total_0",other_var[j]));
      }
      
      
    }
    
    //multisims
    /////////////
  } else {

    for(int i=0; i < num_universes; i++){
      for(int j=0; j < num_var; j++){
	h_muon[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_%d",var0[j],i));
	h_leading[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_%d",var0[j],i));
	h_recoil[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_%d",var0[j],i));
      }

      for(int j=0; j < num_other_var; j++){
	h_other[i][j] = (TH1D*)file->Get(Form("h%s_total_%d",other_var[j],i));
      }
    }
    
 } //end of else
  
}//end of grab_files
