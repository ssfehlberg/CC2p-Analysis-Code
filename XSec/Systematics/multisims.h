#include "tools/xsec.h"
#include "tools/plotting.h"
#include "tools/Covariance_Matrix.h"

class Process_Multisims{

 public:
  void Grab_Histograms(const char* input_file, int num_universes, bool genie_unisim, const char* file_type = "");
  void Grab_Iterations();
  void main();

  //XSEc class
  //////////////
  xsec Xsec;
  
  //Stuff for individual files
  /////////////////////////////
  /*static const int number_of_files = 15;
  const char* directory_name;
  const char* directory_name_list[number_of_files] = {"flux_all","reint_all","All_UBGenie",
						"AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
						"NormNCCOH_UBGenie","ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
						"VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","xsr_scc_Fa3_SCC",
						"xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie","RootinoFix_UBGenie"};
  int universes; //1000,1000,500
  std::vector<int> uni = {1000,1000,500,
			  1,1,1,
			  1,1,1,
			  1,1,1,
			  1,2,1};
  bool unisim[number_of_files] = {false,false,false,
			    true,true,true,
			    true,true,true,
			    true,true,true,
			    true,true,true};
  */

  //Stuff for individual files                                                                                                                                                                                               /////////////////////////////
  static const int number_of_files = 12;
  const char* directory_name;
  const char* directory_name_list[number_of_files] = {"AxFFCCQEshape_UBGenie","DecayAngMEC_UBGenie","NormCCCOH_UBGenie",
						      "NormNCCOH_UBGenie","ThetaDelta2NRad_UBGenie","Theta_Delta2Npi_UBGenie",
						      "VecFFCCQEshape_UBGenie","XSecShape_CCMEC_UBGenie","xsr_scc_Fa3_SCC",
						      "xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie","RootinoFix_UBGenie"};
  int universes; //1000,1000,500                                                                                                                                                                                            
  std::vector<int> uni = {1,1,1,
                          1,1,1,
                          1,1,1,
                          1,2,1};
  bool unisim[number_of_files] = {true,true,true,
				  true,true,true,
				  true,true,true,
				  true,true,true};
  /*
  
  //THIS IS FOR TESTING ONLY
  static const int number_of_files = 3;
  const char* directory_name;
  const char* directory_name_list[number_of_files] = {"xsr_scc_Fa3_SCC","xsr_scc_Fv3_SCC","RPA_CCQE_UBGenie"};
  int universes;
  std::vector<int> uni = {1,1,2};
  bool unisim[number_of_files] = {true, true,true};
  */
  

 //CV Files we will grab CV plots from
  /////////////////////////////////
  TFile* file_overlay = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_overlay_wgt.root "));
  TFile* file_eff = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_mc_eff.root"));
  TFile* file_bnb = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_bnb.root "));
  TFile* file_ext = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_ext.root "));
  TFile* file_dirt = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/root_files/pelee/Run_all/histograms_pelee_xsec_dirt_wgt.root "));

  /////////////////////////////////
  //Stuff Needed for the Overlay Plots
  ////////////////////////////////
  static const int num_channels = 2;
  const char* channel[num_channels] = {"_total", "_cc2p0pi"};

  ///////////////////////////
  //Stuff for the particle plots
  //////////////////////////////
  static const int num_var = 3;
  const char* var0[num_var] = {"_mom","_theta","_phi"};
  const char* var[num_var] = {"_mom","_costheta","_phi"};
  const char* var_titles[num_var] = {"Momentum (GeV/c)","cos(#theta) (Rad.)","#phi (Rad.)"};

  static const int num_particles = 3;
  const char* particles[num_particles] = {"_muon","_leading","_recoil"};
  const char* particles_titles[num_particles] = {"Muon","Leading Proton","Recoiling Proton"};

  //Muon
  ///////////////
  TH1D* h_bnb_muon[num_var]; //CV BNB
  TH1D* h_ext_muon[num_var]; //CV EXT
  TH1D* h_dirt_muon[num_var]; //CV Dirt
  TH1D* h_overlay_muon[num_var][num_channels]; //CV Overlay for total and cc2p    
  TH1D* h_muon_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_muon_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_muon_num_CV[num_var]; //CV muon numberator of efficiency
  TH1D* h_muon_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_muon_denom_CV[num_var]; //CV muon denominator of the efficiency
  TH1D* h_muon_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_muon_matrices_CV[num_var]; //the CV muon matrices
  TH2D* h_muon_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_muon_2d[num_var];
  
  TH1D* h_muon_xsec[1000][num_var]; //cross-section in each universe
  double muon_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_muon[num_var] = {false,true,false};
  int muon_iter[num_var] = {2,2,1}; //number of iterations used in the unfolding. Really dumbly defined

  //Multisim drawing stuff
  TCanvas* canv_muon[num_var]; 
  TLegend* legend_muon[num_var];
  double ymax_muon[num_var] = {2000,2000,2000};
  double ymin_muon[num_var] = {0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_muon[num_var];
  TCanvas* canv_2D_muon[num_var];

  //Error drawing stuff
  TH1D* h_muon_error[num_var];
  TCanvas* canv_muon_error[num_var];

  //fractional error stuff
  TH1D* h_muon_fractional_error[num_var];
  TCanvas* canv_muon_fractional_error[num_var];


  //Leading
  ///////////////
  TH1D* h_bnb_leading[num_var]; //CV BNB
  TH1D* h_ext_leading[num_var]; //CV EXT
  TH1D* h_dirt_leading[num_var]; //CV Dirt
  TH1D* h_overlay_leading[num_var][num_channels]; //CV Overlay for total and cc2p    
  TH1D* h_leading_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_leading_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_leading_num_CV[num_var]; //CV leading numberator of efficiency
  TH1D* h_leading_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_leading_denom_CV[num_var]; //CV leading denominator of the efficiency
  TH1D* h_leading_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_leading_matrices_CV[num_var]; //the CV leading matrices
  TH2D* h_leading_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_leading_2d[num_var];
  
  TH1D* h_leading_xsec[1000][num_var]; //cross-section in each universe                                                                                                                                                              
  double leading_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_leading[num_var] = {false,true,false};
  int leading_iter[num_var] = {1,2,1}; //number of iterations used in the unfolding. Really dumbly defined                                                                                                                           
  
  //Multisim drawing stuff                                                                                                                                                                                                           
  TCanvas* canv_leading[num_var];
  TLegend* legend_leading[num_var];
  double ymax_leading[num_var] = {2000,2000,2000};
  double ymin_leading[num_var] = {0,0,0};
  
  //matrix drawing stuff
  TH2D* h_2D_leading[num_var];
  TCanvas* canv_2D_leading[num_var];

  //Error drawing stuff
  TH1D* h_leading_error[num_var];
  TCanvas* canv_leading_error[num_var];

  //fractional error stuff
  TH1D* h_leading_fractional_error[num_var];
  TCanvas* canv_leading_fractional_error[num_var];


  //Recoil
  ///////////////
  TH1D* h_bnb_recoil[num_var]; //CV BNB
  TH1D* h_ext_recoil[num_var]; //CV EXT
  TH1D* h_dirt_recoil[num_var]; //CV Dirt
  TH1D* h_overlay_recoil[num_var][num_channels]; //CV Overlay for total and cc2p    
  TH1D* h_recoil_total[1000][num_var]; //total MC contribution in a universe
  TH1D* h_recoil_cc2p[1000][num_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_recoil_num_CV[num_var]; //CV recoil numberator of efficiency
  TH1D* h_recoil_num[1000][num_var]; //variations of numberator of efficiency
  TH1D* h_recoil_denom_CV[num_var]; //CV recoil denominator of the efficiency
  TH1D* h_recoil_denom[1000][num_var]; //variations of denominator of efficiency
  
  TH2D* h_recoil_matrices_CV[num_var]; //the CV recoil matrices
  TH2D* h_recoil_matrices[1000][num_var]; //the variation matrices
  TCanvas* canv_recoil_2d[num_var];
  
  TH1D* h_recoil_xsec[1000][num_var]; //cross-section in each universe                                                                                                                                                               
  double recoil_xsec_max[num_var] = {17E-38,12E-38,2E-38};
  bool flip_legend_recoil[num_var] = {false,true,false};
  int recoil_iter[num_var] = {2,2,2}; //number of iterations used in the unfolding. Really dumbly defined                                                                                                                            

  //Multisim drawing stuff                                                                                                                                                                                                           
  TCanvas* canv_recoil[num_var];
  TLegend* legend_recoil[num_var];
  double ymax_recoil[num_var] = {2000,2000,2000};
  double ymin_recoil[num_var] = {0,0,0};

  //matrix drawing stuff
  TH2D* h_2D_recoil[num_var];
  TCanvas* canv_2D_recoil[num_var];

  //Error drawing stuff
  TH1D* h_recoil_error[num_var];
  TCanvas* canv_recoil_error[num_var];

  //fractional error stuff
  TH1D* h_recoil_fractional_error[num_var];
  TCanvas* canv_recoil_fractional_error[num_var];

  TH1D* h_recoil[1000][num_var];

  //Other Variables
  /////////////////
  static const int num_other_var = 5;
  const char* other_var[num_other_var] = {"_opening_angle_protons_lab","_opening_angle_mu_both",
                                          "_delta_PT","_delta_alphaT","_delta_phiT"};
  const char* other_var_titles[num_other_var] = {"cos(#gamma_{Lab})","cos(#gamma_{#mu,p_{L} + p_{R}})",
                                                 "#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)"};

  TH1D* h_bnb_other[num_other_var]; //CV BNB
  TH1D* h_ext_other[num_other_var]; //CV EXT
  TH1D* h_dirt_other[num_other_var]; //CV Dirt
  TH1D* h_overlay_other[num_other_var][num_channels]; //CV Overlay for total and cc2p    
  TH1D* h_other_total[1000][num_other_var]; //total MC contribution in a universe
  TH1D* h_other_cc2p[1000][num_other_var]; //CC2p MC contribution in a universe
  
  //cross-section drawing stuff
  TH1D* h_other_num_CV[num_other_var]; //CV other numberator of efficiency
  TH1D* h_other_num[1000][num_other_var]; //variations of numberator of efficiency
  TH1D* h_other_denom_CV[num_other_var]; //CV other denominator of the efficiency
  TH1D* h_other_denom[1000][num_other_var]; //variations of denominator of efficiency
  
  TH2D* h_other_matrices_CV[num_other_var]; //the CV other matrices
  TH2D* h_other_matrices[1000][num_other_var]; //the variation matrices
  TCanvas* canv_other_2d[num_other_var];
  
  TH1D* h_other_xsec[1000][num_other_var]; //cross-section in each universe                                                                                                                                                          
  double other_xsec_max[num_other_var] = {8E-38,8E-38,25E-38,0.18E-38,0.1E-38};//,2E-38,2E-38};                                                                                                                                      
  bool flip_legend_other[num_other_var] = {false,true,false,false,false};//,false,false};                                                                                                                                            
  int other_iter[num_other_var] = {2,2,1,1,2}; //number of iterations used in the unfolding. Really dumbly defined  

  //Multisim drawing stuff
  TCanvas* canv_other[num_other_var];
  TLegend* legend_other[num_other_var];
  double ymax_other[num_other_var] = {2000,2000,2000,2000,2000};
  double ymin_other[num_other_var] = {0,0,0,0,0};
 
  //matrix drawing stuff
  TH2D* h_2D_other[num_other_var];
  TCanvas* canv_2D_other[num_other_var];

  //Error drawing stuff
  TH1D* h_other_error[num_other_var];
  TCanvas* canv_other_error[num_other_var];

  //fractional error stuff
  TH1D* h_other_fractional_error[num_other_var];
  TCanvas* canv_other_fractional_error[num_other_var];
  
};

void Process_Multisims::Grab_Histograms(const char* input_file, int num_universes, bool genie_unisim, const char* file_type = ""){

  std::cout<<"[GRAABBING ALL THE FILES]"<<std::endl;

  ////////////////////////////////////////////////////////////////////////
  //Grab the overlay smearing matrix and efficiency for all the variables
  ///////////////////////////////////////////////////////////////////////
  
  for(int i=0; i < num_var; i++){

    h_bnb_muon[i] = (TH1D*)file_bnb->Get(Form("h_muon%s_bnb",var0[i])); //bnb
    h_ext_muon[i] = (TH1D*)file_ext->Get(Form("h_muon%s_ext",var0[i])); //ext
    h_dirt_muon[i] = (TH1D*)file_dirt->Get(Form("h_muon%s_dirt_wgt",var0[i])); //dirt
    h_overlay_muon[i][0] = (TH1D*)file_overlay->Get(Form("h_muon%s_total",var0[i])); //total CV overlay
    h_overlay_muon[i][1] = (TH1D*)file_overlay->Get(Form("h_muon%s_cc2p0pi",var0[i])); //cc2p0pi CV overlay
    h_muon_matrices_CV[i] =  (TH2D*)file_eff->Get(Form("h_particle_matrices_muon_all%s",var[i])); //CV matrix
    h_muon_num_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_num_muon_all%s",var[i])); //CV numerator 
    h_muon_denom_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_denom_muon_all%s",var[i])); //CV denominator

    h_bnb_leading[i] = (TH1D*)file_bnb->Get(Form("h_leading%s_bnb",var0[i])); //bnb
    h_ext_leading[i] = (TH1D*)file_ext->Get(Form("h_leading%s_ext",var0[i])); //ext
    h_dirt_leading[i] = (TH1D*)file_dirt->Get(Form("h_leading%s_dirt_wgt",var0[i])); //dirt
    h_overlay_leading[i][0] = (TH1D*)file_overlay->Get(Form("h_leading%s_total",var0[i])); //total CV overlay
    h_overlay_leading[i][1] = (TH1D*)file_overlay->Get(Form("h_leading%s_cc2p0pi",var0[i])); //cc2p0pi CV overlay
    h_leading_matrices_CV[i] =  (TH2D*)file_eff->Get(Form("h_particle_matrices_lead_proton%s",var[i])); //CV matrix
    h_leading_num_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_num_lead_proton%s",var[i])); //CV numerator 
    h_leading_denom_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_denom_lead_proton%s",var[i])); //CV denominator

    h_bnb_recoil[i] = (TH1D*)file_bnb->Get(Form("h_recoil%s_bnb",var0[i])); //bnb
    h_ext_recoil[i] = (TH1D*)file_ext->Get(Form("h_recoil%s_ext",var0[i])); //ext
    h_dirt_recoil[i] = (TH1D*)file_dirt->Get(Form("h_recoil%s_dirt_wgt",var0[i])); //dirt
    h_overlay_recoil[i][0] = (TH1D*)file_overlay->Get(Form("h_recoil%s_total",var0[i])); //total CV overlay
    h_overlay_recoil[i][1] = (TH1D*)file_overlay->Get(Form("h_recoil%s_cc2p0pi",var0[i])); //cc2p0pi CV overlay
    h_recoil_matrices_CV[i] =  (TH2D*)file_eff->Get(Form("h_particle_matrices_recoil_proton%s",var[i])); //CV matrix
    h_recoil_num_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_num_recoil_proton%s",var[i])); //CV numerator 
    h_recoil_denom_CV[i] = (TH1D*)file_eff->Get(Form("h_particle_denom_recoil_proton%s",var[i])); //CV denominator

  }

  for(int i=0; i < num_other_var; i++){
    h_bnb_other[i] = (TH1D*)file_bnb->Get(Form("h%s_bnb",other_var[i])); //bnb
    h_ext_other[i] = (TH1D*)file_ext->Get(Form("h%s_ext",other_var[i])); //ext
    h_dirt_other[i] = (TH1D*)file_dirt->Get(Form("h%s_dirt_wgt",other_var[i])); //dirt
    h_overlay_other[i][0] = (TH1D*)file_overlay->Get(Form("h%s_total",other_var[i])); //total CV overlay
    h_overlay_other[i][1] = (TH1D*)file_overlay->Get(Form("h%s_cc2p0pi",other_var[i])); //cc2p0pi CV overlay
    h_other_matrices_CV[i] =  (TH2D*)file_eff->Get(Form("h_other_matrices%s",other_var[i])); //CV matrix
    h_other_num_CV[i] = (TH1D*)file_eff->Get(Form("h_other_eff_num%s",other_var[i])); //CV numerator 
    h_other_denom_CV[i] = (TH1D*)file_eff->Get(Form("h_other_eff_denom%s",other_var[i])); //CV denominator
  }
  

  /////////////////////////////////////////
  //Now we have to grab the universe values
  /////////////////////////////////////////
  TFile* file = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/%s/histograms_pelee_xsec.root",input_file));
  TFile* file_efficiency = new TFile(Form("/Users/ssfehlberg/Research/Thesis/2Proton_Pandora/XSec/Systematics/root_files/%s/histograms_mc_eff.root",input_file));
  
  //Unisims
  ///////////
  if(genie_unisim == true){

    if(strncmp(file_type,"RPA_CCQE_UBGenie",16) == 0){

      //values stored in _0 and _1
      for(int i = 0; i < num_universes; i++){
	for(int j=0; j < num_var; j++){
	  h_muon_total[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_%d",var0[j],i));
	  h_muon_cc2p[i][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_%d",var0[j],i));
	  h_muon_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_%d",var[j],i)); //CV matrix
	  h_muon_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_%d",var[j],i)); //CV numerator 
	  h_muon_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_%d",var[j],i)); //CV denominator

	  h_leading_total[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_%d",var0[j],i));
	  h_leading_cc2p[i][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_%d",var0[j],i));
	  h_leading_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_%d",var[j],i)); //CV matrix
	  h_leading_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_%d",var[j],i)); //CV numerator 
	  h_leading_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_%d",var[j],i)); //CV denominator

	  h_recoil_total[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_%d",var0[j],i));
	  h_recoil_cc2p[i][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_%d",var0[j],i));
	  h_recoil_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_%d",var[j],i)); //CV matrix
	  h_recoil_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_%d",var[j],i)); //CV numerator 
	  h_recoil_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_%d",var[j],i)); //CV denominator
	  
	}

	for(int j=0; j < num_other_var; j++){
	  h_other_total[i][j] = (TH1D*)file->Get(Form("h%s_total_%d",other_var[j],i));
	  h_other_cc2p[i][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_%d",other_var[j],i));
	  h_other_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_%d",other_var[j],i)); //CV matrix
	  h_other_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_%d",other_var[j],i)); //CV numerator 
	  h_other_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_%d",other_var[j],i)); //CV denominator
	}
      }

     
    } else if (strncmp(file_type,"XSecShape_CCMEC_UBGenie",23) == 0){

      //Values are stored in _1
      for(int j=0; j < num_var; j++){
	h_muon_total[0][j] = (TH1D*)file->Get(Form("h_muon%s_total_1",var0[j]));
	h_muon_cc2p[0][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_1",var0[j]));
	h_muon_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_1",var[j])); //CV matrix
	h_muon_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_1",var[j])); //CV numerator 
	h_muon_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_1",var[j])); //CV denominator

	h_leading_total[0][j] = (TH1D*)file->Get(Form("h_leading%s_total_1",var0[j]));
	h_leading_cc2p[0][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_1",var0[j]));
	h_leading_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_1",var[j])); //CV matrix
	h_leading_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_1",var[j])); //CV numerator 
	h_leading_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_1",var[j])); //CV denominator
	
	h_recoil_total[0][j] = (TH1D*)file->Get(Form("h_recoil%s_total_1",var0[j]));
	h_recoil_cc2p[0][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_1",var0[j]));
	h_recoil_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_1",var[j])); //CV matrix
	h_recoil_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_1",var[j])); //CV numerator 
	h_recoil_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_1",var[j])); //CV denominator
	
      }
      
      for(int j=0; j < num_other_var; j++){
	h_other_total[0][j] = (TH1D*)file->Get(Form("h%s_total_1",other_var[j]));
	h_other_cc2p[0][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_1",other_var[j]));
	h_other_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_1",other_var[j])); //CV matrix
	h_other_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_1",other_var[j])); //CV numerator 
	h_other_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_1",other_var[j])); //CV denominator
      }
      
    } else {
    
      //Values are stored in _0      
      for(int j=0; j < num_var; j++){
	h_muon_total[0][j] = (TH1D*)file->Get(Form("h_muon%s_total_0",var0[j]));
	h_muon_cc2p[0][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_0",var0[j]));
	h_muon_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_0",var[j])); //CV matrix
	h_muon_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_0",var[j])); //CV numerator 
	h_muon_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_0",var[j])); //CV denominator

	h_leading_total[0][j] = (TH1D*)file->Get(Form("h_leading%s_total_0",var0[j]));
	h_leading_cc2p[0][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_0",var0[j]));
	h_leading_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_0",var[j])); //CV matrix
	h_leading_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_0",var[j])); //CV numerator 
	h_leading_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_0",var[j])); //CV denominator
	
	h_recoil_total[0][j] = (TH1D*)file->Get(Form("h_recoil%s_total_0",var0[j]));
	h_recoil_cc2p[0][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_0",var0[j]));
	h_recoil_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_0",var[j])); //CV matrix
	h_recoil_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_0",var[j])); //CV numerator 
	h_recoil_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_0",var[j])); //CV denominator
	
      }
      
      for(int j=0; j < num_other_var; j++){
	h_other_total[0][j] = (TH1D*)file->Get(Form("h%s_total_0",other_var[j]));
	h_other_cc2p[0][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_0",other_var[j]));
	h_other_matrices[0][j] =  (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_0",other_var[j])); //CV matrix
	h_other_num[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_0",other_var[j])); //CV numerator 
	h_other_denom[0][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_0",other_var[j])); //CV denominator
      }
      
    }
    
    //multisims
    /////////////
  } else {

    for(int i=0; i < num_universes; i++){
      for(int j=0; j < num_var; j++){
	h_muon_total[i][j] = (TH1D*)file->Get(Form("h_muon%s_total_%d",var0[j],i));
	h_muon_cc2p[i][j] = (TH1D*)file->Get(Form("h_muon%s_cc2p0pi_%d",var0[j],i));
	h_muon_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_muon_all%s_%d",var[j],i)); //CV matrix
	h_muon_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_muon_all%s_%d",var[j],i)); //CV numerator 
	h_muon_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_muon_all%s_%d",var[j],i)); //CV denominator


	h_leading_total[i][j] = (TH1D*)file->Get(Form("h_leading%s_total_%d",var0[j],i));
	h_leading_cc2p[i][j] = (TH1D*)file->Get(Form("h_leading%s_cc2p0pi_%d",var0[j],i));
	h_leading_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_lead_proton%s_%d",var[j],i)); //CV matrix
	h_leading_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_lead_proton%s_%d",var[j],i)); //CV numerator 
	h_leading_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_lead_proton%s_%d",var[j],i)); //CV denominator
	
	h_recoil_total[i][j] = (TH1D*)file->Get(Form("h_recoil%s_total_%d",var0[j],i));
	h_recoil_cc2p[i][j] = (TH1D*)file->Get(Form("h_recoil%s_cc2p0pi_%d",var0[j],i));
	h_recoil_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_particle_matrices_recoil_proton%s_%d",var[j],i)); //CV matrix
	h_recoil_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_num_recoil_proton%s_%d",var[j],i)); //CV numerator 
	h_recoil_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_particle_denom_recoil_proton%s_%d",var[j],i)); //CV denominator
	
      }
      
      for(int j=0; j < num_other_var; j++){
	h_other_total[i][j] = (TH1D*)file->Get(Form("h%s_total_%d",other_var[j],i));
	h_other_cc2p[i][j] = (TH1D*)file->Get(Form("h%s_cc2p0pi_%d",other_var[j],i));
	h_other_matrices[i][j] =  (TH2D*)file_efficiency->Get(Form("h_other_matrices%s_%d",other_var[j],i)); //CV matrix
	h_other_num[i][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_num%s_%d",other_var[j],i)); //CV numerator 
	h_other_denom[i][j] = (TH1D*)file_efficiency->Get(Form("h_other_eff_denom%s_%d",other_var[j],i)); //CV denominator
      }	
    }//end of for loop over universes  
  } //end of else

 std::cout<<"[FINISHED GRAABBING ALL THE FILES]"<<std::endl;

} //end of grab histograms

void Process_Multisims::Grab_Iterations(){

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

void Process_Multisims::main(){

  for(int s = 0; s < number_of_files; s++){

    directory_name = directory_name_list[s];
    universes = uni[s];
    Grab_Histograms(directory_name,universes,unisim[s],directory_name);
  
    ///////////////////////////
    //Now to fucking plot shit
    //////////////////////////

    std::cout<<Form("[Starting to Process the %s Multisims]",directory_name)<<std::endl;
    
    //File to save shit to
    TFile *tfile_mine = new TFile(Form("root_files/%s/systematics.root",directory_name),"RECREATE"); //output root file   
  
    for(int j=0; j < num_var; j++){


      //Creating the cross-section
      //h_muon_xsec[0][j] = CV cross-section using h_muon[0][j] (i.e. dirt = 100%)
      //h_muon_xsec[1][j] = Modified cross-section using h_muon[1][j] (i.e. dirt = 130%)
      //We push the contributions back to h_muon_universes to be used in the plotting



      //Muon
      ////////////////////////////////////
      std::vector<TH1D*> h_muon_universes;

      h_muon_xsec[0][j] =  Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_overlay_muon[j][0],h_overlay_muon[j][1],h_bnb_muon[j],
					      h_muon_matrices_CV[j],muon_iter[j],h_muon_num_CV[j], h_muon_denom_CV[j],
					      muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_CV",var[j]),directory_name,false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
	h_muon_xsec[k+1][j] = Xsec.cross_section(h_ext_muon[j],h_dirt_muon[j],h_muon_total[k][j],h_muon_cc2p[k][j],h_bnb_muon[j],
					       h_muon_matrices[k][j],muon_iter[j],h_muon_num[k][j], h_muon_denom[k][j],
					       muon_xsec_max[j],Form("True Muon %s",var_titles[j]),Form("_muon%s_%d",var[j],k),directory_name,false,true); //CV Montecarlo values
	h_muon_universes.push_back(h_muon_xsec[k+1][j]); //h_muon_xsec[0][j] should be same as h_muon_xsec_CV[j]
      }

      //multisim
      plot_multisims(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]), h_muon_universes, h_muon_xsec[0][j], ymin_muon[j], muon_xsec_max[j], universes, directory_name);

      
      //Covariance shit
      h_2D_muon[j] = make_covariance_matrix(h_muon_xsec[0][j],h_muon_universes, universes, Form("muon_%s",var0[j]),true);
      plot_covariance(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_2D_muon[j],h_2D_muon[j],h_muon_xsec[0][j],directory_name);
      h_2D_muon[j]->Write();
      
      //Getting the Error
      h_muon_error[j] = make_error_histogram(Form("_muon%s",var0[j]),h_muon_xsec[0][j],h_2D_muon[j]);
      plot_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_error[j],h_muon_error[j],directory_name);
      h_muon_error[j]->Write();
      
      //Getting the fractional error
      h_muon_fractional_error[j] = make_fractional_error_histogram(Form("_muon%s",var0[j]),h_muon_xsec[0][j],h_muon_error[j]);
      plot_fractional_error(Form("_muon%s",var0[j]), Form("Muon %s",var_titles[j]),canv_muon_fractional_error[j],h_muon_fractional_error[j],directory_name);
      h_muon_fractional_error[j]->Write();


      //Leading
      ////////////////////////////////////
      std::vector<TH1D*> h_leading_universes;

      h_leading_xsec[0][j] =  Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_overlay_leading[j][0],h_overlay_leading[j][1],h_bnb_leading[j],
					      h_leading_matrices_CV[j],leading_iter[j],h_leading_num_CV[j], h_leading_denom_CV[j],
					      leading_xsec_max[j],Form("True Leading %s",var_titles[j]),Form("_leading%s_CV",var[j]),directory_name,false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
	h_leading_xsec[k+1][j] = Xsec.cross_section(h_ext_leading[j],h_dirt_leading[j],h_leading_total[k][j],h_leading_cc2p[k][j],h_bnb_leading[j],
					       h_leading_matrices[k][j],leading_iter[j],h_leading_num[k][j], h_leading_denom[k][j],
					       leading_xsec_max[j],Form("True Leading %s",var_titles[j]),Form("_leading%s_%d",var[j],k),directory_name,false,true); //CV Montecarlo values
	h_leading_universes.push_back(h_leading_xsec[k+1][j]); //h_leading_xsec[0][j] should be same as h_leading_xsec_CV[j]
      }

      //multisim
      plot_multisims(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]), h_leading_universes, h_leading_xsec[0][j], ymin_leading[j], leading_xsec_max[j], universes, directory_name);

      
      //Covariance shit
      h_2D_leading[j] = make_covariance_matrix(h_leading_xsec[0][j],h_leading_universes, universes, Form("leading_%s",var0[j]),true);
      plot_covariance(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_2D_leading[j],h_2D_leading[j],h_leading_xsec[0][j],directory_name);
      h_2D_leading[j]->Write();
      
      //Getting the Error
      h_leading_error[j] = make_error_histogram(Form("_leading%s",var0[j]),h_leading_xsec[0][j],h_2D_leading[j]);
      plot_error(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_leading_error[j],h_leading_error[j],directory_name);
      h_leading_error[j]->Write();
      
      //Getting the fractional error
      h_leading_fractional_error[j] = make_fractional_error_histogram(Form("_leading%s",var0[j]),h_leading_xsec[0][j],h_leading_error[j]);
      plot_fractional_error(Form("_leading%s",var0[j]), Form("Leading %s",var_titles[j]),canv_leading_fractional_error[j],h_leading_fractional_error[j],directory_name);
      h_leading_fractional_error[j]->Write();


      //Recoil
      ////////////////////////////////////
      std::vector<TH1D*> h_recoil_universes;

      h_recoil_xsec[0][j] =  Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_overlay_recoil[j][0],h_overlay_recoil[j][1],h_bnb_recoil[j],
					      h_recoil_matrices_CV[j],recoil_iter[j],h_recoil_num_CV[j], h_recoil_denom_CV[j],
					      recoil_xsec_max[j],Form("True Recoil %s",var_titles[j]),Form("_recoil%s_CV",var[j]),directory_name,false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
	h_recoil_xsec[k+1][j] = Xsec.cross_section(h_ext_recoil[j],h_dirt_recoil[j],h_recoil_total[k][j],h_recoil_cc2p[k][j],h_bnb_recoil[j],
					       h_recoil_matrices[k][j],recoil_iter[j],h_recoil_num[k][j], h_recoil_denom[k][j],
					       recoil_xsec_max[j],Form("True Recoil %s",var_titles[j]),Form("_recoil%s_%d",var[j],k),directory_name,false,true); //CV Montecarlo values
	h_recoil_universes.push_back(h_recoil_xsec[k+1][j]); //h_recoil_xsec[0][j] should be same as h_recoil_xsec_CV[j]
      }

      //multisim
      plot_multisims(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]), h_recoil_universes, h_recoil_xsec[0][j], ymin_recoil[j], recoil_xsec_max[j], universes, directory_name);

      
      //Covariance shit
      h_2D_recoil[j] = make_covariance_matrix(h_recoil_xsec[0][j],h_recoil_universes, universes, Form("recoil_%s",var0[j]),true);
      plot_covariance(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_2D_recoil[j],h_2D_recoil[j],h_recoil_xsec[0][j],directory_name);
      h_2D_recoil[j]->Write();
      
      //Getting the Error
      h_recoil_error[j] = make_error_histogram(Form("_recoil%s",var0[j]),h_recoil_xsec[0][j],h_2D_recoil[j]);
      plot_error(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_recoil_error[j],h_recoil_error[j],directory_name);
      h_recoil_error[j]->Write();
      
      //Getting the fractional error
      h_recoil_fractional_error[j] = make_fractional_error_histogram(Form("_recoil%s",var0[j]),h_recoil_xsec[0][j],h_recoil_error[j]);
      plot_fractional_error(Form("_recoil%s",var0[j]), Form("Recoil %s",var_titles[j]),canv_recoil_fractional_error[j],h_recoil_fractional_error[j],directory_name);
      h_recoil_fractional_error[j]->Write();
    }
      

    //Other Variables
    //////////////////////
    for(int j=0; j <num_other_var; j++){

      std::vector<TH1D*> h_other_universes;

      h_other_xsec[0][j] =  Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_overlay_other[j][0],h_overlay_other[j][1],h_bnb_other[j],
					       h_other_matrices_CV[j],other_iter[j],h_other_num_CV[j], h_other_denom_CV[j],
					       other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("%s_CV",other_var[j]),directory_name,false,true); //CV Montecarlo values

      for(int k = 0; k < universes; k++){
	h_other_xsec[k+1][j] = Xsec.cross_section(h_ext_other[j],h_dirt_other[j],h_other_total[k][j],h_other_cc2p[k][j],h_bnb_other[j],
						  h_other_matrices[k][j],other_iter[j],h_other_num[k][j], h_other_denom[k][j],
						  other_xsec_max[j],Form("True %s",other_var_titles[j]),Form("%s_%d",other_var[j],k),directory_name,false,true); //CV Montecarlo values
	h_other_universes.push_back(h_other_xsec[k+1][j]); //h_other_xsec[0][j] should be same as h_other_xsec_CV[j]
      }

      //multisim
      plot_multisims(Form("%s",other_var[j]), Form("%s",other_var_titles[j]), h_other_universes, h_other_xsec[0][j], ymin_other[j], other_xsec_max[j], universes, directory_name);

      
      //Covariance shit
      h_2D_other[j] = make_covariance_matrix(h_other_xsec[0][j],h_other_universes, universes, Form("%s",other_var[j]),true);
      plot_covariance(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_2D_other[j],h_2D_other[j],h_other_xsec[0][j],directory_name);
      h_2D_other[j]->Write();
      
      //Getting the Error
      h_other_error[j] = make_error_histogram(Form("%s",other_var[j]),h_other_xsec[0][j],h_2D_other[j]);
      plot_error(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_other_error[j],h_other_error[j],directory_name);
      h_other_error[j]->Write();
      
      //Getting the fractional error
      h_other_fractional_error[j] = make_fractional_error_histogram(Form("%s",other_var[j]),h_other_xsec[0][j],h_other_error[j]);
      plot_fractional_error(Form("%s",other_var[j]), Form("%s",other_var_titles[j]),canv_other_fractional_error[j],h_other_fractional_error[j],directory_name);
      h_other_fractional_error[j]->Write();

    } //end of loop over other var

    tfile_mine->Close();

  } //end of loop over directories
}//end of program
