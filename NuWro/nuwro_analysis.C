#define nuwro_analysis_cxx
#include "nuwro_analysis.h"


void nuwro_analysis::main(){

  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  //FIRST: CHECK WHICH SAMPLE THIS IS: PELEE,FILTERED,OR UNFILTERED & WHAT RUN: JAN, RUN 1,RUN 2,RUN or 3
  /////////////////////////////////////////////////////////////////////////////////////////////////////////
  const char* sample = which_sample();
  std::pair<const char*, const char*> r = which_run();
  const char* run = r.first;
  const char* run_title = r.second;

  ///////////////////////////////////////////////////////
  //SECOND: DEFINE HISTOGRAMS FILES AND GENERAL VARIABLES
  //////////////////////////////////////////////////////
  TFile *f_overlay=new TFile(Form("../root_files/nuwro/%s/histograms_nuwro%s_overlay_wgt.root",run,sample));//overlay histograms.
  TFile *f_bnb=new TFile(Form("../root_files/pelee/%s/histograms_pelee%s_bnb.root",run,sample));//bnb histograms  
  TFile *f_ext=new TFile(Form("../root_files/pelee/%s/histograms_pelee%s_ext.root",run,sample));//extbnb histograms 
  TFile *f_dirt=new TFile(Form("../root_files/pelee/%s/histograms_pelee%s_dirt_wgt.root",run,sample));//dirt histograms    
  TFile* f_eff = new TFile(Form("../root_files/pelee/%s/histograms_mc_eff.root",run)); //All efficiency histograms  
  TFile *f_mom_thresholds =new TFile(Form("../root_files/pelee/%s/efficiency1.root",run)); //All efficiency histograms

  //Define Parameters and color scheme
  ////////////////////////////////////
  Define_Parameters(run, sample);
  
  //Color Scheme
  tolcols::init();
  Color_t colors[] = {0,9031, 9025, 9030, 9024, 9029, kViolet-8, 9028, 9032, 9026, kGray+2};
  Color_t colors_raquel[] = {0,9012,9007,9011,kViolet-8,9010,9001,9026,9032};
  Color_t colors_chi2[11] = {kBlack,kRed,kRed+2,kBlue,kBlue+2,kGreen,kGreen+1,kYellow,kYellow-3,kOrange+8,204};

  /////////////////////////////////////////////////////////////
  //THIRD: MAKE DIRECTORY WITH TODAY'S DATE TO STORE ALL IMAGES
  //////////////////////////////////////////////////////////////
  const char* pathname = Form("images/nuwro%s/%s/",sample,run);
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

  /////////////////////////////////////////
  //GRAB ALL THE HISTOGRAMS FROM THE FILES
  ////////////////////////////////////////

  std::cout<<"Grabbing histograms"<<std::endl;
  Grab_Histograms(f_overlay,f_bnb,f_ext,f_dirt,f_eff,f_mom_thresholds);
  std::cout<<"After Grabbing histograms"<<std::endl;
  
  //////////////////////////////
  //Particle specific plots
  /////////////////////////////
  for(int i = 0; i < num_var; i++){
    for(int k=0; k < num_channels; k++){
      h_muon_overlay_vec.push_back(h_muon_overlay[i][k]);
      h_recoil_overlay_vec.push_back(h_recoil_overlay[i][k]);
      h_leading_overlay_vec.push_back(h_leading_overlay[i][k]);
    }
    for(int k=0; k < num_channels_raquel; k++){
      h_muon_overlay_raquel_vec.push_back(h_muon_overlay_raquel[i][k]);
      h_recoil_overlay_raquel_vec.push_back(h_recoil_overlay_raquel[i][k]);
      h_leading_overlay_raquel_vec.push_back(h_leading_overlay_raquel[i][k]);
      
    }
  
    //muon:mine
    Plot_Histograms(run,pot_num,sample_name,colors, h_muon_overlay_vec, h_muon_overlay0[i][0],h_muon_overlay0[i][1],h_muon_overlay0[i][2],h_muon_stat[i],h_muon_systematic[i],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon[i], h_muon[i],
		    pad_muon[i], pad0_muon[i], legend_muon[i],channel_legend,
		    muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon%s",var[i]), "",false, true,flip_muon[i],0.0,0.19,0.0,muon_bnb_ymin[run_num][i]);
    h_muon_overlay_vec.clear();
       
    //muon:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_muon_overlay_raquel_vec, h_muon_overlay0_raquel[i][0],h_muon_overlay0_raquel[i][1],h_muon_overlay0_raquel[i][2],h_muon_stat[i],h_muon_systematic[i],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon_raquel[i], h_muon_raquel[i],
		     pad_muon_raquel[i], pad0_muon_raquel[i], legend_muon_raquel[i], channel_legend_raquel,
		     muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels_raquel, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon_raquel%s",var[i]), "",false, true,flip_muon[i],0.0,0.19,0.0,muon_bnb_ymin[run_num][i]);
    h_muon_overlay_raquel_vec.clear();
    
    //leading proton:mine
    Plot_Histograms(run,pot_num,sample_name,colors, h_leading_overlay_vec, h_leading_overlay0[i][0],h_leading_overlay0[i][1],h_leading_overlay0[i][2],h_leading_stat[i],h_leading_systematic[i],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading[i], h_leading[i],
		    pad_leading[i], pad0_leading[i], legend_leading[i],channel_legend,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading%s",var[i]), "",false, true, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_vec.clear();

    //leading proton:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_leading_overlay_raquel_vec, h_leading_overlay0_raquel[i][0],h_leading_overlay0_raquel[i][1],h_leading_overlay0_raquel[i][2],h_leading_stat[i],h_leading_systematic[i],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading_raquel[i], h_leading_raquel[i],
		    pad_leading_raquel[i], pad0_leading_raquel[i], legend_leading_raquel[i],channel_legend_raquel,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels_raquel, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading_raquel%s",var[i]), "",false, true, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_raquel_vec.clear();
        
    //recoil proton:mine
    Plot_Histograms(run,pot_num,sample_name,colors, h_recoil_overlay_vec, h_recoil_overlay0[i][0],h_recoil_overlay0[i][1],h_recoil_overlay0[i][2],h_recoil_stat[i],h_recoil_systematic[i],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil[i], h_recoil[i],
		    pad_recoil[i], pad0_recoil[i], legend_recoil[i],channel_legend,
		    recoil_ylim[run_num][i], recoil_ymin[run_num][i], num_channels, Form("Recoil Proton: %s",titles_var[i]), path,  a1,"", Form("_recoil%s",var[i]), "",false, true, flip_recoil[i],0.0105, 0.19, 0.0, 2.7);
    h_recoil_overlay_vec.clear();

    //recoil proton:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_recoil_overlay_raquel_vec, h_recoil_overlay0_raquel[i][0],h_recoil_overlay0_raquel[i][1],h_recoil_overlay0_raquel[i][2],h_recoil_stat[i],h_recoil_systematic[i],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil_raquel[i], h_recoil_raquel[i],
		    pad_recoil_raquel[i], pad0_recoil_raquel[i], legend_recoil_raquel[i], channel_legend_raquel,
		    recoil_ylim[run_num][i], recoil_ymin[run_num][i], num_channels_raquel, Form("Recoil Proton: %s",titles_var[i]), path,  a1,"", Form("_recoil_raquel%s",var[i]), "",false, true,flip_recoil[i],0.0,0.19,0.0,2.7);
    h_recoil_overlay_raquel_vec.clear();
    
  }

  
  //PHYSICS PLOTS!!!
  ////////////////////////////////
   for(int i=0; i < num_phys; i++){
    for(int j=0; j < num_channels; j++){
      h_phys_overlay_vec.push_back(h_phys_overlay[i][j]);
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_phys_overlay_raquel_vec.push_back(h_phys_overlay_raquel[i][j]);
      h_phys_overlay_raquel[i][j]->Draw("hist");
    }

    //mine 
    Plot_Histograms(run,pot_num,sample_name,colors, h_phys_overlay_vec, h_phys_overlay0[i][0],h_phys_overlay0[i][1],h_phys_overlay0[i][2],h_phys_stat[i],h_phys_systematic[i],
                    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys[i], h_phys[i],
                    pad_phys[i], pad0_phys[i], legend_phys[i],channel_legend,
                    phys_ylim[run_num][i], phys_ymin[run_num][i], num_channels, physics_titles[i], path,  a1,"", physics[i], "",false, true, false,0.0, 0.19,0.0,phys_bnb_ymin[run_num][i]);
    h_phys_overlay_vec.clear();
     
    //raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_phys_overlay_raquel_vec, h_phys_overlay0_raquel[i][0],h_phys_overlay0_raquel[i][1],h_phys_overlay0_raquel[i][2],h_phys_stat[i],h_phys_systematic[i],
		    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys_raquel[i], h_phys_raquel[i],
                    pad_phys_raquel[i], pad0_phys_raquel[i], legend_phys_raquel[i],channel_legend_raquel,
                    phys_ylim[run_num][i], phys_ymin[run_num][i], num_channels_raquel, physics_titles[i], path,  a1,"", Form("%s_raquel",physics[i]), "",false, true, false,0.0,0.19,0.0,phys_bnb_ymin[run_num][i]);
		    h_phys_overlay_raquel_vec.clear();
		    
  }
  
   
  //STV PLOTS
  ///////////////////////////////////////////
  for(int i=0; i < num_stv; i++){
    for(int j=0; j < num_channels; j++){   
      h_stv_overlay_vec.push_back(h_stv_overlay[i][j]);
    }
    for(int j=0; j < num_channels_raquel; j++){
      h_stv_overlay_raquel_vec.push_back(h_stv_overlay_raquel[i][j]);
    }

    //mine 
    Plot_Histograms(run,pot_num,sample_name,colors, h_stv_overlay_vec, h_stv_overlay0[i][0],h_stv_overlay0[i][1],h_stv_overlay0[i][2],h_stv_stat[i],h_stv_systematic[i],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv[i], h_stv[i],
                    pad_stv[i], pad0_stv[i], legend_stv[i],channel_legend,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels, stv_titles[i], path, a1,"", stv[i], "",false, true,false,0.0,0.19,0.0,stv_bnb_ymin[run_num][i]);
    h_stv_overlay_vec.clear();
    
    //raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_stv_overlay_raquel_vec, h_stv_overlay0_raquel[i][0],h_stv_overlay0_raquel[i][1],h_stv_overlay0_raquel[i][2],h_stv_stat[i],h_stv_systematic[i],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv_raquel[i], h_stv_raquel[i],
                    pad_stv_raquel[i], pad0_stv_raquel[i], legend_stv_raquel[i],channel_legend_raquel,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels_raquel, stv_titles[i], path, a1,"", Form("%s_raquel",stv[i]), "",false, true,false,0.0,0.19,0.0,stv_bnb_ymin[run_num][i]);
    h_stv_overlay_raquel_vec.clear();
  }

  

} //end of program
