#define analysis_cxx
#include "analysis.h"


void analysis::main(){

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
  TFile *f_overlay=new TFile(Form("../root_files/pelee/%s/histograms_%s_overlay_wgt.root",run,sample));//overlay histograms.
  TFile *f_bnb=new TFile(Form("../root_files/pelee/%s/histograms_%s_bnb.root",run,sample));//bnb histograms  
  TFile *f_ext=new TFile(Form("../root_files/pelee/%s/histograms_%s_ext.root",run,sample));//extbnb histograms 
  TFile *f_dirt=new TFile(Form("../root_files/pelee/%s/histograms_%s_dirt_wgt.root",run,sample));//dirt histograms    
  TFile* f_eff = new TFile(Form("../root_files/pelee/%s/histograms_mc_eff.root",run)); //All efficiency histograms  
  TFile *f_mom_thresholds =new TFile(Form("../root_files/pelee/%s/efficiency1.root",run)); //All efficiency histograms

  //Define Parameters and color scheme
  ////////////////////////////////////
  Define_Parameters(run, sample);
  
  //Color Scheme
  tolcols::init();
  //Color_t colors[] = {0,9031, 9030, 9029, 9028, 9026, 9025, 9024, 9032, kGray+2 ,9027,kViolet+3};  //Black, light pink, light orange, light yellow, olive, mint, light cyan, light blue, light grey, darker grey, olive
  //Color_t colors_raquel[] = {0,9012,9011,9010,9009,9008,9007,kGray+2,9032,9027}; //black, magenta, red, orange, mint, cyan, blue, dark gray, light gray, ccnue     
  Color_t colors[] = {0,9031, 9025, 9030, 9024, 9029, kViolet-8, 9028, 9032, 9026, kGray+2}; //black, //light pink,light cyan,light orange, light blue, light yellow, purple, olive, light gray, mint, darker gray

  //Color_t colors[] = {0,9031, 9030, 9024, kViolet-8, 9029, 9032, 9028, kGray+2, 9026, 9025};
  Color_t colors_raquel[] = {0,9012,9007,9011,kViolet-8,9010, 9032, 9026};     
  Color_t colors_chi2[11] = {kBlack,kRed,kRed+2,kBlue,kBlue+2,kGreen,kGreen+1,kYellow,kYellow-3,kOrange+8,204};

  /////////////////////////////////////////////////////////////
  //THIRD: MAKE DIRECTORY WITH TODAY'S DATE TO STORE ALL IMAGES
  //////////////////////////////////////////////////////////////
  const char* pathname = Form("images/%s/%s/",sample,run);
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
  
  ///////////////////////////////////
  //Plots of the Truth Variables
  ////////////////////////////////////
  /*  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_truth; j++){
      
      canv_truth[i][j] = new TCanvas(Form("C_truth%s%s",truth[j],cut[i]),Form("C_truth%s%s",truth[j],cut[i]),2000,1500);
      h_truth[i][j] = new THStack(Form("h_truth%s%s",truth[j],cut[i]),Form("h_truth%s%s",truth[j],cut[i]));
      h_truth[i][j]->Draw("HIST");
      h_truth[i][j]->SetTitle(Form(" ; %s ; Number of Events",truth_titles[j]));
      h_truth[i][j]->SetMaximum(ylim_truth[run_num][i][j]);
      //h_truth[i][j]->GetXaxis()->SetTitle(Form("%s",truth_titles[j]));

      if(plot_ccnue_truth){
	z = num_channels;
	f = num_channels;
      }else{
	z = num_channels - 1;
	f = num_channels - 1;
      }
      
      for(int k=1; k < z ; k++){
	h_mc[i][j][k]->SetLineColor(colors[k]);
	h_mc[i][j][k]->SetFillColor(colors[k]);
	h_mc[i][j][k]->SetLineWidth(1);
	h_truth[i][j]->Add(h_mc[i][j][k]);
      }
      
      if(plot_total_truth){
      h_mc[i][j][0]->Draw("e1SAME");
      }

      legend_truth[i][j] = new TLegend(0.71, 0.54, 0.899, 0.89);
      for(int k =1; k < z; k++){	
	legend_truth[i][j]->AddEntry(h_mc[i][j][f-k],Form("%s",channel_legend[f-k]),"f");

      }
      legend_truth[i][j]->SetLineWidth(0);
      legend_truth[i][j]->SetFillColor(kWhite);
      legend_truth[i][j]->SetTextSize(0.03);
      legend_truth[i][j]->Draw("same");
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{%s: %s}",cut_titles[i],truth_titles[j]));
      t->DrawLatex(0.195,0.92,Form("%s",pot_num));
      t->DrawLatex(0.82,0.92,Form("%s",sample_name));
      
      canv_truth[i][j]->Print(Form("%s%s%s.png",path.c_str(),truth[j],cut[i]));
      canv_truth[i][j]->Print(Form("%s%s%s.pdf",path.c_str(),truth[j],cut[i]));
      
    }
  }
  
  //////////////////////////
  //Plots of the PFP Stuff:
  //////////////////////////
  for(int i = 0; i < num_group; i ++){
    
      canv_pfp[i] = new TCanvas(Form("C_pfp_%s",group[i]),Form("C_pfp_%s",group[i]),2000,1500);
      h_pfp[i] = new THStack(Form("h_pfp%s",group[i]),Form("h_pfp%s",group[i]));
      h_pfp[i]->Add(h_overlay_pfp[i]);
      h_pfp[i]->Add(h_ext_pfp[i]);
      h_pfp[i]->Draw("hist");
      h_pfp[i]->GetYaxis()->SetTitle("Number of Events");
      h_pfp[i]->GetXaxis()->SetTitle(Form("%s",titles_pfp[i]));
      h_pfp[i]->SetTitle(Form("%s",titles_pfp[i]));
      h_pfp[i]->SetMaximum(ylim_pfp[i]);

      h_ext_pfp[i]->SetFillColor(kRed);
      h_overlay_pfp[i]->SetFillColor(kBlue);

      h_bnb_pfp[i]->Draw("e1SAME");
      h_bnb_pfp[i]->SetLineColor(kBlack);
      h_bnb_pfp[i]->SetLineWidth(1);

      legend_pfp[i] = new TLegend(0.63,0.58,1.0,0.89);
      legend_pfp[i]->AddEntry(h_bnb_pfp[i], "BNB", "lepf");
      legend_pfp[i]->AddEntry(h_ext_pfp[i], "EXTBNB", "f");
      legend_pfp[i]->AddEntry(h_overlay_pfp[i], "Overlay MC", "f");
      legend_pfp[i]->Draw("same");
      
      canv_pfp[i]->Print(Form("%s_%s.png",path.c_str(),group[i]));
      canv_pfp[i]->Print(Form("%s_%s.pdf",path.c_str(),group[i]));
  }
  */
  
  //////////////////////////////////
  //Plot the 2D Histograms:
  //////////////////////////////////
  /*  for(int i =0; i < num_group2d; i++){
    canv_2d[i] = new TCanvas(Form("C_2D_%s",group2d[i]),Form("C_2D_%s",group2d[i]),2000,1500);
    h_overlay2D[i]->Draw("colz"); 
    canv_2d[i]->Print(Form("%s_%s.png",path.c_str(),group2d[i]));
    canv_2d[i]->Print(Form("%s_%s.pdf",path.c_str(),group2d[i]));
  }
  */
  
  /////////////////////////////
  //Plot the Efficiency Stuff:
  //////////////////////////
  efficiency.Plot_Efficiency(path);
  
  ///////////////////////////////////
  //Plots of the Reconstructed Vertex
  ///////////////////////////////////
  vertex.Plot_Vertex(run,colors,colors_raquel,path);

  
  /* 
  //////////////////////////
  //Track plots, such as PID
  //////////////////////////
  for(int i=0; i < num_track; i++){
    for(int k=0; k < track_cut; k++){
      for(int j=0; j < num_particles; j ++){
	h_track_overlay_vec.push_back(h_track_overlay[i][k][j]);
      }
    
      Plot_Histograms(run,pot_num,sample_name,colors_chi2, h_track_overlay_vec, h_track_overlay0[i][k][0],h_track_overlay0[i][k][1],h_track_overlay0[i][k][2],
		      h_track_ext[i][k][0], h_track_ext[i][k][1], h_track_ext[i][k][2],
		      h_track_dirt[i][k][0],h_track_dirt[i][k][1],h_track_dirt[i][k][2],
		      h_track_bnb[i][k][0],h_track_bnb[i][k][1],canv_track[i][k], h_track[i][k],
		      pad_track[i][k], pad0_track[i][k], legend_track[i][k], channel_legend_chi2,
		      ymax_track[i], ymin_track[i], num_particles, Form("%s",titles_track[i]), path, a_track[i],Form(""), Form("%s%s",variable[i],which_track_cut[k]), Form(""),false, true, true, 0.0,0.19,0.5,1.5, xlim_track[i]);
      h_track_overlay_vec.clear();
    }
  }

  */
  
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
    Plot_Histograms(run,pot_num,sample_name,colors, h_muon_overlay_vec, h_muon_overlay0[i][0],h_muon_overlay0[i][1],h_muon_overlay0[i][2],h_muon_stat[i],h_muon_systematic[i],h_covar_muon[i],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon[i], h_muon[i],
		    pad_muon[i], pad0_muon[i], legend_muon[i],channel_legend,
		    muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon%s",var[i]), "",false, true,flip_muon[i],0.0,0.19,0.0,muon_bnb_ymin[run_num][i]);
    h_muon_overlay_vec.clear();
       
    //muon:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_muon_overlay_raquel_vec, h_muon_overlay0_raquel[i][0],h_muon_overlay0_raquel[i][1],h_muon_overlay0_raquel[i][2],h_muon_stat[i],h_muon_systematic[i],h_covar_muon[i],
		    h_muon_ext[i][0], h_muon_ext[i][1], h_muon_ext[i][2],
		    h_muon_dirt[i][0],h_muon_dirt[i][1],h_muon_dirt[i][2],
		    h_muon_bnb[i][0],h_muon_bnb[i][1],canv_muon_raquel[i], h_muon_raquel[i],
		     pad_muon_raquel[i], pad0_muon_raquel[i], legend_muon_raquel[i], channel_legend_raquel,
		     muon_ylim[run_num][i], muon_ymin[run_num][i], num_channels_raquel, Form("Muon: %s",titles_var[i]), path,  a1,"", Form("_muon_raquel%s",var[i]), "",false, true,flip_muon[i],0.0,0.19,0.0,muon_bnb_ymin[run_num][i]);
    h_muon_overlay_raquel_vec.clear();
    
    //leading proton:mine
    Plot_Histograms(run,pot_num,sample_name,colors, h_leading_overlay_vec, h_leading_overlay0[i][0],h_leading_overlay0[i][1],h_leading_overlay0[i][2],h_leading_stat[i],h_leading_systematic[i],h_covar_leading[i],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading[i], h_leading[i],
		    pad_leading[i], pad0_leading[i], legend_leading[i],channel_legend,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading%s",var[i]), "",false, true, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_vec.clear();

    //leading proton:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_leading_overlay_raquel_vec, h_leading_overlay0_raquel[i][0],h_leading_overlay0_raquel[i][1],h_leading_overlay0_raquel[i][2],h_leading_stat[i],h_leading_systematic[i],h_covar_leading[i],
		    h_leading_ext[i][0], h_leading_ext[i][1], h_leading_ext[i][2],
		    h_leading_dirt[i][0],h_leading_dirt[i][1],h_leading_dirt[i][2],
		    h_leading_bnb[i][0],h_leading_bnb[i][1],canv_leading_raquel[i], h_leading_raquel[i],
		    pad_leading_raquel[i], pad0_leading_raquel[i], legend_leading_raquel[i],channel_legend_raquel,
		    leading_ylim[run_num][i], leading_ymin[run_num][i], num_channels_raquel, Form("Leading Proton: %s",titles_var[i]), path,  a1,"", Form("_leading_raquel%s",var[i]), "",false, true, flip_lead[i],0.0,0.19,0.0,2.7);
    h_leading_overlay_raquel_vec.clear();
        
    //recoil proton:mine
    Plot_Histograms(run,pot_num,sample_name,colors, h_recoil_overlay_vec, h_recoil_overlay0[i][0],h_recoil_overlay0[i][1],h_recoil_overlay0[i][2],h_recoil_stat[i],h_recoil_systematic[i],h_covar_recoil[i],
		    h_recoil_ext[i][0], h_recoil_ext[i][1], h_recoil_ext[i][2],
		    h_recoil_dirt[i][0],h_recoil_dirt[i][1],h_recoil_dirt[i][2],
		    h_recoil_bnb[i][0],h_recoil_bnb[i][1],canv_recoil[i], h_recoil[i],
		    pad_recoil[i], pad0_recoil[i], legend_recoil[i],channel_legend,
		    recoil_ylim[run_num][i], recoil_ymin[run_num][i], num_channels, Form("Recoil Proton: %s",titles_var[i]), path,  a1,"", Form("_recoil%s",var[i]), "",false, true, flip_recoil[i],0.0105, 0.19, 0.0, 2.7);
    h_recoil_overlay_vec.clear();

    //recoil proton:raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_recoil_overlay_raquel_vec, h_recoil_overlay0_raquel[i][0],h_recoil_overlay0_raquel[i][1],h_recoil_overlay0_raquel[i][2],h_recoil_stat[i],h_recoil_systematic[i],h_covar_recoil[i],
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
    Plot_Histograms(run,pot_num,sample_name,colors, h_phys_overlay_vec, h_phys_overlay0[i][0],h_phys_overlay0[i][1],h_phys_overlay0[i][2],h_phys_stat[i],h_phys_systematic[i],h_covar_phys[i],
                    h_phys_ext[i][0], h_phys_ext[i][1], h_phys_ext[i][2],
                    h_phys_dirt[i][0],h_phys_dirt[i][1],h_phys_dirt[i][2],
                    h_phys_bnb[i][0],h_phys_bnb[i][1],canv_phys[i], h_phys[i],
                    pad_phys[i], pad0_phys[i], legend_phys[i],channel_legend,
                    phys_ylim[run_num][i], phys_ymin[run_num][i], num_channels, physics_titles[i], path,  a1,"", physics[i], "",false, true, false,0.0, 0.19,0.0,phys_bnb_ymin[run_num][i]);
    h_phys_overlay_vec.clear();
     
    //raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_phys_overlay_raquel_vec, h_phys_overlay0_raquel[i][0],h_phys_overlay0_raquel[i][1],h_phys_overlay0_raquel[i][2],h_phys_stat[i],h_phys_systematic[i],h_covar_phys[i],
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
    Plot_Histograms(run,pot_num,sample_name,colors, h_stv_overlay_vec, h_stv_overlay0[i][0],h_stv_overlay0[i][1],h_stv_overlay0[i][2],h_stv_stat[i],h_stv_systematic[i],h_covar_stv[i],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv[i], h_stv[i],
                    pad_stv[i], pad0_stv[i], legend_stv[i],channel_legend,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels, stv_titles[i], path, a1,"", stv[i], "",false, true,false,0.0,0.19,0.0,stv_bnb_ymin[run_num][i]);
    h_stv_overlay_vec.clear();
    
    //raquel
    Plot_Histograms(run,pot_num,sample_name,colors_raquel, h_stv_overlay_raquel_vec, h_stv_overlay0_raquel[i][0],h_stv_overlay0_raquel[i][1],h_stv_overlay0_raquel[i][2],h_stv_stat[i],h_stv_systematic[i],h_covar_stv[i],
                    h_stv_ext[i][0], h_stv_ext[i][1], h_stv_ext[i][2],
                    h_stv_dirt[i][0],h_stv_dirt[i][1],h_stv_dirt[i][2],
                    h_stv_bnb[i][0],h_stv_bnb[i][1],canv_stv_raquel[i], h_stv_raquel[i],
                    pad_stv_raquel[i], pad0_stv_raquel[i], legend_stv_raquel[i],channel_legend_raquel,
                    stv_ylim[run_num][i], stv_ymin[run_num][i], num_channels_raquel, stv_titles[i], path, a1,"", Form("%s_raquel",stv[i]), "",false, true,false,0.0,0.19,0.0,stv_bnb_ymin[run_num][i]);
    h_stv_overlay_raquel_vec.clear();
  }

  

} //end of program
