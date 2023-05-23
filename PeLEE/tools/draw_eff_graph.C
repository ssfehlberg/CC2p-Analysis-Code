#define draw_eff_graph_cxx

void draw_eff_graph(){

  //Stuff for Drawing
  //////////////////////
  TLatex* t = new TLatex(); //latex
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetHistMinimumZero(kFALSE);
  t->SetNDC();
  t->SetTextAlign(22);
  char const* pot_num="#scale[0.6]{Runs 1+2+3 Accumulated POT: 6.79e+20}";
  char const* sample_name="#scale[0.6]{MicroBooNE In-Progress}";


  //File and number of slices
  //////////////////////////
  TFile *f_eff1=new TFile("root_files/pelee/Run1/histograms_pelee_overlay_wgt.root"); //All efficiency histograms: Run 1
  TFile *f_eff2=new TFile("root_files/pelee/Run2/histograms_pelee_overlay_wgt.root"); //All efficiency histograms: Run 2
  TFile *f_eff3=new TFile("root_files/pelee/Run3/histograms_pelee_overlay_wgt.root"); //All efficiency histograms: Run 3

  std::vector<TFile*> files;
  files.push_back(f_eff1);
  files.push_back(f_eff2);
  files.push_back(f_eff3);
	 
  
  const char* run[3] = {"Run1","Run2","Run3"};
  TGraph* eff_graph[3];
  TGraph* pur_graph[3];
  TCanvas* canv_eff_pur[3];
  TLegend* lg[3];


  const char* eff_cut_label[8] = {"No Cuts","FV","3 PFP","Trk Scr","Vtx","PID","Containment","Reco. Mom Cut"};

  
  //efficiency and purity as function of cuts
  for(int i =0; i < 3; i++){

    eff_graph[i] = (TGraph*)files[i]->Get("eff_graph");
    pur_graph[i] = (TGraph*)files[i]->Get("pur_graph");
    
    canv_eff_pur[i] = new TCanvas(Form("canv_eff_pur%s",run[i]),Form("canv_eff_pur%s",run[i]),2000,1500);
    eff_graph[i]->Draw("alp");
    eff_graph[i]->SetLineColor(kBlue);
    eff_graph[i]->SetMarkerColor(kBlue);
    eff_graph[i]->SetLineWidth(3);
    eff_graph[i]->SetMarkerStyle(20);
    eff_graph[i]->SetTitle("Efficiency and Purity as a Function of Cuts");
    eff_graph[i]->GetYaxis()->SetRangeUser(0.,1.);
    eff_graph[i]->GetYaxis()->SetTitle("Efficiency or Purity");

    TAxis* xax = eff_graph[i]->GetXaxis();
    for(int j = 1; j < xax->GetXmax()-1; j++){
      int bin_index_up = xax->FindBin(j+1);
      int bin_index_low = xax->FindBin(j);
      xax->SetBinLabel(bin_index_low,Form("%s", eff_cut_label[j-1]));
      xax->SetBinLabel(bin_index_up,Form("%s", eff_cut_label[j]));
    }
    xax->Draw("hist");
    
    pur_graph[i]->Draw("SAME lp");
    pur_graph[i]->SetLineColor(kRed);
    pur_graph[i]->SetMarkerColor(kRed);
    pur_graph[i]->SetLineWidth(3);
    pur_graph[i]->SetMarkerStyle(20);
    
    lg[i] = new TLegend(0.65,0.65,0.85,0.85);
    lg[i]->AddEntry(eff_graph[i],"Efficiency","lp");
    lg[i]->AddEntry(pur_graph[i],"Purity","lp");
    lg[i]->Draw("same");
      
    canv_eff_pur[i]->Print(Form("images/pelee/%s/_efficiency_and_purity.png",run[i]));
    canv_eff_pur[i]->Print(Form("images/pelee/%s/_efficiency_and_purity.pdf",run[i]));
  }


}//end of code
