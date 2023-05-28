class plotting_tools{

 public:
  virtual void plot_multisims(const char* variable, const char* title, TCanvas* canv, TLegend* legend, std::vector<TH1D*> hist, TH1D* hist_CV, double ymin, double ymax, int num_universes, const char* directory);
  virtual void plot_covariance(const char* variable, const char* title, TCanvas* canv, TH2D* hist, TH1D* hist_CV, const char* directory);
  virtual void plot_error(const char* variable, const char* title, TCanvas* canv, TH1D* hist_error, const char* directory);
  virtual void plot_fractional_error(const char* variable, const char* title, TCanvas* canv, TH1D* hist_error, const char* directory);
  virtual TH1D* total_error(const char* name,std::vector<TH1D*> hist);
  virtual void plot_detvar_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total); //same as below
  virtual void plot_GENIE_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total); //same as above
  virtual void plot_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total);
  virtual void make_total_covariance(const char* variable, const char* title,std::vector<TH2D*> hist_2D, const char* directory);

}; //end of class definition

void plotting_tools::plot_multisims(const char* variable, const char* title, TCanvas* canv, TLegend* legend, std::vector<TH1D*> hist, TH1D* hist_CV, double ymin, double ymax, int num_universes, const char* directory){

  canv = new TCanvas(Form("canv%s",variable),Form("canv%s",variable),2000,1500);

  hist[0]->Draw("hist");
  hist[0]->SetLineColor(kGreen+1);
  hist[0]->SetLineWidth(3);
  hist[0]->SetTitle(Form("%s",title)); //title
  hist[0]->SetXTitle(Form("%s",title)); //xtitle
  hist[0]->GetXaxis()->SetTitleSize(40); 
  hist[0]->GetXaxis()->SetTitleFont(43);
  hist[0]->GetXaxis()->SetTitleOffset(1.5);
  hist[0]->GetXaxis()->SetLabelFont(43);
  hist[0]->GetXaxis()->SetLabelSize(30);
  hist[0]->SetYTitle("Number of Events"); //Y title
  hist[0]->GetYaxis()->SetTitleSize(40);
  hist[0]->GetYaxis()->SetTitleFont(43);
  hist[0]->GetYaxis()->SetTitleOffset(1.5);
  hist[0]->GetYaxis()->SetLabelFont(43);
  hist[0]->GetYaxis()->SetLabelSize(30);
  hist[0]->SetMaximum(ymax); //max
  hist[0]->SetMinimum(ymin); //min

  //draw other universes
  for(int k=1; k < num_universes; k++){
    hist[k]->Draw("hist same");
    hist[k]->SetLineColor(kGreen+1);
    hist[k]->SetLineWidth(3);
  }

  //Draw CV
  hist_CV->Draw("hist same");
  hist_CV->SetLineColor(kBlack);
  hist_CV->SetLineWidth(3);
  
  legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  legend->AddEntry(hist_CV,Form("Central Value"),"L");
  legend->AddEntry(hist[0],Form("Multisim Universe"),"L");
  legend->Draw("same");
  gStyle->SetOptStat(0);
  canv->Print(Form("images/%s/Multisims/%s.png",directory,variable));
  
}

void plotting_tools::plot_covariance(const char* variable, const char* title, TCanvas* canv, TH2D* hist, TH1D* hist_CV, const char* directory){

  canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  hist->Draw("colz text");
  
  hist->SetTitle(Form("Covariance Matrix: %s",title)); //title
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
    edges[i] = hist_CV->GetBinLowEdge(i+1);
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
  canv->Print(Form("images/%s/Covariance_Matrices/_2D%s.png",directory,variable));
    

}

void plotting_tools::plot_error(const char* variable, const char* title, TCanvas* canv, TH1D* hist_error, const char* directory){

 canv =  new TCanvas(Form("canv_error%s",variable),Form("canv_error%s",variable),2000,1500);
 hist_error->Draw("hist");
 hist_error->SetLineColor(kBlack);
 hist_error->SetLineWidth(2);
  
 hist_error->SetTitle(Form("Error: %s",title)); //title
 hist_error->SetXTitle(Form("%s",title));
 hist_error->GetXaxis()->SetTitleSize(40);
 hist_error->GetXaxis()->SetTitleFont(43);
 hist_error->GetXaxis()->SetTitleOffset(1.5);
 hist_error->GetXaxis()->SetLabelFont(43);
 hist_error->GetXaxis()->SetLabelSize(30);

 hist_error->SetYTitle("Number of Events");
 hist_error->GetYaxis()->SetTitleSize(40);
 hist_error->GetYaxis()->SetTitleFont(43);
 hist_error->GetYaxis()->SetTitleOffset(1.5);
 hist_error->GetYaxis()->SetLabelFont(43);
 hist_error->GetYaxis()->SetLabelSize(30);

 gStyle->SetOptStat(0);
 canv->Print(Form("images/%s/Error_Number_of_Events/_error%s.png",directory,variable));

}

void plotting_tools::plot_fractional_error(const char* variable, const char* title, TCanvas* canv, TH1D* hist_error, const char* directory){

 canv =  new TCanvas(Form("canv_fractional_error%s",variable),Form("canv_fractional_error%s",variable),2000,1500);
 hist_error->Draw("hist text");
 hist_error->SetLineColor(kBlack);
 hist_error->SetLineStyle(10);

 TH1D* h_clone = (TH1D*)hist_error->Clone();
 h_clone->Scale(-1.0);
 h_clone->Draw("same hist");
 h_clone->SetLineColor(kBlack);
 h_clone->SetLineStyle(4);
   
 hist_error->SetTitle(Form("Fractional Uncertainty: %s",title)); //title
 hist_error->SetXTitle(Form("%s",title));
 hist_error->GetXaxis()->SetTitleSize(40);
 hist_error->GetXaxis()->SetTitleFont(43);
 hist_error->GetXaxis()->SetTitleOffset(1.5);
 hist_error->GetXaxis()->SetLabelFont(43);
 hist_error->GetXaxis()->SetLabelSize(30);

 hist_error->SetYTitle("Fractional Uncertainty");
 hist_error->GetYaxis()->SetTitleSize(40);
 hist_error->GetYaxis()->SetTitleFont(43);
 hist_error->GetYaxis()->SetTitleOffset(1.5);
 hist_error->GetYaxis()->SetLabelFont(43);
 hist_error->GetYaxis()->SetLabelSize(30);

 hist_error->SetMinimum(-1.0);
 hist_error->SetMaximum(1.0);

 int nbins = hist_error->GetNbinsX();
 double xlow = hist_error->GetBinLowEdge(1);
 double xhigh = hist_error->GetBinLowEdge(nbins+1);
 TLine* a = new TLine(xlow,0,xhigh,0);
 a->Draw("same");
 a->SetLineColor(kRed);
 a->SetLineWidth(3);

 gStyle->SetOptStat(0);
 canv->Print(Form("images/%s/Fractional_Error/_fractional_error%s.png",directory,variable));

}

TH1D* plotting_tools::total_error(const char* name,std::vector<TH1D*> hist){

  int nbins = hist[0]->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist[0]->GetBinLowEdge(i+1);
  }

  TH1D* h_total_error = new TH1D(Form("hist_fractional_errors%s",name),Form("hist_fractional_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    double sum = 0;

    for(int j=0; j < hist.size(); j++){
      double error = std::pow(hist[j]->GetBinContent(i),2);
      sum += error;
    }

    double total_error = std::sqrt(sum + std::pow(0.02,2)); //applying the POT counting error
    h_total_error->SetBinContent(i,total_error);

  }

  return h_total_error;

}

void plotting_tools::plot_detvar_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total){

  //Stuff for Drawing
  //////////////////////
  Color_t colors[] = {kRed, kOrange+6, kYellow-3, kGreen+2, kCyan, kBlue, kViolet+1, kMagenta, kGray};
  int line_style[] = {1,9,1,9,1,9,1,9,1,9,1,9,1};
  static const int num_files = 9;
  const char* file_titles[num_files] = {"LY Attenuation","LY Down","LY Rayleigh",
					"ThetaXZ","ThetaYZ",
					"X","YZ",
					"Recombination","SCE"};
  //Now to plot
  //////////////
  canv = new TCanvas(Form("canv%s",name),Form("canv%s",name),2000,1500);
  if(flip_legend == true){
    legend = new TLegend(0.115, 0.65, 0.525, 0.89);
  } else {
    legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  }
  legend->SetNColumns(2);
  
  hist[0]->Draw("HIST TEXT00");
  hist[0]->SetLineColor(colors[0]);
  hist[0]->SetLineWidth(3);
  hist[0]->SetLineStyle(line_style[0]);
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
    hist[f]->SetLineStyle(line_style[f]);
    hist[f]->SetMarkerColor(colors[f]);
    legend->AddEntry(hist[f],Form("%s",file_titles[f]),"lpf");
  }
    
  h_total = total_error(Form("%s",name),hist);
  h_total->Draw("HIST TEXT00 SAME");
  h_total->SetLineColor(kBlack);
  h_total->SetLineWidth(5);
  h_total->SetLineStyle(1);
  h_total->SetMarkerColor(kBlack);
  h_total->Write();
  
  legend->AddEntry(h_total,"Total detVar Uncertainty","lpf");
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, In-Progress","C");
  legend->Draw("same");
  gStyle->SetOptStat(0);
  canv->Print(Form("images/detVar/Total/_total_error%s.png",name));
 
}

void plotting_tools::plot_GENIE_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total){

  //Stuff for Drawing
  //////////////////////
  Color_t colors[] = {kRed, kRed, kOrange+6, kOrange+6, kYellow-3, kYellow-3, kGreen+2 ,kGreen+2, kBlue, kBlue, kViolet+1, kViolet+1,kMagenta};
  int line_style[] = {1,9,1,9,1,9,1,9,1,9,1,9,1};
  static const int num_files = 13;
  const char* file_titles[num_files] = {"Multisims", "AxFFCCQE Shape","DecayAngleMEC","NormCCCOH",
					"NormNCCOH","RPA_CCQE",
					"ThetaDelta2NRad","ThetaDelta2Npi",
					"VecFFCCQEshape","XSecShape_CCMEC",
					"SCC: F^{3}_{V}","SCC: F^{3}_{A}","Rootino Fix"};
   
  //Now to plot
  //////////////
  canv = new TCanvas(Form("canv%s",name),Form("canv%s",name),2000,1500);
  if(flip_legend == true){
    legend = new TLegend(0.115, 0.65, 0.525, 0.89);
  } else {
    legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  }
  legend->SetNColumns(2);
  
  hist[0]->Draw("HIST TEXT00");
  hist[0]->SetLineColor(colors[0]);
  hist[0]->SetLineWidth(3);
  hist[0]->SetLineStyle(line_style[0]);
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
    hist[f]->SetLineStyle(line_style[f]);
    hist[f]->SetMarkerColor(colors[f]);
    legend->AddEntry(hist[f],Form("%s",file_titles[f]),"lpf");
  }
    
  h_total = total_error(Form("%s",name),hist);
  h_total->Draw("HIST TEXT00 SAME");
  h_total->SetLineColor(kBlack);
  h_total->SetLineWidth(5);
  h_total->SetLineStyle(1);
  h_total->SetMarkerColor(kBlack);
  h_total->Write();
  
  legend->AddEntry(h_total,"Total GENIE Uncertainty","lpf");
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, In-Progress","C");
  legend->Draw("same");
  gStyle->SetOptStat(0);
  canv->Print(Form("images/GENIE_Total/_total_error%s.png",name));
 
}

void plotting_tools::plot_total_error(const char* name,TCanvas* canv, bool flip_legend,TLegend* legend,std::vector<TH1D*> hist,const char* title, TH1D* h_total){

  //Stuff for Drawing
  /////////////////////
  static const int num_files = 7;
  const char* file_titles[num_files] = {"Cross-Section","Detector Variations","Flux","Reinteractions","Dirt","Iteration","Statistical"};
  Color_t colors[] = {kRed, kOrange+6, kYellow-3, kGreen, kBlue, kViolet+1, kMagenta};

  std::cout<<"Inside of plot_total_error"<<std::endl;
  
  //Now to plot
  /////////////
  canv = new TCanvas(Form("canv%s",name),Form("canv%s",name),2000,1500);
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

   //hist[5]->Draw("HIST TEXT00 SAME");
   
   h_total = total_error(Form("%s",name),hist);
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
   canv->Print(Form("images/Total_Error/_total_error%s.png",name));

}

void plotting_tools::make_total_covariance(const char* variable, const char* title,std::vector<TH2D*> hist_2D, const char* directory){

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
  TH2D* hist = (TH2D*)hist_2D[0]->Clone();
  for(int i=1; i < hist_2D.size(); i++){
    hist->Add(hist_2D[i]);
  }
  hist->Draw("colz text");
  hist->SetTitle(Form("%s",title)); //title
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
  hist->Write();
  gStyle->SetOptStat(0);
  canv->Print(Form("images/%s/Covariance_Matrices/_2D%s.png",directory,variable));
 
}
