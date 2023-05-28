/////////////////////////////////////
//July 11th, 2022
//Author: Samantha Sword-Fehlberg
//Base plotting class
//includes functions to calculate covariance matrices, plot fractional error, and other stuff
//Functioins are all free to allow for multiple class access
///////////////////////////////////////////////

namespace shared{

  void plot_multisims(const char* variable, const char* title, std::vector<TH1D*> hist, TH1D* hist_CV, double ymin, double ymax, int num_universes, const char* directory);
  void make_total_covariance_matrix(const char* variable, const char* title, std::vector<TH2D*> hist_2D, const char* directory);
  void plot_covariance(const char* variable, const char* title, TH2D* hist, const char* directory);
  void plot_error(const char* variable, const char* title, TH1D* hist_error, const char* directory);
  void plot_fractional_error(const char* variable, const char* title, TH1D* hist_error, const char* directory);
  TH1D* total_error(const char* name,std::vector<TH1D*> hist);
  
}; //end of namespace definition

//////////
//[PLOT MULTISIMS]
//plots the multisim curves in terms of number of events
//variable: which variable we are considering
//title: title of plot
//hist: vector of TH1D. Each TH1D represents a single multisim universe prediction
//hist_CV: our central value MC prediction for CC2p events. Follows our predefined signal
// ymin/ymax: limmits of the y-axis
//num_universes: number of multisim universes
//directory: where to save the images
/////////////////////
void shared::plot_multisims(const char* variable, const char* title, std::vector<TH1D*> hist, TH1D* hist_CV, double ymin, double ymax, int num_universes, const char* directory){

  TCanvas* canv = new TCanvas(Form("canv%s",variable),Form("canv%s",variable),2000,1500);

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

  
  TLegend* legend = new TLegend(0.48, 0.65, 0.89, 0.89);
  legend->AddEntry(hist_CV,Form("Central Value"),"L");
  legend->AddEntry(hist[0],Form("Multisim Universe"),"L");
  legend->Draw("same");
  gStyle->SetOptStat(0);
  canv->Print(Form("../images/%s/Multisims/%s.png",directory,variable));
  
}

///////////////////////////
//[make_total_covariance_matrix]
//Adds vector of covariance matrices (hist_2D) to create covariance matrix
//Plot is then saved to specified directory
///////////////////////////
void shared::make_total_covariance_matrix(const char* variable, const char* title, std::vector<TH2D*> hist_2D, const char* directory){

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
  canv->Print(Form("../images/%s/Covariance_Matrices/_2D%s.png",directory,variable));
 
}

////////////////////////////////
//[PLOT COVARIANCE]
//Plots the covariance matrix for a specific sample
//variable: variable we are looking at
//title: title of the plot
//hist: a TH2D that represents the covariance matrix
//hist_CV: TH1D central value of our CC2p signal. Need this to get the bin edges
//directory: where to save the image
//////////////////////////////
void shared::plot_covariance(const char* variable, const char* title, TH2D* hist, const char* directory){

  TCanvas* canv =  new TCanvas(Form("canv_2D%s",variable),Form("canv_2D%s",variable),2000,1500);
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
  canv->Print(Form("../images/%s/Covariance_Matrices/_2D%s.png",directory,variable));
    
}

//////////
//[PLOT ERROR]
//Plots error in terms of number of events
//variable: variaable name
//title: title of plot and x axis
//hist_error: histograms we want to plot
//directory: where to save the plot images
//////////////////
void shared::plot_error(const char* variable, const char* title, TH1D* hist_error, const char* directory){

  TCanvas* canv =  new TCanvas(Form("canv_error%s",variable),Form("canv_error%s",variable),2000,1500);
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
  canv->Print(Form("../images/%s/Error_Number_of_Events/_error%s.png",directory,variable));
}

////////////////////////
//[PLOT FRACTIONAL ERROR]
//Plots the fractional error
//variable: which variable
//title: title of plot
//hist_error: Th1D we wish to plot. NOTE: ASSUMES HIST_ERROR HAS ALREADY BEEN AREA NORMALIZED
//directory: where to save the image
/////////////////////////
void shared::plot_fractional_error(const char* variable, const char* title, TH1D* hist_error, const char* directory){

 TCanvas* canv =  new TCanvas(Form("canv_fractional_error%s",variable),Form("canv_fractional_error%s",variable),2000,1500);
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
 canv->Print(Form("../images/%s/Fractional_Error/_fractional_error%s.png",directory,variable));

}

////////////////////////////////////////////////////////////////
//[TOTAL ERROR]
//give a variable name (name) and a vector of TH1D (hist), calcultes the total
//error, by adding the bin contributions in quadrature
//Produces a new Th1D (h_total_error) which is a histogram of the total error
/////////////////////////////////////////////////////////////////////////
TH1D* shared::total_error(const char* name,std::vector<TH1D*> hist){

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

//////////////////////////////////////
//[]
