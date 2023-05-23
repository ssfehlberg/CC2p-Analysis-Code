#define proton_resolution_cxx

void proton_resolution(){

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
  TFile *f_eff=new TFile(Form("root_files/pelee/Run_all/histograms_mc_eff.root")); //All efficiency histograms
  static const int num_slices = 10;
  const char* slices[num_slices] = {"0.25_0.35","0.35_0.45","0.45_0.55","0.55_0.65","0.65_0.75","0.75_0.85","0.85_0.95","0.95_1.05","1.05_1.15","1.15_1.25"};
  

  //Lead Proton
  ///////////////
  TH1D* h_leading_resolution[num_slices-1];
  TCanvas* canv_leading[num_slices-1];
  TF1* lead[num_slices-1];
  TGraphErrors* lead_resolution = new TGraphErrors(num_slices-1);

  for(int i=0; i < num_slices-1; i++){
    h_leading_resolution[i] = (TH1D*)f_eff->Get(Form("h_leading_resolution_%s",slices[i]));
    canv_leading[i] = new TCanvas(Form("canv_leading_%s",slices[i]),Form("canv_leading_%s",slices[i]),2000,1500);
    //h_leading_resolution[i]->Draw("1e1p");
    h_leading_resolution[i]->Fit("gaus");
    h_leading_resolution[i]->GetXaxis()->SetTitle("Range Resolution (%)");
    h_leading_resolution[i]->GetYaxis()->SetTitle("# of Events");
    lead[i] = (TF1*)h_leading_resolution[i]->GetListOfFunctions()->FindObject("gaus");
    canv_leading[i]->Update();
    
    double mean_lead = lead[i]->GetParameter(1);
    double mean_lead_error = lead[i]->GetParError(1);
    lead_resolution->SetPoint(i,i+1,mean_lead);
    lead_resolution->SetPointError(i,0.5,mean_lead_error);
  }

  TCanvas* canv_lead_graph = new TCanvas("canv_lead_graph","canv_lead_graph",2000,1500);
  lead_resolution->Draw("ap");
  lead_resolution->SetLineColor(kBlack);
  lead_resolution->SetLineWidth(2);
  lead_resolution->GetXaxis()->SetTitle("Reco. Momentum Range (GeV/c)");
  lead_resolution->GetYaxis()->SetTitle("Mean Value (%)");
  lead_resolution->SetTitle("Lead Proton");
  //lead_resolution->SetMinimum(-0.6);
  //lead_resolution->SetMaximum(0.6);  
   const char* eff_cut_label[9] = {"0.35","0.45","0.55","0.65","0.75","0.85","0.95","1.05","1.15"};
  TAxis *xax = lead_resolution->GetXaxis();
  for(int i = 1; i < xax->GetXmax()-1; i++){
    int bin_index_up = xax->FindBin(i+1);
    int bin_index_low = xax->FindBin(i);
    xax->SetBinLabel(bin_index_low,Form("%s", eff_cut_label[i-1]));
    xax->SetBinLabel(bin_index_up,Form("%s", eff_cut_label[i]));
  }
  xax->Draw("hist");
  
  TLine* a = new TLine(0,0,11,0);
  a->Draw("same");
  a->SetLineColor(kRed);
  canv_lead_graph->Update();
  canv_lead_graph->Print("images/pelee/Run_all/lead_resolution.png");

  std::cout<<"-------Starting the Recoil--------"<<std::endl;
  
  //Recoil Proton
  ////////////////
  TH1D* h_recoil_resolution[num_slices-3];
  TCanvas* canv_recoil[num_slices-3];
  TF1* recoil[num_slices-3];
  TGraphErrors* recoil_resolution = new TGraphErrors(num_slices-3);
  
  for(int i=0; i < num_slices-3; i++){
    h_recoil_resolution[i] = (TH1D*)f_eff->Get(Form("h_recoil_resolution_%s",slices[i]));
    canv_recoil[i] = new TCanvas(Form("canv_recoil_%s",slices[i]),Form("canv_recoil_%s",slices[i]),2500,1500);
    h_recoil_resolution[i]->Draw("1e1p");
    h_recoil_resolution[i]->Fit("gaus");
    h_recoil_resolution[i]->GetXaxis()->SetTitle("Range Resolution (%)");
    h_recoil_resolution[i]->GetYaxis()->SetTitle("# of Events");
    recoil[i] = (TF1*)h_recoil_resolution[i]->GetListOfFunctions()->FindObject("gaus");
    canv_recoil[i]->Update();
    
    double mean_recoil = recoil[i]->GetParameter(1);
    double mean_recoil_error = recoil[i]->GetParError(1);
    recoil_resolution->SetPoint(i,i+1,mean_recoil);
    recoil_resolution->SetPointError(i,0.5,mean_recoil_error);

  }

  TCanvas* canv_recoil_graph = new TCanvas("canv_recoil_graph","canv_recoil_graph",2000,1500);
  recoil_resolution->Draw("ap");
  recoil_resolution->SetLineColor(kBlack);
  recoil_resolution->SetLineWidth(2);
  recoil_resolution->GetXaxis()->SetTitle("Reco. Momentum Range (GeV/c)");
  recoil_resolution->GetYaxis()->SetTitle("Mean Value (%)");
  recoil_resolution->SetTitle("Recoil Proton");
  recoil_resolution->SetMinimum(-20);
  recoil_resolution->SetMaximum(20);  
  const char* eff_cut_label0[7] = {"0.35","0.45","0.55","0.65","0.75","0.85","0.95"};
  TAxis *xax0 = recoil_resolution->GetXaxis();
  for(int i = 1; i < xax0->GetXmax()-1; i++){
    int bin_index_up = xax0->FindBin(i+1);
    int bin_index_low = xax0->FindBin(i);
    xax0->SetBinLabel(bin_index_low,Form("%s", eff_cut_label0[i-1]));
    xax0->SetBinLabel(bin_index_up,Form("%s", eff_cut_label0[i]));
  }
  xax0->Draw("hist");
  
  TLine* a0 = new TLine(0,0,8,0);
  a0->Draw("same");
  a0->SetLineColor(kRed);
  canv_recoil_graph->Update();
  canv_recoil_graph->Print("images/pelee/Run_all/recoil_resolution.png");
  
  
}
