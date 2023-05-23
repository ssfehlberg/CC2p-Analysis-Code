class Flux{
 public:
  virtual std::pair<double,double> calculate_flux();
};


////////////////////////////////////////////////////////////////////////
//We're going to calculate our Number of Targest and Neutrino Flux First
//////////////////////////////////////////////////////////////////////

std::pair<double,double> Flux::calculate_flux(){

  double rho_Ar = 1.3836; //density of argon in g*cm^-3                                                                      
  double V = (256.35 - 20.0) * (233.0 - 20.0) * (1036.8 - 20.0); //volume of detector cm: x*y*z                               
  double N_A = 6.023e+23; //avogadro's number                                                                                
  double m_mol = 39.95; //mass of argon in g*mol^-1                                                                                                                                 
  double pot_num = 6.79E+20;
  
  //Neutrino flux stuff
  TFile* flux_file = new TFile("../../root_files/pelee/Run_all/neutrino_flux.root","READ");
  TCanvas* canv_flux = new TCanvas("canv_flux","canv_flux",700,500);
  TH1D* h_flux = (TH1D*)flux_file->Get("hEnumu_cv");
  double scale_factor= 1/(4997.*5e8)/(256.35*233);
  h_flux->Scale(scale_factor);
  h_flux->Draw("HIST");
  
  h_flux->GetXaxis()->SetRangeUser(0.0,4.0);
  h_flux->GetXaxis()->SetTitle("E_{#nu} (GeV)");
  h_flux->GetYaxis()->SetTitle("# #nu/4.183e+7 POT/GeV/cm^{2}");
  h_flux->SetTitle("Incoming Neutrino Flux");
  TF1* g = new TF1("g","gaus");
  TFitResultPtr fit = h_flux->Fit(g,"Q");
  Double_t* parameters = g->GetParameters();
  double ymax = 33.5E-12;//h.GetMaximum()
  
  Double_t mean = parameters[1];
  TLine* mean_line = new TLine(mean,0,mean,ymax);
  mean_line->Draw("SAME");
  mean_line->SetLineColor(2);
  mean_line->SetLineStyle(1);
  
  Double_t sigma = parameters[2];
  double sigma_plus = mean+sigma;
  double sigma_minus = mean-sigma;
  
  TLine* a_plus = new TLine(sigma_plus,0,sigma_plus,ymax);
  a_plus->Draw("SAME");
  a_plus->SetLineColor(4);
  a_plus->SetLineStyle(9);
 
  TLine* a_minus = new TLine(sigma_minus,0,sigma_minus,ymax);
  a_minus->Draw("SAME");
  a_minus->SetLineColor(4);
  a_minus->SetLineStyle(9);
  
  TLegend* legend = new TLegend(0.57,0.57,0.87,0.87);
  legend->AddEntry(h_flux,"Simulated BNB #nu_{#mu} Flux","L");
  legend->AddEntry(mean_line,Form("<E_{#nu}>= %f GeV",mean),"L");
  legend->AddEntry(a_plus,Form("+1#sigma Energy Range= %f GeV",sigma_plus),"L");
  legend->AddEntry(a_minus,Form("-1#sigma Energy Range= %f GeV",sigma_minus),"L");
  legend->SetBorderSize(0);
  legend->Draw("SAME");
  
  TLatex* t;
  t->DrawLatex(0.8,0.93,"#scale[0.6]{MicroBooNE In-Progress}");
  
  canv_flux->Print("images/neutrino_flux.png");
  
  //Save the number of targets and the flux value
  double N_targets = (rho_Ar*V*N_A)/m_mol; //number of target nuclei 
  double flux_value = (h_flux->Integral()) * pot_num;

  std::cout<<Form("Density of 40-Argon: %f g*cm^3",rho_Ar)<<std::endl;
  std::cout<<Form("Molar Mass of 40-Argon: %f g*mol^-1",m_mol)<<std::endl;
  std::cout<<Form("Volume of Detector Used: %E cm^3",V)<<std::endl;
  std::cout<<Form("Estimate for Number of Target Nuclei w. 10 cm FV: %E",N_targets)<<std::endl;
  std::cout<<Form("Neutrino Flux Estimate: %E v cm^2",flux_value)<<std::endl;

  return std::make_pair(N_targets,flux_value);
  
}
