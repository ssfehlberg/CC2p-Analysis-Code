#ifndef xsec_prep_h
#define xsec_prep_h

class xsec_prep{

public:
  virtual void main();
  virtual void Grab_Histograms(TFile* f_bnb,TFile* f_ext,TFile* f_overlay,TFile* f_dirt, TFile* f_eff);
  virtual void Define_Constants();
  virtual void Determine_Bin_Size(TH2D* hist);
  
private:


  //Debug Statements
  //////////////////
  bool _debug = false;

  //Stuff for date and time
  /////////////////////////
  time_t now = time(0);
  tm *ltm = localtime(&now);
  int Day = ltm->tm_mday;
  int Month = ltm->tm_mon + 1;
  int Year = ltm->tm_year + 1900;

  //Migration Matrices
  //////////////////////////////////////////////
  static const int num_particles_matrices_plots = 3;
  const char* particles_matrices_plots[num_particles_matrices_plots] ={"_mom","_costheta","_phi"};
  static const int num_particles_matrices = 5;
  const char* particles_matrices[num_particles_matrices] = {"_muon_all","_muon_contained","_muon_uncontainied","_lead_proton","_recoil_proton"};
  TH2D* h_particle_matrices[num_particles_matrices][num_particles_matrices_plots];
  TH2D* h_particle_matrices_normalized[num_particles_matrices][num_particles_matrices_plots];
  TCanvas* canv_particle_matrices[num_particles_matrices][num_particles_matrices_plots];
  TCanvas* canv1_particle_matrices[num_particles_matrices][num_particles_matrices_plots];
  const char* titles_particle_matrices[num_particles_matrices_plots] = {"Momentum (GeV/c)","cos(#theta)","#phi (Rad.)"};
  const char* titles_particles_matrices[num_particles_matrices] = {"All Muons","Contained Muons","Uncontained Muons","Leading Proton","Recoil Proton"};

  static const int num_other_matrices =9;
  const char* other_matrices[num_other_matrices] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading","_opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  TH2D* h_other_matrices[num_other_matrices];
  TH2D* h_other_matrices_normalized[num_other_matrices];
  TCanvas* canv_other_matrices[num_other_matrices];
  TCanvas* canv1_other_matrices[num_other_matrices];
  const char* titles_other_matrices[num_other_matrices] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})","cos(Opening Angle Muon and Leading)","cos(Opening Angle Muon and Both Protons)","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)", "#delta #phi_{T} (Deg.)","p_{n} (GeV/c)","Neutrino Energy"};

  //Efficiency Plots. Needed to make the Smearing Matrices
  ///////////////////////////////////////////////////////

  //Momentum, theta and phi of particles
  static const int num_particles_eff_plots = 3;
  const char* particles_eff_var[num_particles_eff_plots] = {"_mom","_costheta","_phi"};
  const char* particles_eff_var_titles[num_particles_eff_plots] = {"True Momentum (GeV/c)","True cos(#theta)","True #phi (Rad.)"};
  static const int num_particles_eff = 5;
  const char* particles_eff[num_particles_eff] = {"_muon_all","_muon_contained","_muon_uncontainied","_lead_proton","_recoil_proton"};
  const char* particles_eff_titles[num_particles_eff] = {"All Muons","Contained Muons","Uncontained Muons","Leading Proton","Recoil Proton"};
  TH1D* h_particle_num[num_particles_eff][num_particles_eff_plots];
  TH1D* h_particle_denom[num_particles_eff][num_particles_eff_plots];     
  TCanvas* canv_particle_eff[num_particles_eff][num_particles_eff_plots];
  //double xlim_eff[num_eff] = {0.1,0.1,0.1,0.25,0.25};
  //double xlim_high_eff[num_eff] = {-9999.,-9999.,-9999.,1.2,1.2};
  //TLine* a[num_eff];
  //TLine* a1[num_eff];
    
  //efficiency plots of other variables to determine potential xsec candidates
  static const int num_other_eff = 9;
  const char* other_eff[num_other_eff] = {"_opening_angle_protons_lab","_opening_angle_protons_com","_opening_angle_mu_leading","_opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_pn","_nu_E"};
  const char* other_eff_titles[num_other_eff] = {"cos(#gamma_{Lab})","cos(#gamma_{1#mu2p COM})","cos(Opening Angle Muon and Leading_{Lab})","cos(Opening Angle Muon and Both Protons_{Lab})","#delta P_{T} (GeV/c)","#delta #alpha_{T} (Deg.)","#delta #phi_{T} (Deg.)","True p_{n} (GeV/c)","True Neutrino Energy (GeV)"};
  TH1D* h_other_eff_num[num_other_eff];
  TH1D* h_other_eff_denom[num_other_eff];
  TCanvas* canv_other_eff[num_other_eff];

  //Smearing Matrices
  /////////////////////
  TH2D* h_particle_smearing[num_particles_matrices][num_particles_matrices_plots];
  TCanvas* canv_particle_smearing[num_particles_matrices][num_particles_matrices_plots];
  TH2D* h_other_smearing[num_other_matrices];
  TCanvas* canv_other_smearing[num_other_matrices];

  //Things for the XSec
  //////////////////////////
  static const int num_particles = 5;
  static const int num_particle_plots = 3;
  const char* particles[num_particles] = {"_muon_all","_muon_contained","_muon_uncontained","leading_proton","_recoil_proton"};
  const char* particle_plots[num_particle_plots] = {"_mom","_cos_theta","phi"};
  TH1D* h_particle[num_particles][num_particle_plots];
  TH1D* h_particle_bnb[num_particles][num_particle_plots];
  TH1D* h_particle_ext[num_particles][num_particle_plots];
  TH1D* h_particle_overlay[num_particles][num_particle_plots];
  TH1D* h_particle_dirt[num_particles][num_particle_plots];

  static const int num_plots = 9;
  const char* plots[num_plots] = {"_opening_angle_lab","opening_angle_com","opening_angle_mu_leading","opening_angle_mu_both","_delta_PT","_delta_alphaT","_delta_phiT","_pn","_neutrino_energy"};
  TH1D* h[num_plots];
  TH1D* h_bnb[num_plots];
  TH1D* h_ext[num_plots];
  TH1D* h_overlay[num_plots];
  TH1D* h_dirt[num_plots];

  //List to Write Histograms to
  ////////////////////////////
  std::vector<TH2D*> h_2D_list;
  
};
#endif
#ifdef xsec_prep_cxx

void xsec_prep::Grab_Histograms(TFile* f_bnb,TFile* f_ext,TFile* f_overlay,TFile* f_dirt, TFile* f_eff){

  //Migration Matrices
  /////////////////////
  for(int i=0; i < num_particles_matrices; i++){
    for(int j=0; j < num_particles_matrices_plots; j++){
      h_particle_matrices[i][j] = (TH2D*)f_eff->Get(Form("h_particle_matrices%s%s",particles_matrices[i],particles_matrices_plots[j]));
    }
  }
  
  for(int i=0; i < num_other_matrices; i++){
    h_other_matrices[i] = (TH2D*)f_eff->Get(Form("h_other_matrices%s",other_matrices[i]));
  }

  //Efficiency Plots
  ///////////////////////////
  for(int i=0; i < num_particles_eff; i++){
    for(int j=0; j < num_particles_eff_plots; j++){
      h_particle_num[i][j] = (TH1D*)f_eff->Get(Form("h_particle_num%s%s",particles_eff[i],particles_eff_var[j]));
      h_particle_denom[i][j] = (TH1D*)f_eff->Get(Form("h_particle_denom%s%s",particles_eff[i],particles_eff_var[j]));
    }
  }

  for(int i = 0; i < num_other_eff; i++){
    h_other_eff_num[i] = (TH1D*)f_eff->Get(Form("h_other_eff_num%s",other_eff[i]));
    h_other_eff_denom[i] = (TH1D*)f_eff->Get(Form("h_other_eff_denom%s",other_eff[i]));
  }

} //end of grab histograms

void xsec_prep::Define_Constants(){

  //Calculate the Constants for the XSec
  ///////////////////////////////////////////
  double rho_Ar = 1.3836; //density of argon in g*cm^-3
  double V = (256.35 - 20.0) * (233.0 - 20.0) * (1036.8 - 20.0);//volume of detector cm: x*y*z
  double N_A = 6.023e+23; //avogadro's number
  double m_mol = 39.95; //mass of argon in g*mol^-1
  double N_targets = (rho_Ar*V*N_A)/m_mol; //number of target nuclei
  std::cout<<"Value of Volume: "<<V<<std::endl;
  std::cout<<"Value of N_targets assuming 10cm FV border: "<<N_targets<<std::endl;
  
  double pot = 6.79e+20;
  double flux = pot*(1.16e+11/1.6e+20); //fraction is MCC8 CCInclusive flux/POT
  std::cout<<"Value of Flux: "<<flux<<""<<std::endl;


} //end of define constants


void xsec_prep::Determine_Bin_Size(TH2D* hist){

  int num_bins =  hist->GetXaxis()->GetNbins();
  double first_bin = hist->GetXaxis()->GetBinCenter(1);
  double last_bin =  hist->GetXaxis()->GetBinCenter(num_bins);
  double bin_width = 0.1;//0.1;
  
  if(_debug) std::cout<<"Value of num_bins: "<<num_bins<<std::endl;
  if(_debug) std::cout<<"Value of first Bin: "<<first_bin<<std::endl;
  if(_debug) std::cout<<"Value of last Bin: "<<last_bin<<std::endl;
  
  TH1D* p[num_bins][num_bins];
  TCanvas* canv_bin[num_bins][num_bins];
  TF1* g[num_bins][num_bins];

  for(int i=1; i < num_bins + 2; i++){
    int x_low = i;
    std::vector<std::pair<double,char*>> values;

    for(int j=1; j < num_bins + 2; j++){
      int x_high = j;

      if(x_low < x_high){
	 p[i][j] = hist->ProjectionY(Form("p_%d-%d",x_low,x_high),x_low,x_high);
	 int num_entries = p[i][j]->GetEntries();

	 if(num_entries != 0){
	   canv_bin[i][j] = new TCanvas(Form("canv_bin_%d-%d",x_low,x_high),Form("canv_bin_%d-%d",x_low,x_high),2000,1500);
	   p[i][j]->Draw("hist");
	   p[i][j]->Fit("gaus","Q");
	   g[i][j] = (TF1*)p[i][j]->GetListOfFunctions()->FindObject("gaus");
	   g[i][j]->Draw("SAME");

	   double sigma = g[i][j]->GetParameter(2);
	   double bin_width = 0.05;
	   double xlow  = double(p[i][j]->GetXaxis()->GetBinCenter(x_low))  - bin_width;
	   double xhigh = double(p[i][j]->GetXaxis()->GetBinCenter(x_high)) - bin_width;
	   double bin_size = double(xhigh - xlow);
	   //double value = std::abs((2.0*sigma) - bin_size)/std::abs(bin_size);
	   double value = std::pow((2.0*sigma) -  bin_size,2)/ bin_size;
	   if(value < 1.0){
	   //std::cout<<"Value of Value: "<<value<<" Bin Range: "<<xlow<<"-"<<xhigh<<std::endl;
	   values.push_back(std::make_pair(value,Form("%f-%f",xlow,xhigh)));
	 
	   }
	   std::sort(values.begin(),values.end(),greater());
	     
	   if(_debug) std::cout<<"Value of X_low: "<<x_low<<std::endl;
	   if(_debug) std::cout<<"Value of X_high: "<<x_high<<std::endl;
	   if(_debug) std::cout<<"Number of Entries: "<<num_entries<<std::endl;
	   if(_debug) std::cout<<"Value of Xlow: "<<xlow<<std::endl;
	   if(_debug) std::cout<<"Value of XHigh: "<<xhigh<<std::endl;
	   if(_debug) std::cout<<"Value of Sigma: "<<sigma<<std::endl;
	   if(_debug) std::cout<<"Value of Bin Size: "<<bin_size<<std::endl;
	   if(_debug) std::cout<<"Value of Value: "<<value<<std::endl;

	 }// if entries	 
      } // if xlow
    } // j loop

    if(_debug) std::cout<<"Size of values: "<<values.size()<<std::endl;
    if(values.size() != 0){
      std::cout<<"[LOOKIE HERE] Minimum Difference: "<<values[0].first<<" Bin Range: "<<values[0].second<<std::endl;
      //std::cout<<"[LOOKIE HERE AGAIN] Second Minimum Difference: "<<values[1].first<<" Bin Range: "<<values[1].second<<std::endl;
      //std::cout<<"[LOOKIE HERE LAST TIME] Third Minimum Difference: "<<values[2].first<<" Bin Range: "<<values[2].second<<std::endl;
    } //values size
          
  } // i loop

} //Determine bin size
 
#endif
