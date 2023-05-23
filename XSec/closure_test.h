#include "constants.h"
using namespace Constants;

class closure{

 public:
  virtual void Check_Efficiency(TCanvas* c_eff,TH1D* h_eff, TH1D* h_num, TH1D* h_denom, const char* title, const char* name);
  virtual void Check_Smearing_Matrix(TCanvas* c_2d, TH2D* h_matrix, const char* title,const char* name);
  virtual void Closure_Test(TCanvas* c_eff,TH1D* h_eff,TH1D* h_num, TH1D* h_denom, const char* title, const char* name,
			    TCanvas* c_2d, TH2D* h_matrix,
			    TH1D* h_cc2p,double maximum, bool print_contents = true);
  

};

///////////////////////////////
//Check to make sure the grabbed efficiency is the same as the calculated efficiency
///////////////////////////////
void closure::Check_Efficiency(TCanvas* c_eff,TH1D* h_eff,TH1D* h_num, TH1D* h_denom, const char* title, const char* name){
  
  c_eff = new TCanvas("c_eff","Efficiency",800,600);
  c_eff->SetRightMargin(0.09);
  c_eff->SetLeftMargin(0.15);
  c_eff->SetBottomMargin(0.15);

  //eff from xsec prep
  h_eff->Draw("HIST") ;
  h_eff->SetFillColor(kBlue);
  h_eff->SetTitle(title);
  h_eff->SetMaximum(0.4);

  //sanity check of xsec prep eff
  TH1D* h_num_eff = (TH1D*)h_num->Clone();
  TH1D* h_denom_eff = (TH1D*)h_denom->Clone();
  h_num_eff->Divide(h_num_eff,h_denom_eff,1,1,"cp");
  h_num_eff->Draw("1p SAME");
  h_num_eff->SetLineColor(kMagenta);
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);
  c_eff->Print(Form("images/Closure_Test/efficiency_check%s.png",name));
 

}

//Make Sure the Smearing Matrix is Okay
//////////////////////////////////////
void closure::Check_Smearing_Matrix(TCanvas* c_2d, TH2D* h_matrix, const char* title,const char* name){

  c_2d =  new TCanvas("c_2d","Smearing Matrix",800,600);
  c_2d->SetRightMargin(0.09);
  c_2d->SetLeftMargin(0.15);
  c_2d->SetBottomMargin(0.15);
  TH2D* h_smearing = (TH2D*)h_matrix->Clone();
  h_smearing->Draw("colz text");
  h_smearing->SetTitle(title);
  h_smearing->GetXaxis()->SetTitle(Form("True %s",title));
  h_smearing->GetYaxis()->SetTitle(Form("Reco. %s",title));
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);
  c_2d->Print(Form("images/Closure_Test/smearing_matrix%s.png",name));

}


void closure::Closure_Test(TCanvas* c_eff,TH1D* h_eff,TH1D* h_num, TH1D* h_denom, const char* title, const char* name,
				TCanvas* c_2d, TH2D* h_matrix,
				TH1D* h_cc2p, double maximum, bool print_contents = true){

  Check_Efficiency(c_eff,h_eff,h_num,h_denom, title,name); //checks to make sure the efficiency plot is calculated properly
  Check_Smearing_Matrix(c_2d,h_matrix,title,name); //checks to make sure we are using the correct smearing matrix


  //Now to actually perform the closure test and plot it.
  ///////////////////////////////////////////////////////
  RooUnfoldResponse ClosureResponseMatrix(0,0,h_matrix,"ClosureResponseMatrix");
  ClosureResponseMatrix.UseOverflow(kTRUE); //specify to include contributions from overflow bins. 

  RooUnfoldBayes RooUnfoldBayes_MC (&ClosureResponseMatrix, h_cc2p, 1);
  TH1D* h_unfolded_signal_closure = (TH1D*) RooUnfoldBayes_MC.Hreco();

  //RooUnfoldResponse ResponseMatrix(0,0,h2, "ResponseMatrix");
  //RooUnfoldBayes unfold_data (&ResponseMatrix, h_sig, iter);
  //TH1F* h_uf_sig= (TH1F*) unfold_data.Hreco();
 
  //Now to plot everything:
  /////////////////////////
  TCanvas* c_closure = new TCanvas("c_closure","Closure Test",800,600);
  c_closure->SetRightMargin(0.09);
  c_closure->SetLeftMargin(0.15);
  c_closure->SetBottomMargin(0.15);
  TH1D* h_closure = (TH1D*)h_unfolded_signal_closure->Clone();
  h_closure->Draw("hist");
  
  //divide by efficiency
  h_closure->Divide(h_unfolded_signal_closure,h_eff,1,1,"cp");
 
  //Drawing Effieicny on Same Plot
  h_denom->Draw("SAME 1e1p"); 
  h_denom->SetLineColor(kGreen);
  h_num->Draw("SAME 1e1p");
  h_num->SetLineColor(kRed);

  //Stupid setting parameters
  h_closure->SetMaximum(maximum);
  h_closure->SetTitle(title);
  h_closure->GetXaxis()->SetTitle(title);
  h_closure->GetYaxis()->SetTitle("Number of Events");
  TLegend* legend_closure = new TLegend(0.5, 0.71, 0.799, 0.89);
  legend_closure->AddEntry(h_closure,"1 Iteration Unfolded CC2p Signal","lp");
  legend_closure->AddEntry(h_denom,"All Generated CC2p Events","lp");
  legend_closure->AddEntry(h_num,"All True and Selected CC2p Events","lp");
  legend_closure->Draw("Same");
  gStyle->SetHistMinimumZero(kTRUE);
  gStyle->SetOptStat(0);
  c_closure->Print(Form("images/Closure_Test/closure_test%s.png",name));


  //Checking to make sure the closure test actually worked
  ////////////////////////////////////////////////////////

  if(print_contents){
    std::cout<<"---Bin Contents of Lower Bin For Muon All---"<<std::endl;
    std::cout<<"Signal: "<<h_cc2p->GetBinContent(0)<<std::endl;
    std::cout<<"Numerator: "<<h_num->GetBinContent(0)<<std::endl;
    std::cout<<"Denominator: "<<h_denom->GetBinContent(0)<<std::endl;
    std::cout<<"Smearing: "<<h_matrix->GetBinContent(0)<<std::endl;

    double num_bins = h_cc2p->GetNbinsX() + 1;
    std::cout<<"---Bin Contents of Upper Bin For Muon All---"<<std::endl;
    std::cout<<"Signal: "<<h_cc2p->GetBinContent(num_bins)<<std::endl;
    std::cout<<"Numerator: "<<h_num->GetBinContent(num_bins)<<std::endl;
    std::cout<<"Denominator: "<<h_denom->GetBinContent(num_bins)<<std::endl;
    std::cout<<"Smearing: "<<h_matrix->GetBinContent(num_bins)<<std::endl;


    std::cout<<"Number of Bins: "<<num_bins<<std::endl;
    double sum_particle = 0;
    double sum_num = 0;
    for(int i=0; i < num_bins+1; i++){
      double bin_content_particle = h_cc2p->GetBinContent(i);
      sum_particle += bin_content_particle;
      std::cout<<Form("PARTICLE: Bin %d Content: %f",i,bin_content_particle)<<std::endl;
      
      double bin_content_num = h_num->GetBinContent(i);
      sum_num += bin_content_num;
      std::cout<<Form("NUM: Bin %d Content: %f",i,bin_content_particle)<<std::endl;
    }
    std::cout<<"Sum Particle: "<<sum_particle<<std::endl;
    std::cout<<"Sum Numerator: "<<sum_num<<std::endl;
    std::cout<<"Difference: "<<(sum_particle - sum_num)<<std::endl;
    std::cout<<"Num. Integral: "<<h_num->Integral()<<std::endl;
    std::cout<<"Num. # of Entries: "<<h_num->GetEntries()<<std::endl;
  }
  
} //end of Closure test
