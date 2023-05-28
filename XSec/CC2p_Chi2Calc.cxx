#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TString.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TMath.h>
#include <TLatex.h>
#include <TMatrixD.h>

#include <iostream>
#include <vector>
#include <sstream>
#include <string>

using namespace std;

// -------------------------------------------------------------------------------------------------------------------------------------

TString ToString(double num) {

	std::ostringstream start;
	start << num;
	string start1 = start.str();
	return start1;

}

//----------------------------------------//

TString to_string_with_precision(double a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return TString(out.str());
}

//----------------------------------------//

double FindOneDimHistoMaxValue(TH1D* h){

	int NBins = h->GetXaxis()->GetNbins();
	double HistoMax = -9999.;	

	for (int ibin = 1; ibin<= NBins; ibin++) {

		double LocalMax = h->GetBinContent(ibin);
		if (LocalMax > HistoMax) { HistoMax = LocalMax; }

	}

	return HistoMax;

}

//----------------------------------------//

void CalcChiSquared(TH1D* h_model, TH1D* h_data, TH2D* cov, double &chi, int &ndof, double &pval) {

	// Clone them so we can scale them 

	TH1D* h_model_clone = (TH1D*)h_model->Clone();
	TH1D* h_data_clone  = (TH1D*)h_data->Clone();
	TH2D* h_cov_clone   = (TH2D*)cov->Clone();
	int NBins = h_cov_clone->GetNbinsX();

	// Getting covariance matrix in TMatrix form

	TMatrixD cov_m;
	cov_m.Clear();
	cov_m.ResizeTo(NBins,NBins);

	// loop over rows

	for (int i = 0; i < NBins; i++) {			

		// loop over columns

		for (int j = 0; j < NBins; j++) {

			cov_m[i][j] = h_cov_clone->GetBinContent(i+1, j+1);
 
		}
	
	}

	TMatrixD copy_cov_m = cov_m;

	// Inverting the covariance matrix
	TMatrixD inverse_cov_m = cov_m.Invert();

	// Calculating the chi2 = Summation_ij{ (x_i - mu_j)*E_ij^(-1)*(x_j - mu_j)  }
	// x = data, mu = model, E^(-1) = inverted covariance matrix 

	chi = 0.;
	
	for (int i = 0; i < NBins; i++) {

		//double XWidth = h_data_clone->GetBinWidth(i+1);

		for (int j = 0; j < NBins; j++) {

			//double YWidth = h_data_clone->GetBinWidth(i+1);

			double diffi = h_data_clone->GetBinContent(i+1) - h_model_clone->GetBinContent(i+1);
			double diffj = h_data_clone->GetBinContent(j+1) - h_model_clone->GetBinContent(j+1);
			double LocalChi = diffi * inverse_cov_m[i][j] * diffj; 
			chi += LocalChi;

			/*std::cout<<"Value of diffi: "<<diffi<<std::endl;
			std::cout<<"Value of diffj: "<<diffj<<std::endl;
			std::cout<<"Value of inverse_cov_m: "<<inverse_cov_m[i][j]<<std::endl;
			std::cout<<"Value of LocalChi: "<<LocalChi<<std::endl;*/
		}

	}

	ndof = h_data_clone->GetNbinsX();
	pval = TMath::Prob(chi, ndof);

	delete h_model_clone;
	delete h_data_clone;
	delete h_cov_clone;

}

//----------------------------------------//

void CC2p_Chi2Calc() {

	//----------------------------------------//

	int DecimalAccuracy = 2;
	int FontStyle = 132;
	double TextSize = 0.06;	

	TH1D::SetDefaultSumw2();
	gStyle->SetEndErrorSize(4);		

	gStyle->SetPalette(55); 
	const Int_t NCont = 999; 
	gStyle->SetNumberContours(NCont); 
	gStyle->SetTitleSize(TextSize,"t"); 
	gStyle->SetTitleFont(FontStyle,"t");
	gStyle->SetOptStat(0);	
	gStyle->SetPaintTextFormat("4.2f");		

	//----------------------------------------//

	TFile* fXSec = TFile::Open("CC2p_xsec_results.root");

	//----------------------------------------//	

	vector<TString> PlotNames;

	PlotNames.push_back("opening_angle_protons_lab");
	PlotNames.push_back("opening_angle_mu_both");
	PlotNames.push_back("delta_PT");		

	const int NPlots = PlotNames.size();
	cout << "Number of 1D Plots = " << NPlots << endl;

	//----------------------------------------//

	// MC Samples to compare

	vector<TString> MCSampleBand;

	MCSampleBand.push_back("uboone_tune");	
	MCSampleBand.push_back("empirical");
	MCSampleBand.push_back("nieves");		
	MCSampleBand.push_back("susav2");
	MCSampleBand.push_back("nuisance");
//	MCSampleBand.push_back("nuwro");		

	int NMC = MCSampleBand.size();

	//----------------------------------------//
	
	// Loop over the plots

	for (int iplot = 0; iplot < NPlots; iplot ++) {

		TH1D* DataPlot = (TH1D*)( fXSec->Get("h_data_xsec_" + PlotNames[iplot]) );
		TH2D* Cov = (TH2D*)fXSec->Get("h_2D_total_covariance_matrix_"+PlotNames[iplot]);		
		//DataPlot->Draw();

		//----------------------------------------//						

		TH1D* MCPlot[NMC];
		double Chi2[NMC];
		int Ndof[NMC];
		double pval[NMC];

		// Loop over the MC predictions

		cout << endl;

		for (int igen = 0; igen < NMC; igen++) {			

			TString InteName = "h_" + MCSampleBand[igen] + "_xsec_" + PlotNames[iplot];
			MCPlot[igen] = (TH1D*)( fXSec->Get(InteName) );	
			//MCPlot[igen]->Draw("same hist");	

			CalcChiSquared(MCPlot[igen],DataPlot,Cov,Chi2[igen],Ndof[igen],pval[igen]);
			TString Chi2NdofAlt = " (" + to_string_with_precision(Chi2[igen],1) + "/" + TString(std::to_string(Ndof[igen])) +")";									

			cout << PlotNames[iplot] << " " << MCSampleBand[igen] + Chi2NdofAlt << endl;

		} // End of the loop over the generators

		cout << endl;

	} // End of the loop over the plots

} // End of the program 
