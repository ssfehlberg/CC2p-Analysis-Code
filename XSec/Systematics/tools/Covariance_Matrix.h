TH2D* make_covariance_matrix(TH1D* h_xsec_CV,std::vector<TH1D*> h_xsec_universe, double N_univ, const char* name,bool DIRT = false){

  ///////////////////////////////
  //First let us calculate nu_CV
  //////////////////////////////
  int nbins = h_xsec_CV->GetNbinsX();
  double nu_CV[nbins]; //array of length nbins that will contain the central values   
   for(int reco_bin = 1; reco_bin <nbins+1; reco_bin++){
     nu_CV[reco_bin-1] =  h_xsec_CV->GetBinContent(reco_bin);
  }

  ///////////////////////////////////////////
  //Now we have to get the nu_Universe 
  ///////////////////////////////////////
  int num_universes = N_univ;
  double nu_universe[nbins][num_universes];
  for(int reco_bin = 1; reco_bin < nbins+1; reco_bin++){
    for(int u = 0; u < num_universes; u++){
      nu_universe[reco_bin-1][u] = h_xsec_universe[u]->GetBinContent(reco_bin);  
    }
  }
    
  ////////////////////////////////////
  //Now to Make our Covariance Matrix
  ///////////////////////////////////
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = h_xsec_CV->GetBinLowEdge(i+1);
  }

  TH2D* covariance_matrix = new TH2D(Form("h_2D_covariance%s",name),Form("h_2D_covariance%s",name),nbins,edges,nbins,edges);

  for(int a = 1; a < nbins+1; a++){
    for(int b=1; b < nbins+1; b++){
	double sum = 0;

	for(int u = 0; u < num_universes; u++){
	  double x = nu_CV[a-1] - nu_universe[a-1][u];
	  double y = nu_CV[b-1] - nu_universe[b-1][u];
	  double product = x * y;
	  sum += product;
	}
	double value = std::pow(num_universes,-1) * sum;
	covariance_matrix->SetBinContent(a,b,value);

    } //end of loop over b                                                                                                                                                                
  } //end of loop over a                                                                                                                                                                  

  return covariance_matrix;

} //end of make covariance matrix

TH1D* make_error_histogram(const char* name, TH1D* hist_CV,TH2D* h_covariance){

  int nbins = hist_CV->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist_CV->GetBinLowEdge(i+1);
  }

  TH1D* hist_errors = new TH1D(Form("hist_errors%s",name),Form("hist_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    for(int j=1; j < nbins+1; j++){
      if(i == j){
	double bin_content = h_covariance->GetBinContent(i,j);
	double value = std::pow(bin_content,0.5);
	hist_errors->SetBinContent(i,value);
      }
    }
  }

  return hist_errors;
  
} //end of make_error_histograms

TH1D* make_fractional_error_histogram(const char* name, TH1D* hist_CV,TH1D* hist_error){

  int nbins = hist_CV->GetNbinsX();
  Double_t edges[nbins+1];
  for(int i=0; i < nbins+1; i++){
    edges[i] = hist_CV->GetBinLowEdge(i+1);
  }

  TH1D* hist_fractional_errors = new TH1D(Form("hist_fractional_errors%s",name),Form("hist_fractional_errors%s",name),nbins,edges);

  for(int i=1; i < nbins+1; i++){
    double bin_content_error = hist_error->GetBinContent(i);
    double bin_content_CV = hist_CV->GetBinContent(i);
    double value = bin_content_error / bin_content_CV;
    hist_fractional_errors->SetBinContent(i,value);
  }
  
  return hist_fractional_errors;
}
