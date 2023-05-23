#include <sys/types.h>
#include <sys/stat.h>

const char* which_sample(){

  //Get the correct Sample
  char response;
  const char* Sample;
  std::cout<<"Which Sample is This?"<<std::endl;
  std::cout<<" 0 = pelee \n 1 = pelee with xsec binning"<<std::endl;
  std::cin>>response;

  if(response == '0'){
    Sample = "pelee";
  }else if (response == '1'){
    Sample = "pelee_xsec";
  }else{
    std::cout<<"Invalid Response. Please Type 0, 1, or 2 for pelee, filtered, and unfiltered samples respectively."<<std::endl;
  }

  return Sample;

}

std::pair<const char*, const char*> which_run(){
  
  //Get the Correct Run
  char response1;
  const char* Run;
  const char* Run_Title;
  std::cout<<"Which Run is This?"<<std::endl;
  std::cout<<" 0 = Jan \n 1 = Run 1 \n 2 = Run 2 \n 3 = Run 3 \n 4 = Runs 1+2+3"<<std::endl;
  std::cin>>response1;

  if(response1 == '0'){
    Run = "Jan";
    Run_Title = "January Sample";
  } else if(response1 == '1'){
    Run = "Run1";
    Run_Title = "Run 1";
  } else if(response1 == '2') {
    Run = "Run2";
    Run_Title = "Run 2";
  } else if (response1 == '3'){
    Run = "Run3";
    Run_Title = "Run 3";
  }else if(response1 == '4'){
    Run = "Run_all";
    Run_Title = "Runs 1+2+3";
  }else{
      std::cout<<"Invalid Response. Please Type 0, 1, 2 ,3, 4, or 5 for January, Run 1, Run2 , Run 3, Runs 1+2+3, respectively."<<std::endl;
  }

  return std::make_pair(Run,Run_Title);

}

/******************************************************************************  
 * Checks to see if a directory exists. Note: This method only checks the                                                   
 * existence of the full path AND if path leaf is a dir.                                                                      
 *                                                                                                             
 * @return  >0 if dir exists AND is a dir,                                                                                    
 *           0 if dir does not exist OR exists but not a dir,                                                 
 *          <0 if an error occurred (errno is also set)                                                                       
*****************************************************************************/

int dirExists(const char* const path)
{
    struct stat info;

    int statRC = stat( path, &info );
    if( statRC != 0 )
    {
        if (errno == ENOENT)  { return 0; } // something along the path does not exist                                                                                                                                                     
        if (errno == ENOTDIR) { return 0; } // something in path prefix is not a dir                                                                                                                                                       
        return -1;
    }

    return ( info.st_mode & S_IFDIR ) ? 1 : 0;
}

double calculatePearsonChiSq(TH1D* O, TH1D* E){
    double chisq = 0;
    for (int i = 1; i < O->GetNbinsX()+1; i++){
        double O_i = O->GetBinContent(i);
        double E_i = E->GetBinContent(i);
        if (O_i == 0 && E_i == 0){
	  chisq += 0;

	} else if (E_i == 0){
	  chisq += 0;

	}else{
	  //chisq += std::pow(O_i - E_i,2)/((O_i+E_i)/2);
	  chisq += std::pow(O_i - E_i,2)/(E_i);

        } //end of else
    } //end of for loop

    return chisq;
}


void Make_Stat_Sys_Plot(TH1D* h,TH1D* h_MC_CV,TH1D* h_ext, TH1D* h_dirt,TH1D* h_stat,TH1D* h_systematic, bool up){

  double nbins = h_stat->GetXaxis()->GetNbins();

  for(int i = 1; i < nbins+1; i++){
    double MC_CV = h_MC_CV->GetBinContent(i);
    double EXT = h_ext->GetBinContent(i);
    double dirt = h_dirt->GetBinContent(i);
    double total = MC_CV+EXT+dirt;

    double statistical = h_stat->GetBinContent(i);
    double systematic = h_systematic->GetBinContent(i);
    double value_stat = total * statistical; //taking it from another code
    //double value_stat = std::sqrt(MC_CV+EXT+dirt); //calculating directly
    double value_sys = total * systematic;
    double value = std::sqrt(std::pow(value_stat,2) + std::pow(value_sys,2));
    
    std::cout<<"Value of MC_CV: "<<MC_CV<<std::endl;
    std::cout<<"Value of EXT: "<<EXT<<std::endl;
    std::cout<<"Value of dirt: "<<dirt<<std::endl;
    std::cout<<"Value of total: "<<total<<std::endl;
    std::cout<<"value of statistical: "<<statistical<<std::endl;
    std::cout<<"Value of systematic: "<<systematic<<std::endl;
    std::cout<<"Value of Value_Stat: "<<value_stat<<std::endl;
    std::cout<<"Value of Value_sys: "<<value_sys<<std::endl;
    std::cout<<"Value of Value: "<<value<<std::endl;

    //h->SetBinContent(i,total);
    if(up == true){
      h->SetBinContent(i,total+value);
    } else{
      h->SetBinContent(i,total-value);
    }
  }

}


void Fix_Statistical(TH1D* h,TH1D* h_MC_CV,TH1D* h_stat){

  double nbins = h_stat->GetXaxis()->GetNbins();
  for(int i = 1; i < nbins+1; i++){
    double MC_CV = h_MC_CV->GetBinContent(i);
    double statistical = h_stat->GetBinContent(i);
    std::cout<<"value of statistical: "<<statistical<<std::endl;
    double value = MC_CV * statistical;
    h->SetBinError(i,value);
  }
}

void Fix_Systematic(TH1D* h,TH1D* h_MC_CV,TH1D* h_systematic){

  double nbins = h->GetXaxis()->GetNbins();
  for(int i = 1; i < nbins+1; i++){
    double MC_CV = h_MC_CV->GetBinContent(i);
    double systematic = h_systematic->GetBinContent(i);
    double value = MC_CV * systematic;
    std::cout<<"Value of MC_CV: "<<MC_CV<<std::endl;
    std::cout<<"Value of systematic: "<<systematic<<std::endl;
    std::cout<<"Value of Value: "<<value<<std::endl;
    h->SetBinError(i,value);
  }
}
