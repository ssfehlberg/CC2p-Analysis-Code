#include "constants.h"
using namespace Constants;


class Efficiency{

 public:
  virtual void Plot_Efficiency(string path);

};//


void Efficiency::Plot_Efficiency(string path){


  gStyle->SetEndErrorSize(2);  

  
  /////////////////////////////
  //Plot the Efficiency Stuff:
  //////////////////////////
  for(int i=0; i < num_eff; i++){
    
    h_num0[i] = (TH1D*)h_num[i]->Clone();
    h_denom0[i] = (TH1D*)h_denom[i]->Clone();
    canv_both[i] = new TCanvas(Form("C_both_%s",eff[i]),Form("C_both_%s",eff[i]),2000,1500);
    h_denom0[i]->Draw("hist");
    h_denom0[i]->SetFillColor(kBlue);
    h_denom0[i]->SetTitle(Form("Efficiency(%s); %s ; Number of Entries",titles_eff[i],titles_eff[i]));
    h_denom0[i]->SetMaximum(ylim_eff[run_num][i]);
    h_num0[i]->Draw("histSAME");
    h_num0[i]->SetFillColor(kRed);    
    legend_eff[i] = new TLegend(0.63,0.58,1.0,0.89);
    legend_eff[i]->AddEntry(h_num0[i], "Numerator", "f");
    legend_eff[i]->AddEntry(h_denom0[i], "Denomiator", "f");
    legend_eff[i]->Draw("same");

    std::cout<<Form("Number of Entries in %s Denominator: ",titles_eff[i])<<h_denom0[i]->GetEntries()<<std::endl;
    std::cout<<Form("Number of Entries in %s  Numerator: ",titles_eff[i])<<h_num0[i]->GetEntries()<<std::endl;
    
    canv_both[i]->Print(Form("%s_%s_both.png",path.c_str(),eff[i]));
    canv_both[i]->Print(Form("%s_%s_both.pdf",path.c_str(),eff[i]));

    h_num1[i] = (TH1D*)h_num[i]->Clone();
    h_denom1[i] = (TH1D*)h_denom[i]->Clone();
    canv_eff[i] = new TCanvas(Form("C_eff_%s",eff[i]),Form("C_eff_%s",eff[i]),2000,1500);
    h_num1[i]->Divide(h_num1[i],h_denom1[i],1.0,1.0, "B");
    h_num1[i]->Draw("1e1p");
    h_num1[i]->SetTitle(Form(" ; %s ; Efficiency",titles_eff[i]));
    h_num1[i]->SetLineColor(kViolet);
    h_num1[i]->SetMaximum(0.4);
    h_num1[i]->SetMinimum(0);

    
    t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s}",titles_eff[i]));
    //t->DrawLatex(0.3,0.92,Form("%s",pot_num));
    t->DrawLatex(0.73,0.88,"#scale[0.5]{MicroBooNE 6.79 x 10^{20} POT, Preliminary}");
    a[i] = new TLine(xlim_eff[run_num][i],0,xlim_eff[run_num][i],0.4);
    a[i]->Draw("same");
    a[i]->SetLineColor(kBlack);
    a[i]->SetLineWidth(4);
    a_up[i] = new TLine(xlim_up_eff[run_num][i],0,xlim_up_eff[run_num][i],0.4);
    a_up[i]->Draw("same");
    a_up[i]->SetLineColor(kBlack);
    a_up[i]->SetLineWidth(4);
    canv_eff[i]->Print(Form("%s_%s_eff.png",path.c_str(),eff[i]));
    canv_eff[i]->Print(Form("%s_%s_eff.pdf",path.c_str(),eff[i]));
  }

  //Effieincy of other variables
  for(int i =0; i < num_particles_eff; i++){
    for (int j=0; j <num_particles_eff_plots; j++){
      canv_particle_eff[i][j] = new TCanvas(Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),Form("C_eff%s%s",particles_eff[j],particles_eff_var[i]),2000,1500);

      h_particle_num[i][j]->Divide(h_particle_num[i][j],h_particle_denom[i][j],1.0,1.0, "B");
      h_particle_num[i][j]->Draw("1e1p");
      h_particle_num[i][j]->SetTitle(Form(" ; %s ; Efficiency",particles_eff_var_titles[j]));
      h_particle_num[i][j]->SetLineColor(kViolet);
      h_particle_num[i][j]->SetMaximum(0.4);
      h_particle_num[i][j]->SetMinimum(0);
      t->DrawLatex(0.515,0.97,Form("#scale[1.0]{Efficiency: %s of %s}",particles_eff_var_titles[j],particles_eff_titles[i]));
      t->DrawLatex(0.23,0.92,Form("%s",pot_num));
      t->DrawLatex(0.8,0.92,"#scale[0.5]{MicroBooNE In-Progress}");
      //a[i] = new TLine(xlim_eff[run_num][i],0,xlim_eff[run_num][i],1);
      //a[i]->Draw("same");
      //a[i]->SetLineColor(kBlack);
      //a[i]->SetLineWidth(4);
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff.png",path.c_str(),particles_eff[i],particles_eff_var[j]));
      canv_particle_eff[i][j]->Print(Form("%s%s%s_eff.pdf",path.c_str(),particles_eff[i],particles_eff_var[j]));
 
    }
  }
} //end of plot efficiency
