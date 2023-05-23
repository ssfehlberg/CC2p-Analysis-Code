#include "constants.h"
using namespace Constants;


class Vertex{

 public:
  virtual void Plot_Vertex(const char* run,  Color_t colors[],Color_t colors_raquel[],string path);

 private:
  virtual void Fix_Statistical(TH1D* h_ext,TH1D* h_dirt,TH1D* h_overlay,TH1D* h);
  virtual void plotting(const char* run, char const* pot_num, const char* sample_name,Color_t colors[], std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2, TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1, TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend, int ymax, int ymin, int num_channels, const char* titles, string path,  TLine* a2,const char* titles2 = "", const char* plots="", const char* cut = "", bool plot_total = false, bool plot_ccnue = false, bool flip_legend = false,double pad_lim = 0.0, double pad0_lim = 0.19,double bnb_min = 0.0, double bnb_max = 2.7, Double_t xlim = -9999.0, Double_t xlim_up = +9999.0);

};//

void Vertex::Fix_Statistical(TH1D* h_ext,TH1D* h_dirt,TH1D* h_overlay,TH1D* h){

  double nbins = h_ext->GetXaxis()->GetNbins();
  
  for(int i = 1; i < nbins+1; i++){
    
    double MC = h_overlay->GetBinContent(i);
    double EXT = h_ext->GetBinContent(i);
    double Dirt = h_dirt->GetBinContent(i);
    double value = std::sqrt(MC + EXT + Dirt);
    h->SetBinError(i,value);
  }
}


void Vertex::plotting(const char* run, char const* pot_num, const char* sample_name, Color_t colors[], std::vector<TH1D*> h_overlay, TH1D* h_overlay00,  TH1D* h_overlay1,  TH1D* h_overlay2, TH1D* h_ext,  TH1D* h_ext1,  TH1D* h_ext2,  TH1D* h_dirt,  TH1D* h_dirt1,  TH1D* h_dirt2,  TH1D* h_bnb,  TH1D* h_bnb1, TCanvas* canv, THStack* h, TPad* pad, TPad* pad0, TLegend* legend, std::vector<const char*> channel_legend, int ymax, int ymin, int num_channels, const char* titles, string path,  TLine* a2,const char* titles2 = "", const char* plots="", const char* cut = "", bool plot_total = false, bool plot_ccnue = false, bool flip_legend = false,double pad_lim = 0.0, double pad0_lim = 0.19,double bnb_min = 0.0, double bnb_max = 2.7, Double_t xlim = -9999.0, Double_t xlim_up = +9999.0){

  

  h_ext1 = (TH1D*)h_ext->Clone();
  h_ext2 = (TH1D*)h_ext->Clone();
  h_dirt1 = (TH1D*)h_dirt->Clone();
  h_dirt2 = (TH1D*)h_dirt->Clone();
  h_bnb1 = (TH1D*)h_bnb->Clone();
  h_overlay00 = (TH1D*)h_overlay[0]->Clone();
  h_overlay1 = (TH1D*)h_overlay[0]->Clone();
  h_overlay2 = (TH1D*)h_overlay[0]->Clone();
  
  canv = new TCanvas(Form("C%s%s",plots,cut),Form("C%s%s",plots,cut),2000,1500);
  canv->SetGridx();
  h = new THStack(Form("h%s%s",plots,cut),Form("h%s%s",plots,cut));
   
  //Stacked Histrogram parameters
  h->Draw("HIST");
  h->SetTitle("");
  h->SetMaximum(ymax); 
  h->SetMinimum(ymin);
  
  //MC
  int z;
  int f;
  if(plot_ccnue){
    z = num_channels;
    f = num_channels;
  }else{
    z = num_channels-1;
    f = num_channels-1;
  }
          
  for(int k=1; k < z ; k++){
    if(k%2 != 0){
      h_overlay[k]->SetLineColor(colors[k]);
      h_overlay[k]->SetFillColor(colors[k]);
      h_overlay[k]->SetLineWidth(1);
      h->Add(h_overlay[k]);
      }
  }
  
   for(int k=1; k < z ; k++){
    if(k%2 == 0){
      h_overlay[k]->SetLineColor(colors[k]);
      h_overlay[k]->SetFillColor(colors[k]);
      h_overlay[k]->SetLineWidth(1);
      h->Add(h_overlay[k]);
    }
  }
  
  //Dirt
  h->Add(h_dirt);
  h_dirt->SetFillColor(kOrange-8);
  h_dirt->SetLineColor(kOrange-8);
  h_dirt->SetLineWidth(1);
  
  //EXT
  h->Add(h_ext);
  h_ext->SetFillColor(kViolet-7);
  h_ext->SetFillStyle(3005);
  h_ext->SetLineWidth(1);
  
  //BNB
  h_bnb->Draw("1e1p SAME");
  h_bnb->SetLineColor(kBlack);
  h_bnb->SetLineWidth(4);
  h_bnb->SetMarkerSize(1);

  h->GetXaxis()->SetTitle(Form("%s",titles));
  h->GetXaxis()->SetTitleSize(50); //35
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleOffset(1.7);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(50);

  h->GetYaxis()->SetTitle("No. Events");
  h->GetYaxis()->SetTitleSize(50);
  h->GetYaxis()->SetTitleFont(43); //4 = hevelatica normal 3 = precision
  h->GetYaxis()->SetTitleOffset(1.7);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(50);
  
  //Drawing cut lines if needed

  Double_t YMAX = gPad->GetFrame()->GetY2();
  
  a2 = new TLine(xlim,ymin,xlim,YMAX);  
  a2->Draw("SAME");
  a2->SetLineColor(kBlack);
  a2->SetLineWidth(3);
  
  a3 = new TLine(xlim_up,ymin,xlim_up,YMAX);  
  a3->Draw("SAME");
  a3->SetLineColor(kBlack);
  a3->SetLineWidth(3);
  
  //if you want to plot the total for sanity sake:
  if(plot_total){
    h_overlay[0]->Draw("SAME");
  }

  std::cout<<"Right before statistical"<<std::endl;

  h_overlay1->Add(h_ext1);
  h_overlay1->Add(h_dirt1);
  
  //Make sure to plot the statistical uncertainty
  Fix_Statistical(h_ext1,h_dirt1,h_overlay2,h_overlay1);
  h_overlay1->Draw("e2SAME");
  h_overlay1->SetLineColor(kBlack);
  h_overlay1->SetFillColor(kBlack);
  h_overlay1->SetFillStyle(3004);
  h_overlay1->SetMarkerSize(0);
  h_overlay1->SetLineWidth(1);

  if(flip_legend){
    legend = new TLegend(0.065, 0.52, 0.584, 0.87);
  }else{
    legend = new TLegend(0.37, 0.52, 0.889, 0.87);
  }

  legend->SetNColumns(2);
  legend->SetHeader("MicroBooNE 6.79 x 10^{20} POT, Preliminary","C"); // option "C" allows to center the header
  legend->AddEntry(h_bnb,"BNB Data (Stat.)","lepf");
  legend->AddEntry(h_ext,"EXT Data","f");
  legend->AddEntry(h_overlay1,"MC Stat. Unc.","f");
  legend->AddEntry(h_dirt,"Dirt","f");
  for(int k =1; k < z; k++){	
    legend->AddEntry(h_overlay[k],Form("%s",channel_legend[k]),"f");
  }
  legend->SetLineWidth(2);
  legend->SetLineColor(kGray);
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.03);
  legend->Draw("same");
  t->SetNDC();
  t->SetTextAlign(22);
  t->DrawLatex(0.515,0.95,Form("#scale[1.0]{%s %s}",titles,titles2));

  canv->Print(Form("%s%s%s.png",path.c_str(),plots,cut));
  canv->Print(Form("%s%s%s.pdf",path.c_str(),plots,cut));
   
}

void Vertex::Plot_Vertex(const char* run,  Color_t colors[],Color_t colors_raquel[],string path){

  for(int i = 0; i< num_cuts; i++){
    for(int j = 0; j< num_variables; j++){
      for(int k=0; k < num_channels; k++){
	h_overlay_vec.push_back(h_overlay[i][j][k]);
      }
      for(int k=0; k < num_channels_raquel; k++){
	h_overlay_raquel_vec.push_back(h_overlay_raquel[i][j][k]);
      }

      plotting(run,pot_num,sample_name,colors, h_overlay_vec, h_overlay0[i][j][0],h_overlay0[i][j][1],h_overlay0[i][j][2],
		       h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2],
		       h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],
		       h_bnb[i][j][0],h_bnb[i][j][1],canv[i][j], h[i][j],
		       pad[i][j], pad0[i][j], legend[i][j], channel_legend,
		       ylim[run_num][i][j],ymin[run_num][i][j], num_channels, titles[j], path, a1 ,"", plots[j], cut[i], false, true, false, 0.0, 0.19, 0.5, 1.5,  xmin_vtx[j], xmax_vtx[j]);
      h_overlay_vec.clear();
      
      plotting(run,pot_num,sample_name,colors_raquel, h_overlay_raquel_vec, h_overlay0_raquel[i][j][0],h_overlay0_raquel[i][j][1],h_overlay0_raquel[i][j][2],
	       h_ext[i][j][0], h_ext[i][j][1], h_ext[i][j][2],
	       h_dirt[i][j][0],h_dirt[i][j][1],h_dirt[i][j][2],
	       h_bnb[i][j][0],h_bnb[i][j][1],canv_raquel[i][j], h_raquel[i][j],
	       pad_raquel[i][j], pad0_raquel[i][j], legend_raquel[i][j], channel_legend_raquel,
	       ylim[run_num][i][j],ymin[run_num][i][j], num_channels_raquel, titles[j], path,   a1,"", plots[j], Form("%s_raquel",cut[i]), false, true, false, 0.0, 0.19, 0.5, 1.5, xmin_vtx[j], xmax_vtx[j]);
      h_overlay_raquel_vec.clear();
     
    }
  }

  
}//
