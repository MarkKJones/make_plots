#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TEllipse.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_hist_beam(TString basename) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 //
     TString outputpdf;
 outputpdf="plots/"+basename+"_beam.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_beam_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   static const Int_t nh1=10;
   TString h1name[nh1]={"hxrast","hyrast","hxbeam_0","hxbeam_1","hxbeam_2","hxbeam_3","hybeam_0","hybeam_1","hybeam_2","hybeam_3"};
   TH1F* hist1d[nh1];
  for (Int_t ip=0;ip<nh1;ip++) {
    // cout <<  hname[ip] << endl;
       hist1d[ip] = (TH1F*)fhistroot->Get(h1name[ip]);
       if (!hist1d[ip]) cout << " no hist = " << h1name[ip] << endl;
  }
   static const Int_t nh2=3;
   TString h2name[nh2]={"hxrast_yrast","hxrast_yrast_tr","hxbeam_ybeam"};
   TH1F* hist2d[nh1];
  for (Int_t ip=0;ip<nh2;ip++) {
       hist2d[ip] = (TH1F*)fhistroot->Get(h2name[ip]);
       if (!hist2d[ip]) cout << " no hist = " << h2name[ip] << endl;
  }
  Double_t mean_x[4];
  Double_t mean_y[4];
     Double_t zbpm[4]={-320.82,-224.86,-129.44,0.};
  
  for (Int_t nb=0;nb<4;nb++) {
    mean_x[nb]=hist1d[2+nb]->GetMean();
    mean_y[nb]=hist1d[6+nb]->GetMean();
    cout << zbpm[nb] << " " << mean_x[nb] << " " << mean_y[nb] << endl;
  }
  //
  TCanvas *cbeam = new TCanvas("cbeam","Beam",800,700);
  cbeam->Divide(2,2);
  cbeam->cd(1);
  gPad->SetLogz();
  gPad->SetGridy();
  gPad->SetGridx();
  hist2d[0]->Draw("colz");
  TEllipse *yEllipse= new TEllipse(-mean_x[3],-mean_y[3],.1,.1);
  yEllipse->Draw("same");
  yEllipse->SetLineColor(2);
  yEllipse->SetFillStyle(0);
  yEllipse->SetLineWidth(2);
  TLine *yline= new TLine(-.3,-mean_y[3],.3,-mean_y[3]);
  yline->Draw();
  yline->SetLineColor(2);
  yline->SetLineWidth(2);
  TLine *xline= new TLine(-mean_x[3],-.3,-mean_x[3],.3);
  xline->Draw();
  xline->SetLineColor(2);
  xline->SetLineWidth(2);
  cbeam->cd(3);
  //gPad->SetLogz();
  gStyle->SetOptStat(1111);
  hist1d[5]->Draw();
  cbeam->cd(4);
  //gPad->SetLogz();
  gStyle->SetOptStat(1111);
  hist1d[9]->Draw();
  cbeam->Print(outputpdf+"(");
  //
  TGraph *grbeamx = new TGraph(4,zbpm,mean_x);
  TGraph *grbeamy = new TGraph(4,zbpm,mean_y);
  grbeamx->SetTitle("; Z BPM (cm); X BPM (cm) (+X beam right) ");
  grbeamx->SetMarkerStyle(20);
  grbeamx->SetMarkerSize(3);
  grbeamy->SetMarkerStyle(21);
  grbeamy->SetMarkerSize(3);
  grbeamy->SetTitle("; Z BPM (cm); Y BPM (cm) (+Y up) ");
 TCanvas *cbeam2 = new TCanvas("cbeam2","Beam plots",800,700);
  cbeam2->Divide(1,2);
  cbeam2->cd(1);
  grbeamx->GetXaxis()->SetLimits(-400,50);
  grbeamx->SetMinimum(-.5);
  grbeamx->SetMaximum(.5);
  grbeamx->Draw("AP");
  cbeam2->cd(2);
  grbeamy->GetXaxis()->SetLimits(-400,50);
 grbeamy->Draw("AP");
  grbeamy->SetMinimum(-.5);
  grbeamy->SetMaximum(.5);
   cbeam2->Print(outputpdf+")");

  //
}
