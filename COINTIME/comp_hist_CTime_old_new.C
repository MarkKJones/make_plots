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
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void comp_hist_CTime_old_new(Int_t nr1) {
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
     outputpdf=Form("plots/comp_%d_CTime_old_new.pdf",nr1);
  TString inputroot[2];
   TFile *fhistroot[2];
   const char* lname[2]={"new code","old code"};
   inputroot[0]=Form("hist/coin_replay_production_%d_new_CTime_hist.root",nr1);
   inputroot[1]=Form("hist/coin_replay_production_%d_old_CTime_hist.root",nr1);
     cout << " infile root = " << inputroot[0] << endl;
     cout << " infile root = " << inputroot[1] << endl;
   fhistroot[0] =  new TFile(inputroot[0]);
   fhistroot[1] =  new TFile(inputroot[1]);
   static const Int_t nhist=3;
   static const Int_t nhist2d=4;
   //  TString hname[nhist]={"hepicoinTime_ROC2"};
   TString hname[nhist]={"hpbetatrack","hhbetatrack","hepcoinTime_ROC2"};
   TString hname2d[nhist2d]={"hepcoinTime_hxptar","hepcoinTime_pxptar","hepcoinTime_hdelta","hepcoinTime_pdelta"};
    TH1F *hist[2][nhist];
    TH2F *hist2d[2][nhist2d];
  for (Int_t nf=0;nf<2;nf++) {
  for (Int_t ip=0;ip<nhist;ip++) {
    // cout <<  hname[ip] << endl;
       hist[nf][ip] = (TH1F*)fhistroot[nf]->Get(hname[ip]);
       if (!hist[nf][ip]) cout << " no hist = " << hname[ip] << endl;
       hist[nf][ip]->SetTitle(Form(" Run %d %s",nr1,lname[nf]));
  }
   for (Int_t ip=0;ip<nhist2d;ip++) {
    // cout <<  hname[ip] << endl;
       hist2d[nf][ip] = (TH2F*)fhistroot[nf]->Get(hname2d[ip]);
       if (!hist2d[nf][ip]) cout << " no hist = " << hname2d[ip] << endl;
       hist2d[nf][ip]->SetTitle(Form(" Run %d %s",nr1,lname[nf]));
  }
 }
  //
  TCanvas *cPlot[nhist];
     TF1* fW[2][nhist];
    Double_t sigma[2][nhist];
    Double_t Max;
 for (Int_t ip=0;ip<nhist;ip++) {
    cPlot[ip] = new TCanvas(Form("can_%d",ip),Form("can_%d",ip),700,700);
    Max = hist[0][ip]->GetBinCenter(hist[0][ip]->GetMaximumBin());
cPlot[ip]->Divide(1,2);
  cPlot[ip]->cd(1);
  //  hist[0][ip]->Draw();
  fW[0][ip] = new TF1(Form("fW_0_%d",ip),"gaus",Max*.90,Max*1.1);
  hist[0][ip]->Fit(Form("fW_0_%d",ip),"QR");
  sigma[0][ip]=fW[0][ip]->GetParameter(2);
  if (ip==0) cout << " sigma = " << sigma[0][ip] << " time = " << 200./30.*sigma[0][ip]<< endl;
  if (ip==1) cout << " sigma = " << sigma[0][ip] << " time = " << 300./30.*sigma[0][ip]<< endl;
cPlot[ip]->cd(2);
    Max = hist[1][ip]->GetBinCenter(hist[1][ip]->GetMaximumBin());
  fW[1][ip] = new TF1(Form("fW_0_%d",ip),"gaus",Max*.90,Max*1.1);
  hist[1][ip]->Fit(Form("fW_0_%d",ip),"QR");
  sigma[1][ip]=fW[1][ip]->GetParameter(2);
  if (ip==0) cout << " sigma = " << sigma[1][ip] << " time = " << 200./30.*sigma[1][ip]<< endl;
  if (ip==1) cout << " sigma = " << sigma[1][ip] << " time = " << 300./30.*sigma[1][ip]<< endl;
    // hist[1][ip]->Draw();
  if (ip==0) cPlot[ip]->Print(outputpdf+"(");
   if (ip!=0) cPlot[ip]->Print(outputpdf);
 }
   //
  TCanvas *cPlot2d[nhist2d];
  for (Int_t ip=0;ip<nhist2d;ip++) {
    cPlot2d[ip] = new TCanvas(Form("can2d_%d",ip),Form("can2d_%d",ip),700,700);
  cPlot2d[ip]->Divide(1,2);
  cPlot2d[ip]->cd(1);
  hist2d[0][ip]->Draw("colz");
  cPlot2d[ip]->cd(2);
  hist2d[1][ip]->Draw("colz");
  if (ip==nhist2d-1) cPlot2d[ip]->Print(outputpdf+")");
  if (ip<nhist2d-1) cPlot2d[ip]->Print(outputpdf);
  }
  
   //
}
