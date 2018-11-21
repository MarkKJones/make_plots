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

void plot_hist_CTime_planes(TString basename) {
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
 outputpdf="plots/"+basename+"_CTime_planes.pdf";
 TString  inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_CTime_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
  static const Int_t plnum=4;
 static const Int_t iside=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 const char* sidename[iside]={"Neg","Pos"};
 const char* sname[iside]={"neg","pos"};
 static const Int_t npad[plnum]={16,10,16,10};
 static const Int_t npadshms[plnum]={13,13,14,21};
 TH1F *htdc1d_plane[plnum];
 TH2F *htdc2d_plane[plnum];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
       htdc1d_plane[ipl] = (TH1F*)fhistroot->Get(Form("hhs%s_pre40_diff",plname[ipl]));
       if (!htdc1d_plane[ipl]) cout << " no " << Form("hhs%s_pre40_diff",plname[ipl]) << endl;
 }
 for (Int_t ipl=0;ipl<3;ipl++) {
       htdc2d_plane[ipl] = (TH2F*)fhistroot->Get(Form("hhs1x_s%s_pre40_diff",plname[ipl+1]));
       if (!htdc2d_plane[ipl]) cout << " no " << Form("hhs1x_s%s_pre40_diff",plname[ipl+1]) << endl;
 }
 //
 TCanvas *can2d = new TCanvas("can2d","can2d",700,700);
 can2d->Divide(2,2);
 for (Int_t ipl=0;ipl<3;ipl++) {
   can2d->cd(ipl+1);
   htdc2d_plane[ipl]->Draw("colz");
 }
 can2d->Print(outputpdf+"(");
  //
 TCanvas *can1d = new TCanvas("can1d","can1d",700,700);
 can1d->Divide(2,2);
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   can1d->cd(ipl+1);
   htdc1d_plane[ipl]->Draw();
 }
 can1d->Print(outputpdf+")");
 
 //
}
