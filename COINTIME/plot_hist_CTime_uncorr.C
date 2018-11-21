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

void plot_hist_CTime_uncorr(TString basename) {
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
 outputpdf="plots/"+basename+"_CTime_uncorr.pdf";
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
TH1F *htdc_uncorr[plnum][iside][16];
 TH2F *htdc_uncorr_pm[plnum][16];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
    // cout <<  hname[ip] << endl;
       htdc_uncorr_pm[ipl][ipad] = (TH2F*)fhistroot->Get(Form("tdc_uncorr_pm_%s_pad_%d",plname[ipl],ipad));
       if (!htdc_uncorr_pm[ipl][ipad]) cout << " no 2dhist = " << ipl << " " << ipad << endl;
 for (Int_t is=0;is<iside;is++) {
   htdc_uncorr[ipl][is][ipad] = (TH1F*)fhistroot->Get(Form("tdc_uncorr_%s_%s_pad_%d",plname[ipl],sidename[is],ipad));
   if (!htdc_uncorr[ipl][is][ipad]) cout << " no 1dhist = " << ipl << " " << is << " " << ipad << endl;
 } }}
   //
   TCanvas *can2d[4];
   Int_t npl=0;
   Int_t padnum=0;
   for (Int_t ip=0;ip<4;ip++) {
     can2d[ip] = new TCanvas(Form("can2d_hms_1x_%d",ip),Form("can2d_hms_1x_%d",ip),700,700);
     can2d[ip]->Divide(3,4);
     Int_t cn=1;
     for (Int_t ii=0;ii<4;ii++) {
     can2d[ip]->cd(cn++);
     htdc_uncorr_pm[npl][padnum]->Draw("colz");
     can2d[ip]->cd(cn++);
     htdc_uncorr[npl][0][padnum]->Draw();
     can2d[ip]->cd(cn++);
     htdc_uncorr[npl][1][padnum++]->Draw();
   }   
     if (ip==0) can2d[ip]->Print(outputpdf+"(");
     if (ip!=0) can2d[ip]->Print(outputpdf);
     //if (ip==nhist2d-1) can2d[ip]->Print(outputpdf+")");
   }
   //
  //
   TCanvas *can2d_hms_1y[4];
   npl=1;
   padnum=0;
   for (Int_t ip=0;ip<2;ip++) {
     can2d_hms_1y[ip] = new TCanvas(Form("can2d_hms_1y_%d",ip),Form("can2d_hms_1y_%d",ip),700,700);
     can2d_hms_1y[ip]->Divide(3,5);
     Int_t cn=1;
     for (Int_t ii=0;ii<5;ii++) {
     can2d_hms_1y[ip]->cd(cn++);
     htdc_uncorr_pm[npl][padnum]->Draw("colz");
     can2d_hms_1y[ip]->cd(cn++);
     htdc_uncorr[npl][0][padnum]->Draw();
     can2d_hms_1y[ip]->cd(cn++);
     htdc_uncorr[npl][1][padnum++]->Draw();
   }   
   can2d_hms_1y[ip]->Print(outputpdf);
      }
   //
  //
   TCanvas *can2d_hms_2y[4];
   npl=1;
   padnum=0;
   for (Int_t ip=0;ip<2;ip++) {
     can2d_hms_2y[ip] = new TCanvas(Form("can2d_hms_2y_%d",ip),Form("can2d_hms_2y_%d",ip),700,700);
     can2d_hms_2y[ip]->Divide(3,5);
     Int_t cn=1;
     for (Int_t ii=0;ii<5;ii++) {
     can2d_hms_2y[ip]->cd(cn++);
     htdc_uncorr_pm[npl][padnum]->Draw("colz");
     can2d_hms_2y[ip]->cd(cn++);
     htdc_uncorr[npl][0][padnum]->Draw();
     can2d_hms_2y[ip]->cd(cn++);
     htdc_uncorr[npl][1][padnum++]->Draw();
   }   
   can2d_hms_2y[ip]->Print(outputpdf);
      }
   //
   TCanvas *can2d_hms_2x[4];
   npl=2;
   padnum=0;
   for (Int_t ip=0;ip<4;ip++) {
     can2d_hms_2x[ip] = new TCanvas(Form("can2d_hms_2x_%d",ip),Form("can2d_hms_2x_%d",ip),700,700);
     can2d_hms_2x[ip]->Divide(3,4);
     Int_t cn=1;
     for (Int_t ii=0;ii<4;ii++) {
     can2d_hms_2x[ip]->cd(cn++);
     htdc_uncorr_pm[npl][padnum]->Draw("colz");
     can2d_hms_2x[ip]->cd(cn++);
     htdc_uncorr[npl][0][padnum]->Draw();
     can2d_hms_2x[ip]->cd(cn++);
     htdc_uncorr[npl][1][padnum++]->Draw();
   }   
     if (ip<3) can2d_hms_2x[ip]->Print(outputpdf);
    if (ip==3) can2d_hms_2x[ip]->Print(outputpdf+")");
   }
   //
 }
