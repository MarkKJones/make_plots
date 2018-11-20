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

void plot_hist_CTime(TString basename) {
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
 outputpdf="plots/"+basename+"_CTime.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   static const Int_t nhist2d=8;
   TString hname2d[nhist2d]={"hptrig4_ptrig1_ROC1_cut","hptrig4_ptrigdiff_ROC1_cut","hpxfp_ptrigdiff_ROC1_cut","hpyfp_ptrigdiff_ROC1_cut","hhxfp_ptrigdiff_ROC1_cut","hhyfp_ptrigdiff_ROC1_cut","hptrig14diff_pstarttimediff_ROC1_cut","hptrig13diff_pstarttimediff_ROC1_cut"};
   TH2F *hist2d[nhist2d];
   for (Int_t ip=0;ip<nhist2d;ip++) {
    // cout <<  hname[ip] << endl;
       hist2d[ip] = (TH2F*)fhistroot->Get(hname2d[ip]);
       if (!hist2d[ip]) cout << " no hist = " << hname2d[ip] << endl;
  }
   //
   TCanvas *can2d[nhist2d];
   for (Int_t ip=0;ip<nhist2d;ip++) {
     can2d[ip] = new TCanvas(Form("can2d_%d",ip),Form("can2d_%d",ip),700,700);
     can2d[ip]->Divide(1,1);
     can2d[ip]->cd(1);
     hist2d[ip]->Draw("colz");
     if (ip==0) can2d[ip]->Print(outputpdf+"(");
     if (ip>0&&ip<nhist2d-1) can2d[ip]->Print(outputpdf);
     if (ip==nhist2d-1) can2d[ip]->Print(outputpdf+")");
   }   
   //
}
