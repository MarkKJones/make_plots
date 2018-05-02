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

void comp_hms_cer(TString basename, TString basename1) {
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
 outputpdf="plots/"+basename+"_hms_cer.pdf";
  TString inputroot;
   TFile *fhistroot[2];
     inputroot="hist/"+basename+"_cer_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"_cer_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
   TString hname[6]={"cer1_pamp","cer2_pamp","cer1_pint","cer2_pint","cer1_pamp_nocut","cer2_pamp_nocut"};
   TH1F *hist[2][6];
  for (Int_t ifn=0;ifn<2;ifn++) {
  for (Int_t ip=0;ip<6;ip++) {
    // cout <<  hname[ip] << endl;
       hist[ifn][ip] = (TH1F*)fhistroot[ifn]->Get(hname[ip]);
       if (!hist[ifn][ip]) cout << " no hist = " << hname[ip] << endl;
  }}
  // 
  TCanvas *cplot = new TCanvas("cplot"," Cer1 ampl",700,700);
  cplot->Divide(1,2);
  cplot->cd(1);
  hist[0][4]->Draw();
  hist[0][4]->SetAxisRange(0.,50.);
  hist[0][4]->SetTitle("Run 3933 Cer 1 HV = 2590 ");
  cplot->cd(2);
  hist[1][4]->Draw();
  hist[1][4]->SetAxisRange(0.,50.);
  hist[1][4]->SetTitle("Run 3930 Cer 1 HV = 2507 ");
  //
  TCanvas *cplot2 = new TCanvas("cplot2"," Cer2 ampl",700,700);
  cplot2->Divide(1,2);
  cplot2->cd(1);
  hist[0][5]->Draw();
  hist[0][5]->SetAxisRange(0.,50.);
  hist[0][5]->SetTitle("Run 3933 Cer 2 HV = 2590 ");
  cplot2->cd(2);
  hist[1][5]->Draw();
  hist[1][5]->SetAxisRange(0.,50.);
  hist[1][5]->SetTitle("Run 3930 Cer 2 HV = 2573 ");
  //
  TCanvas *cplot3 = new TCanvas("cplot3"," Cer1 int",700,700);
  cplot3->Divide(1,2);
  cplot3->cd(1);
  hist[0][2]->Draw();
  hist[0][2]->SetAxisRange(0.,50.);
  cplot3->cd(2);
  hist[1][2]->Draw();
   hist[1][2]->SetAxisRange(0.,50.);
 //
  TCanvas *cplot4 = new TCanvas("cplot4"," Cer2 int",700,700);
  cplot4->Divide(1,2);
  cplot4->cd(1);
  hist[0][3]->Draw();
  hist[0][3]->SetAxisRange(0.,50.);
  cplot4->cd(2);
  hist[1][3]->Draw();
  hist[1][3]->SetAxisRange(0.,50.);

//
}
