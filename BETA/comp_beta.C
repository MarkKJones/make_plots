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

void comp_beta(TString basename, TString basename1) {
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
 outputpdf="plots/"+basename+"_hms_beta.pdf";
  TString inputroot;
   TFile *fhistroot[2];
     inputroot="hist/"+basename+"hms_beta_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"hms_beta_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
   TString hname[4]={"hbetanotrack","hbetatrack","hstarttime","hfptime"};
   TH1F *hist[2][4];
  for (Int_t ifn=0;ifn<2;ifn++) {
  for (Int_t ip=0;ip<4;ip++) {
    // cout <<  hname[ip] << endl;
       hist[ifn][ip] = (TH1F*)fhistroot[ifn]->Get(hname[ip]);
       if (!hist[ifn][ip]) cout << " no hist = " << hname[ip] << endl;
  }}
  // 
  TCanvas *cplot = new TCanvas("cplot","HMS hodo",700,700);
  cplot->Divide(2,2);
  cplot->cd(1);
gPad->SetLogy();
  hist[0][0]->Draw();
  hist[1][0]->Draw("same");
  hist[1][0]->SetLineColor(2);
  cplot->cd(2);
gPad->SetLogy();
  hist[0][1]->Draw();
  hist[1][1]->Draw("same");
  hist[1][1]->SetLineColor(2);
  cplot->cd(3);
gPad->SetLogy();
  hist[0][2]->Draw();
  hist[1][2]->Draw("same");
  hist[1][2]->SetLineColor(2);
  cout << " Mean for first file = " << hist[0][2]->GetMean() << endl;
  cout << " Mean for 2nd file = " << hist[1][2]->GetMean() << endl;
  cplot->cd(4);
gPad->SetLogy();
  hist[0][3]->Draw();
  hist[1][3]->Draw("same");
  hist[1][3]->SetLineColor(2);
  // 


//
}
