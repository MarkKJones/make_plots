#include <TSystem.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TLine.h>
#include <TObjArray.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void comp_hist_hms_sieve(TString basename="",TString basename2=""){
  //
   if (basename=="") {
     cout << " Input the basename of the root file" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot[2];
   inputroot[0]="hist/"+basename+"_hms_sieve_hist.root";
   inputroot[1]="hist/"+basename2+"_hms_sieve_hist.root";
   TString pdffile;
  TFile *fhistroot[2];
   fhistroot[0] =  new TFile(inputroot[0]);
   fhistroot[1] =  new TFile(inputroot[1]);
    TH1F *hys[2][3];
   TH1F *hxs[2][3];
   TH1F *hyptar[2][3];
   TH1F *hxptar[2][3];
   for (int nf=0;nf<2;nf++) {
   for (int nd=0;nd<3;nd++) {
        hys[nf][nd] = (TH1F*)fhistroot[nf]->Get(Form("hys_%d",nd));
        hxs[nf][nd] = (TH1F*)fhistroot[nf]->Get(Form("hxs_%d",nd));
        hyptar[nf][nd] = (TH1F*)fhistroot[nf]->Get(Form("hyptar_%d",nd));
        hxptar[nf][nd] = (TH1F*)fhistroot[nf]->Get(Form("hxptar_%d",nd));
  }
   //
//
      TCanvas *cplot1 = new TCanvas("cplot1","ang plots", 900,700);
      cplot1->Divide(2,2);
      cplot1->cd(1);
      hxptar[0][1]->Draw();
      hxptar[1][1]->Draw("same");
      hxptar[1][1]->SetLineColor(2);
       cplot1->cd(2);
      hyptar[0][1]->Draw();
      hyptar[1][1]->Draw("same");
      hyptar[1][1]->SetLineColor(2);
       cplot1->cd(3);
      hxs[0][1]->Draw();
      hxs[1][1]->Draw("same");
      hxs[1][1]->SetLineColor(2);
       cplot1->cd(4);
      hys[0][1]->Draw();
      hys[1][1]->Draw("same");
      hys[1][1]->SetLineColor(2);
       //
}
