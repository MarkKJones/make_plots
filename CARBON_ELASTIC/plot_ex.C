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
#include "TF1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <TCutG.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_ex(Int_t nrun=2043,Int_t ntot=100000){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   TString basename;
   basename=Form("hms_replay_production_mkj_%d_%d",nrun,ntot);
   if (nrun > 3000) basename=Form("shms_replay_production_mkj_%d_%d",nrun,ntot);
   inputroot= "hist/"+basename+"_helium_elastic_hist.root";
   const UInt_t nplots=3;
   TString hname[nplots]={"hExHe","hExHe_scincut","hExHe_scincut1"};
   TFile *fhistroot;
   fhistroot =  new TFile(inputroot);
   TH1F* fhist[nplots];
    for (UInt_t nh=0;nh<nplots;nh++) {
     fhist[nh] = (TH1F*)fhistroot->Get(hname[nh]);
     if (!fhist[nh]) cout << " No hist = " << hname[nh] << endl;
   }
    TF1* fPeak;
    TF1* fBg;
    TF1* fSum;
    fPeak = new TF1("fPeak","gaus");
    fPeak->SetParameter(0,20);
    fPeak->SetParameter(1,-2);
    fPeak->SetParameter(2,-2);
    fBg = new TF1("fBg","pol1");
    fBg->SetParameter(0,5.77);
    fBg->SetParameter(1,0.14);
    fSum = new TF1("fSum","fPeak+fBg");
    TCanvas *cHMS;
    TLegend *lHMS;
     cHMS = new TCanvas("cHMS","HMS", 1000,700);
     cHMS->Divide(1,2);
     cHMS->cd(0);
     fhist[0]->Draw();
     fhist[0]->Fit("fSum","","",-8,12);
     cHMS->cd(1);
     fhist[1]->Draw();
     fhist[1]->Fit("fSum","","",-8,12);
     //     Double_t nPeak = fPeak->Integral(fhist[0]->FindBin(-4.),fhist[0]->FindBin(4.));
     Double_t exhi=6.;
     Double_t exlo=-4.;
     Double_t exbin=100;
     Double_t bg_con = fBg->GetParameter(0);
     Double_t bg_line = fBg->GetParameter(1);
     Double_t nBg = bg_con*(exhi-exlo)+bg_line*(exhi*exhi-exlo*exlo)/2.;
     Double_t nPeak = 0.;
     cout << nPeak << " " <<  nBg << " " << nPeak-nBg<< endl;
}
