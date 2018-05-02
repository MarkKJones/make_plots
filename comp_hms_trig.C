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

void comp_hms_trig(TString basename, TString basename1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
  static const Int_t adcnum=9;
 const char* adcname[adcnum]={
   "hASUM_adc",
   "hBSUM_adc",
   "hCSUM_adc",
   "hDSUM_adc",
   "hPSHWR_adc",
   "hSHWR_adc",
 "hAER_adc",
   "hCER_adc",
   "hFADC_TREF_ROC1_adc"
 };
 static const Int_t tdcnum=53;
 const char* tdcname[tdcnum]={
   "h1X_tdc",
   "h1Y_tdc",
   "h2X_tdc",
   "h2Y_tdc",
   "h1T_tdc",
   "h2T_tdc",
   "hT1_tdc",
   "hASUM_tdc",
   "hBSUM_tdc",
   "hCSUM_tdc",
   "hDSUM_tdc",
   "hPRLO_tdc",
   "hPRHI_tdc",
   "hSHWR_tdc",
   "hEDTM_tdc",
   "hCER_tdc",
   "hT2_tdc",
    "hDCREF1_tdc",
   "hDCREF2_tdc",
   "hDCREF3_tdc",
   "hDCREF4_tdc",
   "hTRIG1_tdc",
   "hTRIG2_tdc",
   "hTRIG3_tdc",
   "hTRIG4_tdc",
   "hTRIG5_tdc",
   "hTRIG6_tdc",
   "pTRIG1_tdc",
   "pTRIG2_tdc",
   "pTRIG3_tdc",
   "pTRIG4_tdc",
   "pTRIG5_tdc",
   "pTRIG6_tdc",
   "pSTOF_tdc",
   "pEL_LO_LO_tdc",
   "pEL_LO_tdc",
   "pEL_HI_tdc",
   "pEL_REAL_tdc",
   "pEL_CLEAN_tdc",
   "hSTOF_tdc",
   "hEL_LO_LO_tdc",
   "hEL_LO_tdc",
   "hEL_HI_tdc",
   "hEL_REAL_tdc",
   "hEL_CLEAN_tdc",
   "pPRE40_tdc",
   "pPRE100_tdc",
   "pPRE150_tdc",
   "pPRE200_tdc",
   "hPRE40_tdc",
   "hPRE100_tdc",
   "hPRE150_tdc",
   "hPRE200_tdc"
};
//
     TString outputpdf;
 outputpdf="plots/"+basename+"_hms_trig.pdf";
  TString inputroot;
   TFile *fhistroot[2];
     inputroot="hist/"+basename+"_trig_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
     inputroot="hist/"+basename1+"_trig_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[1] =  new TFile(inputroot);
  TH1F* tdchist[2][tdcnum];
  for (Int_t ifn=0;ifn<2;ifn++) {
  for (Int_t ip=0;ip<tdcnum;ip++) {
       TString hname= Form("tdchist_%s",tdcname[ip]);
       tdchist[ifn][ip] = (TH1F*)fhistroot[ifn]->Get(hname);
       if (!tdchist[ifn][ip]) cout << " no hist = " << hname << endl;
  }}
  // 
  const static Int_t numcan=8;
  Int_t nplots=tdcnum/numcan+1;
  Int_t plcnt=0;
  cout << " nplots = " << nplots << endl;
    TCanvas *cplot[numcan];
     for (Int_t nh=0;nh<numcan;nh++) {
       cplot[nh] = new TCanvas(Form("cplot_%d",nh),Form("plot_%d",nh), 700,700);
       cplot[nh]->Divide(2,int(nplots/2));
        for (Int_t nl=0;nl<nplots;nl++) {
	  cplot[nh]->cd(nl+1);
	  if (plcnt<tdcnum) {
            gPad->SetLogy();
            if (tdchist[0][plcnt]) tdchist[0][plcnt]->Draw();
            if (tdchist[1][plcnt]) tdchist[1][plcnt]->Draw("same");
            if (tdchist[1][plcnt]) tdchist[1][plcnt]->SetLineColor(2);
	    plcnt++;
	  }
	}
     }
//
}
