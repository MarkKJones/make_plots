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

void plot_hms_trig(TString basename,  Bool_t setcuts = kFALSE) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
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
 outputpdf="plots/"+basename+"_hms_trig.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_trig_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
      TCutG *cutg[tdcnum];
      Double_t tdcmin[tdcnum];
      Double_t tdcmax[tdcnum];
  TH1F* tdchist[tdcnum];
  for (Int_t ip=0;ip<tdcnum;ip++) {
       TString hname= Form("tdchist_%s",tdcname[ip]);
       tdchist[ip] = (TH1F*)fhistroot->Get(hname);
       if (!tdchist[ip]) cout << " no hist = " << hname << endl;
   }
   //
  if (setcuts) {
    TCanvas *ctrig;
     ctrig = new TCanvas("ctrig","trig", 700,700);
     for (Int_t nh=0;nh<tdcnum;nh++) {
     ctrig->Divide(1,1);
     gPad->SetLogy();
     if (tdchist[nh]) tdchist[nh]->Draw();
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
  gPad->Update();
     TString cutname=Form("cut_%d",nh);
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg[nh]=(TCutG*)(tempg->Clone());
      	  cutg[nh]->SetName(cutname);
	  Int_t npts=cutg[nh]->GetN();
          Double_t ydummy;
	  cutg[nh]->GetPoint(0,tdcmin[nh],ydummy);
	  cutg[nh]->GetPoint(1,tdcmax[nh],ydummy);
    }
     
   }
     
     for (Int_t nh=0;nh<tdcnum;nh++) {
       if (nh==0) cout <<  "hms_TDCTimeWindowMin = ";
		    cout << tdcmin[nh] << "," ;	
		    if (nh!=0&&nh%10==0) cout << endl;
     }
     for (Int_t nh=0;nh<tdcnum;nh++) {
       if (nh==0) cout <<  "hms_TDCTimeWindowMax = ";
		    cout << tdcmax[nh] << "," ;	
		    if (nh!=0&&nh%10==0) cout << endl;
     }
  }
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
            if (tdchist[plcnt]) tdchist[plcnt++]->Draw();
	  }
	}
     }
//
}
