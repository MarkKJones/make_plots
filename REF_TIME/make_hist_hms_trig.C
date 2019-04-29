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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_hms_trig(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_trig_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
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
 Double_t tdcdata[tdcnum];
 Double_t tdcrawdata[tdcnum];
 Double_t tdcmultdata[tdcnum];
 for (Int_t ip=0;ip<tdcnum;ip++) {
   tsimc->SetBranchAddress(Form("T.shms.%sTime",tdcname[ip]),&tdcdata[ip]) ;
   tsimc->SetBranchAddress(Form("T.shms.%sTimeRaw",tdcname[ip]),&tdcrawdata[ip]) ;
   tsimc->SetBranchAddress(Form("T.shms.%sMultiplicity",tdcname[ip]),&tdcmultdata[ip]) ;
 }
 Double_t adcdata[tdcnum];
 Double_t adcrawdata[tdcnum];
 Double_t adcmultdata[tdcnum];
 for (Int_t ip=0;ip<adcnum;ip++) {
   tsimc->SetBranchAddress(Form("T.shms.%sPulseTime",adcname[ip]),&adcdata[ip]) ;
   tsimc->SetBranchAddress(Form("T.shms.%sPulseTimeRaw",adcname[ip]),&adcrawdata[ip]) ;
   tsimc->SetBranchAddress(Form("T.shms.%sMultiplicity",adcname[ip]),&adcmultdata[ip]) ;
 }
   // Define histograms
 TH1F *tdchist[tdcnum];
 TH1F *tdcrawhist[tdcnum];
 TH1F *tdchistmult1[tdcnum];
 for (Int_t ip=0;ip<tdcnum;ip++) {
   tdchist[ip]= new TH1F(Form("tdchist_%s",tdcname[ip]),Form("; %s (ns) ; Counts ",tdcname[ip]),1000,0,2000);
   HList.Add( tdchist[ip]);
   tdcrawhist[ip]= new TH1F(Form("tdcrawhist_%s",tdcname[ip]),Form(";Raw  %s (ns) ; Counts ",tdcname[ip]),1000,0,2000);
   HList.Add( tdcrawhist[ip]);
   tdchistmult1[ip]= new TH1F(Form("tdchistmult1_%s",tdcname[ip]),Form("; %s (ns) Mult = 1; Counts ",tdcname[ip]),1000,0,2000);
   HList.Add( tdchistmult1[ip]);
 }
 TH1F *adchist[tdcnum];
 TH1F *adcrawhist[tdcnum];
 TH1F *adchistmult1[tdcnum];
 for (Int_t ip=0;ip<adcnum;ip++) {
   adchist[ip]= new TH1F(Form("adchist_%s",adcname[ip]),Form("; %s (ns) ; Counts ",adcname[ip]),1000,0,20000);
   HList.Add( adchist[ip]);
   adcrawhist[ip]= new TH1F(Form("adcrawhist_%s",adcname[ip]),Form(";Raw  %s (ns) ; Counts ",adcname[ip]),1000,0,20000);
   HList.Add( adchist[ip]);
     adchistmult1[ip]= new TH1F(Form("adchistmult1_%s",adcname[ip]),Form("; %s (ns) Mult = 1; Counts ",adcname[ip]),1000,0,20000);
   HList.Add( adchistmult1[ip]);
 }
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
                 for (Int_t ip=0;ip<tdcnum;ip++) {
		   tdchist[ip]->Fill(tdcdata[ip]);
		   tdcrawhist[ip]->Fill(tdcrawdata[ip]);
		   if (tdcmultdata[ip]==1) tdchistmult1[ip]->Fill(tdcdata[ip]);
		 }
                 for (Int_t ip=0;ip<adcnum;ip++) {
		   adchist[ip]->Fill(adcdata[ip]);
		   adcrawhist[ip]->Fill(adcrawdata[ip]);
		   if (adcmultdata[ip]==1) adchistmult1[ip]->Fill(adcdata[ip]);
		 }
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
