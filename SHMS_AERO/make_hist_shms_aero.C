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

void make_hist_shms_aero(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_aero_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t npad=7;
 static const Int_t iside=2;
 const char* sidename[iside]={"Neg","Pos"};
 const char* sname[iside]={"neg","pos"};
 //
 Double_t adctdcdiff[iside][npad];
 Double_t adcmult[iside][npad];
 Double_t npe[iside][npad];
 for (Int_t is=0;is<iside;is++) {
   tsimc->SetBranchAddress(Form("P.aero.good%sAdcTdcDiffTime",sidename[is]),&adctdcdiff[is]) ;
   tsimc->SetBranchAddress(Form("P.aero.good%sAdcMult",sidename[is]),&adcmult[is]) ;
   tsimc->SetBranchAddress(Form("P.aero.%sNpe",sname[is]),&npe[is]) ;
 }
 Double_t xAtAero;
 Double_t yAtAero; 
 Double_t ntr; 
 tsimc->SetBranchAddress("P.aero.xAtAero",&xAtAero);
 tsimc->SetBranchAddress("P.aero.yAtAero",&yAtAero);
 tsimc->SetBranchAddress("P.dc.ntrack",&ntr);
 //
 TH1F *hAero_npe[iside][npad];
 TH1F *hAero_mult[iside][npad];
  TH1F *hAero_adctdcdiff[iside][npad];
TH2F *hAero_adctdcdiff_x[iside][npad];
 TH2F *hAero_adctdcdiff_y[iside][npad];
TH2F *hAero_npe_x[iside][npad];
 TH2F *hAero_npe_y[iside][npad];
 Double_t lo[iside]={150,-50};
 Double_t hi[iside]={250,50};
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad;ipad++) {
   hAero_npe[is][ipad]= new TH1F(Form("hAero_npe_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; Npe; Counts",sidename[is],ipad),80,0,10.);
   hAero_mult[is][ipad]= new TH1F(Form("hAero_mult_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; Mult; Counts",sidename[is],ipad),10,0,10.);
    hAero_npe_x[is][ipad]= new TH2F(Form("hAero_npe_x_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; Npe; X",sidename[is],ipad),80,0,10.,100,-60,60);
    hAero_npe_y[is][ipad]= new TH2F(Form("hAero_npe_y_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; Npe; Y",sidename[is],ipad),80,0,10.,100,-40,40);
 hAero_adctdcdiff[is][ipad]= new TH1F(Form("hAero_tdiff_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; ADc TDC diff time; Counts",sidename[is],ipad),100,lo[is],hi[is]);
   hAero_adctdcdiff_x[is][ipad]= new TH2F(Form("hAeroX_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; ADc TDC diff time; X ",sidename[is],ipad),100,lo[is],hi[is],100,-60,60);
   hAero_adctdcdiff_y[is][ipad]= new TH2F(Form("hAeroY_%s_pad_%d",sidename[is],ipad),Form("%s pad_%d; ADc TDC diff time; Y ",sidename[is],ipad),100,lo[is],hi[is],100,-40,40);
   HList.Add(hAero_npe[is][ipad]);
   HList.Add(hAero_mult[is][ipad]);
   HList.Add(hAero_npe_x[is][ipad]);
   HList.Add(hAero_npe_y[is][ipad]);
   HList.Add(hAero_adctdcdiff[is][ipad]);
    HList.Add(hAero_adctdcdiff_x[is][ipad]);
    HList.Add(hAero_adctdcdiff_y[is][ipad]);
}}  
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad;ipad++) {
   if (adcmult[is][ipad] == 1) hAero_adctdcdiff[is][ipad]->Fill(adctdcdiff[is][ipad]);
   if (adcmult[is][ipad] >= 1) hAero_npe[is][ipad]->Fill(npe[is][ipad]);
   hAero_mult[is][ipad]->Fill(adcmult[is][ipad]);
   if (ntr>0 && adcmult[is][ipad] >=1 ) {
   hAero_adctdcdiff_x[is][ipad]->Fill(adctdcdiff[is][ipad],xAtAero);
  hAero_adctdcdiff_y[is][ipad]->Fill(adctdcdiff[is][ipad],yAtAero);
   hAero_npe_x[is][ipad]->Fill(npe[is][ipad],xAtAero);
   hAero_npe_y[is][ipad]->Fill(npe[is][ipad],yAtAero);
   }
  }}  
		
	}
TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
