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
#include <math.h>
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
using namespace std;

void make_hist_raster(TString basename="",Int_t nrun=3288){
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
   outputhist= "hist/"+basename+"_raster_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t  rastx;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&rastx);
 Double_t  rasty;
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&rasty);
 Double_t  rastxaRawAdc;
   tsimc->SetBranchAddress("P.rb.raster.frxaRawAdc",&rastxaRawAdc);
 Double_t  rastyaRawAdc;
   tsimc->SetBranchAddress("P.rb.raster.fryaRawAdc",&rastyaRawAdc);
 Double_t  rastxbRawAdc;
   tsimc->SetBranchAddress("P.rb.raster.frxbRawAdc",&rastxbRawAdc);
 Double_t  rastybRawAdc;
   tsimc->SetBranchAddress("P.rb.raster.frybRawAdc",&rastybRawAdc);
   //
     TH1F *hxrast = new TH1F("hxrast",Form("Run %d ; Xrast (cm) (+X beam right);Counts",nrun), 300, -.3,.3);
    HList.Add(hxrast);
   TH1F *hyrast = new TH1F("hyrast",Form("Run %d ; Yrast (cm) (+Y up) ;Counts",nrun), 300, -.3,.3);
    HList.Add(hyrast);
     TH1F *hrastxaRawAdc = new TH1F("hrastxaRawAdc",Form("Run %d ; A_X raw  ;Counts",nrun), 1000,0.,100000.);
    HList.Add(hrastxaRawAdc);
     TH1F *hrastyaRawAdc = new TH1F("hrastyaRawAdc",Form("Run %d ; A_Y raw  ;Counts",nrun), 1000,0.,100000.);
    HList.Add(hrastyaRawAdc);
     TH1F *hrastxbRawAdc = new TH1F("hrastxbRawAdc",Form("Run %d ; B_X raw  ;Counts",nrun), 1000,0.,100000.);
    HList.Add(hrastxbRawAdc);
     TH1F *hrastybRawAdc = new TH1F("hrastybRawAdc",Form("Run %d ; B_Y raw  ;Counts",nrun), 1000,0.,100000.);
    HList.Add(hrastybRawAdc);
     TH1F *hrastxaVolts = new TH1F("hrastxaVolts",Form("Run %d ; A_X raw (V) ;Counts",nrun), 500, 0. ,1.);
    HList.Add(hrastxaVolts);
     TH1F *hrastyaVolts = new TH1F("hrastyaVolts",Form("Run %d ; A_Y raw (V) ;Counts",nrun), 500, 0. ,1.);
    HList.Add(hrastyaVolts);
     TH1F *hrastxbVolts = new TH1F("hrastxbVolts",Form("Run %d ; B_X raw (V) ;Counts",nrun), 500, 0. ,1.);
    HList.Add(hrastxbVolts);
     TH1F *hrastybVolts = new TH1F("hrastybVolts",Form("Run %d ; B_Y raw (V) ;Counts",nrun), 500, 0. ,1.);
    HList.Add(hrastybVolts);
 // loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hrastxaRawAdc->Fill(rastxaRawAdc);
		hrastxbRawAdc->Fill(rastxbRawAdc);
		hrastyaRawAdc->Fill(rastyaRawAdc);
		hrastybRawAdc->Fill(rastybRawAdc);
		hrastxaVolts->Fill(rastxaRawAdc*(1./4096.)/25.);
		hrastxbVolts->Fill(rastxbRawAdc*(1./4096.)/25.);
		hrastyaVolts->Fill(rastyaRawAdc*(1./4096.)/25.);
		hrastybVolts->Fill(rastybRawAdc*(1./4096.)/25.);
	}
//
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
}
