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

void make_hist_shms_hgcer_ll(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_hgcer_ll_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t npeSum;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&npeSum) ;
 Double_t Amp[4];
   tsimc->SetBranchAddress("P.hgcer.goodAdcPulseAmp",&Amp) ;
   //
   TString temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
   temp=Form("Run %d ; Pulse Amp cer0  ; Counts",nrun);
   TH1F *hcer0PulseAmp= new TH1F("hcer0PulseAmp",temp,50,0,500.);
   HList.Add(hcer0PulseAmp);
   temp=Form("Run %d ; Pulse Amp cer1  ; Counts",nrun);
   TH1F *hcer1PulseAmp= new TH1F("hcer1PulseAmp",temp,50,0,500.);
  HList.Add(hcer1PulseAmp);
  temp=Form("Run %d ; Pulse Amp cer2  ; Counts",nrun);
   TH1F *hcer2PulseAmp= new TH1F("hcer2PulseAmp",temp,50,0,500.);
   HList.Add(hcer2PulseAmp);
   temp=Form("Run %d ; Pulse Amp cer3  ; Counts",nrun);
   TH1F *hcer3PulseAmp= new TH1F("hcer3PulseAmp",temp,50,0,500.);
  HList.Add(hcer3PulseAmp);
   //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		hcernpeSum->Fill(npeSum);
		hcer0PulseAmp->Fill(Amp[0]);
		hcer1PulseAmp->Fill(Amp[1]);
		hcer2PulseAmp->Fill(Amp[2]);
		hcer3PulseAmp->Fill(Amp[3]);
	}




}
