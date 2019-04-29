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

void make_hist_shms_ngcer_positron(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_ngcer_positron_hist.root";
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
   tsimc->SetBranchAddress("P.ngcer.npeSum",&npeSum) ;
 Double_t xAtCer;
   tsimc->SetBranchAddress("P.ngcer.xAtCer",&xAtCer) ;
 Double_t yAtCer;
   tsimc->SetBranchAddress("P.ngcer.yAtCer",&yAtCer) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntrack) ;
 Double_t xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp) ;
 Double_t etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm) ;
   //
   TString temp=Form("Run %d ; Etottracknorm  ; Counts",nrun);
   TH1F *hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   temp=Form("Run %d ; NpeSUm  ; Counts",nrun);
   TH1F *hcernpeSum = new TH1F("hcernpeSum",temp,160,0,40.);
   HList.Add(hcernpeSum);
   temp=Form("Run %d ; NpeSUm  ; Etottracknorm",nrun);
   TH2F *hcernpeSum_cal = new TH2F("hcernpeSum_cal",temp,160,0,40.,100,0,2.0);
   HList.Add(hcernpeSum_cal);
   temp=Form("Run %d ; NpeSUm  ; X at Cer",nrun);
   TH2F *hcernpeSum_cerX = new TH2F("hcernpeSum_cerX",temp,160,0,40.,12,-60,60);
   HList.Add(hcernpeSum_cerX);
   temp=Form("Run %d ; NpeSUm  ; Xfp",nrun);
   TH2F *hcernpeSum_xfp = new TH2F("hcernpeSum_xfp",temp,160,0,40.,60,-60,60);
   HList.Add(hcernpeSum_xfp);
   temp=Form("Run %d E/p>0.85; NpeSUm  ; Xfp",nrun);
   TH2F *hcernpeSum_xfp_cut = new TH2F("hcernpeSum_xfp_cut",temp,160,0,40.,60,-60,60);
   HList.Add(hcernpeSum_xfp_cut);
   temp=Form("Run %d ; Y at CEr ; X at Cer",nrun);
   TH2F *hcerXY = new TH2F("hcerXY",temp,160,-40,40.,100,-60,60.0);
   HList.Add(hcerXY);
   TH2F *hcerXYnpe = new TH2F("hcerXYnpe",temp,160,-40,40.,100,-60,60.0);
   HList.Add(hcerXYnpe);

//
   Int_t hit_index=0;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		if (ntrack>0) {
		hetotnorm->Fill(etottracknorm);
		hcernpeSum_cal->Fill(npeSum,etottracknorm);
		hcernpeSum_xfp->Fill(npeSum,xfp);
		if (etottracknorm>0.85) hcernpeSum_xfp_cut->Fill(npeSum,xfp);
		hcernpeSum->Fill(npeSum);
	}
	}
 TFile hsimc(outputhist,"recreate");
        HList.Write();
//
}
