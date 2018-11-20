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

void make_hist_track_eff(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_track_eff_hist.root";
 TObjArray HList(0);
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t  GoodScinHit;
   tsimc->SetBranchAddress("P.hod.goodscinhit",&GoodScinHit);
 Double_t  TRIG5_ROC2_rawtime;
   tsimc->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw",&TRIG5_ROC2_rawtime);
 Double_t  PreShEplane;
   tsimc->SetBranchAddress("P.cal.pr.eplane",&PreShEplane);
 Double_t  Etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&Etotnorm);
 Double_t  BetaNoTrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&BetaNoTrack);
 Double_t  gind;
   tsimc->SetBranchAddress("P.gtr.index",&gind);
//
 Double_t PreShCut = 0.1;
 Int_t ncut=7;
   TH1F *pPreShEplane = new TH1F("pPreShEplane","; PreShower Energy",100,0,2.);
   TH1F *pHodGood = new TH1F("pHodGood",";Hod Good Flag",10,0,10);
   TH1F *pHodGoodTrack = new TH1F("pHodGoodTrack",";Hod Good Flag (ntrack>0)",10,0,10);
   TH1F *pHodGoodCut[7];
   TH1F *pHodGoodTrackCut[7];
   for (int ii=0 ; ii<ncut ;ii++) {
     pHodGoodCut[ii] = new TH1F(Form("pHodGoodCut_%d",ii),Form(";Hod Good Flag (PrShr > %5.2f",PreShCut+(ii-1)*.1),10,0,10);
     pHodGoodTrackCut[ii] = new TH1F(Form("pHodGoodTrackCut_%d",ii),Form(";Hod Good Flag (ntrack>0 && PrShr > %5.2f ",PreShCut+(ii-1)*.1),10,0,10);
   }
// 
Long64_t nentries = tsimc->GetEntries();
 
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (BetaNoTrack>0.5&&BetaNoTrack<1.4 && TRIG5_ROC2_rawtime>0) {
	  pPreShEplane->Fill(PreShEplane);
	  pHodGood->Fill(GoodScinHit);
	  if (gind>-1) pHodGoodTrack->Fill(GoodScinHit);
          for (int ii=0 ; ii<ncut ;ii++) {
	    if (PreShEplane > PreShCut+(ii-1)*.1) {
	       pHodGoodCut[ii]->Fill(GoodScinHit);
	       if (gind>-1) pHodGoodTrackCut[ii]->Fill(GoodScinHit);
	    }
	  }
		}
	}
	  //
	    cout << " Track efficiency  = " << pHodGoodTrack->Integral(2,2)/pHodGood->Integral(2,2) << endl;
          for (int ii=0 ; ii<ncut ;ii++) {
	    cout << " Preshower cut = " << PreShCut+(ii-1)*.1 << endl;
	cout << " Track efficiency (Preshower cut) = " << pHodGoodTrackCut[ii]->Integral(2,2)/pHodGoodCut[ii]->Integral(2,2) << endl;
	  }
}
