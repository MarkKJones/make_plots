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
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
using namespace std;

void make_hist_check_hodo_time_cuts(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"_check_hodo_time_cuts_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  betanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&betanotrack);
 Double_t  ntrack;
   tsimc->SetBranchAddress("H.dc.ntrack",&ntrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress("H.hod.beta",&betatrack);
 Double_t  hodTrackXPos;
   tsimc->SetBranchAddress("H.hod.1x.TrackXPos",&hodTrackXPos);
 Double_t  hodTrackYPos;
   tsimc->SetBranchAddress("H.hod.1x.TrackYPos",&hodTrackYPos);
 Double_t  starttime;
   tsimc->SetBranchAddress("H.hod.starttime",&starttime);
 Double_t  TimeHist_Hits;
 tsimc->SetBranchAddress("H.hod.TimeHist_Hits",&TimeHist_Hits);
 Double_t  TimeHist_Sigma;
 tsimc->SetBranchAddress("H.hod.TimeHist_Sigma",&TimeHist_Sigma);
 Double_t  TimeHist_Mean;
 tsimc->SetBranchAddress("H.hod.TimeHist_Mean",&TimeHist_Mean);
   Double_t  GoodPosTofCorr[4][16];
   Double_t  GoodNegTofCorr[4][16];
   TString plname[4]={"1x","1y","2x","2y"};
   Int_t plnum[4]={16,10,16,10};
   for (Int_t np=0;np<4;np++) {
     TString name = "H.hod."+plname[np]+".GoodPosTdcTimeTOFCorr";
     tsimc->SetBranchAddress(name,&GoodPosTofCorr[np]);
     name = "H.hod."+plname[np]+".GoodNegTdcTimeTOFCorr";
     tsimc->SetBranchAddress(name,&GoodNegTofCorr[np]);
   }
   // Define histograms
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack);
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 80, 0,40.0);
    HList.Add(hstarttime);
    TH1F *hsigma = new TH1F("hsigma",Form("Run %d ; Time Hist SIgma;Counts",nrun), 80, 0.,4.0);
    HList.Add(hsigma);
    TH1F *hghits = new TH1F("hghits",Form("Run %d ; Time Hist Hits;Counts",nrun), 40, 0.,40.0);
    HList.Add(hghits);
    TH1F *htrackxpos = new TH1F("htrackxpos",Form("Run %d ; Track Xpos at Hod 1x;Counts",nrun), 100, -60.,60.0);
    HList.Add(htrackxpos);
    TH1F *htrackypos = new TH1F("htrackypos",Form("Run %d ; Track Ypos at Hod 1x;Counts",nrun), 100, -60.,60.0);
    HList.Add(htrackypos);
    TH1F *hghits_betazero = new TH1F("hghits_betazero",Form("Run %d ; Time Hist Hits (beta=0);Counts",nrun), 40, 0.,40.0);
    HList.Add(hghits_betazero);
    TH2F *hsigma_betatrack = new TH2F("hsigma_betatrack",Form("Run %d ; Time Hist SIgma; Beta track",nrun), 80, 0.,4.0,200, -1.,2.0);
    TH2F *hstarttime_betatrack = new TH2F("hstarttime_betatrack",Form("Run %d ; Starttime; Beta track",nrun), 80, 0.,40.0,200, -1.,2.0);
    HList.Add(hsigma_betatrack);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 600, -1.,2.0);
  TH1F *hTime = new TH1F("hTime","",400,0,200);
  // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
   for (Int_t np=0;np<4;np++) {
   for (Int_t npad=0;npad<plnum[np];npad++) {
     if (GoodPosTofCorr[np][npad]<1000 || GoodNegTofCorr[np][npad]<1000) {
       hTime->Fill(GoodPosTofCorr[np][npad]);
       hTime->Fill(GoodNegTofCorr[np][npad]);
     }
   }
   } 		
   if (ntrack>0) {
   hbetanotrack->Fill(betanotrack);
   hbetatrack->Fill(betatrack);
   hstarttime->Fill(starttime);
   hsigma->Fill(TimeHist_Sigma);
   hghits->Fill(TimeHist_Hits);
   if (betatrack==0) hghits_betazero->Fill(TimeHist_Hits);
   if (betatrack==0) htrackxpos->Fill(hodTrackXPos);
   if (betatrack==0) htrackypos->Fill(hodTrackYPos);
   hsigma_betatrack->Fill(TimeHist_Sigma,betatrack);
   hstarttime_betatrack->Fill(starttime,betatrack);
   }
   hTime->Reset();
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
