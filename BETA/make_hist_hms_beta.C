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

void make_hist_hms_beta(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"hms_beta_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  betanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&betanotrack);
 Double_t  betachisqnotrack;
   tsimc->SetBranchAddress("H.hod.betachisqnotrack",&betachisqnotrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress("H.hod.beta",&betatrack);
 Double_t  etotnorm;
 tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  s1x_x;
   tsimc->SetBranchAddress("H.hod.1x.TrackXPos",&s1x_x);
 Double_t  s1x_y;
   tsimc->SetBranchAddress("H.hod.1x.TrackYPos",&s1x_y);
 Double_t  s2x_x;
   tsimc->SetBranchAddress("H.hod.2x.TrackXPos",&s2x_x);
 Double_t  s2x_y;
   tsimc->SetBranchAddress("H.hod.2x.TrackYPos",&s2x_y);
 Double_t  ntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&ntr);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  starttime;
   tsimc->SetBranchAddress("H.hod.starttime",&starttime);
 Double_t  fptime;
    tsimc->SetBranchAddress("H.hod.fpHitsTime",&fptime);
 Double_t  TimeHist_Hits;
   tsimc->SetBranchAddress("H.hod.TimeHist_Hits",&TimeHist_Hits);
 Double_t  TimeHist_Sigma;
   tsimc->SetBranchAddress("H.hod.TimeHist_Sigma",&TimeHist_Sigma);
 Double_t  TimeHist_Peak;
   tsimc->SetBranchAddress("H.hod.TimeHist_Peak",&TimeHist_Peak);
   // Define histograms
   TH1F *hTimeHist_Hits = new TH1F("hTimeHist_Hits",Form("Run %d ; Time Histogram Nhits;Counts",nrun), 20,0,20);
    HList.Add(hTimeHist_Hits);
   TH1F *hTimeHist_Sigma = new TH1F("hTimeHist_Sigma",Form("Run %d ; Time Histogram Sigma;Counts",nrun), 80,0,20);
    HList.Add(hTimeHist_Sigma);
   TH1F *hTimeHist_Peak = new TH1F("hTimeHist_Hits_Peak",Form("Run %d ; Time Histogram Peak (ns);Counts",nrun), 120,0,60);
    HList.Add(hTimeHist_Peak);
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack);
    TH1F *hbetachisqnotrack = new TH1F("hbetachisqnotrack",Form("Run %d ; Beta chisq notrack;Counts",nrun), 15, -5.,10.0);
    HList.Add(hbetachisqnotrack);
    TH2F *hbetanotrack_delta = new TH2F("hbetanotrack_delta",Form("Run %d ; Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hbetanotrack_delta);
    TH2F *hbetanotrack_hits = new TH2F("hbetanotrack_hits",Form("Run %d ; Beta notrack;Time Histogram Nhits",nrun), 200, .5,1.5,20,0,20);
    HList.Add(hbetanotrack_hits);
    TH2F *hbetatrack_hits = new TH2F("hbetatrack_hits",Form("Run %d ; Beta notrack;Time Histogram Nhits",nrun), 200, -.1,1.5,20,0,20);
    HList.Add(hbetatrack_hits);
    TH2F *hs1x_x_y = new TH2F("hs1x_x_y",Form("Run %d ; Y at S1x ; x at S1x",nrun),100,-40,40,100,-20,20);
    HList.Add(hs1x_x_y);
    TH2F *hs2x_x_y = new TH2F("hs2x_x_y",Form("Run %d ; Y at S2x ; x at S2x",nrun),100,-40,40,100,-20,20);
    HList.Add(hs2x_x_y);
    TH2F *hbetanotrack_xfp = new TH2F("hbetanotrack_xfp",Form("Run %d ; Beta notrack;X_fp",nrun), 200, .5,1.5,100,-40,40);
    HList.Add(hbetanotrack_xfp);
    TH2F *hbetanotrack_yfp = new TH2F("hbetanotrack_yfp",Form("Run %d ; Beta notrack;Y_fp",nrun), 200, .5,1.5,100,-20,20);
    HList.Add(hbetanotrack_yfp);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack);
    TH1F *hbetatrack_LowSigma = new TH1F("hbetatrack_LowSigma",Form("Run %d ; Beta track _LowSigma;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack_LowSigma);
    TH1F *hbetatrack_HighSigma = new TH1F("hbetatrack_HighSigma",Form("Run %d ; Beta track _HighSigma;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack_HighSigma);
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime);
    TH1F *hfptime = new TH1F("hfptime",Form("Run %d ; Fp time;Counts",nrun), 280, -10.,60.0);
    HList.Add(hfptime);
    TH1F *hstarttime_betazero = new TH1F("hstarttime_betazero",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime_betazero);
    TH1F *hfptime_betazero = new TH1F("hfptime_betazero",Form("Run %d ; Fp time;Counts",nrun), 280, -10.,60.0);
    HList.Add(hfptime_betazero);
 // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (1==1) {
		  hbetanotrack_hits->Fill(betanotrack,TimeHist_Hits);
		  if (ntr >0)  hbetatrack_hits->Fill(betatrack,TimeHist_Hits);
		  hTimeHist_Sigma->Fill(TimeHist_Sigma);
		  hTimeHist_Hits->Fill(TimeHist_Hits);
		  hTimeHist_Peak->Fill(TimeHist_Peak);
		if (ntr >0 ) hbetatrack->Fill(betatrack);
		if (ntr >0 && TimeHist_Sigma<0.6) hbetatrack_LowSigma->Fill(betatrack);
		if (ntr >0 && TimeHist_Sigma>=0.6) hbetatrack_HighSigma->Fill(betatrack);
                hbetanotrack_xfp->Fill(betanotrack,xfp);
                hbetanotrack_yfp->Fill(betanotrack,yfp);
		if (ntr>0) {
		hbetanotrack->Fill(betanotrack);
                hbetanotrack_delta->Fill(betanotrack,delta);
                hs1x_x_y->Fill(s1x_x,s1x_y);
                hs2x_x_y->Fill(s2x_x,s2x_y);
		hstarttime->Fill(starttime);
		hfptime->Fill(fptime);
                if (betatrack==0) hstarttime_betazero->Fill(starttime);
                if (betatrack==0) hfptime_betazero->Fill(fptime);
                if (betanotrack==0) hbetachisqnotrack->Fill(betachisqnotrack);
		}
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
