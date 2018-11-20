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

void make_hist_shms_beta(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"shms_beta_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  betanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
 Double_t  etottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etottracknorm);
 Double_t  betatrack;
   tsimc->SetBranchAddress("P.hod.beta",&betatrack);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  starttime;
   tsimc->SetBranchAddress("P.hod.starttime",&starttime);
   // Define histograms
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack);
    TH1F *hetottracknorm = new TH1F("hetottracknorm",Form("Run %d ; Etot tracknorm;Counts",nrun), 300, 0.,2.0);
    HList.Add(hetottracknorm);
    TH2F *hbetanotrack_delta = new TH2F("hbetanotrack_delta",Form("Run %d ; Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hbetanotrack_delta);
    TH2F *hbetatrack_delta = new TH2F("hbetatrack_delta",Form("Run %d ; Beta track;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hbetatrack_delta);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack);
    static const Int_t nxpcut=8;
    TH1F *hbeta_xfpcut[nxpcut];
    for (Int_t nx=0;nx<nxpcut;nx++) {
       hbeta_xfpcut[nx] = new TH1F(Form("hbeta_xfpcut_%d",nx),Form("Run %d Cut %d; Beta track;Counts",nrun,nx), 200, 0.5,1.5);
    HList.Add(hbeta_xfpcut[nx]);
    }
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime);
 // loop over entries
    Double_t xlow=-40.;
    Double_t xcent;
    Double_t xfpstep=80./nxpcut;
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hetottracknorm->Fill(etottracknorm);
                if (etottracknorm> .7) {
		hbetanotrack->Fill(betanotrack);
		hbetatrack->Fill(betatrack);
		hstarttime->Fill(starttime);
                hbetanotrack_delta->Fill(betanotrack,delta);
                hbetatrack_delta->Fill(betatrack,delta);
                   for (Int_t nx=0;nx<nxpcut;nx++) {
		     xcent = xlow +nx*xfpstep+xfpstep/2.;
		     if (TMath::Abs(xfp-xcent)<xfpstep) hbeta_xfpcut[nx]->Fill(betatrack); 
		   }
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
