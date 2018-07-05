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

void make_hist_beta(TString basename="",Int_t nrun=3288,Double_t pfac=0.985){
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
   outputhist= "hist/"+basename+"_beta_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  betanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&betanotrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress("H.hod.beta",&betatrack);
 Double_t  starttime;
   tsimc->SetBranchAddress("H.hod.starttime",&starttime);
 Double_t  pbetanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&pbetanotrack);
 Double_t  pbetatrack;
   tsimc->SetBranchAddress("P.hod.beta",&pbetatrack);
 Double_t  pstarttime;
   tsimc->SetBranchAddress("P.hod.starttime",&pstarttime);
   // Define histograms
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetatrack);
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime);
    TH1F *hpbetanotrack = new TH1F("hpbetanotrack",Form("Run %d ;SHMS Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetanotrack);
    TH1F *hpbetatrack = new TH1F("hpbetatrack",Form("Run %d ;SHMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrack);
    TH1F *hpstarttime = new TH1F("hpstarttime",Form("Run %d ;SHMS Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hpstarttime);
  // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hbetanotrack->Fill(betanotrack);
		hbetatrack->Fill(betatrack);
		hstarttime->Fill(starttime);
		hpbetanotrack->Fill(pbetanotrack);
		hpbetatrack->Fill(pbetatrack);
		hpstarttime->Fill(pstarttime);
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
