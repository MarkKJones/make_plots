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

void make_hist_vcs2(TString basename="",Int_t nrun=3288,Bool_t vcs_flag=kFALSE){
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
   if (vcs_flag) {
      inputroot="ROOTfiles_vcs/"+basename+".root";
    } else {
      inputroot="ROOTfiles/"+basename+".root";
   }
   TString outputhist;
   outputhist= "hist/"+basename+"_vcs_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot);
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t ctime;
tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&ctime);
 Double_t W;
tsimc->SetBranchAddress("P.kin.primary.W",&W);
 Double_t Q2;
tsimc->SetBranchAddress("P.kin.primary.Q2",&Q2);
 Double_t hdelta;
tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t hbeta;
tsimc->SetBranchAddress("H.hod.beta",&hbeta);
 Double_t emiss;
tsimc->SetBranchAddress("H.kin.secondary.emiss",&emiss);
 Double_t pmiss;
tsimc->SetBranchAddress("H.kin.secondary.pmiss",&pmiss);
 Double_t MM2;
//
    TH1F *hCtime = new TH1F("hCtime",Form("Run %d ; Ctime (ns);Counts",nrun), 800, 20,50);
    HList.Add(hCtime);
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW);
    TH1F *hBeta = new TH1F("hBeta",Form("Run %d ; Beta;Counts",nrun), 300, -.1,1.5);
    HList.Add(hBeta);
    TH2F *hBeta_Ctime = new TH2F("hBeta_Ctime",Form("Run %d ; Beta;Ctime",nrun), 300, -.1,1.5,400,30,50);
    HList.Add(hBeta_Ctime);
    TH1F *hQ2 = new TH1F("hQ2",Form("Run %d ; Q2 (Gev2);Counts",nrun), 300, .0,.5);
    HList.Add(hQ2);
    TH1F *hW_pi0 = new TH1F("hW_pi0",Form("Run %d ; W _pi0 (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW_pi0);
    TH1F *hW_vcs = new TH1F("hW_vcs",Form("Run %d ; W _vcs (Gev);Counts",nrun), 300, .9,1.5);
    HList.Add(hW_vcs);
    TH1F *hMM2 = new TH1F("hMM2",Form("Run %d ; MM2 (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2);
    TH2F *hMM2_W = new TH2F("hMM2_W",Form("Run %d ; MM2 (Gev2);W",nrun), 100, -.05,.05,100,.9,1.5);
    HList.Add(hMM2);
    TH1F *hMM2_ran = new TH1F("hMM2_ran",Form("Run %d ; MM2 random (Gev2);Counts",nrun), 100, -.05,.05);
    HList.Add(hMM2_ran);
    //
    Double_t Ctime_ran_cent=32.89;
    Double_t Ctime_ran_sig=0.5;
    Double_t Ctime_cent=38.85;
    Double_t Ctime_sig=0.4;
    Double_t pi_cent=0.02;
    Double_t pi_sig =0.01;
    Double_t vcs_cent=0.002;
    Double_t vcs_sig =0.004;
 Long64_t nentries = tsimc->GetEntries();
 // nentries=200000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hCtime->Fill(ctime);
		if (abs(hdelta)<8) {
		  hBeta->Fill(hbeta);
		  if (abs(hbeta-.7)<.1) {
		  hBeta_Ctime->Fill(hbeta,ctime);
		if (abs(ctime-Ctime_cent)<2*Ctime_sig) {
		  hW->Fill(W);
		  if (abs(MM2-pi_cent)<pi_sig) hW_pi0->Fill(W);
		  if (abs(MM2-vcs_cent)<vcs_sig) hW_vcs->Fill(W);
		  hQ2->Fill(Q2);
		  MM2=emiss*emiss-pmiss*pmiss;
		  hMM2->Fill(MM2); 
		  hMM2_W->Fill(MM2,W); 
		}
		if (abs(ctime-Ctime_ran_cent)<6*Ctime_ran_sig) {
		  MM2=emiss*emiss-pmiss*pmiss;
		  hMM2_ran->Fill(MM2); 
		}
		  }
		}
	}
//
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
} 
