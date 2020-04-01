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
#include <TCutG.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_shms(Int_t nrun=2043,Int_t ntot=1000000){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   TString basename;
   basename=Form("shms_replay_production_default_0%d_%d",nrun,ntot);
   inputroot="polhe_rootfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_polhe_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  etotnorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etotnorm);
 Double_t  ys;
   tsimc->SetBranchAddress("P.extcor.ysieve",&ys);
 Double_t  xs;
   tsimc->SetBranchAddress("P.extcor.xsieve",&xs);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("P.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("P.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("P.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("P.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  eprime;
   tsimc->SetBranchAddress("P.gtr.p",&eprime);
 Double_t  eth;
   tsimc->SetBranchAddress("P.kin.scat_ang_rad",&eth);
 Double_t  yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
 Double_t  ntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntr);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
    //
   // Define histograms
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_sieve ; X_sieve",nrun), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys);
	  TH2F *hxpfp_ypfp = new TH2F("hxpfp_ypfp", Form("Run %d ; Yp_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.05,0.05);
          HList.Add(hxpfp_ypfp);
	  TH2F *hxfp_xpfp = new TH2F("hxfp_xpfp", Form("Run %d ; X_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.05,0.05);
          HList.Add(hxpfp_ypfp);
	  TH1F *hdelta = new TH1F("hdelta", Form("Run %d ; Delta ; Counts",nrun), 360,-15,30);
          HList.Add(hdelta);
	  TH1F *hdelta1 = new TH1F("hdelta1", Form("Run %d ; Delta ; Counts",nrun), 50,5,20);
          HList.Add(hdelta1);
	  TH2F *hdelta_the = new TH2F("hdelta_the", Form("Run %d ; Delta ; Theta",nrun), 50,5,20,80,6,10);
          HList.Add(hdelta_the);
	  TH2F *hdelta_reactz = new TH2F("hdelta_reactz", Form("Run %d ; Delta ; Reactz",nrun), 100,5,20,100,-25,25);
          HList.Add(hdelta_reactz);
	  TH2F *hex_the = new TH2F("hex_the", Form("Run %d ; Ex ; Theta",nrun), 500,-10,20,80,6,10);
          HList.Add(hex_the);
	  TH1F *hthe = new TH1F("hthe", Form("Run %d ; Scat ang (deg) ; Counts",nrun), 80,6,10);
          HList.Add(hthe);
	  TH1F *hEx = new TH1F("hEx", Form("Run %d ; Ex ; Counts",nrun), 90,-10.,20);
          HList.Add(hEx);
// loop over entries
	  Double_t Eb=2.18151;
	  Double_t Ex=0;
          Double_t Mc=12*.9315;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (etotnorm>.6 && ntr>0 && abs(delta+5)<1) {
	Double_t epr=eprime;
		  Double_t A = 0.5;
		  Double_t B = epr-Mc-Eb;
		  Double_t C = epr*Eb*(TMath::Cos(eth)-1)+Mc*(Eb-epr);
                  Ex=-10000;
		  if ((B*B-4*A*C)>0) {
 		    Ex = (1./2./A)*(-B-TMath::Sqrt(B*B-4*A*C));
		  }
		  hEx->Fill(Ex*1000.);
                hxs_ys->Fill(ys,xs);
		hxpfp_ypfp->Fill(ypfp,xpfp);
                hdelta->Fill(delta);
                hdelta1->Fill(delta);
                hdelta_the->Fill(delta,eth*57.3);
                hdelta_reactz->Fill(delta,reactz);
                hex_the->Fill(Ex*1000,eth*57.3);
                hthe->Fill(eth*57.3);
		}
	}
	//
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
