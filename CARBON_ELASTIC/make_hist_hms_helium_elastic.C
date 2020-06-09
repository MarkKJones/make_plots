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

void make_hist_hms_helium_elastic(Int_t nrun=2043,Int_t ntot=100000){
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
   basename=Form("hms_replay_production_mkj_%d_%d",nrun,ntot);
   inputroot="online-ROOTfiles/"+basename+".root";
    cout << "root file = " << inputroot << endl;
  TString outputhist;
   outputhist= "hist/"+basename+"_helium_elastic_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  etotnorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etotnorm);
 Double_t  ys;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ys);
 Double_t  xs;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xs);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  eprime;
   tsimc->SetBranchAddress("H.gtr.p",&eprime);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  W;
   tsimc->SetBranchAddress("H.kin.W",&W);
 Double_t  eth;
   tsimc->SetBranchAddress("H.kin.scat_ang_rad",&eth);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  ntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&ntr);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
    //
      TCutG* xpfp_xfp_cut[11];
      TCutG* ypfp_yfp_cut[11];
      Bool_t cuttest=0;
      Int_t ind_xpfp_xfp=-1;
      Int_t ind_ypfp_yfp=-1;
      Int_t ntot_xcut=0;
      Int_t ntot_ycut=0;
    // Define histograms
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_sieve ; X_sieve",nrun), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys);
	  TH2F *hxpfp_ypfp = new TH2F("hxpfp_ypfp", Form("Run %d ; Yp_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.1,0.1);
          HList.Add(hxpfp_ypfp);
	  TH2F *hxfp_xpfp = new TH2F("hxfp_xpfp", Form("Run %d ; X_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.05,0.05);
          HList.Add(hxpfp_ypfp);
	  TH1F *hdelta = new TH1F("hdelta", Form("Run %d ; Delta ; Counts",nrun), 360,-15,30);
          HList.Add(hdelta);
	  TH1F *hdelta1 = new TH1F("hdelta1", Form("Run %d ; Delta ; Counts",nrun), 25,0,2.5);
          HList.Add(hdelta1);
	  TH2F *hdelta_the = new TH2F("hdelta_the", Form("Run %d ; Delta ; Theta",nrun), 50,-2.5,2.5,100,9,13);
          HList.Add(hdelta_the);
	  TH2F *hex_the = new TH2F("hex_the", Form("Run %d ; Ex ; Theta",nrun), 90,-10,20,100,9,13);
          HList.Add(hex_the);
	  TH2F *hex_delta = new TH2F("hex_delta", Form("Run %d ; Ex ; Delta",nrun), 90,-10,20,50,-2.5,2.5);
          HList.Add(hex_delta);
	  TH2F *hex_xtar = new TH2F("hex_xtar", Form("Run %d ; Ex ; Xtar",nrun), 90,-10,20,25,-1,1);
          HList.Add(hex_xtar);
	  TH2F *hex_reactz = new TH2F("hex_reactz", Form("Run %d ; Ex ; reactz",nrun), 90,-10,20,100,-25,25);
          HList.Add(hex_reactz);
	  TH2F *hdelta_reactz = new TH2F("hdelta_reactz", Form("Run %d ; Delta ; reactz",nrun), 80,-10,10,100,-25,25);
          HList.Add(hdelta_reactz);
	  TH1F *hreactz = new TH1F("hreactz", Form("Run %d ; Z (cm) ; Counts",nrun), 100,-25,25);
          HList.Add(hreactz);
	  TH1F *hreactz_small = new TH1F("hreactz_small", Form("Run %d ; Z (cm) ; Counts",nrun), 100,-25,25);
          HList.Add(hreactz_small);
	  TH1F *hthe = new TH1F("hthe", Form("Run %d ; Scat ang (deg) ; Counts",nrun), 360,9,13);
          HList.Add(hthe);
	  TH1F *hExHe = new TH1F("hExHe", Form("Run %d ; Ex Helium ; Counts",nrun), 600,-10.,20);
          HList.Add(hExHe);
	  TH1F *hExN = new TH1F("hExN", Form("Run %d ; Ex Nitrogen ; Counts",nrun), 60,-10.,20);
          HList.Add(hExHe);
// loop over entries
	  Double_t Eb=2.181*.997;
	  Double_t ExHe=0;
          Double_t MHe=2.80793;
	  Double_t ExN=0;
          Double_t MN=14.0067;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (etotnorm>.6 && ntr>0 ) {
		  Double_t h_entr = 3.536;
		  Double_t v_entr = 9.003;
		  Bool_t small_coll = kTRUE;
		  if (ys > h_entr) small_coll = kFALSE;
		  if (xs > v_entr) small_coll = kFALSE;
		  if (abs(xs) >  (-v_entr/h_entr*abs(ys)+3*v_entr/2))   small_coll = kFALSE;
		  if (small_coll) hreactz_small->Fill(reactz);
		  hreactz->Fill(reactz);
		}
		if (etotnorm>.6 && ntr>0 && abs(reactz)<10.) {
	          Double_t epr=eprime*.997;
		  Double_t A = 0.5;
		  Double_t B = epr-MHe-Eb;
		  Double_t C = epr*Eb*(TMath::Cos(eth)-1)+MHe*(Eb-epr);
                  ExHe=-10000;
                  ExN=-10000;
		  if ((B*B-4*A*C)>0) {
 		    ExHe = (1./2./A)*(-B-TMath::Sqrt(B*B-4*A*C));
		  }
		  B = epr-MN-Eb;
		  C = epr*Eb*(TMath::Cos(eth)-1)+MN*(Eb-epr);
		  if ((B*B-4*A*C)>0) {
 		    ExN = (1./2./A)*(-B-TMath::Sqrt(B*B-4*A*C));
		  }
		hExHe->Fill(ExHe*1000.);
		hExN->Fill(ExN*1000.);
                hxs_ys->Fill(ys,xs);
		hxpfp_ypfp->Fill(ypfp,xpfp);
                hdelta->Fill(delta);
                hdelta1->Fill(delta);
                hdelta_the->Fill(delta,eth*57.3);
                hex_the->Fill(ExHe*1000,eth*57.3);
                hex_delta->Fill(ExHe*1000,delta);
                hex_xtar->Fill(ExHe*1000,xtar);
                hex_reactz->Fill(ExHe*1000,reactz);
                hdelta_reactz->Fill(delta,reactz);
                hthe->Fill(eth*57.3);
		}
	}
	//
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
