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

void make_hist_shms_xs_ys(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_xs_ys_hist.root";
   TObjArray HList(0);
   //
  TString outCutFile;
  outCutFile=Form("cuts/ytar_delta_%d_cut.root", nrun);
    TFile fcut(outCutFile);
       TCutG* delta_ytar_cut[3];
   for (Int_t nc=0;nc<3;nc++) {
      delta_ytar_cut[nc] = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
      if (delta_ytar_cut[nc]) {
      Int_t npt = delta_ytar_cut[nc]->GetN();
      cout << "delta_ytar_cut = " << nc << " npts = " << npt << endl;
      }
   }
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  sumhgnpe;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
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
 Double_t  yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
 Double_t  ysieve;
   tsimc->SetBranchAddress("P.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("P.extcor.xsieve",&xsieve);

   // Define histograms
	TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etotnorm ; Counts",nrun),100,0.,2.);
        HList.Add(hetot);
	TH1F *hngsum = new TH1F("hngsum",Form("Run %d ; NG Npe SUM ; Counts",nrun),100,0.,40.);
	HList.Add(hngsum);
	TH2F *hYtarDelta = new TH2F("hYtarDelta",Form("Run %d ; Ytar ; Delta",nrun),100,-10.,10.,100,-15.,25.);
	HList.Add(hYtarDelta);
	TH2F *hZtarDelta = new TH2F("hZtarDelta",Form("Run %d ; Ztar ; Delta",nrun),100,-15.,15.,100,-15.,25.);
	HList.Add(hZtarDelta);
	TH2F *hXsYs[3];
	TH1F *hYptar[3];
	static const Int_t ndtot=5;
	TH2F *hXsYs_delta[3][ndtot];
	//TH2F *hXsYs_calc_delta[3][ndtot];
	Double_t delta_low[ndtot]={-15.,-10,-5.,0,5};
	Double_t delta_hi[ndtot]={-10.,-5,0.,5,24};
        for (Int_t nc=0;nc<3;nc++) {
	  hXsYs[nc] = new TH2F(Form("hXsYs_foil_%d",nc),Form("Run %d ; Ys ; Xs",nrun),90,-9,9,60,-15,15);
          HList.Add(hXsYs[nc]);
	  hYptar[nc] = new TH1F(Form("hYptar_foil_%d",nc),Form("Run %d Foil %d; Yptar ",nrun,nc),100,-.06,.06);
          HList.Add(hYptar[nc]);
           for (Int_t nd=0;nd<ndtot;nd++) {
	     hXsYs_delta[nc][nd] = new TH2F(Form("hXsYs_foil_%d_delta_%d",nc,nd),Form("Delta reg %d Run %d ; Ys ; Xs",nd,nrun),90,-9,9,90,-18,18);
	     HList.Add(hXsYs_delta[nc][nd]);
	     //	     hXsYs_calc_delta[nc][nd] = new TH2F(Form("hXsYs_calc_foil_%d_delta_%d",nc,nd),Form("Delta reg %d Run %d ; Ys ; Xs",nd,nrun),90,-9,9,90,-18,18);
	     //     HList.Add(hXsYs_calc_delta[nc][nd]);
	   }
	}
// loop over entries
	Double_t ztar[3]={-10,0,10};
	Double_t ytar_cent[3];
	Double_t ang=30.065;
		      Double_t ys_calc= 0.; 
		      Double_t xs_calc= 0.; 
                  for (Int_t nc=0;nc<3;nc++) {
		    ytar_cent[nc] = -ztar[nc]*TMath::Sin(ang*3.14159/180);
		  }	
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (sumnpe > 6.) hetot->Fill(etracknorm);
		if (etracknorm>.8) hngsum->Fill(sumnpe);
		if (etracknorm>.8 && sumnpe > 6. && delta>-15 && delta<24) {
		  hYtarDelta->Fill(ytar,delta);
		  hZtarDelta->Fill(reactz,delta);
                  for (Int_t nc=0;nc<3;nc++) {
		    if (delta_ytar_cut[nc] && delta_ytar_cut[nc]->IsInside(ytar,delta))  {
	 		      ys_calc= ytar_cent[nc] + yptar*240.; 
		      xs_calc= xptar*240.; 
		      hXsYs[nc]->Fill(ysieve,xsieve);
		      hYptar[nc]->Fill(yptar);
                      for (Int_t nd=0;nd<ndtot;nd++) {
			if (delta>=delta_low[nd] && delta<delta_hi[nd]) {
                               hXsYs_delta[nc][nd]->Fill(ysieve,xsieve);
			       //     hXsYs_calc_delta[nc][nd]->Fill(ys_calc,xs_calc);
			}}
		    }
                  }
		  
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
