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

void make_hist_shms_ztar_sieve(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_ztar_sieve_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  sumhgnpe;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&sumhgnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
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
   TH1F *hztar[9];
   TH1F *hytar[9];
   TH1F *hyptar[9];
   TH1F *hyptar_cent_foil;
   TH1F *hxptar_cent_foil;
   TH2F *hytar_yptar;
   TH2F *hyfp_yxfp_cent_foil;
   TH2F *hys_xs_cent_foil;
   TH2F *hyfp_yxfp_cent_foil_ypcut[9];
   TH2F *hztar_yptar;
   cout << " nrun = " << nrun << endl;
   Double_t yp_cutlo[9]={-0.0294451,-0.0236552,-0.0171921,-0.0100557,-0.00359254,0.00313991,0.00946841,0.0152583,0.0215868};
   Double_t yp_cuthi[9]={-0.0236552,-0.0171921,-0.00992104,-0.00359254,0.00192807,0.00946841,0.0152583,0.0214522,0.0289925};
   hyptar_cent_foil = new TH1F("hyptar_cent_foil", Form("Run %d Cetner Foil; Yp_tar ; Counts",nrun), 140,-.035,.035);
   HList.Add(hyptar_cent_foil);
   hxptar_cent_foil = new TH1F("hxptar_cent_foil", Form("Run %d Cetner Foil; Xp_tar ; Counts",nrun), 140,-.05,.05);
   HList.Add(hyptar_cent_foil);
   hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Y_tar ; Yp_tar",nrun), 190,-4.,4.,60,-.035,.035);
   HList.Add(hytar_yptar);
   hztar_yptar = new TH2F("hztar_yptar", Form("Run %d ; Yp_tar ; Z_tar",nrun), 60,-.035,.035,120,-15.,15.);
   HList.Add(hztar_yptar);
	  hyfp_yxfp_cent_foil = new TH2F("hyfp_yxfp_cent_foil", Form("Run %d ; Xfp; Yfp",nrun), 160,-40,40,80,-20,20);	  
          HList.Add(hyfp_yxfp_cent_foil);
	  hys_xs_cent_foil = new TH2F("hys_xs_cent_foil", Form("Run %d ; Ys; Xs",nrun), 160,-10,10,80,-15,15);	  
          HList.Add(hys_xs_cent_foil);
	for (Int_t iz = 0; iz < 9; iz++) {
	  hztar[iz] = new TH1F(Form("hztar_%d",iz), Form("Run %d ; Z_tar %5.3f < YP < %5.3f ; Counts",nrun,yp_cutlo[iz],yp_cuthi[iz]), 80,-4.,4.);	  
          HList.Add(hztar[iz]);
	  hytar[iz] = new TH1F(Form("hytar_%d",iz), Form("Run %d ; Y_tar %5.3f < YP < %5.3f ; Counts",nrun,yp_cutlo[iz],yp_cuthi[iz]), 80,-4.,4.);	  
          HList.Add(hytar[iz]);
	  hyptar[iz] = new TH1F(Form("hyptar_%d",iz), Form("Run %d ; Yp_tar %5.3f < YP < %5.3f ; Counts",nrun,yp_cutlo[iz],yp_cuthi[iz]), 80,-.03,.03);	  
          HList.Add(hyptar[iz]);
	  hyfp_yxfp_cent_foil_ypcut[iz] = new TH2F(Form("hyfp_yxfp_cent_foil_ypcut_%d",iz), Form("Run %d  %5.3f < YP < %5.3f; Xfp ; Yfp",nrun,yp_cutlo[iz],yp_cuthi[iz]), 160,-40,40,80,-20,20);	  
          HList.Add(hyfp_yxfp_cent_foil_ypcut[iz]);
 	}
	TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etotnorm ; Counts",nrun),100,0.,2.);
	TH1F *hngsum = new TH1F("hngsum",Form("Run %d ; NG Npe SUM ; Counts",nrun),100,0.,40.);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hetot->Fill(etracknorm);
		hngsum->Fill(sumhgnpe);
		if (etracknorm>.8 && sumhgnpe > 2. && delta>-10 && delta<15) {
		  hztar_yptar->Fill(yptar,reactz);
		  hytar_yptar->Fill(ytar,yptar);
		  if (reactz <5 && reactz>-5) {
                    if (abs(xsieve) < 1) hyptar_cent_foil->Fill(yptar);
                    if (abs(ysieve) < 1) hxptar_cent_foil->Fill(xptar);
		    hyfp_yxfp_cent_foil->Fill(xfp,yfp);
		    hys_xs_cent_foil->Fill(ysieve,xsieve);
	           for (int iz = 0; iz < 9; iz++) {
		     if (yptar>yp_cutlo[iz] && yptar <=yp_cuthi[iz]) {
                         hztar[iz]->Fill(reactz);
                         hytar[iz]->Fill(ytar);
		         hyptar[iz]->Fill(yptar);
                    hyfp_yxfp_cent_foil_ypcut[iz]->Fill(xfp,yfp);
		     }
		   }
		  }
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
