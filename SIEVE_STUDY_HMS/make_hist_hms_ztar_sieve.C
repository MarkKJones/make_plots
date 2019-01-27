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

void make_hist_hms_ztar_sieve(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hms_ztar_sieve_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);

   // Define histograms
   
   TH1F *hztar[3][7];
   TH1F *hyptar[3][7];
   TH1F *hyptar_foil[3];
   TH1F *hys_foil[3];
   TH1F *hxs_foil[3];
   TH1F *hys_xscent_foil[3];
   TH1F *hys2_xs2cent_foil[3];
   TH2F *hytar_yptar;
   TH2F *hyfp_yxfp_foil[3];
   TH2F *hys_xs_foil[3];
   TH2F *hys2_xs2_foil[3];
   TH2F *hyfp_yxfp_foil_ypcut[3][7];
   TH2F *hztar_yptar;
   cout << " nrun = " << nrun << endl;
   Double_t yp_cutlo[3][7]={-0.0297144,-0.0213662,-0.0128833,-0.00359254,0.00623683,0.0155276,0.0250877
  ,-0.0314526,-0.0241815,-0.0144868,-0.00513834,0.00542198,0.0135586,0.0239458
,-0.0297144,-0.0213662,-0.0128833,-0.00359254,0.00623683,0.0155276,0.0250877};
   Double_t yp_cuthi[3][7]={-0.0240592,-0.0154416,-0.00709341,0.00354385,0.0133732,0.0223947,0.029935
     ,-0.0245278,-0.0155255,-0.00617706,0.00438326,0.0133855,0.0232533,0.03139
     ,-0.0240592,-0.0154416,-0.00709341,0.00354385,0.0133732,0.0223947,0.029935};
 	for (Int_t nf = 0; nf < 3; nf++) {
          hyptar_foil[nf] = new TH1F(Form("hyptar_foil_%d",nf), Form("Run %d Foil %d; Yp_tar ; Counts",nrun,nf), 140,-.045,.045);
          HList.Add(hyptar_foil[nf]);
          hys_foil[nf] = new TH1F(Form("hys_foil_%d",nf), Form("Run %d Foil %d; Ys ; Counts",nrun,nf), 140,-10.,10.);
          HList.Add(hys_foil[nf]);
          hxs_foil[nf] = new TH1F(Form("hxs_foil_%d",nf), Form("Run %d Foil %d; Xs ; Counts",nrun,nf), 140,-15.,15.);
          HList.Add(hxs_foil[nf]);
          hys_xscent_foil[nf] = new TH1F(Form("hys_xscent_foil_%d",nf), Form("Run %d XS cent Foil %d; Ys ; Counts",nrun,nf), 140,-10.,10.);
          HList.Add(hys_xscent_foil[nf]);
          hys2_xs2cent_foil[nf] = new TH1F(Form("hys2_xs2cent_foil_%d",nf), Form("Run %d XS_new cent Foil %d; Ys_new ; Counts",nrun,nf), 140,-10.,10.);
          HList.Add(hys2_xs2cent_foil[nf]);
          hys_xs_foil[nf] = new TH2F(Form("hys_xs_foil_%d",nf), Form("Run %d Foil %d; Ys ; Xs",nrun,nf), 70,-10.,10., 70,-15.,15.);
          HList.Add(hys_xs_foil[nf]);
          hys2_xs2_foil[nf] = new TH2F(Form("hys2_xs2_foil_%d",nf), Form("Run %d Foil %d; Ys_new ; Xs_new",nrun,nf), 70,-10.,10., 70,-15.,15.);
          HList.Add(hys2_xs2_foil[nf]);
	  hyfp_yxfp_foil[nf] = new TH2F(Form("hyfp_yxfp_foil_%d",nf), Form("Run %d Foil %d; Xfp; Yfp",nrun,nf), 160,-40,40,80,-20,20);	  
          HList.Add(hyfp_yxfp_foil[nf]);
	}
   hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Y_tar ; Yp_tar",nrun), 190,-4.,4.,60,-.045,.045);
   HList.Add(hytar_yptar);
   hztar_yptar = new TH2F("hztar_yptar", Form("Run %d ; Yp_tar ; Z_tar",nrun), 60,-.045,.045,120,-15.,15.);
   HList.Add(hztar_yptar);
	for (Int_t iz = 0; iz < 7; iz++) {
 	for (Int_t nf = 0; nf < 3; nf++) {
	  hztar[nf][iz] = new TH1F(Form("hztar_%d_%d",nf,iz), Form("Run %d foil %d ; Z_tar %5.3f < YP < %5.3f ; Counts",nrun,nf,yp_cutlo[nf][iz],yp_cuthi[nf][iz]), 150,-15.,15.);	  
          HList.Add(hztar[nf][iz]);
	  hyptar[nf][iz] = new TH1F(Form("hyptar_%d_%d",nf,iz), Form("Run %d foil %d ; Yp_tar %5.3f < YP < %5.3f ; Counts",nrun,nf,yp_cutlo[nf][iz],yp_cuthi[nf][iz]), 80,-.045,.045);	  
          HList.Add(hyptar[nf][iz]);
	  hyfp_yxfp_foil_ypcut[nf][iz] = new TH2F(Form("hyfp_yxfp_foil_%d_ypcut_%d",nf,iz), Form("Run %d Foil %d %5.3f < YP < %5.3f; Xfp ; Yfp",nrun,nf,yp_cutlo[nf][iz],yp_cuthi[nf][iz]), 160,-40,40,80,-20,20);	  
          HList.Add(hyfp_yxfp_foil_ypcut[nf][iz]);
	}
 	}
	TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etotnorm ; Counts",nrun),100,0.,2.);
	TH1F *hetot_good = new TH1F("hetot_good",Form("Run %d ; Etotnorm (good) ; Counts",nrun),100,0.,2.);
	TH1F *hngsum = new TH1F("hngsum",Form("Run %d ; Npe SUM ; Counts",nrun),100,0.,40.);
// loop over entries
Long64_t nentries = tsimc->GetEntries();
 Double_t zlo[3]={-15.,-2.5,5};
 Double_t zhi[3]={-5., 2.5,15.};
 //nentries=200000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hetot->Fill(etracknorm);
		hngsum->Fill(sumnpe);
				if (etracknorm>.7 && sumnpe > 5.&& delta>-8 && delta<8) {
		//		if (delta>-8 && delta<8) {
				  //		hztar_yptar->Fill(yptar,reactz);
		hetot_good->Fill(etracknorm);
		hytar_yptar->Fill(ytar,yptar);
				 		hztar_yptar->Fill(yptar,reactz);
		Double_t slope_lo = 239.281 ;
		Double_t slope_hi = 238.095;
		Double_t int_lo = -3.01667;
		  Double_t int_hi = -0.11;
		  //Int_t nf=1;
		  	for (Int_t nf = 0; nf < 3; nf++) {
		  		if (reactz <= zhi[nf] && reactz> zlo[nf]) {
		  //if (TMath::Abs(reactz-(slope_hi*yptar+int_hi))<2.5  ) {
		     Double_t ys_calc= yptar*168.;
		     Double_t xs_calc= xptar*168.;
                     hyptar_foil[nf]->Fill(yptar);
                     hys_foil[nf]->Fill(ysieve);
                     if (xsieve>-1 && xsieve < 1.5) hys_xscent_foil[nf]->Fill(ysieve);
                     if (xsieve>-1 && xsieve < 1.5) hys2_xs2cent_foil[nf]->Fill(ys_calc);
                     hxs_foil[nf]->Fill(xsieve);
                     hys_xs_foil[nf]->Fill(ysieve,xsieve);
                     hys2_xs2_foil[nf]->Fill(ys_calc,xs_calc);
		    hyfp_yxfp_foil[nf]->Fill(xfp,yfp);
	        for (int iz = 0; iz < 7; iz++) {
		  if (yptar>yp_cutlo[nf][iz] && yptar <=yp_cuthi[nf][iz]) {
		    hztar[nf][iz]->Fill(reactz);
		    hyptar[nf][iz]->Fill(yptar);
                    hyfp_yxfp_foil_ypcut[nf][iz]->Fill(xfp,yfp);
		    }
		}
		  }
			}
				}
	}
	//
	cout << " Integral = " << hetot_good->Integral() << endl;
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
