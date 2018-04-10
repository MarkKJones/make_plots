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
#include "TF1.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCutG.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TCutG.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_shms_sieve(TString basename="",Int_t nrun=2043, Int_t cutflag=0){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 Double_t cer_cut;
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
 Double_t  cer_adc[4];
 if (nrun>=1721) {
   cout << " 3pass " << endl;
    tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
    tsimc->SetBranchAddress("P.ngcer.goodAdcPulseInt",cer_adc);
 }
 if (nrun<1721) {
   cout << " Onepass  nrun= " << nrun << endl;
    tsimc->SetBranchAddress("P.hgcer.npeSum",&sumnpe);
    tsimc->SetBranchAddress("P.hgcer.goodAdcPulseInt",cer_adc);
 }
  Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
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
 Double_t  ztar;
   tsimc->SetBranchAddress("P.react.z",&ztar);
   // Define histograms
   TH1F *hcer_adc[4];
   for (Int_t n=0;n<4;n++) {
     hcer_adc[n] = new TH1F(Form("hcer_adc_%d",n), Form("Run %d ; Cer_%d Adc ; Counts",nrun,n+1),500,0,100);
   HList.Add(hcer_adc[n]);
   }
   TH2F *hetot_npe = new TH2F("hetot_npe", Form("Run %d ; Etrack norm ; Cer Sum",nrun), 200,0,2.,120,0,30);
   HList.Add(hetot_npe);
   TH2F *hetot_delta = new TH2F("hetot_delta", Form("Run %d ; Etrack norm ; delta",nrun), 200,0,2.,120,-20,30);
   HList.Add(hetot_npe);
   TH2F *hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Ytar ; Yptar",nrun), 200,-10.,10,120,-.05,.05);
   HList.Add(hytar_yptar);
   TH2F *hztar_yptar = new TH2F("hztar_yptar", Form("Run %d ; Ztar ; Yptar",nrun), 150,-15.,15,120,-.05,.05);
   HList.Add(hztar_yptar);
   TH1F *hztar = new TH1F("hztar", Form("Run %d ; Ztar ; Counts",nrun), 600,-15.,15.);
   HList.Add(hztar);
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_s ; X_s",nrun), 200,-10,10,200,-20,20);
   HList.Add(hxs_ys);
   TH2F *hxfp_yfp = new TH2F("hxfp_yfp", Form("Run %d ; X_fp ; Y_fp",nrun), 200,-40,40,200,-40,40);
   HList.Add(hxfp_yfp);
   //
     TCutG* ytar_yp_foil[3];
     TH2F *hxfp_yfp_foil[3];
     TH2F *hxs_ys_foil[3];
     Int_t nfoils;
     Double_t ztar_array[3];
     Double_t xs_lo=-1.;
     Double_t xs_hi=+1.;
     Int_t nys=9;
     TH2F *hxfp_yfp_foil1_xscent_yscut[9];
     TH1F *hys_xscent_foil[3];
     Double_t ys_lo[9]={-8.,-6.,-4.5,-2.8,-0.7,0.3,2.1,3.9,5.3};
     Double_t ys_hi[9]={-6.,-4.5,-2.8,-1.0,0.3,2.1,3.9,5.3,7.};
     Double_t theta_cent;
     cer_cut = 0.5;
     if (nrun==1703) {
          cer_cut =-1;
          nfoils = 3;
          ztar_array[0]=-10;
          ztar_array[1]=0;
          ztar_array[2]=10;
          theta_cent=25.*TMath::Pi()/180.;
     }
     if (nrun==1805) {
        nfoils = 2;
          ztar_array[0]=+5;
          ztar_array[1]=-5;
          theta_cent=15.*TMath::Pi()/180.;
     }
     if (nrun==2258) cer_cut =3.;
   if (cutflag==1) {
    TString outputcut="hist/"+basename+"_cut.root";
    TFile *fcut= new TFile(outputcut);
    for (Int_t nc=0;nc<nfoils;nc++) {
      cout << " Looking for cuts " << nc << endl;
      fcut->cd();
      ytar_yp_foil[nc] = (TCutG*)gROOT->FindObject(Form("ytar_yptar_cut_%d",nc));
      Int_t npt =ytar_yp_foil[nc]->GetN();
      cout << " npts = " << npt << endl;
      fsimc->cd();
    }
   }
    for (Int_t nc=0;nc<3;nc++) {
      hxfp_yfp_foil[nc]= new TH2F(Form("hxfp_yfp_foil_%d",nc), Form("Run %d  foil %d; X_fp ; Y_fp",nrun,nc), 200,-40,40,200,-40,40);
      HList.Add(hxfp_yfp_foil[nc]);
      hxs_ys_foil[nc] = new TH2F(Form("hxs_ys_foil_%d",nc), Form("Run %d  foil %d; Y_s ; X_s",nrun,nc), 200,-10,10,200,-20,20);
   HList.Add(hxs_ys_foil[nc]);
      hys_xscent_foil[nc] = new TH1F(Form("hys_xscent_foil_%d",nc), Form("Run %d  Ysieve foil %d xs cent cut; Y_s ; Counts",nrun,nc), 80,-10,10);
   HList.Add(hys_xscent_foil[nc]);
    }   
    for (Int_t nc=0;nc<nys;nc++) {
      hxfp_yfp_foil1_xscent_yscut[nc]= new TH2F(Form("hxfp_yfp_foil1_xscent_yscut_%d",nc), Form("Run %d  foil 1 xs cent %4.2f > ys > %4.2f ; X_fp ; Y_fp",nrun,ys_hi[nc],ys_lo[nc]), 25,-10,15,80,-10,10);
      HList.Add(hxfp_yfp_foil1_xscent_yscut[nc]);
    }

   //
   if (nrun==1796 || nrun==1778 || nrun==2258) {
     nfoils = 1;
   }
   TH2F *hxs_ys_foil1_delta[7];
   Double_t dstart=-7.5;
   Double_t dstep=5.;
   Double_t dcent=dstart;
    for (Int_t nc=0;nc<7;nc++) {
      hxs_ys_foil1_delta[nc] = new TH2F(Form("hxs_ys_foil1_delta_%d",nc), Form("Run %d foil1 %4.2f > delta > %4.2f; Y_s ; X_s",nrun,dcent-dstep/2,dcent+dstep/2), 200,-10,10,200,-15,15);
      dcent+=dstep;
   HList.Add(hxs_ys_foil1_delta[nc]);
    }
   //
// loop over entries
Long64_t nentries = tsimc->GetEntries();
//
 Double_t xs;
 Double_t ys;
	// loop data and apply cuts
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entries = " << i << endl;
		hetot_npe->Fill(etracknorm,sumnpe);
		hetot_delta->Fill(etracknorm,delta);
                   for (Int_t n=0;n<4;n++) {
		     hcer_adc[n]->Fill(cer_adc[n]);
		       }
		if (sumnpe>cer_cut&&etracknorm > .8 && delta>-10 && delta< 24) {
		  xs=xptar*253.;
		  ys=(-0.019*(delta+0.0))+yptar*253.-40.*0.00052*(delta+0.0);
		   hxfp_yfp->Fill(xfp,yfp);
                   hztar->Fill(ztar);
                   hytar_yptar->Fill(ytar,yptar);
                   hztar_yptar->Fill(ztar,yptar);
		   hxs_ys->Fill(ys,xs);
                   Int_t hit_foil=-1;
                   if (cutflag==1) {
	             for (Int_t nc=0;nc<nfoils;nc++) {
		       if (ytar_yp_foil[nc]->IsInside(ytar,yptar)) {
			 hit_foil=nc;
                         ztar=ztar_array[nc];
			 ys = ys-ztar*TMath::Sin(theta_cent)-ztar*TMath::Cos(theta_cent)*yptar;
			 if ( TMath::Abs(xs)< 1) {
                            hxfp_yfp_foil[nc]->Fill(xfp,yfp);
			    hys_xscent_foil[nc]->Fill(ys);
 			 }
		         hxs_ys_foil[nc]->Fill(ys,xs);
		       }
		     }
                           if (TMath::Abs(xs)< 1 && ytar_yp_foil[1]->IsInside(ytar,yptar)) {
	                    for (Int_t ns=0;ns<nys;ns++) {
			      if (ys < ys_hi[ns] && ys >= ys_lo[ns]) hxfp_yfp_foil1_xscent_yscut[ns]->Fill(xfp,yfp);
			    }
			    }
		   } else {
		     hit_foil=1;
		         hxs_ys_foil[1]->Fill(ys,xs);
			 if ( TMath::Abs(xs)< 1) {
                            hxfp_yfp_foil[1]->Fill(xfp,yfp);
			    hys_xscent_foil[1]->Fill(ys);
	                    for (Int_t ns=0;ns<nys;ns++) {
		      if (ys < ys_hi[ns] && ys >= ys_lo[ns]) hxfp_yfp_foil1_xscent_yscut[ns]->Fill(xfp,yfp);
			    }
			 }
		   }
	          
		    if (hit_foil==1) { 
                         Double_t delta_cent=dstart;
			 for (Int_t nd=0;nd<7;nd++) {
			   if(TMath::Abs(delta-delta_cent)<dstep/2.) hxs_ys_foil1_delta[nd]->Fill(ys,xs);
			      delta_cent+=dstep;
			 }
		    }
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
	//
}
