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

void make_hist_hms_sieve(TString basename="",Int_t nrun=2043, Int_t nfoils=3){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_hms_sieve_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  beamx;
   tsimc->SetBranchAddress("H.rb.raster.fr_xbpm_tar",&beamx);
 Double_t  beamy;
   tsimc->SetBranchAddress("H.rb.raster.fr_ybpm_tar",&beamy);
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
   TH1F *hztar[3];
   TH1F *hytar[3];
   TH1F *hxtar[3];
   TH1F *hytar_cent[3];
   TH1F *hztar_cent[3];
   TH1F *hys[3];
   TH1F *hxs[3];
   TH1F *hyptar[3];
   TH1F *hxptar[3];
   TH2F *hxptar_yptar[3];
   TH2F *hxs_ys[3];
   cout << " nrun = " << nrun << endl;
   Double_t ztar_cutlo[3]={7.5,-2.5,-12.5};
   Double_t ztar_cuthi[3]={12.5,2.5,-7.5};
   if (nfoils == 2) {
     ztar_cutlo[0]=2.5;
     ztar_cuthi[0]=7.5;
     ztar_cutlo[1]=-7.5;
     ztar_cuthi[1]=-2.5;
     ztar_cutlo[2]=0.;
     ztar_cuthi[2]=0.;
   }
     cout << " make histo " << endl;
   TH2F* hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Y_tar ; Yp_tar",nrun), 80,-6.,6.,60,-.045,.045);
   HList.Add(hytar_yptar);
   TH1F* hetracknorm = new TH1F("hetracknorm", Form("Run %d ; etracknorm ; Counts",nrun), 500,0.,2.0);
   HList.Add(hetracknorm);
   TH1F* hnpesum = new TH1F("hnpesum", Form("Run %d ; npesum ; Counts",nrun), 200,0.,20.0);
   HList.Add(hnpesum);
   TH2F* hetracknorm_npesum = new TH2F("etracknorm_hnpesum", Form("Run %d ; etracknorm; npesum ",nrun), 500,0.,2.0, 200,0.,20.0);
   HList.Add(hnpesum);
   TH1F* hbeamy = new TH1F("hbeamy", Form("Run %d ; Beamy ; Counts",nrun), 500,-.5,.5);
   HList.Add(hbeamy);
   TH1F* hbeamx = new TH1F("hbeamx", Form("Run %d ; Beamx ; Counts",nrun), 500,-.5,.5);
   HList.Add(hbeamx);
   TH1F* hreacty = new TH1F("hreacty", Form("Run %d ; Reacty ; Counts",nrun), 500,-.5,.5);
   HList.Add(hreacty);
   TH1F* hreactx = new TH1F("hreactx", Form("Run %d ; Reactx ; Counts",nrun), 500,-.5,.5);
   HList.Add(hreactx);
   TH1F* hreactz = new TH1F("hreactz", Form("Run %d ; Reactz ; Counts",nrun), 260,-13.,13.);
   HList.Add(hreactz);
    TH2F* hztar_yptar = new TH2F("hztar_yptar", Form("Run %d ; Yp_tar ; Z_tar",nrun), 60,-.045,.045,60,-15.,15.);
   HList.Add(hztar_yptar);
	for (Int_t iz = 0; iz < 3; iz++) {
	  hztar[iz] = new TH1F(Form("hztar_%d",iz), Form("Run %d ; Z_tar ; Counts",nrun), 80,ztar_cutlo[iz],ztar_cuthi[iz]);	  
          HList.Add(hztar[iz]);
	  hytar[iz] = new TH1F(Form("hytar_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; Y_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-5,5);	  
         HList.Add(hytar[iz]);
	  hxtar[iz] = new TH1F(Form("hxtar_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; X_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-1,1);	  
         HList.Add(hxtar[iz]);
	  hytar_cent[iz] = new TH1F(Form("hytar_cent_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f Center Sieve hole; Y_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-5,5);	  
         HList.Add(hytar_cent[iz]);
	  hztar_cent[iz] = new TH1F(Form("hztar_cent_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f Center Sieve hole; Z_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-15,15);	  
         HList.Add(hztar_cent[iz]);
 	  hyptar[iz] = new TH1F(Form("hyptar_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; Yp_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-.04,.04);	  
          HList.Add(hyptar[iz]);
 	  hxptar[iz] = new TH1F(Form("hxptar_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; Xp_tar ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 100,-.1,.1);	  
          HList.Add(hxptar[iz]);
	  hys[iz] = new TH1F(Form("hys_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; Y_s ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 200,-10,10);	  
          HList.Add(hys[iz]);
	  hxs[iz] = new TH1F(Form("hxs_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; X_s ; Counts",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 300,-15,15);	  
          HList.Add(hxs[iz]);
	   hxptar_yptar[iz] = new TH2F(Form("hxptar_yptar_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; Xp_tar ; Yp_tar",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 80,-0.1,.1,60,-.045,.045);
          HList.Add(hxptar_yptar[iz]);
 	   hxs_ys[iz] = new TH2F(Form("hxs_ys_%d",iz), Form("Run %d %5.3f < Ztar < %5.3f; X_sieve ; Y_sieve",nrun,ztar_cutlo[iz],ztar_cuthi[iz]), 150,-15,15.,160,-8,8);
          HList.Add(hxs_ys[iz]);
	}
// loop over entries
	Double_t ysmean=-0.15;
	Double_t yssig = 0.2;
	Double_t xsmean=-0.52;
	Double_t xssig = 5.;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hetracknorm->Fill(etracknorm);
		hnpesum->Fill(sumnpe);
		hetracknorm_npesum->Fill(etracknorm,sumnpe);
		if (sumnpe>2 && etracknorm>.8 && TMath::Abs(delta)<8) {
		  hbeamx->Fill(beamx);
		  hbeamy->Fill(beamy);
		  hreactx->Fill(reactx);
		  hreacty->Fill(reacty);
		  hreactz->Fill(reactz);
		  hztar_yptar->Fill(yptar,reactz);
  		  hytar_yptar->Fill(ytar,yptar);
	        for (int iz = 0; iz < 3; iz++) {
		  if (reactz>ztar_cutlo[iz] && reactz <=ztar_cuthi[iz]) {
		    hztar[iz]->Fill(reactz);
		    hytar[iz]->Fill(ytar);
		    hxtar[iz]->Fill(xtar);
		    if (TMath::Abs(ysieve-ysmean) < 3*yssig && TMath::Abs(xsieve-xsmean) < 3*xssig) hytar_cent[iz]->Fill(ytar);
		    if (TMath::Abs(ysieve-ysmean) < 3*yssig && TMath::Abs(xsieve-xsmean) < 3*xssig) hztar_cent[iz]->Fill(reactz);
		    hyptar[iz]->Fill(yptar);
		    hxptar[iz]->Fill(xptar);
		    hys[iz]->Fill(ysieve);
		    hxs[iz]->Fill(xsieve);
		    hxptar_yptar[iz]->Fill(xptar,yptar);
		    hxs_ys[iz]->Fill(xsieve,ysieve);
		  }
		}
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
