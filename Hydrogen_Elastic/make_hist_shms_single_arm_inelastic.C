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

void make_hist_shms_single_arm_inelastic(TString basename="",Int_t nrun=3288,Double_t pfac=1.0){
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
   outputhist= "hist/"+basename+"_shms_single_arm_inelastic_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  beamx;
  tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&beamx);
 Double_t  beamy;
   tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&beamy);
 Double_t  rastx;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&rastx);
 Double_t  rasty;
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&rasty);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
   Double_t  npeSum;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&npeSum);
 Double_t  DeltaDp;
   tsimc->SetBranchAddress("P.extcor.delta_dp",&DeltaDp);
 Double_t  e_xtar;
   tsimc->SetBranchAddress("P.gtr.x",&e_xtar);
 Double_t  gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etracknorm);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("P.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&e_delta);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&e_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&e_xpfp);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 2.,3.0);
    HList.Add(hW);
    TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; Delta ;Counts",nrun), 100, -15.,25.0);
    HList.Add(hdelta);
    TH2F *hWdelta = new TH2F("hWdelta",Form("Run %d ; Delta  ; W",nrun), 50, -15.,25.0, 50, .8,1.2);
    HList.Add(hWdelta);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.6);
    HList.Add(hetot);
    TH2F *hetot_npe = new TH2F("hetot_npe",Form("Run %d ; Etrack norm; NG npe",nrun), 120, 0.0,1.6,80,0,40);
    HList.Add(hetot_npe);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
      TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
  // loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
 		if (gindex>-1 )  {
                   hetot->Fill(etracknorm);
                   hetot_npe->Fill(etracknorm,npeSum);
		}
		if (gindex>-1 && e_delta > -10. && e_delta < 20. && npeSum>3. && etracknorm>.8) { 
		hW->Fill(W);
		  if (1==1) {
		  hxfp->Fill(e_xfp);		  
		  hyfp->Fill(e_yfp);		  
		  hxpfp->Fill(e_xpfp);		  
		  hypfp->Fill(e_ypfp);		  
		  hxptar->Fill(e_xptar);		  
		  hyptar->Fill(e_yptar);		  
		  hdelta->Fill(e_delta);
		  }
		//
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
