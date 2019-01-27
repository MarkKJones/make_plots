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

void make_hist_shms_ep_elastic(TString basename="",Int_t nrun=3288,Double_t pfac=0.981){
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
   outputhist= "hist/"+basename+"_shms_ep_elastic_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
 Double_t  h_ytar;
   tsimc->SetBranchAddress("H.gtr.y",&h_ytar);
 Double_t  gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex);
 Double_t  h_bpmtarx;
   tsimc->SetBranchAddress("H.rb.raster.fr_xbpm_tar",&h_bpmtarx);
 Double_t  h_bpmtary;
   tsimc->SetBranchAddress("H.rb.raster.fr_ybpm_tar",&h_bpmtary);
 Double_t  e_bpmtarx;
   tsimc->SetBranchAddress("P.rb.raster.fr_xbpm_tar",&e_bpmtarx);
 Double_t  e_bpmtary;
   tsimc->SetBranchAddress("P.rb.raster.fr_ybpm_tar",&e_bpmtary);
 Double_t  e_rasterx;
   tsimc->SetBranchAddress("P.rb.raster.fr_xa",&e_rasterx);
 Double_t  e_rastery;
   tsimc->SetBranchAddress("P.rb.raster.fr_ya",&e_rastery);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("P.react.z",&e_reactz);
 Double_t  e_beamy;
   tsimc->SetBranchAddress("P.react.y",&e_beamy);
 Double_t  e_beamx;
   tsimc->SetBranchAddress("P.react.x",&e_beamx);
 Double_t  h_beamy;
   tsimc->SetBranchAddress("H.react.y",&h_beamy);
 Double_t  h_beamx;
   tsimc->SetBranchAddress("H.react.x",&h_beamx);
 Double_t  h_reactz;
   tsimc->SetBranchAddress("H.react.z",&h_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&e_delta);
 Double_t  e_mom;
   tsimc->SetBranchAddress("P.gtr.p",&e_mom);
 Double_t  ptrig6;
   tsimc->SetBranchAddress("T.coin.pTRIG6_TdcTime",&ptrig6);
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
   tsimc->SetBranchAddress("P.kin.primary.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.primary.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 300, 0.9,1.5);
    HList.Add(hW);
    TH1F *h_ereactz = new TH1F("h_ereactz",Form("Run %d ; SHMS react_z (cm) ;Counts",nrun), 300, -15.,15.);
    HList.Add(h_ereactz);
    TH1F *h_ereactx = new TH1F("h_ereactx",Form("Run %d ; SHMS react_x (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(h_ereactx);
    TH1F *h_eytar = new TH1F("h_eytar",Form("Run %d ; SHMS Ytar (cm);Counts",nrun), 300, -5,5.);
    HList.Add(h_eytar);
    TH1F *h_hytar = new TH1F("h_hytar",Form("Run %d ; HMS Ytar (cm);Counts",nrun), 300, -5,5.);
    HList.Add(h_hytar);
    TH1F *h_ereacty = new TH1F("h_ereacty",Form("Run %d ; SHMS react_y (cm) ;Counts",nrun), 300, -.5,.5);
    HList.Add(h_ereacty);
    TH1F *h_eEpicsx = new TH1F("h_eEpicsx",Form("Run %d ; SHMS EPICS X BPM at target (cm) ;Counts",nrun), 100, -.3,.3);
    HList.Add(h_eEpicsx);
    TH1F *h_eEpicsy = new TH1F("h_eEpicsy",Form("Run %d ; SHMS EPICS Y BPM at target (cm)  ;Counts",nrun), 100, -.3,.3);
    HList.Add(h_eEpicsy);
    TH1F *h_eRasterx = new TH1F("h_eRasterx",Form("Run %d ; SHMS Raster X (cm) ;Counts  ",nrun), 300, -.5,.5);
    HList.Add(h_eRasterx);
    TH1F *h_eRastery = new TH1F("h_eRastery",Form("Run %d ; SHMS Raster Y (cm) ;Counts  ",nrun), 300, -.5,.5);
    HList.Add(h_eRastery);
    TH1F *h_hreactz = new TH1F("h_hreactz",Form("Run %d ; HMS react_z (cm) ;Counts",nrun), 300, -15.,15.);
    HList.Add(h_hreactz);
    TH1F *h_hreactx = new TH1F("h_hreactx",Form("Run %d ; HMS react_x (cm) ;Counts",nrun), 300, -.5,.5);
    HList.Add(h_hreactx);
    TH1F *h_hreacty = new TH1F("h_hreacty",Form("Run %d ; HMS react_y(cm)  ;Counts",nrun), 300, -.5,.5);
    HList.Add(h_hreacty);
    TH1F *h_diffreactz = new TH1F("h_diffreactz",Form("Run %d ; HMS react_z -  SHMS react_z (cm) ;Counts",nrun), 300, -5.,5.);
    HList.Add(h_diffreactz);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; Xptar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
     TH1F *hetracknorm = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetracknorm);
     TH2F *hWXfp = new TH2F("hWXfp",Form("Run %d ; W (GeV) ; SHMS  Xfp ",nrun), 100, 0.8,1.2, 80,-40,40);
    HList.Add(hWXfp);
        TH2F *hWYfp = new TH2F("hWYfp",Form("Run %d ; W (GeV) ; SHMS  Yfp ",nrun), 100, 0.8,1.2, 80,-40,40);
   HList.Add(hWYfp);
    TH2F *hWXpfp = new TH2F("hWXpfp",Form("Run %d ; W (GeV) ; SHMS  Xpfp ",nrun), 100, 0.8,1.2, 80,-.05,.05);
   HList.Add(hWXpfp);
        TH2F *hWYpfp = new TH2F("hWYpfp",Form("Run %d ; W (GeV) ; SHMS Ypfp ",nrun), 100, 0.8,1.2, 80,-.02,.02);
   HList.Add(hWYpfp);
  // loop over entries
    Double_t th_cent=7.5;
  Double_t Mp = .93827;
   Double_t Ei=2.221;
   Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (gindex>-1 ) {
		  /*		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
		hW->Fill(W);
		Double_t Enew=p_spec*(1+e_delta/100.);
		Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		hWcalc->Fill(W_calc);*/
		  if (TMath::Abs(W-.910)<.05) hxptar->Fill(e_xptar);
		  hetracknorm->Fill(etracknorm);
		  if (etracknorm>1) hW->Fill(W);
		if (W<1.075 && etracknorm>1) {
		  h_eytar->Fill(e_ytar);
		  h_hytar->Fill(h_ytar);
		  h_ereactz->Fill(e_reactz);
		  h_ereactx->Fill(e_beamx);
		  h_ereacty->Fill(e_beamy);
		  h_eEpicsx->Fill(e_bpmtarx);
		  h_eEpicsy->Fill(e_bpmtary);
		  h_eRasterx->Fill(e_rasterx);
		  h_eRastery->Fill(e_rastery);
		  h_hreactz->Fill(h_reactz);
		  h_hreactx->Fill(h_beamx);
		  h_hreacty->Fill(h_beamy);
		  h_diffreactz->Fill(h_reactz-e_reactz);
		hWXfp->Fill(W,e_xfp);
		hWYfp->Fill(W,e_yfp);
		hWXpfp->Fill(W,e_xpfp);
		hWYpfp->Fill(W,e_ypfp);
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
