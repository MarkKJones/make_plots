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

void make_hist_shms_single_arm(TString basename="",Int_t nrun=3288,Double_t pfac=1.0){
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
   outputhist= "hist/"+basename+"_shms_single_arm_hist.root";
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
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 100, 0.8,1.2);
    HList.Add(hW);
    TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; Delta ;Counts",nrun), 100, -15.,25.0);
    HList.Add(hdelta);
    TH2F *hWdelta = new TH2F("hWdelta",Form("Run %d ; Delta  ; W",nrun), 50, -15.,25.0, 50, .8,1.2);
    HList.Add(hWdelta);
    TH1F *hDeltaDp = new TH1F("hDeltaDp",Form("Run %d ; DeltaDP;Counts",nrun), 300, -1.,1.);
     HList.Add(hDeltaDp);
    TH1F *hxtar = new TH1F("hxtar",Form("Run %d ; Xtar (cm);Counts",nrun), 300, -.3,.3);
    HList.Add(hxtar);
    TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar (cm);Counts",nrun), 300, -10.,10.);
    HList.Add(hztar);
    Int_t nzcut=3;
    Double_t zcut_cent[3]={-4.,0,4.};
    TH2F *hW_xtar[3];
     TH2F *hW_xptar[3];
     for (Int_t nz=0;nz<3;nz++) {
       hW_xtar[nz] = new TH2F(Form("hW_xtar_%d",nz),Form("Run %d Zcut_cent=%3.1f;  W (GeV);Xtar (cm)",nrun,zcut_cent[nz]), 100, 0.90,1.0,100, -.3,.3);
    HList.Add(hW_xtar[nz]);
    hW_xptar[nz] = new TH2F(Form("hW_xptar_%d",nz),Form("Run %d Zcut_cent=%f5.2;  W (GeV);Xptar (rad)",nrun,zcut_cent[nz]), 100, 0.90,1.0,100, -.1,.1);
    HList.Add(hW_xptar[nz]);
 }
   TH2F *hW_DeltaDp = new TH2F("hW_DeltaDp",Form("Run %d ; W (GeV) ; DeltaDP ",nrun), 100, 0.90,1.0 ,100, -.3,.3);
    HList.Add(hW_DeltaDp);
    TH2F *hxtar_DeltaDp = new TH2F("hxtar_DeltaDp",Form("Run %d ; Xtar (cm) ; DeltaDP ",nrun), 100, -0.3,.3 ,100, -.3,.3);
    HList.Add(hxtar_DeltaDp);
    TH2F *hztar_xtar = new TH2F("hztar_xtar",Form("Run %d ;  Ztar (cm) ;Xtar (cm)",nrun), 300, -10.,10. ,100, -.3,.3);
    HList.Add(hztar_xtar);
    TH2F *hztar_yptar = new TH2F("hztar_yptar",Form("Run %d ;  Ztar (cm) ;Yptar (rad)",nrun), 300, -10.,10. ,100, -.04,.04);
    HList.Add(hztar_yptar);
    TH2F *hztar_xptar = new TH2F("hztar_xptar",Form("Run %d ;  Ztar (cm) ;Xptar (rad)",nrun), 300, -10.,10. ,100, -.1,.1);
    HList.Add(hztar_xptar);
    TH2F *hxptar_xtar = new TH2F("hxptar_xtar",Form("Run %d ;  Xptar (rad) ;Xtar (cm)",nrun), 300, -.1,.1 ,100, -.3,.3);
    HList.Add(hxptar_xtar);
    TH2F *hztar_xbeam = new TH2F("hztar_xbeam",Form("Run %d ;  Ztar (cm) ;Xbeam (cm)",nrun), 300, -10.,10. ,100, -.3,.3);
    HList.Add(hztar_xbeam);
     TH1F *hxrast = new TH1F("hxrast",Form("Run %d ; Xrast (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hxrast);
   TH1F *hyrast = new TH1F("hyrast",Form("Run %d ; Yrast (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hyrast);
    TH1F *hxbeam = new TH1F("hxbeam",Form("Run %d ; Xbeam (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hxbeam);
    TH1F *hybeam = new TH1F("hybeam",Form("Run %d ; Ybeam (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(hybeam);
    TH1F *hWcalc = new TH1F("hWcalc",Form("Run %d Fac = %5.3f; W (GeV);Counts",nrun,pfac),  300, 0.9,1.5);
    HList.Add(hWcalc);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.6);
    HList.Add(hetot);
    TH1F *hDeltaDiff = new TH1F("hDeltaDiff",Form("Run %d  Fac = %5.3f; DeltaDiff ;Counts",nrun,pfac), 300, -5,5);
    HList.Add(hDeltaDiff);
    TH2F *hDeltaDiffXfp = new TH2F("hDeltaDiffXfp",Form("Run %d  Fac = %5.3f; DeltaDiff  ; Xfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffYfp = new TH2F("hDeltaDiffYfp",Form("Run %d  Fac = %5.3f; DeltaDiff  ; Yfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffXpfp = new TH2F("hDeltaDiffXpfp",Form("Run %d  Fac = %5.3f; DeltaDiff  ; Xpfp ",nrun,pfac), 300, -1,2, 80,-.05,.05);
    TH2F *hDeltaDiffYpfp = new TH2F("hDeltaDiffYpfp",Form("Run %d  Fac = %5.3f; DeltaDiff  ; Ypfp ",nrun,pfac), 300, -1,2, 80,-.02,.02);
    HList.Add(hDeltaDiffXfp);
    HList.Add(hDeltaDiffXpfp);
    HList.Add(hDeltaDiffYfp);
    HList.Add(hDeltaDiffYpfp);
    TH2F *hWXfp = new TH2F("hWXfp",Form("Run %d ; W (GeV) ; Xfp ",nrun), 100, 0.8,1.2, 80,-40,40);
    TH2F *hWYfp = new TH2F("hWYfp",Form("Run %d ; W (GeV) ; Yfp ",nrun), 100, 0.8,1.2, 80,-40,40);
    TH2F *hWXpfp = new TH2F("hWXpfp",Form("Run %d ; W (GeV) ; Xpfp ",nrun), 100, 0.8,1.2, 80,-.05,.05);
    TH2F *hWYpfp = new TH2F("hWYpfp",Form("Run %d ; W (GeV) ; Ypfp ",nrun), 100, 0.8,1.2, 80,-.02,.02);
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
    Double_t th_cent=19.995;
  Double_t Mp = .93827;
   Double_t Ei=10.587;
    Double_t p_spec=6.3;
  if (nrun >=4784 && nrun >=4786) {
       th_cent = 19.995;
       p_spec = 6.3;
   }
    if (nrun ==4789 ) {
       th_cent = 18.;
       p_spec = 6.842;
   }
    if (nrun ==4791 ) {
       th_cent = 16.25;
       p_spec = 6.842;
   }
    if (nrun >=4793 && nrun >=4798 ) {
       th_cent = 19.85;
       p_spec = 6.842;
   }
    Double_t scalfac=1.0;
   Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                hxrast->Fill(rastx);
                hyrast->Fill(rasty);
                hxbeam->Fill(beamx);
                hybeam->Fill(beamy);
 		if (gindex>-1 )  hetot->Fill(etracknorm);
		if (gindex>-1 && e_delta > -10. && e_delta < 30. && npeSum>2. && etracknorm>.5) { 
		hW->Fill(W);
		  if ( W> 0.850 &&W < 1.075) {
		  hxfp->Fill(e_xfp,scalfac);		  
		  hyfp->Fill(e_yfp,scalfac);		  
		  hxpfp->Fill(e_xpfp,scalfac);		  
		  hypfp->Fill(e_ypfp,scalfac);		  
		  hxptar->Fill(e_xptar,scalfac);		  
		  hyptar->Fill(e_yptar,scalfac);		  
		  hdelta->Fill(e_delta,scalfac);
               hxtar->Fill(e_xtar);
                hDeltaDp->Fill(DeltaDp);
                hztar->Fill(e_reactz);
		      for (Int_t nz=0;nz<3;nz++) {
			if (TMath::Abs(e_reactz-zcut_cent[nz]) < .5 ) {
			  if ( TMath::Abs(e_xptar) < .005) hW_xtar[nz]->Fill(W,e_xtar);
                           if ( TMath::Abs(rasty) < .05)   hW_xptar[nz]->Fill(W,e_xptar);
			}
		      }
                        hW_DeltaDp->Fill(W,DeltaDp);
                        hxtar_DeltaDp->Fill(e_xtar,DeltaDp);
                        hxptar_xtar->Fill(e_xptar,e_xtar);
                hztar_xtar->Fill(e_reactz,e_xtar);
                hztar_yptar->Fill(e_reactz,e_yptar);
                hztar_xptar->Fill(e_reactz,e_xptar);
                hztar_xbeam->Fill(e_reactz,rastx+beamx);
		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
		Double_t Enew=p_spec*(1+e_delta/100.);
		Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		hWcalc->Fill(W_calc);
		hWdelta->Fill(e_delta,W);
		if (W<1.075) {
		hDeltaDiff->Fill(delta_diff);
		hDeltaDiffXfp->Fill(delta_diff,e_xfp);
		hDeltaDiffYfp->Fill(delta_diff,e_yfp);
		hDeltaDiffXpfp->Fill(delta_diff,e_xpfp);
		hDeltaDiffYpfp->Fill(delta_diff,e_ypfp);
		hWXfp->Fill(W,e_xfp);
		hWYfp->Fill(W,e_yfp);
		hWXpfp->Fill(W,e_xpfp);
		hWYpfp->Fill(W,e_ypfp);
		}
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
