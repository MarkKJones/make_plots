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
  Double_t  p_delta;
   tsimc->SetBranchAddress("H.gtr.dp",&p_delta);
Double_t  e_mom;
   tsimc->SetBranchAddress("P.gtr.p",&e_mom);
Double_t  p_mom;
   tsimc->SetBranchAddress("H.gtr.p",&p_mom);
 Double_t  ptrig6;
   tsimc->SetBranchAddress("T.coin.pTRIG6_TdcTime",&ptrig6);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&e_xptar);
 Double_t  p_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&p_yptar);
 Double_t  p_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&p_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&e_xpfp);
 Double_t  p_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&p_yfp);
 Double_t  p_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&p_ypfp);
   Double_t  p_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&p_xfp);
 Double_t  p_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&p_xpfp);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.primary.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.primary.scat_ang_rad",&ThScat);
 Double_t  emiss;
   tsimc->SetBranchAddress("H.kin.secondary.emiss",&emiss);
 Double_t  pmiss;
   tsimc->SetBranchAddress("H.kin.secondary.pmiss",&pmiss);
 Double_t  pmissx;
   tsimc->SetBranchAddress("H.kin.secondary.pmiss_x",&pmissx);
 Double_t  pmissy;
   tsimc->SetBranchAddress("H.kin.secondary.pmiss_y",&pmissy);
 Double_t  pmissz;
   tsimc->SetBranchAddress("H.kin.secondary.pmiss_z",&pmissz);
   // Define histograms
    TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV);Counts",nrun), 200, -.05,.1);
    HList.Add(hEmiss);
    TH1F *hPmiss = new TH1F("hPmiss",Form("Run %d ; Pmiss (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmiss);
    TH1F *hMmiss2 = new TH1F("hMmiss2",Form("Run %d ; Mmiss2 (GeV2);Counts",nrun), 400, -.05,.05);
    HList.Add(hMmiss2);
    TH1F *hMmiss2_el = new TH1F("hMmiss2_el",Form("Run %d ; Mmiss2 (GeV2) W el;Counts",nrun), 100, -.05,.05);
    HList.Add(hMmiss2_el);
    TH1F *hPmissx = new TH1F("hPmissx",Form("Run %d ; Pmissx (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissx);
    TH1F *hPmissy = new TH1F("hPmissy",Form("Run %d ; Pmissy (GeV);Counts",nrun), 200, -.1,.1);
     HList.Add(hPmissy);
     TH1F *hPmissz = new TH1F("hPmissz",Form("Run %d ; Pmissz (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissz);
    TH1F *h_pdelta = new TH1F("h_pdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(h_pdelta);
    TH1F *h_edelta = new TH1F("h_edelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-10.,20.);
    HList.Add(h_edelta);
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 125, 0.8,1.3);
    HList.Add(hW);
    TH1F *hW_2 = new TH1F("hW_2",Form("Run %d ; W (GeV) (p_delta < 10);Counts",nrun), 250, 0.8,1.3);
    HList.Add(hW_2);
    TH1F *hW_3 = new TH1F("hW_3",Form("Run %d ; W (GeV) (p_delta < 8);Counts",nrun), 250, 0.8,1.3);
    HList.Add(hW_3);
      TH1F *h_pxfp = new TH1F("h_pxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(h_pxfp);
    TH1F *h_pyfp = new TH1F("h_pyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(h_pyfp);
    TH1F *h_pxpfp = new TH1F("h_pxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(h_pxpfp);
    TH1F *h_pypfp = new TH1F("h_pypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(h_pypfp);
    TH1F *h_exfp = new TH1F("h_exfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -50.,50.);
    HList.Add(h_exfp);
    TH1F *h_eyfp = new TH1F("h_eyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(h_eyfp);
   TH2F *h_eyfp_xfp = new TH2F("h_eyfp_xfp",Form("Run %d ; SHMS Y_fp; X_fp",nrun), 100, -20.,20., 100, -50.,50.);
    HList.Add(h_eyfp_xfp);
   TH2F *h_expfp_xfp = new TH2F("h_expfp_xfp",Form("Run %d ; SHMS Xp_fp; X_fp",nrun), 100, -.1,.1, 100, -50.,50.);
    HList.Add(h_expfp_xfp);
   TH2F *h_eypfp_xfp = new TH2F("h_eypfp_xfp",Form("Run %d ; SHMS Yp_fp; X_fp",nrun), 100, -.05,.05, 100, -50.,50.);
    HList.Add(h_eypfp_xfp);
   TH1F *h_expfp = new TH1F("h_expfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.12,.12);
    HList.Add(h_expfp);
    TH1F *h_eypfp = new TH1F("h_eypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(h_eypfp);
   TH1F *h_pxptar = new TH1F("h_pxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(h_pxptar);
    TH1F *h_pyptar = new TH1F("h_pyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.04,.04);
    HList.Add(h_pyptar);
    TH1F *h_exptar = new TH1F("h_exptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(h_exptar);
    TH1F *h_eyptar = new TH1F("h_eyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.04,.04);
    HList.Add(h_eyptar);
    TH1F *h_ereactz = new TH1F("h_ereactz",Form("Run %d ; SHMS react_z (cm) ;Counts",nrun), 300, -15.,15.);
    HList.Add(h_ereactz);
    TH1F *h_ereactx = new TH1F("h_ereactx",Form("Run %d ; SHMS react_x (cm);Counts",nrun), 300, -.5,.5);
    HList.Add(h_ereactx);
    TH1F *h_eytar = new TH1F("h_eytar",Form("Run %d ; SHMS Ytar (cm);Counts",nrun), 100, -5,5.);
    HList.Add(h_eytar);
    TH1F *h_hytar = new TH1F("h_hytar",Form("Run %d ; HMS Ytar (cm);Counts",nrun), 100, -5,5.);
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
    TH1F *hprot_mom_calc = new TH1F("hprot_mom_calc",Form("Run %d ; (pcalc-p)/p ;Counts",nrun), 100,-0.1,.1 );
    HList.Add(hprot_mom_calc);
  // loop over entries
    Double_t th_cent=25.6;
  Double_t Mp = .93827;
   Double_t Ei=10.600;
   Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (1==1 ) {
		  /*		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
		hW->Fill(W);
		Double_t Enew=p_spec*(1+e_delta/100.);
		Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		hWcalc->Fill(W_calc);*/
		  hetracknorm->Fill(etracknorm);
		  if (e_delta>-10. && e_delta<22. && p_delta>-8. && p_delta<8.) {
         	if (p_delta>-10. && p_delta<10.) hW_2->Fill(W);
         	if (p_delta>-8. && p_delta<8.) hW_3->Fill(W);
		hW->Fill(W);
		      Double_t theta_hms = TMath::ACos((cos_ts + p_yptar*sin_ts) / TMath::Sqrt( 1. + p_xptar*p_xptar + p_yptar * p_yptar ));
		      Double_t pcalc=2*Mp*Ei*(Ei+Mp)*cos(theta_hms)/(Mp*Mp+2*Mp*Ei+Ei*Ei*sin(theta_hms)*sin(theta_hms));
		if (W<1.05 && W>0.85) {
		      hprot_mom_calc->Fill((pcalc-p_mom)/p_mom);
                  h_exfp->Fill(e_xfp);
                  h_expfp->Fill(e_xpfp);
                  h_eyfp->Fill(e_yfp);
		  h_eyfp_xfp->Fill(e_yfp,e_xfp);
		  h_expfp_xfp->Fill(e_xpfp,e_xfp);
		  h_eypfp_xfp->Fill(e_ypfp,e_xfp);
                 h_eypfp->Fill(e_ypfp);
                  h_pxfp->Fill(p_xfp);
                  h_pxpfp->Fill(p_xpfp);
                  h_pyfp->Fill(p_yfp);
                  h_pypfp->Fill(p_ypfp);
                  h_exptar->Fill(e_xptar);
                  h_eyptar->Fill(e_yptar);
                  h_pxptar->Fill(p_xptar);
                  h_pyptar->Fill(p_yptar);
                  h_edelta->Fill(e_delta);
                  h_pdelta->Fill(p_delta);
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
		  hEmiss->Fill(emiss);		  
		  hPmiss->Fill(pmiss);		  
		  if (W<1.075) hMmiss2_el->Fill(emiss*emiss-pmiss*pmiss);	
		  hPmissx->Fill(pmissx);		  
		  hPmissy->Fill(pmissy);		  
		  hPmissz->Fill(pmissz);
		}
		//
		  }
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
