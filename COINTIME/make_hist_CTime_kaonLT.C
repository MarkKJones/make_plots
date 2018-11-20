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

void make_hist_CTime_kaonLT(TString basename="",Int_t nrun=1267){
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
   outputhist= "hist/"+basename+"_CTime_kaonLT_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Int_t  nhits_RF;
   tsimc->SetBranchAddress("Ndata.D.pRF_ROC2",&nhits_RF);
 Double_t  RFtime[10];
   tsimc->SetBranchAddress("D.pRF_ROC2",&RFtime);
   Double_t etot;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&etot);
   Double_t aerosum;
   tsimc->SetBranchAddress("P.aero.npeSum",&aerosum);
  Double_t hgsum;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&hgsum);
 Double_t  CoinTime_RAW_ROC1;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC1",&CoinTime_RAW_ROC1);
 Double_t  pTRIG1_time_ROC2;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime",&pTRIG1_time_ROC2);
 Double_t  pTRIG4_time_ROC2;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime",&pTRIG4_time_ROC2);
 Double_t  CoinTime_RAW_ROC2;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&CoinTime_RAW_ROC2);
 Double_t  epCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC1",&epCoinTime_ROC1);
 Double_t  epCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&epCoinTime_ROC2);
 Double_t  eKCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.eKCoinTime_ROC1",&eKCoinTime_ROC1);
 Double_t  eKCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.eKCoinTime_ROC2",&eKCoinTime_ROC2);
 Double_t  epiCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC1",&epiCoinTime_ROC1);
 Double_t  epiCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&epiCoinTime_ROC2);
 Double_t  pbetanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&pbetanotrack);
 Double_t  pbetatrack;
   tsimc->SetBranchAddress("P.hod.beta",&pbetatrack);
 Double_t  hbetanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&hbetanotrack);
 Double_t  hbetatrack;
   tsimc->SetBranchAddress("H.hod.beta",&hbetatrack);
 Double_t  pdelta;
   tsimc->SetBranchAddress("P.gtr.dp",&pdelta);
 Double_t  p_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&p_xfp);
 Double_t  p_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&p_xpfp);
 Double_t  p_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&p_yfp);
 Double_t  p_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&p_ypfp);
 Double_t  h_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&h_xfp);
 Double_t  h_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&h_xpfp);
 Double_t  h_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&h_yfp);
 Double_t  h_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&h_ypfp);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  pntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&pntr);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  hstarttime;
   tsimc->SetBranchAddress("H.hod.starttime",&hstarttime);
 Double_t  pstarttime;
   tsimc->SetBranchAddress("P.hod.starttime",&pstarttime);
   // Define histograms
    TH1F *hhbetanotrack = new TH1F("hhbetanotrack",Form("Run %d ; HMS Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetanotrack);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; SHMS etot;Counts",nrun), 600, 0.,2.0);
    HList.Add(hhbetanotrack);
    TH2F *hhbetanotrack_delta = new TH2F("hhbetanotrack_delta",Form("Run %d ; HMS Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,10);
    HList.Add(hhbetanotrack_delta);
    TH2F *hhbetatrack_delta = new TH2F("hhbetatrack_delta",Form("Run %d ; HMS Beta track;Delta",nrun), 200, .5,1.5,100,-10,10);
    HList.Add(hhbetatrack_delta);
    TH1F *hhbetatrack = new TH1F("hhbetatrack",Form("Run %d ; HMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetatrack);
    TH1F *hpbetanotrack = new TH1F("hpbetanotrack",Form("Run %d ; SHMS Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetanotrack);
    TH2F *hpbetanotrack_delta = new TH2F("hpbetanotrack_delta",Form("Run %d ; SHMS Beta notrack;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetanotrack_delta);
    TH2F *hpbetatrack_delta = new TH2F("hpbetatrack_delta",Form("Run %d ; SHMS Beta track;Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetatrack_delta);
    TH1F *hpbetatrack = new TH1F("hpbetatrack",Form("Run %d ; SHMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrack);
    TH1F *hhstarttime = new TH1F("hhstarttime",Form("Run %d ; HMS Starttime;Counts",nrun), 280, -10.,80.0);
    HList.Add(hhstarttime);
    TH1F *hpstarttime = new TH1F("hpstarttime",Form("Run %d ; SHMS Starttime;Counts",nrun), 280, -10.,80.0);
    HList.Add(hpstarttime);
    TH2F *hpstarttime_pdelta = new TH2F("hpstarttime_pdelta",Form("Run %d ; SHMS Starttime; SHMS delta",nrun), 280, -10.,80.0,100,-10,20);
    HList.Add(hpstarttime_pdelta);
    TH2F *hhstarttime_hdelta = new TH2F("hhstarttime_pdelta",Form("Run %d ; HMS Starttime; HMS delta",nrun), 280, -10.,80.0,100,-10,10);
    HList.Add(hhstarttime_hdelta);
    TH1F *hCoinTime_RAW_ROC1 = new TH1F("hcoinTime_RAW_ROC1",Form("Run %d ; CoinTime_RAW_ROC1  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC1);
    TH2F *hCoinTime_RAW_ROC1_pdelta = new TH2F("hcoinTime_RAW_ROC1_pdelta",Form("Run %d ; CoinTime_RAW_ROC1  ; SHMS delta",nrun),600, -50.,100.0,100,-10,20);
    HList.Add(hCoinTime_RAW_ROC1_pdelta);
    TH1F *hCoinTime_RAW_ROC2 = new TH1F("hcoinTime_RAW_ROC2",Form("Run %d ; CoinTime_RAW_ROC2  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC2);
    TH1F *hepCoinTime_ROC1 = new TH1F("hepcoinTime_ROC1",Form("Run %d ; epCoinTime_ROC1  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepCoinTime_ROC1);
    TH1F *hepCoinTime_ROC2 = new TH1F("hepcoinTime_ROC2",Form("Run %d ; epCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepCoinTime_ROC2);
    TH2F *hepiCoinTime_ROC2_pdelta = new TH2F("hepicoinTime_ROC2_pdelta",Form("Run %d ; epiCoinTime_ROC2  ; SHMS delta",nrun),100,40.,50.0,60,-10,20);
    HList.Add(hepiCoinTime_ROC2_pdelta);
    TH2F *hepiCoinTime_ROC2_pxfp = new TH2F("hepicoinTime_ROC2_pxfp",Form("Run %d ; epiCoinTime_ROC2  ; SHMS xfp",nrun),40,40.,50.0,60,-50,50);
    TH2F *hepiCoinTime_ROC2_pyfp = new TH2F("hepicoinTime_ROC2_pyfp",Form("Run %d ; epiCoinTime_ROC2  ; SHMS yfp",nrun),40,40.,50.0,60,-30,30);
    HList.Add(hepiCoinTime_ROC2_pyfp);
    TH2F *hepiCoinTime_ROC2_pxpfp = new TH2F("hepicoinTime_ROC2_pxpfp",Form("Run %d ; epiCoinTime_ROC2  ; SHMS xpfp",nrun),40,40.,50.0,60,-.15,.15);
    HList.Add(hepiCoinTime_ROC2_pxpfp);
    TH2F *hepiCoinTime_ROC2_pypfp = new TH2F("hepicoinTime_ROC2_pypfp",Form("Run %d ; epiCoinTime_ROC2  ; SHMS ypfp",nrun),40,40.,50.0,60,-.05,.05);
    HList.Add(hepiCoinTime_ROC2_pypfp);
    TH2F *hRawCoinTime_ROC2_pxfp = new TH2F("hRawcoinTime_ROC2_pxfp",Form("Run %d ; RawCoinTime_ROC2  ; SHMS xfp",nrun),40,40.,50.0,60,-50,50);
    TH2F *hRawCoinTime_ROC2_pyfp = new TH2F("hRawcoinTime_ROC2_pyfp",Form("Run %d ; RawCoinTime_ROC2  ; SHMS yfp",nrun),40,40.,50.0,60,-30,30);
    HList.Add(hRawCoinTime_ROC2_pyfp);
    TH2F *hRawCoinTime_ROC2_pxpfp = new TH2F("hRawcoinTime_ROC2_pxpfp",Form("Run %d ; RawCoinTime_ROC2  ; SHMS xpfp",nrun),40,40.,50.0,60,-.15,.15);
    HList.Add(hRawCoinTime_ROC2_pxpfp);
    TH2F *hRawCoinTime_ROC2_pypfp = new TH2F("hRawcoinTime_ROC2_pypfp",Form("Run %d ; RawCoinTime_ROC2  ; SHMS ypfp",nrun),40,40.,50.0,60,-.05,.05);
    HList.Add(hRawCoinTime_ROC2_pypfp);
    TH2F *hRawCoinTime_ROC2_hdelta = new TH2F("hRawcoinTime_ROC2_hdelta",Form("Run %d ; RawCoinTime_ROC2  ; HMS delta",nrun),40,40.,50.0,40,-10,10);
    HList.Add(hRawCoinTime_ROC2_hdelta);
    TH2F *hepiCoinTime_ROC2_hdelta = new TH2F("hepicoinTime_ROC2_hdelta",Form("Run %d ; epiCoinTime_ROC2  ; HMS delta",nrun),40,40.,50.0,40,-10,10);
    HList.Add(hepiCoinTime_ROC2_hdelta);
    TH2F *hepiCoinTime_ROC2_hxfp = new TH2F("hepicoinTime_ROC2_hxfp",Form("Run %d ; epiCoinTime_ROC2  ; HMS xfp",nrun),40,40.,50.0,60,-50,50);
    TH2F *hepiCoinTime_ROC2_hyfp = new TH2F("hepicoinTime_ROC2_hyfp",Form("Run %d ; epiCoinTime_ROC2  ; HMS yfp",nrun),40,40.,50.0,60,-30,30);
    HList.Add(hepiCoinTime_ROC2_hyfp);
    TH2F *hepiCoinTime_ROC2_hxpfp = new TH2F("hepicoinTime_ROC2_hxpfp",Form("Run %d ; epiCoinTime_ROC2  ; HMS xpfp",nrun),40,40.,50.0,60,-.15,.15);
    HList.Add(hepiCoinTime_ROC2_hxpfp);
    TH2F *hepiCoinTime_ROC2_hypfp = new TH2F("hepicoinTime_ROC2_hypfp",Form("Run %d ; epiCoinTime_ROC2  ; HMS ypfp",nrun),40,40.,50.0,60,-.05,.05);
    HList.Add(hepiCoinTime_ROC2_hypfp);
    TH1F *hepiCoinTime_ROC2 = new TH1F("hepicoinTime_ROC2",Form("Run %d ; epiCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(hepiCoinTime_ROC2);
    TH1F *heKCoinTime_ROC2 = new TH1F("heKcoinTime_ROC2",Form("Run %d ; eKCoinTime_ROC2  ;Counts",nrun),800,0.,100.0);
    HList.Add(heKCoinTime_ROC2);
    TH2F *heKCoinTime_hstarttime = new TH2F("heKcoinTime_hstarttime",Form("Run %d ; eKCoinTime_ROC2  ; HMS startitime",nrun),160,20.,60.0,160,20.,60.0);
    HList.Add(heKCoinTime_hstarttime);
     TH2F *heKCoinTime_pstarttime = new TH2F("heKcoinTime_pstarttime",Form("Run %d ; eKCoinTime_ROC2  ; SHMS startitime",nrun),160,20.,60.0,160,20.,60.0);
    HList.Add(heKCoinTime_pstarttime);
    TH2F *hRawCoinTime_hstarttime = new TH2F("hRawcoinTime_hstarttime",Form("Run %d ; Raw CoinTime_ROC2  ; HMS startitime",nrun),160,20.,60.0,160,20.,60.0);
    HList.Add(hRawCoinTime_hstarttime);
    TH2F *hTrig4_hstarttime = new TH2F("hTrig4_hstarttime",Form("Run %d ; Trig4  ; HMS starttime",nrun),200,200.,600.0,160,20.,60.0);
    HList.Add(hTrig4_hstarttime);
    TH2F *hTrig4_pstarttime = new TH2F("hTrig4_pstarttime",Form("Run %d ; Trig4  ; SHMS starttime",nrun),200,200.,600.0,160,20.,60.0);
    HList.Add(hTrig4_pstarttime);
     TH2F *hRawCoinTime_pstarttime = new TH2F("hRawcoinTime_pstarttime",Form("Run %d ; Raw CoinTime_ROC2  ; SHMS startitime",nrun),160,20.,60.0,160,20.,60.0);
    HList.Add(hRawCoinTime_pstarttime);
    TH2F *hRawCoinTime_pbeta = new TH2F("hRawcoinTime_pbeta",Form("Run %d ; Raw CoinTime_ROC2  ; SHMS beta",nrun),160,20.,60.0,80,.8,1.2);
    HList.Add(hRawCoinTime_pbeta);
    TH2F *hRawCoinTime_hbeta = new TH2F("hRawcoinTime_hbeta",Form("Run %d ; Raw CoinTime_ROC2  ; HMS beta",nrun),160,20.,60.0,80,.8,1.2);
    HList.Add(hRawCoinTime_hbeta);
    TH2F *hRawCoinTime_calc = new TH2F("hRawcoinTime_calc",Form("Run %d ; Raw CoinTime_ROC2  ; calc",nrun),160,20.,60.0,160,20.,60.);
    HList.Add(hRawCoinTime_calc);
    TH1F *hRawCoinTime_calc_diff = new TH1F("hRawcoinTime_calc_diff",Form("Run %d ; Raw CoinTime_ROC2-calc",nrun),160,-10.,10.0);
    HList.Add(hRawCoinTime_calc_diff);
    TH1F *hRawCoinTimecalc = new TH1F("hRawcoinTimecalc",Form("Run %d ; Raw Calc",nrun),160,20.,60.0);
    HList.Add(hRawCoinTimecalc);
    TH2F *hTrig14diff_starttimediff = new TH2F("hTrig14diff_starttimediff",Form("Run %d ; Trig4   - trig1; H-SHMSMS starttime",nrun),200,-20.,20.0,160,-30.,-15.);
    HList.Add(hTrig14diff_starttimediff);
 // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                hetot->Fill(etot);
		if (etot>0.05 && aerosum<2.) {
		hhbetanotrack->Fill(hbetanotrack);
		hhbetatrack->Fill(hbetatrack);
		hhstarttime->Fill(hstarttime);
		hpbetanotrack->Fill(pbetanotrack);
		hpbetatrack->Fill(pbetatrack);
		hpstarttime->Fill(pstarttime);
                hhbetanotrack_delta->Fill(hbetanotrack,hdelta);
                hpbetanotrack_delta->Fill(pbetanotrack,pdelta);
                hhbetatrack_delta->Fill(hbetatrack,hdelta);
                hpbetatrack_delta->Fill(pbetatrack,pdelta);
		hCoinTime_RAW_ROC1->Fill(CoinTime_RAW_ROC1);
		hCoinTime_RAW_ROC1_pdelta->Fill(CoinTime_RAW_ROC1,pdelta);
  		hCoinTime_RAW_ROC2->Fill(CoinTime_RAW_ROC2);
		hepCoinTime_ROC1->Fill(epCoinTime_ROC1);
 		hepCoinTime_ROC2->Fill(epCoinTime_ROC2);
 		hepiCoinTime_ROC2_hdelta->Fill(epiCoinTime_ROC2,hdelta);
 		hepiCoinTime_ROC2_pdelta->Fill(epiCoinTime_ROC2,pdelta);
 		hepiCoinTime_ROC2_pxfp->Fill(epiCoinTime_ROC2,p_xfp);
 		hepiCoinTime_ROC2_pxpfp->Fill(epiCoinTime_ROC2,p_xpfp);
 		hepiCoinTime_ROC2_pyfp->Fill(epiCoinTime_ROC2,p_yfp);
 		hepiCoinTime_ROC2_pypfp->Fill(epiCoinTime_ROC2,p_ypfp);
 		hepiCoinTime_ROC2_hxfp->Fill(epiCoinTime_ROC2,h_xfp);
 		hepiCoinTime_ROC2_hxpfp->Fill(epiCoinTime_ROC2,h_xpfp);
 		hepiCoinTime_ROC2_hyfp->Fill(epiCoinTime_ROC2,h_yfp);
 		hepiCoinTime_ROC2_hypfp->Fill(epiCoinTime_ROC2,h_ypfp);
 		hepiCoinTime_ROC2_hdelta->Fill(epiCoinTime_ROC2,hdelta);
 		hepiCoinTime_ROC2->Fill(epiCoinTime_ROC2);
 		heKCoinTime_ROC2->Fill(eKCoinTime_ROC2);
		heKCoinTime_hstarttime->Fill(eKCoinTime_ROC2,hstarttime);
		heKCoinTime_pstarttime->Fill(eKCoinTime_ROC2,pstarttime);
		hRawCoinTime_hstarttime->Fill(CoinTime_RAW_ROC2,hstarttime);
		hRawCoinTime_pstarttime->Fill(CoinTime_RAW_ROC2,pstarttime);
		hRawCoinTime_pbeta->Fill(CoinTime_RAW_ROC2,pbetatrack);
		hRawCoinTime_hbeta->Fill(CoinTime_RAW_ROC2,hbetatrack);
		hRawCoinTime_calc->Fill(CoinTime_RAW_ROC2,-pTRIG4_time_ROC2+pTRIG1_time_ROC2-hstarttime+pstarttime);
		hRawCoinTime_calc_diff->Fill(CoinTime_RAW_ROC2-(-pTRIG4_time_ROC2+pTRIG1_time_ROC2-hstarttime+pstarttime));
		hRawCoinTimecalc->Fill((-pTRIG4_time_ROC2+pTRIG1_time_ROC2-hstarttime+pstarttime));
		hTrig14diff_starttimediff->Fill(pTRIG4_time_ROC2-pTRIG1_time_ROC2,hstarttime-pstarttime);
		hTrig4_hstarttime->Fill(pTRIG4_time_ROC2,hstarttime);
		hTrig4_pstarttime->Fill(pTRIG4_time_ROC2,pstarttime);
	       hpstarttime_pdelta->Fill(pstarttime,pdelta);
	       hhstarttime_hdelta->Fill(hstarttime,hdelta);
		  }
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
