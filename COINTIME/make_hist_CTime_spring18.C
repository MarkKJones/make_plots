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

void make_hist_CTime_spring18(TString basename="",Int_t nrun=1267,Double_t cpeak=45){
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
   outputhist= "hist/"+basename+"_CTime_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  pTRIG1_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw",&pTRIG1_ROC1_tdcTimeRaw);
 Double_t  pTRIG4_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw",&pTRIG4_ROC1_tdcTimeRaw);
 Double_t  pTRIG1_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&pTRIG1_ROC2_tdcTimeRaw);
 Double_t  pTRIG4_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw",&pTRIG4_ROC2_tdcTimeRaw);
 Double_t  pTRIG1_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime",&pTRIG1_ROC1_tdcTime);
 Double_t  pTRIG4_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime",&pTRIG4_ROC1_tdcTime);
 Double_t  pTRIG1_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime",&pTRIG1_ROC2_tdcTime);
 Double_t  pTRIG4_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime",&pTRIG4_ROC2_tdcTime);
 Double_t  CoinTime_RAW_ROC1;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC1",&CoinTime_RAW_ROC1);
 Double_t  CoinTime_RAW_ROC2;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&CoinTime_RAW_ROC2);
 Double_t  epCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC1",&epCoinTime_ROC1);
 Double_t  epCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&epCoinTime_ROC2);
 Double_t  had_coinCorr_Positron;
   tsimc->SetBranchAddress("CTime.had_coinCorr_Positron",&had_coinCorr_Positron);
 Double_t  elec_coinCorr;
   tsimc->SetBranchAddress("CTime.elec_coinCorr",&elec_coinCorr);
 Double_t  eposCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.ePositronCoinTime_ROC1",&eposCoinTime_ROC1);
 Double_t  eposCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.ePositronCoinTime_ROC2",&eposCoinTime_ROC2);
 Double_t  ePiCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC1",&ePiCoinTime_ROC1);
 Double_t  ePiCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&ePiCoinTime_ROC2);
 Double_t  ePiCoinTime_TRIG1;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_TRIG1",&ePiCoinTime_TRIG1);
 Double_t  ePiCoinTime_TRIG4;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_TRIG4",&ePiCoinTime_TRIG4);
 Double_t  betanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
 Double_t  hbetanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&hbetanotrack);
 Double_t  pbetatrack;
   tsimc->SetBranchAddress("P.hod.beta",&pbetatrack);
 Double_t  hbetatrack;
   tsimc->SetBranchAddress("H.hod.beta",&hbetatrack);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
 Double_t  hytar;
   tsimc->SetBranchAddress("H.gtr.y",&hytar);
 Double_t  hxptar;
   tsimc->SetBranchAddress("H.gtr.th",&hxptar);
 Double_t  pxfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&pxfp);
 Double_t  pyfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&pyfp);
 Double_t  hxfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&hxfp);
 Double_t  hyfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&hyfp);
 Double_t  pntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&pntr);
 Double_t  pcer;
   tsimc->SetBranchAddress("P.hgcer.npeSum",&pcer);
 Double_t  paero;
   tsimc->SetBranchAddress("P.aero.npeSum",&paero);
 Double_t  hcer;
   tsimc->SetBranchAddress("H.cer.npeSum",&hcer);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  hetottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&hetottracknorm);
 Double_t  petottracknorm;
   tsimc->SetBranchAddress("P.cal.etottracknorm",&petottracknorm);
 Double_t  pstarttime;
   tsimc->SetBranchAddress("P.hod.starttime",&pstarttime);
 Double_t  hstarttime;
   tsimc->SetBranchAddress("H.hod.starttime",&hstarttime);
   // Define histograms
    TH1F *hpbetatrack = new TH1F("hpbetatrack",Form("Run %d ;SHMS  Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrack);
    TH1F *hhbetatrack = new TH1F("hhbetatrack",Form("Run %d ;HMS  Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetatrack);
    TH1F *hpbetanotrack = new TH1F("hpbetanotrack",Form("Run %d ;SHMS  Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetanotrack);
    TH1F *hhbetanotrack = new TH1F("hhbetanotrack",Form("Run %d ;HMS  Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetanotrack);
    TH2F *hhbetatrack_delta = new TH2F("hhbetatrack_delta",Form("Run %d ; HMS Beta ; HMS Delta",nrun), 200, .5,1.5,100,-10,10);
    HList.Add(hhbetatrack_delta);
    TH2F *hpbetatrack_delta = new TH2F("hpbetatrack_delta",Form("Run %d ; SHMS Beta ; SHMS Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetatrack_delta);
    TH2F *hpbetanotrack_delta = new TH2F("hpbetanotrack_delta",Form("Run %d ; SHMS Beta No track ; SHMS Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetanotrack_delta);
    TH1F *hpbetatrackpi = new TH1F("hpbetatrackpi",Form("Run %d ; SHMS Beta track pion;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrackpi);
    TH1F *hpbetatrackp = new TH1F("hpbetatrackp",Form("Run %d ; SHMS Beta track proton ;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrackp);
    TH1F *hpstarttime = new TH1F("hpstarttime",Form("Run %d ; SHMS Starttime;Counts",nrun), 280, 0.,100.0);
    HList.Add(hpstarttime);
    TH1F *hhstarttime = new TH1F("hhstarttime",Form("Run %d ; HMS Starttime;Counts",nrun), 280, 0.,100.0);
    HList.Add(hhstarttime);
    TH2F *hhstarttime_hdelta = new TH2F("hhstarttime_hdelta",Form("Run %d ; HMS Starttime; HMS delta",nrun), 280, -10.,60.0,100,-10,10);
    HList.Add(hhstarttime_hdelta);
    TH2F *hpstarttime_pdelta = new TH2F("hpstarttime_pdelta",Form("Run %d ; SHMS Starttime; SHMS delta",nrun), 280, -10.,60.0,100,-10,20);
    HList.Add(hpstarttime_pdelta);
    //
    TH1F *hCoinTime_RAW_ROC1 = new TH1F("hcoinTime_RAW_ROC1",Form("Run %d ; CoinTime_RAW_ROC1  ;Counts",nrun),100, -1000.,1000.0);
    HList.Add(hCoinTime_RAW_ROC1);
    TH1F *hCoinTime_RAW_ROC1_calc = new TH1F("hcoinTime_RAW_ROC1_calc",Form("Run %d ; CoinTime_RAW_ROC1_cla  ;Counts",nrun),400, -200.,+200.0);
    HList.Add(hCoinTime_RAW_ROC1_calc);
    TH1F *hCoinTime_RAW_ROC2 = new TH1F("hcoinTime_RAW_ROC2",Form("Run %d ; CoinTime_RAW_ROC2  ;Counts",nrun),100, -1000.,1000);
   TH1F *hCoinTime_RAW_ROC2_calc = new TH1F("hcoinTime_RAW_ROC2_calc",Form("Run %d ; CoinTime_RAW_ROC2_calc  ;Counts",nrun),400, -200.,200.0);
   TH2F *hCoinTime_RAW_ROC2_calc_trig1 = new TH2F("hcoinTime_RAW_ROC2_calc_trig1",Form("Run %d ; CoinTime_RAW_ROC2_calc  ; Raw trig1",nrun),100, -1000.,1000.0,100,0,1000);
   TH2F *hCoinTime_RAW_ROC2_calc_trig4 = new TH2F("hcoinTime_RAW_ROC2_calc_trig4",Form("Run %d ; CoinTime_RAW_ROC2_calc  ; Raw trig4 ",nrun),100, -1000.,1000.0,100,0,1000);
    HList.Add(hCoinTime_RAW_ROC2);
     HList.Add(hCoinTime_RAW_ROC2_calc);
     //
    TH2F *hCoinTime_RAW_ROC1_pdelta = new TH2F("hcoinTime_RAW_ROC1_pdelta",Form("Run %d ; CoinTime_RAW_ROC1  ; SHMS delta",nrun),600, -50.,100.0,100,-10,20);
    HList.Add(hCoinTime_RAW_ROC1_pdelta);
   TH1F *hepiCoinTime_ROC1 = new TH1F("hepicoinTime_ROC1",Form("Run %d ; epiCoinTimeW_ROC1  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_ROC1);
    TH1F *hepiCoinTime_ROC2 = new TH1F("hepicoinTime_ROC2",Form("Run %d ; epiCoinTimeW_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_ROC2);
   TH1F *hepiCoinTime_TRIG1 = new TH1F("hepicoinTime_TRIG1",Form("Run %d ; epiCoinTimeW_TRIG1  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_TRIG1);
    TH1F *hepiCoinTime_TRIG4 = new TH1F("hepicoinTime_TRIG4",Form("Run %d ; epiCoinTimeW_TRIG4  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_TRIG4);
    Double_t ctlo,cthi;
      ctlo=30.;
      cthi=50.;
    TH2F *hepiCoinTime_pdelta = new TH2F("hepicoinTime_pdelta",Form("Run %d ; epiCoinTimeW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(hepiCoinTime_pdelta);
    TH2F *hepiCoinTime_hdelta = new TH2F("hepicoinTime_hdelta",Form("Run %d ; epiCoinTimeW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(hepiCoinTime_hdelta);
    TH2F *hepiCoinTime_pxptar = new TH2F("hepicoinTime_pxptar",Form("Run %d ; epiCoinTimeW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepiCoinTime_pxptar);
    TH2F *hepiCoinTime_hxptar = new TH2F("hepicoinTime_hxptar",Form("Run %d ; epiCoinTimeW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepiCoinTime_hxptar);

    TH1F *hepCoinTime_ROC1 = new TH1F("hepcoinTime_ROC1",Form("Run %d ; epCoinTimeW_ROC1  ;Counts",nrun),800, 0.,100.0);
    TH1F *h_had_coinCorr_Positron = new TH1F("h_had_coinCorr_Positron",Form("Run %d ; had_coinCorr_Positron ;Counts",nrun),200, 0.,100.0);
    HList.Add(h_had_coinCorr_Positron);
    TH1F *h_elec_coinCorr = new TH1F("h_elec_coinCorr",Form("Run %d ; elec_coinCorr ;Counts",nrun),200, 0.,100.0);
    HList.Add(h_elec_coinCorr);
    TH1F *hepCoinTime_ROC2 = new TH1F("hepcoinTime_ROC2",Form("Run %d ; epCoinTimeW_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepCoinTime_ROC2);
      TH1F *heposCoinTime_ROC1 = new TH1F("heposcoinTime_ROC1",Form("Run %d ; ePositronCoinTime_ROC1  ;Counts",nrun),800, 0.,100.0);
    HList.Add(heposCoinTime_ROC1);
    TH2F *heposCoinTime_ROC2_Hetot = new TH2F("heposcoinTime_ROC2_Hetot",Form("Run %d ; ePositronCoinTime_ROC2  ;HMS Etot",nrun),120, 30.,60.0,120,0,1.5);
    HList.Add(heposCoinTime_ROC2_Hetot);
    TH2F *heposCoinTime_ROC2_Hcer = new TH2F("heposcoinTime_ROC2_Hcer",Form("Run %d ; ePositronCoinTime_ROC2  ;HMS Cer",nrun),120, 30.,60.0,80,0,20);
    HList.Add(heposCoinTime_ROC2_Hcer);
    TH2F *heposCoinTime_ROC2_Pcer = new TH2F("heposcoinTime_ROC2_Pcer",Form("Run %d ; ePositronCoinTime_ROC2  ;SHMS Cer",nrun),120, 30.,60.0,200,0,50);
    HList.Add(heposCoinTime_ROC2_Pcer);
    TH2F *heposCoinTime_ROC2_Petot = new TH2F("heposcoinTime_ROC2_Petot",Form("Run %d ; ePositronCoinTime_ROC2  ;SHMS Etot",nrun),120, 30.,60.0,120,0,1.5);
    HList.Add(heposCoinTime_ROC2_Petot);
    TH2F *heposCoinTime_ROC2_Hxfp = new TH2F("heposcoinTime_ROC2_Hxfp",Form("Run %d ; ePositronCoinTime_ROC2  ;HMS Xfp",nrun),120, 30.,60.0,80,-40,40);
    TH2F *heposCoinTime_ROC2_Hyfp = new TH2F("heposcoinTime_ROC2_Hyfp",Form("Run %d ; ePositronCoinTime_ROC2  ;HMS Yfp",nrun),120, 30.,60.0,80,-40,40);
    TH2F *heposCoinTime_ROC2_Hstime = new TH2F("heposcoinTime_ROC2_Hstime",Form("Run %d ; ePositronCoinTime_ROC2  ;HMS Starttime",nrun),120, 30.,60.0,80,0,80);
    TH2F *heposCoinTime_ROC2_Pstime = new TH2F("heposcoinTime_ROC2_Pstime",Form("Run %d ; ePositronCoinTime_ROC2  ;SHMS Starttime",nrun),120, 30.,60.0,100,20,120);
      TH1F *heposCoinTime_ROC1_cut = new TH1F("heposcoinTime_ROC1_cut",Form("Run %d ; ePositronCoinTime_ROC1 (cut) ;Counts",nrun),800, 0.,100.0);
    HList.Add(heposCoinTime_ROC1_cut);
    TH1F *heposCoinTime_ROC2 = new TH1F("heposcoinTime_ROC2",Form("Run %d ; ep=PositronCoinTime_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(heposCoinTime_ROC2);
   TH2F *hepCoinTime_pdelta = new TH2F("hepcoinTime_pdelta",Form("Run %d ; epCoinTimeW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(hepCoinTime_pdelta);
    TH2F *hepCoinTime_hdelta = new TH2F("hepcoinTime_hdelta",Form("Run %d ; epCoinTimeW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(hepCoinTime_hdelta);
    TH2F *hepCoinTime_pxptar = new TH2F("hepcoinTime_pxptar",Form("Run %d ; epCoinTimeW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepCoinTime_pxptar);
    TH2F *hepCoinTime_hxptar = new TH2F("hepcoinTime_hxptar",Form("Run %d ; epCoinTimeW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepCoinTime_hxptar);
    TH2F *hepCoinTime_pytar = new TH2F("hepcoinTime_pytar",Form("Run %d ; epCoinTimeW_ROC2  ;pytar",nrun),120, ctlo,cthi,80,-4.,4.);
    HList.Add(hepCoinTime_pytar);

    TH2F *heCoinTimeRaw_pdelta = new TH2F("heCoinTimeRaw_pdelta",Form("Run %d ; eCoinTimeRawW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(heCoinTimeRaw_pdelta);
    TH2F *heCoinTimeRaw_hdelta = new TH2F("heCoinTimeRaw_hdelta",Form("Run %d ; eCoinTimeRawW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(heCoinTimeRaw_hdelta);
    TH2F *heCoinTimeRaw_pxptar = new TH2F("heCoinTimeRaw_pxptar",Form("Run %d ; eCoinTimeRawW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(heCoinTimeRaw_pxptar);
    TH2F *heCoinTimeRaw_hxptar = new TH2F("heCoinTimeRaw_hxptar",Form("Run %d ; eCoinTimeRawW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(heCoinTimeRaw_hxptar);
// loop over entries
    Long64_t nentries = tsimc->GetEntries();
    //
    // nentries=50000;
	for (int ie = 0; ie < nentries; ie++) {
      		tsimc->GetEntry(ie);
                if (ie%50000==0) cout << " Entry = " << ie << endl;
		Bool_t selcut=pntr>0 && hntr >0;
		if (nrun == 7217) selcut = pntr>0 && hntr >0;
		if (nrun == 7220) selcut = pntr>0 && hntr >0;
		if (selcut && TMath::Abs(hdelta)<8. && delta>-10 && delta < 22.) {
		Double_t coin_ROC1_calc =  (pTRIG1_ROC2_tdcTime + pstarttime) - (pTRIG1_ROC1_tdcTime + hstarttime);
		Double_t coin_ROC2_calc =  (pTRIG4_ROC2_tdcTime + pstarttime) - (pTRIG4_ROC2_tdcTime + hstarttime);
		hpbetanotrack->Fill(betanotrack);
		hhbetanotrack->Fill(hbetanotrack);
 		hpbetatrack->Fill(pbetatrack);
		hhbetatrack->Fill(hbetatrack);
 		hpbetatrack_delta->Fill(pbetatrack,delta);
		hhbetatrack_delta->Fill(hbetatrack,hdelta);
		hpstarttime->Fill(pstarttime);
 		hhstarttime->Fill(hstarttime);
                h_had_coinCorr_Positron->Fill(had_coinCorr_Positron);
		h_elec_coinCorr->Fill(elec_coinCorr);
 		hCoinTime_RAW_ROC2->Fill(CoinTime_RAW_ROC2);
 		hCoinTime_RAW_ROC1->Fill(CoinTime_RAW_ROC1);
 		hCoinTime_RAW_ROC2_calc->Fill(coin_ROC2_calc);
  		hCoinTime_RAW_ROC2_calc_trig1->Fill(coin_ROC2_calc,pTRIG1_ROC2_tdcTimeRaw*.09776 );
 		hCoinTime_RAW_ROC2_calc_trig4->Fill(coin_ROC2_calc,pTRIG4_ROC2_tdcTimeRaw*.09776 );
		hCoinTime_RAW_ROC1_calc->Fill(coin_ROC1_calc);
		if (paero<2) {
  		hpbetatrackp->Fill(pbetatrack);
                 hepCoinTime_ROC1->Fill(epCoinTime_ROC1);
 		hepCoinTime_ROC2->Fill(epCoinTime_ROC2);
 		hepCoinTime_pdelta->Fill(epCoinTime_ROC2,delta);
 		hepCoinTime_hdelta->Fill(epCoinTime_ROC2,hdelta);
 		hepCoinTime_pxptar->Fill(epCoinTime_ROC2,xptar);
 		hepCoinTime_pytar->Fill(epCoinTime_ROC2,ytar);
 		hepCoinTime_hxptar->Fill(epCoinTime_ROC2,hxptar);
		}
                heposCoinTime_ROC1->Fill(eposCoinTime_ROC1);
 		heposCoinTime_ROC2->Fill(eposCoinTime_ROC2);
                heposCoinTime_ROC2_Hetot->Fill(eposCoinTime_ROC2,hetottracknorm);
                heposCoinTime_ROC2_Petot->Fill(eposCoinTime_ROC2,petottracknorm);
                heposCoinTime_ROC2_Hcer->Fill(eposCoinTime_ROC2,hcer);
                 heposCoinTime_ROC2_Hstime->Fill(eposCoinTime_ROC2,hstarttime);
                 heposCoinTime_ROC2_Pstime->Fill(eposCoinTime_ROC2,pstarttime);
                 heposCoinTime_ROC2_Hxfp->Fill(eposCoinTime_ROC2,hxfp);
                 heposCoinTime_ROC2_Hyfp->Fill(eposCoinTime_ROC2,hyfp);
                heposCoinTime_ROC2_Pcer->Fill(eposCoinTime_ROC2,pcer);
		if (hcer > 0.5 && hetottracknorm >0.8 ) heposCoinTime_ROC1_cut->Fill(eposCoinTime_ROC1);
		if (paero>=2) {
  		hpbetatrackpi->Fill(pbetatrack);
                   hepiCoinTime_ROC1->Fill(ePiCoinTime_ROC1);
 		 hepiCoinTime_ROC2->Fill(ePiCoinTime_ROC2);
                   hepiCoinTime_TRIG1->Fill(ePiCoinTime_TRIG1);
 		 hepiCoinTime_TRIG4->Fill(ePiCoinTime_TRIG4);
 		 hepiCoinTime_pdelta->Fill(ePiCoinTime_ROC2,delta);
 		 hepiCoinTime_hdelta->Fill(ePiCoinTime_ROC2,hdelta);
 		hepiCoinTime_pxptar->Fill(ePiCoinTime_ROC2,xptar);
 		hepiCoinTime_hxptar->Fill(ePiCoinTime_ROC2,hxptar);
		}
		heCoinTimeRaw_pdelta->Fill(CoinTime_RAW_ROC2,delta);
 		heCoinTimeRaw_hdelta->Fill(CoinTime_RAW_ROC2,hdelta);
 		heCoinTimeRaw_pxptar->Fill(CoinTime_RAW_ROC2,xptar);
 		heCoinTimeRaw_hxptar->Fill(CoinTime_RAW_ROC2,hxptar);
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
 }
