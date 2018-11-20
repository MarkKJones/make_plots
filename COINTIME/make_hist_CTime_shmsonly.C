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

void make_hist_CTime_shmsonly(TString basename="",Int_t nrun=1267,Double_t cpeak=45){
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
 Double_t  pTRIG6_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTimeRaw",&pTRIG6_ROC1_tdcTimeRaw);
 Double_t  pRF_tdcTime;
   tsimc->SetBranchAddress("T.coin.pRF_tdcTime",&pRF_tdcTime);
 Double_t  pT2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pT2_tdcTime",&pT2_tdcTime);
 Double_t  pT2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pT2_tdcTimeRaw",&pT2_tdcTimeRaw);
 Double_t  pTRIG6_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime",&pTRIG6_ROC1_tdcTime);
 Double_t  pTRIG6tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_tdcTime",&pTRIG6tdcTime);
 Double_t  pTRIG4_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime",&pTRIG4_ROC1_tdcTime);
 Double_t  pTRIG4tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_tdcTimeRaw",&pTRIG4tdcTimeRaw);
 Double_t  pTRIG4_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw",&pTRIG4_ROC1_tdcTimeRaw);
 Double_t  pTRIG4tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_tdcTime",&pTRIG4tdcTime);
 Double_t  pTRIG3_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcTime",&pTRIG3_ROC1_tdcTime);
 Double_t  pTRIG3tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG3_tdcTime",&pTRIG3tdcTime);
 Double_t  pTRIG1_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime",&pTRIG1_ROC1_tdcTime);
 Double_t  pTRIG1tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_tdcTime",&pTRIG1tdcTime);
 Double_t  pTRIG1_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw",&pTRIG1_ROC1_tdcTimeRaw);
 Double_t  pTRIG1tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_tdcTimeRaw",&pTRIG1tdcTimeRaw);
 Double_t  hTRIG1_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.hTRIG1_ROC1_tdcTime",&hTRIG1_ROC1_tdcTime);
 Double_t  CoinTime_RAW_ROC1;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC1",&CoinTime_RAW_ROC1);
 Double_t  CoinTime_RAW_ROC2;
   tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&CoinTime_RAW_ROC2);
 Double_t  epCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC1",&epCoinTime_ROC1);
 Double_t  epCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.epCoinTime_ROC2",&epCoinTime_ROC2);
 Double_t  ePiCoinTime_ROC1;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC1",&ePiCoinTime_ROC1);
 Double_t  ePiCoinTime_ROC2;
   tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&ePiCoinTime_ROC2);
 Double_t  betanotrack;
   tsimc->SetBranchAddress("P.hod.betanotrack",&betanotrack);
 Double_t  hbetanotrack;
   tsimc->SetBranchAddress("H.hod.betanotrack",&hbetanotrack);
 Double_t  pbetatrack;
   tsimc->SetBranchAddress("P.hod.beta",&pbetatrack);
 Double_t  ps1xfptime;
   tsimc->SetBranchAddress("P.hod.1x.fptime",&ps1xfptime);
 Double_t  hs1xfptime;
   tsimc->SetBranchAddress("H.hod.1x.fptime",&hs1xfptime);
 Double_t  ps2xfptime;
   tsimc->SetBranchAddress("P.hod.2x.fptime",&ps2xfptime);
 Double_t  hs2xfptime;
   tsimc->SetBranchAddress("H.hod.2x.fptime",&hs2xfptime);
 Double_t  ps1yfptime;
   tsimc->SetBranchAddress("P.hod.1y.fptime",&ps1yfptime);
 Double_t  hs1yfptime;
   tsimc->SetBranchAddress("H.hod.1y.fptime",&hs1yfptime);
 Double_t  ps2yfptime;
   tsimc->SetBranchAddress("P.hod.2y.fptime",&ps2yfptime);
 Double_t  hs2yfptime;
   tsimc->SetBranchAddress("H.hod.2y.fptime",&hs2yfptime);
 Double_t  hbetatrack;
   tsimc->SetBranchAddress("H.hod.beta",&hbetatrack);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
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
 Double_t  paero;
   tsimc->SetBranchAddress("P.aero.npeSum",&paero);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  pstarttime;
   tsimc->SetBranchAddress("P.hod.starttime",&pstarttime);
 Double_t  hstarttime;
   tsimc->SetBranchAddress("H.hod.starttime",&hstarttime);
   // Define histograms
    TH1F *haero = new TH1F("haero",Form("Run %d ; Aero ;Counts",nrun), 100, -1.,50.0);
    HList.Add(haero);
    TH1F *hhs1xfptime= new TH1F("hhs1xfptime",Form("Run %d ;HMS s1X fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hhs1xfptime);
    TH1F *hhs1yfptime= new TH1F("hhs1yfptime",Form("Run %d ;HMS s1Y fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hhs1yfptime);
    TH1F *hhs2xfptime= new TH1F("hhs2xfptime",Form("Run %d ;HMS s2X fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hhs2xfptime);
    TH1F *hhs2yfptime= new TH1F("hhs2yfptime",Form("Run %d ;HMS s2Y fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hhs2yfptime);
    TH1F *hps1xfptime= new TH1F("hps1xfptime",Form("Run %d ;SHMS s1X fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hps1xfptime);
    TH1F *hps1yfptime= new TH1F("hps1yfptime",Form("Run %d ;SHMS s1Y fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hps1yfptime);
    TH1F *hps2xfptime= new TH1F("hps2xfptime",Form("Run %d ;SHMS s2X fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hps2xfptime);
    TH1F *hps2yfptime= new TH1F("hps2yfptime",Form("Run %d ;SHMS s2Y fp time ;Counts",nrun), 200,0.,50.);
    HList.Add(hps2yfptime);
    TH1F *haero_cut = new TH1F("haero_cut",Form("Run %d ; Aero (track) ;Counts",nrun), 100, -1.,50.0);

     TH2F *hps1xfptime_trig14diff = new TH2F("hps1xfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s1x fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps1xfptime_trig14diff);
     TH2F *hps2xfptime_trig14diff = new TH2F("hps2xfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s2X fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps2xfptime_trig14diff);
     TH2F *hps1yfptime_trig14diff = new TH2F("hps1yfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s1Y fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps1yfptime_trig14diff);
     TH2F *hps2yfptime_trig14diff = new TH2F("hps2yfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s2Y fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps2yfptime_trig14diff);

     TH2F *hhs1xfptime_trig14diffraw = new TH2F("hhs1xfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s1x fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs1xfptime_trig14diffraw);
     TH2F *hhs2xfptime_trig14diffraw = new TH2F("hhs2xfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s2X fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs2xfptime_trig14diffraw);
     TH2F *hhs1yfptime_trig14diffraw = new TH2F("hhs1yfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s1Y fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs1yfptime_trig14diffraw);
     TH2F *hhs2yfptime_trig14diffraw = new TH2F("hhs2yfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s2Y fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs2yfptime_trig14diffraw);

     TH2F *hps1xfptime_trig14diffraw = new TH2F("hps1xfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s1x fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps1xfptime_trig14diffraw);
     TH2F *hps2xfptime_trig14diffraw = new TH2F("hps2xfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s2X fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps2xfptime_trig14diffraw);
     TH2F *hps1yfptime_trig14diffraw = new TH2F("hps1yfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s1Y fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps1yfptime_trig14diffraw);
     TH2F *hps2yfptime_trig14diffraw = new TH2F("hps2yfptime_trig14diffraw",Form("Run %d ; raw ptrig4_roc1 - ptrig1_roc1; SHMS s2Y fp time",nrun),100, -20.,20.0, 100,20.,60.);
    HList.Add(hps2yfptime_trig14diffraw);

     TH2F *hhs1xfptime_trig14diff = new TH2F("hhs1xfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s1x fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs1xfptime_trig14diff);
     TH2F *hhs2xfptime_trig14diff = new TH2F("hhs2xfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s2X fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs2xfptime_trig14diff);
     TH2F *hhs1yfptime_trig14diff = new TH2F("hhs1yfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s1Y fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs1yfptime_trig14diff);
     TH2F *hhs2yfptime_trig14diff = new TH2F("hhs2yfptime_trig14diff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1; SHMS s2Y fp time",nrun),100, -20.,20.0, 100,0.,40.);
    HList.Add(hhs2yfptime_trig14diff);

   HList.Add(haero_cut);
    TH1F *hpbetanotrack = new TH1F("hpbetanotrack",Form("Run %d ;SHMS  Beta notrack;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetanotrack);
    TH1F *hbetanotrack_cut = new TH1F("hbetanotrack_cut",Form("Run %d ; Beta notrack (track)0000000;Counts",nrun), 600, -1.,2.0);
    HList.Add(hbetanotrack_cut);
    TH2F *hhbetatrack_delta = new TH2F("hhbetatrack_delta",Form("Run %d ; HMS Beta ; HMS Delta",nrun), 200, .5,1.5,100,-10,10);
    HList.Add(hhbetatrack_delta);
    TH2F *hpbetatrack_delta = new TH2F("hpbetatrack_delta",Form("Run %d ; SHMS Beta ; SHMS Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetatrack_delta);
    TH2F *hpbetanotrack_delta = new TH2F("hpbetanotrack_delta",Form("Run %d ; SHMS Beta No track ; SHMS Delta",nrun), 200, .5,1.5,100,-10,20);
    HList.Add(hpbetanotrack_delta);
    TH1F *hpbetatrack = new TH1F("hpbetatrack",Form("Run %d ; SHMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hpbetatrack);
    TH1F *hhbetatrack = new TH1F("hhbetatrack",Form("Run %d ; HMS Beta track;Counts",nrun), 600, -1.,2.0);
    HList.Add(hhbetatrack);
    TH1F *hpstarttime = new TH1F("hpstarttime",Form("Run %d ; SHMS Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hpstarttime);
    TH1F *hs1xdiff = new TH1F("hs1xdiff",Form("Run %d ; S1X-H1X fptimes ;Counts",nrun), 280, -60.,60.0);
    HList.Add(hs1xdiff);
    TH1F *hhstarttime = new TH1F("hhstarttime",Form("Run %d ; HMS Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hhstarttime);
    TH2F *hhstarttime_hdelta = new TH2F("hhstarttime_hdelta",Form("Run %d ; HMS Starttime; HMS delta",nrun), 280, -10.,60.0,100,-10,10);
    HList.Add(hhstarttime_hdelta);
    TH2F *hpstarttime_pdelta = new TH2F("hpstarttime_pdelta",Form("Run %d ; SHMS Starttime; SHMS delta",nrun), 280, -10.,60.0,100,-10,20);
    HList.Add(hpstarttime_pdelta);
    TH1F *hCoinTime_RAW_ROC1 = new TH1F("hcoinTime_RAW_ROC1",Form("Run %d ; CoinTime_RAW_ROC1  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC1);
    TH2F *hCoinTime_RAW_ROC1_pdelta = new TH2F("hcoinTime_RAW_ROC1_pdelta",Form("Run %d ; CoinTime_RAW_ROC1  ; SHMS delta",nrun),600, -50.,100.0,100,-10,20);
    HList.Add(hCoinTime_RAW_ROC1_pdelta);
    TH1F *hCoinTime_RAW_ROC2 = new TH1F("hcoinTime_RAW_ROC2",Form("Run %d ; CoinTime_RAW_ROC2  ;Counts",nrun),600, -50.,100.0);
    HList.Add(hCoinTime_RAW_ROC2);
    TH1F *hepiCoinTime_ROC1 = new TH1F("hepicoinTime_ROC1",Form("Run %d ; epiCoinTimeW_ROC1  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_ROC1);
    TH1F *hepiCoinTime_ROC2 = new TH1F("hepicoinTime_ROC2",Form("Run %d ; epiCoinTimeW_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTime_ROC2);
    Double_t ctlo,cthi;
    if (nrun<4000) {
      ctlo=-10.;
      cthi=20.;
    } else {
      ctlo=20.;
      cthi=50.;
    }
    TH2F *hepiCoinTime_pdelta = new TH2F("hepicoinTime_pdelta",Form("Run %d ; epiCoinTimeW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(hepiCoinTime_pdelta);
    TH2F *hepiCoinTime_hdelta = new TH2F("hepicoinTime_hdelta",Form("Run %d ; epiCoinTimeW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(hepiCoinTime_hdelta);
    TH2F *hepiCoinTime_pxptar = new TH2F("hepicoinTime_pxptar",Form("Run %d ; epiCoinTimeW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepiCoinTime_pxptar);
    TH2F *hepiCoinTime_hxptar = new TH2F("hepicoinTime_hxptar",Form("Run %d ; epiCoinTimeW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepiCoinTime_hxptar);

    TH1F *hepCoinTime_ROC1 = new TH1F("hepcoinTime_ROC1",Form("Run %d ; epCoinTimeW_ROC1  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepCoinTime_ROC1);
    TH1F *hepCoinTime_ROC2 = new TH1F("hepcoinTime_ROC2",Form("Run %d ; epCoinTimeW_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepCoinTime_ROC2);
    TH2F *hepCoinTime_pdelta = new TH2F("hepcoinTime_pdelta",Form("Run %d ; epCoinTimeW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(hepCoinTime_pdelta);
    TH2F *hepCoinTime_hdelta = new TH2F("hepcoinTime_hdelta",Form("Run %d ; epCoinTimeW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(hepCoinTime_hdelta);
    TH2F *hepCoinTime_pxptar = new TH2F("hepcoinTime_pxptar",Form("Run %d ; epCoinTimeW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepCoinTime_pxptar);
    TH2F *hepCoinTime_hxptar = new TH2F("hepcoinTime_hxptar",Form("Run %d ; epCoinTimeW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(hepCoinTime_hxptar);

    TH2F *heCoinTimeRaw_pdelta = new TH2F("heCoinTimeRaw_pdelta",Form("Run %d ; eCoinTimeRawW_ROC2  ;pdelta",nrun),120, ctlo,cthi,50,-10,20);
    HList.Add(heCoinTimeRaw_pdelta);
    TH2F *heCoinTimeRaw_hdelta = new TH2F("heCoinTimeRaw_hdelta",Form("Run %d ; eCoinTimeRawW_ROC2  ;hdelta",nrun),120, ctlo,cthi,50,-10,10);
    HList.Add(heCoinTimeRaw_hdelta);
    TH2F *heCoinTimeRaw_pxptar = new TH2F("heCoinTimeRaw_pxptar",Form("Run %d ; eCoinTimeRawW_ROC2  ;pxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(heCoinTimeRaw_pxptar);
    TH2F *heCoinTimeRaw_hxptar = new TH2F("heCoinTimeRaw_hxptar",Form("Run %d ; eCoinTimeRawW_ROC2  ;hxptar",nrun),120, ctlo,cthi,50,-.06,.06);
    HList.Add(heCoinTimeRaw_hxptar);

    TH1F *hCTcalc = new TH1F("hCTcalc",Form("Run %d ; CT Calc  ;Counts",nrun),4000, -500.,500.0);
    HList.Add(hCTcalc);
    TH1F *hCTcalc_all = new TH1F("hCTcalc_all",Form("Run %d ; CT Calc  ;Counts",nrun),400, -50.,50.0);
    HList.Add(hCTcalc_all);
    TH1F *hCTcalc_all2 = new TH1F("hCTcalc_all2",Form("Run %d ; CT Calc 2  ;Counts",nrun),400, -50.,50.0);
    HList.Add(hCTcalc_all2);
    TH1F *hepiCoinTimecut = new TH1F("hepicoinTimecut",Form("Run %d ; epiCoinTimeW_ROC2  ;Counts",nrun),800, 0.,100.0);
    HList.Add(hepiCoinTimecut);
    TH1F *hptrig1_ROC1 = new TH1F("hptrig1_ROC1",Form("Run %d ; ptrig1_ROC1 (SHMS trig in HMS ROC) ;Counts",nrun),1200, 200.,500.0);
    HList.Add(hptrig1_ROC1);
    TH1F *hptrig4_ROC1 = new TH1F("hptrig4_ROC1",Form("Run %d ; ptrig4_ROC1 (HMS trig in HMS ROC)  ;Counts",nrun),1200, 200.,500.0);
    HList.Add(hptrig4_ROC1);
    TH1F *hptrig1_ROC2 ;
    TH1F *hptrig4_ROC2 ;
    if (nrun<4000) {
       hptrig1_ROC2   = new TH1F("hptrig1_ROC2",Form("Run %d ; ptrig1_ROC2 (SHMS trig in SHMS ROC) ;Counts",nrun),40., 268.,278.0);
       hptrig4_ROC2 = new TH1F("hptrig4_ROC2",Form("Run %d ; ptrig4_ROC2 (HMS trig in SHMS ROC)  ;Counts",nrun),1200, 100.,400.0);
    } else {
       hptrig1_ROC2   = new TH1F("hptrig1_ROC2",Form("Run %d ; ptrig1_ROC2 (SHMS trig in SHMS ROC) ;Counts",nrun),40., 388.,398.0);
       hptrig4_ROC2 = new TH1F("hptrig4_ROC2",Form("Run %d ; ptrig4_ROC2 (HMS trig in SHMS ROC)  ;Counts",nrun),800, 300.,500.0);
    }
    HList.Add(hptrig1_ROC2);
    HList.Add(hptrig4_ROC2);
   TH1F *hpt2raw_ROC2 = new TH1F("hpt2raw_ROC2",Form("Run %d ; pt2raw_ROC2 (SHMS trig in SHMS ROC) ;Counts",nrun),1000, 0.,10000.0);
    HList.Add(hpt2raw_ROC2);
   TH1F *hptrig1raw_ROC2 = new TH1F("hptrig1raw_ROC2",Form("Run %d ; ptrig1raw_ROC2 (SHMS trig in SHMS ROC) ;Counts",nrun),1000, 0.,10000.0);
    HList.Add(hptrig1raw_ROC2);
   TH1F *hptrig1_t2_diff = new TH1F("hptrig1_t2_diff",Form("Run %d ; ptrig1 - pt2 raw ;Counts",nrun),2000, -1000.,1000.0);
    HList.Add(hptrig1_t2_diff);
    TH1F *hptrig4raw_ROC2 = new TH1F("hptrig4raw_ROC2",Form("Run %d ; ptrig4raw_ROC2 (HMS trig in SHMS ROC)  ;Counts",nrun),1000, 0.,10000.0);
    HList.Add(hptrig4raw_ROC2);
    TH2F *hptrig4_ptrig1_ROC1 = new TH2F("hptrig4_ptrig1_ROC1",Form("Run %d ; ptrig4_ROC1 ; ptrig1_ROC1 ",nrun),200, -50.,300.0,200, 0.,200.0);
    HList.Add(hptrig4_ptrig1_ROC1);
    TH2F *hptrig4_ptrig1cut = new TH2F("hptrig4_ptrig1cut",Form("Run %d ; ptrig4_ROC2 ; ptrig1_ROC2 ",nrun),200, 350.,450.0,200, 350.,450.0);
    HList.Add(hptrig4_ptrig1cut);
    
    //     TH2F *hptrig4_ptrigdiff_ROC1_cut = new TH2F("hptrig4_ptrigdiff_ROC1_cut",Form("Run %d ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ; T4-t1 ",nrun),200, 260.,360.0,200, -20.,20.0);
       TH2F *hptrig4_ptrigdiffcut = new TH2F("hptrig4_ptrigdiffcut",Form("Run %d ; ptrig4_ROC2 (HMS trig in SHMS ROC)  ; T4-t1 ",nrun),200, 100.,500.0,200, -20.,20.0);
    HList.Add(hptrig4_ptrigdiffcut);
       TH2F *hptrig1_ptrigdiffcut = new TH2F("hptrig1_ptrigdiffcut",Form("Run %d ; ptrig1_ROC2 (SHMS trig in SHMS ROC)  ; T4-t1 ROC2 ",nrun),1200, 200.,500.0,200, -20.,20.0);
    HList.Add(hptrig1_ptrigdiffcut);
    TH2F *hs1xdiff_ptrigdiff_ROC1_cut;
    if (nrun<4000) {
       hs1xdiff_ptrigdiff_ROC1_cut   = new TH2F("hs1xdiff_ptrigdiff_ROC1_cut",Form("Run %d ; S1x-H1x fptime  diff ; T4-t1 ",nrun),200, 0.,30.0,200, -20.,20.0);
    } else {
       hs1xdiff_ptrigdiff_ROC1_cut = new TH2F("hs1xdiff_ptrigdiff_ROC1_cut",Form("Run %d ; S1x-H1x fptime  diff ; T4-t1 ",nrun),200, 0.,30.0,200, -20.,20.0);
    }
    HList.Add(hs1xdiff_ptrigdiff_ROC1_cut);
    TH2F *hpxfp_ptrigdiff_ROC1_cut = new TH2F("hpxfp_ptrigdiff_ROC1_cut",Form("Run %d ; SHMS Xfp ; T4-t1 ",nrun),40,-40.,40.0,200, -20.,20.0);
    HList.Add(hpxfp_ptrigdiff_ROC1_cut);
    TH2F *hpyfp_ptrigdiff_ROC1_cut = new TH2F("hpyfp_ptrigdiff_ROC1_cut",Form("Run %d ; SHMS Yfp ; T4-t1 ",nrun),40,-40.,40.0,200, -20.,20.0);
    HList.Add(hpyfp_ptrigdiff_ROC1_cut);
   TH2F *hhxfp_ptrigdiff_ROC1_cut = new TH2F("hhxfp_ptrigdiff_ROC1_cut",Form("Run %d ; HMS Xfp ; T4-t1 ",nrun),40,-40.,40.0,200, -20.,20.0);
    HList.Add(hhxfp_ptrigdiff_ROC1_cut);
    TH2F *hhyfp_ptrigdiff_ROC1_cut = new TH2F("hhyfp_ptrigdiff_ROC1_cut",Form("Run %d ; HMS Yfp ; T4-t1 ",nrun),40,-40.,40.0,200, -20.,20.0);
    HList.Add(hhyfp_ptrigdiff_ROC1_cut);
    //    TH2F *hptrig4_ptrig1_ROC1_cut = new TH2F("hptrig4_ptrig1_ROC1_cut",Form("Run %d  ; ptrig4_ROC1  ; ptrig1_ROC1  ",nrun),200, 260.,360.0,200, 260.,360.0);
         TH2F *hptrig4_ptrig1_ROC1_cut = new TH2F("hptrig4_ptrig1_ROC1_cut",Form("Run %d  ; ptrig4_ROC1  ; ptrig1_ROC1  ",nrun),200, 160.,260.0,200, 150.,250.0);
    HList.Add(hptrig4_ptrig1_ROC1_cut);
    TH2F *hptrig4_htrig1_ROC1_cut = new TH2F("hptrig4_htrig1_ROC1_cut",Form("Run %d ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ; htrig1_ROC1 (HMS trig in HMS ROC) ",nrun),200, 260.,360.0,200, 0.,500.0);
    HList.Add(hptrig4_htrig1_ROC1_cut);
    TH2F *hptrig4_hstarttime_ROC1_cut = new TH2F("hptrig4_hstarttime_ROC1_cut",Form("Run %d SHMS set time ; ptrig4_ROC1 (HMS trig in SHMS ROC)  ; HMS startime ",nrun),200, 260.,360.0,200, 0.,60.0);
    HList.Add(hptrig4_hstarttime_ROC1_cut);
    TH2F *hptrig1_pstarttime_ROC1_cut = new TH2F("hptrig1_pstarttime_ROC1_cut",Form("Run %d SHMS set time ; ptrig1_ROC1 (SHMS trig in SHMS ROC)  ; SHMS startime ",nrun),200, 260.,360.0,200, 0.,60.0);
    HList.Add(hptrig1_pstarttime_ROC1_cut);
    TH2F *hptrig14diff_starttimediff_ROC1_cut ;
    TH2F *hptrig14diff_starttimediff_cut1 ;
    TH2F *hptrig13diff_starttimediff_ROC1_cut ;
    TH2F *hepiCoinTime_starttimediff;
    TH2F *hepiCoinTime_ptrig14diff;
    TH2F *hepiCoinTime_starttimediff_roc1;
    TH2F *hepiCoinTime_ptrig14diff_roc1;
    TH2F *hCTcalc_ptrig14diff;
     if (nrun<4000) {
    hptrig14diff_starttimediff_ROC1_cut= new TH2F("hptrig14diff_pstarttimediff_ROC1_cut",Form("Run %d ; ptrig4-ptrig1  ; SHMS startime-HMS starttime ",nrun),200, -10.,30.0,200, -5.,20.0);
    hptrig13diff_starttimediff_ROC1_cut= new TH2F("hptrig13diff_pstarttimediff_ROC1_cut",Form("Run %d ; ptrig3-ptrig1  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, -5.,20.0);
    hepiCoinTime_starttimediff= new TH2F("hepicoinTime_starttimediff",Form("Run %d ; epiCoinTimeW_ROC2  ;Starttime diff",nrun),800, 0.,100.0,200, -5.,20.0);
    hepiCoinTime_ptrig14diff= new TH2F("hepicoinTime_ptrig14diff",Form("Run %d ; epiCoinTimeW_ROC2  ;ST4-T1 ROC2",nrun),800, 0.,100.0,200, -10.,30.0);
     hCTcalc_ptrig14diff= new TH2F("hCTcalc_ptrig14diff",Form("Run %d ; CT_calc ;T4-T1 ROC2",nrun),800, 0.,100.0,200, -20.,20.0);
          hptrig14diff_starttimediff_cut1= new TH2F("hptrig14diff_pstarttimediff_cut1",Form("Run %d ; ptrig4-ptrig1 roc2  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, -5.,20.0);
    hepiCoinTime_starttimediff_roc1= new TH2F("hepicoinTime_starttimediff_roc1",Form("Run %d ; epiCoinTimeW_ROC1  ;Starttime diff",nrun),800, 0.,100.0,200, -5.,20.0);
      hepiCoinTime_ptrig14diff_roc1= new TH2F("hepicoinTime_ptrig14diff_roc1",Form("Run %d ; epiCoinTimeW_ROC1  ;T4-T1 ROC1",nrun),800, 0.,100.0,200, -20.,20.0);
   }
    if (nrun>5300) {
    hepiCoinTime_starttimediff= new TH2F("hepicoinTime_starttimediff",Form("Run %d ; epiCoinTimeW_ROC2  ;Starttime diff",nrun),800, 0.,100.0,200, 15.,40.0);
    hepiCoinTime_starttimediff_roc1= new TH2F("hepicoinTime_starttimediff_roc1",Form("Run %d ; epiCoinTimeW_ROC1  ;Starttime diff",nrun),800, 0.,100.0,200, 15.,40.0);
      //    hptrig14diff_starttimediff_ROC1_cut= new TH2F("hptrig14diff_pstarttimediff_ROC1_cut",Form("Run %d ; ptrig4-ptrig1  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, 5.,30.0);
          hptrig14diff_starttimediff_ROC1_cut= new TH2F("hptrig14diff_pstarttimediff_ROC1_cut",Form("Run %d ; ptrig4-ptrig1  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, 15.,40.0);
          hptrig14diff_starttimediff_cut1= new TH2F("hptrig14diff_pstarttimediff_cut1",Form("Run %d ; ptrig4-ptrig1 roc2  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, 15.,40.0);
    hptrig13diff_starttimediff_ROC1_cut= new TH2F("hptrig13diff_pstarttimediff_ROC1_cut",Form("Run %d ; ptrig3-ptrig1  ; SHMS startime-HMS starttime ",nrun),200, -20.,20.0,200, 15.,40.0);
     hepiCoinTime_ptrig14diff= new TH2F("hepicoinTime_ptrig14diff",Form("Run %d ; epiCoinTimeW_ROC2  ;T4-T1 ROC2",nrun),800, 0.,100.0,200, -20.,20.0);
     hCTcalc_ptrig14diff= new TH2F("hCTcalc_ptrig14diff",Form("Run %d ; CT_calc ;T4-T1 ROC2",nrun),800, 0.,100.0,200, -20.,20.0);
     hepiCoinTime_ptrig14diff_roc1= new TH2F("hepicoinTime_ptrig14diff_roc1",Form("Run %d ; epiCoinTimeW_ROC1  ;T4-T1 ROC1",nrun),800, 0.,100.0,200, -20.,20.0);
   }
    HList.Add(hepiCoinTime_starttimediff);
    HList.Add(hepiCoinTime_ptrig14diff);
    HList.Add(hptrig14diff_starttimediff_ROC1_cut);
    HList.Add(hptrig13diff_starttimediff_ROC1_cut);
    TH2F *hptrig6shms_ptrig6hms = new TH2F("hptrig6shms_ptrig6hms",Form("Run %d ; coin trig in SHMS ROC  ; coin trig in HMS ROC ",nrun),200, -50.,300.0,200, 0.,200.0);
    HList.Add(hptrig6shms_ptrig6hms);
    TH1F *hptrig14diff = new TH1F("hptrigd14iff",Form("Run %d ; ptrig4_roc1 - ptrig1_roc1;",nrun),200, -20.,20.0);
    HList.Add(hptrig14diff);
    TH1F *hptrig14diffroc2 = new TH1F("hptrigd14iffroc2",Form("Run %d ; ptrig4_roc2 - ptrig1_roc2;",nrun),200, -20.,20.0);
    HList.Add(hptrig14diffroc2);
 // loop over entries
    Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		Bool_t selcut=pntr>0 && hntr >0 && paero>2.;
		Double_t CT_calc=(pTRIG1tdcTimeRaw-pTRIG4tdcTimeRaw)*.0976+(pstarttime-hstarttime);
		Double_t CT_calc2=(pTRIG1tdcTimeRaw-pTRIG4tdcTimeRaw)*.0976+(-pstarttime+hstarttime);
		hpbetanotrack->Fill(betanotrack);
		haero->Fill(paero);
		hpt2raw_ROC2->Fill(pT2_tdcTimeRaw);
 		hptrig1raw_ROC2->Fill(pTRIG1tdcTimeRaw);
		hptrig1_t2_diff->Fill(pTRIG1tdcTimeRaw-pT2_tdcTimeRaw);
 		hptrig4raw_ROC2->Fill(pTRIG4tdcTimeRaw);
		if (pntr>0 && hntr >0 && paero>2.) hCTcalc_all->Fill(CT_calc);
		if (pntr>0 && hntr >0 && paero>2.) hCTcalc_all2->Fill(CT_calc2);
		if (pntr>0 && hntr >0) haero_cut->Fill(paero);
		if (nrun==3288) selcut=pntr>0 && hntr >0;
		if (selcut) {
		hps1xfptime->Fill(ps1xfptime);
		hps2xfptime->Fill(ps2xfptime);
		hps1yfptime->Fill(ps1yfptime);
		hps2yfptime->Fill(ps2yfptime);
		hhs1xfptime->Fill(hs1xfptime);
		hhs2xfptime->Fill(hs2xfptime);
		hhs1yfptime->Fill(hs1yfptime);
		hhs2yfptime->Fill(hs2yfptime);

		hps1xfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,ps1xfptime);
		hps2xfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,ps2xfptime);
		hps1yfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,ps1yfptime);
		hps2yfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,ps2yfptime);
		hhs1xfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,hs1xfptime);
		hhs2xfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,hs2xfptime);
		hhs1yfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,hs1yfptime);
		hhs2yfptime_trig14diff->Fill(pTRIG4tdcTime-pTRIG1tdcTime,hs2yfptime);

		hps1xfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,ps1xfptime);
		hps2xfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,ps2xfptime);
		hps1yfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,ps1yfptime);
		hps2yfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,ps2yfptime);
		hhs1xfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,hs1xfptime);
		hhs2xfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,hs2xfptime);
		hhs1yfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,hs1yfptime);
		hhs2yfptime_trig14diffraw->Fill((pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw)*.0976,hs2yfptime);
		hpbetatrack->Fill(pbetatrack);
		hhbetatrack->Fill(hbetatrack);
		hpstarttime->Fill(pstarttime);
		hhstarttime->Fill(hstarttime);
		hCoinTime_RAW_ROC1->Fill(CoinTime_RAW_ROC1);
		hCoinTime_RAW_ROC1_pdelta->Fill(CoinTime_RAW_ROC1,delta);
 		hptrig4_ptrig1_ROC1->Fill(pTRIG4_ROC1_tdcTime,pTRIG1_ROC1_tdcTime);
		//		if (TMath::Abs(ePiCoinTime_ROC1-cpeak)<4) {
		  hs1xdiff->Fill(ps1xfptime-hs1xfptime);
		  hs1xdiff_ptrigdiff_ROC1_cut->Fill(ps1xfptime-hs1xfptime,pTRIG4_ROC1_tdcTime-pTRIG1_ROC1_tdcTime);
		  hptrig1_ptrigdiffcut->Fill(pTRIG1tdcTime,pTRIG4tdcTime-pTRIG1tdcTime);
		hptrig1_ROC1->Fill(pTRIG1_ROC1_tdcTime);
 		hptrig4_ROC1->Fill(pTRIG4_ROC1_tdcTime);
		hptrig1_ROC2->Fill(pTRIG1tdcTime);
 		hptrig4_ROC2->Fill(pTRIG4tdcTime);
		hptrig4_ptrig1_ROC1_cut->Fill(pTRIG4_ROC1_tdcTime,pTRIG1_ROC1_tdcTime);
		  hptrig4_ptrig1cut->Fill(pTRIG4tdcTime,pTRIG1tdcTime);
		  hptrig4_ptrigdiffcut->Fill(pTRIG4_ROC1_tdcTime,pTRIG4tdcTime-pTRIG1tdcTime);
		  hptrig4_hstarttime_ROC1_cut->Fill(pTRIG4_ROC1_tdcTime,hstarttime);
		  hptrig1_pstarttime_ROC1_cut->Fill(pTRIG1_ROC1_tdcTime,pstarttime);
		  hptrig14diff_starttimediff_ROC1_cut->Fill(pTRIG4tdcTime-pTRIG1tdcTime,pstarttime-hstarttime);
		  if (TMath::Abs(pTRIG1tdcTime-391.9)< 0.1) {
                        hptrig14diff_starttimediff_cut1->Fill(pTRIG4tdcTime-pTRIG1tdcTime,pstarttime-hstarttime);
		  }
		  if (nrun==5387) hptrig13diff_starttimediff_ROC1_cut->Fill(pTRIG3_ROC1_tdcTime-pTRIG1_ROC1_tdcTime,pstarttime-hstarttime);
		  hptrig4_htrig1_ROC1_cut->Fill(pTRIG4_ROC1_tdcTime,hTRIG1_ROC1_tdcTime);	
		  hptrig14diff->Fill(0.09766*(pTRIG4_ROC1_tdcTimeRaw-pTRIG1_ROC1_tdcTimeRaw));
		  hptrig14diffroc2->Fill(0.09766*(pTRIG4tdcTimeRaw-pTRIG1tdcTimeRaw));
		  if (TMath::Abs(pTRIG4tdcTime-pTRIG1tdcTime-5.75)<0.5 && TMath::Abs(pstarttime-hstarttime-18)<2) {
		  hpxfp_ptrigdiff_ROC1_cut->Fill(pxfp,pTRIG4tdcTime-pTRIG1tdcTime);
		  hpyfp_ptrigdiff_ROC1_cut->Fill(pyfp,pTRIG4tdcTime-pTRIG1tdcTime);
		  hhxfp_ptrigdiff_ROC1_cut->Fill(hxfp,pTRIG4tdcTime-pTRIG1tdcTime);
		  hhyfp_ptrigdiff_ROC1_cut->Fill(hyfp,pTRIG4tdcTime-pTRIG1tdcTime);
		  }
		  //}
 		hCoinTime_RAW_ROC2->Fill(CoinTime_RAW_ROC2);
		hepCoinTime_ROC1->Fill(epCoinTime_ROC1);
 		hepCoinTime_ROC2->Fill(epCoinTime_ROC2);
		hepiCoinTime_ROC1->Fill(ePiCoinTime_ROC1);
 		hepiCoinTime_ROC2->Fill(ePiCoinTime_ROC2);
 		hepiCoinTime_pdelta->Fill(ePiCoinTime_ROC2,delta);
 		hepiCoinTime_hdelta->Fill(ePiCoinTime_ROC2,hdelta);
 		hepiCoinTime_pxptar->Fill(ePiCoinTime_ROC2,xptar);
 		hepiCoinTime_hxptar->Fill(ePiCoinTime_ROC2,hxptar);
		heCoinTimeRaw_pdelta->Fill(CoinTime_RAW_ROC2,delta);
 		heCoinTimeRaw_hdelta->Fill(CoinTime_RAW_ROC2,hdelta);
 		heCoinTimeRaw_pxptar->Fill(CoinTime_RAW_ROC2,xptar);
 		heCoinTimeRaw_hxptar->Fill(CoinTime_RAW_ROC2,hxptar);
 		hepCoinTime_pdelta->Fill(epCoinTime_ROC2,delta);
 		hepCoinTime_hdelta->Fill(epCoinTime_ROC2,hdelta);
 		hepCoinTime_pxptar->Fill(epCoinTime_ROC2,xptar);
 		hepCoinTime_hxptar->Fill(epCoinTime_ROC2,hxptar);
		hCTcalc->Fill(CT_calc);
                hepiCoinTime_starttimediff->Fill(ePiCoinTime_ROC2,pstarttime-hstarttime);
		hepiCoinTime_ptrig14diff->Fill(ePiCoinTime_ROC2,pTRIG4tdcTime-pTRIG1tdcTime);
		hCTcalc_ptrig14diff->Fill(CT_calc,pTRIG4tdcTime-pTRIG1tdcTime);
                hepiCoinTime_starttimediff_roc1->Fill(ePiCoinTime_ROC2,pstarttime-hstarttime);
		hepiCoinTime_ptrig14diff_roc1->Fill(ePiCoinTime_ROC1,pTRIG4_ROC1_tdcTime-pTRIG1_ROC1_tdcTime);
                if (TMath::Abs(pTRIG1tdcTime-391.9)< 0.1) hepiCoinTimecut->Fill(ePiCoinTime_ROC2);
               hhbetatrack_delta->Fill(hbetatrack,hdelta);
               hpbetatrack_delta->Fill(pbetatrack,delta);
               hpbetanotrack_delta->Fill(betanotrack,delta);
	       hptrig6shms_ptrig6hms->Fill(pTRIG6tdcTime,pTRIG6_ROC1_tdcTime);
	       hpstarttime_pdelta->Fill(pstarttime,delta);
	       hhstarttime_hdelta->Fill(hstarttime,hdelta);
		}
		//		
	}
	//
	cout << " Ratio of tracked event/total events =  " << haero_cut->Integral()/haero->Integral() << endl;
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
