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

void make_hist_sidis(TString basename="",Int_t nrun=3288){
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
   outputhist= "hist/"+basename+"_sidis_hist.root";
 TObjArray HList(0);
 //
    
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t Ctime_raw_roc2;
 tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&Ctime_raw_roc2);
 Double_t t2_mult;
 tsimc->SetBranchAddress("T.coin.pT2_tdcMultiplicity",&t2_mult);
 Double_t t2;
 tsimc->SetBranchAddress("T.coin.pT2_tdcTimeRaw",&t2);
 Double_t trig1;
 tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&trig1);
 Double_t trig4;
 tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw",&trig4);
 Double_t  hstart_time;
 tsimc->SetBranchAddress("H.hod.starttime", &hstart_time);
 Double_t  pstart_time;
 tsimc->SetBranchAddress("P.hod.starttime", &pstart_time);
 Double_t  pTimeHist_Sigma;
 tsimc->SetBranchAddress("P.hod.TimeHist_Sigma", &pTimeHist_Sigma);
 Double_t  pTimeHist_Hits;
 tsimc->SetBranchAddress("P.hod.TimeHist_Hits", &pTimeHist_Hits);
 Double_t  pTimeHist_Peak;
 tsimc->SetBranchAddress("P.hod.TimeHist_Peak", &pTimeHist_Peak);

Double_t  etottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etottracknorm);
Double_t  etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm);
 Double_t  ghindex;
   tsimc->SetBranchAddress("H.gtr.index",&ghindex);
 Double_t  gpindex;
   tsimc->SetBranchAddress("P.gtr.index",&gpindex);
 Double_t  e_delta;
   tsimc->SetBranchAddress("H.gtr.dp",&e_delta);
 Double_t  p_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&p_delta);
 Double_t  e_mom;
   tsimc->SetBranchAddress("H.gtr.p",&e_mom);
 Double_t  p_mom;
   tsimc->SetBranchAddress("P.gtr.p",&p_mom);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&e_xptar);
 Double_t  p_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&p_yptar);
 Double_t  p_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&p_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&e_xpfp);
 Double_t  p_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&p_yfp);
 Double_t  p_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&p_ypfp);
   Double_t  p_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&p_xfp);
 Double_t  p_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&p_xpfp);
   //
//
 static const Int_t plnum=4;
 static const Int_t iside=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 const char* sidename[iside]={"Neg","Pos"};
 const char* sidename2[iside]={"neg","pos"};
 static const Int_t npad[plnum]={13,13,14,21};
 //
 Double_t plhits[plnum];
 Int_t pladchits[plnum][iside];
 Int_t pltdchits[plnum][iside];
 Double_t pltdcpad[plnum][iside][100];
 Double_t pladcpad[plnum][iside][100];
 Double_t tof_corr[plnum][iside][21];
 Double_t tw_uncorr[plnum][iside][21];
 Double_t pulseamp[plnum][iside][21];
 Double_t tdctime[plnum][iside][30];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.nhits",plname[ipl]),&plhits[ipl]) ;
  for (Int_t is=0;is<iside;is++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sTdcTimeUnCorr",plname[ipl],sidename[is]),&tw_uncorr[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Good%sTdcTimeTOFCorr",plname[ipl],sidename[is]),&tof_corr[ipl][is]) ;
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%sAdcCounter",plname[ipl],sidename2[is]),&pladchits[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.%sAdcCounter",plname[ipl],sidename2[is]),&pladcpad[ipl][is]) ;
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.%sTdcCounter",plname[ipl],sidename2[is]),&pltdchits[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.%sTdcTime",plname[ipl],sidename2[is]),&tdctime[ipl][is]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.%sTdcCounter",plname[ipl],sidename2[is]),&pltdcpad[ipl][is]) ;
     tsimc->SetBranchAddress(Form("P.hod.%s.Good%sAdcPulseAmp",plname[ipl],sidename[is]),&pulseamp[ipl][is]) ;
 }}  

// define histograms
/*
     TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
    TH1F *pxfp = new TH1F("pxfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -50.,50.);
    HList.Add(pxfp);
    TH1F *pyfp = new TH1F("pyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(pyfp);
    TH1F *pxpfp = new TH1F("pxpfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.12,.12);
    HList.Add(pxpfp);
    TH1F *pypfp = new TH1F("pypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(pypfp);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
    TH1F *pxptar = new TH1F("pxptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pxptar);
    TH1F *pyptar = new TH1F("pyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pyptar);
    TH1F *pdelta = new TH1F("pdelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(pdelta);
    TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(hdelta);
*/
//
 TH1F *hPlane_hits[plnum];
 TH1F *hPlane_padhits[plnum];
 TH1F *hPlane_rhits[plnum][iside];
 TH1F *hPlane_cut_rhits[plnum][iside];
 TH1F *hPlane_cut_hits[plnum];
 TH1F *hPlane_cut_padhits[plnum];
 TH1F *hTdcTime[plnum][iside][21];
 TH1F *hTdcTime_cut[plnum][iside][21];
 TH1F *hTdcGTime[plnum][iside][21];
 TH1F *hTdcGTime_cut[plnum][iside][21];
for (Int_t ipl=0;ipl<plnum;ipl++) {
   hPlane_hits[ipl]= new TH1F(Form("Nhits_%s",plname[ipl]),Form("%s ; Nhits; Counts",plname[ipl]),10,0,10);
   HList.Add(hPlane_hits[ipl]);
   hPlane_padhits[ipl]= new TH1F(Form("padhits_%s",plname[ipl]),Form("%s ; Paddle; Counts",plname[ipl]),npad[ipl],0,npad[ipl]);
   HList.Add(hPlane_padhits[ipl]);
   hPlane_cut_hits[ipl]= new TH1F(Form("Nhits_cut_%s",plname[ipl]),Form("%s ; Nhits; Counts",plname[ipl]),10,0,10);
   HList.Add(hPlane_cut_hits[ipl]);
   hPlane_cut_padhits[ipl]= new TH1F(Form("padhits_cut_%s",plname[ipl]),Form("%s ; Paddle; Counts",plname[ipl]),npad[ipl],0,npad[ipl]);
   HList.Add(hPlane_cut_padhits[ipl]);
 for (Int_t is=0;is<iside;is++) {
   hPlane_rhits[ipl][is]= new TH1F(Form("Rhits_%s_%s",plname[ipl],sidename[is]),Form("%s %s ; raw hits; Counts",plname[ipl],sidename[is]),10,0,10);
   HList.Add(hPlane_rhits[ipl][is]);
   hPlane_cut_rhits[ipl][is]= new TH1F(Form("Rhits_cut_%s_%s",plname[ipl],sidename[is]),Form("%s %s ; raw hits; Counts",plname[ipl],sidename[is]),10,0,10);
   HList.Add(hPlane_cut_rhits[ipl][is]);
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   hTdcTime_cut[ipl][is][ipad]= new TH1F(Form("tdctime_cut_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Tdc Time; Counts",plname[ipl],sidename[is],ipad),1200,-4000,8000.);
   HList.Add(hTdcTime_cut[ipl][is][ipad]);
   hTdcTime[ipl][is][ipad]= new TH1F(Form("tdctime_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Tdc Time; Counts",plname[ipl],sidename[is],ipad),1200,-4000,8000.);
   HList.Add(hTdcTime[ipl][is][ipad]);
   hTdcGTime_cut[ipl][is][ipad]= new TH1F(Form("tdcgtime_cut_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Tdc good Time; Counts",plname[ipl],sidename[is],ipad),1200,-400,800.);
   HList.Add(hTdcGTime_cut[ipl][is][ipad]);
   hTdcGTime[ipl][is][ipad]= new TH1F(Form("tdcgtime_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form("%s %spad_%d; Tdc good Time; Counts",plname[ipl],sidename[is],ipad),1200,-400,800.);
   HList.Add(hTdcGTime[ipl][is][ipad]);
   //cout << " added " << ipl << " " << is << " " << ipad << endl;
 }}  
 }
//
    TH1F *hEtotTrack = new TH1F("hEtotTrack",Form("Run %d ; EtotTrack;Counts",nrun), 200, 0.0,3.0);
    HList.Add(hEtotTrack);
    TH1F *hEtotTrack_gdelta = new TH1F("hEtotTrack_gdelta ",Form("Run %d ; EtotTrack (delta cut);Counts",nrun), 200, 0.0,3.0);
    HList.Add(hEtotTrack_gdelta );
    TH1F *hEtot = new TH1F("hEtot",Form("Run %d ; Etot;Counts",nrun), 200, 0.0,3.0);
    HList.Add(hEtot);
    TH1F *hCtimeRaw = new TH1F("hCtimeRaw",Form("Run %d ; Ctime raw;Counts",nrun), 4000, -1000.0,1000.0);
    HList.Add(hCtimeRaw);
    TH2F *hCtimeRaw_trig1 = new TH2F("hCtimeRaw_trig1",Form("Run %d ; Ctime raw;Trig 1 raw",nrun), 4000, -1000.0,1000.0,3000,0,8000);
    HList.Add(hCtimeRaw_trig1);
   TH2F *hCtimeRaw_trig4 = new TH2F("hCtimeRaw_trig4",Form("Run %d ; Ctime raw;Trig 4 raw",nrun), 4000, -1000.0,1000.0,3000,0,8000);
    HList.Add(hCtimeRaw_trig4);
    TH2F *hT2_trig1 = new TH2F("hT2_trig1",Form("Run %d ; T2 raw;Trig 1 raw",nrun), 4000, 0.0,8000.0,3000,0,8000);
    HList.Add(hT2_trig1);
    TH2F *hT2_trig4 = new TH2F("hT2_trig4",Form("Run %d ; T2 raw;Trig 4 raw",nrun), 4000, 0.0,8000.0,3000,0,8000);
    HList.Add(hT2_trig4);
   TH2F *hCtimeRaw_t2 = new TH2F("hCtimeRaw_t2",Form("Run %d ; Ctime raw;T2 raw",nrun), 4000, -1000.0,1000.0,1000,0,8000);
    HList.Add(hCtimeRaw_t2);
   TH2F *hCtimeRaw_t2_mult1 = new TH2F("hCtimeRaw_t2_mult1",Form("Run %d ; Ctime raw;T2 raw",nrun), 4000, -1000.0,1000.0,1000,0,8000);
    HList.Add(hCtimeRaw_t2_mult1);
   TH2F *hCtimeRaw_t2_mult2 = new TH2F("hCtimeRaw_t2_mult2",Form("Run %d ; Ctime raw;T2 raw",nrun), 4000, -1000.0,1000.0,1000,0,8000);
    HList.Add(hCtimeRaw_t2_mult2);
    TH1F *hHstarttime = new TH1F("hHstarttime",Form("Run %d ; Ctime raw;Counts",nrun), 400, -100.0,100.0);
    HList.Add(hHstarttime);
    TH1F *hPstarttime = new TH1F("hPstarttime",Form("Run %d ; SHMS starttime;Counts",nrun), 600, -100.0,200.0);
    HList.Add(hPstarttime);
    TH1F *hPTimeHist_Sigma = new TH1F("hPTimeHist_Sigma",Form("Run %d ; SHMS TimeHist_Sigma;Counts",nrun), 60, 0.0,10.0);
    HList.Add(hPTimeHist_Sigma);
    TH1F *hPTimeHist_Peak = new TH1F("hPTimeHist_Peak",Form("Run %d ; SHMS TimeHist_Peak;Counts",nrun), 200, 0.0,100.0);
    HList.Add(hPTimeHist_Peak);
    TH1F *hPTimeHist_Hits = new TH1F("hPTimeHist_Hits",Form("Run %d ; SHMS TimeHist_Hits;Counts",nrun), 10, 0.0,10.0);
    HList.Add(hPTimeHist_Hits);
    TH1F *hPTimeHist_cut_Sigma = new TH1F("hPTimeHist_cut_Sigma",Form("Run %d ; SHMS TimeHist_cut_Sigma;Counts",nrun), 60, 0.0,10.0);
    HList.Add(hPTimeHist_cut_Sigma);
    TH1F *hPTimeHist_cut_Peak = new TH1F("hPTimeHist_cut_Peak",Form("Run %d ; SHMS TimeHist_cut_Peak;Counts",nrun), 200, 0.0,100.0);
    HList.Add(hPTimeHist_cut_Peak);
    TH1F *hPTimeHist_cut_Hits = new TH1F("hPTimeHist_cut_Hits",Form("Run %d ; SHMS TimeHist_cut_Hits;Counts",nrun), 10, 0.0,10.0);
    HList.Add(hPTimeHist_cut_Hits);
    TH2F *hP_Hstarttime = new TH2F("hP_Hstarttime",Form("Run %d ; SHMS starttime;HMS start time",nrun), 600, 0.0,150.0,200,0,100.);
    HList.Add(hP_Hstarttime);
    


// loop
 Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (ghindex > -1) hEtotTrack->Fill(etottracknorm);
		if (ghindex > -1 && abs(e_delta)<=8 ) hEtotTrack_gdelta->Fill(etottracknorm);
		 hEtot->Fill(etotnorm);
		 hCtimeRaw->Fill(Ctime_raw_roc2);
		 hCtimeRaw_trig1->Fill(Ctime_raw_roc2,trig1);
		 hCtimeRaw_trig4->Fill(Ctime_raw_roc2,trig4);
		 if (Ctime_raw_roc2 < -200 ) hT2_trig1->Fill(t2,trig1);
		 if (Ctime_raw_roc2 < -200 ) hT2_trig4->Fill(t2,trig4);
		 hCtimeRaw_t2->Fill(Ctime_raw_roc2,t2);
		 if (t2_mult==1) hCtimeRaw_t2_mult1->Fill(Ctime_raw_roc2,t2);
		 if (t2_mult>1) hCtimeRaw_t2_mult2->Fill(Ctime_raw_roc2,t2);
		 hHstarttime->Fill(hstart_time);
		 hPstarttime->Fill(pstart_time);
		 hP_Hstarttime->Fill(pstart_time,hstart_time);
		 hPTimeHist_Sigma->Fill(pTimeHist_Sigma);
		 hPTimeHist_Hits->Fill(pTimeHist_Hits);
		 hPTimeHist_Peak->Fill(pTimeHist_Peak);
   for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t is=0;is<iside;is++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   hTdcGTime[ipl][is][ipad]->Fill(tw_uncorr[ipl][is][ipad]);
   if (pstart_time ==32)hTdcGTime_cut[ipl][is][ipad]->Fill(tw_uncorr[ipl][is][ipad]);
    }}}
 for (Int_t ipl=0;ipl<plnum;ipl++) {
                hPlane_hits[ipl]->Fill(plhits[ipl]);
               if (pstart_time ==32)hPlane_cut_hits[ipl]->Fill(plhits[ipl]);
 for (Int_t is=0;is<iside;is++) {
                hPlane_rhits[ipl][is]->Fill(pltdchits[ipl][is]);
               if (pstart_time ==32)hPlane_cut_rhits[ipl][is]->Fill(pltdchits[ipl][is]);
   for (Int_t nh=0;nh<pltdchits[ipl][is];nh++) {
     hPlane_padhits[ipl]->Fill(pltdcpad[ipl][is][nh]);
     if (pstart_time ==32)hPlane_cut_padhits[ipl]->Fill(pltdcpad[ipl][is][nh]);
     Int_t pad = pltdcpad[ipl][is][nh]-1;
   if (pad>=npad[ipl])  cout << " event = " << i << " hit " << nh << " " << ipl << " " << is << " " << pad << endl;
    if (pad<npad[ipl]) hTdcTime[ipl][is][pad]->Fill(tdctime[ipl][is][nh]+2000);
    if (pstart_time ==32) hTdcTime_cut[ipl][is][pad]->Fill(tdctime[ipl][is][nh]+2000);
   }
 }}
		 if (pstart_time ==32) {
		 hPTimeHist_cut_Peak->Fill(pTimeHist_Peak);
		 hPTimeHist_cut_Sigma->Fill(pTimeHist_Sigma);
		 hPTimeHist_cut_Hits->Fill(pTimeHist_Hits);
		 }
	}
	//
}
