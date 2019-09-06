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

void make_hist_shms_hodo(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_hodo_hist.root";
 TObjArray HList(0);
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t plnum=4;
 static const Int_t clmax=10;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 Double_t ctimeraw;
 tsimc->SetBranchAddress("CTime.CoinTime_RAW_ROC2",&ctimeraw) ;
 Double_t Epictime;
 tsimc->SetBranchAddress("CTime.ePiCoinTime_ROC2",&Epictime) ;
 Double_t trig1;
 tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime",&trig1) ;
 Double_t trig4;
 tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime",&trig4) ;
 Double_t gevtyp;
 tsimc->SetBranchAddress("g.evtyp",&gevtyp) ;
 Double_t gindex;
 tsimc->SetBranchAddress("P.gtr.index",&gindex) ;
 Double_t etot;
 tsimc->SetBranchAddress("P.cal.etot",&etot) ;
 Double_t goodscinhit;
 tsimc->SetBranchAddress("P.hod.goodscinhit",&goodscinhit) ;
 Double_t hstarttime;
 tsimc->SetBranchAddress("H.hod.starttime",&hstarttime) ;
 Double_t pstarttime;
 tsimc->SetBranchAddress("P.hod.starttime",&pstarttime) ;
 Int_t ncl[plnum];
 Double_t trackdiffpos[plnum];
 Double_t ScinTransPos[plnum];
 Double_t ScinLongPos[plnum];
 Double_t trackxpos[plnum];
 Double_t trackypos[plnum];
 Double_t pos[plnum][clmax];
 Double_t size[plnum][clmax];
 Double_t flag[plnum][clmax];
 Double_t usedflag[plnum][clmax];
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("P.hod.%s.DiffDisTrack",plname[ip]),&trackdiffpos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.ScintTranversePos",plname[ip]),&ScinTransPos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.ScintLongPos",plname[ip]),&ScinLongPos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.TrackXPos",plname[ip]),&trackxpos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.TrackYPos",plname[ip]),&trackypos[ip]) ;
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.Clus.Pos",plname[ip]),&ncl[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Pos",plname[ip]),&pos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Size",plname[ip]),&size[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Flag",plname[ip]),&flag[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.UsedFlag",plname[ip]),&usedflag[ip]) ;
 }
 //
 TH1F*  hgoodscinhit = new TH1F("hgoodscinhit",";Good Scin Hit ;",10,0,10);
HList.Add(hgoodscinhit);
 TH1F*  hCtimeRaw = new TH1F("hCtimeRaw",";Raw Coin Time  ;",8000,-2000,2000);
 HList.Add(hCtimeRaw);
 TH1F*  hCalcCtimeRaw = new TH1F("hCalcCtimeRaw","; Calc  Coin Time  ;",8000,-2000,2000);
 HList.Add(hCalcCtimeRaw);
 TH1F*  hEpiCtime = new TH1F("hEpiCtime",";Epi Coin Time  ;",400,0,100);
 HList.Add(hEpiCtime);
 TH1F*  hgoodscinhit_track = new TH1F("hgoodscinhit_track",";Good Scin Hit ;",10,0,10);
 HList.Add(hgoodscinhit_track);
 TH1F*  hXposgoodscinhit = new TH1F("hXposgoodscinhit",";Scint X pos (good scint) ;",24,-60,60);
 HList.Add(hXposgoodscinhit);
 TH1F*  hXposgoodscinhit_track = new TH1F("hXposgoodscinhit_track",";Scint X pos (good scint_track) ;",24,-60,60);
 HList.Add(hXposgoodscinhit_track);
 TH1F*  hXposgoodscinhit_ratio = new TH1F("hXposgoodscinhit_ratio",";Scint X pos (good scint_track) ;",24,-60,60);
 TH1F*  hXposgoodscinhit_yposcut = new TH1F("hXposgoodscinhit_yposcut",";Scint X pos  abs(ypos)< 30 (good scint) ;",24,-60,60);
 HList.Add(hXposgoodscinhit_yposcut);
 TH1F*  hXposgoodscinhit_yposcut_track = new TH1F("hXposgoodscinhit_yposcut_track",";Scint X pos  abs(ypos)< 30 (good scint_track) ;",24,-60,60);
 HList.Add(hXposgoodscinhit_track);
 TH1F*  hXposgoodscinhit_yposcut_ratio = new TH1F("hXposgoodscinhit_yposcut_ratio",";Scint X pos  abs(ypos)< 30 (good scint_track) ;",24,-60,60);
 HList.Add(hXposgoodscinhit_ratio);
 TH1F*  hYposgoodscinhit = new TH1F("hYposgoodscinhit",";Scint Y pos abs(xpos)< 30 (good scint) ;",24,-60,60);
 HList.Add(hYposgoodscinhit);
 TH1F*  hYposgoodscinhit_track = new TH1F("hYposgoodscinhit_track",";Scint Y pos abs(xpos)< 30 (good scint_track) ;",24,-60,60);
 HList.Add(hYposgoodscinhit_track);
 TH1F*  hYposgoodscinhit_ratio = new TH1F("hYposgoodscinhit_ratio",";Scint Y pos abs(xpos)< 30 (good scint_track) ;",24,-60,60);
 HList.Add(hYposgoodscinhit_ratio);
  TH1F*  hxdiff_nclusone = new TH1F("hxdiff_nclusone","One cluster in plane; X1-X2 (cm)  ;",80,-20,20);
 TH1F*  hxdiff_nclusone_2 = new TH1F("hxdiff_nclusone_2","No shift One cluster in plane; X1-X2 (cm)  ;",80,-20,20);
 TH1F*  hydiff_nclusone = new TH1F("hydiff_nclusone","One cluster in plane; Y1-Y2 (cm)  ;",80,-20,20);
 TH2F*  hxdiff_xpos_nclusone = new TH2F("hxdiff_xpos_nclusone","One cluster in plane; X1-X2 (cm)  ; X1",120,-30,30,120,-30,30);
 TH2F*  hydiff_ypos_nclusone = new TH2F("hydiff_ypos_nclusone","One cluster in plane; Y1-Y2 (cm)  ; Y1",100,-40,40,100,-40,40);
 TH1F*  hclus[plnum];
 TH1F*  hpos[plnum];
 TH1F*  hTrackDiffpos[plnum];
 TH1F*  hTrackScinXDiffpos[plnum];
 TH1F*  hTrackScinYDiffpos[plnum];
 TH1F*  hTrackXpos_goodscin[plnum];
 TH1F*  hTrackYpos_goodscin[plnum];
 TH2F*  hTrackXYpos[plnum];
 TH2F*  hScinXYpos[plnum];
 TH2F*  hTrackXYpos_goodscin[plnum];
 TH2F*  hScinXYpos_goodscin[plnum];
 for (Int_t i=0;i<plnum;i++) {
   hclus[i] = new TH1F(Form("Hclus_%d",i),Form("; Pl %d Number of Cluster",i),10,0,10);
   hpos[i] = new TH1F(Form("Hpos_%d",i),Form("; Pl %d Position",i),100,-60,60);
   hTrackDiffpos[i] = new TH1F(Form("HTrackDiffpos_%d",i),Form("; Pl %s Diff track - hodo Position",plname[i]),100,-10,10);
 HList.Add(hTrackDiffpos[i]);
   hTrackScinXDiffpos[i] = new TH1F(Form("HTrackScinXDiffpos_%d",i),Form("; Pl %s Diff track - hodo X Position",plname[i]),100,-10,10);
 HList.Add(hTrackScinXDiffpos[i]);
   hTrackScinYDiffpos[i] = new TH1F(Form("HTrackScinYDiffpos_%d",i),Form("; Pl %s Diff track - hodo Y Position",plname[i]),100,-10,10);
 HList.Add(hTrackScinYDiffpos[i]);
   hTrackXYpos[i] = new TH2F(Form("HTrackXYpos_%d",i),Form("; Pl %s Y Position; X Position",plname[i]),100,-80,80,100,-60,60);
 HList.Add(hTrackXYpos[i]);
   hTrackXYpos_goodscin[i] = new TH2F(Form("HTrackXYpos_goodscin_%d",i),Form("goodscin; Pl %s Y Position; X Position",plname[i]),100,-80,80,100,-60,60);
 HList.Add(hTrackXYpos_goodscin[i]);
   hScinXYpos[i] = new TH2F(Form("HScinXYpos_%d",i),Form("; Pl %s Scin Y Position; Scin X Position",plname[i]),100,-80,80,100,-60,60);
 HList.Add(hScinXYpos[i]);
   hScinXYpos_goodscin[i] = new TH2F(Form("HScinXYpos_goodscin_%d",i),Form("goodscin; Pl %s Scin Y Position; Scin X Position",plname[i]),100,-80,80,100,-60,60);
 HList.Add(hScinXYpos_goodscin[i]);
   hTrackXpos_goodscin[i] = new TH1F(Form("HTrackXpos_goodscin_%d",i),Form("goodscin; Pl %s X Position; ",plname[i]),160,-80,80);
 HList.Add(hTrackXpos_goodscin[i]);
   hTrackYpos_goodscin[i] = new TH1F(Form("HTrackYpos_goodscin_%d",i),Form("goodscin; Pl %s Y Position; ",plname[i]),160,-80,80);
 HList.Add(hTrackYpos_goodscin[i]);
 }
 //
 Double_t xpfp_to_xfp=0.0018;//*0.66;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
                if (gevtyp == 4 && etot !=0) {
		hgoodscinhit->Fill(goodscinhit);
		hCtimeRaw->Fill(ctimeraw);
		hEpiCtime->Fill(Epictime);
		Double_t calc_cointime=(trig1-trig4)+(pstarttime-hstarttime);
		hCalcCtimeRaw->Fill(calc_cointime);
		if (gindex>-1) hgoodscinhit_track->Fill(goodscinhit);
		if (goodscinhit==1) hXposgoodscinhit->Fill(ScinLongPos[1]);
		if (goodscinhit==1 && gindex>-1) hXposgoodscinhit_track->Fill(ScinLongPos[1]);
		if (goodscinhit==1 && abs(ScinLongPos[0])<30) hXposgoodscinhit_yposcut->Fill(ScinLongPos[1]);
		if (goodscinhit==1 && gindex>-1 && abs(ScinLongPos[0])<30) hXposgoodscinhit_yposcut_track->Fill(ScinLongPos[1]);
		if (goodscinhit==1 && abs(ScinLongPos[1])<30) hYposgoodscinhit->Fill(ScinLongPos[0]);
		if (goodscinhit==1 && gindex>-1 && abs(ScinLongPos[1])<30) hYposgoodscinhit_track->Fill(ScinLongPos[0]);
                  for (Int_t i=0;i<plnum;i++) {
		    hTrackDiffpos[i]->Fill(trackdiffpos[i]);
		    hTrackXYpos[i]->Fill(trackypos[i],trackxpos[i]);
		    if ( i == 0 || i==2) hScinXYpos[i]->Fill(ScinLongPos[i],ScinTransPos[i]);
		    if ( i == 1 || i==3) hScinXYpos[i]->Fill(ScinTransPos[i],ScinLongPos[i]);
                    if ( i == 0 || i==2) {
		      if (trackxpos[i]!=ScinTransPos[i])hTrackScinXDiffpos[i]->Fill(trackxpos[i]-ScinTransPos[i]);
		      if (trackypos[i]!=ScinLongPos[i])hTrackScinYDiffpos[i]->Fill(trackypos[i]-ScinLongPos[i]);
		    }
                    if ( i == 1 || i==3) {
		      if (trackxpos[i]!=ScinLongPos[i])hTrackScinXDiffpos[i]->Fill(trackxpos[i]-ScinLongPos[i]);
		      if (trackypos[i]!=ScinTransPos[i])hTrackScinYDiffpos[i]->Fill(trackypos[i]-ScinTransPos[i]);
		    }
		  }
		if (goodscinhit==1) {
                  for (Int_t i=0;i<plnum;i++) {
		    hTrackXYpos_goodscin[i]->Fill(trackypos[i],trackxpos[i]);
		    hTrackXpos_goodscin[i]->Fill(trackxpos[i]);
		    hTrackYpos_goodscin[i]->Fill(trackypos[i]);
		    if ( i == 0 || i==2) hScinXYpos_goodscin[i]->Fill(ScinLongPos[i],ScinTransPos[i]);
		    if ( i == 1 || i==3) hScinXYpos_goodscin[i]->Fill(ScinTransPos[i],ScinLongPos[i]);
                    hclus[i]->Fill(ncl[i]);
		    if (ncl[i]==1) hpos[i]->Fill(pos[i][0]);
                   }
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_nclusone->Fill(pos[0][0]*(1+xpfp_to_xfp*200)-pos[2][0]);
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_nclusone_2->Fill(pos[0][0]-pos[2][0]);
		  if (ncl[1]==1&&ncl[3]==1) hydiff_nclusone->Fill(pos[1][0]-pos[3][0]);
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_xpos_nclusone->Fill(pos[0][0]*(1+xpfp_to_xfp*200)-pos[2][0],pos[0][0]);
		  if (ncl[1]==1&&ncl[3]==1) hydiff_ypos_nclusone->Fill(pos[1][0]-pos[3][0],pos[1][0]);
		}
		}
	}
//
	hXposgoodscinhit_track->Sumw2();
	hXposgoodscinhit->Sumw2();
	hXposgoodscinhit_ratio->Divide(hXposgoodscinhit_track,hXposgoodscinhit);
	hXposgoodscinhit_track->Sumw2();
	hXposgoodscinhit_yposcut->Sumw2();
	hXposgoodscinhit_yposcut_ratio->Divide(hXposgoodscinhit_yposcut_track,hXposgoodscinhit_yposcut);
	hYposgoodscinhit_track->Sumw2();
	hYposgoodscinhit->Sumw2();
	hYposgoodscinhit_ratio->Divide(hYposgoodscinhit_track,hYposgoodscinhit);
	cout << " tracking eff = " << hgoodscinhit_track->Integral(2,2) << " " << hgoodscinhit->Integral(2,2) << " " << float(hgoodscinhit_track->Integral(2,2))/float(hgoodscinhit->Integral(2,2)) << endl;
}
