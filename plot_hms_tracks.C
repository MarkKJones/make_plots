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

void plot_hms_tracks(TString basename="",Int_t nrun=2043){
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
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 Double_t  edtm_time;
   tsimc->SetBranchAddress("T.hms.hEDTM_tdcTime",&edtm_time);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm);
 Double_t  W;
   tsimc->SetBranchAddress("H.kin.W",&W);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  beta;
   tsimc->SetBranchAddress("H.hod.beta",&beta);
 Double_t  stime;
   tsimc->SetBranchAddress("H.hod.starttime",&stime);
 Double_t  track_index;
   tsimc->SetBranchAddress("H.gtr.index",&track_index);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
   Int_t ch1_numtime[6];
   Int_t ch2_numtime[6];
   TString clab[6]={"u1","u2","v1","v2","x1","x2"};
   Double_t ch1_time[6][100];
   Double_t ch2_time[6][100];
   Double_t ch1_wire[6][100];
   Double_t ch2_wire[6][100];
   TString temp;
   for (Int_t i=0;i<6;i++) {
             temp="Ndata.H.dc.1"+clab[i]+".time";
	     tsimc->SetBranchAddress(temp,&ch1_numtime[i]);
             temp="H.dc.1"+clab[i]+".time";
	     tsimc->SetBranchAddress(temp,&ch1_time[i][0]);
             temp="H.dc.1"+clab[i]+".wirenum";
	     tsimc->SetBranchAddress(temp,&ch1_wire[i][0]);
             temp="Ndata.H.dc.2"+clab[i]+".time";
	     tsimc->SetBranchAddress(temp,&ch2_numtime[i]);
             temp="H.dc.2"+clab[i]+".time";
	     tsimc->SetBranchAddress(temp,&ch2_time[i][0]);
             temp="H.dc.2"+clab[i]+".wirenum";
	     tsimc->SetBranchAddress(temp,&ch2_wire[i][0]);
	   }
   //
   TH1F *hetracknorm = new TH1F("hetracknorm", Form(" Run %d ; Track Cal E/p ; Counts",nrun), 200, 0.0,1.5);
   TH1F *hW = new TH1F("hW", Form(" Run %d ; W (GeV) ; Counts",nrun), 800, 0.8,1.5);
   TH1F *hbeta = new TH1F("hbeta", Form(" Run %d ; Beta ; Counts",nrun), 200, -.1,1.5);
   TH1F *hdelta = new TH1F("hdelta", Form(" Run %d ; Delta ; Counts",nrun), 200, -10.,15.);
   TH2F *hxfp_yfp = new TH2F("hxfp_yfp", Form(" Run %d ; X_fp ; Y_fp",nrun), 200, -50.,50., 200, -50.,50.);
   TH2F *hxfp_ypfp = new TH2F("hxfp_ypfp", Form(" Run %d ; X_fp ; Yp_fp",nrun), 200, -50.,50., 200, -.03,.03);
   TH2F *hxpfp_ypfp = new TH2F("hxpfp_ypfp", Form(" Run %d ; Xp_fp ; Yp_fp",nrun), 200, -.05,.05, 200, -.03,.03);
   TH2F *hxfp_beta = new TH2F("hxfp_beta", Form(" Run %d ;Beta ; X_fp",nrun), 200, -.1,1.5, 200, -50.,50.);
   TH2F *hstime_beta = new TH2F("hstime_beta", Form(" Run %d ;Beta ; Starttime",nrun), 200, -.1,1.5, 200, -50.,50.);
   TH2F *hyfp_beta = new TH2F("hyfp_beta", Form(" Run %d ;Beta ; Y_fp",nrun), 200, -.1,1.5, 200, -50.,50.);
   TH2F *hytar_beta = new TH2F("hytar_beta", Form(" Run %d ;Beta ; Y_tar",nrun), 200, -.1,1.5, 200, -5.,5.);
     TH1F *hch1time[6];
     TH1F *hch2time[6];
     TH2F *hch1time_wire[6];
     TH2F *hch2time_wire[6];
   for (Int_t j=0;j<6;j++) {
     temp= Form(" Run %d ; Ch1 ",nrun)+clab[j]+" Time (ns) ; Counts";
     hch1time[j] = new TH1F(Form("hch1time_%d",j+1),temp, 350,-50,250);
     temp= Form(" Run %d ; Ch1 ",nrun)+clab[j]+" Time (ns) ; Wire";
     hch1time_wire[j] = new TH2F(Form("hch1time_wire_%d",j+1),temp, 350,-50,250,100,0,100);
    temp= Form(" Run %d ; Ch2 ",nrun)+clab[j]+" Time (ns) ; Counts";
     hch2time[j] = new TH1F(Form("hch2time_%d",j+1),temp, 350,-50,250);
     temp= Form(" Run %d ; Ch2 ",nrun)+clab[j]+" Time (ns) ; Wire";
     hch2time_wire[j] = new TH2F(Form("hch2time_wire_%d",j+1),temp, 350,-50,250,100,0,100);
   }
   //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (edtm_time==0) {
		  hetracknorm->Fill(etracknorm);
		  hW->Fill(W);
		  hdelta->Fill(delta);
		  hxfp_yfp->Fill(xfp,yfp);
		  hxfp_ypfp->Fill(xfp,ypfp);
		  hxpfp_ypfp->Fill(xpfp,ypfp);
		  if (etracknorm > .7) hxfp_beta->Fill(beta,xfp);
		  if (etracknorm > .7) hyfp_beta->Fill(beta,yfp);
		  hbeta->Fill(beta);
		  hytar_beta->Fill(beta,ytar);
		  hstime_beta->Fill(beta,stime);
                  for (Int_t j=0;j<6;j++) {
		    if (ch1_numtime[j]>0) hch1time[j]->Fill(ch1_time[j][0]);
		    if (ch1_numtime[j]>0) hch1time_wire[j]->Fill(ch1_time[j][0],ch1_wire[j][0]);
		    if (ch2_numtime[j]>0) hch2time[j]->Fill(ch2_time[j][0]);
		    if (ch2_numtime[j]>0) hch2time_wire[j]->Fill(ch2_time[j][0],ch2_wire[j][0]);
		  }
		}
	}
   // plots
	//
//
}
