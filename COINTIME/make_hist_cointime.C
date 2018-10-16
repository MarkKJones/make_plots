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

void make_hist_cointime(TString basename="",Int_t nrun=2043){
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
   inputroot="Online_ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_coin_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
 //
   // branches for cointime calc
    Double_t TcoinpTRIG1_ROC1_tdcTimeRaw;
    Double_t TcoinpTRIG4_ROC1_tdcTimeRaw;
    Double_t TcoinpTRIG1_ROC2_tdcTimeRaw;
    Double_t TcoinpTRIG4_ROC2_tdcTimeRaw;
    Double_t TcoinpTRIG6_ROC2_tdcTimeRaw;
    Double_t TcoinpTRIG1_ROC1_tdcTime;
    Double_t TcoinpTRIG4_ROC1_tdcTime;
    Double_t TcoinpTRIG1_ROC2_tdcTime;
    Double_t TcoinpTRIG4_ROC2_tdcTime;
    
    Double_t TcoinhTRIG1_ROC1_tdcTimeRaw;
    Double_t TcoinhTRIG4_ROC1_tdcTimeRaw;
    Double_t TcoinhTRIG1_ROC2_tdcTimeRaw;
    Double_t TcoinhTRIG4_ROC2_tdcTimeRaw;
    //Hodoscope start time 
    Double_t HhodStartTime;
    Double_t PhodStartTime;
		    Double_t PhodoStartTimeMean=30.;
		    Double_t HhodoStartTimeMean=30.;
    tsimc->SetBranchAddress("H.hod.starttime", &HhodStartTime);    
    tsimc->SetBranchAddress("P.hod.starttime", &PhodStartTime);
    tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw", &TcoinpTRIG1_ROC1_tdcTimeRaw);
    tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw", &TcoinpTRIG4_ROC1_tdcTimeRaw);
    tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw", &TcoinpTRIG1_ROC2_tdcTimeRaw);
    tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw", &TcoinpTRIG4_ROC2_tdcTimeRaw);
    tsimc->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw", &TcoinpTRIG6_ROC2_tdcTimeRaw);
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime", &TcoinpTRIG1_ROC1_tdcTime);
    tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime", &TcoinpTRIG4_ROC1_tdcTime);
    tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime", &TcoinpTRIG1_ROC2_tdcTime);
    tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime", &TcoinpTRIG4_ROC2_tdcTime);
    //
    Double_t = epCT;
    tsimc->SetBranchAddress("CT.epCoinTime_ROC2", &epCT);
    Double_t = epiCT;
    tsimc->SetBranchAddress("CT.epiCoinTime_ROC2", &epiCT);
    Double_t = ekCT;
    tsimc->SetBranchAddress("CT.ekCoinTime_ROC2", &ekCT);
    Double_t = eposCT;
    tsimc->SetBranchAddress("CT.ePositronCoinTime_ROC2", &eposCT);
    //
    Double_t pOffset = 11.5;  // in ns
    Double_t hOffset = 335.;
    Double_t speedOfLight = 29.9792; // in cm/ns
    
    Double_t SHMScentralPathLen = 18.1*100;  // SHMS Target to focal plane path length converted to cm
    Double_t SHMSpathLength = 18.1*100;      // For now assume that it's same as SHMScentralPathLen
    Double_t HMScentralPathLen = 22*100;     // HMS Target to focal plane path length converted to cm
    Double_t DeltaHMSpathLength;    
    Double_t DeltaSHMSpathLength;    
    Double_t SHMSpartMass=0.510998/1000.0; 
    Double_t HMSpartMass= 938.27208/1000.0;

// Define branches
   Double_t  hms_etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&hms_etracknorm);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
 Double_t  gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex);
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
  Double_t  pcaletot, paeronpe, pcalepr, pcernpe, pngcernpe, hcalepr, hcaletot, hcernpe;
  tsimc->SetBranchAddress("P.cal.eprtracknorm", &pcalepr);
  tsimc->SetBranchAddress("P.cal.etottracknorm", &pcaletot); 
  tsimc->SetBranchAddress("P.hgcer.npeSum", &pcernpe);                                                 
  tsimc->SetBranchAddress("P.ngcer.npeSum", &pngcernpe);                                                 
  tsimc->SetBranchAddress("P.aero.npeSum", &paeronpe);                                                 
  tsimc->SetBranchAddress("H.cal.eprtracknorm", &hcalepr);                                            
  tsimc->SetBranchAddress("H.cal.etracknorm", &hcaletot);                                          
  tsimc->SetBranchAddress("H.cer.npeSum", &hcernpe); 
   // Define histograms
  TH1D *h1_epi_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected epi Coin Time","SHMS ROC2 Corrected epi Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_epi_PcointimeROC2_noPID    = new TH1D("SHMS ROC2 Corrected epi Coin Time NO PID","SHMS ROC2 Corrected epi Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ep_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected ep Coin Time","SHMS ROC2 Corrected ep Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ep_PcointimeROC2_noPID    = new TH1D("SHMS ROC2 Corrected ep Coin Time NO PID","SHMS ROC2 Corrected ep Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ek_PcointimeROC2    = new TH1D("SHMS ROC2 Corrected ek Coin Time","SHMS ROC2 Corrected ek Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h1_ek_PcointimeROC2_noPID    = new TH1D("SHMS ROC2 Corrected ek Coin Time NO PID","SHMS ROC2 Corrected ek Coin Time (no PID cut); cointime [ns]",       480, -24, 24); 
   TH2F *hROC1_Trig1_4raw =  new TH2F("hROC1_Trig1_4raw","; ROC1 trig1 Rawtime; ROC1 Trig4 Rawtime",500,1000,4500,500,1000,4500);
   TH2F *hROC1_Trig1_4 =  new TH2F("hROC1_Trig1_4","; ROC1 trig1 Time; ROC1 Trig4 time",200,60,140,400,-200,600);
   TH2F *hROC1_Trig1_4diff =  new TH2F("hROC1_Trig1_4diff","; ROC1 trig1 Time; ROC1 Trig4 -1 time",200,0,200,400,-30,30);
   TH2F *hStart_trigdiff =  new TH2F("hStart_trigdiff","; SHMS -HMS Startime; HMS - SHMS Trig time",200,-50,50,300,-100,60);
   TH2F *hStart_coin =  new TH2F("hStart_coin","; Coin time; SHMS -HMS Startime",200,-50,50,300,-100,60);
   TH2F *hxfp_yfp =  new TH2F("hxfp_yfp","; Xfp; Yfp",200,-40,40,300,-40,40);
   TH1F *hSHMS_Start = new TH1F("hSHMS_Start","; SHMS starttime ", 120 ,0,60);
   TH1F *hHMS_Start = new TH1F("hHMS_Start","; HMS starttime ", 120 ,0,60);
  // loop over entries
   Double_t cos_ts=TMath::Cos(14./180*3.14159);
   Double_t sin_ts=TMath::Sin(14./180*3.14159);
   Double_t Mp = .93827;
   Double_t Ei=10.6;
   Double_t p_spec=8.7*.985;
Long64_t nentries = tsimc->GetEntries();
// nentries=150000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		//
                if (hms_etracknorm>.6 && TcoinpTRIG1_ROC2_tdcTimeRaw>0) {
		hSHMS_Start->Fill(PhodStartTime);
		hHMS_Start->Fill(HhodStartTime);
		}
                if (hms_etracknorm>.6 && TcoinpTRIG6_ROC2_tdcTimeRaw>0) {
		hROC1_Trig1_4raw->Fill(TcoinpTRIG1_ROC1_tdcTimeRaw,TcoinpTRIG4_ROC1_tdcTimeRaw);
		hROC1_Trig1_4->Fill(TcoinpTRIG1_ROC1_tdcTime,TcoinpTRIG4_ROC1_tdcTime);
		hROC1_Trig1_4diff->Fill(TcoinpTRIG1_ROC1_tdcTime,0.1*(TcoinpTRIG4_ROC1_tdcTimeRaw-TcoinpTRIG1_ROC1_tdcTimeRaw));
		hStart_trigdiff->Fill(PhodStartTime-HhodStartTime,(TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime));
		hStart_coin->Fill((TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime)-(PhodfpHitsTime-HhodfpHitsTime),(TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime));
		if ( PhodfpHitsTime-HhodfpHitsTime > 32 && PhodStartTime-HhodStartTime < 42 && (TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime) > 35 && (TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime) <50) {
                  hxfp_yfp->Fill(e_xfp,e_yfp);
		  
		}
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
