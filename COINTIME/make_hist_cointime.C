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

void make_hist_cointime(TString basename="",Int_t nrun=2043,Int_t shms_polarity=1){
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
   inputroot="pb_ROOTfiles/"+basename+".root";
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
    Double_t HhodfpHitsTime;
    Double_t PhodfpHitsTime;
		    Double_t PhodoStartTimeMean=30.;
		    Double_t HhodoStartTimeMean=30.;
    tsimc->SetBranchAddress("H.hod.starttime", &HhodStartTime);    
    tsimc->SetBranchAddress("P.hod.starttime", &PhodStartTime);
    tsimc->SetBranchAddress("H.hod.fpHitsTime", &HhodfpHitsTime);    
    tsimc->SetBranchAddress("P.hod.fpHitsTime", &PhodfpHitsTime);
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
    Double_t pOffset = 11.5;  // in ns
    Double_t hOffset = 335.;
    Double_t SHMScentralPathLen = 18.1*100;  // SHMS Target to focal plane path length converted to cm
    Double_t SHMSpathLength = 18.1*100;      // For now assume that it's same as SHMScentralPathLen
    Double_t HMScentralPathLen = 22*100;     // HMS Target to focal plane path length converted to cm
    Double_t DeltaHMSpathLength;    
    Double_t DeltaSHMSpathLength;    
    Double_t SHMSpartMass=0.510998/1000.0; 
    Double_t HMSpartMass= 938.27208/1000.0;
  Double_t speedofLight = 29.9792; // in cm/ns                                   

// Define branches
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
 Double_t  gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("P.react.z",&e_reactz);
 Double_t  pdelta;
   tsimc->SetBranchAddress("P.gtr.dp",&pdelta);
 Double_t  hdelta;
   tsimc->SetBranchAddress("H.gtr.dp",&hdelta);
 Double_t  e_beta;
   tsimc->SetBranchAddress("P.gtr.beta",&e_beta);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&e_xptar);
 Double_t  HgtrTh;
   tsimc->SetBranchAddress("H.gtr.th",&HgtrTh);
 Double_t  PgtrTh;
   tsimc->SetBranchAddress("P.gtr.th",&PgtrTh);
 Double_t  HgtrPh;
   tsimc->SetBranchAddress("H.gtr.ph",&HgtrPh);
 Double_t  HgtrX;
   tsimc->SetBranchAddress("H.gtr.x",&HgtrX);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&e_yfp);
 Double_t  PgtrP;
   tsimc->SetBranchAddress("P.gtr.p",&PgtrP);
 Double_t  HgtrP;
   tsimc->SetBranchAddress("H.gtr.p",&HgtrP);
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
  tsimc->SetBranchAddress("H.cal.etottracknorm", &hcaletot);                                          
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
   TH1F *hetrack = new TH1F("hetrack","; HMS etottrack ", 120 ,0,1.5);
  //
  TH1D *h2_epi_PcointimeROC2    = new TH1D("h2_epi_PcointimeROC2","SHMS ROC2 Corrected epi Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_epi_PcointimeROC2_noPID    = new TH1D("h2_epi_PcointimeROC2_noPID","SHMS ROC2 Corrected epi Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ep_PcointimeROC2    = new TH1D("h2_ep_PcointimeROC2","SHMS ROC2 Corrected ep Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ep_PcointimeROC2_noPID    = new TH1D("h2_ep_PcointimeROC2_noPID","SHMS ROC2 Corrected ep Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ek_PcointimeROC2    = new TH1D("h2_ek_PcointimeROC2","SHMS ROC2 Corrected ek Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ek_PcointimeROC2_noPID    = new TH1D("h2_ek_PcointimeROC2_noPID","SHMS ROC2 Corrected ek Coin Time (no PID cut); cointime [ns]",       480, -24, 24); 
  // loop over entries
   Double_t cos_ts=TMath::Cos(14./180*3.14159);
   Double_t sin_ts=TMath::Sin(14./180*3.14159);
   Double_t Mp = .93827;
   Double_t Ei=10.6;
   Double_t p_spec=8.7*.985;
  Bool_t epievent_cut, epevent_cut, ekevent_cut, positron_cut, event_cut, hpdelta_cut;
Long64_t nentries = tsimc->GetEntries();
// nentries=150000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		//
		hetrack->Fill(hcaletot);
                hpdelta_cut = hdelta > -10 && hdelta < 10 && pdelta > -10 && pdelta < 20 ;
	       Double_t hms_elec = hcernpe > 1.0 && hcaletot > 0.6 && hcaletot < 2.0 && hcalepr > 0.2;
                event_cut = hpdelta_cut &&  hms_elec;
                epievent_cut =hms_elec && hpdelta_cut; // default 
                if (shms_polarity == 1) {
                 epievent_cut = pcalepr*PgtrP > 0.02 &&pcernpe >0.5 &&  paeronpe >= 2. ;
                 if (PgtrP>2.65) epievent_cut =pcalepr*PgtrP > 0.02 && pcaletot > 0.0 && pcaletot < 0.6&& paeronpe >= 2.&&  pcernpe >0.5 ;
                }
                if (shms_polarity == 2) {
                 epievent_cut =  pcalepr*PgtrP > 0.02 && pcernpe <=0.5 && paeronpe >= 2.;
                 if (PgtrP>=2.65 && PgtrP<4.2) epievent_cut = pcalepr*PgtrP > 0.02 && pcaletot > 0.0 && pcaletot < 0.6&& pngcernpe <=0.5&& paeronpe >= 2.&&  pcernpe >0.5 ; 
                 if (PgtrP>=4.2) epievent_cut = pcalepr*PgtrP > 0.02 && pcaletot > 0.0 && pcaletot < 0.6 && paeronpe >= 2.; 
                }
    //    epievent_cut = paeronpe > 0. && pcernpe > 0. && hcaletot > 0.6 && pcaletot > 0.6 && hpdelta_cut ;
                epevent_cut  =  pcalepr*PgtrP > 0.02 && paeronpe <= 1.0 ;
                ekevent_cut  =  pcalepr*PgtrP < 0.3 && pcalepr*PgtrP > 0.02 && pcernpe <0.5 && paeronpe >= 2.0 ;
                if (PgtrP<=3.0) ekevent_cut  = pcalepr*PgtrP > 0.02 && pcernpe <0.5 && paeronpe <=1 ;

                if (!event_cut) continue;
                if (hcaletot>.6 && TcoinpTRIG1_ROC2_tdcTimeRaw>0) {
		hSHMS_Start->Fill(PhodStartTime);
		hHMS_Start->Fill(HhodStartTime);
		}
                if (hcaletot>.6 && TcoinpTRIG6_ROC2_tdcTimeRaw>0) {
		hROC1_Trig1_4raw->Fill(TcoinpTRIG1_ROC1_tdcTimeRaw,TcoinpTRIG4_ROC1_tdcTimeRaw);
		hROC1_Trig1_4->Fill(TcoinpTRIG1_ROC1_tdcTime,TcoinpTRIG4_ROC1_tdcTime);
		hROC1_Trig1_4diff->Fill(TcoinpTRIG1_ROC1_tdcTime,0.1*(TcoinpTRIG4_ROC1_tdcTimeRaw-TcoinpTRIG1_ROC1_tdcTimeRaw));
		hStart_trigdiff->Fill(PhodStartTime-HhodStartTime,(TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime));
		hStart_coin->Fill((TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime)-(PhodStartTime-HhodStartTime),(TcoinpTRIG4_ROC1_tdcTime-TcoinpTRIG1_ROC1_tdcTime));
    Double_t PgtrBetaCalcP = PgtrP/sqrt(PgtrP*PgtrP + ( 0.93827231* 0.93827231));        
    Double_t PgtrBetaCalcK = PgtrP/sqrt(PgtrP*PgtrP + (0.497648*0.497648));        
    Double_t PgtrBetaCalcPi = PgtrP/sqrt(PgtrP*PgtrP + (0.1395704*0.1395704));        
        Double_t DeltaHMSpathLength = 12.462*HgtrTh + 0.1138*HgtrTh*HgtrX - 0.0154*HgtrX - 72.292*HgtrTh*HgtrTh - 0.0000544*HgtrX*HgtrX - 116.52*HgtrPh*HgtrPh;    
        Double_t DeltaSHMSpathLength= +0.11*atan2(PgtrTh,1)*1000;
	Double_t HgtrBetaCalc = HgtrP/sqrt(HgtrP*HgtrP + HMSpartMass*HMSpartMass);          
	Double_t HMScoinCorr = (HMScentralPathLen + DeltaHMSpathLength)/speedofLight*(1./HgtrBetaCalc-1)  - HhodfpHitsTime;     
        Double_t SHMScoinCorrP = (SHMScentralPathLen+ DeltaSHMSpathLength)/speedofLight*( 1./PgtrBetaCalcP-1) - PhodfpHitsTime; 
	HMScoinCorr = (HMScentralPathLen + DeltaHMSpathLength)/speedofLight*(1./HgtrBetaCalc-1)  - HhodStartTime;     
        SHMScoinCorrP = (SHMScentralPathLen+ DeltaSHMSpathLength)/speedofLight*( 1./PgtrBetaCalcP-1) - PhodStartTime; 
    Double_t coinP= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrP) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
    Double_t SHMScoinCorrK =   (SHMScentralPathLen+ DeltaSHMSpathLength)/speedofLight*( 1./PgtrBetaCalcK-1) - PhodfpHitsTime; 
    Double_t coinK= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrK) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
    Double_t SHMScoinCorrPi =  (SHMScentralPathLen+ DeltaSHMSpathLength)/speedofLight*( 1./PgtrBetaCalcPi-1) - PhodfpHitsTime; 
    Double_t coinPi= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrPi) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
	h2_epi_PcointimeROC2_noPID->Fill(coinPi);    
	h2_ep_PcointimeROC2_noPID->Fill(coinP);    
	h2_ek_PcointimeROC2_noPID->Fill(coinK);    
                if (epievent_cut) h2_epi_PcointimeROC2->Fill(coinPi);
                if (epevent_cut) h2_ep_PcointimeROC2->Fill(coinP);
                if (ekevent_cut) h2_ek_PcointimeROC2->Fill(coinK);
		  
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
