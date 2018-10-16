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

void make_hist_cointime2(TString basename="",Int_t nrun=2043,Int_t shms_polarity=1){
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
    
    //Hodoscope start time 
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
    Double_t epCT;
    tsimc->SetBranchAddress("CT.epCoinTime_ROC2", &epCT);
    Double_t epiCT;
    tsimc->SetBranchAddress("CT.ePiCoinTime_ROC2", &epiCT);
    Double_t ekCT;
    tsimc->SetBranchAddress("CT.eKCoinTime_ROC2", &ekCT);
    Double_t eposCT;
    tsimc->SetBranchAddress("CT.ePositronCoinTime_ROC2", &eposCT);
    //

// Define branches
  Double_t  pcaletot, paeronpe, pcalepr, pcernpe, pngcernpe, hcalepr, hcaletot, hcernpe;
  Double_t HgtrX, HgtrTh, HgtrY, HgtrPh, hdelta, PgtrX, PgtrP, HgtrP, PgtrTh, PgtrY, PgtrPh,pbeta,hbeta, pdelta;
  Double_t HhodStartTime, HhodfpHitsTime, PhodStartTime, PhodfpHitsTime;
  tsimc->SetBranchAddress("P.hod.starttime", &PhodStartTime);                                               
  tsimc->SetBranchAddress("P.hod.fpHitsTime", &PhodfpHitsTime);                                             
  tsimc->SetBranchAddress("H.hod.starttime", &HhodStartTime);                                               
  tsimc->SetBranchAddress("H.hod.fpHitsTime", &HhodfpHitsTime); 
  tsimc->SetBranchAddress("P.gtr.p", &PgtrP); 
  tsimc->SetBranchAddress("H.gtr.p", &HgtrP); 
  tsimc->SetBranchAddress("H.gtr.th", &HgtrTh); 
  tsimc->SetBranchAddress("H.gtr.th", &HgtrPh); 
  tsimc->SetBranchAddress("P.gtr.th", &PgtrTh); 
  tsimc->SetBranchAddress("P.gtr.beta", &pbeta);   
  tsimc->SetBranchAddress("H.gtr.beta", &hbeta); 
  tsimc->SetBranchAddress("H.gtr.dp", &hdelta);   
  tsimc->SetBranchAddress("P.gtr.dp", &pdelta);  
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
  //
  TH1D *h2_epi_PcointimeROC2    = new TH1D("h2_epi_PcointimeROC2","SHMS ROC2 Corrected epi Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_epi_PcointimeROC2_noPID    = new TH1D("h2_epi_PcointimeROC2_noPID","SHMS ROC2 Corrected epi Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ep_PcointimeROC2    = new TH1D("h2_ep_PcointimeROC2","SHMS ROC2 Corrected ep Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ep_PcointimeROC2_noPID    = new TH1D("h2_ep_PcointimeROC2_noPID","SHMS ROC2 Corrected ep Coin Time (no PID cut) ; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ek_PcointimeROC2    = new TH1D("h2_ek_PcointimeROC2","SHMS ROC2 Corrected ek Coin Time; cointime [ns]",       480, -24, 24); 
  TH1D *h2_ek_PcointimeROC2_noPID    = new TH1D("h2_ek_PcointimeROC2_noPID","SHMS ROC2 Corrected ek Coin Time (no PID cut); cointime [ns]",       480, -24, 24); 
   // loop over entries
Long64_t nentries = tsimc->GetEntries();
  Bool_t epievent_cut, epevent_cut, ekevent_cut, positron_cut, event_cut, hpdelta_cut;
  Double_t pOffset = -3.34; //9.5 + 10;  // in ns                                  
  Double_t speedOfLight = 29.9792; // in cm/ns                                   
  Double_t SHMScentralPathLen = 18.1*100;  // SHMS Target to focal plane path length converted to cm  
  Double_t HMScentralPathLen = 22*100;     // HMS Target to focal plane path length converted to cm
  Double_t HMSpartMass = 0.000510998; // electron mass in GeV/c^2
  HgtrX=0;

// nentries=150000;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		//
                hpdelta_cut = hdelta > -10 && hdelta < 10 && pdelta > -10 && pdelta < 20 ;
                Double_t hms_elec = hcernpe > 1.0 && hcaletot > 0.6 && hcaletot < 2.0 && hcalepr > 0.2;


                event_cut = hpdelta_cut &&  hms_elec;

                if (!event_cut) continue;

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
                epevent_cut  =  pcalepr*PgtrP > 0.02 && paeronpe < 0.5 && paeronpe <= 1.0 ;
                ekevent_cut  =  pcalepr*PgtrP < 0.3 && pcalepr*PgtrP > 0.02 && pcernpe <0.5 && paeronpe >= 2.0 ;
                if (PgtrP<=3.0) ekevent_cut  = pcalepr*PgtrP > 0.02 && pcernpe <0.5 && paeronpe <=1 ;
		Double_t toffset=-10.;
                h1_epi_PcointimeROC2_noPID->Fill(epiCT+toffset);
                h1_ep_PcointimeROC2_noPID->Fill(epCT+toffset);
                h1_ek_PcointimeROC2_noPID->Fill(ekCT+toffset);
                if (epievent_cut) h1_epi_PcointimeROC2->Fill(epiCT+toffset);
                if (epevent_cut) h1_ep_PcointimeROC2->Fill(epCT+toffset);
                if (ekevent_cut) h1_ek_PcointimeROC2->Fill(ekCT+toffset);
    Double_t PgtrBetaCalcP = PgtrP/sqrt(PgtrP*PgtrP + ( 0.93827231* 0.93827231));        
    Double_t PgtrBetaCalcK = PgtrP/sqrt(PgtrP*PgtrP + (0.497648*0.497648));        
    Double_t PgtrBetaCalcPi = PgtrP/sqrt(PgtrP*PgtrP + (0.1395704*0.1395704));        
        Double_t DeltaHMSpathLength = 12.462*HgtrTh + 0.1138*HgtrTh*HgtrX - 0.0154*HgtrX - 72.292*HgtrTh*HgtrTh - 0.0000544*HgtrX*HgtrX - 116.52*HgtrPh*HgtrPh;    
        Double_t DeltaSHMSpathLength= +0.11*atan2(PgtrTh,1)*1000;
	Double_t HgtrBetaCalc = HgtrP/sqrt(HgtrP*HgtrP + HMSpartMass*HMSpartMass);          
	Double_t HMScoinCorr = (HMScentralPathLen + DeltaHMSpathLength)/speedOfLight*(1./HgtrBetaCalc-1)  - HhodfpHitsTime;     Double_t SHMScoinCorrP = (SHMScentralPathLen+ DeltaSHMSpathLength)/speedOfLight*( 1./PgtrBetaCalcP-1) - PhodfpHitsTime; 
    Double_t coinP= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrP) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
    Double_t SHMScoinCorrK =   (SHMScentralPathLen+ DeltaSHMSpathLength)/speedOfLight*( 1./PgtrBetaCalcK-1) - PhodfpHitsTime; 
    Double_t coinK= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrK) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
    Double_t SHMScoinCorrPi =  (SHMScentralPathLen+ DeltaSHMSpathLength)/speedOfLight*( 1./PgtrBetaCalcPi-1) - PhodfpHitsTime; 
    Double_t coinPi= (TcoinpTRIG1_ROC2_tdcTimeRaw*0.1 - SHMScoinCorrPi) - (TcoinpTRIG4_ROC2_tdcTimeRaw*0.1 - HMScoinCorr) - pOffset;
	h2_epi_PcointimeROC2_noPID->Fill(coinPi);    
	h2_ep_PcointimeROC2_noPID->Fill(coinP);    
	h2_ek_PcointimeROC2_noPID->Fill(coinK);    
                if (epievent_cut) h2_epi_PcointimeROC2->Fill(coinPi);
                if (epevent_cut) h2_ep_PcointimeROC2->Fill(coinP);
                if (ekevent_cut) h2_ek_PcointimeROC2->Fill(coinK);

	}
		//		
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
