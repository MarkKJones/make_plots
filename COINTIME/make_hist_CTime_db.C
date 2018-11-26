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

void make_hist_CTime_db(TString basename="",Int_t nrun=1267,Double_t cpeak=45){
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
   outputhist= "hist/"+basename+"_CTime_db.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Int_t hhodoref1_mult;
 Double_t hhodoref1_tdcTimeRaw[10];
 tsimc->SetBranchAddress("Ndata.D.hhodoref1",&hhodoref1_mult);
 tsimc->SetBranchAddress("D.hhodoref1",&hhodoref1_tdcTimeRaw);
 Int_t hstof_mult;
 Double_t hstof_tdcTimeRaw[10];
 tsimc->SetBranchAddress("Ndata.D.hstof",&hstof_mult);
 tsimc->SetBranchAddress("D.hstof",&hstof_tdcTimeRaw);
 Int_t htrigref1_mult;
 Double_t htrigref1_tdcTimeRaw[10];
 tsimc->SetBranchAddress("Ndata.D.htrigref1",&htrigref1_mult);
 tsimc->SetBranchAddress("D.htrigref1",&htrigref1_tdcTimeRaw);
 Int_t helclean_mult;
 Double_t helclean_tdcTimeRaw[10];
 tsimc->SetBranchAddress("Ndata.D.helclean",&helclean_mult);
 tsimc->SetBranchAddress("D.helclean",&helclean_tdcTimeRaw);
 Int_t helreal_mult;
 Double_t helreal_tdcTimeRaw[10];
 tsimc->SetBranchAddress("Ndata.D.helreal",&helreal_mult);
 tsimc->SetBranchAddress("D.helreal",&helreal_tdcTimeRaw);
  Double_t  pTRIG6_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime",&pTRIG6_ROC1_tdcTime);
 Double_t  pTRIG6_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTime",&pTRIG6_ROC2_tdcTime);
   Double_t  pTRIG6_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTimeRaw",&pTRIG6_ROC1_tdcTimeRaw);
 Double_t  pTRIG6_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTimeRaw",&pTRIG6_ROC2_tdcTimeRaw);
 Double_t  pTRIG4_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTime",&pTRIG4_ROC1_tdcTime);
 Double_t  pTRIG4_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTimeRaw",&pTRIG4_ROC2_tdcTimeRaw);
 Double_t  pTRIG4_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC1_tdcTimeRaw",&pTRIG4_ROC1_tdcTimeRaw);
 Double_t  pTRIG4_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime",&pTRIG4_ROC2_tdcTime);
 Double_t  pTRIG5_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcTime",&pTRIG5_ROC1_tdcTime);
 Double_t  pTRIG5_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTimeRaw",&pTRIG5_ROC2_tdcTimeRaw);
 Double_t  pTRIG5_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG5_ROC1_tdcTimeRaw",&pTRIG5_ROC1_tdcTimeRaw);
 Double_t  pTRIG5_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG5_ROC2_tdcTime",&pTRIG5_ROC2_tdcTime);
 Double_t  hTRIG1_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.hTRIG1_ROC2_tdcTime",&hTRIG1_ROC2_tdcTime);
 Double_t  pTRIG3_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG3_ROC1_tdcTime",&pTRIG3_ROC1_tdcTime);
 Double_t  pTRIG3_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG3_ROC2_tdcTime",&pTRIG3_ROC2_tdcTime);
 Double_t  pTRIG1_ROC1_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTime",&pTRIG1_ROC1_tdcTime);
 Double_t  pTRIG1_ROC2_tdcTime;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime",&pTRIG1_ROC2_tdcTime);
 Double_t  pTRIG1_ROC1_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC1_tdcTimeRaw",&pTRIG1_ROC1_tdcTimeRaw);
 Double_t  pTRIG1_ROC2_tdcTimeRaw;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTimeRaw",&pTRIG1_ROC2_tdcTimeRaw);

			 //
    TH1F *h_hodoref_mult= new TH1F("h_hodoref_mult",Form("Run %d ;HMS Hodo reftime mult ;Counts",nrun), 10,0.,10.);
    TH1F *h_hstof_mult= new TH1F("h_hstof_mult",Form("Run %d ;HMS STOF mult ;Counts",nrun), 10,0.,10.);
    TH1F *h_htrigref1_mult= new TH1F("h_htrigref1_mult",Form("Run %d ;HMS TRIGREF1 mult ;Counts",nrun), 10,0.,10.);
    TH1F *h_helclean_mult= new TH1F("h_helclean_mult",Form("Run %d ;HMS ELCLEAN mult ;Counts",nrun), 10,0.,10.);
    TH1F *h_helreal_mult= new TH1F("h_helreal_mult",Form("Run %d ;HMS ELREAL mult ;Counts",nrun), 10,0.,10.);
    HList.Add(h_hodoref_mult);
    HList.Add(h_hstof_mult);
			 TH1F *h_hhodo_tdcraw[4];
			 TH1F *h_htrigref1_tdcraw[4];
			 TH1F *h_helclean_tdcraw[4];
			 TH1F *h_helreal_tdcraw[4];
			 TH1F *h_hstof_tdcraw[4];
			 for (Int_t ih=0;ih<4;ih++)         {
			   h_hstof_tdcraw[ih] = new TH1F(Form("h_hstof_tdcraw_%d",ih),Form("HMS stof hit = %d; STOF raw time hit=%d;Counts",ih,ih),1000,0,10000);
			   h_hhodo_tdcraw[ih] = new TH1F(Form("h_hhodo_tdcraw_%d",ih),Form("HMS Hodo ref hit = %d; Hodo ref raw time hit=%d;Counts",ih,ih),1000,0,10000);
			   h_htrigref1_tdcraw[ih] = new TH1F(Form("h_htrigref1_tdcraw_%d",ih),Form("HMS Trig ref hit = %d; Trig ref raw time hit=%d;Counts",ih,ih),1000,0,10000);
			   h_helclean_tdcraw[ih] = new TH1F(Form("h_helclean_tdcraw_%d",ih),Form("HMS Trig ref hit = %d; HMS ELCLEAN raw time hit=%d;Counts",ih,ih),1000,0,10000);
			   h_helreal_tdcraw[ih] = new TH1F(Form("h_helreal_tdcraw_%d",ih),Form("HMS Trig ref hit = %d HMS ELREAL raw time hit=%d;Counts",ih,ih),1000,0,10000);
			 }
			 h_hhodo_tdcraw_mult1_trig4 = new TH1F(Form("h_hhodo_tdcraw_mult1_trig4_%d",1),Form("HMS Hodo ref hit = %d mult = 1 Trig 4; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_hhodo_tdcraw_mult1_trig6 = new TH1F(Form("h_hhodo_tdcraw_mult1_trig6_%d",1),Form("HMS Hodo ref hit = %d mult = 1 Trig 6; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_hhodo_tdcraw_mult2_trig5_hit1 = new TH1F(Form("h_hhodo_tdcraw_mult2_trig5_%d",1),Form("HMS Hodo ref hit = %d mult = 2 TRIG 5; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_hhodo_tdcraw_mult2_trig5_hit2 = new TH1F(Form("h_hhodo_tdcraw_mult2_trig5_%d",2),Form("HMS Hodo ref hit = %d mult = 2 TRIG 5; Hodo ref raw time hit=%d;Counts",2,2),1000,0,10000);
      			 h_hhodo_tdcraw_mult2_trig6_hit1 = new TH1F(Form("h_hhodo_tdcraw_mult2_trig6_%d",1),Form("HMS Hodo ref hit = %d mult = 2 TRIG 6; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_hhodo_tdcraw_mult2_trig6_hit2 = new TH1F(Form("h_hhodo_tdcraw_mult2_trig6_%d",2),Form("HMS Hodo ref hit = %d mult = 2 TRIG 6; Hodo ref raw time hit=%d;Counts",2,2),1000,0,10000);
			 h_hhodo_tdcraw_mult2_elclean = new TH1F(Form("h_hhodo_tdcraw_elclean_mult2_hit_%d",1),Form("Elclean HMS Hodo ref hit = %d mult = 2; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_hhodo_tdcraw_mult2_elreal = new TH1F(Form("h_hhodo_tdcraw_elreal_mult2_hit_%d",1),Form("Elreal HMS Hodo ref hit = %d mult = 2; Hodo ref raw time hit=%d;Counts",1,1),1000,0,10000);
			 h_htrigref1_tdcraw_mult1 = new TH1F(Form("h_htrigref1_tdcraw_mult_%d",1),Form("HMS Trig ref hit = %d mult = 1; Trig 1 ref raw time hit=%d;Counts",1,1),1000,0,10000);
    Long64_t nentries = tsimc->GetEntries();
	for (int ie = 0; ie < nentries; ie++) {
      		tsimc->GetEntry(ie);
                if (ie%50000==0) cout << " Entry = " << ie << endl;
		h_hodoref_mult->Fill(hhodoref1_mult);
		h_htrigref1_mult->Fill(htrigref1_mult);
 		h_hstof_mult->Fill(hstof_mult);
 		h_helreal_mult->Fill(helreal_mult);
 		h_helclean_mult->Fill(helclean_mult);
                if (hhodoref1_mult==1&&pTRIG4_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult1_trig4->Fill(hhodoref1_tdcTimeRaw[0]);
                if (hhodoref1_mult==1&&pTRIG6_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult1_trig6->Fill(hhodoref1_tdcTimeRaw[0]);
                if (hhodoref1_mult==2&&pTRIG5_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult2_trig5_hit1->Fill(hhodoref1_tdcTimeRaw[0]);
                if (hhodoref1_mult==2&&pTRIG5_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult2_trig5_hit2->Fill(hhodoref1_tdcTimeRaw[1]);
                if (hhodoref1_mult==2&&pTRIG6_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult2_trig6_hit1->Fill(hhodoref1_tdcTimeRaw[0]);
                if (hhodoref1_mult==2&&pTRIG6_ROC2_tdcTimeRaw>0) h_hhodo_tdcraw_mult2_trig6_hit2->Fill(hhodoref1_tdcTimeRaw[1]);
                if (hhodoref1_mult>=1&&helreal_mult==1) h_hhodo_tdcraw_mult2_elreal->Fill(hhodoref1_tdcTimeRaw[0]);
                if (hhodoref1_mult>=1&&helclean_mult==1) h_hhodo_tdcraw_mult2_elclean->Fill(hhodoref1_tdcTimeRaw[0]);
                if (htrigref1_mult==1) h_htrigref1_tdcraw_mult1->Fill(htrigref1_tdcTimeRaw[0]);
                for (Int_t ih=0;ih<htrigref1_mult;ih++){
		  h_htrigref1_tdcraw[ih]->Fill(htrigref1_tdcTimeRaw[ih]);
		}
                for (Int_t ih=0;ih<hhodoref1_mult;ih++){
		  h_hhodo_tdcraw[ih]->Fill(hhodoref1_tdcTimeRaw[ih]);
		}
                for (Int_t ih=0;ih<hstof_mult;ih++){
		  h_hstof_tdcraw[ih]->Fill(hstof_tdcTimeRaw[ih]);
		}
                for (Int_t ih=0;ih<helclean_mult;ih++){
		  h_helclean_tdcraw[ih]->Fill(helclean_tdcTimeRaw[ih]);
		}
                for (Int_t ih=0;ih<helreal_mult;ih++){
		  h_helreal_tdcraw[ih]->Fill(helreal_tdcTimeRaw[ih]);
		}
	}
 TFile hsimc(outputhist,"recreate");
	HList.Write();
			 //
			 }
