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

void make_hist_sidis_coin(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_sidis_coin_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  evtyp;
   tsimc->SetBranchAddress("g.evtyp",&evtyp);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
 Double_t  pindex;
   tsimc->SetBranchAddress("P.gtr.index",&pindex);
 Double_t  hindex;
   tsimc->SetBranchAddress("H.gtr.index",&hindex);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("P.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&e_delta);
 Double_t  ptrig6;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC2_tdcTime",&ptrig6);
 Double_t  ptrig4;
   tsimc->SetBranchAddress("T.coin.pTRIG4_ROC2_tdcTime",&ptrig4);
 Double_t  ptrig1;
   tsimc->SetBranchAddress("T.coin.pTRIG1_ROC2_tdcTime",&ptrig1);
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
 Double_t  W;
   tsimc->SetBranchAddress("H.kin.primary.W",&W);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm);
 Double_t  npesum;
   tsimc->SetBranchAddress("H.cer.npeSum",&npesum);
 Double_t  Qsq;
   tsimc->SetBranchAddress("H.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("H.kin.primary.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 125, 0.7,2.0);
    HList.Add(hW);
    TH2F *hetotnorm_npe = new TH2F("hetotnorm_npe",Form("Run %d ; Etot norm;Npe",nrun), 120, 0.0,1.2,50,0,25);
    HList.Add(hetotnorm_npe);
    TH1F *hetotnorm = new TH1F("hetotnorm",Form("Run %d ; Etot norm; Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetotnorm);
    TH1F *hetotnorm_coin = new TH1F("hetotnorm_coin",Form("Run %d ; Etot norm; Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetotnorm_coin);
    TH1F *hetotnorm_hms = new TH1F("hetotnorm_hms",Form("Run %d ; Etot norm; Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetotnorm_hms);
    TH2F *hptrig4_ptrig6 = new TH2F("hptrig6_ptrig4",Form("Run %d ; Ptrig 4 ; Ptrig 6",nrun), 250, 0. ,500.,250,0,500);
    HList.Add(hptrig4_ptrig6);
  // loop over entries
    Double_t ch;
    Double_t coin_lt;
    Double_t hms_lt;
    Double_t hms_ps;
    if (nrun==3572){
      ch=4.175;
      coin_lt=53.1341/100.;
      hms_lt=49.7630/100.;
      hms_ps=257;
    }
    if (nrun==3573){
      ch=4.578;
      coin_lt=74.2566/100.;
      hms_lt=71.7261/100.;
      hms_ps=257;
    }
    if (nrun==3650){
      ch=13.125;
      coin_lt=1.0/100.;
      hms_lt=99.5/100.;
      hms_ps=9;
    }
    if (nrun==3648){
      ch=5.840;
      coin_lt=1.0/100.;
      hms_lt=99.5/100.;
      hms_ps=3;
    }
    if (nrun==3492){
      ch=26.382;
      coin_lt=90.98/100.;
      hms_lt=88.19/100.;
      hms_ps=257;
    }
    //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		hetotnorm_npe->Fill(etotnorm,npesum);
		if (etotnorm > 0.8) hetotnorm->Fill(etotnorm);
		hptrig4_ptrig6->Fill(ptrig4,ptrig6);
		if (etotnorm > 0.8 && evtyp>=4&&evtyp<=7) hetotnorm_coin->Fill(etotnorm);
		if (etotnorm > 0.8 && evtyp==2) hetotnorm_hms->Fill(etotnorm);
		if (pindex>-1&&hindex>-1 ) {
		  hW->Fill(W);
		//
		}
		//		
	}
	//
	cout << " Run = " << nrun << endl;
	cout << " Integral coin = " << hetotnorm_coin->Integral() << endl;
	cout << " Integral hms = " << hetotnorm_hms->Integral() << endl;
	cout << " Yield coin = " << hetotnorm_coin->Integral()/coin_lt/ch << endl;
	cout << " Yield hms = " << hetotnorm_hms->Integral()/hms_lt/ch*hms_ps << endl;
	cout << " Yield  = " << hetotnorm_coin->Integral()/coin_lt/ch+hetotnorm_hms->Integral()/hms_lt/ch*hms_ps << endl;
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
