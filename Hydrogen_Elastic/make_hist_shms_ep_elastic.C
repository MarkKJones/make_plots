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

void make_hist_shms_ep_elastic(TString basename="",Int_t nrun=3288,Double_t pfac=0.981){
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
   outputhist= "hist/"+basename+"_shms_ep_elastic_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  e_ytar;
   tsimc->SetBranchAddress("P.gtr.y",&e_ytar);
 Double_t  gindex;
   tsimc->SetBranchAddress("P.gtr.index",&gindex);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("P.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&e_delta);
 Double_t  ptrig6;
   tsimc->SetBranchAddress("T.coin.pTRIG6_tdctime",&ptrig6);
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
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 300, 0.9,1.5);
    HList.Add(hW);
    TH1F *hWcalc = new TH1F("hWcalc",Form("Run %d Fac = %5.3f; W (GeV);Counts",nrun,pfac),  300, 0.9,1.5);
    HList.Add(hWcalc);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetot);
    TH1F *hDeltaDiff = new TH1F("hDeltaDiff",Form("Run %d  Fac = %5.3f; DeltaDiff ;Counts",nrun,pfac), 300, -5,5);
    HList.Add(hDeltaDiff);
    TH2F *hDeltaDiffXfp = new TH2F("hDeltaDiffXfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Xfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffYfp = new TH2F("hDeltaDiffYfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Yfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffXpfp = new TH2F("hDeltaDiffXpfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Xpfp ",nrun,pfac), 300, -1,2, 80,-.05,.05);
    TH2F *hDeltaDiffYpfp = new TH2F("hDeltaDiffYpfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Ypfp ",nrun,pfac), 300, -1,2, 80,-.02,.02);
    HList.Add(hDeltaDiffXfp);
    HList.Add(hDeltaDiffXpfp);
    HList.Add(hDeltaDiffYfp);
    HList.Add(hDeltaDiffYpfp);
    TH2F *hWXfp = new TH2F("hWXfp",Form("Run %d ; W (GeV) ; Xfp ",nrun), 100, 0.8,1.2, 80,-40,40);
    TH2F *hWYfp = new TH2F("hWYfp",Form("Run %d ; W (GeV) ; Yfp ",nrun), 100, 0.8,1.2, 80,-40,40);
    TH2F *hWXpfp = new TH2F("hWXpfp",Form("Run %d ; W (GeV) ; Xpfp ",nrun), 100, 0.8,1.2, 80,-.05,.05);
    TH2F *hWYpfp = new TH2F("hWYpfp",Form("Run %d ; W (GeV) ; Ypfp ",nrun), 100, 0.8,1.2, 80,-.02,.02);
  // loop over entries
    Double_t th_cent=7.5;
  Double_t Mp = .93827;
   Double_t Ei=2.221;
   Double_t p_spec=2.177*pfac;
   Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (gindex>-1 ) {
		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
		hW->Fill(W);
		Double_t Enew=p_spec*(1+e_delta/100.);
		Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		hWcalc->Fill(W_calc);
		if (W<1.075) {
		hDeltaDiff->Fill(delta_diff);
		hDeltaDiffXfp->Fill(delta_diff,e_xfp);
		hDeltaDiffYfp->Fill(delta_diff,e_yfp);
		hDeltaDiffXpfp->Fill(delta_diff,e_xpfp);
		hDeltaDiffYpfp->Fill(delta_diff,e_ypfp);
		hWXfp->Fill(W,e_xfp);
		hWYfp->Fill(W,e_yfp);
		hWXpfp->Fill(W,e_xpfp);
		hWYpfp->Fill(W,e_ypfp);
		}
		//
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
