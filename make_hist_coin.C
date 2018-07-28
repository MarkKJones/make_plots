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

void make_hist_coin(TString basename="",Int_t nrun=3288,Double_t pfac=0.983){
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
 Double_t  e_mom;
   tsimc->SetBranchAddress("P.gtr.p",&e_mom);
 Double_t  p_mom;
   tsimc->SetBranchAddress("H.gtr.p",&p_mom);
 Double_t  ptrig6;
   tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime",&ptrig6);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&e_xptar);
 Double_t  h_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&h_yptar);
 Double_t  h_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&h_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&e_xpfp);
 Double_t  emiss;
   tsimc->SetBranchAddress("H.kin.secondary.emiss_nuc",&emiss);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.primary.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.primary.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 300, 0.9,1.5);
    HList.Add(hW);
    TH1F *hms_momdiff = new TH1F("hms_momdiff",Form("Run %d ; P(theta)-P ;Counts",nrun), 300, -.1,.1);
    HList.Add(hms_momdiff);
    TH1F *hWcalc = new TH1F("hWcalc",Form("Run %d Fac = %5.3f; W (GeV);Counts",nrun,pfac),  300, 0.9,1.5);
    HList.Add(hWcalc);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetot);
    TH1F *hThetaDiff = new TH1F("hThetaDiff",Form("Run %d ; Electron Theta- Theta_calc (mr);Counts",nrun), 120, -20,20 );
    HList.Add(hThetaDiff);
    TH2F *hThetaDiffTheta = new TH2F("hThetaDiffTheta",Form("Run %d ; Electron Theta- Theta_calc (mr);Theta",nrun), 120, -20,20,100,-30.,30 );
    HList.Add(hThetaDiffTheta);
    TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV) ;Counts",nrun), 300, -.3,.3);
    HList.Add(hEmiss);
    TH1F *hDeltaDiff = new TH1F("hDeltaDiff",Form("Run %d  Fac = %5.3f; DeltaDiff ;Counts",nrun,pfac), 300, -5,5);
    HList.Add(hDeltaDiff);
    TH2F *hEmissW = new TH2F("hEmissW",Form("Run %d ; Emiss (GeV) ;W (GeV)",nrun), 40, -.02,.02, 75, 0.9,1.2);
    HList.Add(hEmissW);
    TH2F *hEmissXfp = new TH2F("hEmissXfp",Form("Run %d ; Emiss (GeV) ; Xfp ",nrun), 40, -.02,.02, 44,-22,22);
    HList.Add(hEmissXfp);
    TH2F *hEmissYfp = new TH2F("hEmissYfp",Form("Run %d ; Emiss (GeV) ; Yfp ",nrun), 40, -.02,.02, 44,-22,22);
    HList.Add(hEmissYfp);
    TH2F *hEmissXpfp = new TH2F("hEmissXpfp",Form("Run %d ; Emiss (GeV) ; Xpfp ",nrun), 40, -.02,.02, 80,-.06,.06);
    HList.Add(hEmissXpfp);
    TH2F *hEmissYpfp = new TH2F("hEmissYpfp",Form("Run %d ; Emiss (GeV) ; Ypfp ",nrun), 40, -.02,.02, 60,-.03,.03);
    HList.Add(hEmissYpfp);
    TH2F *hDeltaDiffXfp = new TH2F("hDeltaDiffXfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Xfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffYfp = new TH2F("hDeltaDiffYfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Yfp ",nrun,pfac), 300, -1,2, 80,-40,40);
    TH2F *hDeltaDiffXpfp = new TH2F("hDeltaDiffXpfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Xpfp ",nrun,pfac), 300, -1,2, 80,-.05,.05);
    TH2F *hDeltaDiffYpfp = new TH2F("hDeltaDiffYpfp",Form("Run %d  Fac = %5.3f; DeltaDiff (%) ; Ypfp ",nrun,pfac), 300, -1,2, 80,-.02,.02);
    HList.Add(hDeltaDiffXfp);
    HList.Add(hDeltaDiffXpfp);
    HList.Add(hDeltaDiffYfp);
    HList.Add(hDeltaDiffYpfp);
    TH2F *hWXfp = new TH2F("hWXfp",Form("Run %d ; W (GeV) ; Xfp ",nrun), 50, 0.9,1., 44,-22,22);
    TH2F *hWYfp = new TH2F("hWYfp",Form("Run %d ; W (GeV) ; Yfp ",nrun), 50, 0.9,1., 44,-22,22);
    TH2F *hWXpfp = new TH2F("hWXpfp",Form("Run %d ; W (GeV) ; Xpfp ",nrun), 50, 0.9,1., 80,-.06,.06);
    TH2F *hWYpfp = new TH2F("hWYpfp",Form("Run %d ; W (GeV) ; Ypfp ",nrun), 50, 0.9,1., 60,-.03,.03);
    HList.Add(hWXfp);
    HList.Add(hWXpfp);
    HList.Add(hWYfp);
    HList.Add(hWYpfp);
  // loop over entries
    Double_t th_cent=12.2;
  Double_t Mp = .93827;
   Double_t Ei=10.587;
   Double_t p_spec=8.7*pfac;
   if (nrun==3288) th_cent=12.2;
    if (nrun==3371) th_cent=13.935;
    if (nrun==3373) th_cent=9.93;
    if (nrun==3374) th_cent=9.93;
    if (nrun>=3375&&nrun<=3379) th_cent=8.5;
    Double_t htheta_lab;
    if (nrun ==1711) {
      htheta_lab = -53.33;
      th_cent=25.01;
      Ei=2.221;
      p_spec=1.816*pfac;
   }
    if (nrun ==1713) {
      htheta_lab = -47.9;
      th_cent=30.0;
      Ei=2.221;
      p_spec=1.816*pfac;
   }
    if (nrun ==1718) {
      htheta_lab = -47.9;
      th_cent=30.0;
      Ei=2.221;
      p_spec=1.686*pfac;
   }
    if (nrun ==1719) {
      htheta_lab = -47.9;
      th_cent=30.0;
      Ei=2.221;
      p_spec=1.555*pfac;
   }
    if (nrun ==1900) {
      th_cent=27.615;
      Ei=6.430;
      p_spec=3.609*pfac;
   }
    if (nrun ==1854 ||nrun ==1863||nrun ==1864 ) {
      th_cent=27.615;
      Ei=6.430;
      p_spec=3.609*pfac;
   }
   Double_t cos_ts=TMath::Cos(th_cent/180*3.14159);
   Double_t sin_ts=TMath::Sin(th_cent/180*3.14159);
   Double_t cos_hms=TMath::Cos(htheta_lab/180*3.14159);
   Double_t sin_hms=TMath::Sin(htheta_lab/180*3.14159);
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (gindex>-1 && ptrig6 > 50 ) {
		  //		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t theta_hms = TMath::ACos((cos_hms - h_yptar * sin_hms) / TMath::Sqrt( 1. + h_xptar*h_xptar + h_yptar * h_yptar )); // polar 			scattering angle relative to the beam line //rad
		Double_t h_mom_calc= 2*Mp*Ei*(Ei+Mp)*TMath::Cos(theta_hms);
		h_mom_calc=h_mom_calc/(Mp*Mp+2*Mp*Ei+Ei*Ei*TMath::Sin(theta_hms)*TMath::Sin(theta_hms));
		hms_momdiff->Fill(h_mom_calc-p_mom);
		Double_t e_calc = Ei+Mp-TMath::Sqrt(p_mom*p_mom+Mp*Mp);
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
		hW->Fill(W);
		Double_t Enew=p_spec*(1+e_delta/100.);
		Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		Double_t theta_calc = TMath::ACos((p_mom*p_mom-e_mom*e_mom-Ei*Ei)/(-2*Ei*e_mom));
		hThetaDiff->Fill((ThScat-theta_calc)*1000);
		hThetaDiffTheta->Fill((ThScat-theta_calc)*1000,1000*(ThScat-th_cent*3.1459/180));
		hWcalc->Fill(W_calc);
                Double_t Em=Ei+Mp-e_mom-TMath::Sqrt(p_mom*p_mom+Mp*Mp);
		hEmiss->Fill(Em);
		hEmissW->Fill(Em,W);
		hEmissXfp->Fill(Em,e_xfp);
		hEmissYfp->Fill(Em,e_yfp);
		hEmissXpfp->Fill(Em,e_xpfp);
		hEmissYpfp->Fill(Em,e_ypfp);
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
