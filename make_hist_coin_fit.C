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

void make_hist_coin(TString basename="",Int_t nrun=2043){
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
 const static Int_t npar=7;
 const static Int_t nparold=7;
 const static Int_t nfit_max=100000;
 Double_t etemp; 
Int_t nfit=0;
  TVectorD b_delta(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
  //
  // Double_t coeff[npar]={2.25778,54.7339,244.862,-17.8256,15.8021,172.152,1282.58,-466.21,-5878.52,5929.85};
  //Double_t coeff[npar]=(1.02696,58.2212,258.529,341.727,-593.522,589.651,-17.9675,19.717,110.911,1277.94,-181.51,-6090.1,5924.33};
  Double_t coeff[nparold]={1677.59,20.3899,42.4909,2400.09,-18.1085,-254.928,-1821.43};
  //Double_t coeff[npar]={-10.0326,44.1946,69.2675,-957.418,-146.041,60.8183,-2.56417,5.80568,238.447,305.949,-217.468,3243.73,344.064};

  //   Int_t xfpexpon[npar] ={0,0,1,1,0,0,1,0,2,0,1,1,2};
  // Int_t xpfpexpon[npar]={0,0,0,0,1,1,0,1,0,2,1,2,1};
  // Int_t yfpexpon[npar] ={1,0,1,0,1,0,0,0,0,0,0,0,0};
  // Int_t ypfpexpon[npar]={0,1,0,1,0,1,0,0,0,0,0,0,0};
  Int_t xfpexponfit[npar] ={0,0,0,0,1,1,1};
  Int_t xpfpexponfit[npar]={2,1,0,0,0,1,2};
  Int_t yfpexponfit[npar] ={0,0,0,0,0,0,0};
  Int_t ypfpexponfit[npar]={0,0,1,2,0,0,0};
  Int_t xfpexpon[npar] ={0,0,0,0,1,1,1};
  Int_t xpfpexpon[npar]={2,1,0,0,0,1,2};
  Int_t yfpexpon[npar] ={0,0,0,0,0,0,0};
  Int_t ypfpexpon[npar]={0,0,1,2,0,0,0};

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
 Double_t  emiss;
   tsimc->SetBranchAddress("H.kin.secondary.emiss_nuc",&emiss);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.primary.W",&W);
 Double_t  Qsq;
   tsimc->SetBranchAddress("P.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("P.kin.primary.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 250, 0.7,1.2);
    TH1F *hW_calc = new TH1F("hWcalc",Form("Run %d ; W_calc (GeV);Counts",nrun), 250, 0.7,1.2);
    HList.Add(hW);
    TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV) ;Counts",nrun), 300, -.3,.3);
    TH1F *hDeltaDiff = new TH1F("hDeltaDiff",Form("Run %d ; DeltaDiff ;Counts",nrun), 300, -1,2);
    TH1F *hDeltaCorrDiff = new TH1F("hDeltaCorrDiff",Form("Run %d ; DeltaDiff ;Counts",nrun), 300, -1,2);
    TH2F *hEmissW = new TH2F("hEmissW",Form("Run %d ; Emiss (GeV) ;W (GeV)",nrun), 300, -.3,.3, 250, 0.7,1.2);
    TH2F *hEmissXfp = new TH2F("hEmissXfp",Form("Run %d ; Emiss (GeV) ; Xfp ",nrun), 300, -.3,.3, 80,-10,10);
    HList.Add(hEmissXfp);
    TH2F *hEmissYfp = new TH2F("hEmissYfp",Form("Run %d ; Emiss (GeV) ; Yfp ",nrun), 300, -.3,.3, 80,-10,10);
    HList.Add(hEmissYfp);
    TH2F *hEmissXpfp = new TH2F("hEmissXpfp",Form("Run %d ; Emiss (GeV) ; Xpfp ",nrun), 300, -.3,.3, 80,-.05,.05);
    HList.Add(hEmissXpfp);
    TH2F *hEmissYpfp = new TH2F("hEmissYpfp",Form("Run %d ; Emiss (GeV) ; Ypfp ",nrun), 300, -.3,.3, 80,-.02,.02);
    HList.Add(hEmissYpfp);
    TH2F *hDeltaDiffXfp = new TH2F("hDeltaDiffXfp",Form("Run %d ; DeltaDiff (%) ; Xfp ",nrun), 300, -1,2, 80,-10,10);
    TH2F *hDeltaDiffYfp = new TH2F("hDeltaDiffYfp",Form("Run %d ; DeltaDiff (%) ; Yfp ",nrun), 300, -1,2, 80,-10,10);
    TH2F *hDeltaDiffXpfp = new TH2F("hDeltaDiffXpfp",Form("Run %d ; DeltaDiff (%) ; Xpfp ",nrun), 300, -1,2, 80,-.05,.05);
    TH2F *hDeltaDiffYpfp = new TH2F("hDeltaDiffYpfp",Form("Run %d ; DeltaDiff (%) ; Ypfp ",nrun), 300, -1,2, 80,-.02,.02);
    TH2F *hDeltaCOrrDiffXfp = new TH2F("hDeltaCOrrDiffXfp",Form("Run %d ; DeltaCOrrDiff (%) ; Xfp ",nrun), 300, -1,2, 80,-10,10);
    TH2F *hDeltaCOrrDiffYfp = new TH2F("hDeltaCOrrDiffYfp",Form("Run %d ; DeltaCOrrDiff (%) ; Yfp ",nrun), 300, -1,2, 80,-10,10);
    TH2F *hDeltaCOrrDiffXpfp = new TH2F("hDeltaCOrrDiffXpfp",Form("Run %d ; DeltaCOrrDiff (%) ; Xpfp ",nrun), 300, -1,2, 80,-.05,.05);
    TH2F *hDeltaCOrrDiffYpfp = new TH2F("hDeltaCOrrDiffYpfp",Form("Run %d ; DeltaCOrrDiff (%) ; Ypfp ",nrun), 300, -1,2, 80,-.02,.02);
    TH2F *hWXfp = new TH2F("hWXfp",Form("Run %d ; W (GeV) ; Xfp ",nrun), 300, 0.7,1.2, 80,-10,10);
    TH2F *hWYfp = new TH2F("hWYfp",Form("Run %d ; W (GeV) ; Yfp ",nrun), 300, 0.7,1.2, 80,-10,10);
    TH2F *hWXpfp = new TH2F("hWXpfp",Form("Run %d ; W (GeV) ; Xpfp ",nrun), 300, 0.7,1.2, 80,-.05,.05);
    TH2F *hWYpfp = new TH2F("hWYpfp",Form("Run %d ; W (GeV) ; Ypfp ",nrun), 300, 0.7,1.2, 80,-.02,.02);
  // loop over entries
   Double_t cos_ts=TMath::Cos(12.2/180*3.14159);
   Double_t sin_ts=TMath::Sin(12.2/180*3.14159);
   Double_t Mp = .93827;
   Double_t Ei=10.6;
   Double_t p_spec=8.7*.985;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (gindex>-1) {
		Double_t e_calc = Mp*Ei/(Mp + 2*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.));
		Double_t delta_diff = 100*(e_calc - p_spec)/p_spec - e_delta;
                Double_t delta_corr = 0.;
                  for( Int_t icoeff_fit=0; icoeff_fit<nparold; icoeff_fit++ ){
        	     delta_corr+= coeff[icoeff_fit]*
	               pow( e_xfp / 100.0, xfpexpon[icoeff_fit] ) * 
	               pow( e_xpfp, xpfpexpon[icoeff_fit] ) *
		       pow( e_yfp/100., yfpexpon[icoeff_fit] ) * 
	               pow( e_ypfp, ypfpexpon[icoeff_fit] ) ;
		  }
		  Double_t Enew = p_spec*(1.+(e_delta-delta_corr)/100.);
		  Double_t W_calc= TMath::Sqrt(-4*Enew*Ei*TMath::Sin(ThScat/2.)*TMath::Sin(ThScat/2.)+Mp*Mp+2*Mp*(Ei-Enew));
		hW->Fill(W);
		hW_calc->Fill(W_calc);
		hEmiss->Fill(emiss);
		hEmissW->Fill(emiss,W);
		hEmissXfp->Fill(emiss,e_xfp);
		hEmissYfp->Fill(emiss,e_yfp);
		hEmissXpfp->Fill(emiss,e_xpfp);
		hEmissYpfp->Fill(emiss,e_ypfp);
		hDeltaDiff->Fill(delta_diff);
		hDeltaCorrDiff->Fill(delta_diff-delta_corr);
		hDeltaDiffXfp->Fill(delta_diff,e_xfp);
		hDeltaDiffYfp->Fill(delta_diff,e_yfp);
		hDeltaDiffXpfp->Fill(delta_diff,e_xpfp);
		hDeltaDiffYpfp->Fill(delta_diff,e_ypfp);
		hDeltaCOrrDiffXfp->Fill(delta_diff-delta_corr,e_xfp);
		hDeltaCOrrDiffYfp->Fill(delta_diff-delta_corr,e_yfp);
		hDeltaCOrrDiffXpfp->Fill(delta_diff-delta_corr,e_xpfp);
		hDeltaCOrrDiffYpfp->Fill(delta_diff-delta_corr,e_ypfp);
		hWXfp->Fill(W,e_xfp);
		hWYfp->Fill(W,e_yfp);
		hWXpfp->Fill(W,e_xpfp);
		hWYpfp->Fill(W,e_ypfp);
		//
		if ( TMath::Abs(delta_diff-delta_corr) < .5 && nfit<(nfit_max-100) && TMath::Abs(e_xfp) < 20. && TMath::Abs(e_xpfp) < .06 && TMath::Abs(e_yfp) < 20. && TMath::Abs(e_ypfp) < .06 ) {
		  if (nfit%1000==0) cout << " nfit = " << nfit << endl;
                  for( Int_t icoeff_fit=0; icoeff_fit<npar; icoeff_fit++ ){
        	     etemp= 
	               pow( e_xfp / 100.0, xfpexponfit[icoeff_fit] ) * 
	               pow( e_xpfp, xpfpexponfit[icoeff_fit] ) *
		       pow( e_yfp/100., yfpexponfit[icoeff_fit] ) * 
	               pow( e_ypfp, ypfpexponfit[icoeff_fit] ) ;
                         lambda[icoeff_fit][nfit] = etemp;
	                   b_delta[icoeff_fit] += (delta_diff) * etemp;
 	           } // for icoeff_fit loop
	           nfit++;
		}
		}
		//		
	}
	//
	cout << " doing fit = " << nfit << endl;
 for(Int_t i=0; i<npar; i++){
    for(Int_t j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
 }
   for( Int_t ifit=0; ifit<nfit; ifit++){
     if( ifit % 10000 == 0 ) cout << " at number = " << ifit << endl;
    for( Int_t ipar=0; ipar<npar; ipar++){
      for( Int_t jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
       }
    }
  }
   TDecompSVD Ay_svd(Ay);
  bool ok;
  ok = Ay_svd.Solve( b_delta );
  cout << "Delta solution ok = " << ok << endl;
  b_delta.Print();
                  for( Int_t icoeff_fit=0; icoeff_fit<npar; icoeff_fit++ ){
		    cout << b_delta[icoeff_fit] << "," ;
		  }
	//
		  cout << endl;
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
