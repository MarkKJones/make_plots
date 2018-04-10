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

void make_hist_hms_sieve(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hms_sieve_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etracknorm",&etracknorm);
 Double_t  ytar;
   tsimc->SetBranchAddress("H.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("H.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("H.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("H.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("H.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("H.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("H.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp);

   // Define histograms
   TH1F *hztar[5];
   TH1F *hyptar[5];
   TH2F *hytar_yptar;
   TH2F *hztar_yptar;
   cout << " nrun = " << nrun << endl;
   Double_t yp_cutlo[5]={-.025,-0.0125,-0.005,0.005,0.0125};
   Double_t yp_cuthi[5]={-0.0125,-0.005,0.005,0.0125,0.025};
     cout << " make histo " << endl;
   hytar_yptar = new TH2F("hytar_yptar", Form("Run %d ; Y_tar ; Yp_tar",nrun), 80,-1.,1.,60,-.035,.035);
   HList.Add(hytar_yptar);
   hztar_yptar = new TH2F("hztar_yptar", Form("Run %d ; Yp_tar ; Z_tar",nrun), 60,-.035,.035,60,-3.,3.);
   HList.Add(hztar_yptar);
	for (Int_t iz = 0; iz < 5; iz++) {
	  hztar[iz] = new TH1F(Form("hztar_%d",iz), Form("Run %d ; Z_tar %5.3f < YP < %5.3f ; Counts",nrun,yp_cutlo[iz],yp_cuthi[iz]), 80,-4.,4.);	  
          HList.Add(hztar[iz]);
	  hyptar[iz] = new TH1F(Form("hyptar_%d",iz), Form("Run %d ; Yp_tar %5.3f < YP < %5.3f ; Counts",nrun,yp_cutlo[iz],yp_cuthi[iz]), 80,-.03,.03);	  
          HList.Add(hyptar[iz]);
 	}
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (sumnpe>2 && etracknorm>.6) {
		  hztar_yptar->Fill(yptar,reactz);
		hytar_yptar->Fill(ytar,yptar);
	        for (int iz = 0; iz < 5; iz++) {
		  if (yptar>yp_cutlo[iz] && yptar <=yp_cuthi[iz]) hztar[iz]->Fill(reactz);
		  if (yptar>yp_cutlo[iz] && yptar <=yp_cuthi[iz]&&reactz<4.&&reactz>-4.) hyptar[iz]->Fill(yptar);
		}
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
