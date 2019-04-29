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

void make_hist_shms_cal(TString basename="",Int_t nrun=1267,Double_t cpeak=45){
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
   outputhist= "hist/"+basename+"_shms_cal_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t nblk=224;
 Double_t adcTdcDiffTime[nblk];
 Double_t adcMult[nblk];
 tsimc->SetBranchAddress("P.cal.fly.goodAdcMult",adcMult);
 tsimc->SetBranchAddress("P.cal.fly.goodAdcTdcDiffTime",adcTdcDiffTime);
 //
 TH1F *hAdcTdcDiffTime[nblk];
 for (Int_t ip=0;ip<nblk;ip++) {
   hAdcTdcDiffTime[ip] = new TH1F(Form("pcal_%d",ip),Form(" Block %d ; Adctdcdifftime (ns); counts",ip),100,-100,100);
   HList.Add(hAdcTdcDiffTime[ip]);
 }
 //
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		 for (Int_t ip=0;ip<nblk;ip++) {
		   if (adcMult[ip]==1) hAdcTdcDiffTime[ip]->Fill(adcTdcDiffTime[ip]);
		 }
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
 //
}
