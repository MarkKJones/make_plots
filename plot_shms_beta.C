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

void plot_shms_beta(TString basename="",Int_t nrun=2043){
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
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 TString hlab[4]={"1x","1y","2x","2y"};
 Int_t hnumpad[4]={13,13,14,21};
 Double_t  hodo_neg_tdc[4][21];
 Double_t  hodo_pos_tdc[4][21];
 TString temp;
 for (Int_t i=0;i<4;i++) {
   temp="P.hod."+hlab[i]+".GoodNegTdcTimeUnCorr";
 tsimc->SetBranchAddress(temp,&hodo_neg_tdc[i]) ;
 temp="P.hod."+hlab[i]+".GoodPosTdcTimeUnCorr";
 tsimc->SetBranchAddress(temp,&hodo_pos_tdc[i]) ;
 }
 Double_t beta;
 tsimc->SetBranchAddress("P.hod.beta",&beta) ;
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
 Double_t  track_index;
   tsimc->SetBranchAddress("P.gtr.index",&track_index);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
  //
   TH1F *hetracknorm = new TH1F("hetracknorm", Form(" Run %d ; Track Cal E/p ; Counts",nrun), 200, 0.0,1.5);
   TH2F *hetracknorm_xfp = new TH2F("hetracknorm_xfp", Form(" Run %d ; Track Cal E/p ; Xfp",nrun), 200, 0.0,1.5,200,-45,45);
   TH2F *hetracknorm_yfp = new TH2F("hetracknorm_yfp", Form(" Run %d ; Track Cal E/p ; Yfp",nrun), 200, 0.0,1.5,200,-45,45);
   TH1F *hbeta = new TH1F("hbeta", Form(" Run %d ; Beta ; Counts",nrun), 200, -0.1,1.5);
   TH2F *hbeta_xfp = new TH2F("hbeta_xfp", Form(" Run %d ; Beta ; Xfp",nrun), 200, -0.1,1.5,200,-45,45);
   TH2F *hbeta_yfp = new TH2F("hbeta_yfp", Form(" Run %d ; Beta ; Yfp",nrun), 200, -0.1,1.5,200,-45,45);
   //// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (track_index>-1) {
                  hbeta->Fill(beta);
                  hbeta_xfp->Fill(beta,xfp);
                  hbeta_yfp->Fill(beta,yfp);
		  hetracknorm->Fill(etracknorm);
		  hetracknorm_xfp->Fill(etracknorm,xfp);
		  hetracknorm_yfp->Fill(etracknorm,yfp);
		}
	}
	//
}
