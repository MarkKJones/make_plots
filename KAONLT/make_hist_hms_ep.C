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

void make_hist_hms_ep(TString basename="",Int_t nrun=3288,Double_t pfac=0.981){
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
   outputhist= "hist/"+basename+"_hms_ep_elastic_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  etottracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etottracknorm);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("H.gtr.y",&e_ytar);
 Double_t  gindex;
   tsimc->SetBranchAddress("H.gtr.index",&gindex);
 Double_t  gpindex;
   tsimc->SetBranchAddress("P.gtr.index",&gpindex);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("H.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("H.gtr.dp",&e_delta);
 Double_t  p_delta;
   tsimc->SetBranchAddress("P.gtr.dp",&p_delta);
 Double_t  e_mom;
   tsimc->SetBranchAddress("H.gtr.p",&e_mom);
 Double_t  p_mom;
   tsimc->SetBranchAddress("P.gtr.p",&p_mom);
   // Double_t  ptrig6;
   //tsimc->SetBranchAddress("T.coin.pTRIG6_ROC1_tdcTime",&ptrig6);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&e_xptar);
 Double_t  p_yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&p_yptar);
 Double_t  p_xptar;
   tsimc->SetBranchAddress("P.gtr.th",&p_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&e_xpfp);
 Double_t  p_yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&p_yfp);
 Double_t  p_ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&p_ypfp);
   Double_t  p_xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&p_xfp);
 Double_t  p_xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&p_xpfp);
 Double_t  W;
   tsimc->SetBranchAddress("H.kin.primary.W",&W);
 Double_t  emiss;
   tsimc->SetBranchAddress("P.kin.secondary.emiss",&emiss);
 Double_t  pmiss;
   tsimc->SetBranchAddress("P.kin.secondary.pmiss",&pmiss);
 Double_t  pmissx;
   tsimc->SetBranchAddress("P.kin.secondary.pmiss_x",&pmissx);
 Double_t  pmissy;
   tsimc->SetBranchAddress("P.kin.secondary.pmiss_y",&pmissy);
 Double_t  pmissz;
   tsimc->SetBranchAddress("P.kin.secondary.pmiss_z",&pmissz);
 Double_t  Qsq;
   tsimc->SetBranchAddress("H.kin.primary.Q2",&Qsq);
 Double_t  ThScat;
   tsimc->SetBranchAddress("H.kin.primary.scat_ang_rad",&ThScat);
   // Define histograms
    TH1F *hW = new TH1F("hW",Form("Run %d ; W (GeV);Counts",nrun), 200, 0.9,3.0);
    HList.Add(hW);
     TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; HMS X_fp;Counts",nrun), 100, -40.,40.);
    HList.Add(hxfp);
    TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; HMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(hyfp);
    TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; HMS Xp_fp;Counts",nrun), 100, -.1,.1);
    HList.Add(hxpfp);
    TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; HMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(hypfp);
    TH1F *pxfp = new TH1F("pxfp",Form("Run %d ; SHMS X_fp;Counts",nrun), 100, -50.,50.);
    HList.Add(pxfp);
    TH1F *pyfp = new TH1F("pyfp",Form("Run %d ; SHMS Y_fp;Counts",nrun), 100, -20.,20.);
    HList.Add(pyfp);
    TH1F *pxpfp = new TH1F("pxpfp",Form("Run %d ; SHMS Xp_fp;Counts",nrun), 100, -.12,.12);
    HList.Add(pxpfp);
    TH1F *pypfp = new TH1F("pypfp",Form("Run %d ; SHMS Yp_fp;Counts",nrun), 100, -.05,.05);
    HList.Add(pypfp);
    TH1F *hxptar = new TH1F("hxptar",Form("Run %d ; HMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hxptar);
    TH1F *hyptar = new TH1F("hyptar",Form("Run %d ; HMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(hyptar);
    TH1F *pxptar = new TH1F("pxptar",Form("Run %d ; SHMS Xp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pxptar);
    TH1F *pyptar = new TH1F("pyptar",Form("Run %d ; SHMS Yp_tar;Counts",nrun), 100, -.1,.1);
    HList.Add(pyptar);
    TH1F *pdelta = new TH1F("pdelta",Form("Run %d ; SHMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(pdelta);
    TH1F *hdelta = new TH1F("hdelta",Form("Run %d ; HMS Delta;Counts",nrun), 100,-10.,10.);
    HList.Add(hdelta);
    TH1F *hEmiss = new TH1F("hEmiss",Form("Run %d ; Emiss (GeV);Counts",nrun), 200, -.05,.1);
    HList.Add(hEmiss);
    TH1F *hPmiss = new TH1F("hPmiss",Form("Run %d ; Pmiss (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmiss);
    TH1F *hMmiss2 = new TH1F("hMmiss2",Form("Run %d ; Mmiss2 (GeV2);Counts",nrun), 400, 0.,1.5);
    HList.Add(hMmiss2);
    TH1F *hPmissx = new TH1F("hPmissx",Form("Run %d ; Pmissx (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissx);
    TH1F *hPmissy = new TH1F("hPmissy",Form("Run %d ; Pmissy (GeV);Counts",nrun), 200, -.1,.1);
     HList.Add(hPmissy);
     TH1F *hPmissz = new TH1F("hPmissz",Form("Run %d ; Pmissz (GeV);Counts",nrun), 200, -.1,.1);
    HList.Add(hPmissz);
     TH1F *hetottracknorm = new TH1F("hetot",Form("Run %d ; Etottrack norm;Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetottracknorm);
 Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (gindex>-1 && gpindex>-1 ) {
		  hetottracknorm->Fill(etottracknorm);
		  //		  if (etottracknorm > 0.8 ) {
		  		  if (1==1 ) {
                  if (	TMath::Abs(e_delta) < 8 && p_delta>-10 && p_delta<22) {
		    hW->Fill(W);
		  hxfp->Fill(e_xfp);		  
		  hyfp->Fill(e_yfp);		  
		  hxpfp->Fill(e_xpfp);		  
		  hypfp->Fill(e_ypfp);		  
		  pxfp->Fill(p_xfp);		  
		  pyfp->Fill(p_yfp);		  
		  pxpfp->Fill(p_xpfp);		  
		  pypfp->Fill(p_ypfp);		  
		  hxptar->Fill(e_xptar);		  
		  pxptar->Fill(p_xptar);		  
		  hyptar->Fill(e_yptar);		  
		  pyptar->Fill(p_yptar);		  
		  hdelta->Fill(e_delta);		  
		  pdelta->Fill(p_delta);		  
		  hEmiss->Fill(emiss);		  
		  hPmiss->Fill(pmiss);		  
		  hMmiss2->Fill(emiss*emiss-pmiss*pmiss);	
		  hPmissx->Fill(pmissx);		  
		  hPmissy->Fill(pmissy);		  
		  hPmissz->Fill(pmissz);
		  }		  
				  } // etot cut
		//
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
