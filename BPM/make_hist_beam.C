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

void make_hist_beam(TString basename="",Int_t nrun=3288){
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
   outputhist= "hist/"+basename+"_beam_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  beamx;
  tsimc->SetBranchAddress("H.rb.raster.fr_xcent",&beamx);
 Double_t  beamy;
   tsimc->SetBranchAddress("H.rb.raster.fr_ycent",&beamy);
 Double_t  rastx;
   tsimc->SetBranchAddress("H.rb.raster.fr_xa",&rastx);
 Double_t  rasty;
   tsimc->SetBranchAddress("H.rb.raster.fr_ya",&rasty);
 Double_t  e_ytar;
   tsimc->SetBranchAddress("H.gtr.y",&e_ytar);
 Double_t  e_xtar;
   tsimc->SetBranchAddress("H.gtr.x",&e_xtar);
 Double_t  gindex;
   tsimc->SetBranchAddress("H.gtr.index",&gindex);
 Double_t  e_reactz;
   tsimc->SetBranchAddress("H.react.z",&e_reactz);
 Double_t  e_delta;
   tsimc->SetBranchAddress("H.gtr.dp",&e_delta);
 Double_t  e_yptar;
   tsimc->SetBranchAddress("H.gtr.ph",&e_yptar);
 Double_t  e_xptar;
   tsimc->SetBranchAddress("H.gtr.th",&e_xptar);
 Double_t  e_yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&e_yfp);
 Double_t  e_ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&e_ypfp);
   Double_t  e_xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&e_xfp);
 Double_t  e_xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&e_xpfp);
  // Define histograms
    TH1F *hxtar = new TH1F("hxtar",Form("Run %d ; Xtar (cm);Counts",nrun), 300, -.3,.3);
    HList.Add(hxtar);
    TH1F *hztar = new TH1F("hztar",Form("Run %d ; Ztar (cm);Counts",nrun), 300, -10.,10.);
    HList.Add(hztar);
    TH2F *hxrast_yrast = new TH2F("hxrast_yrast",Form("Run %d ;  Xrast (cm) (+X beam right); Yrast (cm) (+Y up) ",nrun), 100, -.3,.3 ,100, -.3,.3);
    HList.Add(hxrast_yrast);
    TH2F *hxbeam_ybeam = new TH2F("hxbeam_ybeam",Form("Run %d ;  Xbeam (cm) (+X beam right); Ybeam (cm) (+Y up)",nrun), 100, -.3,.3 ,100, -.3,.3);
    HList.Add(hxbeam_ybeam);
    TH2F *hxrast_yrast_tr = new TH2F("hxrast_yrast_tr",Form("Run %d Track in spec;  Xrast (cm) (+X beam right) ; Yrast (cm) (+Y up) ",nrun), 100, -.3,.3 ,100, -.3,.3);
    HList.Add(hxrast_yrast_tr);
    TH2F *hztar_xtar = new TH2F("hztar_xtar",Form("Run %d ;  Ztar (cm) ;Xtar (cm) (+X down )",nrun), 300, -10.,10. ,100, -.3,.3);
    HList.Add(hztar_xtar);
    TH2F *hztar_yptar = new TH2F("hztar_yptar",Form("Run %d ;  Ztar (cm) ;Yptar (rad)",nrun), 300, -10.,10. ,100, -.04,.04);
    HList.Add(hztar_yptar);
    TH2F *hztar_xptar = new TH2F("hztar_xptar",Form("Run %d ;  Ztar (cm) ;Xptar (rad)",nrun), 300, -10.,10. ,100, -.1,.1);
    HList.Add(hztar_xptar);
    TH2F *hxptar_xtar = new TH2F("hxptar_xtar",Form("Run %d ;  Xptar (rad) ;Xtar (cm) (+X down )",nrun), 300, -.1,.1 ,100, -.3,.3);
    HList.Add(hxptar_xtar);
    TH2F *hztar_xrastbeam = new TH2F("hztar_xrastbeam",Form("Run %d ;  Ztar (cm) ;Xrast+Xbeam (cm) (+X beam right)",nrun), 300, -10.,10. ,100, -.3,.3);
    HList.Add(hztar_xrastbeam);
     TH1F *hxrast = new TH1F("hxrast",Form("Run %d ; Xrast (cm) (+X beam right);Counts",nrun), 300, -.3,.3);
    HList.Add(hxrast);
   TH1F *hyrast = new TH1F("hyrast",Form("Run %d ; Yrast (cm) (+Y up) ;Counts",nrun), 300, -.3,.3);
    HList.Add(hyrast);
    TH1F *hxbeam = new TH1F("hxbeam",Form("Run %d ; Xbeam (cm) (+X beam right);Counts",nrun), 300, -.3,.3);
    HList.Add(hxbeam);
    TH1F *hybeam = new TH1F("hybeam",Form("Run %d ; Ybeam (cm) (+Y up);Counts",nrun), 300, -.3,.3);
    HList.Add(hybeam);
    TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etrack norm;Counts",nrun), 120, 0.0,1.2);
    HList.Add(hetot);
  // loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
                hxrast_yrast->Fill(rastx,rasty);
                hxrast->Fill(rastx);
                hyrast->Fill(rasty);
                hxbeam->Fill(beamx);
                hxbeam_ybeam->Fill(beamx,beamy);
                hybeam->Fill(beamy);
		if (gindex>-1  ) { 
                hxrast_yrast_tr->Fill(rastx,rasty);
                hxtar->Fill(e_xtar);
                hztar->Fill(e_reactz);
                        hxptar_xtar->Fill(e_xptar,e_xtar);
                hztar_xtar->Fill(e_reactz,e_xtar);
                hztar_yptar->Fill(e_reactz,e_yptar);
                hztar_xptar->Fill(e_reactz,e_xptar);
                hztar_xrastbeam->Fill(e_reactz,rastx+beamx);
		//
		}
		//		
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
