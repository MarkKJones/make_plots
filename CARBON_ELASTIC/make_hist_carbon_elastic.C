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
#include <TCutG.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_carbon_elastic(Int_t nrun=2043,Int_t flag=1){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   TString basename;
   basename=Form("shms_replay_matrixopt_%d_orig",nrun);
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_carbon_elastic_hist.root";
 TObjArray HList(0);
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm);
 Double_t  ys;
   tsimc->SetBranchAddress("P.extcor.ysieve",&ys);
 Double_t  xs;
   tsimc->SetBranchAddress("P.extcor.xsieve",&xs);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
 Double_t  xtar;
   tsimc->SetBranchAddress("P.gtr.x",&xtar);
 Double_t  reactx;
   tsimc->SetBranchAddress("P.react.x",&reactx);
 Double_t  reacty;
   tsimc->SetBranchAddress("P.react.y",&reacty);
 Double_t  reactz;
   tsimc->SetBranchAddress("P.react.z",&reactz);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  eprime;
   tsimc->SetBranchAddress("P.gtr.p",&eprime);
 Double_t  eth;
   tsimc->SetBranchAddress("P.kin.scat_ang_rad",&eth);
 Double_t  yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
 Double_t  ntr;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntr);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
   //
    vector<double> save_calc_delta;
    vector<double> save_delta;
    vector<double> save_nx;
    vector<double> save_ny;
    vector<double> save_nd;
    vector<double> save_xsieve;
    vector<double> save_ysieve;
    vector<double> save_ex;
    vector<double> save_xfp;
    vector<double> save_yfp;
    vector<double> save_ypfp;
    vector<double> save_xpfp;
    vector<double> save_xptar_true;
    vector<double> save_ytar_true;
    vector<double> save_yptar_true;
    vector<double> save_xptar;
    vector<double> save_ytar;
    vector<double> save_yptar;
  //
    Int_t drun[19]={1657,1658,1659,1660,1661,1662,1663,1664,1665,1667,1668,1669,1670,1671,1672,1673,1674,1675,1676};
    Double_t nd_run;
    for (Int_t ndd=0;ndd<19;ndd++) {
      if (drun[ndd]==nrun) nd_run=ndd;
    }
    //
      TCutG* xpfp_xfp_cut[11];
      TCutG* ypfp_yfp_cut[11];
      Bool_t cuttest=0;
      Int_t ind_xpfp_xfp=-1;
      Int_t ind_ypfp_yfp=-1;
      Int_t ntot_xcut=0;
      Int_t ntot_ycut=0;
   if (flag==1) {
     TString outputcut=Form("../shms_carbon/hist/%d_cut.root",nrun);
    TFile fcut(outputcut);
    for (Int_t nc=0;nc<11;nc++) {
      xpfp_xfp_cut[nc] = (TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",nc));
      if (xpfp_xfp_cut[nc]) {
      Int_t npt = xpfp_xfp_cut[nc]->GetN();
      ntot_xcut++;
      cout << " xpfp v xfp cut = " << nc << " npts = " << npt << " ntot = " << ntot_xcut<< endl;
      }
    }
      cout << " get ypfp_v yfp cuts " << endl;
    for (Int_t nc=0;nc<11;nc++) {
      ypfp_yfp_cut[nc] = (TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",nc));
      if (ypfp_yfp_cut[nc]) {
      Int_t npt = ypfp_yfp_cut[nc]->GetN();
      ntot_ycut++;
      cout  << " ypfp v yfp cut = " << nc << " npts = " << npt << endl;
      }
    }
   }
   // Define histograms
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_sieve ; X_sieve",nrun), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys);
	  TH2F *hxpfp_ypfp = new TH2F("hxpfp_ypfp", Form("Run %d ; Yp_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.05,0.05);
          HList.Add(hxpfp_ypfp);
	  TH2F *hxfp_xpfp = new TH2F("hxfp_xpfp", Form("Run %d ; X_fp ; Xp_fp",nrun), 100,-.03,0.03,120,-.05,0.05);
          HList.Add(hxpfp_ypfp);
	  TH1F *hdelta = new TH1F("hdelta", Form("Run %d ; Delta ; Counts",nrun), 360,-15,30);
          HList.Add(hdelta);
	  TH1F *hthe = new TH1F("hthe", Form("Run %d ; Scat ang (deg) ; Counts",nrun), 360,5,10);
          HList.Add(hthe);
	  TH1F *hEx = new TH1F("hEx", Form("Run %d ; Ex ; Counts",nrun), 400,-50.,50);
          HList.Add(hEx);
	  TH2F *hxs_ys_bin[11][11];
	  TH1F *hxptar[11][11];
	  TH1F *hyptar[11][11];
	  TH1F *hxptardiff[11][11];
	  TH1F *hyptardiff[11][11];
	  TH1F *hytar[11][11];
	  for (Int_t ntx=0;ntx<11;ntx++) {
	  for (Int_t nty=0;nty<11;nty++) {
            hxs_ys_bin[ntx][nty]  = new TH2F(Form("hxs_ys_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d ; Y_sieve ; X_sieve",nrun,ntx,nty), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys_bin[ntx][nty]);
	  hxptar[ntx][nty]  = new TH1F(Form("hxptar_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d; Xp_tar ; Counts",nrun,ntx,nty), 100,-.06,.06);
          HList.Add(hxptar[ntx][nty]);
	  hyptar[ntx][nty]  = new TH1F(Form("hyptar_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d; Yp_tar ; Counts",nrun,ntx,nty), 100,-.04,.04);
          HList.Add(hyptar[ntx][nty]);
	  hxptardiff[ntx][nty]  = new TH1F(Form("hxptardiff_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d; Xp_tar ; Counts",nrun,ntx,nty), 100,-.06,.06);
          HList.Add(hxptardiff[ntx][nty]);
	  hyptardiff[ntx][nty]  = new TH1F(Form("hyptardiff_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d; Yp_tar ; Counts",nrun,ntx,nty), 100,-.04,.04);
          HList.Add(hyptardiff[ntx][nty]);
	  hytar[ntx][nty]  = new TH1F(Form("hytar_%d_%d",ntx,nty), Form("Run %d nx = %d ny = %d; Yp_tar ; Counts",nrun,ntx,nty), 100,-1.,1.);
          HList.Add(hytar[ntx][nty]);
	  }}
// loop over entries
	  Double_t Eb=2.221;
	  Double_t Ex=0;
          Double_t Mc=12*.9315;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (etotnorm>1. && ntr>0 ) {
	  cuttest = 0;
	if (flag==1) {
	for (Int_t nx=0;nx<11;nx++) {
	  if (xpfp_xfp_cut[nx] && xpfp_xfp_cut[nx]->IsInside(xfp,xpfp)) {
	   ind_xpfp_xfp = nx;
	     for (Int_t ny=0;ny<11;ny++) {
	      if (ypfp_yfp_cut[ny] && ypfp_yfp_cut[ny]->IsInside(yfp,ypfp)) {
                cuttest = 1;
	        ind_ypfp_yfp = ny;
              }
	     }
	  }
	}
	}
	Double_t epr=eprime;//*1.018;
		  Double_t A = 0.5;
		  Double_t B = epr-Mc-Eb;
		  Double_t C = epr*Eb*(TMath::Cos(eth)-1)+Mc*(Eb-epr);
                  Ex=999./1000.;
		  if ((B*B-4*A*C)>0) {
 		    Ex = (1./2./A)*(-B-TMath::Sqrt(B*B-4*A*C));
		  }
		  hEx->Fill(Ex*1000.);
                hxs_ys->Fill(ys,xs);
		hxpfp_ypfp->Fill(ypfp,xpfp);
                hdelta->Fill(delta);
                hthe->Fill(eth*57.3);
                if (cuttest && Ex*1000<20 &&  Ex*1000>-20) {
		 hxs_ys_bin[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(ys,xs);
		 hxptar[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(xptar);
		 hyptar[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(yptar);
		 hytar[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(ytar);
		 save_nx.push_back(float(ind_xpfp_xfp));
		 save_ny.push_back(float(ind_ypfp_yfp));
		 save_nd.push_back(nd_run);
	save_xfp.push_back(xfp);
	save_xpfp.push_back(xpfp);
	save_yfp.push_back(yfp);
	save_ypfp.push_back(ypfp);
	save_xptar.push_back(xptar);
	save_ytar.push_back(ytar);
	save_yptar.push_back(yptar);
	Double_t xptar_true = (5-ind_xpfp_xfp)*2.5/253.;
	Double_t yptar_true = ( (5-ind_ypfp_yfp)*1.64+(0.019+40.*.01*0.052)*delta-(0.00019+40*.01*.00052)*delta*delta)/253.;
		 hxptardiff[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(xptar-xptar_true);
		 hyptardiff[ind_xpfp_xfp][ind_ypfp_yfp]->Fill(yptar-yptar_true);
	save_xptar_true.push_back(xptar_true);
	save_ytar_true.push_back(0.);
	save_yptar_true.push_back(yptar_true);
	save_delta.push_back(delta);
	save_ex.push_back(Ex);		 
		}
		}
	}
	//
	if (flag==1) {
	    TString outputroot=Form("hist/"+basename+"_%d_tree.root",nrun);
   TFile hroot(outputroot,"recreate");
   TNtuple ntuple("fit","Select data","nx:ny:nd:xfp:xpfp:yfp:ypfp:xptar:yptar:ytar:xptarT:yptarT:ytarT:delta:ex");
   for (UInt_t n=0;n<save_delta.size();n++) {
     ntuple.Fill(save_nx[n],save_ny[n],save_nd[n],save_xfp[n],save_xpfp[n],save_yfp[n],save_ypfp[n],save_xptar[n],save_yptar[n],save_ytar[n],save_xptar_true[n],save_yptar_true[n],save_ytar_true[n],save_delta[n],save_ex[n]);
   }
   hroot.Write();
     }
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
