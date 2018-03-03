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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_tdcref(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000001);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="../../ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 const Int_t nfadcmax=10;
 Int_t nfadchits;
 Double_t fadc_ref[nfadcmax];
 tsimc->SetBranchAddress("Ndata.D.pfadc1",&nfadchits);
 tsimc->SetBranchAddress("D.pfadc1",&fadc_ref);
 //
//
 const Int_t nhfadcmax=10;
 Int_t nhfadchits;
 Double_t hfadc_ref[nhfadcmax];
 tsimc->SetBranchAddress("Ndata.D.hfadc1",&nhfadchits);
 tsimc->SetBranchAddress("D.hfadc1",&hfadc_ref);
 //
 const Int_t maxpl=21;
 Int_t nhodopad[4];
 Double_t hodo_neg_tdc[4][maxpl];
 Double_t hodo_pos_tdc[4][maxpl];
 tsimc->SetBranchAddress("Ndata.P.hod.1x.GoodNegTdcTimeUnCorr",&nhodopad[0]);
 tsimc->SetBranchAddress("P.hod.1x.GoodNegTdcTimeUnCorr",&hodo_neg_tdc[0]);
 tsimc->SetBranchAddress("Ndata.P.hod.1y.GoodNegTdcTimeUnCorr",&nhodopad[1]);
 tsimc->SetBranchAddress("P.hod.1y.GoodNegTdcTimeUnCorr",&hodo_neg_tdc[1]);
 tsimc->SetBranchAddress("Ndata.P.hod.2x.GoodNegTdcTimeUnCorr",&nhodopad[2]);
 tsimc->SetBranchAddress("P.hod.2x.GoodNegTdcTimeUnCorr",&hodo_neg_tdc[2]);
 tsimc->SetBranchAddress("Ndata.P.hod.2y.GoodNegTdcTimeUnCorr",&nhodopad[3]);
 tsimc->SetBranchAddress("P.hod.2y.GoodNegTdcTimeUnCorr",&hodo_neg_tdc[3]);
 tsimc->SetBranchAddress("P.hod.1x.GoodPosTdcTimeUnCorr",&hodo_pos_tdc[0]);
 tsimc->SetBranchAddress("P.hod.1y.GoodPosTdcTimeUnCorr",&hodo_pos_tdc[1]);
 tsimc->SetBranchAddress("P.hod.2x.GoodPosTdcTimeUnCorr",&hodo_pos_tdc[2]);
 tsimc->SetBranchAddress("P.hod.2y.GoodPosTdcTimeUnCorr",&hodo_pos_tdc[3]);

//
 const Int_t nreftdc=10;
 const Int_t nmultmax=6;
 Int_t ntdchits[nreftdc];
 Double_t tdc_ref[nreftdc][nmultmax];
 for (Int_t i=0;i<nreftdc;i++) {
   TString name=Form("Ndata.D.pdcref%d",i+1);
  tsimc->SetBranchAddress(name, &ntdchits[i]);
   name=Form("D.pdcref%d",i+1);
  tsimc->SetBranchAddress(name, &tdc_ref[i]);
 }
Double_t tdc_coin_ref[nreftdc];
 for (Int_t i=0;i<nreftdc;i++) {
  TString name=Form("T.coin.pDCREF%d_tdcTime",i+1);
  tsimc->SetBranchAddress(name, &tdc_coin_ref[i]);
 }
//
 const Int_t nhdcreftdc=4;
 Int_t nhdctdchits[nhdcreftdc];
 Double_t hdctdc_ref[nhdcreftdc][nmultmax];
 for (Int_t i=0;i<nhdcreftdc;i++) {
   TString name=Form("Ndata.D.hdcref%d",i+1);
  tsimc->SetBranchAddress(name, &nhdctdchits[i]);
   name=Form("D.hdcref%d",i+1);
  tsimc->SetBranchAddress(name, &hdctdc_ref[i]);
 }
Double_t hdctdc_coin_ref[nreftdc];
 for (Int_t i=0;i<nhdcreftdc;i++) {
  TString name=Form("T.coin.hDCREF%d_tdcTime",i+1);
  tsimc->SetBranchAddress(name, &hdctdc_coin_ref[i]);
 }
 //
 Double_t edtm_time;
 tsimc->SetBranchAddress("T.coin.pEDTM_tdcTime",&edtm_time);
 Double_t edtm_rawtime;
 tsimc->SetBranchAddress("T.coin.pEDTM_tdcTimeRaw",&edtm_rawtime);
 Double_t edtm_mult;
 tsimc->SetBranchAddress("T.coin.pEDTM_tdcMultiplicity",&edtm_mult);
 Double_t W;
 tsimc->SetBranchAddress("H.kin.primary.W",&W);
 Double_t hodo_starttime;
 tsimc->SetBranchAddress("P.hod.starttime",&hodo_starttime);
 //
 const UInt_t NPLANES=6;
   const TString plane_names[NPLANES] = {"x1", "x2", "v1", "v2", "u1", "u2"};
 TString pl_names[NPLANES] = {"x1", "x2", "v1", "v2", "u1", "u2"};
 Int_t ndchits[NPLANES];
 Double_t dc1_tdc[NPLANES][nmultmax];
 for (UInt_t i=0;i<NPLANES;i++) {
   TString name="Ndata.P.dc.1"+plane_names[i]+".time";
  tsimc->SetBranchAddress(name, &ndchits[i]);
   name="P.dc.1"+plane_names[i]+".time";
  tsimc->SetBranchAddress(name, &dc1_tdc[i]);
 }
 // SHMS
 const Int_t nrefhodotdc=2;
 const Int_t nmulthodomax=10;
 Int_t nhodotdchits[nrefhodotdc];
 Double_t hodotdc_ref[nrefhodotdc][nmulthodomax];
 for (Int_t i=0;i<nrefhodotdc;i++) {
   TString name=Form("Ndata.D.phodo%d",i+1);
  tsimc->SetBranchAddress(name, &nhodotdchits[i]);
   name=Form("D.phodo%d",i+1);
  tsimc->SetBranchAddress(name, &hodotdc_ref[i]);
 }
 //
 // HMS 
 const Int_t nrefhhodotdc=1;
 const Int_t nmulthhodomax=10;
 Int_t nhhodotdchits[nrefhhodotdc];
 Double_t hhodotdc_ref[nrefhhodotdc][nmulthhodomax];
 for (Int_t i=0;i<nrefhhodotdc;i++) {
   TString name=Form("Ndata.D.hhodo%d",i+1);
  tsimc->SetBranchAddress(name, &nhhodotdchits[i]);
   name=Form("D.hhodo%d",i+1);
  tsimc->SetBranchAddress(name, &hhodotdc_ref[i]);
 }
 //
 const Int_t nreftrigtdc=2;
 const Int_t nmulttrigmax=10;
 Int_t ntrigtdchits[nreftrigtdc];
 Double_t trigtdc_ref[nreftrigtdc][nmulttrigmax];
 for (Int_t i=0;i<nreftrigtdc;i++) {
   TString name=Form("Ndata.D.ptrig%d",i+1);
  tsimc->SetBranchAddress(name, &ntrigtdchits[i]);
   name=Form("D.ptrig%d",i+1);
  tsimc->SetBranchAddress(name, &trigtdc_ref[i]);
 }

 //
 // define histograms
 TH1F *hW = new TH1F("hW","; W (GeV)",175,.700,1.050);
 TH1F *hStarttime = new TH1F("hStarttime","; SHMS Hodo startime (ns)",300,-50,100);
 TH1F *hStarttime_mult1 = new TH1F("hStarttime_mult1",";  SHMS Hodo startime (ns) Ref mult=1 ",300,-50,100);
 TH1F *hStarttime_mult2 = new TH1F("hStarttime_mult2",";  SHMS Hodo startime (ns) Ref mult=2 ",300,-50,100);
 //
 TH1F *hhodonegtdc[4][maxpl];
 TH1F *hhodopostdc[4][maxpl];
 for (Int_t j=0;j<4;j++) {
 for (Int_t i=0;i<maxpl;i++) {
   hhodonegtdc[j][i] = new TH1F(Form("hhodonegtdc_pl%d_pad%d",j+1,i+1),Form(";SHMS Plane %d Pad %d Neg TDC uncorr; counts",j+1,i+1),300,-50,100);
   hhodopostdc[j][i] = new TH1F(Form("hhodopostdc_pl%d_pad%d",j+1,i+1),Form(";SHMS Plane %d Pad %d Pos TDC uncorr; counts",j+1,i+1),300,-50,100);
 } 
 }
 //
 TH1F *hmult[nreftdc];
 TH1F *hdccoin_ref_tdc[nreftdc];
 TH1F *hreftdc1[nreftdc][nmultmax];
 TH1F *hreftdc2[nreftdc][nmultmax];
 TH1F *hreftdc3[nreftdc][nmultmax];
 for (Int_t i=0;i<nreftdc;i++) {
   hmult[i] = new TH1F(Form("hmult_%d",i),Form(";SHMS Mult DCRef_%d; counts",i+1),5,0,5);
   hdccoin_ref_tdc[i] = new TH1F(Form("hdccoin_ref_tdc_%d",i),Form(";Coin SHMS DCRef_%d; counts",i+1),200,0,2000);
  for (Int_t j=0;j<nmultmax;j++) {
    hreftdc1[i][j] = new TH1F(Form("hreftdc_ref_%d_mult_%d_hit_1",i,j),Form(";SHMS DCRef_%d Time (ns) for Mult = %d, hit = 1; counts",i+1,j+1),200,0,2000);
    hreftdc2[i][j] = new TH1F(Form("hreftdc_ref_%d_mult_%d_hit_2",i,j),Form(";SHMS DCRef_%d Time (ns) for Mult = %d, hit = 2; counts",i+1,j+1),200,0,2000);
    hreftdc3[i][j] = new TH1F(Form("hreftdc_ref_%d_mult_%d_hit_3",i,j),Form(";SHMS DCRef_%d Time (ns) for Mult = %d, hit = 3; counts",i+1,j+1),200,0,2000);
  }
 }
 //
 //
 TH1F *hhdcmult[nhdcreftdc];
 TH1F *hhdccoin_ref_tdc[nhdcreftdc];
 TH1F *hrefhdctdc1[nhdcreftdc][nmultmax];
 TH1F *hrefhdctdc2[nhdcreftdc][nmultmax];
 TH1F *hrefhdctdc3[nhdcreftdc][nmultmax];
 for (Int_t i=0;i<nhdcreftdc;i++) {
   hhdcmult[i] = new TH1F(Form("hhdcmult_%d",i),Form(";HMS Mult DCRef_%d; counts",i+1),5,0,5);
   hhdccoin_ref_tdc[i] = new TH1F(Form("hddccoin_ref_tdc_%d",i),Form(";HMS Coin DCRef_%d; counts",i+1),200,0,2000);
  for (Int_t j=0;j<nmultmax;j++) {
    hrefhdctdc1[i][j] = new TH1F(Form("hrefhdctdc_ref_%d_mult_%d_hit_1",i,j),Form(";HMS DCRef_%d Time (ns) for Mult = %d, hit = 1; counts",i+1,j+1),200,0,2000);
    hrefhdctdc2[i][j] = new TH1F(Form("hrefhdctdc_ref_%d_mult_%d_hit_2",i,j),Form(";HMS DCRef_%d Time (ns) for Mult = %d, hit = 2; counts",i+1,j+1),200,0,2000);
    hrefhdctdc3[i][j] = new TH1F(Form("hrefhdctdc_ref_%d_mult_%d_hit_3",i,j),Form(";HMS DCRef_%d Time (ns) for Mult = %d, hit = 3; counts",i+1,j+1),200,0,2000);
  }
 }
 //
 TH1F *hdctdc[NPLANES][3];
 for (UInt_t i=0;i<NPLANES;i++) {
 for (Int_t j=0;j<3;j++) {
   TString name="hdctdc_1"+plane_names[i]+Form("_ mult_%d",j+1);
   TString n2=";DC "+plane_names[i]+Form(" Ref mult = %d",j+1);
   hdctdc[i][j] = new TH1F(name,n2,600,-100,500);
 }
 }
 //
 TH1F  *hfadcmult = new TH1F("hfadcmult",";SHMS Mult FADC reftime; counts",10,0,10);
 TH1F *hfadcreftdc[nfadcmax];
 for (Int_t j=0;j<nfadcmax;j++) {
     hfadcreftdc[j] = new TH1F(Form("hfadctdc_mult4_hit_%d",j+1),Form(";SHMS FADC Ref Time (ns)  mult = 4 hit = %d",j+1),3000,0,3000);
 }
 //
 TH1F  *hhfadcmult = new TH1F("hhfadcmult",";HMS Mult FADC reftime; counts",10,0,10);
 TH1F *hhfadcreftdc[nhfadcmax];
 for (Int_t j=0;j<nhfadcmax;j++) {
     hhfadcreftdc[j] = new TH1F(Form("hhfadctdc_mult4_hit_%d",j+1),Form(";HMS FADC Ref Time (ns)  mult = 4 hit = %d",j+1),3000,0,3000);
 }

 //
 TH1F *hmulthodo[nrefhodotdc];
 TH1F *hhodoreftdc1[nrefhodotdc][nmulthodomax];
 TH1F *hhodoreftdc2[nrefhodotdc][nmulthodomax];
 TH1F *hhodoreftdc3[nrefhodotdc][nmulthodomax];
 for (Int_t i=0;i<nrefhodotdc;i++) {
   hmulthodo[i] = new TH1F(Form("hmulthodo_%d",i),Form(";SHMS Mult Hodo Ref_%d; counts",i+1),5,0,5);
  for (Int_t j=0;j<nmulthodomax;j++) {
    hhodoreftdc1[i][j] = new TH1F(Form("hhodoreftdc_ref_%d_mult_%d_hit_1",i,j),Form(";SHMS Hodo Ref_%d Time (ns) for Mult = %d, hit = 1; counts",i+1,j+1),200,0,1000);
    hhodoreftdc2[i][j] = new TH1F(Form("hhodoreftdc_ref_%d_mult_%d_hit_2",i,j),Form(";SHMS Hodo Ref_%d Time (ns) for Mult = %d, hit = 2; counts",i+1,j+1),200,0,1000);
    hhodoreftdc3[i][j] = new TH1F(Form("hhodoreftdc_ref_%d_mult_%d_hit_3",i,j),Form(";SHMS Hodo Ref_%d Time (ns) for Mult = %d, hit = 3; counts",i+1,j+1),200,0,1000);
  }
 }
 //
 // HMS
 TH1F *hmulthhodo[nrefhhodotdc];
 TH1F *hhhodoreftdc1[nrefhhodotdc][nmulthhodomax];
 TH1F *hhhodoreftdc2[nrefhhodotdc][nmulthhodomax];
 TH1F *hhhodoreftdc3[nrefhhodotdc][nmulthhodomax];
 for (Int_t i=0;i<nrefhhodotdc;i++) {
   hmulthhodo[i] = new TH1F(Form("hmulthhodo_%d",i),Form(";HMS Mult hodo Ref_%d; counts",i+1),5,0,5);
  for (Int_t j=0;j<nmulthhodomax;j++) {
    hhhodoreftdc1[i][j] = new TH1F(Form("hhhodoreftdc_ref_%d_mult_%d_hit_1",i,j),Form(";HMS hodo Ref_%d Time (ns) for Mult = %d, hit = 1; counts",i+1,j+1),200,0,1000);
    hhhodoreftdc2[i][j] = new TH1F(Form("hhhodoreftdc_ref_%d_mult_%d_hit_2",i,j),Form(";HMS hodo Ref_%d Time (ns) for Mult = %d, hit = 2; counts",i+1,j+1),200,0,1000);
    hhhodoreftdc3[i][j] = new TH1F(Form("hhhodoreftdc_ref_%d_mult_%d_hit_3",i,j),Form(";HMS hodo Ref_%d Time (ns) for Mult = %d, hit = 3; counts",i+1,j+1),200,0,1000);
  }
 }
 //
 TH1F *hmulttrig[nreftrigtdc];
 TH1F *htrigreftdc1[nreftrigtdc][nmulttrigmax];
 TH1F *htrigreftdc2[nreftrigtdc][nmulttrigmax];
 TH1F *htrigreftdc3[nreftrigtdc][nmulttrigmax];
 for (Int_t i=0;i<nreftrigtdc;i++) {
   hmulttrig[i] = new TH1F(Form("hmulttrig_%d",i),Form(";Mult Trig Ref_%d; counts",i+1),10,0,10);
  for (Int_t j=0;j<nmulttrigmax;j++) {
    htrigreftdc1[i][j] = new TH1F(Form("htrigreftdc_ref_%d_mult_%d_hit_1",i,j),Form(";Trig Ref_%d Time (ns) for Mult = %d, hit = 1; counts",i+1,j+1),200,0,1000);
    htrigreftdc2[i][j] = new TH1F(Form("htrigreftdc_ref_%d_mult_%d_hit_2",i,j),Form(";Trig Ref_%d Time (ns) for Mult = %d, hit = 2; counts",i+1,j+1),200,0,1000);
    htrigreftdc3[i][j] = new TH1F(Form("htrigreftdc_ref_%d_mult_%d_hit_3",i,j),Form(";Trig Ref_%d Time (ns) for Mult = %d, hit = 3; counts",i+1,j+1),200,0,1000);
  }
 }
 // loop over data
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (edtm_mult==0) {
		  // SHMS DC ref time
                 for (Int_t i=0;i<nreftdc;i++) {
		   hdccoin_ref_tdc[i]->Fill(tdc_coin_ref[i]);
		   hmult[i]->Fill(ntdchits[i]);
                    for (Int_t j=0;j<ntdchits[i];j++) {
		    if (j==0) hreftdc1[i][ntdchits[i]-1]->Fill(tdc_ref[i][j]*.1);
		    if (j==1) hreftdc2[i][ntdchits[i]-1]->Fill(tdc_ref[i][j]*.1);
		    if (j==2) hreftdc3[i][ntdchits[i]-1]->Fill(tdc_ref[i][j]*.1);
                    }
                  }
		  // HMS DC ref time
                 for (Int_t i=0;i<nhdcreftdc;i++) {
		   hhdccoin_ref_tdc[i]->Fill(hdctdc_coin_ref[i]);
		   hhdcmult[i]->Fill(nhdctdchits[i]);
                    for (Int_t j=0;j<nhdctdchits[i];j++) {
		    if (j==0) hrefhdctdc1[i][nhdctdchits[i]-1]->Fill(hdctdc_ref[i][j]*.1);
		    if (j==1) hrefhdctdc2[i][nhdctdchits[i]-1]->Fill(hdctdc_ref[i][j]*.1);
		    if (j==2) hrefhdctdc3[i][nhdctdchits[i]-1]->Fill(hdctdc_ref[i][j]*.1);
                    }
                  }
		  // SHMS Hodo ref time
                 for (Int_t i=0;i<nrefhodotdc;i++) {
		   hmulthodo[i]->Fill(nhodotdchits[i]);
                    for (Int_t j=0;j<nhodotdchits[i];j++) {
		    if (j==0) hhodoreftdc1[i][nhodotdchits[i]-1]->Fill(hodotdc_ref[i][j]*.1);
		    if (j==1) hhodoreftdc2[i][nhodotdchits[i]-1]->Fill(hodotdc_ref[i][j]*.1);
		    if (j==2) hhodoreftdc3[i][nhodotdchits[i]-1]->Fill(hodotdc_ref[i][j]*.1);
                    }
                  }
		 //
		  // HMS Hodo ref time
                 for (Int_t i=0;i<nrefhhodotdc;i++) {
		   hmulthhodo[i]->Fill(nhhodotdchits[i]);
                    for (Int_t j=0;j<nhhodotdchits[i];j++) {
		    if (j==0) hhhodoreftdc1[i][nhhodotdchits[i]-1]->Fill(hhodotdc_ref[i][j]*.1);
		    if (j==1) hhhodoreftdc2[i][nhhodotdchits[i]-1]->Fill(hhodotdc_ref[i][j]*.1);
		    if (j==2) hhhodoreftdc3[i][nhhodotdchits[i]-1]->Fill(hhodotdc_ref[i][j]*.1);
                    }
                  }
		 //
		  for (Int_t i=0;i<nreftrigtdc;i++) {
		   hmulttrig[i]->Fill(ntrigtdchits[i]);
		    for (Int_t j=0;j<ntrigtdchits[i];j++) {
		    if (j==0) htrigreftdc1[i][ntrigtdchits[i]-1]->Fill(trigtdc_ref[i][j]*.1);
		    if (j==1) htrigreftdc2[i][ntrigtdchits[i]-1]->Fill(trigtdc_ref[i][j]*.1);
		    if (j==2) htrigreftdc3[i][ntrigtdchits[i]-1]->Fill(trigtdc_ref[i][j]*.1);
		    }
		  }
		 // SHMS fadc
		 hfadcmult->Fill(nfadchits);
		   for (Int_t j=0;j<nfadchits;j++) {
		      if (nfadchits==4) hfadcreftdc[j]->Fill(fadc_ref[j]*.0625);
		   }
		 // HMS fadc
		 hhfadcmult->Fill(nhfadchits);
		   for (Int_t j=0;j<nhfadchits;j++) {
		      if (nhfadchits==4) hhfadcreftdc[j]->Fill(hfadc_ref[j]*.0625);
		   }
		 // dc time
		 for (UInt_t i=0;i<NPLANES;i++) {
		   for (Int_t j=0;j<ndchits[i];j++) {
		     if( ntdchits[0] == 1) hdctdc[i][0]->Fill(dc1_tdc[i][j]);
		     if( ntdchits[0] == 2) hdctdc[i][1]->Fill(dc1_tdc[i][j]);
		     if( ntdchits[0] == 3) hdctdc[i][2]->Fill(dc1_tdc[i][j]);
		 }
		 }
		 //
                 for (Int_t j=0;j<4;j++) {
                 for (Int_t i=0;i<nhodopad[j];i++) {
		   hhodonegtdc[j][i]->Fill(hodo_neg_tdc[j][i]);
		   hhodopostdc[j][i]->Fill(hodo_pos_tdc[j][i]);
                  } 
                   }
		 hW->Fill(W);
		 hStarttime->Fill(hodo_starttime);
		 if (nhodotdchits[0]==1) hStarttime_mult1->Fill(hodo_starttime);
		 if (nhodotdchits[0]==2) hStarttime_mult2->Fill(hodo_starttime);
		} // if edtmmult
		//
	}
	// SHMS DC
 TCanvas *c1[nreftdc];
for (Int_t nref=0;nref<nreftdc;nref++) {
  c1[nref] = new TCanvas(Form("c1_%d",nref),Form("SHMS DC ref_%d",nref+1), 900,700);
c1[nref]->Divide(2,4);
c1[nref]->cd(1);
 gPad->SetLogy();
 hmult[nref]->Draw();
c1[nref]->cd(2);
 gPad->SetLogy();
 hreftdc1[nref][0]->Draw();
c1[nref]->cd(3);
 gPad->SetLogy();
 hreftdc1[nref][1]->Draw();
c1[nref]->cd(4);
 gPad->SetLogy();
 hreftdc1[nref][2]->Draw();
c1[nref]->cd(5);
 gPad->SetLogy();
 hreftdc2[nref][1]->Draw();
c1[nref]->cd(6);
 gPad->SetLogy();
 hreftdc2[nref][2]->Draw();
c1[nref]->cd(7);
 gPad->SetLogy();
 hreftdc3[nref][2]->Draw();
c1[nref]->cd(8);
 gPad->SetLogy();
 hdccoin_ref_tdc[nref]->Draw();
 outputpdf= "plots/"+basename+Form("_SHMS_ref_%d",nref+1)+".pdf";
 c1[nref]->Print(outputpdf);
 c1[nref]->Close();
 }
	// HMS DC
 TCanvas *chdc[nhdcreftdc];
for (Int_t nref=0;nref<nhdcreftdc;nref++) {
  chdc[nref] = new TCanvas(Form("chdc_%d",nref),Form("HMS DC ref_%d",nref+1), 900,700);
chdc[nref]->Divide(2,4);
chdc[nref]->cd(1);
 gPad->SetLogy();
 hhdcmult[nref]->Draw();
chdc[nref]->cd(2);
 gPad->SetLogy();
 hrefhdctdc1[nref][0]->Draw();
chdc[nref]->cd(3);
 gPad->SetLogy();
 hrefhdctdc1[nref][1]->Draw();
chdc[nref]->cd(4);
 gPad->SetLogy();
 hrefhdctdc1[nref][2]->Draw();
chdc[nref]->cd(5);
 gPad->SetLogy();
 hrefhdctdc2[nref][1]->Draw();
chdc[nref]->cd(6);
 gPad->SetLogy();
 hrefhdctdc2[nref][2]->Draw();
chdc[nref]->cd(7);
 gPad->SetLogy();
 hrefhdctdc3[nref][2]->Draw();
chdc[nref]->cd(8);
 gPad->SetLogy();
 hhdccoin_ref_tdc[nref]->Draw();
 outputpdf= "plots/"+basename+Form("_HMS_dc_ref_%d",nref+1)+".pdf";
 chdc[nref]->Print(outputpdf);
 //chdc[nref]->Close();
 }
	// SHMS
 TCanvas *chodo[nrefhodotdc];
for (Int_t nref=0;nref<nrefhodotdc;nref++) {
  chodo[nref] = new TCanvas(Form("chodo_%d",nref),Form("SHMS Hodo ref_%d",nref+1), 900,700);
chodo[nref]->Divide(2,4);
chodo[nref]->cd(1);
 gPad->SetLogy();
 hmulthodo[nref]->Draw();
chodo[nref]->cd(2);
 gPad->SetLogy();
 hhodoreftdc1[nref][0]->Draw();
chodo[nref]->cd(3);
 gPad->SetLogy();
 hhodoreftdc1[nref][1]->Draw();
chodo[nref]->cd(4);
 gPad->SetLogy();
 hhodoreftdc1[nref][2]->Draw();
chodo[nref]->cd(5);
 gPad->SetLogy();
 hhodoreftdc2[nref][1]->Draw();
chodo[nref]->cd(6);
 gPad->SetLogy();
 hhodoreftdc2[nref][2]->Draw();
chodo[nref]->cd(7);
 gPad->SetLogy();
 hhodoreftdc3[nref][2]->Draw();
 outputpdf= "plots/"+basename+Form("_SHMS_hodoref_%d",nref+1)+".pdf";
 chodo[nref]->Print(outputpdf);
 }
	// HMS
 TCanvas *chhodo[nrefhhodotdc];
for (Int_t nref=0;nref<nrefhhodotdc;nref++) {
  chhodo[nref] = new TCanvas(Form("chhodo_%d",nref),Form("HMS hodo ref_%d",nref+1), 900,700);
chhodo[nref]->Divide(2,4);
chhodo[nref]->cd(1);
 gPad->SetLogy();
 hmulthhodo[nref]->Draw();
chhodo[nref]->cd(2);
 gPad->SetLogy();
 hhhodoreftdc1[nref][0]->Draw();
chhodo[nref]->cd(3);
 gPad->SetLogy();
 hhhodoreftdc1[nref][1]->Draw();
chhodo[nref]->cd(4);
 gPad->SetLogy();
 hhhodoreftdc1[nref][2]->Draw();
chhodo[nref]->cd(5);
 gPad->SetLogy();
 hhhodoreftdc2[nref][1]->Draw();
chhodo[nref]->cd(6);
 gPad->SetLogy();
 hhhodoreftdc2[nref][2]->Draw();
chhodo[nref]->cd(7);
 gPad->SetLogy();
 hhhodoreftdc3[nref][2]->Draw();
 outputpdf= "plots/"+basename+Form("_HMS_hodoref_%d",nref+1)+".pdf";
 chhodo[nref]->Print(outputpdf);
 }
	//
 TCanvas *ctrig[nreftrigtdc];
for (Int_t nref=0;nref<nreftrigtdc;nref++) {
  ctrig[nref] = new TCanvas(Form("ctrig_%d",nref),Form("Trig ref_%d",nref+1), 900,700);
ctrig[nref]->Divide(2,4);
ctrig[nref]->cd(1);
 gPad->SetLogy();
 hmulttrig[nref]->Draw();
ctrig[nref]->cd(2);
 gPad->SetLogy();
 htrigreftdc1[nref][0]->Draw();
ctrig[nref]->cd(3);
 gPad->SetLogy();
 htrigreftdc1[nref][1]->Draw();
ctrig[nref]->cd(4);
 gPad->SetLogy();
 htrigreftdc1[nref][2]->Draw();
ctrig[nref]->cd(5);
 gPad->SetLogy();
 htrigreftdc2[nref][1]->Draw();
ctrig[nref]->cd(6);
 gPad->SetLogy();
 htrigreftdc2[nref][2]->Draw();
ctrig[nref]->cd(7);
 gPad->SetLogy();
 htrigreftdc3[nref][2]->Draw();
 outputpdf= "plots/"+basename+Form("_SHMS_trigref_%d",nref+1)+".pdf";
 ctrig[nref]->Print(outputpdf);
 }
	//
 TCanvas *cdcplane[3];
for (Int_t nref=0;nref<3;nref++) {
  cdcplane[nref] = new TCanvas(Form("cdcplane_%d",nref),Form("DC DC1 Ref mult = %d",nref+1), 900,700);
cdcplane[nref]->Divide(2,3);
 for (UInt_t npl=0;npl<NPLANES;npl++) {
 cdcplane[nref]->cd(npl+1);
 gPad->SetLogy();
 hdctdc[npl][nref]->Draw();
 }
 outputpdf= "plots/"+basename+Form("_SHMS_dc1_tdc_mult_%d",nref+1)+".pdf";
 cdcplane[nref]->Print(outputpdf);
 }
	//
 TCanvas *cfadcref;
 cfadcref = new TCanvas("cfadcref","SHMS FADC Ref mult = 4", 900,700);
cfadcref->Divide(2,3);
 cfadcref->cd(1);
 hfadcmult->Draw();
 for (UInt_t npl=1;npl<5;npl++) {
 cfadcref->cd(npl+1);
 gPad->SetLogy();
 hfadcreftdc[npl-1]->Draw();
 }
 outputpdf= "plots/"+basename+Form("_SHMS_fadc_mult_%d",4)+".pdf";
 cfadcref->Print(outputpdf);
 //
	// HMS
 TCanvas *chfadcref;
 chfadcref = new TCanvas("chfadcref","HMS FADC Ref mult = 4", 900,700);
chfadcref->Divide(2,3);
 chfadcref->cd(1);
 hhfadcmult->Draw();
 for (UInt_t npl=1;npl<5;npl++) {
 chfadcref->cd(npl+1);
 gPad->SetLogy();
 hhfadcreftdc[npl-1]->Draw();
 }
 outputpdf= "plots/"+basename+Form("_HMS_fadc_mult_%d",4)+".pdf";
 chfadcref->Print(outputpdf);
 //
 TCanvas *cstime;
 cstime = new TCanvas("cstime","Start time", 900,700);
cstime->Divide(1,3);
cstime->cd(1);
 hStarttime->Draw();
cstime->cd(2);
 hStarttime_mult1->Draw();
cstime->cd(3);
 hStarttime_mult2->Draw();
 outputpdf= "plots/"+basename+"_SHMS_starttime"+".pdf";
 cstime->Print(outputpdf);
	//
	//
// end of file
}
