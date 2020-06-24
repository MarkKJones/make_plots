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
#include <TCutG.h>
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
#include <iostream>
#include <fstream>
using namespace std;

void make_hist_hms_optics(Int_t nrun=1813,Bool_t CutYtarFlag=kFALSE,Bool_t cflag=kFALSE,Bool_t newfit=kFALSE){
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
 //  Get info for that optics run
 TString OpticsFile = "HMS_OPTICS_A1N/list_of_optics_run.dat";
   ifstream file_optics(OpticsFile.Data());
 TString opticsline;
  TString OpticsID="";
  Int_t RunNum=0.;
  Double_t CentAngle=0.;
  Int_t SieveFlag=1;
  Int_t NumFoil=0;
  TString temp;
  if (file_optics.is_open()) {
    //
    cout << " Open file = " << OpticsFile << endl;
    while (RunNum!=nrun  ) {
      temp.ReadToDelim(file_optics,',');
      cout << temp << endl;
      if (temp.Atoi() == nrun) {
	RunNum = temp.Atoi();
      } else {
	temp.ReadLine(file_optics);
      }
    }
    if (RunNum==nrun) {
      temp.ReadToDelim(file_optics,',');
      OpticsID = temp;
      temp.ReadToDelim(file_optics,',');
      CentAngle = temp.Atof();
      temp.ReadToDelim(file_optics,',');
      NumFoil = temp.Atoi();
      temp.ReadToDelim(file_optics);
      SieveFlag = temp.Atoi();
    }
  } else {
    cout << " No file = " << OpticsFile << endl;    
  }
  cout << RunNum << " " << OpticsID << " " << CentAngle << " " << NumFoil << " " << SieveFlag << endl;
  if (NumFoil==0) return;
 //
   TString inputroot;
   TString outputhist;
   if ( newfit) {
   inputroot="hms_optics_newfit/hms_replay_matrixopt_"+OpticsID+"_-5.root";
   outputhist= "hist/Optics_"+OpticsID+"_hms_newfit5_hist.root";
   } else {
   inputroot="hms_optics_rootfiles/hms_replay_matrixopt_"+OpticsID+"_-1.root";
   outputhist= "hist/Optics_"+OpticsID+"_hms_hist.root";
   }
   cout << " input root = " << inputroot << endl;
 TObjArray HList(0);
 //
	static const Int_t ndelcut=5;
	Double_t delcut[ndelcut+1]={-10.,-7.5,-5,0,5,10};
 //
  TString YtarDeltaCutFile;
  TFile *fYtarDeltaCut;
  vector <TCutG*> ytar_delta_cut;
  if (CutYtarFlag) {
    YtarDeltaCutFile="cuts/ytar_delta_"+OpticsID+"_cut.root";
    fYtarDeltaCut = new TFile(YtarDeltaCutFile);
    cout << " Cut file = " << YtarDeltaCutFile << endl;
   for (Int_t nc=0;nc<NumFoil;nc++) {
    fYtarDeltaCut->cd();
      TCutG* tempcut = (TCutG*)gROOT->FindObject(Form("delta_vs_ytar_cut_foil%d",nc));
      if (tempcut) {
      Int_t npt = tempcut->GetN();
      cout << "hYtarDelta_cut = " << nc << " npts = " << npt << endl;
      ytar_delta_cut.push_back(tempcut);
      } else {
      cout << " No hYtarDelta_cut = " << nc << endl;
      }
   }
  }
 //
  TString outCutFile;
  TFile *fcut;
  vector<vector<vector<TCutG*> > > ypfp_ypfp_cut;
  vector<vector<vector<Int_t> > > ypfp_ypfp_cut_flag;
  ypfp_ypfp_cut.resize(NumFoil);
  ypfp_ypfp_cut_flag.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
          ypfp_ypfp_cut[nf].resize(ndelcut);
          ypfp_ypfp_cut_flag[nf].resize(ndelcut);
   }
  if (cflag) {
  outCutFile="cuts/YpFpYFp_"+OpticsID+"_cut.root";
    fcut = new TFile(outCutFile);
    cout << " Cut file = " << outCutFile << endl;
    fcut->cd();
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	for  (Int_t nd=0;nd<ndelcut;nd++) {
        for (Int_t nc=0;nc<11;nc++) {
	  TCutG* tempg  = (TCutG*)gROOT->FindObject(Form("hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d",nc,nf,nd));
	  if (tempg)  {
	    //cout << "hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
	  ypfp_ypfp_cut[nf][nd].push_back(tempg);
      } else {
	    //cout << " No hYpFpYFp_cut = " << nc << " " << nf << " " << nd << endl;
     	  ypfp_ypfp_cut[nf][nd].push_back(tempg);
      }
	}}}
  }
//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("H.cer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("H.cal.etottracknorm",&etracknorm);
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
 Double_t  ysieve;
   tsimc->SetBranchAddress("H.extcor.ysieve",&ysieve);
 Double_t  xsieve;
   tsimc->SetBranchAddress("H.extcor.xsieve",&xsieve);

   // Define histograms
	TH1F *hetot = new TH1F("hetot",Form("Run %d ; Etotnorm ; Counts",nrun),100,0.,2.);
        HList.Add(hetot);
	TH1F *hngsum = new TH1F("hngsum",Form("Run %d ; NG Npe SUM ; Counts",nrun),100,0.,40.);
	HList.Add(hngsum);
	TH1F *hytar = new TH1F("hytar",Form("Run %d ; Ytar; Counts",nrun),500,-10.,10.);
	HList.Add(hytar);
	TH2F *hXptarDelta = new TH2F("hXptarDelta",Form("Run %d ; Xptar ; Delta",nrun),120,-.06,.06,100,-15.,25.);
	HList.Add(hXptarDelta);
	TH2F *hYptarDelta = new TH2F("hYptarDelta",Form("Run %d ; Yptar ; Delta",nrun),120,-.04,.04,100,-15.,25.);
	HList.Add(hYptarDelta);
	TH2F *hYtarDelta = new TH2F("hYtarDelta",Form("Run %d ; Ytar ; Delta",nrun),100,-10.,10.,100,-15.,25.);
	HList.Add(hYtarDelta);
	//
	TH2F *hYpFpXFp = new TH2F("hYpFpXFp",Form("Run %d ; Ypfp ; Xfp",nrun),100,-.04,.04,100,-40.,40.);
	HList.Add(hYpFpXFp);
	TH2F *hYFpXFp = new TH2F("hYFpXFp",Form("Run %d ; Yfp ; Xfp",nrun),100,-40,40,100,-40.,40.);
	HList.Add(hYFpXFp);
	TH2F *hXpFpXFp = new TH2F("hXpFpXFp",Form("Run %d ; Xpfp ; Xfp",nrun),100,-.1,.1,100,-50.,50.);
	HList.Add(hXpFpXFp);
	TH2F *hYtarYptar = new TH2F("hYtarYptar",Form("Run %d ; Yptar ; Ytar",nrun),100,-.05,.05,100,-10.,10.);
	HList.Add(hYtarYptar);
	TH2F *hZtarDelta = new TH2F("hZtarDelta",Form("Run %d ; Ztar ; Delta",nrun),140,-35.,25.,100,-15.,25.);
	HList.Add(hZtarDelta);
	//
	vector <TH2F*> hYsDelta;
	hYsDelta.resize(NumFoil);
	vector <TH2F*> hXsDelta;
	hXsDelta.resize(NumFoil);
	vector <TH2F*> hYpFpYFp;
	hYpFpYFp.resize(NumFoil);
	vector<vector<vector<TH2F*> > > hYsXs_DelCut_YpYfpCut;
	vector<vector<vector<TH1F*> > > hXs_DelCut_YpYfpCut;
	vector<vector<TH2F*> > hYsXs_DelCut;
	vector<vector<TH2F*> > hYpFpYFp_DelCut;
	vector<vector<TH2F*> > hXpFpXFp_DelCut;
	cout << " setup DelCut 2d" << endl;
	hYsXs_DelCut.resize(NumFoil);
	hYsXs_DelCut_YpYfpCut.resize(NumFoil);
	hXs_DelCut_YpYfpCut.resize(NumFoil);
	hYpFpYFp_DelCut.resize(NumFoil);
	hXpFpXFp_DelCut.resize(NumFoil);
	for  (Int_t nf=0;nf<NumFoil;nf++) {
	  hYsXs_DelCut[nf].resize(ndelcut);
	  hYsXs_DelCut_YpYfpCut[nf].resize(ndelcut);
	  hXs_DelCut_YpYfpCut[nf].resize(ndelcut);
	  hYpFpYFp_DelCut[nf].resize(ndelcut);
	  hXpFpXFp_DelCut[nf].resize(ndelcut);
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	  hYsXs_DelCut_YpYfpCut[nf][nd].resize(11);
	  hXs_DelCut_YpYfpCut[nf][nd].resize(11);
	}
	}
	cout << " finish setup DelCut 2d" << endl;
	for  (Int_t nc=0;nc<NumFoil;nc++) {
	 hYsDelta[nc] = new TH2F(Form("hYsDelta_Foil_%d",nc),Form("Run %d Foil %d; Ys ; Delta",nc,nrun),100,-12,12,50,-15.,25.);
	HList.Add(hYsDelta[nc]);
	 hXsDelta[nc] = new TH2F(Form("hXsDelta_Foil_%d",nc),Form("Run %d Foil %d; Xs ; Delta",nc,nrun),100,-15,15,50,-15.,25.);
	HList.Add(hXsDelta[nc]);
	  hYpFpYFp[nc] = new TH2F(Form("hYpFpYFp_%d",nc),Form("Run %d Foil %d; Ypfp ; Yfp",nrun,nc),100,-.045,.045,100,-30.,30.);
	HList.Add(hYpFpYFp[nc]);
	for  (Int_t nd=0;nd<ndelcut;nd++) {
	 hYsXs_DelCut[nc][nd]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Ys ; Xs",nrun,nc,(delcut[nd+1]+delcut[nd])/2),50,-12,12,100,-15.,15.);
	HList.Add(hYsXs_DelCut[nc][nd]);
	for  (Int_t ny=0;ny<11;ny++) {
	  hYsXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH2F(Form("hYsXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d DelCut %3.1f Ys=%d; Ys ; Xs",nrun,nc,(delcut[nd+1]+delcut[nd])/2,ny),100,-12,12,100,-15.,15.);
	HList.Add(hYsXs_DelCut_YpYfpCut[nc][nd][ny]);
	  hXs_DelCut_YpYfpCut[nc][nd][ny]  = new TH1F(Form("hXs_Foil_%d_DelCut_%d_FpCut_%d",nc,nd,ny),Form("Run %d Foil %d DelCut %3.1f Ys=%d; Xs",nrun,nc,(delcut[nd+1]+delcut[nd])/2,ny),100,-15.,15.);
	HList.Add(hXs_DelCut_YpYfpCut[nc][nd][ny]);
	}
	 hYpFpYFp_DelCut[nc][nd]  = new TH2F(Form("hYpFpYFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Ypfp ; Yfp",nrun,nc,(delcut[nd+1]+delcut[nd])/2),75,-.045,.045,150,-35.,35.);
	HList.Add(hYpFpYFp_DelCut[nc][nd]);
	hXpFpXFp_DelCut[nc][nd]= new TH2F(Form("hXpFpXFp_%d_DelCut_%d",nc,nd),Form("Run %d Foil %d DelCut %3.1f; Xpfp ; Xfp",nrun,nc,(delcut[nd+1]+delcut[nd])/2),100,-.1,.1,100,-50.,50.);
	HList.Add(hXpFpXFp_DelCut[nc][nd]);
	}
        }	  
	//
// loop over entries
Long64_t nentries = tsimc->GetEntries();
 cout << " start loop " << nentries << endl;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%50000==0) cout << " Entry = " << i << endl;
		if (sumnpe > 6.) hetot->Fill(etracknorm);
		if (etracknorm>.8) hngsum->Fill(sumnpe);
		if (sumnpe > 6. && delta>-15 && delta<24) {
		  hytar->Fill(ytar);
		  hXptarDelta->Fill(xptar,delta);
		  hYptarDelta->Fill(yptar,delta);
		  hYtarDelta->Fill(ytar,delta);
		  hYtarYptar->Fill(yptar,ytar);
		  hXpFpXFp->Fill(xpfp,xfp);
		  hYFpXFp->Fill(yfp,xfp); 
		  hYtarYptar->Fill(yptar,ytar);
		  hZtarDelta->Fill(reactz,delta);
	          for  (UInt_t nc=0;nc<ytar_delta_cut.size();nc++) {
		       if (ytar_delta_cut[nc]->IsInside(ytar,delta))	{ 
		       hYsDelta[nc]->Fill(ysieve,delta);
		       hXsDelta[nc]->Fill(xsieve,delta);
		       hYpFpYFp[nc]->Fill(ypfp,yfp);
                            for  (UInt_t nd=0;nd<ndelcut;nd++) {
		             if ( delta >=delcut[nd] && delta <delcut[nd+1]) {
                               hYsXs_DelCut[nc][nd]->Fill(ysieve,xsieve); 
		               hYpFpYFp_DelCut[nc][nd]->Fill(ypfp,yfp);
		               hXpFpXFp_DelCut[nc][nd]->Fill(xpfp,xfp);
                               for  (UInt_t ny=0;ny<11;ny++) {
			        if (cflag && ypfp_ypfp_cut[nc][nd][ny] && ypfp_ypfp_cut[nc][nd][ny]->IsInside(ypfp,yfp)) {
				hYsXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(ysieve,xsieve);
				hXs_DelCut_YpYfpCut[nc][nd][ny]->Fill(xsieve);
			        }
			       }
			     }			     
		            }
			 }
		  }
		}
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
