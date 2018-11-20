#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TCutG.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_hist_shms_ztar_sieve(TString basename,TString label) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
 //
 outputpdf="plots/"+basename+"_shms_ztar_sieve_hist.pdf";
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_shms_ztar_sieve_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   TH2F *fhist[2];
       fhist[0] = (TH2F*)fhistroot->Get("hytar_yptar");
       fhist[1] = (TH2F*)fhistroot->Get("hztar_yptar");
       TH1F *zhist[9];
       TH1F *yphist[9];
       for (Int_t n=0;n<9;n++) {
	 zhist[n] = (TH1F*)fhistroot->Get(Form("hztar_%d",n));
	 yphist[n] = (TH1F*)fhistroot->Get(Form("hyptar_%d",n));
       }
       TCanvas *ctime2;
       Double_t ztar_cent[9];
       Double_t ztar_error[9];
      ctime2 = new TCanvas("ctime2","ztar", 700,700);
     ctime2->Divide(1,1);
     Double_t xlo;
     Double_t xhi;
 	for (Int_t iz = 0; iz < 9; iz++) {
	  ctime2->cd(1);
	    cout << iz << endl;
	  zhist[iz]->Draw();
	  TCutG *cutg;     
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
      gPad->Update();
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg=(TCutG*)(tempg->Clone());
      		  cutg->SetName("cuttemp");
     		  //cutg[nf][nh]->Print();
		  Int_t npts=cutg->GetN();
		  cutg->RemovePoint(npts-1);
		  Double_t ydum;
		  cutg->GetPoint(0,xlo,ydum);
		  cutg->GetPoint(1,xhi,ydum);
		  cout << " xlo " << xlo << " xhi " << xhi << endl;
     }
	  zhist[iz]->Fit("gaus","","",xlo,xhi);
	  TF1 *fit = zhist[iz]->GetFunction("gaus");

	  ztar_cent[iz] = fit->GetParameter(1);
	  ztar_error[iz] = fit->GetParError(1);
	}
       TCanvas *ctime3;
       Double_t yptar_cent[9];
       Double_t yptar_error[9];
      ctime3 = new TCanvas("ctime3","yptar", 700,700);
     ctime3->Divide(1,1);
 	for (Int_t iz = 0; iz < 9; iz++) {
	  ctime3->cd(1);
	    cout << iz << endl;
	  yphist[iz]->Draw();
	  TCutG *cutg;     
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
      gPad->Update();
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg=(TCutG*)(tempg->Clone());
      		  cutg->SetName("cuttemp");
     		  //cutg[nf][nh]->Print();
		  Int_t npts=cutg->GetN();
		  cutg->RemovePoint(npts-1);
		  Double_t ydum;
		  cutg->GetPoint(0,xlo,ydum);
		  cutg->GetPoint(1,xhi,ydum);
		  cout << " xlo " << xlo << " xhi " << xhi << endl;
     }
	  yphist[iz]->Fit("gaus","","",xlo,xhi);
	  TF1 *fit = yphist[iz]->GetFunction("gaus");
	  yptar_cent[iz] = fit->GetParameter(1);
	  yptar_error[iz] = fit->GetParError(1);
	}
        TCanvas *ctime4;
      ctime4 = new TCanvas("ctime4","yptar", 700,700);
     ctime4->Divide(1,1);
     ctime4->cd(1);
     TGraphErrors *gr;
     gr = new TGraphErrors(9,yptar_cent,ztar_cent,yptar_error,ztar_error);
gr->Draw("ALP");
 gr->SetTitle(label);
 gr->GetYaxis()->SetTitle( "Z tar (cm) ");
 gr->GetXaxis()->SetTitle( "Yp tar (rad) ");
 TF1 *linfit;
 linfit = new TF1("linfit","pol1",-1.,1.);
 gr->Fit("linfit");
     TCanvas *ctime;
     ctime = new TCanvas("ctime","2d", 700,700);
     ctime->Divide(1,1);
     fhist[1]->Draw("colz");
        linfit->Draw("same");
}
