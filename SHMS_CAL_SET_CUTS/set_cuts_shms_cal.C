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

void set_cuts_shms_cal(TString basename) {
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
 outputpdf="plots/"+basename+"_set_cuts_shms_cal.pdf";
   TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_shms_cal_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   //
 static const Int_t nblk=224;
      TCutG *cutg[nblk];
      Double_t tdcmin[nblk];
      Double_t tdcmax[nblk];
  TH1F* tdchist[nblk];
  for (Int_t ip=0;ip<nblk;ip++) {
       TString hname= Form("pcal_%d",ip);
       tdchist[ip] = (TH1F*)fhistroot->Get(hname);
       if (!tdchist[ip]) cout << " no hist = " << hname << endl;
   }
  //
   TCanvas *ctdc;
     ctdc = new TCanvas("ctdc","tdc", 700,700);
     Int_t cont=0;
     for  (Int_t nh=0;nh<2;nh++) {
       cont=0;
       while(cont==0) {
     ctdc->Divide(1,1);
     gPad->SetLogy();
        if (tdchist[nh]) tdchist[nh]->Draw();
     TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
  gPad->Update();
     TString cutname=Form("cut_%d",nh);
    if (!tempg) cout << " no cut" << endl;
    if (tempg)		  {
          cutg[nh]=(TCutG*)(tempg->Clone());
      	  cutg[nh]->SetName(cutname);
	  Int_t npts=cutg[nh]->GetN();
          Double_t ydummy;
	  cutg[nh]->GetPoint(0,tdcmin[nh],ydummy);
	  cutg[nh]->GetPoint(1,tdcmax[nh],ydummy);
    }
 	   TLine *tlow = new TLine(tdcmin[nh],0,tdcmin[nh],tdchist[nh]->GetMaximum()*.9);
	   TLine *thigh = new TLine(tdcmax[nh],0,tdcmax[nh],tdchist[nh]->GetMaximum()*.9);
           tlow->Draw();
           thigh->Draw();
           tlow->SetLineColor(2);
           thigh->SetLineColor(2);
         gPad->Update();
	 cout << " Are gates ok (0=NO,1=YES) " <<endl;
	 cin >> cont;
       }
     }
    //
  const static Int_t numcan=14;
   TCanvas *cplot[numcan];
   Int_t plcnt=0;
     for (Int_t nh=0;nh<numcan;nh++) {
       cplot[nh] = new TCanvas(Form("cplot_%d",nh),Form("plot_%d",nh), 700,700);
       cplot[nh]->Divide(4,4);
        for (Int_t nl=0;nl<16;nl++) {
	  cplot[nh]->cd(nl+1);
	  if (tdchist[plcnt]) {
           tdchist[plcnt]->Draw();
	   TLine *tlow = new TLine(tdcmin[plcnt],0,tdcmin[plcnt],tdchist[plcnt]->GetMaximum()*.9);
	   TLine *thigh = new TLine(tdcmax[plcnt],0,tdcmax[plcnt],tdchist[plcnt]->GetMaximum()*.9);
           tlow->Draw();
           thigh->Draw();
           tlow->SetLineColor(2);
           thigh->SetLineColor(2);
	   plcnt++;
	  }
	}
        if (nh==0) cplot[nh]->Print(outputpdf+"(");
        if (nh==numcan-1) cplot[nh]->Print(outputpdf+")");
	if (nh>0&&nh<numcan-1) cplot[nh]->Print(outputpdf);
     }

     //
 

}
