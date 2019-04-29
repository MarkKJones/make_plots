#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
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

void plot_hist_shms_hodo_raw(TString basename) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
     TString outputpdf;
  TString inputroot;
 static const Int_t nftot=1;
   TFile *fhistroot[nftot];
     inputroot="hist/"+basename+"_shms_hodo_raw_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot[0] =  new TFile(inputroot);
  static const Int_t plnum=4;
 static const Int_t ns=2;
 static const Int_t clmax=200;
 static const Int_t adcsigtypes=2;
 static const Int_t tdcsigtypes=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 static const Int_t npad[plnum]={13,13,14,21};
 const char* sname[ns]={"neg","pos"};
 const char* adcname[adcsigtypes]={"AdcCounter","AdcPulseTime"};
 const char* tdcname[tdcsigtypes]={"TdcCounter","TdcTime"};
 TH1F* h_adcpulsetimediff[plnum][ns][16][nftot];
TH1F* h_adcpulsetime[plnum][ns][16][nftot];
TH1F* h_tdctime[plnum][ns][16][nftot];
  for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t is=0;is<ns;is++) {
 for (Int_t ipad=0;ipad<npad[ip];ipad++) {
 for (Int_t ifn=0;ifn<nftot;ifn++) {
   TString hname = Form("h_adcpulsetimediff_%d_%d_%d",ip,is,ipad);
   h_adcpulsetimediff[ip][is][ipad][ifn] = (TH1F*)fhistroot[ifn]->Get(hname);
   hname= Form("h_adcpulsetime_%d_%d_%d",ip,is,ipad);
   h_adcpulsetime[ip][is][ipad][ifn]  = (TH1F*)fhistroot[ifn]->Get(hname);
   hname = Form("h_tdctime_%d_%d_%d",ip,is,ipad);
   h_tdctime[ip][is][ipad][ifn] =  (TH1F*)fhistroot[ifn]->Get(hname);
 }}}}
 //     
  TCanvas* can[plnum][ns][3];
  const char* tname[3]={"AdcTdcDiffTime","AdcPulseTime","TdcTime"};
  TLine *minline;
  TLine *maxline;
  Double_t gmin[plnum]={70,70,70,70};
  Double_t gmax[plnum]={130,130,130,130};
  for (Int_t ic=2;ic<3;ic++) {
   for (Int_t ip=0;ip<plnum;ip++) {
   for (Int_t is=0;is<ns;is++) {
     can[ip][is][ic] = new TCanvas(Form("can_%s_%s_%d",plname[ip],sname[is],ic),Form("can_%s_%s_%d",plname[ip],sname[is],ic), 1000,700);
     if (ip==0 || ip==1 || ip==2) can[ip][is][ic]->Divide(2,7);
     if (ip==3) can[ip][is][ic]->Divide(2,8);
        for (Int_t ipad=0;ipad<npad[ip];ipad++) {
	  can[ip][is][ic]->cd(ipad+1);
	  if (ic==0) {
	    h_adcpulsetimediff[ip][is][ipad][0]->Draw();	    
	  }
	  if (ic==1) {
	    h_adcpulsetime[ip][is][ipad][0]->Draw();	    
	  }
	  if (ic==2) {
	    h_tdctime[ip][is][ipad][0]->Draw();	    
	    minline = new TLine(gmin[ip],0.,gmin[ip],h_tdctime[ip][is][ipad][0]->GetMaximum());
	    maxline = new TLine(gmax[ip],0.,gmax[ip],h_tdctime[ip][is][ipad][0]->GetMaximum());
				minline->Draw();
				maxline->Draw();
				minline->SetLineColor(2);
				maxline->SetLineColor(2);
	  }
        }
    outputpdf=Form("plots/comp_shms_hodo_raw_%s_%s_%s.pdf",plname[ip],sname[is],tname[ic]);
    can[ip][is][ic]->Print(outputpdf);
   }
   }  
  }
  //
}
