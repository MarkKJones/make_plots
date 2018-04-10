#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCutG.h>
#include <TMath.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <iomanip>
#include <fstream>
using namespace std;
void set_ytar_yfp_cuts(TString basename,Int_t flag=0) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 //
     TString inputroot;
    inputroot="hist/"+basename+"_hist.root";
      TFile *fhistroot;
     fhistroot =  new TFile(inputroot);
      TH2F *fhist;
      fhist = (TH2F*)fhistroot->Get("hytar_yptar");
      //
   TCanvas *ccompcut; 
   ccompcut= new TCanvas("ccompcut","cut",900,500);
   ccompcut->Divide(1,1);
   ccompcut->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();    
      fhist->SetTitle(basename);
      fhist->SetMinimum(10);
       fhist->Draw("colz");
   //
    TCutG *cutg;
    TCutG *tmpg,*mycutg;
	TCutG *t ;  
    TString outputhist;
    outputhist="hist/"+basename+"_cut.root";
    TFile f(outputhist,"UPDATE");
    Int_t n=0;
        t=(TCutG*)gROOT->FindObject(Form("ytar_yp_cut_%d",n));
	while(t) {
	  cout << " use file = " << outputhist << endl;
                 ccompcut->cd(1);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
                 t=(TCutG*)gROOT->FindObject(Form("ytar_yp_cut_%d",n));
    ccompcut->Update();
		 }
	//
	// ccompcut->WaitPrimitive();
     Int_t nc=0;
       while (nc!=-1) {
        t=(TCutG*)gROOT->FindObject(Form("ytar_yp_cut_%d",nc));
                 ccompcut->cd(1);
	if (t) t->Draw("same");
	if (t) {
             cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	if (nc!=-1) {
    cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    ccompcut->Update();
    //cout << cutg->GetN() << endl;
    tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
    //cout << tmpg->GetN() << endl;
    mycutg=(TCutG*)(tmpg->Clone(Form("ytar_yp_cut_%d",nc)));
    //cout << mycutg->GetN() <<endl;
    nc++;
    mycutg->Write("",TObject::kOverwrite);
    mycutg->Print();
    mycutg->Draw();
       cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	}
	//
	n=0;
	//	 while ((TCutG*)f.Get(Form("mycutg_%d",n))) {
	gDirectory->ls("m");
        t=(TCutG*)gROOT->FindObject(Form("ytar_yp_cut_%d",n));
	while(t) {
		 t->Draw();
		 t->SetLineColor(2);
		 n++;
        t=(TCutG*)gROOT->FindObject(Form("ytar_yp_cut_%d",n));
    ccompcut->Update();
		 }
       //
}
