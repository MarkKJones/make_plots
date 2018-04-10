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
void plot_fp(TString basename,Int_t flag=0) {
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
     TH2F *fhist[6];
     TH1F *fdelta[10];
    TString hname[6] = {"h_xfp_yfp","h_xfp_ypfp","h_xfp_xpfp","h_yfp_xpfp","h_yfp_ypfp","h_xpfp_ypfp"};
     fhistroot =  new TFile(inputroot);
     fhistroot->ls("d");
     for (Int_t nc=0;nc<6;nc++) {
       fhist[nc] = (TH2F*)fhistroot->Get(hname[nc]);
       cout << inputroot << " " << nc << " " << hname[nc] << endl;
     }
    //
   if (flag==0) {
   TCanvas *ccomp[6];
   for (Int_t nc=0;nc<6;nc++) {
   TString xlabel= hname[nc];
   ccomp[nc] = new TCanvas(Form("ccomp_%d",nc),xlabel,700,800);
   ccomp[nc]->Divide(1,1);
   ccomp[nc]->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
   //gPad->SetLogz();
   cout << " plot = " << hname[nc] << endl;
      fhist[nc]->SetTitle(basename);
      fhist[nc]->SetMinimum(10);
       fhist[nc]->Draw("colz");
   }
   }
   if (flag==2) {
    TString outputhist;
   outputhist="hist/"+basename+"_cut.root";
   TFile f(outputhist,"UPDATE");
   TCanvas *ccompcut; 
   ccompcut= new TCanvas("ccompcut","cut",900,500);
   ccompcut->Divide(1,1);
   ccompcut->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();    
      fhist[2]->SetTitle(basename);
      fhist[2]->SetMinimum(10);
       fhist[2]->Draw("colz");
	TCutG *t ;
        Int_t n=0;  
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
	while(t) {
              t->Write("",TObject::kOverwrite);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
	//  ccompcut->Update();
		 }
   TCanvas *ccompcut1; 
   ccompcut1= new TCanvas("ccompcut1","cut1",900,500);
   ccompcut1->Divide(1,1);
   ccompcut1->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();
    
      fhist[4]->SetTitle(basename);
      fhist[4]->SetMinimum(10);
       fhist[4]->Draw("colz");
        n=0;  
        t=(TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",n));
	while(t) {
              t->Write("",TObject::kOverwrite);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
        t=(TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",n));
        ccompcut->Update();
		 }
    }
   //
   if (flag==1) {
   TCanvas *ccompcut; 
   ccompcut= new TCanvas("ccompcut","cut",900,500);
   ccompcut->Divide(1,1);
   ccompcut->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();    
      fhist[2]->SetTitle(basename);
      fhist[2]->SetMinimum(10);
       fhist[2]->Draw("colz");
    //
    TCutG* cutg;
    TCutG *tmpg,*mycutg;
	TCutG *t ;  
    TString outputhist;
   outputhist="hist/"+basename+"_cut.root";
    TFile f(outputhist,"UPDATE");
    Int_t n=0;
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
	while(t) {
                 ccompcut->cd(1);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
                 t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
    ccompcut->Update();
		 }
	//
	// ccompcut->WaitPrimitive();
     Int_t nc=0;
       while (nc!=-1) {
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",nc));
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
    mycutg=(TCutG*)(tmpg->Clone(Form("xpfp_xfp_cut_%d",nc)));
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
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
	while(t) {
		 t->Draw();
		 t->SetLineColor(2);
		 n++;
        t=(TCutG*)gROOT->FindObject(Form("xpfp_xfp_cut_%d",n));
    ccompcut->Update();
		 }
	//
   TCanvas *ccompcut1; 
   ccompcut1= new TCanvas("ccompcut1","cut1",900,500);
   ccompcut1->Divide(1,1);
   ccompcut1->cd(1);
   gPad->SetGridx();
   gPad->SetGridy();
    gPad->SetLogz();
    
      fhist[4]->SetTitle(basename);
      fhist[4]->SetMinimum(10);
       fhist[4]->Draw("colz");
       n=0;
        t=(TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",n));
	while(t) {
                 ccompcut1->cd(1);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
                 t=(TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",n));
    ccompcut1->Update();
		 }
    nc=0;
        while (nc!=-1) {
        t=(TCutG*)gROOT->FindObject(Form("ypfp_yfp_cut_%d",nc));
                 ccompcut1->cd(1);
	if (t) t->Draw("same");
	if (t) {
             cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	if (nc!=-1) {
    cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    ccompcut1->Update();
    tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
    mycutg=(TCutG*)(tmpg->Clone(Form("ypfp_yfp_cut_%d",nc)));
    nc++;
    mycutg->Write("",TObject::kOverwrite);
    mycutg->Print();
    mycutg->Draw("same");
       cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	}
	 //
   }
      //
}


