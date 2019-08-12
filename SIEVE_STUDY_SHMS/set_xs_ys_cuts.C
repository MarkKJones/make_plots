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


void set_xs_ys_cuts(TString basename="",Int_t nrun=1813,Int_t nfoil=3,Int_t ndelta=3) {
  
  gROOT->Reset();
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(11);
  gStyle->SetTitleOffset(1.,"Y");
  gStyle->SetTitleOffset(.7,"X");
  gStyle->SetLabelSize(0.04,"XY");
  gStyle->SetTitleSize(0.06,"XY");
  gStyle->SetPadLeftMargin(0.14);
  
   static const Int_t nftot=3;
  static const Int_t ndtot=3;
  static const Int_t nxtot=11;
  static const Int_t nytot=11;
  TString inputroot;
   TFile *fhistroot;
     inputroot="hist/"+basename+"_shms_xs_ys_hist.root";
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);

	TH2F *hXsYs_delta[3][ndtot];
        for (Int_t nc=0;nc<nfoil;nc++) {
           for (Int_t nd=0;nd<ndelta;nd++) {
	     hXsYs_delta[nc][nd]=  (TH2F*)fhistroot->Get(Form("hXsYs_foil_%d_delta_%d",nc,nd));
	   }
	}

 
  TString hname_cut[nftot][ndtot][nxtot][nytot];
  for (Int_t nf=0;nf<nfoil;nf++) {
   for (Int_t nd=0;nd<ndelta;nd++) {
  for (Int_t nx=0;nx<nxtot;nx++) {
 for (Int_t ny=0;ny<nytot;ny++) {
   hname_cut[nf][nd][nx][ny]= Form("xs_%d_ys_%d_cut_foil%d_delta_%d",nx,ny,nf,nd);
 }}}}
 

  
  
  //Set Name of ROOTfile containing polygon Cuts
  TString outCutFile;
  outCutFile=Form("cuts/xs_ys_%d_cut.root", nrun);

  //------------------ Flag 1 ----> Set the Polygon Cut ----------------------

  

    TCanvas *histView_Cut; 
    histView_Cut= new TCanvas("histView_Cut","cut",900,500);
    histView_Cut->Divide(1,1);

    //Loop over focal plane variables
	histView_Cut->cd(1);
	gPad->SetGridx();
	gPad->SetGridy();
	//gPad->SetLogz();    
	
	
    cout << "outCutFile =  " << outCutFile << endl;
    TFile f(outCutFile,"UPDATE");
   
    TCutG*t;


	
	Int_t nloop=0;
        Int_t nx=0;
        Int_t ny=0;
        Int_t nd=0;
        Int_t nf=0;
	while( nloop!=-1) {
	fhistroot->cd();
     hXsYs_delta[nf][nd]->Draw("colz");
                  histView_Cut->Update();
	  f.cd();
    t=(TCutG*)gROOT->FindObject(hname_cut[nf][nd][nx][ny]);
        if(!t){cout << "No cut  =  " <<hname_cut[nf][nd][nx][ny]  << endl;}
 	if(t) {
	  cout << " draw cut = " <<hname_cut[nf][nd][nx][ny]<< endl;
	  f.cd();
 		 t->Draw("same");
		 t->SetLineColor(2);
                  histView_Cut->Update();
	  cout << " nx,ny,nd,nf = " << nx << " " << ny << " " << nd << " " << nf << endl;
          cout <<" nloop new cut (or 0 to look, 1 incr, 2 to set, -1 quit, -10 delete cut) "  << endl;
          cin >> nloop ;
	}
	  if (nloop==1) {
	    nx++;
	    if (nx==nxtot) {
	      nx=0;
	      ny++;
	    }
	    if (ny==nytot) {
	      nx=0;
	      ny=0;
	      nd++;
	    }
	    if (nd==ndelta) {
	      nx=0;
	      ny=0;
	      nd=0;
	      nf++;
	    }
	  }
          if(nloop==-10) {
	      f.cd();
	      f.Delete(hname_cut[nf][nd][nx][ny]);
  	      cout << nx << " delete cut = " <<hname_cut[nf][nd][nx][ny]  << endl;
          }
	  
	  if (nloop!=-1 && nloop!=-10) {
                 t=(TCutG*)gROOT->FindObject(hname_cut[nf][nd][nx][ny]);
         	if(t) {
	         cout << " draw cut = " <<hname_cut[nf][nd][nx][ny]  << endl;
	         f.cd();
 		 t->Draw("same");
		 t->SetLineColor(2);
                  histView_Cut->Update();
	        }
		if (nloop==2) {
            cout << " set cut = "  << nx << " " << ny << " " << nd << " " << nf << endl;
  	fhistroot->cd();
        TCutG*cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
          histView_Cut->Update();
          TCutG*tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
	  TCutG*mycutg=(TCutG*)(tmpg->Clone(hname_cut[nf][nd][nx][ny]));
	  f.cd();
	  mycutg->Write("",TObject::kOverwrite);
	  mycutg->Print();
          mycutg->Draw();
          histView_Cut->Update();
		}
	  }
	} // while

   

}


