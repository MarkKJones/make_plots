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
#include <TEllipse.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_shms_raster(Int_t nrun=9266,Int_t nevtot=-1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelOffset(-.22,"Y");
 gStyle->SetLabelSize(0.005,"XY");
 gStyle->SetTitleSize(0.06,"xy");
 gStyle->SetPadLeftMargin(0.2);
 gStyle->SetPadRightMargin(0.14);
 //
     TString outputpdf;
     outputpdf=Form("plots/run_%d_raster.pdf",nrun);
  TString inputroot;
   TFile *fhistroot;

   inputroot=Form("online-ROOTfiles/shms_replay_production_mkj_%d_%d.root",nrun,nevtot);
     cout << " infile root = " << inputroot << endl;
   fhistroot =  new TFile(inputroot);
   static const Int_t nh=2; 
   TString h1name[nh]={"pFRAraw_XvsY_notrack","pFRBraw_XvsY_notrack"};
   TString h2name[nh]={"pFRApos_XvsY_notrack","pFRBpos_XvsY_notrack"};
   TString hrebinname[nh]={"pFRAraw_XvsY_r","pFRBraw_XvsY_r"};
    TH2F* hist2d[nh];
    TH2F* histpos2d[nh];
    TH2* histrebin[nh];
    TH2F* histscale[nh];
    Double_t lxrange[nh]={0,0};
    Double_t lyrange[nh]={0,0};
    Double_t lzrange[nh]={.1,.1};
    Double_t hxrange[nh]={1.0,1.0};
    Double_t hyrange[nh]={1.0,1.0};
    Double_t hzrange[nh]={2.2,2.2};
    Double_t nev=172000.; // run 9474
    //nev=226000.; // 9477
    //nev=258000.; // 9478
    //nev=227000.; // 9479
    Double_t depth=0.;
      Double_t xsize=.55/2;
      Double_t ysize=.66/2;
      if (nrun==9474) {
      xsize = (.81-.43)/2;
      ysize = (.867-.45)/2;
      nev=172000.;
      }
      if (nrun==9480) {
      xsize = (.87-.37)/2;
      ysize = (.948-.37)/2;
      nev=516000.;
      }
      if (nrun==9481) {
      xsize = (.85-.39)/2;
      ysize = (.918-.4)/2;
      nev=466000.;
      }
      if (nrun==9482) {
      xsize = (.89-.34)/2;
      ysize = (.968-.34)/2;
      nev=1062696.;
      }
      if (nrun==9479) {
      xsize = (.6)/2;
      ysize = (.6)/2;
      nev=227000.;
      }
      if (nrun==9485) {
      xsize = (.93-.31)/2;
      ysize = (.999-.33)/2;
      nev=451000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9490) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=459000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9491) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=680000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9492) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=400000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9493) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=400000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9494) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=500000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9496) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=400000.;
      lzrange[0]=.6;
      lzrange[1]=.6;
      }
      if (nrun==9497) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=512000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9498) {
      xsize = (.88-.36)/2;
      ysize = (.97-.35)/2;
      nev=695000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9499) {
      xsize = (.92-.33)/2;
      ysize = (.97-.35)/2;
      nev=521000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9504) {
      xsize = (.92-.33)/2;
      ysize = (.97-.35)/2;
      nev=1100000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9505) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=1024000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9505) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=1024000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
      if (nrun==9506) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
      if (nrun==9507) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=500000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9508) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=623000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9509) {
      xsize = (.84-.4 )/2;
      ysize = (.85-.4 )/2;
      nev=511000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9510) {
      xsize = (1.-.3 )/2;
      ysize = (1.-.3 )/2;
      nev=511000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9511) {
      xsize = (.86-.4 )/2;
      ysize = (.89-.4 )/2;
      nev=516000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9514) {
      xsize = (.86-.4 )/2;
      ysize = (.89-.45 )/2;
      nev=619000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9515) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=619000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9516) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=479000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9540) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
       if (nrun==9541) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9542) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9543) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9544) {
      xsize = (.82-.45 )/2;
      ysize = (.85-.47 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9545) {
      xsize = (.82-.45 )/2;
      ysize = (.85-.47 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9546) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.4;
      lzrange[1]=.4;
      }
        if (nrun==9547) {
      xsize = (.89-.43 )/2;
      ysize = (.81-.4 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9548) {
      xsize = (.704-.608 )/2;
      ysize = (.632-.584 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9550) {
      xsize = (.704-.608 )/2;
      ysize = (.632-.584 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9551) {
      xsize = (.808-.456 )/2;
      ysize = (.840-.480 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9552) {
      xsize = (.848-.416 )/2;
      ysize = (.864-.456 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9552) {
      xsize = (.848-.416 )/2;
      ysize = (.864-.456 )/2;
      nev=300000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9569) {
      xsize = (.704-.560 )/2;
      ysize = (.744-.576 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
         if (nrun==9570) {
      xsize = (.736-.536 )/2;
      ysize = (.784-.536 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
         if (nrun==9571) {
      xsize = (.704-.560 )/2;
      ysize = (.752-.568 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9572) {
      xsize = (.688-.576 )/2;
      ysize = (.720-.600 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9574) {
      xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9576) {
      xsize = (.752-.520 )/2;
      ysize = (.808-.512 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9577) {
      xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9578) {
      xsize = (.656-.616 )/2;
      ysize = (.800-.512 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9579) {
      xsize = (.656-.616 )/2;
      ysize = (.680-.640 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9580) {
      xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9581) {
       xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9582) {
       xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9583) {
       xsize = (.768-.496 )/2;
      ysize = (.832-.488 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9585) {
       xsize = (.76-.504 )/2;
      ysize = (.824-.496 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9585) {
       xsize = (.76-.504 )/2;
      ysize = (.824-.496 )/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9586) {
       xsize = 0.24/2;
      ysize = 0.312/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9623) {
       xsize = 0.416/2;
      ysize = 0.424/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9624) {
       xsize = 0.384/2;
      ysize = 0.408/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9625) {
       xsize = 0.4/2;
      ysize = 0.424/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9626) {
       xsize = 0.416/2;
      ysize = 0.424/2;
      nev=940000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9627) {
       xsize = 0.4/2;
      ysize = 0.424/2;
      nev=540000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9628) {
       xsize = 0.384/2;
      ysize = 0.408/2;
      nev=689000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9629) {
       xsize = 0.384/2;
      ysize = 0.408/2;
      nev=741000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9630) {
       xsize = 0.352/2;
      ysize = 0.392/2;
      nev=609000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9631) {
       xsize = 0.368/2;
      ysize = 0.384/2;
      nev=1600000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
          if (nrun==9632) {
       xsize = 0.352/2;
      ysize = 0.376/2;
      nev=554000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
         if (nrun==9633) {
       xsize = 0.352/2;
      ysize = 0.392/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
         if (nrun==9634) {
       xsize = 0.352/2;
      ysize = 0.392/2;
      nev=510000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
         if (nrun==9635) {
       xsize = 0.336/2;
      ysize = 0.392/2;
      nev=542000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9639) {
       xsize = 0.432/2;
      ysize = 0.4/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9656) {
       xsize = 0.4/2;
      ysize = 0.408/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9657) {
       xsize = 0.4/2;
      ysize = 0.408/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9658) {
       xsize = 0.4/2;
      ysize = 0.408/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9659) {
       xsize = 0.464/2;
      ysize = 0.504/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9660) {
       xsize = 0.336/2;
      ysize = 0.376/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9661) {
       xsize = 0.336/2;
      ysize = 0.376/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9662) {
       xsize = 0.336/2;
      ysize = 0.376/2;
      nev=690000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
        if (nrun==9667) {
       xsize = 0.16/2;
      ysize = 0.176/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9668) {
       xsize = 0.16/2;
      ysize = 0.176/2;
      nev=530000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9669) {
       xsize = 0.16/2;
      ysize = 0.184/2;
      nev=514000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9671) {
       xsize = 0.176/2;
      ysize = 0.192/2;
      nev=514000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9672) {
       xsize = 0.16/2;
      ysize = 0.184/2;
      nev=872000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9673) {
       xsize = 0.16/2;
      ysize = 0.184/2;
      nev=510000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9744) {
       xsize = 0.144/2;
      ysize = 0.168/2;
      nev=500000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9746) {
       xsize = 0.144/2;
      ysize = 0.184/2;
      nev=574000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9747) {
       xsize = 0.144/2;
      ysize = 0.184/2;
      nev=535000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9855) {
       xsize = 0.448/2;
      ysize = 0.424/2;
      nev=517000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9856) {
       xsize = 0.448/2;
      ysize = 0.424/2;
      nev=517000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9855) {
       xsize = 0.448/2;
      ysize = 0.424/2;
      nev=517000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==9857) {
       xsize = 0.448/2;
      ysize = 0.424/2;
      nev=517000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==10095) {
       xsize = 0.368/2;
      ysize = 0.376/2;
      nev=391000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
       if (nrun==10152) {
       xsize = 0.368/2;
      ysize = 0.376/2;
      nev=510000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
      if (nrun==10387) {
       xsize = 0.368/2;
      ysize = 0.376/2;
      nev=7787.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
      if (nrun==10381) {
       xsize = 0.368/2;
      ysize = 0.376/2;
      nev=10000.;
      lzrange[0]=.1;
      lzrange[1]=.1;
      }
	  Double_t area = 3.14159*xsize*ysize*(1-depth);
      Double_t numbin=250.;
      Double_t xaxis_size=1.0;
      Double_t bin_size = xaxis_size/numbin;
      Int_t rebin_size=2;
      Double_t cnts_bin = nev/area*bin_size*bin_size*rebin_size*rebin_size;
      cout << " cnts_bin = " << cnts_bin << endl;
for (Int_t ip=0;ip<nh;ip++) {
    // cout <<  hname[ip] << endl;
       hist2d[ip] = (TH2F*)fhistroot->Get(h1name[ip]);
       histpos2d[ip] = (TH2F*)fhistroot->Get(h2name[ip]);
       if (!hist2d[ip]) cout << " no hist = " << h1name[ip] << endl;
       hist2d[ip]->GetXaxis()->SetRangeUser(lxrange[ip],hxrange[ip]);
       hist2d[ip]->GetYaxis()->SetRangeUser(lyrange[ip],hyrange[ip]);
       histrebin[ip] = hist2d[ip]->Rebin2D(rebin_size,rebin_size,hrebinname[ip]);
       histrebin[ip]->GetXaxis()->SetRangeUser(lxrange[ip],hxrange[ip]);
       histrebin[ip]->GetYaxis()->SetRangeUser(lyrange[ip],hyrange[ip]);
       histrebin[ip]->Scale(1./cnts_bin);
       histrebin[ip]->SetMinimum(lzrange[ip]);
       histrebin[ip]->SetMaximum(hzrange[ip]);
  }
  //
 static const Int_t num_rebin = numbin/rebin_size;
 cout << " num_rebin = " << num_rebin << endl;
 TH1D* HistProjX[2][num_rebin];
 TH1D* HistProjY[2][num_rebin];
 Int_t first_xproj[2]={0,0};
 Int_t last_xproj[2]={0,0};
 Int_t first_yproj[2]={0,0};
 Int_t last_yproj[2]={0,0};
 Double_t ymin=0.0;
 Double_t ymax = 1.0;
 Double_t xmin=0.0;
 Double_t xmax = 1.0;
 Double_t xbin_size= (xmax-xmin)/num_rebin;
 Double_t ybin_size= (ymax-ymin)/num_rebin;
 for (Int_t nr=0;nr<2;nr++) {
 for (Int_t ip=1;ip<num_rebin+1;ip++) {
   HistProjX[nr][ip]=histrebin[nr]->ProjectionX(Form("FRA_ProjX_%d_%d",nr,ip),ip,ip);
   Double_t xbinval=ymin+ip*ybin_size;
   if (nr==0) HistProjX[nr][ip]->SetTitle(Form("FRA Y = %4.3f",xbinval));
   if (nr==1) HistProjX[nr][ip]->SetTitle(Form("FRB Y = %4.3f",xbinval));
   HistProjX[nr][ip]->GetYaxis()->SetRangeUser(lzrange[nr],hzrange[nr]);		   
   HistProjX[nr][ip]->GetYaxis()->SetTitle("Ratio ");
   HistProjX[nr][ip]->GetYaxis()->SetLabelSize(0.05);
   HistProjX[nr][ip]->GetYaxis()->SetLabelOffset(0.025);
   //  cout << ip << " integral = " << HistProjX[nr][ip]->Integral() << endl;
   if (HistProjX[nr][ip]->Integral() > 0 && first_xproj[nr]==0) first_xproj[nr]=ip;
   if (HistProjX[nr][ip]->Integral() == 0 && first_xproj[nr]!=0 && last_xproj[nr]==0) last_xproj[nr]=ip-1;
   HistProjY[nr][ip]=histrebin[nr]->ProjectionY(Form("FRA_ProjY_%d_%d",nr,ip),ip,ip);
   Double_t ybinval=xmin+ip*xbin_size;
   if (nr==0) HistProjY[nr][ip]->SetTitle(Form("FRA X = %4.3f",ybinval));
   if (nr==1) HistProjY[nr][ip]->SetTitle(Form("FRB X = %4.3f",ybinval));
   HistProjY[nr][ip]->GetYaxis()->SetRangeUser(lzrange[0],hzrange[0]);		   
   HistProjY[nr][ip]->GetYaxis()->SetTitle("Ratio ");
   HistProjY[nr][ip]->GetYaxis()->SetLabelSize(0.05);
   HistProjY[nr][ip]->GetYaxis()->SetLabelOffset(0.025);
   //   cout << ip << " integral = " << HistProjY[nr][ip]->Integral() << endl;
   if (HistProjY[nr][ip]->Integral() > 0 && first_yproj[nr]==0) first_yproj[nr]=ip;
   if (HistProjY[nr][ip]->Integral() == 0 && first_yproj[nr]!=0 && last_yproj[nr]==0) last_yproj[nr]=ip-1;
 }
 }
 Double_t xproj_xmin[2];
 Double_t xproj_xmax[2];
 Double_t yproj_xmin[2];
 Double_t yproj_xmax[2];
 for (Int_t nr=0;nr<2;nr++) {
   cout << " nr = " << nr << endl;
 xproj_xmin[nr]= (xmin+ first_yproj[nr]*xbin_size)*0.95;
 xproj_xmax[nr]= (xmin+ last_yproj[nr]*xbin_size)*1.05;
 yproj_xmin[nr]= (ymin+ first_xproj[nr]*ybin_size)*0.95;
 yproj_xmax[nr]= (ymin+ last_xproj[nr]*ybin_size)*1.05;
 cout << first_xproj[nr] << " " << last_xproj[nr] << " " <<last_xproj[nr]-first_xproj[nr] << " " << first_yproj[nr] << " " << last_yproj[nr] << " " <<last_yproj[nr]-first_yproj[nr] << endl;
 cout << " x = " << (xmin+ first_yproj[nr]*xbin_size) << " " << (xmin+ last_yproj[nr]*xbin_size) << " " << (xmin+ last_yproj[nr]*xbin_size)-(xmin+ first_yproj[nr]*xbin_size)<< endl;;
 cout << " y = " << (ymin+ first_xproj[nr]*ybin_size) << " " << (ymin+ last_xproj[nr]*ybin_size)<< " " << (ymin+ last_xproj[nr]*ybin_size)-(ymin+ first_xproj[nr]*ybin_size) << endl;;
 }
 //
 ofstream ofile;
 ofile.open("raster.txt");
 for (Int_t nx=0;nx<num_rebin;nx++) {
   Double_t xbinval=xmin+nx*xbin_size - (.93-.31);
   xbinval= xbinval*2./(.93-.31)*2;
  for (Int_t ny=0;ny<num_rebin;ny++) {
    Double_t ybinval=ymin+ny*ybin_size - (1-.33);
     ybinval= ybinval*4./(1-.33);
     Int_t gbin=histrebin[0]->GetBin(nx,ny);
   if (histrebin[0]->GetBinContent(gbin) >0) ofile<< xbinval << " " << ybinval << " " << histrebin[0]->GetBinContent(gbin) << endl;
  }}
 ofile.close();
//
  TCanvas *craster = new TCanvas("craster","Raster FA",700,700);
  craster->Divide(1,1);
      craster->cd(1);
      //hist2d[0]->Draw("colz");
      histrebin[0]->GetXaxis()->SetRangeUser(xproj_xmin[0],xproj_xmax[0]);
      histrebin[0]->GetYaxis()->SetRangeUser(yproj_xmin[0],yproj_xmax[0]);
       histrebin[0]->Draw("colz");
       craster->Print(outputpdf+"(");
//
  TCanvas *craster4 = new TCanvas("craster4","Raster FB",700,700);
  craster4->Divide(1,1);
      craster4->cd(1);
      //hist2d[1]->Draw("colz");
      histrebin[1]->GetXaxis()->SetRangeUser(xproj_xmin[1],xproj_xmax[1]);
      histrebin[1]->GetYaxis()->SetRangeUser(yproj_xmin[1],yproj_xmax[1]);
       histrebin[1]->Draw("colz");
          craster4->Print(outputpdf);
//
  TCanvas *craster2 = new TCanvas("craster2","Proj XA",700,700);
  Int_t npad=1+int(sqrt(last_xproj[0]-first_xproj[0]));
  cout << npad << endl;
  craster2->Divide(npad,npad);
 for (Int_t ip=first_xproj[0];ip<last_xproj[0]+1;ip++) {
   craster2->cd(ip-first_xproj[0]+1);
   HistProjX[0][ip]->GetXaxis()->SetRangeUser(xproj_xmin[0],xproj_xmax[0]);
    HistProjX[0][ip]->Draw();
    }  
     craster2->Print(outputpdf);
//
  TCanvas *craster22 = new TCanvas("craster22","Proj XB",700,700);
  npad=1+int(sqrt(last_xproj[1]-first_xproj[1]));
  cout << npad << endl;
  craster22->Divide(npad,npad);
 for (Int_t ip=first_xproj[1];ip<last_xproj[1]+1;ip++) {
   craster22->cd(ip-first_xproj[1]+1);
     HistProjX[1][ip]->GetXaxis()->SetRangeUser(xproj_xmin[1],xproj_xmax[1]);
     HistProjX[1][ip]->Draw();
    }  
    craster22->Print(outputpdf);
//
  TCanvas *craster3 = new TCanvas("craster3","Proj YA",700,700);
  Int_t npady=1+int(sqrt(last_yproj[0]-first_yproj[0]));
  cout << npady << endl;
  craster3->Divide(npady,npady);
 for (Int_t ip=first_yproj[0];ip<last_yproj[0]+1;ip++) {
   craster3->cd(ip-first_yproj[0]+1);
   HistProjY[0][ip]->GetXaxis()->SetRangeUser(yproj_xmin[0],yproj_xmax[0]);
    HistProjY[0][ip]->Draw();
    }  
     craster3->Print(outputpdf);
//
  TCanvas *craster33 = new TCanvas("craster33","Proj YB",700,700);
  npady=1+int(sqrt(last_yproj[1]-first_yproj[1]));
  cout << npady << endl;
  craster33->Divide(npady,npady);
 for (Int_t ip=first_yproj[1];ip<last_yproj[1]+1;ip++) {
   craster33->cd(ip-first_yproj[1]+1);
    HistProjY[1][ip]->GetXaxis()->SetRangeUser(yproj_xmin[1],yproj_xmax[1]);
   HistProjY[1][ip]->Draw();
    }  
    craster33->Print(outputpdf);
       //
  TCanvas *craster1 = new TCanvas("craster1","Raster FB",700,700);
  craster1->Divide(1,2);
      craster1->cd(1);
      histpos2d[0]->Draw("colz");
      craster1->cd(2);
         histpos2d[1]->Draw("colz");
         craster1->Print(outputpdf+")");
  //
}
