#include <TString.h>
#include <TLegend.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TLine.h>
#include <TCutG.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TMarker.h>
#include <TProfile.h>
#include <TObjArray.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void set_xfp_cuts(Int_t nrun,Int_t dflag=1) {
 gROOT->Reset();
 gStyle->SetOptStat(0);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.14);
 TString basename = Form("hist/%d_hist_flag4.root",nrun);
     TString outputpdf;
 //
     outputpdf=Form("plots/%d_xfp_cuts.pdf",nrun);
     //
  
  string delta_cut_filename = Form("delta_cuts/%d.dat",nrun);

  TString inputroot;
   TFile *fhistroot;
   TString hname;
   fhistroot =  new TFile(basename);
     static const Int_t ntypes=5;
     TString namtype[ntypes] = {"hypfp_cut","hyfp_cut","hxpfp_cut","hxfp_cut","hdelta_cut"};
     static const Int_t xfpcut_num=11;
     static const Int_t yfpcut_num=11;
     Double_t xlo[yfpcut_num][xfpcut_num];
     Double_t xmid[yfpcut_num][xfpcut_num];
     Double_t xhi[yfpcut_num][xfpcut_num];
       for(Int_t nx=0;nx<xfpcut_num;nx++) {
	 for(Int_t ny=0;ny<yfpcut_num;ny++) {
	   xlo[ny][nx]=0.0;
	   xmid[ny][nx]=0.0;
	   xhi[ny][nx]=0.0;
       }}
     //
     ifstream delta_cut_file;
     string s;
     Int_t nyt,nxt;
     float v1,v2,v3;
     delta_cut_file.open(delta_cut_filename.c_str(),ios::in);
     if (delta_cut_file.is_open()) {
     while( getline(delta_cut_file,s)) {
       //cout << " scan = " << s << endl;
       sscanf(s.c_str(),"%d %d %f %f %f",&nxt,&nyt,&v1,&v2,&v3);
       xlo[nyt][nxt]=v1;
       xmid[nyt][nxt]=v2;
       xhi[nyt][nxt]=v3;
        cout << nyt << " " << nxt  << " " << xlo[nyt][nxt] << endl;
     }
     } else {
       cout << " no file = " << delta_cut_filename << endl;
     }
     //
     Int_t ntot_xcut=0;
     Int_t ntot_ycut=0;
     TH1F* h_cut[ntypes][yfpcut_num][xfpcut_num];
     TH2F* h_2dcut[3][yfpcut_num][xfpcut_num];
     for(Int_t ny=0;ny<yfpcut_num;ny++) {
     for(Int_t nx=0;nx<xfpcut_num;nx++) {
       for(Int_t nt=0;nt<ntypes;nt++) {
       hname= namtype[nt]+Form("_%d_%d",nx,ny);
       h_cut[nt][ny][nx] = (TH1F*)fhistroot->Get(hname);
       if (!h_cut[nt][ny][nx]) cout << " no hist = " << hname << endl;       
       }
       hname =  "hxpfp_ypfp_cut";
       hname= hname+Form("_%d_%d",nx,ny);
       h_2dcut[0][ny][nx] = (TH2F*)fhistroot->Get(hname);
       if (!h_2dcut[0][ny][nx]) cout << " no hist = " << hname << endl;       
       hname =  "hyfp_ypfp_cut";
       hname= hname+Form("_%d_%d",nx,ny);
       h_2dcut[1][ny][nx] = (TH2F*)fhistroot->Get(hname);
       if (!h_2dcut[1][ny][nx]) cout << " no hist = " << hname << endl;       
       hname =  "hxfp_xpfp_cut";
       hname= hname+Form("_%d_%d",nx,ny);
       h_2dcut[2][ny][nx] = (TH2F*)fhistroot->Get(hname);
       if (!h_2dcut[2][ny][nx]) cout << " no hist = " << hname << endl;       
     }}
   //
     if (dflag==1) {
     Int_t flag = 1;
     Int_t ans=1;
     Double_t ydummy;
     TLine *tl;
     TLine *tl2;
     TLine *tl3;
     TMarker *marker;
       TCanvas *cplot = new TCanvas("cplot","xfp",700,500);
       cplot->Divide(1,1);
       cplot->ToggleToolBar();
       cplot->ToggleEventStatus();
       for(Int_t ny=0;ny<yfpcut_num;ny++) {
       for(Int_t nx=0;nx<xfpcut_num;nx++) {
	 if (h_cut[3][ny][nx] && h_cut[3][ny][nx]->Integral()>20) {
	   Double_t max = h_cut[3][ny][nx]->GetBinCenter(h_cut[3][ny][nx]->GetMaximumBin());
	   
           cplot->cd(1);
	   h_cut[3][ny][nx]->SetAxisRange(max-3,max+10,"x");
           h_cut[3][ny][nx]->Draw();
	     tl = new TLine(xlo[ny][nx],0,xlo[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl->Draw();
	     tl->SetLineColor(2);
	     tl2 = new TLine(xmid[ny][nx],0,xmid[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl2->Draw();
	     tl2->SetLineColor(2);
	     tl3 = new TLine(xhi[ny][nx],0,xhi[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl3->Draw();
	     tl3->SetLineColor(2);
       cplot->Update();
       ans=-1;
       cout << "set =10 to skip, Set new gates = -1 , Gates are ok = 1 " << xlo[ny][nx] << " " <<xmid[ny][nx] << " " << xhi[ny][nx] << " " << endl;
	     cin >> ans ;
	   while (ans==-1) {
            cplot->cd(1);
	    //  TCutG *tempg = (TCutG*) gPad->WaitPrimitive("CUTG","CutG");
            cout << " set marker 1 with mouse " << endl;
            gPad->Modified();
            gPad->Update();
	    marker = (TMarker*) gPad->WaitPrimitive("TMarker","Marker");
            gPad->Modified();
            gPad->Update();
              if (marker) {
	       /*tempg->GetPoint(0,xlo[ny][nx],ydummy);
             tempg->GetPoint(1,xmid[ny][nx],ydummy);
             tempg->GetPoint(2,xhi[ny][nx],ydummy);
	       */
	     xlo[ny][nx]=marker->GetX();
	     tl = new TLine(xlo[ny][nx],0,xlo[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl->Draw();
	     tl->SetLineColor(2);
	     }
	     delete marker;
            gPad->Modified();
            gPad->Update();
            cout << " set marker 2 with mouse " << endl;
 	    marker = (TMarker*) gPad->WaitPrimitive("TMarker","Marker");
            gPad->Modified();
            gPad->Update();
             if (marker) {
             cplot->Update();
	     xmid[ny][nx]=marker->GetX();
 	     tl2 = new TLine(xmid[ny][nx],0,xmid[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl2->Draw();
	     tl2->SetLineColor(2);
	     } else {
	       cout << " No marker 2" << endl;
	     }
	     delete marker;
            gPad->Modified();
            gPad->Update();
            cout << " set marker 3 with mouse " << endl;
 	    marker = (TMarker*) gPad->WaitPrimitive("TMarker","Marker");
            gPad->Modified();
            gPad->Update();
             if (marker) {
	     xhi[ny][nx]=marker->GetX();
	     tl3 = new TLine(xhi[ny][nx],0,xhi[ny][nx],h_cut[3][ny][nx]->GetMaximum());
	     tl3->Draw();
	     tl3->SetLineColor(2);
             cplot->Update();
	     } else {
	       cout << " No marker 3" << endl;
	     }
	     delete marker;
            gPad->Modified();
            gPad->Update();
	     cout << "Set new gates = -1 , Gates are ok = 1 " << " " << xlo[ny][nx]<< " " << xmid[ny][nx]<< " " << xhi[ny][nx]<< endl;
	     cin >> ans ;

	   }
	   if (ans==10) break;
	 }
       } // ny loop
       //
       }
     }
     if (dflag==0) {
   TCanvas *cresid[xfpcut_num];
       for(Int_t nx=0;nx<xfpcut_num;nx++) {
        cresid[nx] = new TCanvas(Form("cresid_%d",nx),Form("xfp_%d",nx),1000,700);
        cresid[nx]->Divide(4,3);
       Double_t cnt=1;
       for(Int_t nny=0;nny<yfpcut_num;nny++) {
         cresid[nx]->cd(cnt++);
         if (h_cut[3][nny][nx]) {
	   Double_t max = h_cut[3][nny][nx]->GetBinCenter(h_cut[3][nny][nx]->GetMaximumBin());
	     h_cut[3][nny][nx]->SetAxisRange(max-2,max+2,"x");
	     h_cut[3][nny][nx]->Draw();
	     TLine *ttl = new TLine(xlo[nny][nx],0,xlo[nny][nx],h_cut[3][nny][nx]->GetMaximum());
	     ttl->Draw();
	     ttl->SetLineColor(2);
	     TLine *ttl2 = new TLine(xmid[nny][nx],0,xmid[nny][nx],h_cut[3][nny][nx]->GetMaximum());
	     ttl2->Draw();
	     ttl2->SetLineColor(2);
	     TLine *ttl3 = new TLine(xhi[nny][nx],0,xhi[nny][nx],h_cut[3][nny][nx]->GetMaximum());
	     ttl3->Draw();
	     ttl3->SetLineColor(2);
              gPad->Update();
         }
       }
    if (nx==0) cresid[nx]->Print(outputpdf+"(");
    if (nx>0&&nx<10) cresid[nx]->Print(outputpdf);
    if (nx==10) cresid[nx]->Print(outputpdf+")");
       }
     }
       //
     if (dflag==1) {
       delta_cut_filename = Form("delta_cuts/new_%d.dat",nrun);
     ofstream delta_cut_outfile;
     delta_cut_outfile.open(delta_cut_filename.c_str(),ios::out);
       
       for(Int_t nx=0;nx<xfpcut_num;nx++) {
       for(Int_t ny=0;ny<yfpcut_num;ny++) {
	 delta_cut_outfile << nx << " " << ny << " " << xlo[ny][nx] << " " << xmid[ny][nx] << " " << xhi[ny][nx] << endl;
       }}
       delta_cut_outfile.close();
     }
       //
}
