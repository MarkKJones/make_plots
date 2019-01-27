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
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_hms_trackeff(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1000011);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_trackeff_hist.root";
 TObjArray HList(0);
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t npesum;
    tsimc->SetBranchAddress("H.cer.npeSum",&npesum);
 Double_t etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("H.dc.ntrack",&ntrack) ;
 Double_t xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp) ;
 Double_t xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp) ;
 Double_t yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp) ;
 Double_t goodscinhit;
   tsimc->SetBranchAddress("H.hod.goodscinhit",&goodscinhit) ;
 Double_t ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp) ;
 Double_t evtype;
   tsimc->SetBranchAddress("g.evtyp",&evtype) ;
// Define branches
 static const Int_t plnum=12;
 const char* plname[plnum]={"u1","u2","x1","x2","v1","v2","v2","v1","x2","x1","u2","u1"};
 const char* chname[plnum]={"1","1","1","1","1","1","2","2","2","2","2","2"};
 static const Int_t nwires[plnum]={107, 107, 79, 79, 107, 107,107, 107, 79, 79, 107, 107};
 Int_t nhits[plnum];
 Double_t resid[plnum];
 Double_t exclresid[plnum];
 Double_t dist[plnum][107];
 Double_t wire[plnum][107];
 Double_t time[plnum][107];
 tsimc->SetBranchAddress("H.dc.residual",resid) ;
 tsimc->SetBranchAddress("H.dc.residualExclPlane",exclresid) ;
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("Ndata.H.dc.%s%s.dist",chname[ip],plname[ip]),&nhits[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.dist",chname[ip],plname[ip]),&dist[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.time",chname[ip],plname[ip]),&time[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.wirenum",chname[ip],plname[ip]),&wire[ip]) ;
 }
   //
 TH1F* hetotnorm;
 TH1F* hscinshould;
 TH1F* hscindid;
 TH1F* hscinshould_ev4;
 TH1F* hscindid_ev4;
 TH1F* hscinshould_ev2;
 TH1F* hscindid_ev2;
 TH2F* hetotnorm_npesum;
 TH1F* hdist[plnum];
 TH1F* htime[plnum];
 TH1F* hdist_ev4[plnum];
 TH1F* htime_ev4[plnum];
 TH1F* hdist_ev2[plnum];
 TH1F* htime_ev2[plnum];
 TString temp;
   for (Int_t ip=0;ip<plnum;ip++) {
 
   temp=Form("Run %d ; Plane  %s%s Dist ; Counts",nrun,chname[ip],plname[ip]);
   hdist[ip] = new TH1F(Form("dist_%s%s",chname[ip],plname[ip]),temp,120,0.,.6);
   HList.Add(hdist[ip]);
 
   temp=Form("Run %d ; Plane  %s%s Time ; Counts",nrun,chname[ip],plname[ip]);
   htime[ip] = new TH1F(Form("time_%s%s",chname[ip],plname[ip]),temp,200,-50,250);
   HList.Add(htime[ip]);
 
   temp=Form("Run %d ; Plane  %s%s Dist_ev4 ; Counts",nrun,chname[ip],plname[ip]);
   hdist_ev4[ip] = new TH1F(Form("dist_%s%s_ev4",chname[ip],plname[ip]),temp,120,0.,.6);
   HList.Add(hdist_ev4[ip]);
 
   temp=Form("Run %d ; Plane  %s%s Time_ev4 ; Counts",nrun,chname[ip],plname[ip]);
   htime_ev4[ip] = new TH1F(Form("time_%s%s_ev4",chname[ip],plname[ip]),temp,200,-50,250);
   HList.Add(htime_ev4[ip]);
 
   temp=Form("Run %d ; Plane  %s%s Dist_ev2 ; Counts",nrun,chname[ip],plname[ip]);
   hdist_ev2[ip] = new TH1F(Form("dist_%s%s_ev2",chname[ip],plname[ip]),temp,120,0.,.6);
   HList.Add(hdist_ev2[ip]);
 
   temp=Form("Run %d ; Plane  %s%s Time_ev2 ; Counts",nrun,chname[ip],plname[ip]);
   htime_ev2[ip] = new TH1F(Form("time_%s%s_ev2",chname[ip],plname[ip]),temp,200,-50,250);
   HList.Add(htime_ev2[ip]);

   }

   temp=Form("Run %d ; Etotnorm  ; Counts",nrun);
   hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   temp=Form("Run %d ; Etotnorm  ; NG npeSum",nrun);
   hetotnorm_npesum = new TH2F("hetotnorm_npesum",temp,150,0,1.5,250,0,50);
   HList.Add(hetotnorm_npesum);
   temp=Form("Run %d ; Scin Should ; Counts",nrun);
   hscinshould = new TH1F("hscinshould",temp,5,0,5);
   HList.Add(hscinshould);
   temp=Form("Run %d ; Scin Did ; Counts",nrun);
   hscindid = new TH1F("hscindid",temp,5,0,5);
   HList.Add(hscindid);

   temp=Form("Run %d ; Scin Should_ev4 ; Counts",nrun);
   hscinshould_ev4 = new TH1F("hscinshould_ev4",temp,5,0,5);
   HList.Add(hscinshould_ev4);
   temp=Form("Run %d ; Scin Did_ev4 ; Counts",nrun);
   hscindid_ev4 = new TH1F("hscindid_ev4",temp,5,0,5);
   HList.Add(hscindid_ev4);
   temp=Form("Run %d ; Xfp ev4 ; Counts",nrun);
   hxfp_ev4 = new TH1F("hxfp_ev4",temp,100,-50,50);
   HList.Add(hxfp_ev4);

   temp=Form("Run %d ; Scin Should_ev2 ; Counts",nrun);
   hscinshould_ev2 = new TH1F("hscinshould_ev2",temp,5,0,5);
   HList.Add(hscinshould_ev2);
   temp=Form("Run %d ; Scin Did_ev2 ; Counts",nrun);
   hscindid_ev2 = new TH1F("hscindid_ev2",temp,5,0,5);
   HList.Add(hscindid_ev2);
   temp=Form("Run %d ; Xfp ev2 ; Counts",nrun);
   hxfp_ev2 = new TH1F("hxfp_ev2",temp,100,-50,50);
   HList.Add(hxfp_ev2);
// loop over entries
 Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
 		hetotnorm->Fill(etotnorm);
		hetotnorm_npesum->Fill(etotnorm,npesum);
		if (etotnorm>.6 && npesum>2.) {
		  hscinshould->Fill(goodscinhit);
		  if (ntrack>0) hscindid->Fill(goodscinhit);
		  if (evtype==4) {
		  hscinshould_ev4->Fill(goodscinhit);
		  if (ntrack>0) hscindid_ev4->Fill(goodscinhit);
		  if (goodscinhit==1) hxfp_ev4->Fill(xfp);
		  }
		  if (evtype==2) {
		  hscinshould_ev2->Fill(goodscinhit);
		  if (ntrack>0) hscindid_ev2->Fill(goodscinhit);
		  if (goodscinhit==1) hxfp_ev2->Fill(xfp);
		  }
                 for (Int_t ip=0;ip<plnum;ip++) {
		   for (Int_t id=0;id<nhits[ip];id++) {
                   hdist[ip]->Fill(dist[ip][id]);
                   htime[ip]->Fill(time[ip][id]);
		  if (evtype==2) {
                   hdist_ev2[ip]->Fill(dist[ip][id]);
                   htime_ev2[ip]->Fill(time[ip][id]);
		  }
		  if (evtype==4) {
                   hdist_ev4[ip]->Fill(dist[ip][id]);
                   htime_ev4[ip]->Fill(time[ip][id]);
		  }
		   }		 }		  
 		} // if pid /event selection
	}
	//
	cout << "All  Scin should = " << hscinshould->Integral(1,1) << " Scin did = " << hscindid->Integral(1,1) << " rat = "  << hscindid->Integral(1,1)/hscinshould->Integral(1,1) << endl;
	cout << "Ev4  Scin should = " << hscinshould_ev4->Integral(1,1) << " Scin did = " << hscindid_ev4->Integral(1,1) << " rat = "  << hscindid_ev4->Integral(1,1)/hscinshould_ev4->Integral(1,1) << endl;
	cout << "Ev2  Scin should = " << hscinshould_ev2->Integral(1,1) << " Scin did = " << hscindid_ev2->Integral(1,1) << " rat = "  << hscindid_ev2->Integral(1,1)/hscinshould_ev2->Integral(1,1) << endl;
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
