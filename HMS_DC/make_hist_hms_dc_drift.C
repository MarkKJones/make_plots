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

void make_hist_hms_dc_drift(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hms_dc_drift_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
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
 Double_t rawtime[plnum][107];
 tsimc->SetBranchAddress("H.dc.residual",resid) ;
 tsimc->SetBranchAddress("H.dc.residualExclPlane",exclresid) ;
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("Ndata.H.dc.%s%s.dist",chname[ip],plname[ip]),&nhits[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.dist",chname[ip],plname[ip]),&dist[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.time",chname[ip],plname[ip]),&time[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.rawtdc",chname[ip],plname[ip]),&rawtime[ip]) ;
   tsimc->SetBranchAddress(Form("H.dc.%s%s.wirenum",chname[ip],plname[ip]),&wire[ip]) ;
 }
 Double_t ntrack;
   tsimc->SetBranchAddress("H.dc.ntrack",&ntrack) ;
 Double_t etotnorm;
   tsimc->SetBranchAddress("H.cal.etotnorm",&etotnorm) ;
 Double_t starttime;
   tsimc->SetBranchAddress("H.hod.starttime",&starttime) ;
 Double_t evtype;
   tsimc->SetBranchAddress("g.evtyp",&evtype) ;
 Double_t trig4_rawtime;
   tsimc->SetBranchAddress("T.hms.pTRIG4_tdcTimeRaw",&trig4_rawtime) ;
 Double_t hSTOF_rawtime;
   tsimc->SetBranchAddress("T.hms.hSTOF_tdcTimeRaw",&hSTOF_rawtime) ;
 Double_t hPRHI_rawtime;
   tsimc->SetBranchAddress("T.hms.hPRHI_tdcTimeRaw",&hPRHI_rawtime) ;
 Double_t hPRLO_rawtime;
   tsimc->SetBranchAddress("T.hms.hPRLO_tdcTimeRaw",&hPRLO_rawtime) ;
 Double_t hELREAL_rawtime;
   tsimc->SetBranchAddress("T.hms.hEL_REAL_tdcTimeRaw",&hELREAL_rawtime) ;
 Double_t hELHI_rawtime;
   tsimc->SetBranchAddress("T.hms.hEL_HI_tdcTimeRaw",&hELHI_rawtime) ;
 Double_t hELLO_rawtime;
   tsimc->SetBranchAddress("T.hms.hEL_LO_tdcTimeRaw",&hELLO_rawtime) ;
 Double_t hELCLEAN_rawtime;
   tsimc->SetBranchAddress("T.hms.hEL_CLEAN_tdcTimeRaw",&hELCLEAN_rawtime) ;
 Double_t hDCREF3_rawtime;
   tsimc->SetBranchAddress("T.hms.hDCREF3_tdcTimeRaw",&hDCREF3_rawtime) ;
 Double_t npesum;
 if (nrun>=1721) {
   cout << " 3pass " << endl;
    tsimc->SetBranchAddress("H.cer.npeSum",&npesum);
 }
 if (nrun<1721) {
   cout << " Onepass  nrun= " << nrun << endl;
    tsimc->SetBranchAddress("H.cer.npeSum",&npesum);
 }
 Double_t xfp;
   tsimc->SetBranchAddress("H.dc.x_fp",&xfp) ;
 Double_t xpfp;
   tsimc->SetBranchAddress("H.dc.xp_fp",&xpfp) ;
 Double_t yfp;
   tsimc->SetBranchAddress("H.dc.y_fp",&yfp) ;
 Double_t ypfp;
   tsimc->SetBranchAddress("H.dc.yp_fp",&ypfp) ;
   // Define histograms
 TH1F* hetotnorm;
 TH2F* hetotnorm_npesum;
 TH1F* hwire[plnum];
 TH1F* hdist[plnum];
 TH1F* htime[plnum];
 TH1F* hrawtime[plnum];
 TH1F* hresid[plnum];
 TH1F* hexclresid[plnum];
 TH2F* hwire_resid[plnum];
 TH2F* hwire_time[plnum];
 TH2F* hwire_dist[plnum];
 TH2F* htime_dist[plnum];
 TH2F* htime_trig4_rawtime[plnum];
 TH2F* htime_ELREAL_rawtime[plnum];
 TH2F* htime_ELREAL_rawtime2[plnum];
 TH2F* htime_ELHI_rawtime[plnum];
 TH2F* htime_ELLO_rawtime[plnum];
 TH2F* htime_ELCLEAN_rawtime[plnum];
 TH2F* htime_DCREF3_rawtime[plnum];
 TH2F* hstarttime_resid[plnum];
 TString temp;
 TH1F* hstarttime;
   temp=Form("Run %d ; Hod Starttime ; Counts",nrun);
   hstarttime = new TH1F("hstarttime",temp,200,0,100);
   HList.Add(hstarttime);
   temp=Form("Run %d ; Etotnorm  ; Counts",nrun);
   hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   HList.Add(hetotnorm);
   if (nrun > 1721) {
   temp=Form("Run %d ; Etotnorm  ; NG npeSum",nrun);
   hetotnorm_npesum = new TH2F("hetotnorm_npesum",temp,150,0,1.5,250,0,50);
   } else {
   temp=Form("Run %d ; Etotnorm  ; HG npeSum",nrun);
   hetotnorm_npesum = new TH2F("hetotnorm_npesum",temp,150,0,1.5,60,0,15);
   }
   HList.Add(hetotnorm_npesum);
   for (Int_t ip=0;ip<plnum;ip++) {
   temp=Form("Run %d ; Plane  %s%s Resid ; Counts",nrun,chname[ip],plname[ip]);
   hresid[ip] = new TH1F(Form("resid_%s%s",chname[ip],plname[ip]),temp,200,-.5,.5);
   HList.Add(hresid[ip]);
   temp=Form("Run %d ; Plane  %s%s Excl Resid ; Counts",nrun,chname[ip],plname[ip]);
   hexclresid[ip] = new TH1F(Form("exclresid_%s%s",chname[ip],plname[ip]),temp,200,-.5,.5);
   HList.Add(hexclresid[ip]);
   temp=Form("Run %d ; Plane  %s%s Wire  ; Counts",nrun,chname[ip],plname[ip]);
   hwire[ip] = new TH1F(Form("wire_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip]);
   HList.Add(hwire[ip]);
   temp=Form("Run %d ; Plane  %s%s Dist ; Counts",nrun,chname[ip],plname[ip]);
   hdist[ip] = new TH1F(Form("dist_%s%s",chname[ip],plname[ip]),temp,120,0.,.6);
   HList.Add(hdist[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; Counts",nrun,chname[ip],plname[ip]);
   htime[ip] = new TH1F(Form("time_%s%s",chname[ip],plname[ip]),temp,400,-150,250);
   HList.Add(htime[ip]);
   temp=Form("Run %d ; Plane  %s%s Raw Time ; Counts",nrun,chname[ip],plname[ip]);
   hrawtime[ip] = new TH1F(Form("rawtime_%s%s",chname[ip],plname[ip]),temp,250,-15000,-10000);
   HList.Add(hrawtime[ip]);
  temp=Form("Run %d ; Plane  %s%s Wire ; Resid",nrun,chname[ip],plname[ip]);
   hwire_resid[ip] = new TH2F(Form("wire_resid_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],200,-.2,.2);
   HList.Add(hwire_resid[ip]);
  temp=Form("Run %d ; Plane  %s%s Start time ; Resid",nrun,chname[ip],plname[ip]);
  hstarttime_resid[ip] = new TH2F(Form("starttime_resid_%s%s",chname[ip],plname[ip]),temp,140,0,70,200,-.2,.2);
   HList.Add(hstarttime_resid[ip]);
   temp=Form("Run %d ; Plane  %s%s Wire ; Time",nrun,chname[ip],plname[ip]);
   hwire_time[ip] = new TH2F(Form("wire_time_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],350,-150,200);
   HList.Add(hwire_time[ip]);
   temp=Form("Run %d ; Plane  %s%s Wire ; Dist",nrun,chname[ip],plname[ip]);
   hwire_dist[ip] = new TH2F(Form("wire_dist_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],120,0,.6);
   temp=Form("Run %d ; Plane  %s%s Time ; Dist",nrun,chname[ip],plname[ip]);
   htime_dist[ip] = new TH2F(Form("time_dist_%s%s",chname[ip],plname[ip]),temp,200,-50,250,120,0,.6);
   HList.Add(htime_dist[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; Trig4 raw time ",nrun,chname[ip],plname[ip]);
   htime_trig4_rawtime[ip] = new TH2F(Form("time_trig4_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,3500,6500);
   HList.Add(htime_trig4_rawtime[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; ELREAL raw time (ELHI==0) ",nrun,chname[ip],plname[ip]);
   htime_ELREAL_rawtime[ip] = new TH2F(Form("time_ELREAL_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,3000,5000);
   HList.Add(htime_ELREAL_rawtime[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; ELREAL raw time (ELLO==0)) ",nrun,chname[ip],plname[ip]);
   htime_ELREAL_rawtime2[ip] = new TH2F(Form("time_ELREAL_rawtime2_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,3000,5000);
   HList.Add(htime_ELREAL_rawtime2[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; ELHI raw time ",nrun,chname[ip],plname[ip]);
   htime_ELHI_rawtime[ip] = new TH2F(Form("time_ELHI_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,2000,5000);
   HList.Add(htime_ELHI_rawtime[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; ELLO raw time ",nrun,chname[ip],plname[ip]);
   htime_ELLO_rawtime[ip] = new TH2F(Form("time_ELLO_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,2000,5000);
   HList.Add(htime_ELLO_rawtime[ip]);
    temp=Form("Run %d ; Plane  %s%s Time ; ELCLEAN raw time ",nrun,chname[ip],plname[ip]);
   htime_ELCLEAN_rawtime[ip] = new TH2F(Form("time_ELCLEAN_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,2000,5000);
   HList.Add(htime_ELCLEAN_rawtime[ip]);
   temp=Form("Run %d ; Plane  %s%s Time ; DCREF3 raw time ",nrun,chname[ip],plname[ip]);
   htime_DCREF3_rawtime[ip] = new TH2F(Form("time_DCREF3_rawtime_%s%s",chname[ip],plname[ip]),temp,400,-500,250,200,15000,18000);
   HList.Add(htime_DCREF3_rawtime[ip]);
 }
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		hetotnorm->Fill(etotnorm);
		hetotnorm_npesum->Fill(etotnorm,npesum);
	     		if (evtype==2 &&  starttime>0 && starttime < 100) {
                  hstarttime->Fill(starttime);
                  for (Int_t ip=0;ip<plnum;ip++) {
                  hresid[ip]->Fill(resid[ip]);
                  hexclresid[ip]->Fill(exclresid[ip]);
		  hstarttime_resid[ip]->Fill(starttime,resid[ip]);
  		   for (Int_t id=0;id<nhits[ip];id++) {
		   if (id==0) hwire_resid[ip]->Fill(wire[ip][id],resid[ip]);
                   hwire[ip]->Fill(wire[ip][id]);
                   hdist[ip]->Fill(dist[ip][id]);
                   htime[ip]->Fill(time[ip][id]);
                   hrawtime[ip]->Fill(rawtime[ip][id]);
                   hwire_time[ip]->Fill(wire[ip][id],time[ip][id]);
                   hwire_dist[ip]->Fill(wire[ip][id],dist[ip][id]);
                   htime_dist[ip]->Fill(time[ip][id],dist[ip][id]);
                    htime_trig4_rawtime[ip]->Fill(time[ip][id],trig4_rawtime);
		    if( hELHI_rawtime==0)  htime_ELREAL_rawtime[ip]->Fill(time[ip][id],hELREAL_rawtime);
		    if( hELLO_rawtime==0)  htime_ELREAL_rawtime2[ip]->Fill(time[ip][id],hELREAL_rawtime);
                   htime_ELHI_rawtime[ip]->Fill(time[ip][id],hELHI_rawtime);
                   htime_ELLO_rawtime[ip]->Fill(time[ip][id],hELLO_rawtime);
                   htime_ELCLEAN_rawtime[ip]->Fill(time[ip][id],hELCLEAN_rawtime);
                   htime_DCREF3_rawtime[ip]->Fill(time[ip][id],hDCREF3_rawtime);
		   }
	          }
	      	}// if pid /event selection
	}
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
