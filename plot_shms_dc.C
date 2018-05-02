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

void plot_shms_dc(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_hist.root";
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
 Double_t dist[plnum][107];
 Double_t wire[plnum][107];
 Double_t time[plnum][107];
 tsimc->SetBranchAddress("P.dc.residual",resid) ;
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("Ndata.P.dc.%s%s.dist",chname[ip],plname[ip]),&nhits[ip]) ;
   tsimc->SetBranchAddress(Form("P.dc.%s%s.dist",chname[ip],plname[ip]),&dist[ip]) ;
   tsimc->SetBranchAddress(Form("P.dc.%s%s.time",chname[ip],plname[ip]),&time[ip]) ;
   tsimc->SetBranchAddress(Form("P.dc.%s%s.wirenum",chname[ip],plname[ip]),&wire[ip]) ;
 }
 Double_t etotnorm;
   tsimc->SetBranchAddress("P.cal.etotnorm",&etotnorm) ;
 Double_t sp1_id;
   tsimc->SetBranchAddress("P.dc.sp1_id",&sp1_id) ;
 Double_t ntrack;
   tsimc->SetBranchAddress("P.dc.ntrack",&ntrack) ;
 Double_t ndchit;
   tsimc->SetBranchAddress("P.dc.nhit",&ndchit) ;
 Double_t sp2_id;
   tsimc->SetBranchAddress("P.dc.sp2_id",&sp2_id) ;
 Double_t xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp) ;
 Double_t xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp) ;
 Double_t yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp) ;
 Double_t ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp) ;
 Int_t nsp[2];
 Double_t stub_x[2][100];
 Double_t stub_xp[2][100];
 Double_t stub_y[2][100];
 Double_t stub_yp[2][100];
 const char* chname2[2]={"1","2"};
 for (Int_t ich=0;ich<2;ich++) {
   tsimc->SetBranchAddress(Form("Ndata.P.dc.Ch%s.stub_x",chname2[ich]),&nsp[ich]) ;
   tsimc->SetBranchAddress(Form("P.dc.Ch%s.stub_x",chname2[ich]),&stub_x[ich]) ;
   tsimc->SetBranchAddress(Form("P.dc.Ch%s.stub_xp",chname2[ich]),&stub_xp[ich]) ;
   tsimc->SetBranchAddress(Form("P.dc.Ch%s.stub_y",chname2[ich]),&stub_y[ich]) ;
   tsimc->SetBranchAddress(Form("P.dc.Ch%s.stub_yp",chname2[ich]),&stub_yp[ich]) ;
 }
   // Define histograms
 TH1F* hetotnorm;
 TH1F* hndchit;
 TH2F* hstub1_xp_diff;
 TH1F* hstubx_diff;
 TH1F* hstubxp_diff;
 TH1F* hstuby_diff;
 TH1F* hstubyp_diff;
 TH1F* hwire[plnum];
 TH1F* hdist[plnum];
 TH1F* htime[plnum];
 TH1F* hresid[plnum];
 TH2F* hwire_time[plnum];
 TH2F* hwire_dist[plnum];
 TH2F* hwire_resid[plnum];
 TH2F* hstub1_xp_x;
 TH2F* hstub2_xp_x;
 TH2F* hstubxpdiff_xfp;
 TH2F* hstubxdiff_xfp;
 TH2F* htime_dist[plnum];
 TString temp;
   temp=Form("Run %d ; Etotnorm  ; Counts",nrun);
   hetotnorm = new TH1F("hetotnorm",temp,100,0,2.0);
   temp=Form("Run %d ; Number DC hits  ; Counts",nrun);
   hndchit = new TH1F("hndchit",temp,20,0,20.);
   temp=Form("Run %d ; Track Stub1 Xp -K*X  ; Track Stub1 X",nrun);
   hstub1_xp_diff = new TH2F("hstub1_xp_diff",temp,100,-1.,1.,100,-50.,50.);
   temp=Form("Run %d ; Track Stub1 Xp  ; Track Stub1 X",nrun);
   hstub1_xp_x = new TH2F("hstub1_xp_x",temp,100,-2.,2.,100,-50.,50.);
   temp=Form("Run %d ; Track Stub2 Xp  ; Track Stub2 X",nrun);
   hstub2_xp_x = new TH2F("hstub2_xp_x",temp,100,-2.,2.,100,-50.,50.);
   temp=Form("Run %d ; Stub Xp diff ; Xfp",nrun);
   hstubxpdiff_xfp = new TH2F("hstubxpdiff_Xfp",temp,100,-1.2,1.2,100,-50.,50.);
   temp=Form("Run %d ; Stub X diff ; Xfp",nrun);
   hstubxdiff_xfp = new TH2F("hstubxdiff_Xfp",temp,100,-15.,15.,100,-50.,50.);
 temp=Form("Run %d ; Stub X diff ; Counts",nrun);
   hstubx_diff = new TH1F("hstubx_diff",temp,200,-15.,15.);
   temp=Form("Run %d ; Stub Y diff ; Counts",nrun);
   hstuby_diff = new TH1F("hstuby_diff",temp,200,-10.,10.);
   temp=Form("Run %d ; Stub Xp diff ; Counts",nrun);
   hstubxp_diff = new TH1F("hstubxp_diff",temp,200,-1.,1.);
   temp=Form("Run %d ; Stub Yp diff ; Counts",nrun);
   hstubyp_diff = new TH1F("hstubyp_diff",temp,200,-1.,1.);
   for (Int_t ip=0;ip<plnum;ip++) {
   temp=Form("Run %d ; Plane  %s%s Resid ; Counts",nrun,chname[ip],plname[ip]);
   hresid[ip] = new TH1F(Form("resid_%s%s",chname[ip],plname[ip]),temp,200,-.5,.5);
   temp=Form("Run %d ; Plane  %s%s Wire  ; Counts",nrun,chname[ip],plname[ip]);
   hwire[ip] = new TH1F(Form("wire_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip]);
   temp=Form("Run %d ; Plane  %s%s Dist ; Counts",nrun,chname[ip],plname[ip]);
   hdist[ip] = new TH1F(Form("dist_%s%s",chname[ip],plname[ip]),temp,120,0.,.6);
   temp=Form("Run %d ; Plane  %s%s Dist ; Counts",nrun,chname[ip],plname[ip]);
   htime[ip] = new TH1F(Form("time_%s%s",chname[ip],plname[ip]),temp,250,-50,250);
  temp=Form("Run %d ; Plane  %s%s Wire ; Resid",nrun,chname[ip],plname[ip]);
   hwire_resid[ip] = new TH2F(Form("wire_resid_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],200,-.2,.2);
   temp=Form("Run %d ; Plane  %s%s Wire ; Time",nrun,chname[ip],plname[ip]);
   hwire_time[ip] = new TH2F(Form("wire_time_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],125,-50,200);
   temp=Form("Run %d ; Plane  %s%s Wire ; Dist",nrun,chname[ip],plname[ip]);
   hwire_dist[ip] = new TH2F(Form("wire_dist_%s%s",chname[ip],plname[ip]),temp,nwires[ip],0,nwires[ip],120,0,.6);
   temp=Form("Run %d ; Plane  %s%s Time ; Dist",nrun,chname[ip],plname[ip]);
   htime_dist[ip] = new TH2F(Form("time_dist_%s%s",chname[ip],plname[ip]),temp,200,-50,250,120,0,.6);
 }
// loop over entries
 Int_t nch1;
 Int_t nch2;
 Double_t stub_xpdiff;
 Double_t stub_xdiff;
 Double_t stub_ypdiff;
 Double_t stub_ydiff;
Long64_t nentries = tsimc->GetEntries();
 if (nentries>500000) nentries=500000.;
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                nch1=0;
                nch2=0;
                stub_ydiff=-9999;
                stub_ypdiff=-9999;
                stub_xdiff=-9999;
                stub_xpdiff=-9999;
                for (Int_t ip=0;ip<plnum;ip++) {
		  if (ip<6&&nhits[ip]==1) nch1++;
		  if (ip<5&&nhits[ip]==1) nch2++;
		}
		hetotnorm->Fill(etotnorm);
                 for (Int_t ip=0;ip<plnum;ip++) {
                  hresid[ip]->Fill(resid[ip]);
		   for (Int_t id=0;id<nhits[ip];id++) {
		   if (id==0) hwire_resid[ip]->Fill(wire[ip][id],resid[ip]);
                   hwire[ip]->Fill(wire[ip][id]);
                   hdist[ip]->Fill(dist[ip][id]);
                   htime[ip]->Fill(time[ip][id]);
                   hwire_time[ip]->Fill(wire[ip][id],time[ip][id]);
                   hwire_dist[ip]->Fill(wire[ip][id],dist[ip][id]);
                   htime_dist[ip]->Fill(time[ip][id],dist[ip][id]);
		   }
	          }
		if (nch1>4&&nch2>4&&etotnorm>.6&& ntrack>0) {
		  if (sp1_id <= nsp[0] && sp2_id <= nsp[1]) {
                     stub_ydiff=stub_y[0][int(sp1_id)]-stub_y[1][int(sp2_id)];
                     stub_ypdiff=stub_yp[0][int(sp1_id)]-stub_yp[1][int(sp2_id)];
                     stub_xdiff=stub_x[0][int(sp1_id)]-stub_x[1][int(sp2_id)];
                      stub_xpdiff=stub_xp[0][int(sp1_id)]-stub_xp[1][int(sp2_id)];
		    hstubx_diff->Fill(stub_xdiff);
		    hstuby_diff->Fill(stub_ydiff);
		    hstubxp_diff->Fill(stub_xpdiff);
		    hstubyp_diff->Fill(stub_ypdiff);
		  }
                  hndchit->Fill(ndchit);
                  hstub1_xp_diff->Fill(stub_xp[0][int(sp1_id)]-.0018*stub_x[0][int(sp1_id)],stub_x[0][int(sp1_id)]);
                 hstub1_xp_x->Fill(stub_xp[0][int(sp1_id)],stub_x[0][int(sp1_id)]);
		 if (ntrack ==1) hstubxpdiff_xfp->Fill(stub_xpdiff,xfp);
                  hstub2_xp_x->Fill(stub_xp[1][int(sp2_id)],stub_x[1][int(sp2_id)]);
		  hstubxdiff_xfp->Fill(stub_xdiff,xfp);
 		} // if pid /event selection
	}
//
}
