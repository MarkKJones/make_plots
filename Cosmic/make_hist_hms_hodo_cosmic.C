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

void make_hist_hms_hodo_cosmic(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
   static const char* spec="H";
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1111);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= "hist/"+basename+"_hms_hodo_cosmic_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 static const Int_t plnum=4;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 static const Int_t npad[plnum]={16,10,16,10};
 Double_t plnhits[plnum];
 Double_t negtime[plnum][16];
 Double_t postime[plnum][16];
 Double_t negtime_uncorr[plnum][16];
 Double_t postime_uncorr[plnum][16];
 Double_t negtime_corr[plnum][16];
 Double_t postime_corr[plnum][16];
 Double_t negamp[plnum][16];
 Double_t posamp[plnum][16];
 Int_t plcross[plnum]={1,0,3,2};
 Double_t space[plnum]={-7.5,7.5,-7.5,7.5};
 Double_t pos_pad1[plnum]={33.75,-56.25,33.75,-56.25};
 Double_t  fptime[plnum];
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("%s.hod.%s.nhits",spec,plname[ip]),&plnhits[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodNegTdcTimeWalkCorr",spec,plname[ip]),&negtime[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodPosTdcTimeWalkCorr",spec,plname[ip]),&postime[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodNegTdcTimeUnCorr",spec,plname[ip]),&negtime_uncorr[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodPosTdcTimeUnCorr",spec,plname[ip]),&postime_uncorr[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodNegTdcTimeCorr",spec,plname[ip]),&negtime_corr[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodPosTdcTimeCorr",spec,plname[ip]),&postime_corr[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodNegAdcPulseAmp",spec,plname[ip]),&negamp[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.GoodPosAdcPulseAmp",spec,plname[ip]),&posamp[ip]) ;
   tsimc->SetBranchAddress(Form("%s.hod.%s.fptime",spec,plname[ip]),&fptime[ip]) ;
 }
 Double_t  ntr;
 tsimc->SetBranchAddress(Form("%s.dc.ntrack",spec),&ntr);
 Double_t  xfp;
   tsimc->SetBranchAddress(Form("%s.dc.x_fp",spec),&xfp);
 Double_t  yfp;
   tsimc->SetBranchAddress(Form("%s.dc.y_fp",spec),&yfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress(Form("%s.dc.xp_fp",spec),&xpfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress(Form("%s.dc.yp_fp",spec),&ypfp);
  Double_t  betanotrack;
   tsimc->SetBranchAddress(Form("%s.hod.betanotrack",spec),&betanotrack);
 Double_t  betatrack;
   tsimc->SetBranchAddress(Form("%s.hod.beta",spec),&betatrack);
 Double_t  delta;
   tsimc->SetBranchAddress(Form("%s.gtr.dp",spec),&delta);
 Double_t  starttime;
   tsimc->SetBranchAddress(Form("%s.hod.starttime",spec),&starttime);
   // Define histograms
  TH1F *hxfp = new TH1F("hxfp",Form("Run %d ; Xfp;Counts",nrun), 160,-40,40);
   HList.Add(hxfp);
   TH1F *hyfp = new TH1F("hyfp",Form("Run %d ; Yfp;Counts",nrun), 160,-40,40);
   HList.Add(hyfp);
   TH1F *hxpfp = new TH1F("hxpfp",Form("Run %d ; Xpfp;Counts",nrun), 160,-.1,.1);
   HList.Add(hxpfp);
   TH1F *hypfp = new TH1F("hypfp",Form("Run %d ; Ypfp;Counts",nrun), 160,-.1,.1);
   HList.Add(hypfp);
    TH1F *hbetanotrack = new TH1F("hbetanotrack",Form("Run %d ; Beta notrack;Counts",nrun), 300, -1.5,1.5);
    HList.Add(hbetanotrack);
    TH2F *hbetanotrack_xfp = new TH2F("hbetanotrack_xfp",Form("Run %d ; Beta notrack;Xfp",nrun), 200, -1.5,1.5,100,-60.,60.);
    HList.Add(hbetanotrack_xfp);
    TH2F *hbetanotrack_yfp = new TH2F("hbetanotrack_yfp",Form("Run %d ; Beta notrack;Yfp",nrun), 200, -1.5,1.5,100,-60.,60.);
    HList.Add(hbetanotrack_yfp);
    TH1F *hbetatrack = new TH1F("hbetatrack",Form("Run %d ; Beta track;Counts",nrun), 300, -1.5,1.5);
    HList.Add(hbetatrack);
    TH1F *hstarttime = new TH1F("hstarttime",Form("Run %d ; Starttime;Counts",nrun), 280, -10.,60.0);
    HList.Add(hstarttime);
 //
 TH1F* hfptime[plnum];
 TH1F* neghist[plnum][16];
 TH1F* poshist[plnum][16];
 TH1F* negadchist[plnum][16];
 TH1F* posadchist[plnum][16];
 TH1F* diffhist[plnum][16];
 TH1F* diffadchist[plnum][16];
 TH1F* diffhistcent[plnum][16];
 TH1F* difftdcuncorrcent[plnum][16];
 TH1F* difftdccorrcent[plnum][16];
 TH1F* diffadchistcent[plnum][16];
 TH1F* diffposhistcent[plnum][16];
 TH1F* diffneghistcent[plnum][16];
 TH1F* hdiff_hit2_posdiff[plnum][16];
 TH1F* hdiff_hit2_negdiff[plnum][16];
 for (Int_t ip=0;ip<plnum;ip++) {
   hfptime[ip]= new TH1F(Form("hfptime_%s",plname[ip]),Form("Plane %s ; FPTime ns; Counts",plname[ip]),200,0,100);
   HList.Add(hfptime[ip]);
 for (Int_t ipd=0;ipd<npad[ip];ipd++) {
   neghist[ip][ipd] = new TH1F(Form("hist_%s_neg_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d ;Neg Time ns; COunts",plname[ip],ipd+1),200,-100,100);
   poshist[ip][ipd] = new TH1F(Form("hist_%s_pos_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d ;Pos Time ns; COunts",plname[ip],ipd+1),200,-100,100);
   negadchist[ip][ipd] = new TH1F(Form("hist_%s_adcneg_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d ;Neg ADc amp; COunts",plname[ip],ipd+1),75,0,300);
   posadchist[ip][ipd] = new TH1F(Form("hist_%s_adcpos_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d ;Pos ADc amp; COunts",plname[ip],ipd+1),75,0,300);
   diffhist[ip][ipd] = new TH1F(Form("hist_%s_diff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Pos-Neg Time ns; COunts",plname[ip],ipd+1),200,-100,100);
   diffadchist[ip][ipd] = new TH1F(Form("hist_%s_adcdiff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Pos-Neg Adc amp; COunts",plname[ip],ipd+1),500,-100,100);
   diffhistcent[ip][ipd] = new TH1F(Form("histcent_%s_diff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; TW Corr Pos-Neg Time ns; COunts",plname[ip],ipd+1),60,-100,100);
   difftdcuncorrcent[ip][ipd] = new TH1F(Form("histcent_%s_tdcuncorrdiff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; Uncorr Pos-Neg Time ns; COunts",plname[ip],ipd+1),60,-100,100);
   difftdccorrcent[ip][ipd] = new TH1F(Form("histcent_%s_tdccorrdiff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; Corr Pos-Neg Time ns; COunts",plname[ip],ipd+1),60,-100,100);
   diffadchistcent[ip][ipd] = new TH1F(Form("histcent_%s_adcdiff_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; Pos-Neg Amp; COunts",plname[ip],ipd+1),25,-50,50);
   diffposhistcent[ip][ipd] = new TH1F(Form("histcent_%s_diffpos_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; Pos- Pos(cross pad) ns; COunts",plname[ip],ipd+1),100,-100,100);
   diffneghistcent[ip][ipd] = new TH1F(Form("histcent_%s_diffneg_pad_%d",plname[ip],ipd+1),Form("Plane %s Pad %d Center cross pad ; Neg- Pos(cross pad) ns; COunts",plname[ip],ipd+1),100,-100,100);
   hdiff_hit2_posdiff[ip][ipd] = new TH1F(Form("hdiff_%s_diffpos_2hit_pad_%d",plname[ip],ipd+1),Form("Plane %s  ; Time diff pos pad %d - pad %d ns; COunts",plname[ip],ipd+1,ipd+2),60,-100,100);
   hdiff_hit2_negdiff[ip][ipd] = new TH1F(Form("hdiff_%s_diffneg_2hit_pad_%d",plname[ip],ipd+1),Form("Plane %s  ; Time diff neg pad %d - pad %d ns; COunts",plname[ip],ipd+1,ipd+2),60,-100,100);
   HList.Add(diffhistcent[ip][ipd]);
    HList.Add(difftdcuncorrcent[ip][ipd]);
  HList.Add(difftdccorrcent[ip][ipd]);
   HList.Add(negadchist[ip][ipd]);
   HList.Add(posadchist[ip][ipd]);
   HList.Add(diffadchistcent[ip][ipd]);
   HList.Add(diffposhistcent[ip][ipd]);
   HList.Add(diffneghistcent[ip][ipd]);
   HList.Add(hdiff_hit2_posdiff[ip][ipd]);
   HList.Add(hdiff_hit2_negdiff[ip][ipd]);
 } }
 //
 //
 TH2D *hdtime_dist[plnum][16];
 for (Int_t ip=0;ip<plnum;ip++) {
 for (Int_t ipd=0;ipd<npad[ip];ipd++) {
   hdtime_dist[ip][ipd] = new  TH2D(Form("hist_%s_pad_%d",plname[ip],ipd+1),Form("Plane %s  Pad %d; Dist (cm)  ; Diff time ns",plname[ip],ipd+1),80,-40,40,100,-20,30);
   HList.Add(hdtime_dist[ip][ipd]);
 }
 }
 //
// loop over entries
 Int_t pl_padnum[plnum];
 Double_t pl_padpostime[plnum];
 Double_t pl_padnegtime[plnum];
 Double_t dcut[4]={7.,7.,7.,7.};
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
if (plnhits[0]==1&&plnhits[1]==1&&plnhits[2]==1&&plnhits[3]==1) {
  //  if (TMath::Abs(xfp) < 40) {
    hxfp->Fill(xfp);
    hyfp->Fill(yfp);
    hxpfp->Fill(xpfp);
    hypfp->Fill(ypfp);
		hbetanotrack->Fill(betanotrack);
		hbetatrack->Fill(betatrack);
		hstarttime->Fill(starttime);
                hbetanotrack_xfp->Fill(betanotrack,xfp);
                hbetanotrack_yfp->Fill(betanotrack,yfp);
		//}		  
 for (Int_t ip=0;ip<plnum;ip++) {
   hfptime[ip]->Fill(fptime[ip]);
 for (Int_t ipd=0;ipd<npad[ip];ipd++) {
   if (postime[ip][ipd] < 3000 && negtime[ip][ipd] < 3000) {
       pl_padnum[ip]=ipd;
       pl_padpostime[ip]=postime[ip][ipd];
       pl_padnegtime[ip]=negtime[ip][ipd];
   }
 }}
  for (Int_t ip=0;ip<plnum;ip++) {
   for (Int_t ipd=0;ipd<npad[ip];ipd++) {
    if (negtime[ip][ipd] < 3000) neghist[ip][ipd]->Fill(negtime[ip][ipd]);
    if (postime[ip][ipd] < 3000) poshist[ip][ipd]->Fill(postime[ip][ipd]);
    if (postime[ip][ipd] < 3000 && negtime[ip][ipd] < 3000) {
        diffhist[ip][ipd]->Fill(postime[ip][ipd]-negtime[ip][ipd]);
        Double_t dist=pos_pad1[ip]+space[ip]*pl_padnum[plcross[ip]];
        hdtime_dist[ip][ipd]->Fill(dist,postime[ip][ipd]-negtime[ip][ipd]);
		  if (TMath::Abs(dist) < dcut[ip]) {
           if (negamp[ip][ipd] >0) negadchist[ip][ipd]->Fill(negamp[ip][ipd]);
           if (posamp[ip][ipd] >0) posadchist[ip][ipd]->Fill(posamp[ip][ipd]);
	    if (posamp[ip][ipd] >0 &&negamp[ip][ipd] >0 )diffadchist[ip][ipd]->Fill(posamp[ip][ipd]-negamp[ip][ipd]);
	    diffhistcent[ip][ipd]->Fill(postime[ip][ipd]-negtime[ip][ipd]);
	    difftdcuncorrcent[ip][ipd]->Fill(postime_uncorr[ip][ipd]-negtime_uncorr[ip][ipd]);
	    difftdccorrcent[ip][ipd]->Fill(postime_corr[ip][ipd]-negtime_corr[ip][ipd]);
	    diffadchistcent[ip][ipd]->Fill(posamp[ip][ipd]-negamp[ip][ipd]);
	    diffposhistcent[ip][ipd]->Fill(postime[ip][ipd]-pl_padpostime[plcross[ip]]);
	    diffneghistcent[ip][ipd]->Fill(negtime[ip][ipd]-pl_padpostime[plcross[ip]]);
	     }
    }
   }
  } 
 }	
//
 Int_t pl_hit2_padnum[4][2];
 Double_t pl_hit2_padpostime[4][2];
 Double_t pl_hit2_padnegtime[4][2];
 for (Int_t ip=0;ip<plnum;ip++) {
 if (plnhits[ip]==2) {
   //   cout << " plane = " << ip << " nits " << plnhits[ip] << endl;
   Int_t cnt=0;
   for (Int_t ipd=0;ipd<npad[ip];ipd++) {
       if (postime[ip][ipd] < 3000 && negtime[ip][ipd] < 3000 && cnt <2) {
       pl_hit2_padnum[ip][cnt]=ipd;
       pl_hit2_padpostime[ip][cnt]=postime[ip][ipd];
       pl_hit2_padnegtime[ip][cnt++]=negtime[ip][ipd];
       }
   }
   // cout << " cnt = " << cnt << " " << pl_hit2_padnum[ip][0] << " " << pl_hit2_padnum[ip][1] << endl;
   if ( TMath::Abs(pl_hit2_padnum[ip][0]-pl_hit2_padnum[ip][1])==1 && cnt==2) {
     Int_t padlo = pl_hit2_padnum[ip][0];
     if (pl_hit2_padnum[ip][1]<padlo) padlo = pl_hit2_padnum[ip][1];
      hdiff_hit2_posdiff[ip][padlo]->Fill(pl_hit2_padpostime[ip][0]-pl_hit2_padpostime[ip][1]);
      hdiff_hit2_negdiff[ip][padlo]->Fill(pl_hit2_padnegtime[ip][0]-pl_hit2_padnegtime[ip][1]);
    }
 }
 }
//
 } //nentries
	//
 TFile hsimc(outputhist,"recreate");
	HList.Write();
//
}
