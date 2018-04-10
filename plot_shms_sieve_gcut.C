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
#include <TCutG.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void plot_shms_sieve(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(1);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="../../ROOTfiles/"+basename+".root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
// Define branches
 Double_t  sumnpe;
   tsimc->SetBranchAddress("P.ngcer.npeSum",&sumnpe);
 Double_t  etracknorm;
   tsimc->SetBranchAddress("P.cal.etracknorm",&etracknorm);
 Double_t  W;
   tsimc->SetBranchAddress("P.kin.W",&W);
 Double_t  ytar;
   tsimc->SetBranchAddress("P.gtr.y",&ytar);
 Double_t  delta;
   tsimc->SetBranchAddress("P.gtr.dp",&delta);
 Double_t  yptar;
   tsimc->SetBranchAddress("P.gtr.ph",&yptar);
 Double_t  xptar;
   tsimc->SetBranchAddress("P.gtr.th",&xptar);
 Double_t  yfp;
   tsimc->SetBranchAddress("P.dc.y_fp",&yfp);
 Double_t  ypfp;
   tsimc->SetBranchAddress("P.dc.yp_fp",&ypfp);
 Double_t  xfp;
   tsimc->SetBranchAddress("P.dc.x_fp",&xfp);
 Double_t  xpfp;
   tsimc->SetBranchAddress("P.dc.xp_fp",&xpfp);
   // Define histograms
   TH2F *hxfp_yfp_foil[10];
   TH1F *hys_foil[10];
   Double_t xptar_cut_lo[10]={-0.05,-0.04,-0.032,-0.022,-0.015,-0.005,0.006,0.016,0.026,0.035};
   Double_t xptar_cut_hi[10]={-0.04,-0.032,-0.022,-0.015,-0.005,0.006,0.016,0.026,0.035,0.045};
   for (Int_t i=0;i<10;i++) {
    hxfp_yfp_foil[i] = new TH2F(Form("hxfp_yfp_%d",i+1), Form("Run %d Xptar cut %5.3f to %5.3f ; X_fp ; Y_fp",nrun,xptar_cut_lo[i],xptar_cut_hi[i]), 200,-40,40,200,-40,40);
    hys_foil[i] = new TH1F(Form("hys_%d",i+1), Form("Run %d Cut %d; Y_sieve ; Counts",nrun,i+1), 200, -10,10);
   }
   TH2F *hxfp_yfp = new TH2F("hxfp_yfp", Form("Run %d ; X_fp ; Y_fp",nrun), 200,-40,40,200,-40,40);
   TH2F *hxs_ys = new TH2F("hxs_ys", Form("Run %d ; Y_s ; X_s",nrun), 200,-10,10,200,-20,20);
   TH1F *hxptar = new TH1F("hxptar", Form("Run %d ; Xp_tar ; Counts",nrun), 200, -.05,.05);
   TH1F *hyptar = new TH1F("hyptar", Form("Run %d ; Yp_tar ; Counts",nrun), 200, -.05,.05);
   Double_t xs;
   Double_t ys;
// loop over entries
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe>10&&etracknorm > .6 ) {
		  xs=xptar*253.;
		  ys=(-0.019*(delta+0.0))+yptar*253.-40.*0.00052*(delta+0.0);
                  hxs_ys->Fill(ys,xs);
 		}
	}
	//
	/* Set cuts of xs versus ys
    TCutG* cutg;
    TCutG *tmpg,*mycutg;
	TCutG *t ;  
    TString outputhist;
   outputhist="hist/"+basename+"_cut.root";
   TFile f(outputhist,"UPDATE");
TCanvas *c3 = new TCanvas("c3", "Sieve", 900,800);
c3->Divide(1,1);
c3->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
 hxs_ys->Draw("colz");
   Int_t n=0;
        t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",n));
	while(t) {
                 c3->cd(1);
		 t->Draw("same");
		 t->SetLineColor(2);
		 n++;
                 t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",n));
                 c3->Update();
		 cout << n << endl;
		 }
     Int_t nc=0;
       while (nc!=-1) {
        t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",nc));
                 c3->cd(1);
	if (t) t->Draw("same");
	if (t) {
             cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	if (nc!=-1) {
    cutg=(TCutG*)gPad->WaitPrimitive("CUTG","CutG");
    c3->Update();
    //cout << cutg->GetN() << endl;
    tmpg= (TCutG*)gROOT->GetListOfSpecials()->FindObject("CUTG");
    //cout << tmpg->GetN() << endl;
    mycutg=(TCutG*)(tmpg->Clone(Form("xs_ys_cut_%d",nc)));
    //cout << mycutg->GetN() <<endl;
    nc++;
    mycutg->Write("",TObject::kOverwrite);
    mycutg->Print();
    mycutg->Draw();
       cout <<"enter nc, present nc = " << nc << endl;
             cin >> nc ;
	}
	}
	//	 while ((TCutG*)f.Get(Form("mycutg_%d",n))) {
	gDirectory->ls("m");
        t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",n));
	while(t) {
		 t->Draw();
		 t->SetLineColor(2);
		 n++;
        t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",n));
        c3->Update();
		 }
	//
    outputpdf="plots/"+basename+"_shms_sieve_cuts.pdf";
c3->Print(outputpdf);
	*/
//
	// loop data and apply cuts
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
		if (sumnpe>10&&etracknorm > .6 ) {
		  hxfp_yfp->Fill(xfp,yfp);
		  hxptar->Fill(xptar);
		  hyptar->Fill(yptar);
		  xs=xptar*253.;
		  ys=(-0.019*(delta+0.0))+yptar*253.-40.*0.00052*(delta+0.0);
                  hxs_ys->Fill(ys,xs);
                  for (Int_t i=0;i<10;i++) {
		    t=(TCutG*)gROOT->FindObject(Form("xs_ys_cut_%d",n));
		    if (t->InsidePolygon(xs,ys)) ys_foil[i]->Fill(ys);
                    if (xptar>xptar_cut_lo[i]&&xptar<=xptar_cut_hi[i]) {
		      hxfp_yfp_foil[i]->Fill(xfp,yfp);
		    }
                  }
		}
	}
	//
 gStyle->SetOptStat(0);
	// plot data
TCanvas *c = new TCanvas("c", "focal plane", 900,800);
c->Divide(1,1);
c->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
hxfp_yfp->Draw("colz");
    outputpdf="plots/"+basename+"_shms_sieve.pdf(";
c->Print(outputpdf);
//
//

 TCanvas *cfoil[10];
 for (Int_t i=0;i<10;i++) {
  cfoil[i] = new TCanvas(Form("c_%d",i), Form("Focal plane xptar %d",i), 900,800);
cfoil[i]->Divide(1,1);
cfoil[i]->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
hxfp_yfp_foil[i]->Draw("colz");
    outputpdf="plots/"+basename+"_shms_sieve.pdf";
cfoil[i]->Print(outputpdf);
}
//
TCanvas *c1 = new TCanvas("c1", "angles tar", 900,800);
c1->Divide(1,2);
c1->cd(1);
 gPad->SetGridy();
 gPad->SetGridx();
hxptar->Draw();
c1->cd(2);
 gPad->SetGridy();
 gPad->SetGridx();
hyptar->Draw();
    outputpdf="plots/"+basename+"_shms_sieve.pdf)";
c1->Print(outputpdf);

 //
}
