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

void make_hist_shms_hodo(TString basename="",Int_t nrun=2043){
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
   outputhist= "hist/"+basename+"_shms_hodo_hist.root";
 TObjArray HList(0);
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t plnum=4;
 static const Int_t clmax=10;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 Double_t goodscinhit;
 tsimc->SetBranchAddress("P.hod.goodscinhit",&goodscinhit) ;
 Int_t ncl[plnum];
 Double_t pos[plnum][clmax];
 Double_t size[plnum][clmax];
 Double_t flag[plnum][clmax];
 Double_t usedflag[plnum][clmax];
 for (Int_t ip=0;ip<plnum;ip++) {
   tsimc->SetBranchAddress(Form("Ndata.P.hod.%s.Clus.Pos",plname[ip]),&ncl[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Pos",plname[ip]),&pos[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Size",plname[ip]),&size[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.Flag",plname[ip]),&flag[ip]) ;
   tsimc->SetBranchAddress(Form("P.hod.%s.Clus.UsedFlag",plname[ip]),&usedflag[ip]) ;
 }
 //
 TH1F*  hgoodscinhit = new TH1F("hgoodscinhit",";Good Scin Hit ;",10,0,10);
 HList.Add(hgoodscinhit);
 TH1F*  hxdiff_nclusone = new TH1F("hxdiff_nclusone","One cluster in plane; X1-X2 (cm)  ;",80,-20,20);
 TH1F*  hxdiff_nclusone_2 = new TH1F("hxdiff_nclusone_2","No shift One cluster in plane; X1-X2 (cm)  ;",80,-20,20);
 TH1F*  hydiff_nclusone = new TH1F("hydiff_nclusone","One cluster in plane; Y1-Y2 (cm)  ;",80,-20,20);
 TH2F*  hxdiff_xpos_nclusone = new TH2F("hxdiff_xpos_nclusone","One cluster in plane; X1-X2 (cm)  ; X1",120,-30,30,120,-30,30);
 TH2F*  hydiff_ypos_nclusone = new TH2F("hydiff_ypos_nclusone","One cluster in plane; Y1-Y2 (cm)  ; Y1",100,-40,40,100,-40,40);
 TH1F*  hclus[plnum];
 TH1F*  hpos[plnum];
 for (Int_t i=0;i<plnum;i++) {
   hclus[i] = new TH1F(Form("Hclus_%d",i),Form("; Pl %d Number of Cluster",i),10,0,10);
   hpos[i] = new TH1F(Form("Hpos_%d",i),Form("; Pl %d Position",i),100,-60,60);
 }
 //
 Double_t xpfp_to_xfp=0.0018;//*0.66;
Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (i%10000==0) cout << " Entry = " << i << endl;
		hgoodscinhit->Fill(goodscinhit);
		if (goodscinhit==1) {
                  for (Int_t i=0;i<plnum;i++) {
                    hclus[i]->Fill(ncl[i]);
		    if (ncl[i]==1) hpos[i]->Fill(pos[i][0]);
                   }
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_nclusone->Fill(pos[0][0]*(1+xpfp_to_xfp*200)-pos[2][0]);
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_nclusone_2->Fill(pos[0][0]-pos[2][0]);
		  if (ncl[1]==1&&ncl[3]==1) hydiff_nclusone->Fill(pos[1][0]-pos[3][0]);
		  if (ncl[0]==1&&ncl[2]==1) hxdiff_xpos_nclusone->Fill(pos[0][0]*(1+xpfp_to_xfp*200)-pos[2][0],pos[0][0]);
		  if (ncl[1]==1&&ncl[3]==1) hydiff_ypos_nclusone->Fill(pos[1][0]-pos[3][0],pos[1][0]);
		}
	}
//
}
