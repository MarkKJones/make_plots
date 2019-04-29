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
#include <TRandom3.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
using namespace std;

void make_hist_hms_hodo(TString basename="",Int_t nrun=1267,Double_t cpeak=45){
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
   outputhist= "hist/"+basename+"_hms_hodo_hist.root";
 TObjArray HList(0);
 //

//
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("T");
//
 static const Int_t plnum=4;
 static const Int_t iside=2;
 const char* plname[plnum]={"1x","1y","2x","2y"};
 const char* sidename[iside]={"Neg","Pos"};
 const char* sname[iside]={"neg","pos"};
 static const Int_t npad[plnum]={16,10,16,10};
 static const Int_t npadshms[plnum]={13,13,14,21};
 //
 Double_t tdc_uncorr[plnum][iside][16];
 Double_t hmsadc_time[plnum][iside][16];
 Double_t hmstdc_tofcorr[plnum][iside][16];
 Int_t hmshodo_tdccounter_nhits[plnum][iside];
 Double_t hmshodo_tdccounter[plnum][iside][100];
 Double_t hmshodo_nhits[plnum];
 Double_t hmshodo_Xpos[plnum];
  Double_t hmshodo_Ypos[plnum];
Double_t shmstdc_uncorr[plnum][iside][21];
Double_t shmstdc_tofcorr[plnum][iside][21];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
   tsimc->SetBranchAddress(Form("H.hod.%s.nhits",plname[ipl]),&hmshodo_nhits[ipl]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.TrackXPos",plname[ipl]),&hmshodo_Xpos[ipl]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.TrackYPos",plname[ipl]),&hmshodo_Ypos[ipl]) ;
 for (Int_t is=0;is<iside;is++) {
   tsimc->SetBranchAddress(Form("Ndata.H.hod.%s.%sTdcCounter",plname[ipl],sname[is]),&hmshodo_tdccounter_nhits[ipl][is]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.%sTdcCounter",plname[ipl],sname[is]),&hmshodo_tdccounter[ipl][is]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.Good%sTdcTimeUnCorr",plname[ipl],sidename[is]),&tdc_uncorr[ipl][is]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.Good%sAdcPulseTime",plname[ipl],sidename[is]),&hmsadc_time[ipl][is]) ;
   tsimc->SetBranchAddress(Form("H.hod.%s.Good%sTdcTimeTOFCorr",plname[ipl],sidename[is]),&hmstdc_tofcorr[ipl][is]) ;
 }}
  Double_t  hcer;
   tsimc->SetBranchAddress("H.cer.npeSum",&hcer);
 Double_t  hntr;
   tsimc->SetBranchAddress("H.dc.ntrack",&hntr);
 Double_t  hstarttime;
   tsimc->SetBranchAddress("H.hod.starttime",&hstarttime);
 //
    //
 TH1F *htdc_uncorr[plnum][iside][16];
 TH2F *htdc_uncorr_pm[plnum][16];
 for (Int_t ipl=0;ipl<plnum;ipl++) {
 for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
   htdc_uncorr_pm[ipl][ipad]= new TH2F(Form("tdc_uncorr_pm_%s_pad_%d",plname[ipl],ipad),Form(" %s pad_%d; HMS TDC POS ; HMS TDC NEG",plname[ipl],ipad),240,20,80,240,20,80.);
   HList.Add(htdc_uncorr_pm[ipl][ipad]);
 for (Int_t is=0;is<iside;is++) {
   htdc_uncorr[ipl][is][ipad]= new TH1F(Form("tdc_uncorr_%s_%s_pad_%d",plname[ipl],sidename[is],ipad),Form(" %s %s pad_%d; HMS TDc uncorr;counts",plname[ipl],sidename[is],ipad),240,20,80.);
   HList.Add(htdc_uncorr[ipl][is][ipad]);
 }}}
 // 
    Long64_t nentries = tsimc->GetEntries();
    //
    //  nentries=50000;
	for (int ie = 0; ie < nentries; ie++) {
      		tsimc->GetEntry(ie);
                if (ie%50000==0) cout << " Entry = " << ie << endl;
		Bool_t selcut=hntr >0 && hcer>2. ;
                for (Int_t ipl=0;ipl<plnum;ipl++) {
                  for (Int_t ipad=0;ipad<npad[ipl];ipad++) {
  		     if (tdc_uncorr[ipl][0][ipad]<1000&&tdc_uncorr[ipl][1][ipad]<1000 ) htdc_uncorr_pm[ipl][ipad]->Fill(tdc_uncorr[ipl][1][ipad],tdc_uncorr[ipl][0][ipad]);		    
                    for (Int_t is=0;is<iside;is++) {
		    htdc_uncorr[ipl][is][ipad]->Fill(tdc_uncorr[ipl][is][ipad]);
		   }
		}
		}
	}

 TFile hsimc(outputhist,"recreate");
	HList.Write();

}
