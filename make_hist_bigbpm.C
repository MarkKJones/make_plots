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
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include<math.h>
using namespace std;

void make_hist_bigbpm(TString basename="",Int_t nrun=2043){
   if (basename=="") {
     cout << " Input the basename of the root file (assumed to be in worksim)" << endl;
     cin >> basename;
   }
gStyle->SetPalette(1,0);
 gStyle->SetOptStat(11111);
 gStyle->SetOptFit(11);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.06,"XY");
 gStyle->SetPadLeftMargin(0.12);
   TString inputroot;
   inputroot="ROOTfiles/"+basename+".root";
   TString outputhist;
   outputhist= basename+"_hist.root";
 TObjArray HList(0);
     TString outputpdf;
    outputpdf=basename+".pdf";
 TString htitle=basename;
 TPaveLabel *title = new TPaveLabel(.15,.90,0.95,.99,htitle,"ndc");
 //
TFile *fsimc = new TFile(inputroot); 
TTree *tsimc = (TTree*) fsimc->Get("E");
// Define branches
 Double_t  ibcm1;
   tsimc->SetBranchAddress("ibcm1",&ibcm1);
 Double_t  HBcur;
   tsimc->SetBranchAddress("ecSHB_I_coarse",&HBcur);
 Double_t  Q1cur;
   tsimc->SetBranchAddress("ecSQ1_I_coarse",&Q1cur);
 Double_t  Q2cur;
   tsimc->SetBranchAddress("ecSQ2_I_coarse",&Q2cur);
 Double_t  Q3cur;
   tsimc->SetBranchAddress("ecSQ3_I_coarse",&Q3cur);
 Double_t  DIcur;
   tsimc->SetBranchAddress("ecSDI_I_coarse",&DIcur);
 Double_t  SAXpos;
   tsimc->SetBranchAddress("IPM3H07A.XPOS",&SAXpos);
 Double_t  SAXraw;
   tsimc->SetBranchAddress("IPM3H07A.XRAW",&SAXraw);
 Double_t  SAYpos;
   tsimc->SetBranchAddress("IPM3H07A.YPOS",&SAYpos);
 Double_t  bpm8Xpos;
   tsimc->SetBranchAddress("IPM3H08.XPOS",&bpm8Xpos);
 Double_t  bpm8Ypos;
   tsimc->SetBranchAddress("IPM3H08.YPOS",&bpm8Ypos);
 Double_t  bpm9Xpos;
   tsimc->SetBranchAddress("IPM3H09.XPOS",&bpm9Xpos);
 Double_t  bpm9Ypos;
   tsimc->SetBranchAddress("IPM3H09.YPOS",&bpm9Ypos);
 Double_t  bpm9Xrot;
   tsimc->SetBranchAddress("IPM3H09.XROT",&bpm9Xrot);
 Double_t  bpm9Yrot;
   tsimc->SetBranchAddress("IPM3H09.YROT",&bpm9Yrot);
 Double_t  bpm9Xp;
   tsimc->SetBranchAddress("IPM3H09.FSUM",&bpm9Xp);
 Double_t  bpm9Xm;
   tsimc->SetBranchAddress("IPM3H09.TSUM",&bpm9Xm);
 Double_t  bpm9Yp;
   tsimc->SetBranchAddress("IPM3H09.SSUM",&bpm9Yp);
 Double_t  bpm9Ym;
   tsimc->SetBranchAddress("IPM3H09.RSUM",&bpm9Ym);
 Double_t  bpm8Xrot;
   tsimc->SetBranchAddress("IPM3H08.XROT",&bpm8Xrot);
 Double_t  bpm8Yrot;
   tsimc->SetBranchAddress("IPM3H08.YROT",&bpm8Yrot);
 Double_t  bpm8Xp;
   tsimc->SetBranchAddress("IPM3H08.FSUM",&bpm8Xp);
 Double_t  bpm8Xm;
   tsimc->SetBranchAddress("IPM3H08.TSUM",&bpm8Xm);
 Double_t  bpm8Yp;
   tsimc->SetBranchAddress("IPM3H08.SSUM",&bpm8Yp);
 Double_t  bpm8Ym;
   tsimc->SetBranchAddress("IPM3H08.RSUM",&bpm8Ym);
 Double_t  SAYraw;
   tsimc->SetBranchAddress("IPM3H07A.YRAW",&SAYraw);
 Double_t  SBXpos;
   tsimc->SetBranchAddress("IPM3H07B.XPOS",&SBXpos);
 Double_t  SBXraw;
   tsimc->SetBranchAddress("IPM3H07B.XRAW",&SBXraw);
 Double_t  SBYpos;
   tsimc->SetBranchAddress("IPM3H07B.YPOS",&SBYpos);
 Double_t  SBYraw;
   tsimc->SetBranchAddress("IPM3H07B.YRAW",&SBYraw);
 Double_t  SCXpos;
   tsimc->SetBranchAddress("IPM3H07C.XPOS",&SCXpos);
 Double_t  SCXraw;
   tsimc->SetBranchAddress("IPM3H07C.XRAW",&SCXraw);
 Double_t  SCYpos;
   tsimc->SetBranchAddress("IPM3H07C.YPOS",&SCYpos);
 Double_t  SCYraw;
   tsimc->SetBranchAddress("IPM3H07C.YRAW",&SCYraw);
   Double_t Xtar_AB;
   Double_t Ytar_AB;
   Double_t Xtar_AC;
   Double_t Ytar_AC;
   Double_t Xtar_BC;
   Double_t Ytar_BC;
  // Define histograms
   TH1F *hSAXraw = new TH1F("hSAXraw",Form("Run %d ; SAX Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSAYraw = new TH1F("hSAYraw",Form("Run %d ; SAY Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSBXraw = new TH1F("hSBXraw",Form("Run %d ; SBX Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSBYraw = new TH1F("hSBYraw",Form("Run %d ; SBY Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSCXraw = new TH1F("hSCXraw",Form("Run %d ; SCX Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSCYraw = new TH1F("hSCYraw",Form("Run %d ; SCY Raw (mm); Counts",nrun),100,-5,5);
   TH1F *hSAXpos = new TH1F("hSAXpos",Form("Run %d ; SAX Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hSAYpos = new TH1F("hSAYpos",Form("Run %d ; SAY Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hSBXpos = new TH1F("hSBXpos",Form("Run %d ; SBX Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hSBYpos = new TH1F("hSBYpos",Form("Run %d ; SBY Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hSCXpos = new TH1F("hSCXpos",Form("Run %d ; SCX Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hSCYpos = new TH1F("hSCYpos",Form("Run %d ; SCY Pos (mm); Counts",nrun),100,-5,5);
   TH1F *hXtar_AB = new TH1F("hXtar_AB",Form("Run %d ; X target (use AB) (mm); Counts",nrun),100,-5,5);
   TH1F *hYtar_AB = new TH1F("hYtar_AB",Form("Run %d ; Y target (use AB) (mm); Counts",nrun),100,-5,5);
   TH1F *hXtar_AC = new TH1F("hXtar_AC",Form("Run %d ; X target (use AC) (mm); Counts",nrun),100,-5,5);
   TH1F *hYtar_AC = new TH1F("hYtar_AC",Form("Run %d ; Y target (use AC) (mm); Counts",nrun),100,-5,5);
   TH1F *hXtar_BC = new TH1F("hXtar_BC",Form("Run %d ; X target (use BC) (mm); Counts",nrun),100,-5,5);
   TH1F *hYtar_BC = new TH1F("hYtar_BC",Form("Run %d ; Y target (use BC) (mm); Counts",nrun),100,-5,5);
   //
   vector <Double_t > vHBcur;
   vector <Double_t > v8xpos;
   vector <Double_t > v8ypos;
   vector <Double_t > v8xcalc;
   vector <Double_t > v8ycalc;
   vector <Double_t > vTxpos;
   vector <Double_t > vTypos;
   vector <Double_t > v9xpos;
   vector <Double_t > v9ypos;
   vector <Double_t > v9xsum;
   vector <Double_t > v9ysum;
   vector <Double_t > v9xrot;
   vector <Double_t > v9yrot;
   vector <Double_t > v8xsum;
   vector <Double_t > v8ysum;
   vector <Double_t > v8xp;
   vector <Double_t > v8xm;
   vector <Double_t > v8yp;
   vector <Double_t > v8ym;
   vector <Double_t > v8xrot;
   vector <Double_t > v8yrot;
   vector <Double_t > v9xproj;
   vector <Double_t > v9yproj;
   // loop through data
   Double_t SA_zpos=370.82*10.; //mm
   Double_t SB_zpos=224.96*10.; //mm
   Double_t SC_zpos=129.3*10.; //mm
   Double_t SAX_off = -0.14;
   Double_t SAY_off = -0.17;
   Double_t SBX_off = -0.05;
   Double_t SBY_off = +0.47;
   Double_t SCX_off = -0.84;
   Double_t SCY_off = +0.43;
   Double_t SAX_scale = 1.03;
   Double_t SAY_scale = 0.98;
   Double_t SBX_scale = 1.23;
   Double_t SBY_scale = 1.21;
   Double_t SCX_scale = 0.94;
   Double_t SCY_scale = 0.86;
   Double_t bpm9Xp_off=27800.;
   Double_t bpm9Xm_off=27800.;
   Double_t bpm9Yp_off=25600.;
   Double_t bpm9Ym_off=25600.;
   Double_t bpm8Xp_off=18800.;
   Double_t bpm8Xm_off=18700.;
   Double_t bpm8Yp_off=18400.;
   Double_t bpm8Ym_off=18600.;
   Double_t thet=-45./180.*3.14;
   Double_t bpm8xprime;
   Double_t bpm8yprime;
   Double_t bpm9xprime;
   Double_t bpm9yprime;

Long64_t nentries = tsimc->GetEntries();
	for (int i = 0; i < nentries; i++) {
      		tsimc->GetEntry(i);
                if (ibcm1>0.1) {
		  hSAXraw->Fill(SAXraw);
		  hSAYraw->Fill(SAYraw);
		  hSBXraw->Fill(SBXraw);
		  hSBYraw->Fill(SBYraw);
		  hSCXraw->Fill(SCXraw);
		  hSCYraw->Fill(SCYraw);
		  SAXpos=SAX_off-SAX_scale*SAXraw;
		  SAYpos=SAY_off+SAY_scale*SAYraw;
		  SBXpos=SBX_off-SBX_scale*SBXraw;
		  SBYpos=SBY_off+SBY_scale*SBYraw;
		  SCXpos=SCX_off-SCX_scale*SCXraw;
		  SCYpos=SCY_off+SCY_scale*SCYraw;
		  hSAXpos->Fill(SAXpos);
		  hSAYpos->Fill(SAYpos);
		  hSBXpos->Fill(SBXpos);
		  hSBYpos->Fill(SBYpos);
		  hSCXpos->Fill(SCXpos);
		  hSCYpos->Fill(SCYpos);
		  Xtar_AB = SAXpos - (SAXpos-SBXpos)/(SA_zpos-SB_zpos)*SA_zpos;
		  Ytar_AB = SAYpos - (SAYpos-SBYpos)/(SA_zpos-SB_zpos)*SA_zpos;
		  Xtar_AC = SAXpos - (SAXpos-SCXpos)/(SA_zpos-SC_zpos)*SA_zpos;
		  Ytar_AC = SAYpos - (SAYpos-SCYpos)/(SA_zpos-SC_zpos)*SA_zpos;
		  Xtar_BC = SBXpos - (SBXpos-SCXpos)/(SB_zpos-SC_zpos)*SB_zpos;
		  Ytar_BC = SBYpos - (SBYpos-SCYpos)/(SB_zpos-SC_zpos)*SB_zpos;
                 if (HBcur!=0) {
                  vHBcur.push_back(HBcur);
                  v8xpos.push_back(bpm8Xpos);
                  v8xp.push_back(bpm8Xp);
                  v8xm.push_back(bpm8Xm);
                  v8yp.push_back(bpm8Yp);
                  v8ym.push_back(bpm8Ym);
                  v8ypos.push_back(bpm8Ypos);
                  v9xpos.push_back(bpm9Xpos);
                  v9ypos.push_back(bpm9Ypos);
                  v9xrot.push_back(bpm9Xrot);
                  v9yrot.push_back(bpm9Yrot);
                  v8xrot.push_back(bpm8Xrot);
                  v8yrot.push_back(bpm8Yrot);
                  v9xproj.push_back(Xtar_AB/SA_zpos*(SA_zpos+25900));
                  v9yproj.push_back(Ytar_AB/SA_zpos*(SA_zpos+25900));
		  bpm8xprime=200*(bpm8Xp-bpm8Xm+bpm8Xp_off-bpm8Xm_off)/(bpm8Xp+bpm8Xm+bpm8Xp_off+bpm8Xm_off);
		  bpm8xprime=200*(bpm8Xp-bpm8Yp+bpm8Xp_off-bpm8Xm_off)/(bpm8Xp+bpm8Ym+bpm8Xp_off+bpm8Xm_off);
		  v8xsum.push_back(bpm8xprime);
                  bpm8yprime=200*(bpm8Yp-bpm8Ym+bpm8Yp_off-bpm8Ym_off)/(bpm8Yp+bpm8Ym+bpm8Yp_off+bpm8Ym_off);
                  bpm8yprime=200*(bpm8Ym-bpm8Xm+bpm8Yp_off-bpm8Ym_off)/(bpm8Ym+bpm8Xm+bpm8Yp_off+bpm8Ym_off);
                  v8ysum.push_back(bpm8yprime);
		  //
		  bpm9xprime=200*(bpm9Xp-bpm9Xm+bpm9Xp_off-bpm9Xm_off)/(bpm9Xp+bpm9Xm+bpm9Xp_off+bpm9Xm_off);
		  bpm9xprime=200*(bpm9Xp-bpm9Yp+bpm9Xp_off-bpm9Xm_off)/(bpm9Xp+bpm9Ym+bpm9Xp_off+bpm9Xm_off);
		  v9xsum.push_back(bpm9xprime);
                  bpm9yprime=200*(bpm9Yp-bpm9Ym+bpm9Yp_off-bpm9Ym_off)/(bpm9Yp+bpm9Ym+bpm9Yp_off+bpm9Ym_off);
                  bpm9yprime=200*(bpm9Ym-bpm9Xm+bpm9Yp_off-bpm9Ym_off)/(bpm9Ym+bpm9Xm+bpm9Yp_off+bpm9Ym_off);
                  v9ysum.push_back(bpm9yprime);
		  //
                  vTxpos.push_back(Xtar_AB);
                  vTypos.push_back(Ytar_AB);
		  v8xcalc.push_back(bpm8xprime*cos(thet)+bpm8yprime*sin(thet));
		  v8ycalc.push_back(-bpm8xprime*sin(thet)+bpm8yprime*cos(thet));
		  }
		  hXtar_AB->Fill(Xtar_AB);
		  hYtar_AB->Fill(Ytar_AB);
		  hXtar_AC->Fill(Xtar_AC);
		  hYtar_AC->Fill(Ytar_AC);
		  hXtar_BC->Fill(Xtar_BC);
		  hYtar_BC->Fill(Ytar_BC);
		} else {
                    cout << bpm8Xrot << " "  << bpm8Xp << " "<< bpm8Xm << " " << bpm8Yrot<< " " << bpm8Yp << " " << bpm8Ym << endl;
		  }
	}
	// plot
          outputpdf="plots/"+basename+"_beam.pdf";
	// plot
TCanvas *c8x = new TCanvas("c8x", "bpm 8Xp ", 900,800);
  TGraph *gr_8xp_xm = new TGraph(vHBcur.size(),&(v8xp[0]),&(v8xm[0]));
  TGraph *gr_8xp_ym = new TGraph(vHBcur.size(),&(v8xp[0]),&(v8ym[0]));
  TGraph *gr_8xp_yp = new TGraph(vHBcur.size(),&(v8xp[0]),&(v8yp[0]));
  TGraph *gr_9xpos_8xpos = new TGraph(v9xpos.size(),&(v9xpos[0]),&(v8xpos[0]));
  TGraph *gr_9ypos_8ypos = new TGraph(v9ypos.size(),&(v9ypos[0]),&(v8ypos[0]));
  TGraph *gr_8xpos_9xproj = new TGraph(v8xpos.size(),&(v8xpos[0]),&(v9xproj[0]));
  TGraph *gr_8ypos_9yproj = new TGraph(v8ypos.size(),&(v8ypos[0]),&(v9yproj[0]));
  TGraph *gr_9xpos_9xproj = new TGraph(v9xpos.size(),&(v9xpos[0]),&(v9xproj[0]));
  TGraph *gr_9ypos_9yproj = new TGraph(v9ypos.size(),&(v9ypos[0]),&(v9yproj[0]));
  TGraph *gr_9xsum_9xproj = new TGraph(v9xsum.size(),&(v9xsum[0]),&(v9xproj[0]));
  TGraph *gr_9ysum_9yproj = new TGraph(v9ysum.size(),&(v9ysum[0]),&(v9yproj[0]));
c8x->Divide(2,3);
/* c8x->cd(1);
gr_8xp_xm->Draw();
gr_8xp_xm->SetTitle(" 8xp versus 8xm");
 c8x->cd(2);
gr_8xp_ym->Draw();
gr_8xp_ym->SetTitle(" 8xp versus 8ym");
 c8x->cd(3);
gr_8xp_yp->Draw();
gr_8xp_yp->SetTitle(" 8xp versus 8yp");
*/
 c8x->cd(1);
gr_8xpos_9xproj->Draw();
gr_8xpos_9xproj->SetTitle(" 9x Proj versus 8x pos");
 c8x->cd(2);
gr_8ypos_9yproj->Draw();
gr_8ypos_9yproj->SetTitle(" 9y proj versus 8y pos");
/* c8x->cd(3);
gr_9xpos_8xpos->Draw();
gr_9xpos_8xpos->SetTitle(" 8x POS versus 9X POS");
 c8x->cd(4);
gr_9ypos_8ypos->Draw();
gr_9ypos_8ypos->SetTitle(" 8y POS versus 9y POS");
*/
 c8x->cd(3);
gr_9xsum_9xproj->Draw();
gr_9xsum_9xproj->SetTitle(" 9x Proj versus 9x sum");
 c8x->cd(4);
gr_9ysum_9yproj->Draw();
gr_9ysum_9yproj->SetTitle(" 9y proj versus 9y sum");
 c8x->cd(5);
gr_9xpos_9xproj->Draw();
gr_9xpos_9xproj->SetTitle(" 9x Proj versus 9x pos");
 c8x->cd(6);
gr_9ypos_9yproj->Draw();
gr_9ypos_9yproj->SetTitle(" 9y proj versus 9y pos");
	// plot
TCanvas *c1 = new TCanvas("c1", "bpm 8 vs HB", 900,800);
  TGraph *gr_hbcur_8xpos = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8xpos[0]));
  TGraph *gr_hbcur_8ypos = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8ypos[0]));
  TGraph *gr_hbcur_8xcalc= new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8xcalc[0]));
  TGraph *gr_hbcur_8ycalc= new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8ycalc[0]));
  TGraph *gr_hbcur_8xrot = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8xrot[0]));
  TGraph *gr_hbcur_8yrot = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8yrot[0]));
  TGraph *gr_hbcur_8xsum = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8xsum[0]));
  TGraph *gr_hbcur_8ysum = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v8ysum[0]));
c1->Divide(2,4);
 c1->cd(1);
gr_hbcur_8xpos->Draw();
gr_hbcur_8xpos->SetTitle(" 8x POS versus HB cur");
 c1->cd(2);
gr_hbcur_8ypos->Draw();
gr_hbcur_8ypos->SetTitle(" 8y POS versus HB cur");
 c1->cd(3);
gr_hbcur_8xrot->Draw();
gr_hbcur_8xrot->SetTitle(" 8x rot versus HB cur");
 c1->cd(4);
gr_hbcur_8yrot->Draw();
gr_hbcur_8yrot->SetTitle(" 8y rot versus HB cur");
 c1->cd(5);
gr_hbcur_8xcalc->Draw();
gr_hbcur_8xcalc->SetTitle(" 8x calc versus HB cur");
 c1->cd(6);
gr_hbcur_8ycalc->Draw();
gr_hbcur_8ycalc->SetTitle(" 8y calc versus HB cur");
 c1->cd(7);
gr_hbcur_8xsum->Draw();
gr_hbcur_8xsum->SetTitle(" 8x sum versus HB cur");
 c1->cd(8);
gr_hbcur_8ysum->Draw();
gr_hbcur_8ysum->SetTitle(" 8y sum versus HB cur");
//
TCanvas *c2 = new TCanvas("c2", "bpm 9 vs HB", 900,800);
  TGraph *gr_hbcur_9xsum = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9xsum[0]));
  TGraph *gr_hbcur_9ysum = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9ysum[0]));
  TGraph *gr_hbcur_9x = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9xpos[0]));
  TGraph *gr_hbcur_9y = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9ypos[0]));
  TGraph *gr_hbcur_9xproj = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9xproj[0]));
  TGraph *gr_hbcur_9yproj = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9yproj[0]));
  TGraph *gr_hbcur_9xrot = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9xrot[0]));
  TGraph *gr_hbcur_9yrot = new TGraph(vHBcur.size(),&(vHBcur[0]),&(v9yrot[0]));
c2->Divide(2,3);
 c2->cd(1);
gr_hbcur_9xsum->Draw();
gr_hbcur_9xsum->SetTitle(" 9x SUM versus HB cur");
 c2->cd(2);
gr_hbcur_9ysum->Draw();
gr_hbcur_9ysum->SetTitle(" 9y SUM versus HB cur");
 c2->cd(3);
gr_hbcur_9xproj->Draw();
gr_hbcur_9xproj->SetTitle(" 9x proj versus HB cur");
 c2->cd(4);
gr_hbcur_9yproj->Draw();
gr_hbcur_9yproj->SetTitle(" 9y proj versus HB cur");
 c2->cd(5);
gr_hbcur_9xrot->Draw();
gr_hbcur_9xrot->SetTitle(" 9x rot versus HB cur");
 c2->cd(6);
gr_hbcur_9yrot->Draw();
gr_hbcur_9yrot->SetTitle(" 9y rot versus HB cur");
c2->Print(outputpdf);

//
}
