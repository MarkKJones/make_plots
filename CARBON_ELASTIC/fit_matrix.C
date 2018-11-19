#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <TSystem.h>
#include <TFile.h>
#include <TNtuple.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TMath.h>
#include <TBox.h>
#include <TPolyLine.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>
void fit_matrix() {
  gROOT->Reset();
  gStyle->SetOptStat(0);
 gStyle->SetTitleOffset(1.,"Y");
 gStyle->SetTitleOffset(.7,"X");
 gStyle->SetLabelSize(0.04,"XY");
 gStyle->SetTitleSize(0.05,"XY");
 gStyle->SetPadLeftMargin(0.17);
 // CHain files together
    TString inputroot;
    TChain *fdelta;
    fdelta = new TChain("fit");
    for (Int_t nr=1657;nr<1676;nr++) {
      TString inputroot=Form("shms_replay_matrixopt_%d_orig_%d_tree.root",nr,nr);
      if (gSystem->FindFile("hist",inputroot)) {
		cout << " file = " << inputroot << endl;
		//	TFile fch = new TFile("hist/"+inputroot);
       	fdelta->Add(inputroot);	       
      } else {
	cout << " No file = " << inputroot << endl;
      }
    }
    fdelta->Merge("hist/chain_carbon.root");
    //
    TFile *fsimc = new TFile("hist/chain_carbon.root"); 
    TTree *tsimc = (TTree*) fsimc->Get("fit");
 //
    Float_t xfp,xpfp,yfp,ypfp,delta,ex;
    Float_t xptar,yptar,ytar;
    Float_t xptarT,yptarT,ytarT;
    Float_t nx,ny,nd;
    tsimc->SetBranchAddress("nx",&nx);
    tsimc->SetBranchAddress("ny",&ny);
    tsimc->SetBranchAddress("nd",&nd);
    tsimc->SetBranchAddress("xfp",&xfp);
    tsimc->SetBranchAddress("xpfp",&xpfp);
    tsimc->SetBranchAddress("yfp",&yfp);
    tsimc->SetBranchAddress("ypfp",&ypfp);
    tsimc->SetBranchAddress("delta",&delta);
    tsimc->SetBranchAddress("ex",&ex);
    tsimc->SetBranchAddress("xptar",&xptar);
    tsimc->SetBranchAddress("yptar",&yptar);
    tsimc->SetBranchAddress("ytar",&ytar);
    tsimc->SetBranchAddress("xptarT",&xptarT);
    tsimc->SetBranchAddress("yptarT",&yptarT);
    tsimc->SetBranchAddress("ytarT",&ytarT);
    //
    TString outputhist;
    TObjArray HList(0);
    TH2F *hxs_ys_old[18];
    TH2F *hxs_ys_new[18];
    
    TH1F *hxptardiffold[18][11][11];
    TH1F *hxptardiffnew[18][11][11];
    TH1F *hyptardiffold[18][11][11];
    TH1F *hyptardiffnew[18][11][11];
    TH1F *hytardiffold[18][11][11];
    TH1F *hytardiffnew[18][11][11];
    for (int nnd=0;nnd<18;nnd++) {
      hxs_ys_old[nnd] = new TH2F(Form("hxs_ys_%d_old",nnd), Form("Old fit Delta run = %d; Y_sieve ; X_sieve",-12+2*nnd), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys_old[nnd]);
      hxs_ys_new[nnd] = new TH2F(Form("hxs_ys_%d_new",nnd), Form("New fit Delta run = %d; Y_sieve ; X_sieve",-12+2*nnd), 100,-10,10.,120,-15,15);
          HList.Add(hxs_ys_new[nnd]);
    for (int nnx=0;nnx<11;nnx++) {
    for (int nny=0;nny<11;nny++) {
     
     hxptardiffnew[nnd][nnx][nny] = new TH1F(Form("hxptardiffnew_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Xptar -Xptar_true",nnd,nnx,nny),100,-.02,.02);
    HList.Add(hxptardiffnew[nnd][nnx][nny]);
      hxptardiffold[nnd][nnx][nny] = new TH1F(Form("hxptardiffold_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Xptar -Xptar_true",nnd,nnx,nny),100,-0.02,.02);
    HList.Add(hxptardiffold[nnd][nnx][nny]);
      hyptardiffnew[nnd][nnx][nny] = new TH1F(Form("hyptardiffnew_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Yptar -Yptar_true",nnd,nnx,nny),100,-.02,.02);
    HList.Add(hyptardiffnew[nnd][nnx][nny]);
      hyptardiffold[nnd][nnx][nny] = new TH1F(Form("hyptardiffold_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Yptar -Yptar_true",nnd,nnx,nny),100,-0.02,.02);
    HList.Add(hyptardiffold[nnd][nnx][nny]);
      hytardiffnew[nnd][nnx][nny] = new TH1F(Form("hytardiffnew_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Ytar -Ytar_true",nnd,nnx,nny),100,-1.,1.);
    HList.Add(hytardiffnew[nnd][nnx][nny]);
      hytardiffold[nnd][nnx][nny] = new TH1F(Form("hytardiffold_%d_%d_%d",nnd,nnx,nny),Form(" Nd=%d Nx=%d Ny=%d; Ytar -Ytar_true",nnd,nnx,nny),100,-1.,1.);
    HList.Add(hytardiffold[nnd][nnx][nny]);
    }
    }
    }
    //
  //
		     //
   //
      //
    //string newcoeffsfilename="newfit.dat";
   //    string oldcoeffsfilename="shms-2017-optimized_2plow_delta_cosy_quads_p18.dat";
    string oldcoeffsfilename="CARBON_ELASTIC/shms-2017-optimized_delta_newfit3.dat";
    string newcoeffsfilename="CARBON_ELASTIC/newfit.dat";
  ifstream ifile(oldcoeffsfilename.c_str());
  if(!ifile.is_open()) {
    cout << "error opening reconstruction coefficient file = " << oldcoeffsfilename << endl;
    return;    
  }
   ofstream newcoeffsfile(newcoeffsfilename.c_str());
  vector<double> xptarcoeffs_old;
  vector<double> yptarcoeffs_old;
  vector<double> ytarcoeffs_old;
  vector<double> deltacoeffs_old;
  vector<int> xfpexpon_old;
  vector<int> xpfpexpon_old;
  vector<int> yfpexpon_old;
  vector<int> ypfpexpon_old;
  vector<int> xtarexpon_old;

  vector<double> xptarcoeffs_fit;
  vector<double> yptarcoeffs_fit;
  vector<double> ytarcoeffs_fit;
  vector<double> deltacoeffs_fit;
  vector<int> xfpexpon_fit;
  vector<int> xpfpexpon_fit;
  vector<int> yfpexpon_fit;
  vector<int> ypfpexpon_fit;
  vector<int> xtarexpon_fit;

  vector<double> xptarcoeffs_new;
  vector<double> yptarcoeffs_new;
  vector<double> ytarcoeffs_new;
  vector<double> deltacoeffs_new;
  vector<int> xfpexpon_new;
  vector<int> xpfpexpon_new;
  vector<int> yfpexpon_new;
  vector<int> ypfpexpon_new;
  vector<int> xtarexpon_new;
  vector<double> xtartrue,ytartrue,xptartrue,yptartrue;
  vector<double> ytarold,xptarold,yptarold,deltaold;
  vector<double> xfptrue,yfptrue,xpfptrue,ypfptrue;
  vector<int> nxtrue,nytrue,ndtrue;
  TString currentline;
  int num_recon_terms_old;
  int num_recon_terms_fit;
  num_recon_terms_old = 0;
  num_recon_terms_fit = 0;
  int nfit=0,npar,nfit_max=1000000,npar_final=0,max_order=5,norder,nmax_delta=40000;
  string line="!";
  int good=1;
  while(good && line[0]=='!') {
    good = getline(ifile,line).good();
     cout << line << endl;
  } 
  while(good && line.compare(0,4," ---")!=0) {
    good = getline(ifile,line).good();
     cout << line << endl;
  }
  line=" ";
  good = getline(ifile,line).good();
  Double_t c1,c2,c3,c4;
  Int_t e1,e2,e3,e4,e5;
  while(good && line.compare(0,4," ---")!=0) {
    sscanf(line.c_str()," %le %le %le %le %1d%1d%1d%1d%1d",&c1,&c2,&c3,&c4,&e1,&e2,&e3,&e4,&e5);
    xptarcoeffs_old.push_back(c1);
    ytarcoeffs_old.push_back(c2);
    yptarcoeffs_old.push_back(c3);
    deltacoeffs_old.push_back(c4);
     xfpexpon_old.push_back(e1);
    xpfpexpon_old.push_back(e2);
    yfpexpon_old.push_back(e3);
    ypfpexpon_old.push_back(e4);
    xtarexpon_old.push_back(e5);       
    num_recon_terms_old++;
    good = getline(ifile,line).good();
     cout << line << endl;
    norder= e1+e2+e3+e4;
    if (norder <= max_order && e5==0) {
    xptarcoeffs_fit.push_back(c1);
    ytarcoeffs_fit.push_back(c2);
    yptarcoeffs_fit.push_back(c3);
    deltacoeffs_fit.push_back(c4);
    xfpexpon_fit.push_back(e1);
    xpfpexpon_fit.push_back(e2);
    yfpexpon_fit.push_back(e3);
    ypfpexpon_fit.push_back(e4);
    xtarexpon_fit.push_back(e5);
    num_recon_terms_fit++;
    }
  }
  cout << "num recon terms in OLD matrix = " << num_recon_terms_old << endl;
  cout << "num recon terms in fit matrix = " << num_recon_terms_fit << endl;
   npar= num_recon_terms_fit ;
    //
  TVectorD b_xptar(npar);
  TVectorD b_yptar(npar);
  TVectorD b_ytar(npar);
  TMatrixD lambda(npar,nfit_max);
  TMatrixD Ay(npar,npar);
   //
  Int_t counts[18][11][11];
   	for (int ndd= 0; ndd < 18; ndd++) {
   	for (int nxx = 0; nxx < 11; nxx++) {
   	for (int nyy = 0; nyy < 11; nyy++) {
	  counts[ndd][nxx][nyy]=0;
	}}  }
  //
    Long64_t nentries = tsimc->GetEntriesFast();
   cout << " nent = " << nentries << endl;
   Int_t delta_event_count[25];
   //nentries=100;
   	for (int i = 0; i < nentries; i++) {
    		tsimc->GetEntry(i);
    if (i%500000==0) cout << " Entry = " << i << endl;
    //Double_t ytar = 0.0,yptar=0.0,xptar=0.0,del=0.0;
          Double_t etemp;
          Double_t xtar=0; // temp set xtar=o for all events.
          /*for( int icoeffold=0; icoeffold<num_recon_terms_old; icoeffold++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_old[icoeffold] ) * 
	  pow( yfp / 100.0, yfpexpon_old[icoeffold] ) * 
	  pow( xpfp, xpfpexpon_old[icoeffold] ) * 
	  pow( ypfp, ypfpexpon_old[icoeffold] ) * 
	  pow( xtar/100., xtarexpon_old[icoeffold] );
        	ytar += ytarcoeffs_old[icoeffold] * etemp;
        	yptar += yptarcoeffs_old[icoeffold] * etemp;
        	xptar += xptarcoeffs_old[icoeffold] * etemp;
        	del += deltacoeffs_old[icoeffold] * etemp;
	  } // for icoeffold loop
	  */
   
          if (counts[int(nd)][int(nx)][int(ny)] < 1000) { 
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
        	etemp= 
	  pow( xfp / 100.0, xfpexpon_fit[icoeff_fit] ) * 
	  pow( yfp / 100.0, yfpexpon_fit[icoeff_fit] ) * 
	  pow( xpfp, xpfpexpon_fit[icoeff_fit] ) * 
	  pow( ypfp, ypfpexpon_fit[icoeff_fit] ) * 
	  pow( xtar/100., xtarexpon_fit[icoeff_fit] );
 		 if (nfit < nfit_max ) {
              lambda[icoeff_fit][nfit] = etemp;
	      b_xptar[icoeff_fit] += (xptarT) * etemp;
	      b_yptar[icoeff_fit] += (yptarT) * etemp;
	      //	      b_ytar[icoeff_fit] += (0.001/100.) * etemp;
              }
	  } // for icoeff_fit loop
	  if (nfit < nfit_max  ) {
	    //delta_event_count[delta_run_index]++;
	    counts[int(nd)][int(nx)][int(ny)]++;
	      nfit++;
  	    xfptrue.push_back( xfp );
	    yfptrue.push_back( yfp );
	    xpfptrue.push_back( xpfp );
	    ypfptrue.push_back( ypfp );
	    xtartrue.push_back( xtar );
	    xptartrue.push_back( xptarT  );
	    yptartrue.push_back( yptarT  );
	    ytartrue.push_back( ytarT  );
	    xptarold.push_back( xptar  );
	    yptarold.push_back( yptar  );
	    ytarold.push_back( ytar  );
	    deltaold.push_back( delta  );
	    ndtrue.push_back( int(nd)  );
	    nxtrue.push_back( int(nx)  );
	    nytrue.push_back( int(ny)  );
	  } // nfit < nfit_max
	  } // if deltatemp-delta_calc/100
	} // event loop
	cout << " nfit = " << nfit << " max = " << nfit_max<< endl;
 for(int i=0; i<npar; i++){
    for(int j=0; j<npar; j++){
      Ay[i][j] = 0.0;
    }
 }
   for( int ifit=0; ifit<nfit; ifit++){
    if( ifit % 100000 == 0 ) cout << ifit << endl;
    for( int ipar=0; ipar<npar; ipar++){
      for( int jpar=0; jpar<npar; jpar++){
      	Ay[ipar][jpar] += lambda[ipar][ifit] * lambda[jpar][ifit];
       }
    }
  }
   //
   TDecompSVD Ay_svd(Ay);
  bool ok;
 ok = Ay_svd.Solve( b_xptar );
  cout << "Xptar solution ok = " << ok << endl;
 ok = Ay_svd.Solve( b_yptar );
  cout << "Yptar solution ok = " << ok << endl;
  // ok = Ay_svd.Solve( b_ytar );
  //cout << "Ytar solution ok = " << ok << endl;
  //
  for( int ifit=0; ifit<nfit; ifit++){
          Double_t ytarnew = 0.0,yptarnew=0.0,xptarnew=0.0,deltanew=0.0;
	  Double_t etemp;
     for( int ipar=0; ipar<npar; ipar++){
       etemp=lambda[ipar][ifit];
        	xptarnew += b_xptar[ipar] * etemp;
        	yptarnew += b_yptar[ipar] * etemp;
		// 	ytarnew += b_ytar[ipar] * etemp;
    }
     // cout << ndtrue[ifit] << " " << nxtrue[ifit] << " " << nytrue[ifit] << endl;
     Double_t xsold=xptarold[ifit]*253.;
     Double_t xsnew=xptarnew*253.;
     //     Double_t ystemp=ytarold[ifit]-(0.019+40.*.01*.052)*deltaold[ifit]+(0.00019+40*.01*.00052)*deltaold[ifit]*deltaold[ifit];
          Double_t ystemp=-(0.019+40.*.01*.052)*deltaold[ifit]+(0.00019+40*.01*.00052)*deltaold[ifit]*deltaold[ifit];
    Double_t ysold=ystemp+yptarold[ifit]*253.;
     Double_t ysnew=ystemp+yptarnew*253.;
     hxs_ys_old[ndtrue[ifit]]->Fill(ysold,xsold);
     hxs_ys_new[ndtrue[ifit]]->Fill(ysnew,xsnew);
     hxptardiffnew[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(xptarnew-xptartrue[ifit]);
     hxptardiffold[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(xptarold[ifit]-xptartrue[ifit]);
     hyptardiffnew[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(yptarnew-yptartrue[ifit]);
     hyptardiffold[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(yptarold[ifit]-yptartrue[ifit]);
     //     hytardiffnew[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(100*(ytarnew-ytartrue[ifit]));
     //hytardiffold[ndtrue[ifit]][nxtrue[ifit]][nytrue[ifit]]->Fill(100*(ytarold[ifit]-ytartrue[ifit]));
  } // 
  //
  // write out coeff
  char coeffstring[100];
  Double_t tt;
  cout << "writing new coeffs file" << endl;
  newcoeffsfile << " ---------------------------------------------" << endl;
          for( int icoeff_fit=0; icoeff_fit<num_recon_terms_fit; icoeff_fit++ ){
      newcoeffsfile << " ";
      tt=b_xptar[icoeff_fit];
      sprintf( coeffstring, "%16.9g", tt );
      newcoeffsfile << coeffstring; 
      //      newcoeffsfile << " ";
      sprintf( coeffstring, "%16.9g", ytarcoeffs_old[icoeff_fit] );
      newcoeffsfile << coeffstring;
      sprintf( coeffstring, "%16.9g", b_yptar[icoeff_fit] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      sprintf( coeffstring, "%16.9g", deltacoeffs_old[icoeff_fit] );
      //newcoeffsfile << " ";
      newcoeffsfile << coeffstring; 
      newcoeffsfile << " ";
	newcoeffsfile << setw(1) << setprecision(1) << xfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << xpfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << yfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << ypfpexpon_fit[icoeff_fit]; 
	newcoeffsfile << setw(1) << setprecision(1) << xtarexpon_fit[icoeff_fit]; 
      newcoeffsfile << endl;

	  }
  newcoeffsfile << " ---------------------------------------------" << endl;

  newcoeffsfile.close();
  cout << "wrote new coeffs file" << endl;
  //
 TFile hsimc("hist/fit_carbon_result_fit3.root","recreate");
  HList.Write();
  //
}
