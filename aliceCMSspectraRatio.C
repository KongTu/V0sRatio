#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include <stdio.h>
#include <iomanip>
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;

void aliceCMSspectraRatio(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/aliceQuantativeComparison/aliceQuantativeComparison_90_120_vertexReweight_v4_epos_pt.root");

	TH1D* ksSpectra_cms = (TH1D*)file1->Get("ksSpectra_alice");
	TH1D* laSpectra_cms = (TH1D*)file1->Get("laSpectra_alice");

	double ks_ptbins[35] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2,4.6,5.0,5.5,6.0,8.0};
    double ks_pTbinsBound[35] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,33,36,39,42,46,50,55,60,80};
    
    double ks_binwidth[34] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.4,0.4,0.5,0.5,2.0};
    double ks_ptbincenter[34] = { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
    0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
    1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.15, 3.45, 3.75, 4.05, 
    4.4, 4.8, 5.25, 5.75, 7.0 };

    double la_ptbins[21] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.7,4.2,5.0,6.0,8.0};
    double la_ptbincenter[20] = { 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.3, 1.5, 1.7, 
    1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.45, 3.95, 4.6, 5.5, 
    7.0 };
    double la_pTbinsBound[21] = {6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,32,37,42,50,60,80};
    double la_binwidth[20] = {0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.5,0.5,0.6,1.0,2.0};

    double yields_ks1[] = { 1.532, 1.373, 1.22, 1.026, 0.8654, 0.7002, 0.5571, 0.4442, 0.3476, 
	0.2745, 0.2161, 0.1731, 0.1342, 0.1068, 0.0867, 0.06848, 0.0553, 0.04467, 0.03537, 
	0.02874, 0.02089, 0.01407, 0.009662, 0.006504, 0.004485, 0.002855, 0.001833, 0.00112, 7.169E-4, 
	4.67E-4, 2.456E-4, 1.467E-4, 9.063E-5, 2.868E-5 };

	double yields_ks_error1[] = { 0.161, 0.027, 0.015, 0.011, 0.0076, 0.0053, 0.0039, 0.003, 0.0023, 
    0.0019, 0.0015, 0.0013, 0.0011, 9.0E-4, 7.9E-4, 6.7E-4, 5.8E-4, 5.0E-4, 4.3E-4, 
    3.8E-4, 2.2E-4, 1.7E-4, 1.34E-4, 1.06E-4, 8.5E-5, 5.2E-5, 4.0E-5, 3.0E-5, 2.28E-5, 
    1.57E-5, 1.09E-5, 7.1E-6, 5.37E-6, 1.39E-6 };


    double yields_la1[] = { 0.1833, 0.1867, 0.1737, 0.1574, 0.1356, 0.1232, 0.1006, 0.07763, 0.05522, 
    0.03978, 0.02885, 0.02024, 0.01422, 0.009873, 0.005879, 0.002652, 0.001182, 4.378E-4, 1.2E-4, 
    2.194E-5 };

  	double yields_la_error1[] = { 0.0035, 0.0028, 0.0024, 0.0021, 0.0018, 0.0017, 0.001, 8.2E-4, 6.3E-4, 
    5.0E-4, 4.0E-4, 3.2E-4, 2.6E-4, 2.04E-4, 1.08E-4, 6.0E-5, 3.7E-5, 1.62E-5, 6.5E-6, 
    1.62E-6 };

	TH1D* ksSpectra_alice = new TH1D("ksSpectra_new","ksSpectra_new",34,ks_ptbins);
	TH1D* laSpectra_alice = new TH1D("laSpectra_new","laSpectra_new",20,la_ptbins);

	for(int pt = 0; pt < 34; pt++){

		ksSpectra_alice->SetBinContent(pt+1, yields_ks1[pt] );
			ksSpectra_alice->SetBinError(pt+1, yields_ks_error1[pt] );
	}

	for(pt = 0; pt < 20; pt++){

		laSpectra_alice->SetBinContent(pt+1, yields_la1[pt] );
			laSpectra_alice->SetBinError(pt+1, yields_la_error1[pt] );
	}



	TCanvas* c1 = new TCanvas();
	c1->Divide(1,2,0,0);

	c1->cd(1);

	ksSpectra_cms->Divide(ksSpectra_alice);
	ksSpectra_cms->SetMarkerStyle(20);
	ksSpectra_cms->Draw("P");

	c1->cd(2);

	laSpectra_cms->Divide(laSpectra_alice);
	laSpectra_cms->SetMarkerStyle(20);
	laSpectra_cms->Draw("P");



}