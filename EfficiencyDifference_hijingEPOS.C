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

void EfficiencyDifference_hijingEPOS(){

	gStyle->SetErrorX(0);

	double ks_ptbins[30] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[29] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.05,1.15,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  	double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/eposEfficiencyRapidityTable/EPOSeffKongNew2DTable_5M_Nov3_rapidity_v1_28ks_pTbins.root");

	TH1D* ks_Hijing_rpy[5];
	TH1D* la_Hijing_rpy[5];

	TH1D* ks_Epos_rpy[5];
	TH1D* la_Epos_rpy[5];

	stringstream ksName;
	stringstream laName;

	for (int rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_rpy_";
		ksName << rpy+1;

		laName << "la_eff_rpy_";
		laName << rpy+1;

		ks_Hijing_rpy[rpy] = (TH1D*) file1->Get( ksName.str().c_str() );
		la_Hijing_rpy[rpy] = (TH1D*) file1->Get( laName.str().c_str() );

		ks_Epos_rpy[rpy] = (TH1D*) file2->Get( ksName.str().c_str() );
		la_Epos_rpy[rpy] = (TH1D*) file2->Get( laName.str().c_str() );

	}

	TLatex* r6[5];
    r6[0] = new TLatex(5.5,1.6,"-2.87 < y < -1.8");
    r6[1] = new TLatex(5.5,1.6,"-1.8 < y < -0.9");
    r6[2] = new TLatex(5.5,1.6,"-0.9 < y < 0");
    r6[3] = new TLatex(5.5,1.6,"0 < y < 0.93");
    r6[4] = new TLatex(5.5,1.6,"0.93 < y < 1.93");

    TLatex* r1 = new TLatex(1.0,1.6,"beforeSmoothing");
    TLatex* r4 = new TLatex(1.0,1.6,"afterSmoothing");
    TLatex* r2 = new TLatex(4.0,1.6,"K^{0}_{s}");
    TLatex* r3 = new TLatex(4.0,1.6,"#Lambda/#bar{#Lambda}");

	TCanvas* c1 = new TCanvas();
	c1->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c1->cd(rpy+1);

		ks_Hijing_rpy[rpy]->Divide(ks_Epos_rpy[rpy]);
		ks_Hijing_rpy[rpy]->SetMarkerStyle(20);
		ks_Hijing_rpy[rpy]->GetYaxis()->SetRangeUser(0,2.0);
		ks_Hijing_rpy[rpy]->SetStats(kFALSE);
		ks_Hijing_rpy[rpy]->SetTitle("");
		ks_Hijing_rpy[rpy]->SetLineColor(kBlack);
		ks_Hijing_rpy[rpy]->Draw("P");
		r6[rpy]->Draw("same");
		r2->Draw("same");

	}

	TCanvas* c2 = new TCanvas();
	c2->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c2->cd(rpy+1);

		la_Hijing_rpy[rpy]->Divide(la_Epos_rpy[rpy]);
		la_Hijing_rpy[rpy]->SetMarkerStyle(20);
		la_Hijing_rpy[rpy]->GetYaxis()->SetRangeUser(0,2.0);
		la_Hijing_rpy[rpy]->SetStats(kFALSE);
		la_Hijing_rpy[rpy]->SetTitle("");
		la_Hijing_rpy[rpy]->SetLineColor(kBlack);
		la_Hijing_rpy[rpy]->Draw("P");
		r6[rpy]->Draw("same");
		r3->Draw("same");

	}










}