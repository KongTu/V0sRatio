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


void testVTXreweight(){


	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
    TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/HIJINGeffTable_withVTXreweight_18M_Nov5_rapidity_v1_28ks_pTbins.root");

    TH1D* ks_eff_rpy[5];
    TH1D* la_eff_rpy[5];

    TH1D* ks_eff_rpy_new[5];
    TH1D* la_eff_rpy_new[5];

    double ks_eff[5][28];
    double la_eff[5][20];

    double ks_eff_err[5][28];
    double la_eff_err[5][20];

    stringstream ksName;
    stringstream laName;

    for (int rpy = 0; rpy < 5; rpy++){

        ksName.str("");
        laName.str("");

        ksName << "ks_eff_rpy_";
        ksName << rpy+1;

        laName << "la_eff_rpy_";
        laName << rpy+1;

        ks_eff_rpy[rpy] = (TH1D*) file->Get( ksName.str().c_str() );
        la_eff_rpy[rpy] = (TH1D*) file->Get( laName.str().c_str() );

        ksName.str("");
        laName.str("");

        ksName << "ks_eff_rpy_vtx_";
        ksName << rpy+1;

        laName << "la_eff_rpy_vtx_";
        laName << rpy+1;

        ks_eff_rpy_new[rpy] = (TH1D*) file1->Get( ksName.str().c_str() );
        la_eff_rpy_new[rpy] = (TH1D*) file1->Get( laName.str().c_str() );

    }

    TLatex* ratio[5];
	ratio[0] = new TLatex(1.5,0.95,"-2.87 < y < -1.8");
	ratio[1] = new TLatex(1.5,0.95,"-1.8 < y < -0.9");
	ratio[2] = new TLatex(1.5,0.95,"-0.9 < y < 0");
	ratio[3] = new TLatex(1.5,0.95,"0 < y < 0.93");
	ratio[4] = new TLatex(1.5,0.95,"0.93 < y < 1.93");

    TLine* l1 = new TLine(0,1,9.0,1.0);
    l1->SetLineWidth(2);
    l1->SetLineColor(kRed);
    l1->SetLineStyle(2);


    TCanvas* c1 = new TCanvas();
    c1->Divide(2,3,0,0);


    for(rpy = 0; rpy < 5; rpy++){

    	c1->cd(rpy+1);

    	ks_eff_rpy_new[rpy]->Divide( ks_eff_rpy[rpy] );
    	ks_eff_rpy_new[rpy]->SetMarkerStyle(20);
    	ks_eff_rpy_new[rpy]->GetYaxis()->SetRangeUser(0.9,1.40);
        ks_eff_rpy_new[rpy]->GetXaxis()->SetRangeUser(0,9.0);
    	ks_eff_rpy_new[rpy]->SetStats(kFALSE);
    	ks_eff_rpy_new[rpy]->SetLineColor(kBlack);
    	ks_eff_rpy_new[rpy]->SetTitle("K^{0}_{s}");
        ks_eff_rpy_new[rpy]->SetYTitle("vtx/no_vtx");
        ks_eff_rpy_new[rpy]->GetYaxis()->SetTitleSize(0.06);
        ks_eff_rpy_new[rpy]->GetYaxis()->SetTitleOffset(0.7);
        ks_eff_rpy_new[rpy]->GetYaxis()->SetLabelSize(0.07);

        ks_eff_rpy_new[rpy]->GetYaxis()->SetNdivisions(9,5,0);
        ks_eff_rpy_new[rpy]->SetXTitle("pT(GeV/c)");
        ks_eff_rpy_new[rpy]->GetXaxis()->SetTitleSize(0.06);
        ks_eff_rpy_new[rpy]->GetXaxis()->SetTitleOffset(0.7);
        ks_eff_rpy_new[rpy]->GetXaxis()->SetLabelSize(0.07);
    	ks_eff_rpy_new[rpy]->Draw("P");
        ratio[rpy]->SetTextSize(0.1);
    	ratio[rpy]->Draw("same");
        l1->Draw("same");


    }

    TCanvas* c2 = new TCanvas();
    c2->Divide(2,3,0,0);


    for(rpy = 0; rpy < 5; rpy++){

    	c2->cd(rpy+1);

    	la_eff_rpy_new[rpy]->Divide( la_eff_rpy[rpy] );
        la_eff_rpy_new[rpy]->SetMarkerStyle(20);
        la_eff_rpy_new[rpy]->GetYaxis()->SetRangeUser(0.9,1.40);
        la_eff_rpy_new[rpy]->GetXaxis()->SetRangeUser(0,9.0);
        la_eff_rpy_new[rpy]->SetStats(kFALSE);
        la_eff_rpy_new[rpy]->SetLineColor(kBlack);
        la_eff_rpy_new[rpy]->SetTitle("#Lambda/#bar{#Lambda}");
        la_eff_rpy_new[rpy]->SetYTitle("vtx/no_vtx");
        la_eff_rpy_new[rpy]->GetYaxis()->SetTitleSize(0.06);
        la_eff_rpy_new[rpy]->GetYaxis()->SetTitleOffset(0.7);
        la_eff_rpy_new[rpy]->GetYaxis()->SetLabelSize(0.07);
        la_eff_rpy_new[rpy]->GetYaxis()->SetNdivisions(9,5,0);
        la_eff_rpy_new[rpy]->SetXTitle("pT(GeV/c)");
        la_eff_rpy_new[rpy]->GetXaxis()->SetTitleSize(0.06);
        la_eff_rpy_new[rpy]->GetXaxis()->SetTitleOffset(0.7);
        la_eff_rpy_new[rpy]->GetXaxis()->SetLabelSize(0.07);
        la_eff_rpy_new[rpy]->Draw("P");
        ratio[rpy]->SetTextSize(0.1);
        ratio[rpy]->Draw("same");
        l1->Draw("same");


    }


}