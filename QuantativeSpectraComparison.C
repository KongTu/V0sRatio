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

void QuantativeSpectraComparison(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy_new/new8Multbins_vtxReweight_FullStats_v1_pt.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy_new/new8Multbins_noSmooth_FullStats_v1_pt.root");
	TFile* file3 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy_new/new8Multbins_pSmooth_FullStats_v1_pt.root");

	TH1D* ksSpectra_vtx[8];
	TH1D* laSpectra_vtx[8];

	TH1D* ksSpectra_noSmooth[8];
	TH1D* laSpectra_noSmooth[8];

	TH1D* ksSpectra_pSmooth[8];
	TH1D* laSpectra_pSmooth[8];

	stringstream ksName;
	stringstream laName;

	for(int mult = 0; mult < 8; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_vtx_";
		ksName << mult+1;
		laName << "laSpectra_vtx_";
		laName << mult+1;

		ksSpectra_vtx[mult] = (TH1D*)file1->Get( ksName.str().c_str() );
		laSpectra_vtx[mult] = (TH1D*)file1->Get( laName.str().c_str() );

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_noSmooth_";
		ksName << mult+1;
		laName << "laSpectra_noSmooth_";
		laName << mult+1;

		ksSpectra_noSmooth[mult] = (TH1D*)file2->Get( ksName.str().c_str() );
		laSpectra_noSmooth[mult] = (TH1D*)file2->Get( laName.str().c_str() );

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_pSmooth_";
		ksName << mult+1;
		laName << "laSpectra_pSmooth_";
		laName << mult+1;

		ksSpectra_pSmooth[mult] = (TH1D*)file3->Get( ksName.str().c_str() );
		laSpectra_pSmooth[mult] = (TH1D*)file3->Get( laName.str().c_str() );


	}

	TLatex* ratio[8];
	ratio[0] = new TLatex(1.5,1.1,"0 < N^{offline}_{trk} < 35");
	ratio[1] = new TLatex(1.5,1.1,"35 < N^{offline}_{trk} < 60");
	ratio[2] = new TLatex(1.5,1.1,"60 < N^{offline}_{trk} < 90");
	ratio[3] = new TLatex(1.5,1.1,"90 < N^{offline}_{trk} < 120");
	ratio[4] = new TLatex(1.5,1.1,"120 < N^{offline}_{trk} < 150");
	ratio[5] = new TLatex(1.5,1.1,"150 < N^{offline}_{trk} < 185");
	ratio[6] = new TLatex(1.5,1.1,"185 < N^{offline}_{trk} < 220");
	ratio[7] = new TLatex(1.5,1.1,"220 < N^{offline}_{trk}");

 
	TCanvas* c1 = new TCanvas();
	c1->Divide(2,2,0,0);

	c1->cd(1);

	ksSpectra_vtx[0]->Divide( ksSpectra_noSmooth[0] );
	ksSpectra_vtx[0]->SetMarkerStyle(20);
	ksSpectra_vtx[0]->GetYaxis()->SetRangeUser(0.6,1.3);
	ksSpectra_vtx[0]->SetStats(kFALSE);
	ksSpectra_vtx[0]->SetLineColor(kBlack);
	ksSpectra_vtx[0]->SetTitle("K^{0}_{s} vertex/noSmooth");
	ksSpectra_vtx[0]->Draw("P");

	c1->cd(2);

	ksSpectra_pSmooth[0]->Divide( ksSpectra_noSmooth[0] );
	ksSpectra_pSmooth[0]->SetMarkerStyle(20);
	ksSpectra_pSmooth[0]->GetYaxis()->SetRangeUser(0.6,1.3);
	ksSpectra_pSmooth[0]->SetStats(kFALSE);
	ksSpectra_pSmooth[0]->SetLineColor(kBlack);
	ksSpectra_pSmooth[0]->SetTitle("K^{0}_{s} pSmooth/noSmooth");
	ksSpectra_pSmooth[0]->Draw("P");

	c1->cd(3);

	laSpectra_vtx[0]->Divide( laSpectra_noSmooth[0] );
	laSpectra_vtx[0]->SetMarkerStyle(20);
	laSpectra_vtx[0]->GetYaxis()->SetRangeUser(0.6,1.3);
	laSpectra_vtx[0]->SetStats(kFALSE);
	laSpectra_vtx[0]->SetLineColor(kBlack);
	laSpectra_vtx[0]->SetTitle("#Lambda/#bar{#Lambda} vertex/noSmooth");
	laSpectra_vtx[0]->Draw("P");

	c1->cd(4);

	laSpectra_pSmooth[0]->Divide( laSpectra_noSmooth[0] );
	laSpectra_pSmooth[0]->SetMarkerStyle(20);
	laSpectra_pSmooth[0]->GetYaxis()->SetRangeUser(0.6,1.3);
	laSpectra_pSmooth[0]->SetStats(kFALSE);
	laSpectra_pSmooth[0]->SetLineColor(kBlack);
	laSpectra_pSmooth[0]->SetTitle("#Lambda/#bar{#Lambda} pSmooth/noSmooth");
	laSpectra_pSmooth[0]->Draw("P");






}


