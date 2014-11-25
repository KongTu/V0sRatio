#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingMultDifference(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingMultEfficiencyProducer/HIJING_8multBins_3rdYbins_PbPb_6M_v12.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingMultEfficiencyProducer/HIJING_4multBins_3rdYbins_pPb_10M_v13.root");

	TH1D* ks_eff_pPb[4];
	TH1D* la_eff_pPb[4];

	TH1D* ks_eff_PbPb[4];
	TH1D* la_eff_PbPb[4];

	stringstream ksName;
	stringstream laName;

	for(int mult = 0; mult < 4; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_mult_";
		ksName << mult+1;

		laName << "la_eff_mult_";
		laName << mult+1;

		ks_eff_pPb[mult] = (TH1D*)file2->Get( ksName.str().c_str() );
		la_eff_pPb[mult] = (TH1D*)file2->Get( laName.str().c_str() );

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_mult_PbPb";
		ksName << mult+1;

		laName << "la_eff_mult_PbPb";
		laName << mult+1;
		
		ks_eff_PbPb[mult] = (TH1D*)file1->Get( ksName.str().c_str() );
		la_eff_PbPb[mult] = (TH1D*)file1->Get( laName.str().c_str() );
	
	
	}

	TLine* l1 = new TLine(0,1,9.0,1.0);
	l1->SetLineWidth(2);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(2);

	TLatex* ratio[4];

	ratio[0] = new TLatex(1,1.5,"0 < N^{offline}_{trk} < 35");
	ratio[1] = new TLatex(1,1.5,"35 < N^{offline}_{trk} < 60");
	ratio[2] = new TLatex(1,1.5,"60 < N^{offline}_{trk} < 90");
	ratio[3] = new TLatex(1,1.5,"90 < N^{offline}_{trk} < 120");

	TCanvas* c1 = new TCanvas();
	c1->Divide(2,2,0,0);
	for(mult = 0; mult < 4; mult++){

		c1->cd(mult+1);
		la_eff_pPb[mult]->Divide( la_eff_PbPb[mult] );
		la_eff_pPb[mult]->SetMarkerStyle(20);
		la_eff_pPb[mult]->SetLineColor(kBlack);
		la_eff_pPb[mult]->SetStats(kFALSE);
		la_eff_pPb[mult]->GetYaxis()->SetRangeUser(0.8,1.7);
		la_eff_pPb[mult]->SetYTitle("Hijing pPb/PbPb eff ratio");
		la_eff_pPb[mult]->SetXTitle("pT(GeV/c)");
		la_eff_pPb[mult]->SetTitle("#Lambda/#bar{#Lambda}");
		la_eff_pPb[mult]->Draw("P");
		ratio[mult]->Draw("same");
		l1->Draw("same");

	}





}
