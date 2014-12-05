#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingMultDifference(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingMultEfficiencyProducer/HIJING_4multBins_allEta_pPb_PbPb_v14.root");
	
	TH1D* ks_eff[8];
	TH1D* la_eff[8];

	stringstream ksName;
	stringstream laName;

	for(int mult = 0; mult < 8; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ks_eff_mult_";
		ksName << mult+1;

		laName << "la_eff_mult_";
		laName << mult+1;

		ks_eff[mult] = (TH1D*)file1->Get( ksName.str().c_str() );
		la_eff[mult] = (TH1D*)file1->Get( laName.str().c_str() );

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
		ks_eff[mult]->Divide( ks_eff[mult+4] );
		ks_eff[mult]->SetMarkerStyle(20);
		ks_eff[mult]->SetLineColor(kBlack);
		ks_eff[mult]->SetStats(kFALSE);
		ks_eff[mult]->GetYaxis()->SetRangeUser(0.8,1.7);
		ks_eff[mult]->SetYTitle("Hijing pPb/PbPb eff ratio");
		ks_eff[mult]->SetXTitle("pT(GeV/c)");
		//ks_eff[mult]->SetTitle("#Lambda/#bar{#Lambda}");
		ks_eff[mult]->SetTitle("K^{0}_{s}");
		ks_eff[mult]->Draw("P");
		ratio[mult]->Draw("same");
		l1->Draw("same");

	}





}
