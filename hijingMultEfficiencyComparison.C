#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;


void hijingMultEfficiencyComparison(){

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingMultEfficiencyProducer/HIJING_8multBins_5thYbin_6M_v10.root");

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
		laName << mult+1;;

		ks_eff[mult] = (TH1D*)file->Get(ksName.str().c_str());
		la_eff[mult] = (TH1D*)file->Get(laName.str().c_str());
	
	}

	double pTbinsBound[] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};

    double ks_pTbinsBound[29] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[28] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

	TH1D* ks_eff_new = new TH1D("ks_new","ks_new",28,ks_ptbins);
	TH1D* la_eff_new = new TH1D("la_new","la_new",20,ptbins);

	for (mult = 0; mult < 8; mult++){

		ks_eff_new->Add(ks_eff[mult],1);
		la_eff_new->Add(la_eff[mult],1);

	}

	ks_eff_new->Scale(0.125);
	la_eff_new->Scale(0.125);

	TLatex* ratio[8];

	ratio[0] = new TLatex(1,1.4,"0 < N^{offline}_{trk} < 35");
	ratio[1] = new TLatex(1,1.4,"35 < N^{offline}_{trk} < 60");
	ratio[2] = new TLatex(1,1.4,"60 < N^{offline}_{trk} < 90");
	ratio[3] = new TLatex(1,1.4,"90 < N^{offline}_{trk} < 120");
	ratio[4] = new TLatex(1,1.4,"120 < N^{offline}_{trk} < 150");
	ratio[5] = new TLatex(1,1.4,"150 < N^{offline}_{trk} < 185");
	ratio[6] = new TLatex(1,1.4,"185 < N^{offline}_{trk} < 220");
	ratio[7] = new TLatex(1,1.4,"220 < N^{offline}_{trk} < #infty");

	TLatex* r3 = new TLatex(1.36,1.2,"K^{0}_{s}");
    r3->SetTextSize(0.07);
    TLatex* r4 = new TLatex(1.36,1.2,"#Lambda/#bar{#Lambda}");
    r4->SetTextSize(0.07);

	TLine* l1 = new TLine(0,1,9.0,1.0);
	l1->SetLineWidth(2);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(2);

	TCanvas* c1 = new TCanvas();
	c1->Divide(4,2,0,0);

	for(mult = 0; mult < 8; mult++){

		ks_eff[mult]->GetYaxis()->SetLabelSize(0.08);
		ks_eff[mult]->GetXaxis()->SetLabelSize(0.08);
		ks_eff[mult]->SetYTitle("mult/mean");
		ks_eff[mult]->GetYaxis()->SetTitleSize(0.07);
		ks_eff[mult]->GetYaxis()->SetTitleOffset(0.64);

		ks_eff[mult]->SetXTitle("P^{}_{T}(GeV/c)");
		ks_eff[mult]->GetXaxis()->SetTitleSize(0.07);
		ks_eff[mult]->GetXaxis()->SetTitleOffset(1.0);
		
	}
	

	for(mult = 0; mult < 8; mult++){

		c1->cd(mult+1);
		
		ks_eff[mult]->Divide(ks_eff_new);
		ks_eff[mult]->SetMarkerStyle(20);
		
		ks_eff[mult]->GetYaxis()->SetRangeUser(0.5,1.6);
		ks_eff[mult]->SetStats(kFALSE);
		ks_eff[mult]->SetTitle("");
		ks_eff[mult]->SetLineColor(kBlack);
		ks_eff[mult]->Draw("P");
		ratio[mult]->SetTextSize(0.1);
		ratio[mult]->Draw("same");
		l1->Draw("same");

	}
	
	r3->Draw("same");


}