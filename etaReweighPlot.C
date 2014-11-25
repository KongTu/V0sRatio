#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void etaReweighPlot(){

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingMultEfficiencyProducer/HIJING_2multBins_6M_v10_test.root");

	TH1D* ks_eff_mult_1 = (TH1D*)file->Get("ks_eff_mult_1");
	TH1D* ks_eff_mult_2 = (TH1D*)file->Get("ks_eff_mult_2");

	TH1D* la_eff_mult_1 = (TH1D*)file->Get("la_eff_mult_1");
	TH1D* la_eff_mult_2 = (TH1D*)file->Get("la_eff_mult_2");

	TLine* l1 = new TLine(0,1,9.0,1.0);
	l1->SetLineWidth(2);
	l1->SetLineColor(kRed);
	l1->SetLineStyle(2);

	TCanvas* c1 = new TCanvas();
	c1->Divide(2,1,0,0);
	c1->cd(1);
	ks_eff_mult_1->Divide(ks_eff_mult_2);
	ks_eff_mult_1->SetMarkerStyle(20);
	ks_eff_mult_1->GetYaxis()->SetRangeUser(0.9,1.6);
	ks_eff_mult_1->SetLineColor(kBlack);
	ks_eff_mult_1->SetStats(kFALSE);
	ks_eff_mult_1->Draw("P");
	l1->Draw("same");
	c1->cd(2);
	la_eff_mult_1->Divide(la_eff_mult_2);
	la_eff_mult_1->SetMarkerStyle(20);
	la_eff_mult_1->GetYaxis()->SetRangeUser(0.9,1.6);
	la_eff_mult_1->SetLineColor(kBlack);
	la_eff_mult_1->SetStats(kFALSE);
	la_eff_mult_1->Draw("P");
	l1->Draw("same");

}