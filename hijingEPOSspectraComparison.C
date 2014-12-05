#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEPOSspectraComparison(){

	gStyle->SetErrorX(0);

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEPOSspectra/HIJINGspectra_5rpyBins.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEPOSspectra/EPOSspectra_5rpyBins.root");

	TH1D* ksSpectra_reco[5];
	TH1D* laSpectra_reco[5];

	TH1D* ksSpectra_gen[5];
	TH1D* laSpectra_gen[5];

	stringstream ksName;
	stringstream laName;

	for(int rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_reco_";
		ksName << rpy+1;

		laName << "laSpectra_reco_";
		laName << rpy+1;

		ksSpectra_reco[rpy] = (TH1D*)file1->Get(ksName.str().c_str());
		laSpectra_reco[rpy] = (TH1D*)file1->Get(laName.str().c_str());

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_gen_";
		ksName << rpy+1;

		laName << "laSpectra_gen_";
		laName << rpy+1;

		ksSpectra_gen[rpy] = (TH1D*)file1->Get(ksName.str().c_str());
		laSpectra_gen[rpy] = (TH1D*)file1->Get(laName.str().c_str());

	}


	TH1D* ksSpectra_epos_reco[5];
	TH1D* laSpectra_epos_reco[5];

	TH1D* ksSpectra_epos_gen[5];
	TH1D* laSpectra_epos_gen[5];

	stringstream ksName;
	stringstream laName;

	for(int rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_reco_";
		ksName << rpy+1;

		laName << "laSpectra_epos_reco_";
		laName << rpy+1;

		ksSpectra_epos_reco[rpy] = (TH1D*)file2->Get(ksName.str().c_str());
		laSpectra_epos_reco[rpy] = (TH1D*)file2->Get(laName.str().c_str());

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_gen_";
		ksName << rpy+1;

		laName << "laSpectra_epos_gen_";
		laName << rpy+1;

		ksSpectra_epos_gen[rpy] = (TH1D*)file2->Get(ksName.str().c_str());
		laSpectra_epos_gen[rpy] = (TH1D*)file2->Get(laName.str().c_str());

	}

	TLatex* ratio[5];
	ratio[0] = new TLatex(6.5,1000,"-2.87 < y < -1.8");
	ratio[1] = new TLatex(6.5,1000,"-1.8 < y < -0.9");
	ratio[2] = new TLatex(6.5,1000,"-0.9 < y < 0");
	ratio[3] = new TLatex(6.5,1000,"0 < y < 0.93");
	ratio[4] = new TLatex(6.5,1000,"0.93 < y < 1.93");

	TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->SetLineColor(kWhite);
    w1->SetFillColor(0);
    
    w1->AddEntry(ksSpectra_reco[0],"HIJING");
    w1->AddEntry(ksSpectra_epos_reco[0],"EPOS");

	TCanvas* c1[5];

	for(rpy = 0; rpy < 5; rpy++){

		c1[rpy] = new TCanvas();

		c1[rpy]->Divide(2,2,0,0);

		c1[rpy]->cd(1);
		gPad->SetLogy(1);

		ksSpectra_reco[rpy]->Draw("P");
		ksSpectra_reco[rpy]->SetMarkerStyle(20);
		ksSpectra_reco[rpy]->SetMarkerSize(1.3);
		ksSpectra_reco[rpy]->SetMarkerColor(kBlue);
		ksSpectra_reco[rpy]->SetStats(kFALSE);
		ksSpectra_reco[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
		ksSpectra_reco[rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
		ksSpectra_reco[rpy]->SetTitle("K^{0}_{s} RECO");
		ksSpectra_reco[rpy]->GetYaxis()->SetTitleSize(0.06);
		ksSpectra_reco[rpy]->GetYaxis()->SetTitleOffset(0.7);
		ksSpectra_reco[rpy]->GetYaxis()->SetRangeUser(1,100000000);
		ksSpectra_reco[rpy]->GetXaxis()->SetRangeUser(0,9.0);

		ksSpectra_epos_reco[rpy]->SetMarkerStyle(20);
		ksSpectra_epos_reco[rpy]->SetMarkerSize(1.3);
		ksSpectra_epos_reco[rpy]->SetMarkerColor(kRed);
		ksSpectra_epos_reco[rpy]->Scale(3);
		ksSpectra_epos_reco[rpy]->Draw("same");

		c1[rpy]->cd(2);
		gPad->SetLogy(1);

		ksSpectra_gen[rpy]->Draw("P");
		ksSpectra_gen[rpy]->SetMarkerStyle(20);
		ksSpectra_gen[rpy]->SetMarkerSize(1.3);
		ksSpectra_gen[rpy]->SetMarkerColor(kBlue);
		ksSpectra_gen[rpy]->SetStats(kFALSE);
		ksSpectra_gen[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
		ksSpectra_gen[rpy]->SetYTitle("Normalized Yields");
		ksSpectra_gen[rpy]->SetTitle("K^{0}_{s} GEN");

		ksSpectra_gen[rpy]->GetYaxis()->SetRangeUser(1,100000000);
		ksSpectra_gen[rpy]->GetXaxis()->SetRangeUser(0,9.0);
		
		ksSpectra_epos_gen[rpy]->SetMarkerStyle(20);
		ksSpectra_epos_gen[rpy]->SetMarkerSize(1.3);
		ksSpectra_epos_gen[rpy]->SetMarkerColor(kRed);
		ksSpectra_epos_gen[rpy]->Scale(3);
		ksSpectra_epos_gen[rpy]->Draw("same");

		w1->Draw("same");

		c1[rpy]->cd(3);
		gPad->SetLogy(1);

		laSpectra_reco[rpy]->Draw("P");
		laSpectra_reco[rpy]->SetMarkerStyle(20);
		laSpectra_reco[rpy]->SetMarkerSize(1.3);
		laSpectra_reco[rpy]->SetMarkerColor(kBlue);
		laSpectra_reco[rpy]->SetStats(kFALSE);
		laSpectra_reco[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
		laSpectra_reco[rpy]->SetTitle("#Lambda/#bar{#Lambda} RECO");
		
		laSpectra_reco[rpy]->GetXaxis()->SetTitleSize(0.06);
		laSpectra_reco[rpy]->GetXaxis()->SetTitleOffset(0.7);

		laSpectra_reco[rpy]->GetYaxis()->SetRangeUser(1,100000000);
		laSpectra_reco[rpy]->GetXaxis()->SetRangeUser(0,9.0);

		laSpectra_epos_reco[rpy]->SetMarkerStyle(20);
		laSpectra_epos_reco[rpy]->SetMarkerSize(1.3);
		laSpectra_epos_reco[rpy]->SetMarkerColor(kRed);
		laSpectra_epos_reco[rpy]->Scale(3);
		laSpectra_epos_reco[rpy]->Draw("same");

		c1[rpy]->cd(4);
		gPad->SetLogy(1);

		laSpectra_gen[rpy]->Draw("P");
		laSpectra_gen[rpy]->SetMarkerStyle(20);
		laSpectra_gen[rpy]->SetMarkerSize(1.3);
		laSpectra_gen[rpy]->SetMarkerColor(kBlue);
		laSpectra_gen[rpy]->SetStats(kFALSE);
		laSpectra_gen[rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
		laSpectra_gen[rpy]->SetYTitle("Normalized Yields");
		laSpectra_gen[rpy]->SetTitle("#Lambda/#bar{#Lambda} GEN");

		laSpectra_gen[rpy]->GetXaxis()->SetTitleSize(0.06);
		laSpectra_gen[rpy]->GetXaxis()->SetTitleOffset(0.7);

		laSpectra_gen[rpy]->GetYaxis()->SetRangeUser(1,100000000);
		laSpectra_gen[rpy]->GetXaxis()->SetRangeUser(0,9.0);
		
		laSpectra_epos_gen[rpy]->SetMarkerStyle(20);
		laSpectra_epos_gen[rpy]->SetMarkerSize(1.3);
		laSpectra_epos_gen[rpy]->SetMarkerColor(kRed);
		laSpectra_epos_gen[rpy]->Scale(3);
		laSpectra_epos_gen[rpy]->Draw("same");

		ratio[rpy]->Draw("same");

	}	



}