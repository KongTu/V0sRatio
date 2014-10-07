#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void EasyBlastWaveFit(){

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy/HMbins_18pTBins5rpyBins_wSmooth10Eff_FullStats_5Multbins_v1.root");

	TH1D* ratioHist[5][5];
    TH1D* ksSpectra[5][5];
    TH1D* laSpectra[5][5];
    TF1* fit_ks[5][5];
    TF1* fit_la[5][5];
    TF1* fit;

    double ks_Tprime[5][5];
    double la_Tprime[5][5];

	stringstream ratioHistName;
    stringstream ksSpectraName;
    stringstream laSpectraName;

    TCanvas* c1 = new TCanvas();
    gPad->SetLogy(1);

    c1->Print("ks_ExpoFit1.pdf[");

	for(int mult = 0; mult < 5; mult++){

		for(int rpy = 0; rpy < 5; rpy++){

			ratioHistName.str("");
			ratioHistName << "ratioHist_";
			ratioHistName << mult+1;
			ratioHistName << "_";
			ratioHistName << rpy+1;

            ksSpectraName.str("");
            ksSpectraName << "ksSpectra_";
            ksSpectraName << mult+1;
            ksSpectraName << "_";
            ksSpectraName << rpy+1;

            laSpectraName.str("");
            laSpectraName << "laSpectra_";
            laSpectraName << mult+1;
            laSpectraName << "_";
            laSpectraName << rpy+1;

			ratioHist[mult][rpy] = (TH1D*)file->Get( ratioHistName.str().c_str() );

	            ksSpectra[mult][rpy] = (TH1D*)file->Get( ksSpectraName.str().c_str() );
	            	ksSpectra[mult][rpy]->Fit("expo","","",0.9,2.5);
	            	ksSpectra[mult][rpy]->Draw();
	            		c1->Print("ks_ExpoFit1.pdf");
	            	
	            	fit_ks[mult][rpy] = ksSpectra[mult][rpy]->GetFunction("expo");
	            	ks_Tprime[mult][rpy] = fit_ks[mult][rpy]->GetParameter(1);

	            laSpectra[mult][rpy] = (TH1D*)file->Get( laSpectraName.str().c_str() );
	            	laSpectra[mult][rpy]->Fit("expo","","",0.4,4.0);

	            	fit_la[mult][rpy] = laSpectra[mult][rpy]->GetFunction("expo");
	            	la_Tprime[mult][rpy] = fit_la[mult][rpy]->GetParameter(1);

		}
	}

	c1->Print("ks_ExpoFit1.pdf]");

	double ptbins[] = {1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};

	TGraph* g[5];

	for(int i = 0; i < 5; i++){

		g[0] = new TGraph(2);
		g[1] = new TGraph(2);
		g[2] = new TGraph(2);
		g[3] = new TGraph(2);
		g[4] = new TGraph(2);

		double temp_ks = -1.0/ks_Tprime[i][3];
		double temp_la = -1.0/la_Tprime[i][3];

	}
	
	g[0]->SetPoint(0, 0.497, -1.0/ks_Tprime[0][0] );
	g[0]->SetPoint(1, 1.115, -1.0/la_Tprime[0][0] );

	g[0]->SetMarkerStyle(28);
	g[0]->SetMarkerSize(1.0);
	g[0]->SetMarkerColor(kRed);
	g[0]->SetLineColor(kRed);
	g[0]->GetXaxis()->SetLimits(0.0,1.2);
	g[0]->GetYaxis()->SetRangeUser(0,0.7);
	g[0]->GetXaxis()->SetTitle("mass[GeV/c^{2}]");
	g[0]->GetYaxis()->SetTitle("T'[GeV/c]");
	//g[0]->SetTitle("0.93 < y < 1.93");
	g[0]->Draw("APSimpleLine");

	g[1]->SetPoint(0, 0.497, -1.0/ks_Tprime[1][0] );
	g[1]->SetPoint(1, 1.115, -1.0/la_Tprime[1][0] );

	g[1]->SetMarkerStyle(28);
	g[1]->SetMarkerSize(1.0);
	g[1]->SetMarkerColor(kBlue);
	g[1]->SetLineColor(kBlue);
	g[1]->GetXaxis()->SetLimits(0.0,1.2);
	g[1]->GetYaxis()->SetRangeUser(0,0.7);
	g[1]->Draw("PsameSimpleLine");

	g[2]->SetPoint(0, 0.497, -1.0/ks_Tprime[2][0] );
	g[2]->SetPoint(1, 1.115, -1.0/la_Tprime[2][0] );

	g[2]->SetMarkerStyle(28);
	g[2]->SetMarkerSize(1.0);
	g[2]->SetMarkerColor(kGreen);
	g[2]->SetLineColor(kGreen);
	g[2]->GetXaxis()->SetLimits(0.0,1.2);
	g[2]->GetYaxis()->SetRangeUser(0,0.7);
	g[2]->Draw("PsameSimpleLine");

	g[3]->SetPoint(0, 0.497, -1.0/ks_Tprime[3][0] );
	g[3]->SetPoint(1, 1.115, -1.0/la_Tprime[3][0] );

	g[3]->SetMarkerStyle(28);
	g[3]->SetMarkerSize(1.0);
	g[3]->SetMarkerColor(kBlack);
	g[3]->GetXaxis()->SetLimits(0.0,1.2);
	g[3]->GetYaxis()->SetRangeUser(0,0.7);
	g[3]->Draw("PsameSimpleLine");

	g[4]->SetPoint(0, 0.497, -1.0/ks_Tprime[4][0] );
	g[4]->SetPoint(1, 1.115, -1.0/la_Tprime[4][0] );

	g[4]->SetMarkerStyle(28);
	g[4]->SetMarkerSize(1.0);
	g[4]->SetMarkerColor(kYellow);
	g[4]->SetLineColor(kYellow);
	g[4]->GetXaxis()->SetLimits(0.0,1.2);
	g[4]->GetYaxis()->SetRangeUser(0,0.7);
	g[4]->Draw("PsameSimpleLine");
	
	//g[1]->SetMarkerColor(kBlue);
	//g[1]->Draw("Psame");
	//g[2]->Draw("Psame");
	//g[3]->Draw("Psame");
	//g[4]->Draw("Psame");
	









	
/*
	RooMsgService::instance().setGlobalKillBelow(RooFit::ERROR);

	RooRealVar x("x","mT",0,8);
	RooDataHist data("data","dataset", x, input );
	RooPlot* xframe = x.frame();

	data.plotOn( xframe );

	RooRealVar Tprime("Tprime","Tprime",0,-100000,100000);
	//RooRealVar a("a","a",0,-1000,1000);
	//RooRealVar b("b","b",0,-100000,100000);

	//RooPolynomial poly("poly","poly",x,RooArgList(a));
	RooExponential expo("expo","expo",x,Tprime);
	//RooAddPdf sum("sum","sum",RooArgList(expo,poly),RooArgList(Tprime,a));

	x.setRange("cut",-1,3.0);

	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));
	expo.fitTo(data,Range("cut"));


	expo.plotOn(xframe,Name("expo"),NormRange("cut"),LineWidth(0.5),LineColor(kRed));

	xframe->Draw();
*/



}