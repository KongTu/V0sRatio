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

double errorCal_lambdakshort(double yield_la, double yield_ks, double la_err, double ks_err){

    double first = (la_err*la_err)/(4*(yield_ks)*(yield_ks));
    double second = ( (yield_la)*(yield_la)*ks_err*ks_err )/(16*(yield_ks)*(yield_ks)*(yield_ks)*(yield_ks));

    double error = sqrt( first + second );
    return error;

}

void EfficiencySmoothing(){

	gStyle->SetErrorX(0);

	double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_ptbincenter[28] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
  	double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
	//TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/eposEfficiencyRapidityTable/EPOSeffKongNew2DTable_5M_Nov3_rapidity_v1_28ks_pTbins.root");

	TH1D* ks_eff_rpy[5];
	TH1D* la_eff_rpy[5];

	TH1D* ks_eff_rpy_smooth[5];
	TH1D* la_eff_rpy_smooth[5];

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

		ks_eff_rpy_smooth[rpy] = (TH1D*) file->Get( ksName.str().c_str() );
		la_eff_rpy_smooth[rpy] = (TH1D*) file->Get( laName.str().c_str() );

	}

    double ks_no_smooth[5][28];
    double la_no_smooth[5][20];

    double ks_no_smooth_err[5][28];
    double la_no_smooth_err[5][20];

    double ks_smooth[5][28];
    double la_smooth[5][20];

    double ks_smooth_err[5][28];
    double la_smooth_err[5][20];

    double ks_ratio[5][28];
    double la_ratio[5][20];

    double ks_ratio_err[5][28];
    double la_ratio_err[5][20];

    for(rpy = 0; rpy < 5; rpy++){

    	for(int pt = 0; pt < 28; pt++){

    		ks_no_smooth[rpy][pt] = ks_eff_rpy[rpy]->GetBinContent(pt+1);
    			ks_no_smooth_err[rpy][pt] = ks_eff_rpy[rpy]->GetBinError(pt+1);

    	}

    	for(pt = 0; pt < 20; pt++){

    		la_no_smooth[rpy][pt] = la_eff_rpy[rpy]->GetBinContent(pt+1);
    			la_no_smooth_err[rpy][pt] = la_eff_rpy[rpy]->GetBinError(pt+1);
    	}

    	ks_eff_rpy_smooth[rpy]->Smooth(10);
    	la_eff_rpy_smooth[rpy]->Smooth(10);

    	for(pt = 0; pt < 28; pt++){

    		ks_smooth[rpy][pt] = ks_eff_rpy_smooth[rpy]->GetBinContent(pt+1);
    			ks_smooth_err[rpy][pt] = ks_eff_rpy_smooth[rpy]->GetBinError(pt+1);
    	}

    	for(pt = 0; pt < 20; pt++){

    		la_smooth[rpy][pt] = la_eff_rpy_smooth[rpy]->GetBinContent(pt+1);
    			la_smooth_err[rpy][pt] = la_eff_rpy_smooth[rpy]->GetBinError(pt+1);
    	}

    /*
    CAL ratio between smooth/no_smooth:
     */
    	for(pt = 0; pt < 28; pt++){

    		ks_ratio[rpy][pt] = ks_smooth[rpy][pt]/ks_no_smooth[rpy][pt];
    			ks_ratio_err[rpy][pt] = errorCal_lambdakshort(ks_smooth[rpy][pt],ks_no_smooth[rpy][pt],ks_smooth_err[rpy][pt],ks_no_smooth_err[rpy][pt]);

    		 
    	}

    	for(pt = 0; pt < 20; pt++){

    		la_ratio[rpy][pt] = la_smooth[rpy][pt]/la_no_smooth[rpy][pt];
    			la_ratio_err[rpy][pt] = errorCal_lambdakshort(la_smooth[rpy][pt],la_no_smooth[rpy][pt],la_smooth_err[rpy][pt],la_no_smooth_err[rpy][pt]);
    	}

    }

    TGraphErrors* ks_g[5];
    TGraphErrors* la_g[5];

    for(rpy = 0; rpy < 5; rpy++){

    	ks_g[rpy] = new TGraphErrors(28);
    	la_g[rpy] = new TGraphErrors(20);

    	for(pt = 0; pt < 28; pt++){

    		ks_g[rpy]->SetPoint(pt,ks_ptbincenter[pt],ks_ratio[rpy][pt]);
    			ks_g[rpy]->SetPointError(pt,0,ks_ratio_err[rpy][pt]);
    		 
    	}

    	for(pt = 0; pt < 20; pt++){

    		la_g[rpy]->SetPoint(pt,la_ptbincenter[pt],la_ratio[rpy][pt]);
    			la_g[rpy]->SetPointError(pt,0,la_ratio_err[rpy][pt]);

    	}

    }

    TCanvas* ks = new TCanvas();
    ks->Divide(2,3,0,0);

    for(rpy = 0; rpy < 5; rpy++){

    	ks->cd(rpy+1);
    	ks_g[rpy]->SetTitle("K^{0}_{s}");

    	ks_g[rpy]->SetMarkerStyle(20);
    	ks_g[rpy]->SetMarkerSize(1.3);
    	ks_g[rpy]->GetYaxis()->SetRangeUser(0.65,1.6);

    	ks_g[rpy]->Draw("AP");
    }

    return;

    TLatex* r6[5];
    r6[0] = new TLatex(5.5,0.3,"-2.87 < y < -1.8");
    r6[1] = new TLatex(5.5,0.3,"-1.8 < y < -0.9");
    r6[2] = new TLatex(5.5,0.3,"-0.9 < y < 0");
    r6[3] = new TLatex(5.5,0.3,"0 < y < 0.93");
    r6[4] = new TLatex(5.5,0.3,"0.93 < y < 1.93");

    TLatex* r1 = new TLatex(1.0,0.3,"beforeSmoothing");
    TLatex* r4 = new TLatex(1.0,0.3,"afterSmoothing");
    TLatex* r2 = new TLatex(4.0,0.3,"K^{0}_{s}");
    TLatex* r3 = new TLatex(4.0,0.3,"#Lambda/#bar{#Lambda}");


	TCanvas* c1 = new TCanvas();
	c1->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c1->cd(rpy+1);

		ks_eff_rpy[rpy]->GetYaxis()->SetRangeUser(0,2);
		ks_eff_rpy[rpy]->SetStats(kFALSE);
		ks_eff_rpy[rpy]->SetMarkerStyle(20);
		ks_eff_rpy[rpy]->SetMarkerSize(1.1);
		ks_eff_rpy[rpy]->SetMarkerColor(kBlue);
		ks_eff_rpy[rpy]->SetXTitle("pT(GeV/c)");
		ks_eff_rpy[rpy]->SetYTitle("Efficiency");
		ks_eff_rpy[rpy]->SetTitle("");


		ks_eff_rpy[rpy]->Draw("P");
		
		r6[rpy]->Draw("same");
		r1->Draw("same");
		r2->Draw("same");
	}



	//c1->SaveAs("./beforeSmoothing_EPOS_ks_all5rpyBins.pdf");

	TCanvas* c2 = new TCanvas();
	c2->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c2->cd(rpy+1);

		la_eff_rpy[rpy]->GetYaxis()->SetRangeUser(0,5);
		la_eff_rpy[rpy]->SetStats(kFALSE);
		la_eff_rpy[rpy]->SetMarkerStyle(20);
		la_eff_rpy[rpy]->SetMarkerSize(1.1);
		la_eff_rpy[rpy]->SetMarkerColor(kRed);
		la_eff_rpy[rpy]->SetXTitle("pT(GeV/c)");
		la_eff_rpy[rpy]->SetYTitle("Efficiency");
		la_eff_rpy[rpy]->SetTitle("");

		la_eff_rpy[rpy]->Draw("P");
		r6[rpy]->Draw("same");
		r1->Draw("same");
		r3->Draw("same");
	}

	//c2->SaveAs("./beforeSmoothing_EPOS_la_all5rpyBins.pdf");

	TCanvas* c3 = new TCanvas();
	c3->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c3->cd(rpy+1);	

		ks_eff_rpy[rpy]->GetYaxis()->SetRangeUser(0,0.35);
		ks_eff_rpy[rpy]->SetStats(kFALSE);
		ks_eff_rpy[rpy]->SetMarkerStyle(20);
		ks_eff_rpy[rpy]->SetMarkerSize(1.1);
		ks_eff_rpy[rpy]->SetMarkerColor(kBlue);
		ks_eff_rpy[rpy]->SetXTitle("pT(GeV/c)");
		ks_eff_rpy[rpy]->SetYTitle("Efficiency");
		ks_eff_rpy[rpy]->SetTitle("");
		ks_eff_rpy[rpy]->Draw("P");
		r6[rpy]->Draw("same");
		r4->Draw("same");
		r2->Draw("same");
	}

	//c3->SaveAs("./afterSmoothing_EPOS_ks_all5rpyBins.pdf");

	TCanvas* c4 = new TCanvas();
	c4->Divide(2,3,0,0);

	for(rpy = 0; rpy < 5; rpy++){

		c4->cd(rpy+1);

		la_eff_rpy[rpy]->GetYaxis()->SetRangeUser(0,0.35);
		la_eff_rpy[rpy]->SetStats(kFALSE);
		la_eff_rpy[rpy]->SetMarkerStyle(20);
		la_eff_rpy[rpy]->SetMarkerSize(1.1);
		la_eff_rpy[rpy]->SetMarkerColor(kRed);
		la_eff_rpy[rpy]->SetXTitle("pT(GeV/c)");
		la_eff_rpy[rpy]->SetYTitle("Efficiency");
		la_eff_rpy[rpy]->SetTitle("");
		la_eff_rpy[rpy]->Draw("P");
		r6[rpy]->Draw("same");
		r4->Draw("same");
		r3->Draw("same");

	}

	//c4->SaveAs("./afterSmoothing_EPOS_la_all5rpyBins.pdf");






}