#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEposEfficiencyDifference(){

	gStyle->SetErrorX(0);

	TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyTable/effKongNew2DTable_15M_Sep16_v1_12pTbins.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyTable/EPOS_effKongNew2DTable_15M_Sep16_v1_12pTbins.root");

	TH1D* ks_hijing[6];
	TH1D* ks_epos[6];

	TH1D* la_hijing[6];
	TH1D* la_epos[6];

	stringstream ks_hijingName;
	stringstream ks_eposName;
	stringstream la_hijingName;
	stringstream la_eposName;

	for (int eta = 0; eta < 6; eta++){

		ks_hijingName.str("");
		ks_eposName.str("");
		la_hijingName.str("");
		la_eposName.str("");

		ks_hijingName << "ks_eff_Eta_";
		ks_hijingName << eta + 1;

		ks_eposName << "ks_eff_Eta_";
		ks_eposName << eta + 1;

		la_hijingName << "la_eff_Eta_";
		la_hijingName << eta + 1;

		la_eposName << "la_eff_Eta_";
		la_eposName << eta + 1;

			ks_hijing[eta] = (TH1D*)file1->Get( ks_hijingName.str().c_str() );
			ks_epos[eta] = (TH1D*)file2->Get( ks_eposName.str().c_str() );

			la_hijing[eta] = (TH1D*)file1->Get( la_hijingName.str().c_str() );
			la_epos[eta] = (TH1D*) file2->Get( la_eposName.str().c_str() );

	}


	





/**
 * The following part is to calculate the efficiency in different eta bins:
 */


	TLatex* etarange[6];

	etarange[0] = new TLatex(1,0.15,"-2.4 < #eta < -1.6");
	etarange[1] = new TLatex(1,0.15,"-1.6 < #eta < -0.8");
	etarange[2] = new TLatex(1,0.15,"-0.8 < #eta < 0");
	etarange[3] = new TLatex(1,0.15,"0 < #eta < 0.8");
	etarange[4] = new TLatex(1,0.15,"0.8 < #eta < 1.6");
	etarange[5] = new TLatex(1,0.15,"1.6 < #eta < 2.4");


	TLegend *w1 = new TLegend(0.25,0.37,0.5,0.5);
    w1->SetLineColor(kWhite);

    w1->AddEntry(ks_hijing[0],"HIJING");
    w1->AddEntry(ks_epos[0], "EPOS");

	TCanvas* wq1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    wq1->Divide(2,3,0,0);

    for (eta = 0; eta < 6; eta++){

    	wq1->cd(eta+1);
    	gPad->SetTicks();

    	ks_hijing[eta]->SetStats(kFALSE);
        ks_hijing[eta]->GetXaxis()->SetTitleOffset(1.0);
        ks_hijing[eta]->GetYaxis()->SetTitleOffset(1.0);
        ks_hijing[eta]->GetXaxis()->SetTitleSize(0.05);
        ks_hijing[eta]->GetYaxis()->SetTitleSize(0.05);
        ks_hijing[eta]->GetXaxis()->SetLabelSize(0.06);
        ks_hijing[eta]->GetYaxis()->SetLabelSize(0.06);

        ks_hijing[eta]->SetYTitle("K^{0}_{s} Efficiency");
        ks_hijing[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        ks_hijing[eta]->SetMarkerStyle(28);
        ks_hijing[eta]->SetMarkerSize(1.0);
        ks_hijing[eta]->SetMarkerColor(kBlue);
        ks_hijing[eta]->GetYaxis()->SetRangeUser(0,0.3); 

        ks_epos[eta]->SetStats(kFALSE);
        ks_epos[eta]->GetXaxis()->SetTitleOffset(1.0);
        ks_epos[eta]->GetYaxis()->SetTitleOffset(1.0);
        ks_epos[eta]->GetXaxis()->SetTitleSize(0.05);
        ks_epos[eta]->GetYaxis()->SetTitleSize(0.05);
        ks_epos[eta]->GetXaxis()->SetLabelSize(0.06);
        ks_epos[eta]->GetYaxis()->SetLabelSize(0.06);

        ks_epos[eta]->SetYTitle("K^{0}_{s} Efficiency");
        ks_epos[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        ks_epos[eta]->SetMarkerStyle(24);
        ks_epos[eta]->SetMarkerSize(1.0);
        ks_epos[eta]->SetMarkerColor(kRed);
        ks_epos[eta]->GetYaxis()->SetRangeUser(0,0.3);   

        ks_hijing[eta]->Draw("P");
        ks_epos[eta]->Draw("Psame");
        etarange[eta]->Draw("same");

        

    }

    w1->Draw("same");


 	TCanvas* wq2 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    wq2->Divide(2,3,0,0);

    for (eta = 0; eta < 6; eta++){

    	wq2->cd(eta+1);
    	gPad->SetTicks();

    	la_hijing[eta]->SetStats(kFALSE);
        la_hijing[eta]->GetXaxis()->SetTitleOffset(1.0);
        la_hijing[eta]->GetYaxis()->SetTitleOffset(1.0);
        la_hijing[eta]->GetXaxis()->SetTitleSize(0.05);
        la_hijing[eta]->GetYaxis()->SetTitleSize(0.05);
        la_hijing[eta]->GetXaxis()->SetLabelSize(0.06);
        la_hijing[eta]->GetYaxis()->SetLabelSize(0.06);

        la_hijing[eta]->SetYTitle("#Lambda/#bar{#Lambda} Efficiency");
        la_hijing[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        la_hijing[eta]->SetMarkerStyle(28);
        la_hijing[eta]->SetMarkerSize(1.0);
        la_hijing[eta]->SetMarkerColor(kBlue);
        la_hijing[eta]->GetYaxis()->SetRangeUser(0,0.2);    

        la_epos[eta]->SetStats(kFALSE);
        la_epos[eta]->GetXaxis()->SetTitleOffset(1.0);
        la_epos[eta]->GetYaxis()->SetTitleOffset(1.0);
        la_epos[eta]->GetXaxis()->SetTitleSize(0.05);
        la_epos[eta]->GetYaxis()->SetTitleSize(0.05);
        la_epos[eta]->GetXaxis()->SetLabelSize(0.06);
        la_epos[eta]->GetYaxis()->SetLabelSize(0.06);

        la_epos[eta]->SetYTitle("#Lambda/#bar{#Lambda} Efficiency");
        la_epos[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        la_epos[eta]->SetMarkerStyle(24);
        la_epos[eta]->SetMarkerSize(1.0);
        la_epos[eta]->SetMarkerColor(kRed);
        la_epos[eta]->GetYaxis()->SetRangeUser(0,0.2); 

        la_hijing[eta]->Draw("P");
        la_epos[eta]->Draw("Psame");
        
        etarange[eta]->Draw("same");
       

    }

     w1->Draw("same");


 	TCanvas* wq3 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    wq3->Divide(2,3,0,0);



    for (eta = 0; eta < 6; eta++){

    	wq3->cd(eta+1);
    	gPad->SetTicks();

    	la_hijing[eta]->Divide( la_epos[eta] );

    	la_hijing[eta]->SetStats(kFALSE);
        la_hijing[eta]->GetXaxis()->SetTitleOffset(1.0);
        la_hijing[eta]->GetYaxis()->SetTitleOffset(1.0);
        la_hijing[eta]->GetXaxis()->SetTitleSize(0.05);
        la_hijing[eta]->GetYaxis()->SetTitleSize(0.05);
        la_hijing[eta]->GetXaxis()->SetLabelSize(0.06);
        la_hijing[eta]->GetYaxis()->SetLabelSize(0.06);

        la_hijing[eta]->SetYTitle("HIJING/EPOS Efficiency");
        la_hijing[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        la_hijing[eta]->SetMarkerStyle(28);
        la_hijing[eta]->SetMarkerSize(1.0);
        la_hijing[eta]->SetMarkerColor(kBlue);
        la_hijing[eta]->GetYaxis()->SetRangeUser(0,2);    

        la_hijing[eta]->Draw("P");   
        etarange[eta]->Draw("same");
       

    }






}

