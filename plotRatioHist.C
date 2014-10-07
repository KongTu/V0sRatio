#include "fitting.h"

using namespace std;

void plotRatioHist(){

    gStyle->SetErrorX(0);

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy/HMbins_18pTBins5rpyBins_wSmooth10Eff_FullStats_5Multbins_v1.root");

	TH1D* ratioHist[5][5];
    TH1D* ksSpectra[5][5];
    TH1D* laSpectra[5][5];

	stringstream ratioHistName;
    stringstream ksSpectraName;
    stringstream laSpectraName;

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
            laSpectra[mult][rpy] = (TH1D*)file->Get( laSpectraName.str().c_str() );


		}
	}

   	TLatex* ratio[5];
    ratio[0] = new TLatex(1.65,1.3,"-2.87 < y < -1.8");
    ratio[1] = new TLatex(1.65,1.3,"-1.8 < y < -0.9");
    ratio[2] = new TLatex(1.65,1.3,"-0.9 < y < 0");
    ratio[3] = new TLatex(1.65,1.3,"0 < y < 0.93");
    ratio[4] = new TLatex(1.65,1.3,"0.93 < y < 1.93");

	TCanvas* c1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    c1->Divide(3,2,0,0);

    for(int it = 0; it < 5; it++){

    	c1->cd(it+1);
    	gPad->SetTicks();

    	ratioHist[4][it]->SetStats(kFALSE);
        ratioHist[4][it]->GetXaxis()->SetTitleOffset(1.0);
        ratioHist[4][it]->GetYaxis()->SetTitleOffset(1.0);
        ratioHist[4][it]->GetXaxis()->SetTitleSize(0.05);
        ratioHist[4][it]->GetYaxis()->SetTitleSize(0.05);
        ratioHist[4][it]->GetXaxis()->SetLabelSize(0.06);
        ratioHist[4][it]->GetYaxis()->SetLabelSize(0.06);

        ratioHist[4][it]->SetYTitle("#Lambda+#bar{#Lambda}/2K^{0}_{s}");
        ratioHist[4][it]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        ratioHist[4][it]->SetMarkerStyle(24);
        ratioHist[4][it]->SetMarkerSize(1.2);
        ratioHist[4][it]->SetMarkerColor(kBlue);
        ratioHist[4][it]->GetYaxis()->SetRangeUser(0,1.5);
    
        ratioHist[4][it]->Draw("P");
        ratio[it]->SetTextSize(0.06);
        ratio[it]->Draw("same");
   
    }

    TLatex* r1 = new TLatex(5.5,52,"CMS,p-Pb 2013");
    r1->SetTextSize(0.04);
    TLatex* r2 = new TLatex(5.5,15,"-2.87 < y < -1.9");
    r2->SetTextSize(0.04);
    TLatex* r3 = new TLatex(0.8,0.000001,"K^{0}_{s}");
    r3->SetTextSize(0.07);
    TLatex* r4 = new TLatex(0.8,0.000001,"#Lambda+#bar{#Lambda}");
    r4->SetTextSize(0.07);
    TLatex* r5 = new TLatex(1.65,0.14,"#Lambda/K^{0}_{s}");
    r5->SetTextSize(0.07);

    TLatex* r6[5];
    r6[0] = new TLatex(5.5,15,"-2.87 < y < -1.8");
    r6[1] = new TLatex(5.5,15,"-1.8 < y < -0.9");
    r6[2] = new TLatex(5.5,15,"-0.9 < y < 0");
    r6[3] = new TLatex(5.5,15,"0 < y < 0.93");
    r6[4] = new TLatex(5.5,15,"0.93 < y < 1.93");

    TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->SetLineColor(kWhite);
    //w1->AddEntry(ksSpectra[0],"0  < N^{offline}_{trk} < 35");
    //w1->AddEntry(ksSpectra[1],"35  < N^{offline}_{trk} < 60");
    //w1->AddEntry(ksSpectra[2],"60  < N^{offline}_{trk} < 90");
    //w1->AddEntry(ksSpectra[3],"90  < N^{offline}_{trk} < 120");
    //
    w1->AddEntry(laSpectra[0][0],"120 < N^{offline}_{trk} < 150");
    w1->AddEntry(laSpectra[1][0],"150 < N^{offline}_{trk} < 185");
    w1->AddEntry(laSpectra[2][0],"185 < N^{offline}_{trk} < 220");
    w1->AddEntry(laSpectra[3][0],"220 < N^{offline}_{trk} < 260");
    w1->AddEntry(laSpectra[4][0],"     N^{offline}_{trk} > 260");

    TCanvas* p1 = new TCanvas();
    laSpectra[0][rpy]->SetMarkerColor(kYellow-3);
    laSpectra[0][rpy]->SetLineColor(kYellow-3);
    laSpectra[0][rpy]->SetMarkerStyle(22);
    laSpectra[0][rpy]->SetXTitle("m^{}_{T} - m^{}_{0} (GeV/c^{2})");
    //laSpectra[0][rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
    laSpectra[0][rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
    //laSpectra[0][rpy]->SetAxisRange(0,8,"X");
    laSpectra[0][rpy]->SetStats(kFALSE);
    laSpectra[0][rpy]->Scale(16);

    laSpectra[0][rpy]->GetYaxis()->SetRangeUser(0.0000001,1000);
    laSpectra[0][rpy]->GetXaxis()->SetRangeUser(0,9);

    laSpectra[0][rpy]->Draw("P");


    TCanvas* y1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);
    y1->SetLogy(1);

    y1->Divide(3,2,0,0);

    for (rpy = 0; rpy < 5; rpy++){

        y1->cd(rpy+1);
        gPad->SetLogy(1);
        gPad->SetTicks();

        laSpectra[0][rpy]->SetMarkerColor(kYellow-3);
        laSpectra[0][rpy]->SetLineColor(kYellow-3);
        laSpectra[0][rpy]->SetMarkerStyle(22);
        laSpectra[0][rpy]->SetXTitle("m^{}_{T} - m^{}_{0} (GeV/c^{2})");
        //laSpectra[0][rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        laSpectra[0][rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
        //laSpectra[0][rpy]->SetAxisRange(0,8,"X");
        laSpectra[0][rpy]->SetStats(kFALSE);
        laSpectra[0][rpy]->Scale(16);

        laSpectra[0][rpy]->GetYaxis()->SetRangeUser(0.0000001,1000);
        laSpectra[0][rpy]->GetXaxis()->SetRangeUser(0,9);

        laSpectra[1][rpy]->SetMarkerColor(6);
        laSpectra[1][rpy]->SetMarkerStyle(22);
        laSpectra[1][rpy]->SetLineColor(6);
        laSpectra[1][rpy]->Scale(32);

        laSpectra[2][rpy]->SetMarkerColor(7);
        laSpectra[2][rpy]->SetLineColor(7);
        laSpectra[2][rpy]->SetMarkerStyle(22);
        laSpectra[2][rpy]->Scale(64);

        laSpectra[3][rpy]->SetMarkerColor(8);
        laSpectra[3][rpy]->SetLineColor(8);
        laSpectra[3][rpy]->SetMarkerStyle(22);
        laSpectra[3][rpy]->Scale(128);

        laSpectra[4][rpy]->SetMarkerColor(9);
        laSpectra[4][rpy]->SetLineColor(9);
        laSpectra[4][rpy]->SetMarkerStyle(22);
        laSpectra[4][rpy]->Scale(256);

        laSpectra[0][rpy]->Draw("P");
        laSpectra[1][rpy]->Draw("Psame");
        laSpectra[2][rpy]->Draw("Psame");
        laSpectra[3][rpy]->Draw("Psame");
        laSpectra[4][rpy]->Draw("Psame");
        r6[rpy]->Draw("same");

    }

/*
    ksSpectra[0]->SetMarkerColor(1);
    ksSpectra[0]->SetMarkerStyle(22);
    ksSpectra[0]->SetLineColor(1);
    ksSpectra[0]->Scale(1);

    ksSpectra[1]->SetMarkerColor(2);
    ksSpectra[1]->SetLineColor(2);
    ksSpectra[1]->SetMarkerStyle(22);
    ksSpectra[1]->Scale(2);

    ksSpectra[2]->SetMarkerColor(3);
    ksSpectra[2]->SetLineColor(3);
    ksSpectra[2]->SetMarkerStyle(22);
    ksSpectra[2]->Scale(4);

    ksSpectra[3]->SetMarkerColor(4);
    ksSpectra[3]->SetLineColor(4);
    ksSpectra[3]->SetMarkerStyle(22);
    ksSpectra[3]->Scale(8);
*/

    
    //ksSpectra[0][4]->Draw("Psame");
 

    w1->Draw("same");
    r1->Draw("same");
    //r2->Draw("same");
    r4->Draw("same");

}