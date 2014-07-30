#include "fitting.h"

using namespace std;

void plotSpectraRatio(){

	gStyle->SetErrorX(0);

	TFile* file1 = new TFile("./HMbins_pTEffCorr.root");

	TH1D* ksSpectra[9];
	TH1D* laSpectra[9];
	TH1D* ratioHist[9];

	stringstream HistName;

	for (int it = 0; it < 5; it++){

		HistName.str("");
		HistName << "ksSpectra_";
		HistName << it+1;
			ksSpectra[it+4] = (TH1D*)file1->Get( HistName.str().c_str() );

		HistName.str("");
		HistName << "laSpectra_";
		HistName << it+1;
			laSpectra[it+4] = (TH1D*)file1->Get( HistName.str().c_str() );

		HistName.str("");
		HistName << "ratioHist_";
		HistName << it+1;
			ratioHist[it+4] = (TH1D*)file1->Get( HistName.str().c_str() );

	}

	TFile* file2 = new TFile("./MBbins_pTEffCorr.root");

	for (it = 0; it < 4; it++){

		HistName.str("");
		HistName << "ksSpectra_";
		HistName << it+1;
			ksSpectra[it] = (TH1D*)file2->Get( HistName.str().c_str() );

		HistName.str("");
		HistName << "laSpectra_";
		HistName << it+1;
			laSpectra[it] = (TH1D*)file2->Get( HistName.str().c_str() );

		HistName.str("");
		HistName << "ratioHist_";
		HistName << it+1;
			ratioHist[it] = (TH1D*)file2->Get( HistName.str().c_str() );

	}

    TLatex* ratio[9];
    ratio[0] = new TLatex(1.65,0.14,"0 < N^{offline}_{trk} < 35");
    ratio[1] = new TLatex(1.65,0.14,"35 < N^{offline}_{trk} < 60");
    ratio[2] = new TLatex(1.65,0.14,"60 < N^{offline}_{trk} < 90");
    ratio[3] = new TLatex(1.65,0.14,"90 < N^{offline}_{trk} < 120");
    ratio[4] = new TLatex(1.65,0.14,"120 < N^{offline}_{trk} < 150");
    ratio[5] = new TLatex(1.65,0.14,"150 < N^{offline}_{trk} < 185");
    ratio[6] = new TLatex(1.65,0.14,"185 < N^{offline}_{trk} < 220");
    ratio[7] = new TLatex(1.65,0.14,"220 < N^{offline}_{trk} < 260");
    ratio[8] = new TLatex(1.65,0.14,"N^{offline}_{trk} > 260");

    TLatex* YTtile = new TLatex(1.65,1.2,"#Lambda/K^{0}_{s}");
    YTtile->SetTextSize(0.06);

/**
 * plotting ratio in multipanel:
 */
/*
	TCanvas* c1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);
    

        c1->Divide(3,3,0,0);

        for (it = 0; it < 9; it++){

            c1->cd(it+1);
            gPad->SetTicks();

            ratioHist[it]->SetStats(kFALSE);
            ratioHist[it]->GetXaxis()->SetTitleOffset(1.0);
            ratioHist[it]->GetYaxis()->SetTitleOffset(1.0);
            ratioHist[it]->GetXaxis()->SetTitleSize(0.05);
            ratioHist[it]->GetYaxis()->SetTitleSize(0.05);
            ratioHist[it]->GetXaxis()->SetLabelSize(0.06);
            ratioHist[it]->GetYaxis()->SetLabelSize(0.06);

            ratioHist[it]->SetYTitle("#Lambda+#bar{#Lambda}/2K^{0}_{s}");
            ratioHist[it]->SetXTitle("P^{}_{T,V0}(GeV/c)");
            ratioHist[it]->SetMarkerStyle(24);
            ratioHist[it]->SetMarkerSize(1.2);
            ratioHist[it]->SetMarkerColor(kBlue);
            ratioHist[it]->GetYaxis()->SetRangeUser(0,1.5);
        
                ratioHist[it]->Draw("P");
                ratio[it]->SetTextSize(0.06);
                ratio[it]->Draw("same");
                YTtile->Draw("same");

        }
*/


    TCanvas* y1 = new TCanvas();
        gPad->SetTicks();

	TLatex* r1 = new TLatex(5.5,52,"CMS,p-Pb 2013");
		r1->SetTextSize(0.04);
	TLatex* r2 = new TLatex(5.5,15,"-2.4 < #eta < 2.4");
		r2->SetTextSize(0.04);
	TLatex* r3 = new TLatex(0.8,0.000001,"K^{0}_{s}");
		r3->SetTextSize(0.07);
	TLatex* r4 = new TLatex(0.8,0.000001,"#Lambda+#bar{#Lambda}");
		r4->SetTextSize(0.07);
	TLatex* r5 = new TLatex(1.65,0.14,"#Lambda/K^{0}_{s}");
		r5->SetTextSize(0.07);

    TLegend *w1 = new TLegend(0.25,0.4,0.5,0.5);
    w1->SetLineColor(kWhite);
    w1->AddEntry(ksSpectra[0],"0  < N^{offline}_{trk} < 35");
    w1->AddEntry(ksSpectra[1],"35  < N^{offline}_{trk} < 60");
    w1->AddEntry(ksSpectra[2],"60  < N^{offline}_{trk} < 90");
    w1->AddEntry(ksSpectra[3],"90  < N^{offline}_{trk} < 120");
    w1->AddEntry(ksSpectra[4],"120 < N^{offline}_{trk} < 150");
    w1->AddEntry(ksSpectra[5],"150 < N^{offline}_{trk} < 185");
    w1->AddEntry(ksSpectra[6],"185 < N^{offline}_{trk} < 220");
    w1->AddEntry(ksSpectra[7],"220 < N^{offline}_{trk} < 260");
    w1->AddEntry(ksSpectra[8]," 	N^{offline}_{trk} > 260");

    ksSpectra[5]->SetMarkerColor(6);
    ksSpectra[5]->SetMarkerStyle(22);
    ksSpectra[5]->SetLineColor(6);
    ksSpectra[5]->Scale(32);

    ksSpectra[6]->SetMarkerColor(7);
    ksSpectra[6]->SetLineColor(7);
    ksSpectra[6]->SetMarkerStyle(22);
    ksSpectra[6]->Scale(64);

    ksSpectra[7]->SetMarkerColor(8);
    ksSpectra[7]->SetLineColor(8);
    ksSpectra[7]->SetMarkerStyle(22);
    ksSpectra[7]->Scale(128);

    ksSpectra[8]->SetMarkerColor(9);
    ksSpectra[8]->SetLineColor(9);
    ksSpectra[8]->SetMarkerStyle(22);
    ksSpectra[8]->Scale(256);

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

    ksSpectra[4]->SetMarkerColor(kYellow-3);
    ksSpectra[4]->SetLineColor(kYellow-3);
    ksSpectra[4]->SetMarkerStyle(22);
    ksSpectra[4]->SetXTitle("P^{}_{T,V0}(GeV/ck)");
    ksSpectra[4]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}d#eta) [(GeV/c)^{-2}]");
    ksSpectra[4]->SetAxisRange(0,9,"X");
    ksSpectra[4]->SetStats(kFALSE);
    ksSpectra[4]->Scale(16);

    ksSpectra[4]->GetYaxis()->SetRangeUser(0.0000001,1000);
    ksSpectra[4]->GetXaxis()->SetRangeUser(0,9);

    ksSpectra[4]->Draw("P");
    ksSpectra[3]->Draw("Psame");
    ksSpectra[2]->Draw("Psame");
    ksSpectra[1]->Draw("Psame");
    ksSpectra[0]->Draw("Psame");
    ksSpectra[5]->Draw("Psame");
    ksSpectra[6]->Draw("Psame");
    ksSpectra[7]->Draw("Psame");
    ksSpectra[8]->Draw("Psame");

    w1->Draw("same");
    r1->Draw("same");
    r2->Draw("same");
    r3->Draw("same");


    TCanvas* c2 = new TCanvas();

    gPad->SetTicks();

    laSpectra[5]->SetMarkerColor(6);
    laSpectra[5]->SetLineColor(6);
    laSpectra[5]->SetMarkerStyle(22);
    laSpectra[5]->Scale(32);

    laSpectra[6]->SetMarkerColor(7);
    laSpectra[6]->SetLineColor(7);
    laSpectra[6]->SetMarkerStyle(22);
    laSpectra[6]->Scale(64);

    laSpectra[7]->SetMarkerColor(8);
    laSpectra[7]->SetLineColor(8);
    laSpectra[7]->SetMarkerStyle(22);
    laSpectra[7]->Scale(128);

    laSpectra[8]->SetMarkerColor(9);
    laSpectra[8]->SetLineColor(9);
    laSpectra[8]->SetMarkerStyle(22);
    laSpectra[8]->Scale(256);

    laSpectra[0]->SetMarkerColor(1);
    laSpectra[0]->SetLineColor(1);
    laSpectra[0]->SetMarkerStyle(22);
    laSpectra[0]->Scale(1);

    laSpectra[1]->SetMarkerColor(2);
    laSpectra[1]->SetLineColor(2);
    laSpectra[1]->SetMarkerStyle(22);
    laSpectra[1]->Scale(2);

    laSpectra[2]->SetMarkerColor(3);
    laSpectra[2]->SetLineColor(3);
    laSpectra[2]->SetMarkerStyle(22);
    laSpectra[2]->Scale(4);

    laSpectra[3]->SetMarkerColor(4);
    laSpectra[3]->SetLineColor(4);
    laSpectra[3]->SetMarkerStyle(22);
    laSpectra[3]->Scale(8);

    laSpectra[4]->SetMarkerColor(kYellow-3);
    laSpectra[4]->SetLineColor(kYellow-3);
    laSpectra[4]->SetMarkerStyle(22);
    laSpectra[4]->SetXTitle("P^{}_{T,V0}(GeV/c)");
    laSpectra[4]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}d#eta) [(GeV/c)^{-2}]");
    laSpectra[4]->SetAxisRange(0,9,"X");
    laSpectra[4]->SetStats(kFALSE);
    laSpectra[4]->Scale(16);

    laSpectra[4]->GetYaxis()->SetRangeUser(0.0000001,1000);
    laSpectra[4]->GetXaxis()->SetRangeUser(0,9);

    laSpectra[4]->Draw("P");
    laSpectra[3]->Draw("Psame");
    laSpectra[2]->Draw("Psame");
    laSpectra[1]->Draw("Psame");
    laSpectra[0]->Draw("Psame");
    laSpectra[5]->Draw("Psame");
    laSpectra[6]->Draw("Psame");
    laSpectra[7]->Draw("Psame");
    laSpectra[8]->Draw("Psame");
    
    w1->Draw("same");
    r1->Draw("same");
    r2->Draw("same");
    r4->Draw("same");

	TCanvas* c3 = new TCanvas();

	gPad->SetTicks();

	ratioHist[0]->SetMarkerColor(1);
	ratioHist[0]->SetMarkerStyle(22);
	ratioHist[0]->SetAxisRange(0,1.5,"Y");
	ratioHist[0]->SetYTitle("#Lambda+#bar{#Lambda}/2K^{0}_{s}");
	ratioHist[0]->SetXTitle("P^{}_{T,V0}(GeV/c)");
	ratioHist[0]->SetStats(kFALSE);
	ratioHist[0]->Draw("P");

	ratioHist[1]->SetMarkerColor(2);
	ratioHist[1]->SetMarkerStyle(22);
	ratioHist[1]->Draw("Psame");

	ratioHist[2]->SetMarkerColor(3);
	ratioHist[2]->SetMarkerStyle(22);
	ratioHist[2]->Draw("Psame");

	ratioHist[3]->SetMarkerColor(4);
	ratioHist[3]->SetMarkerStyle(22);
	ratioHist[3]->Draw("Psame");

	ratioHist[4]->SetMarkerColor(kYellow-3);
	ratioHist[4]->SetMarkerStyle(22);
	ratioHist[4]->Draw("Psame");

	ratioHist[5]->SetMarkerColor(6);
	ratioHist[5]->SetMarkerStyle(22);
	ratioHist[5]->Draw("Psame");

	ratioHist[6]->SetMarkerColor(7);
	ratioHist[6]->SetMarkerStyle(22);
	ratioHist[6]->Draw("Psame");

	ratioHist[7]->SetMarkerColor(8);
	ratioHist[7]->SetMarkerStyle(22);
	ratioHist[7]->Draw("Psame");

	ratioHist[8]->SetMarkerColor(9);
	ratioHist[8]->SetMarkerStyle(22);
	ratioHist[8]->Draw("Psame");

	TLatex* r11 = new TLatex(6,1.35,"CMS,p-Pb 2013");
		r11->SetTextSize(0.04);
	TLatex* r22 = new TLatex(6,1.25,"-2.4 < #eta < 2.4");
		r22->SetTextSize(0.04);

	w1->Draw("same");
    r11->Draw("same");
    r22->Draw("same");
    r5->Draw("same");




}