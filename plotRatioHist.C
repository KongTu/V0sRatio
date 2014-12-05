#include "fitting.h"

using namespace std;

void plotRatioHist(){

    gStyle->SetErrorX(0);

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy_new/new8Multbins_EPOSvtx_FullStats_v6_26ksbins_pt.root");

	TH1D* ratioHist[9];
    TH1D* ksSpectra[9];
    TH1D* laSpectra[9];

	stringstream ratioHistName;
    stringstream ksSpectraName;
    stringstream laSpectraName;

	for(int mult = 0; mult < 8; mult++){

	
			ratioHistName.str("");
			ratioHistName << "ratio_vtx_";
			ratioHistName << mult+1;
			//ratioHistName << "_";
			//ratioHistName << rpy+1;

            ksSpectraName.str("");
            ksSpectraName << "ksSpectra_vtx_";
            ksSpectraName << mult+1;
            //ksSpectraName << "_";
            //ksSpectraName << rpy+1;

            laSpectraName.str("");
            laSpectraName << "laSpectra_vtx_";
            laSpectraName << mult+1;
            //laSpectraName << "_";
            //laSpectraName << rpy+1;

			ratioHist[mult] = (TH1D*)file->Get( ratioHistName.str().c_str() );
            ksSpectra[mult] = (TH1D*)file->Get( ksSpectraName.str().c_str() );
            laSpectra[mult] = (TH1D*)file->Get( laSpectraName.str().c_str() );


		
	}

   	/*TLatex* ratio[5];
    ratio[0] = new TLatex(1.65,1.3,"-2.87 < y < -1.8");
    ratio[1] = new TLatex(1.65,1.3,"-1.8 < y < -0.9");
    ratio[2] = new TLatex(1.65,1.3,"-0.9 < y < 0");
    ratio[3] = new TLatex(1.65,1.3,"0 < y < 0.93");
    ratio[4] = new TLatex(1.65,1.3,"0.93 < y < 1.93");*/

    TLatex* ratio[8];

    ratio[0] = new TLatex(4,1.2,"0 < N^{offline}_{trk} < 35");
    ratio[1] = new TLatex(4,1.2,"35 < N^{offline}_{trk} < 60");
    ratio[2] = new TLatex(4,1.2,"60 < N^{offline}_{trk} < 90");
    ratio[3] = new TLatex(4,1.2,"90 < N^{offline}_{trk} < 120");
    ratio[4] = new TLatex(4,1.2,"120 < N^{offline}_{trk} < 150");
    ratio[5] = new TLatex(4,1.2,"150 < N^{offline}_{trk} < 185");
    ratio[6] = new TLatex(4,1.2,"185 < N^{offline}_{trk} < 220");
    ratio[7] = new TLatex(4,1.2,"220 < N^{offline}_{trk} < #infty");

	TCanvas* c1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    c1->Divide(4,2,0,0);

    for(int it = 0; it < 8; it++){

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
        ratioHist[it]->SetMarkerStyle(4);
        ratioHist[it]->SetMarkerSize(1.2);
        ratioHist[it]->SetMarkerColor(kBlue);
        ratioHist[it]->GetYaxis()->SetRangeUser(0,1.5);
    
        ratioHist[it]->Draw("P");
        ratio[it]->SetTextSize(0.06);
        ratio[it]->Draw("same");
   
    }

    TLatex* r1 = new TLatex(2.5,0.52,"CMS,p-Pb 2013,#sqrt{S^{}_{NN}} = 5.02 TeV");
    r1->SetTextSize(0.05);
    TLatex* r2 = new TLatex(5.25,0.13,"-2.87 < y < 1.93");
    r2->SetTextSize(0.05);
    TLatex* r3 = new TLatex(0.8,0.00000001,"K^{0}_{s}");
    r3->SetTextSize(0.07);
    TLatex* r4 = new TLatex(0.8,0.00000001,"#Lambda+#bar{#Lambda}");
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
    w1->SetFillColor(0);
    w1->AddEntry(ksSpectra[0],"0   < N^{offline}_{trk} < 35 (#times2^{-7})");
    w1->AddEntry(ksSpectra[1],"35  < N^{offline}_{trk} < 60 (#times2^{-6})");
    w1->AddEntry(ksSpectra[2],"60  < N^{offline}_{trk} < 90 (#times2^{-5})");
    w1->AddEntry(ksSpectra[3],"90  < N^{offline}_{trk} < 120 (#times2^{-4})");
    w1->AddEntry(ksSpectra[4],"120 < N^{offline}_{trk} < 150 (#times2^{-3})");
    w1->AddEntry(ksSpectra[5],"150 < N^{offline}_{trk} < 185 (#times2^{-2})");
    w1->AddEntry(ksSpectra[6],"185 < N^{offline}_{trk} < 220 (#times2^{-1})");
    w1->AddEntry(ksSpectra[7],"220 < N^{offline}_{trk} (#times1)");
    //w1->AddEntry(ksSpectra[8],"     N^{offline}_{trk} > 260");


    /*TCanvas* p1 = new TCanvas();
    laSpectra[0][rpy]->SetMarkerColor(kYellow-3);
    laSpectra[0][rpy]->SetLineColor(kYellow-3);
    laSpectra[0][rpy]->SetMarkerStyle(20);
    laSpectra[0][rpy]->SetXTitle("m^{}_{T} - m^{}_{0} (GeV/c^{2})");
    //laSpectra[0][rpy]->SetXTitle("P^{}_{T,V0}(GeV/c)");
    laSpectra[0][rpy]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
    //laSpectra[0][rpy]->SetAxisRange(0,8,"X");
    laSpectra[0][rpy]->SetStats(kFALSE);
    laSpectra[0][rpy]->Scale(16);

    laSpectra[0][rpy]->GetYaxis()->SetRangeUser(0.0000001,1000);
    laSpectra[0][rpy]->GetXaxis()->SetRangeUser(0,9);

    laSpectra[0][rpy]->Draw("P");*/


    TCanvas* y1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);
    gPad->SetLogy(1);
    gPad->SetTicks();
    y1->SetLogy(1);

    y1->Divide(2,1,0,0);

    //for (rpy = 0; rpy < 5; rpy++){

    y1->cd(1);

    gPad->SetLogy(1);


    ksSpectra[0]->SetMarkerColor(kYellow-2);
    ksSpectra[0]->SetLineColor(kWhite);
    ksSpectra[0]->SetMarkerStyle(20);
    ksSpectra[0]->SetMarkerSize(1.3);
    //ksSpectra[0]->SetXTitle("m^{}_{T} - m^{}_{0} (GeV/c^{2})");
    ksSpectra[0]->SetXTitle("P^{}_{T,V0}(GeV/c)");
    ksSpectra[0]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
    //ksSpectra[0]->SetAxisRange(0,8,"X");
    ksSpectra[0]->SetStats(kFALSE);
    ksSpectra[0]->Scale(0.0078125);

    ksSpectra[0]->GetYaxis()->SetRangeUser(0.000000001,10);
    ksSpectra[0]->GetXaxis()->SetRangeUser(0,9);

    ksSpectra[1]->SetMarkerColor(kBlack);
    ksSpectra[1]->SetMarkerStyle(21);
    ksSpectra[1]->SetMarkerSize(1.3);
    ksSpectra[1]->SetLineColor(kWhite);
    ksSpectra[1]->Scale(0.015625);

    ksSpectra[2]->SetMarkerColor(kRed);
    ksSpectra[2]->SetLineColor(kWhite);
    ksSpectra[2]->SetMarkerStyle(22);
    ksSpectra[2]->SetMarkerSize(1.3);
    ksSpectra[2]->Scale(0.03125);

    ksSpectra[3]->SetMarkerColor(kGreen+4);
    ksSpectra[3]->SetLineColor(kWhite);
    ksSpectra[3]->SetMarkerStyle(23);
    ksSpectra[3]->SetMarkerSize(1.3);
    ksSpectra[3]->Scale(0.0625);

    ksSpectra[4]->SetMarkerColor(kBlue+1);
    ksSpectra[4]->SetLineColor(kWhite);
    ksSpectra[4]->SetMarkerStyle(24);
    ksSpectra[4]->SetMarkerSize(1.3);
    ksSpectra[4]->Scale(0.125);

    ksSpectra[5]->SetMarkerColor(kCyan+4);
    ksSpectra[5]->SetLineColor(kWhite);
    ksSpectra[5]->SetMarkerStyle(25);
    ksSpectra[5]->SetMarkerSize(1.3);
    ksSpectra[5]->Scale(0.25);

    ksSpectra[6]->SetMarkerColor(kRed);
    ksSpectra[6]->SetLineColor(kWhite);
    ksSpectra[6]->SetMarkerStyle(26);
    ksSpectra[6]->SetMarkerSize(1.3);
    ksSpectra[6]->Scale(0.5);

    ksSpectra[7]->SetMarkerColor(kBlue);
    ksSpectra[7]->SetLineColor(kWhite);
    ksSpectra[7]->SetMarkerStyle(27);
    ksSpectra[7]->SetMarkerSize(1.3);
    ksSpectra[7]->Scale(1);
                  
    ksSpectra[0]->Draw("P");
      
    ksSpectra[1]->Draw("Psame");      
    ksSpectra[2]->Draw("Psame");    
    ksSpectra[3]->Draw("Psame");      
    ksSpectra[4]->Draw("Psame");
    ksSpectra[5]->Draw("Psame");
    ksSpectra[6]->Draw("Psame");
    ksSpectra[7]->Draw("Psame");

    w1->Draw("same");
    //r1->Draw("same");
    r3->Draw("same");

    y1->cd(2);

    gPad->SetLogy(1);

    
    laSpectra[0]->SetMarkerColor(kYellow-2);
    laSpectra[0]->SetLineColor(kWhite);
    laSpectra[0]->SetMarkerStyle(20);
    laSpectra[0]->SetMarkerSize(1.3);
    //laSpectra[0]->SetXTitle("m^{}_{T} - m^{}_{0} (GeV/c^{2})");
    laSpectra[0]->SetXTitle("P^{}_{T,V0}(GeV/c)");
    laSpectra[0]->SetYTitle("1/N^{}_{ev}1/(2#PiP^{}_{T}d^{2}N/(dP^{}_{T}dy) [(GeV/c)^{-2}]");
    //laSpectra[0]->SetAxisRange(0,8,"X");
    laSpectra[0]->SetStats(kFALSE);
    laSpectra[0]->Scale(0.0078125);

    laSpectra[0]->GetYaxis()->SetRangeUser(0.000000001,10);
    laSpectra[0]->GetXaxis()->SetRangeUser(0,9);

    laSpectra[1]->SetMarkerColor(kBlack);
    laSpectra[1]->SetMarkerStyle(21);
    laSpectra[1]->SetMarkerSize(1.3);
    laSpectra[1]->SetLineColor(kWhite);
    laSpectra[1]->Scale(0.015625);

    laSpectra[2]->SetMarkerColor(kRed);
    laSpectra[2]->SetLineColor(kWhite);
    laSpectra[2]->SetMarkerStyle(22);
    laSpectra[2]->SetMarkerSize(1.3);
    laSpectra[2]->Scale(0.03125);

    laSpectra[3]->SetMarkerColor(kGreen+4);
    laSpectra[3]->SetLineColor(kWhite);
    laSpectra[3]->SetMarkerStyle(23);
    laSpectra[3]->SetMarkerSize(1.3);
    laSpectra[3]->Scale(0.0625);

    laSpectra[4]->SetMarkerColor(kBlue+1);
    laSpectra[4]->SetLineColor(kWhite);
    laSpectra[4]->SetMarkerStyle(24);
    laSpectra[4]->SetMarkerSize(1.3);
    laSpectra[4]->Scale(0.125);

    laSpectra[5]->SetMarkerColor(kCyan+4);
    laSpectra[5]->SetLineColor(kWhite);
    laSpectra[5]->SetMarkerStyle(25);
    laSpectra[5]->SetMarkerSize(1.3);
    laSpectra[5]->Scale(0.25);

    laSpectra[6]->SetMarkerColor(kRed);
    laSpectra[6]->SetLineColor(kWhite);
    laSpectra[6]->SetMarkerStyle(26);
    laSpectra[6]->SetMarkerSize(1.3);
    laSpectra[6]->Scale(0.5);

    laSpectra[7]->SetMarkerColor(kBlue);
    laSpectra[7]->SetLineColor(kWhite);
    laSpectra[7]->SetMarkerStyle(27);
    laSpectra[7]->SetMarkerSize(1.3);
    laSpectra[7]->Scale(1);
                  
    laSpectra[0]->Draw("P");
      
    laSpectra[1]->Draw("Psame");      
    laSpectra[2]->Draw("Psame");    
    laSpectra[3]->Draw("Psame");      
    laSpectra[4]->Draw("Psame");
    laSpectra[5]->Draw("Psame");
    laSpectra[6]->Draw("Psame");
    laSpectra[7]->Draw("Psame");


    //w1->Draw("same");
    r1->Draw("same");
    r2->Draw("same");
    r4->Draw("same");

}