#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void eposClosureTestRapidity(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/EPOS_TH3D_Sep22_5M_rpyDependent_v1_2014.root");
	//TFile* file = new TFile("~/Desktop/EPOS_TH3D_Sep16_5M_v1_ptDist_2014.root");

	TH3D* ksHist;
	TH3D* laHist;
	TH3D* genksHist;
	TH3D* genlaHist;
	TH3D* xiHist;

	ksHist = (TH3D*)file->Get("ana/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("ana/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("ana/genKS_underlying");
	genlaHist = (TH3D*)file->Get("ana/genLA_underlying");
	xiHist = (TH3D*)file->Get("ana/XiDaughter");

	TH1D* ks_mass[5][20];
	TH1D* la_mass[5][20];
	TH1D* xiHist_mass[5][20];
	
	TH1D* genks_mass[20];
	TH1D* genla_mass[20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    double rpybins[6] = {33,91,145,199,250,314};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;


    for ( int rpy = 0; rpy < 5; rpy++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");
	        xiHistName.str("");

	        ksHistName << "ks_";
	        ksHistName << rpy+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "la_";
	        laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << rpy+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        ks_mass[rpy][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[rpy][pt] = laHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiHist_mass[rpy][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "genks_";
	        ksHistName << pt;

	        laHistName << "genla_";
	        laHistName << pt;

	        genks_mass[pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),1,350,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),1,350,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }

   ks_mass[1][3]->Draw();


	double genksYield[20];
	double genlaYield[20];

	TH1D* genksSpectra = new TH1D("genksSpectra","genK0short pT spectra",15,ptbins_1);
	TH1D* genlaSpectra = new TH1D("genlaSpectra","genLambda pT spectra",15,ptbins_1);


	for(pt = 3; pt < 15; pt++){

		genksYield[pt] = genks_mass[pt]->GetEntries();
			//genksSpectra->SetBinContent( pt+1, genksYield[eta][pt]/(4.8*binwidth_1[pt]) );

		genlaYield[pt] = genla_mass[pt]->GetEntries();
			//genlaSpectra->SetBinContent( pt+1, genlaYield[pt]/(4.8*binwidth_1[pt]) );

	}
	

	double ksYield[5][20];
	double laYield[5][20];


/**
 * Getting efficiency from the table:
 */
    
    TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_withoutXiRemoval_10M_Sep22_rapidity_v1_12pTbins.root");
    //TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_10M.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    double ks_eff[5][20];
    double la_eff[5][20];

    double ks_eff_err[5][20];
    double la_eff_err[5][20];

    for (int i = 0; i < 5; i++){

        for (int r = 3; r < 15; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+1);
            	
            	ks_eff_err[i][r] = hnew1->GetBinError(i+1,r+1);
            
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+1);
            	
            	la_eff_err[i][r] = hnew2->GetBinError(i+1,r+1);

        }
    }

    double ks_pTYield[20];
	double la_pTYield[20];

	double ks_pTYield_err[20];
	double la_pTYield_err[20];

	double temp_ks_err[5];
	double num_ks_err[20];
	double temp_la_err[5];
	double num_la_err[20];

	double temp = 0;


    for (pt = 3; pt < 15; pt++){

    	ks_pTYield[pt] = 0;
    	la_pTYield[pt] = 0;

    	for (rpy = 0; rpy < 5; rpy++){

    		ksYield[rpy][pt] = ks_YieldCal( ks_mass[rpy][pt] );

    		double ksYield_err = sqrt( ksYield[rpy][pt] );
    		temp_ks_err[rpy] = errorCal_num(ksYield[rpy][pt], ksYield_err, ks_eff[rpy][pt], ks_eff_err[rpy][pt]);
    		
    				ksYield[rpy][pt] = ksYield[rpy][pt]/ks_eff[rpy][pt];

    			ks_pTYield[pt] = ksYield[rpy][pt] + ks_pTYield[pt];


    		la_mass[rpy][pt]->Add( xiHist_mass[rpy][pt], -1);
    		laYield[rpy][pt] = la_YieldCal( la_mass[rpy][pt] );

    		double laYield_err = sqrt( laYield[rpy][pt] );
    		temp_la_err[rpy] = errorCal_num(laYield[rpy][pt], laYield_err, la_eff[rpy][pt], la_eff_err[rpy][pt]);
    		
    				laYield[rpy][pt] = laYield[rpy][pt]/la_eff[rpy][pt];

    			la_pTYield[pt] = laYield[rpy][pt] + la_pTYield[pt];

    	}

    	num_ks_err[pt] = errorCal_sum(temp_ks_err[0], temp_ks_err[1], temp_ks_err[2], temp_ks_err[3], temp_ks_err[4], temp_ks_err[5]);
    	num_la_err[pt] = errorCal_sum(temp_la_err[0], temp_la_err[1], temp_la_err[2], temp_la_err[3], temp_la_err[4], temp_la_err[5]);
    }


	TH1D* ks_close = new TH1D("ks_close","ks_close",15,ptbins_1);
	TH1D* la_close = new TH1D("la_close","la_close",15,ptbins_1);

	for (pt = 3; pt < 15; pt++){

		ks_close->SetBinContent(pt+1, ks_pTYield[pt]/genksYield[pt] );

			double error_genks = sqrt( genksYield[pt] );
			double error_ks = errorCal_num( ks_pTYield[pt], num_ks_err[pt], genksYield[pt], error_genks );
		ks_close->SetBinError(pt+1, error_ks);

		la_close->SetBinContent(pt+1, la_pTYield[pt]/genlaYield[pt] );

			double error_genla = sqrt( genlaYield[pt] );
			double error_la = errorCal_num( la_pTYield[pt], num_la_err[pt], genlaYield[pt], error_genla );
		la_close->SetBinError(pt+1, error_la);

	}
	

	/*
	TLatex* etarange[6];

	etarange[0] = new TLatex(2,1.3,"-2.4 < #eta < -1.6");
	etarange[1] = new TLatex(2,1.3,"-1.6 < #eta < -0.8");
	etarange[2] = new TLatex(2,1.3,"-0.8 < #eta < 0");
	etarange[3] = new TLatex(2,1.3,"0 < #eta < 0.8");
	etarange[4] = new TLatex(2,1.3,"0.8 < #eta < 1.6");
	etarange[5] = new TLatex(2,1.3,"1.6 < #eta < 2.4");

	TCanvas* wq1 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    wq1->Divide(2,3,0,0);

    for (eta = 0; eta < 6; eta++){

    	wq1->cd(eta+1);
    	gPad->SetTicks();

    	ks_close[eta]->SetStats(kFALSE);
        ks_close[eta]->GetXaxis()->SetTitleOffset(1.0);
        ks_close[eta]->GetYaxis()->SetTitleOffset(1.0);
        ks_close[eta]->GetXaxis()->SetTitleSize(0.05);
        ks_close[eta]->GetYaxis()->SetTitleSize(0.05);
        ks_close[eta]->GetXaxis()->SetLabelSize(0.06);
        ks_close[eta]->GetYaxis()->SetLabelSize(0.06);

        ks_close[eta]->SetYTitle("K^{0}_{s} RECO/GEN");
        ks_close[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        ks_close[eta]->SetMarkerStyle(2);
        ks_close[eta]->SetMarkerSize(1.2);
        ks_close[eta]->SetMarkerColor(kBlue);
        ks_close[eta]->GetYaxis()->SetRangeUser(0.8,1.4);    

        ks_close[eta]->Draw("P");
        etarange[eta]->Draw("same");

    }

    TCanvas* wq2 = new TCanvas();
    gStyle->SetOptTitle(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetTitleX(0.15);

    wq2->Divide(2,3,0,0);

    for (eta = 0; eta < 6; eta++){

    	wq2->cd(eta+1);
    	gPad->SetTicks();

    	la_close[eta]->SetStats(kFALSE);
        la_close[eta]->GetXaxis()->SetTitleOffset(1.0);
        la_close[eta]->GetYaxis()->SetTitleOffset(1.0);
        la_close[eta]->GetXaxis()->SetTitleSize(0.05);
        la_close[eta]->GetYaxis()->SetTitleSize(0.05);
        la_close[eta]->GetXaxis()->SetLabelSize(0.06);
        la_close[eta]->GetYaxis()->SetLabelSize(0.06);

        la_close[eta]->SetYTitle("#Lambda/#bar{#Lambda} RECO/GEN");
        la_close[eta]->SetXTitle("P^{}_{T,V0}(GeV/c)");
        la_close[eta]->SetMarkerStyle(2);
        la_close[eta]->SetMarkerSize(1.2);
        la_close[eta]->SetMarkerColor(kBlue);
        la_close[eta]->GetYaxis()->SetRangeUser(0.7,1.4);    

        la_close[eta]->Draw("P");
        etarange[eta]->Draw("same");

    }*/
   	cout << "EPOS genLa 1st:  " << genlaYield[3] << endl;
    cout << "EPOS recola 1st: " << la_pTYield[3] << endl;

    cout << "EPOS genLa 2nd:  " << genlaYield[4] << endl;
    cout << "EPOS recola 2nd: " << la_pTYield[4] << endl;

    cout << "EPOS genLa 3rd:  " << genlaYield[5] << endl;
    cout << "EPOS recola 3rd: " << la_pTYield[5] << endl;

	TFile f1("eposClosureTestResult_wErr_withoutXiRemovel_v1_rpy_12pTBins.root","new");

	ks_close->Write();
	la_close->Write();

	
	



}