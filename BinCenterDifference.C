#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void BinCenterDifference(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HIJING_TH3D_Sep10_10M_v2_ptDist_2014.root");

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

	TH1D* ks_mass[6][20];
	TH1D* la_mass[6][20];
	TH1D* xiHist_mass[6][20];
	
	TH1D* genks_mass[20];
	TH1D* genla_mass[20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");
	        xiHistName.str("");

	        ksHistName << "ks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "la_";
	        laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << eta+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        xiHist_mass[eta][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "genks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "genla_";
	        laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        genks_mass[pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),1,6,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }


   	double ksYield[6][20];
	double laYield[6][20];


/**
 * Getting efficiency from the table:
 */
    
    TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyTable/effKongNew2DTable_15M_Sep10_v1_12pTbins.root");
    //TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_10M.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    double ks_eff[6][20];
    double la_eff[6][20];

    double ks_eff_err[6][20];
    double la_eff_err[6][20];

    for (int i = 0; i < 6; i++){

        for (int r = 3; r < 15; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+1);
            	
            	ks_eff_err[i][r] = hnew1->GetBinError(i+1,r+1);
            
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+1);
            	
            	la_eff_err[i][r] = hnew2->GetBinError(i+1,r+1);

        }
    }


/**
 * Obtain the yield and sum over all eta bins: 
 */

 
  	double ks_pTYield[20];
	double la_pTYield[20];


	for (pt = 3; pt < 15; pt++){

		for (eta = 0; eta < 6; eta++){

			ksYield[eta][pt] = ks_YieldCal( ks_mass[eta][pt] );
				ksYield[eta][pt] = ksYield[eta][pt]/ks_eff[eta][pt];

    				ks_pTYield[pt] = ksYield[eta][pt] + ks_pTYield[pt];

    		la_mass[eta][pt]->Add( xiHist_mass[eta][pt], -1);
    		
    		laYield[eta][pt] = la_YieldCal( la_mass[eta][pt] );
					laYield[eta][pt] = laYield[eta][pt]/la_eff[eta][pt];

    				la_pTYield[pt] = laYield[eta][pt] + la_pTYield[pt];
		}
	}

	TH1D* k0short_spectra =  new TH1D("k0short_spectra","k0short_spectra",15,ptbins_1);
	TH1D* lambda_spectra =  new TH1D("lambda_spectra","lambda_spectra",15,ptbins_1);

	double temp_la[12];

	for (pt = 3; pt < 15; pt++){

		temp_la[pt-3] = la_pTYield[pt]/binwidth_1[pt];

		k0short_spectra->SetBinContent(pt+1,ks_pTYield[pt]/binwidth_1[pt]);
		lambda_spectra->SetBinContent(pt+1, temp_la[pt-3]);
	}

	//lambda_spectra->Draw("P");


	TH1D* ptDist[12];

	double ptDistMean[12];

	stringstream ptDistName;

	for( int i = 0; i < 12; i++){

		ptDistName.str("");

		ptDistName << "ana/ptDist";
		ptDistName << i;

		ptDist[i] = (TH1D*)file->Get( ptDistName.str().c_str() );

		ptDistMean[i] = ptDist[i]->GetMean();

		cout << "ptDistMean: " << ptDistMean[i] << endl;

	}

	//cout << "last bin yield:   " << temp_la[11] << endl;

	TCanvas* c1 = new TCanvas();

	double ptDistBinCenter[12] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.3,2.9,3.7,5.1,7.5};


	TGraph* g1 = new TGraph(12, ptDistMean, temp_la );
	TGraph* g2 = new TGraph(12, ptDistBinCenter, temp_la);

	g1->SetMarkerColor(kRed);
	g2->SetMarkerColor(kBlue);

	g1->Draw("AP");
	g2->Draw("P");

	TFile f1("BinCenterDifference_v1.root","new");
	g1->Write();
	g2->Write();
	lambda_spectra->Write();


	












}
