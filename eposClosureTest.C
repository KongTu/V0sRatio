#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void eposClosureTest(){

	gStyle->SetErrorX(0);

	TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/EPOS_TH3D_Sep2_vtx_2014.root");

	TH3D* ksHist;
	TH3D* laHist;
	TH3D* genksHist;
	TH3D* genlaHist;

	ksHist = (TH3D*)file->Get("ana/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("ana/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("ana/genKS_underlying");
	genlaHist = (TH3D*)file->Get("ana/genLA_underlying");

	TH1D* ks_mass[6][20];
	TH1D* la_mass[6][20];
	
	TH1D* genks_mass[6][20];
	TH1D* genla_mass[6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    stringstream ksHistName;
    stringstream laHistName;

    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 3; pt < 15; pt++){

	        ksHistName.str("");
	        laHistName.str("");

	        ksHistName << "ks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        laHistName << "la_";
	        laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );

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

	        genks_mass[eta][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	        genla_mass[eta][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
	   	
	   	}
   }


	double genksYield[6][20];
	double genlaYield[6][20];

	TH1D* genksSpectra = new TH1D("genksSpectra","genK0short pT spectra",15,ptbins_1);
	TH1D* genlaSpectra = new TH1D("genlaSpectra","genLambda pT spectra",15,ptbins_1);

	for (eta = 0; eta < 6; eta++){

		for(pt = 3; pt < 15; pt++){

			genksYield[eta][pt] = genks_mass[eta][pt]->GetEntries();
				//genksSpectra->SetBinContent( pt+1, genksYield[eta][pt]/(4.8*binwidth_1[pt]) );


			genlaYield[eta][pt] = genla_mass[eta][pt]->GetEntries();
				//genlaSpectra->SetBinContent( pt+1, genlaYield[pt]/(4.8*binwidth_1[pt]) );

		}
	}

	double ksYield[6][20];
	double laYield[6][20];

/**
 * Getting efficiency from the table:
 */
    
    TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyTable/effKongNew2DTable_12pTbins.root");
    //TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_10M.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    double ks_eff[6][20];
    double la_eff[6][20];

    for (int i = 0;i < 6;i++){

        for (int r = 3; r < 15; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+1);
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+1);

        }
    }

	for (eta = 0; eta < 6; eta++){

		for (pt = 3; pt < 15; pt++){

			ksYield[eta][pt] = ks_YieldCal( ks_mass[eta][pt] );
				ksYield[eta][pt] = ksYield[eta][pt]/ks_eff[eta][pt];
			laYield[eta][pt] = la_YieldCal( la_mass[eta][pt] );
				laYield[eta][pt] = laYield[eta][pt]/la_eff[eta][pt];
		}
	}

	double ks_pTYield[20];
	double la_pTYield[20];
/*
	for (pt = 3; pt < 15; pt++){

		ks_pTYield[pt] = 0.0;
		la_pTYield[pt] = 0.0;

		for (eta = 0; eta < 6; eta++){

			ks_pTYield[pt] = ksYield[eta][pt] + ks_pTYield[pt];
			la_pTYield[pt] = laYield[eta][pt] + la_pTYield[pt];
		}

	}*/

	TH1D* ksSpectra = new TH1D("ksSpectra","K0short pT spectra",15,ptbins_1);
	TH1D* laSpectra = new TH1D("laSpectra","Lambda pT spectra",15,ptbins_1);
	TH1D* ratioHist = new TH1D("ratioHist","#Lambda/K^{0}_{s}",15,ptbins_1);
	TH1D* genratioHist = new TH1D("genratioHist","#Lambda/K^{0}_{s}",15,ptbins_1);

	TH1D* ks_close[6]; = new TH1D("ks_close","ks_close",15,ptbins_1);
	TH1D* la_close[6]; = new TH1D("la_close","la_close",15,ptbins_1);

	stringstream kscloseName;
	stringstream lacloseName;

	for (eta = 0; eta < 6; eta++){

		kscloseName.str("");
		kscloseName << "ks_close_";
		kscloseName << eta;

		lacloseName.str("");
		lacloseName << "la_close_";
		lacloseName << eta;
		
		ks_close[eta] = new TH1D( kscloseName.str().c_str(), kscloseName.str().c_str(), 15, ptbins_1);
		la_close[eta] = new TH1D( lacloseName.str().c_str(), lacloseName.str().c_str(), 15, ptbins_1);

		for (pt = 3; pt < 15; pt++){

			ks_close[eta]->SetBinContent(pt+1, ksYield[eta][pt]/genksYield[eta][pt] );
				double error_ks = errorCal( ksYield[eta][pt], genksYield[eta][pt] );
			ks_close[eta]->SetBinError(pt+1, error_ks);

			la_close[eta]->SetBinContent(pt+1, laYield[eta][pt]/genlaYield[eta][pt] );
				double error_la = errorCal( laYield[eta][pt], genlaYield[eta][pt] );
			la_close[eta]->SetBinError(pt+1, error_la);

		}
	}




	/*for(pt = 3; pt < 15; pt++){

			ksSpectra->SetBinContent( pt+1, ks_pTYield[pt]/(4.8*binwidth_1[pt]) );
			laSpectra->SetBinContent( pt+1, la_pTYield[pt]/(4.8*binwidth_1[pt]) );

			

			ratioHist->SetBinContent( pt+1, la_pTYield[pt]/(2*ks_pTYield[pt]) );
				double err = errorCal( la_pTYield[pt], ks_pTYield[pt] );
			ratioHist->SetBinError( pt+1, err );

	}*/


	TFile f1("eposClosureTestResult_etaDependent.root","new");
	//genksSpectra->Write();
	//genlaSpectra->Write();
	//ksSpectra->Write();
	//laSpectra->Write();
	for (eta = 0; eta < 6; eta++){

		ks_close[eta]->Write();
		la_close[eta]->Write();
	}
	
	//ratioHist->Write();


	



}