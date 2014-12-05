#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void ALICEquantativeComparison(){

	TH3D* ksHist;
	TH3D* laHist;
	TH3D* xiHist;
	
	TH1D* eventHist;

	TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov6_90_120.root");
    ksHist = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist = (TH1D*)file1->Get("ana/eventNumber");

    double norm = eventHist->GetEntries();

    TH1D* ks_HM[34];
    TH1D* la_HM[20];
    TH1D* xi_HM[20];

    double ks_ptbins[35] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2,4.6,5.0,5.5,6.0,8.0};
    double ks_pTbinsBound[35] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,22,24,26,28,30,33,36,39,42,46,50,55,60,80};
    
    double ks_binwidth[34] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.3,0.3,0.4,0.4,0.5,0.5,2.0};
    double ks_ptbincenter[34] = { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
    0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
    1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.15, 3.45, 3.75, 4.05, 
    4.4, 4.8, 5.25, 5.75, 7.0 };

    double la_ptbins[21] = {0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.2,3.7,4.2,5.0,6.0,8.0};
    double la_ptbincenter[20] = { 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.3, 1.5, 1.7, 
    1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.45, 3.95, 4.6, 5.5, 
    7.0 };
    double la_pTbinsBound[21] = {6,7,8,9,10,11,12,14,16,18,20,22,24,26,28,32,37,42,50,60,80};
    double la_binwidth[20] = {0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.5,0.5,0.6,1.0,2.0};

  
    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

	for (int pt = 0; pt < 34; pt++){

		ksHistName.str("");

		ksHistName << "ks1_";
		ksHistName << pt;

		ks_HM[pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),36,41,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	}

	for ( pt = 0; pt < 20; pt++){

		laHistName.str("");
		xiHistName.str("");

		laHistName << "la1_";
		laHistName << pt;

		xiHistName << "xi1_";
		xiHistName << pt;


		la_HM[pt] = laHist->ProjectionZ( laHistName.str().c_str(),36,41,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
		xi_HM[pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),36,41,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	}


 	//TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/HIJINGaliceQuantativeComparison_vertexReweighted.root");
	TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/eposEfficiencyRapidityTable/EPOSaliceQuantativeComparison_vertexReweighted.root");


    TH1D* hnew1 = (TH1D*)t1->Get("ks_eff");
    TH1D* hnew2 = (TH1D*)t1->Get("la_eff");

    double ks_eff[34];
    double la_eff[20];
    double ks_eff_err[34];
    double la_eff_err[20];

    for (pt = 0; pt < 34; pt++){
        
        ks_eff[pt] = hnew1->GetBinContent(pt+1);
            ks_eff_err[pt] = hnew1->GetBinError(pt+1);
    }

    for (pt = 0; pt < 20; pt++){

        la_eff[pt] = hnew2->GetBinContent(pt+1);
            la_eff_err[pt] = hnew2->GetBinError(pt+1);
    }



    double ks_HM_yield[34];
    double la_HM_yield[20];

    double temp_ks_err[34];
    double temp_la_err[20];


    for (pt = 0; pt < 34; pt++){
    
		ks_HM_yield[pt] = ks_YieldCal( ks_HM[pt] );	
			double ksYield_err = sqrt( ks_HM_yield[pt] );

		temp_ks_err[pt] = errorCal_num( ks_HM_yield[pt], ksYield_err, ks_eff[pt], ks_eff_err[pt] );

		ks_HM_yield[pt] = ks_HM_yield[pt]/ks_eff[pt];

    }
     
    for (pt = 0; pt < 20; pt++){
         
		la_HM[pt]->Add( xi_HM[pt], -1);
		la_HM_yield[pt] = la_YieldCal( la_HM[pt] );

		double laYield_err = sqrt( la_HM_yield[pt] );
		temp_la_err[pt] = errorCal_num( la_HM_yield[pt], laYield_err, la_eff[pt], la_eff_err[pt] );

		la_HM_yield[pt] = la_HM_yield[pt]/la_eff[pt];

    }

/*
**************************************
 */

    TH1D* ksSpectra = new TH1D("ksSpectra_alice","ksSpectra_alice",34,ks_ptbins);;
    TH1D* laSpectra = new TH1D("laSpectra_alice","laSpectra_alice",20,la_ptbins);

	for (pt = 0; pt < 34; pt++){

	double ks_temp = (ks_HM_yield[pt]/ks_binwidth[pt])/(2*3.1415926*ks_ptbincenter[pt]*0.5*norm);

	ksSpectra->SetBinContent(pt+1, ks_temp );
	ksSpectra->SetBinError(pt+1, temp_ks_err[pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*0.5*norm ));

	}

	for (pt = 0; pt < 20; pt++){

	double la_temp = (la_HM_yield[pt]/la_binwidth[pt])/(2*3.1415926*la_ptbincenter[pt]*0.5*norm);

	laSpectra->SetBinContent(pt+1, la_temp );
	laSpectra->SetBinError(pt+1, temp_la_err[pt]/( 2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*0.5*norm ));

	}
        
    TFile f1("aliceQuantativeComparison_90_120_vertexReweight_v4_epos_pt.root","new");

    ksSpectra->Write();
    laSpectra->Write();
       




}

