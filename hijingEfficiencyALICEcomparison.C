#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiencyALICEcomparison(){

	//TFile* file = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_18M_TH3D_Nov5_rpyDependent_vtxReweight_2014.root");
	TFile* file = new TFile("/Users/kongkong/2014Research/ROOT_file/newEPOSsample/EPOS_TH3D_Nov6_rpyDependent_vertexReweight_2014.root");

	TH3D* ksHist;
	TH3D* laHist;
	TH3D* xiHist;

	TH3D* genksHist;
	TH3D* genlaHist;

	ksHist = (TH3D*)file->Get("ana/InvMass_ks_underlying");
	laHist = (TH3D*)file->Get("ana/InvMass_la_underlying");
	genksHist = (TH3D*)file->Get("ana/genKS_underlying");
	genlaHist = (TH3D*)file->Get("ana/genLA_underlying");
	xiHist = (TH3D*)file->Get("ana/XiDaughter");

	TH1D* ks_mass[34];
	TH1D* la_mass[20];
	TH1D* genks_mass[34];
	TH1D* genla_mass[20];
	TH1D* xiHist_mass[20];

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

		ksHistName << "ks_";
		ksHistName << pt;

		ks_mass[pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),36,41,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

		ksHistName.str("");

		ksHistName << "genks_";
		ksHistName << pt;

		genks_mass[pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),36,41,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	}

	for (pt = 0; pt < 20; pt++){

		laHistName.str("");
		xiHistName.str("");

		laHistName << "la_";
		laHistName << pt;

		xiHistName << "xiDau_";
		xiHistName << pt;

		la_mass[pt] = laHist->ProjectionZ( laHistName.str().c_str(),36,41,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
		xiHist_mass[pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),36,41,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

		laHistName.str("");
		laHistName << "genla_";
		laHistName << pt;

		genla_mass[pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),36,41,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

	}
   

   	int genksYield[34];
   	int genlaYield[20];

	for( pt = 0; pt < 34; pt++){

		genksYield[pt] = genks_mass[pt]->GetEntries();
	}

	for (int i = 0; i < 20; i++){

		genlaYield[i] = genla_mass[i]->GetEntries();

	}
   

   	double ksYield[34];
   	double laYield[20];

	for (pt = 0; pt < 34; pt++){

		ksYield[pt] = ks_YieldCal( ks_mass[pt] );

	}

	for ( i = 0; i < 20; i++){

		la_mass[i]->Add( xiHist_mass[i],-1 );
			laYield[i] = la_YieldCal( la_mass[i] );
	}
   

 	TH1D* ks_eff = new TH1D("ks_eff","ks_eff",34,ks_ptbins);
 	TH1D* la_eff = new TH1D("la_eff","la_eff",20,la_ptbins);
 

	for (pt = 0; pt < 34; pt++){

		ks_eff->SetBinContent(pt+1, ksYield[pt]/genksYield[pt]);
		ks_eff->SetBinError(pt+1, errorCal( ksYield[pt], genksYield[pt] ) );
	}

	for (pt = 0; pt < 20; pt++){

		la_eff->SetBinContent(pt+1, laYield[pt]/genlaYield[pt]);
		la_eff->SetBinError(pt+1, errorCal( laYield[pt], genlaYield[pt] ) );
	}


   	TFile f1("EPOSaliceQuantativeComparison_vertexReweighted.root","new");
   		ks_eff->Write();
		la_eff->Write();
 


}