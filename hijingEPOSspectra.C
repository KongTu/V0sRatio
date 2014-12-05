#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEPOSspectra(){

//HIJING:
	//TFile* file = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_18M_TH3D_Nov5_rpyDependent_vtxReweight_2014.root");
//EPOS:
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

	TH1D* ks_mass[5][28];
	TH1D* la_mass[5][20];
	TH1D* genks_mass[5][28];
	TH1D* genla_mass[5][20];
	TH1D* xiHist_mass[5][20];

    double ks_pTbinsBound[29] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[28] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};
    double ks_ptbincenter[28] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

    double la_pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
    double la_binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double rpybins[6] = {6,16,26,35,44,55};
    double rpybinwidth[5] = {1.0,0.97,0.9,0.93,1.0};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for ( int rpy = 0; rpy < 5; rpy++){

    	for (int pt = 0; pt < 28; pt++){

    		ksHistName.str("");

    		ksHistName << "ks_";
	        ksHistName << rpy+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        ks_mass[rpy][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	        ksHistName.str("");

         	ksHistName << "genks_";
	        ksHistName << rpy+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        genks_mass[rpy][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

    	}

	    for (pt = 0; pt < 20; pt++){

	        laHistName.str("");
	        xiHistName.str("");

	        laHistName << "la_";
	        laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << rpy+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        la_mass[rpy][pt] = laHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	        xiHist_mass[rpy][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

	  
	        laHistName.str("");

	        laHistName << "genla_";
          	laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        genla_mass[rpy][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	   	
	   	}
   }

   	int genksYield[5][28];
   	int genlaYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){


   		for( pt = 0; pt < 28; pt++){

   			genksYield[rpy][pt] = genks_mass[rpy][pt]->GetEntries();
   		}
	   	
	   	for (int i = 0; i < 20; i++){

	   		genlaYield[rpy][i] = genla_mass[rpy][i]->GetEntries();

	   	}
   	}	

   	double ksYield[5][28];
   	double laYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){

   		for (pt = 0; pt < 28; pt++){

   			ksYield[rpy][pt] = ks_YieldCal( ks_mass[rpy][pt] );

   		}

	   	for ( i = 0; i < 20; i++){
	
			la_mass[rpy][i]->Add( xiHist_mass[rpy][i],-1 );
	   			laYield[rpy][i] = la_YieldCal( la_mass[rpy][i] );
	   	}
   	}

/*
Filling into spectrum
 */

	TH1D* ksSpectra_reco[5];
	TH1D* laSpectra_reco[5];

	TH1D* ksSpectra_gen[5];
	TH1D* laSpectra_gen[5];

	stringstream ksName;
	stringstream laName;

	for(rpy = 0; rpy < 5; rpy++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_reco_";
		ksName << rpy+1;

		laName << "laSpectra_epos_reco_";
		laName << rpy+1;

		ksSpectra_reco[rpy] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		laSpectra_reco[rpy] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		for(pt = 0; pt < 28; pt++){

			ksSpectra_reco[rpy]->SetBinContent(pt+1, ksYield[rpy][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));
				ksSpectra_reco[rpy]->SetBinError(pt+1, sqrt( ksYield[rpy][pt] )/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));

		}

		for(pt = 0; pt < 20; pt++){

			laSpectra_reco[rpy]->SetBinContent(pt+1, laYield[rpy][pt]/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));
				laSpectra_reco[rpy]->SetBinError(pt+1, sqrt( laYield[rpy][pt] )/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));

		}

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_epos_gen_";
		ksName << rpy+1;

		laName << "laSpectra_epos_gen_";
		laName << rpy+1;

		ksSpectra_gen[rpy] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		laSpectra_gen[rpy] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		for(pt = 0; pt < 28; pt++){

			ksSpectra_gen[rpy]->SetBinContent(pt+1, genksYield[rpy][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));
				ksSpectra_gen[rpy]->SetBinError(pt+1, sqrt( genksYield[rpy][pt] )/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy] ));

		}

		for(pt = 0; pt < 20; pt++){

			laSpectra_gen[rpy]->SetBinContent(pt+1, genlaYield[rpy][pt]/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));
				laSpectra_gen[rpy]->SetBinError(pt+1, sqrt( genlaYield[rpy][pt] )/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy] ));

		}
	}
	

   	TFile f1("EPOSspectra_5rpyBins.root","new");
 
   	for ( rpy = 0; rpy < 5; rpy++){

   		ksSpectra_reco[rpy]->Write();
   		laSpectra_reco[rpy]->Write();
   	}

   	for ( rpy = 0; rpy < 5; rpy++){

   		ksSpectra_gen[rpy]->Write();
   		laSpectra_gen[rpy]->Write();
   	}






}