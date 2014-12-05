#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingMultLowHigh(){

	TFile* file[2];
	file[0] = new TFile("/Users/kongkong/2014Research/ROOT_file/Hijing_pPb_PbPb_comparison/HIJING_PbPb_TH3D_Dec2_eta_0_60_vertexreweighed_2014.root");
	file[1] = new TFile("/Users/kongkong/2014Research/ROOT_file/Hijing_pPb_PbPb_comparison/HIJING_PbPb_TH3D_Dec1_eta_185plus_vertexreweighed_2014.root");

	TH3D* ksHist[2];
	TH3D* laHist[2];
	TH3D* xiHist[2];

	TH3D* genksHist[2];
	TH3D* genlaHist[2];

	for(int mult = 0; mult < 2; mult++){

		ksHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_ks_underlying");
		laHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_la_underlying");
		genksHist[mult] = (TH3D*)file[mult]->Get("ana/genKS_underlying");
		genlaHist[mult] = (TH3D*)file[mult]->Get("ana/genLA_underlying");
		xiHist[mult] = (TH3D*)file[mult]->Get("ana/XiDaughter");

	}

	TH1D* ks_mass[2][28];
	TH1D* la_mass[2][20];
	TH1D* genks_mass[2][28];
	TH1D* genla_mass[2][20];
	TH1D* xiHist_mass[2][20];

    double ks_pTbinsBound[29] = {0,1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[29] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[28] = {0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};
    double ks_ptbincenter[28] = {0.05,0.15,0.25,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

    double la_pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
    double la_binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double etabins[] = {11,19,27,35,43,51,59};
    double rpybins[6] = {6,16,26,35,44,55};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

	for(mult = 0; mult < 2; mult++){

    	for (int pt = 0; pt < 28; pt++){

    		ksHistName.str("");

    		ksHistName << "ks_";
    		ksHistName << mult+1;
    		ksHistName << "_";
	        ksHistName << pt;

	        ks_mass[mult][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),31,39,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	        ksHistName.str("");

         	ksHistName << "genks_";
         	ksHistName << mult+1;
    		ksHistName << "_";
	        ksHistName << pt;

	        genks_mass[mult][pt] = genksHist[mult]->ProjectionZ( ksHistName.str().c_str(),31,39,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

    	}

	    for (pt = 0; pt < 20; pt++){

	        laHistName.str("");
	        xiHistName.str("");

	        laHistName << "la_";
	        laHistName << mult+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << mult+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        la_mass[mult][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),31,39,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	        xiHist_mass[mult][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),31,39,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

	  
	        laHistName.str("");

	        laHistName << "genla_";
	        laHistName << mult+1;
	        laHistName << "_";
	        laHistName << pt;

	        genla_mass[mult][pt] = genlaHist[mult]->ProjectionZ( laHistName.str().c_str(),31,39,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	   	
	   	}
   
   }

   	int genksYield[2][28];
   	int genlaYield[2][20];

   	for(mult = 0; mult < 2; mult++){

   		for( pt = 0; pt < 28; pt++){

   			genksYield[mult][pt] = genks_mass[mult][pt]->GetEntries();
   		}
	   	
	   	for (int i = 0; i < 20; i++){

	   		genlaYield[mult][i] = genla_mass[mult][i]->GetEntries();

	   	}
	   		
	}

	double ksYield[2][28];
   	double laYield[2][20];

   	TCanvas* c1 = new TCanvas();
   	
   	c1->Print("la_fit.pdf[");

   	for(mult = 0; mult < 2; mult++){

   		for (pt = 0; pt < 28; pt++){

   			ksYield[mult][pt] = ks_YieldCal( ks_mass[mult][pt] );
   			
   		}

	   	for ( i = 0; i < 20; i++){
	
			la_mass[mult][i]->Add( xiHist_mass[mult][i],-1 );
	   			laYield[mult][i] = la_YieldCal( la_mass[mult][i] );
	   			c1->Print("la_fit.pdf");
	   	}
   	
   	}

   	c1->Print("la_fit.pdf]");

	TH1D* ks_eff_mult[2];
 	TH1D* la_eff_mult[2];
 	
 	stringstream ks_effName;
 	stringstream la_effName;

 	for(mult = 0; mult < 2; mult++){

		ks_effName.str("");
		ks_effName << "ks_eff_mult_";
		ks_effName << mult + 1;

		ks_eff_mult[mult] = new TH1D(ks_effName.str().c_str(), ks_effName.str().c_str(),28,ks_ptbins);

		la_effName.str("");
		la_effName << "la_eff_mult_";
		la_effName << mult + 1;

		la_eff_mult[mult] = new TH1D(la_effName.str().c_str(), la_effName.str().c_str(),20,la_ptbins);

		for (pt = 0; pt < 28; pt++){

			ks_eff_mult[mult]->SetBinContent(pt+1, ksYield[mult][pt]/genksYield[mult][pt]);
			ks_eff_mult[mult]->SetBinError(pt+1, errorCal( ksYield[mult][pt], genksYield[mult][pt] ) );
		}

		for (pt = 0; pt < 20; pt++){

			la_eff_mult[mult]->SetBinContent(pt+1, laYield[mult][pt]/genlaYield[mult][pt]);
			la_eff_mult[mult]->SetBinError(pt+1, errorCal( laYield[mult][pt], genlaYield[mult][pt] ) );
		}
	 	
	}


   	TFile f1("HIJING_2multBins_3rdEta_PbPb_v19.root","new");

 	for(mult = 0; mult < 2; mult++){
	   	
   		ks_eff_mult[mult]->Write();
   		la_eff_mult[mult]->Write();
   	}

}