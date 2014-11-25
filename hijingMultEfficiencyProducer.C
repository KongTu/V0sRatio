#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingMultEfficiencyProducer(){

//HIJING:
	/*TFile* file[8];
	file[0] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_0_35_2014.root");
	file[1] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_35_60_2014.root");
	file[2] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_60_90_2014.root");
	file[3] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_90_120_2014.root");
	file[4] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_120_150_2014.root");
	file[5] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_150_185_2014.root");
	file[6] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_185_220_2014.root");
	file[7] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov18_220plus_2014.root");
	//file[8] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov17_300_400_2014.root");
	*/

	TFile* file[4];
	file[0] = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_TH3D_Nov7_MultDependent_0_35_2014.root");
	file[1] = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_TH3D_Nov7_MultDependent_35_60_2014.root");
	file[2] = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_TH3D_Nov7_MultDependent_60_90_2014.root");
	file[3] = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_TH3D_Nov7_MultDependent_90_120_2014.root");

	TH3D* ksHist[8];
	TH3D* laHist[8];
	TH3D* xiHist[8];

	TH3D* genksHist[8];
	TH3D* genlaHist[8];

	TH1D* vertexZ[8];

	stringstream vertexName;

	for(int mult = 0; mult < 4; mult++){

		vertexName.str("");
		vertexName << "vertexDistZ_";
		vertexName << mult+1;

		//vertexZ[mult] = new TH1D(vertexName.str().c_str(),vertexName.str().c_str(),100,-15,15);

		ksHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_ks_underlying");
		laHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_la_underlying");
		genksHist[mult] = (TH3D*)file[mult]->Get("ana/genKS_underlying");
		genlaHist[mult] = (TH3D*)file[mult]->Get("ana/genLA_underlying");
		xiHist[mult] = (TH3D*)file[mult]->Get("ana/XiDaughter");
		//vertexZ[mult] = (TH1D*)file[mult]->Get("ana/vertexDistZ");

	}

	TH1D* ks_mass[8][28];
	TH1D* la_mass[8][20];
	TH1D* genks_mass[8][28];
	TH1D* genla_mass[8][20];
	TH1D* xiHist_mass[8][20];

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


	for(mult = 0; mult < 4; mult++){


    	for (int pt = 0; pt < 28; pt++){

    		ksHistName.str("");

    		ksHistName << "ks_";
    		ksHistName << mult+1;
    		ksHistName << "_";
	        ksHistName << pt;

	        ks_mass[mult][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),27,35,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	        ksHistName.str("");

         	ksHistName << "genks_";
         	ksHistName << mult+1;
    		ksHistName << "_";
	        ksHistName << pt;

	        genks_mass[mult][pt] = genksHist[mult]->ProjectionZ( ksHistName.str().c_str(),27,35,ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

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

	        la_mass[mult][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),27,35,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	        xiHist_mass[mult][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),27,35,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

	  
	        laHistName.str("");

	        laHistName << "genla_";
	        laHistName << mult+1;
	        laHistName << "_";
	        laHistName << pt;

	        genla_mass[mult][pt] = genlaHist[mult]->ProjectionZ( laHistName.str().c_str(),27,35,la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	   	
	   	}
   
   }

   	int genksYield[8][28];
   	int genlaYield[8][20];

   	for(mult = 0; mult < 4; mult++){

   		for( pt = 0; pt < 28; pt++){

   			genksYield[mult][pt] = genks_mass[mult][pt]->GetEntries();
   		}
	   	
	   	for (int i = 0; i < 20; i++){

	   		genlaYield[mult][i] = genla_mass[mult][i]->GetEntries();

	   	}
	   		
	}

/*
GEN spectra distribution;
 */

	TH1D* ksSpectra_gen[8];
	TH1D* laSpectra_gen[8];

	stringstream ksName;
	stringstream laName;

	for(mult = 0; mult < 4; mult++){

		ksName.str("");
		laName.str("");

		ksName << "ksSpectra_hijing_gen_";
		ksName << mult+1;

		laName << "laSpectra_hijing_gen_";
		laName << mult+1;

		ksSpectra_gen[mult] = new TH1D(ksName.str().c_str(),ksName.str().c_str(),28,ks_ptbins);
		laSpectra_gen[mult] = new TH1D(laName.str().c_str(),laName.str().c_str(),20,la_ptbins);

		for(pt = 0; pt < 28; pt++){

			ksSpectra_gen[mult]->SetBinContent(pt+1, genksYield[mult][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*4.8 ));
				ksSpectra_gen[mult]->SetBinError(pt+1, sqrt( genksYield[mult][pt] )/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*4.8 ));

		}

		for(pt = 0; pt < 20; pt++){

			laSpectra_gen[mult]->SetBinContent(pt+1, genlaYield[mult][pt]/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*4.8 ));
				laSpectra_gen[mult]->SetBinError(pt+1, sqrt( genlaYield[mult][pt] )/(2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*4.8));

		}
	}

   	double ksYield[8][28];
   	double laYield[8][20];

   	TCanvas* c1 = new TCanvas();
   	

   	for(mult = 0; mult < 4; mult++){

   		for (pt = 0; pt < 28; pt++){

   			ksYield[mult][pt] = ks_YieldCal( ks_mass[mult][pt] );
   			
   		}

	   	for ( i = 0; i < 20; i++){
	
			la_mass[mult][i]->Add( xiHist_mass[mult][i],-1 );
	   			laYield[mult][i] = la_YieldCal( la_mass[mult][i] );
	   	}
   	
   	}

  
/**
 * Efficiency in different mults bins:
 */

 	TH1D* ks_eff_mult[8];
 	TH1D* la_eff_mult[8];
 	
 	stringstream ks_effName;
 	stringstream la_effName;

 	for(mult = 0; mult < 4; mult++){

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


   	TFile f1("HIJING_4multBins_3rdYbins_pPb_10M_v13.root","new");

 	for(mult = 0; mult < 4; mult++){
	   	
   		ks_eff_mult[mult]->Write();
   		la_eff_mult[mult]->Write();
   	}

   	/*for(mult = 0; mult < 8; mult++){

   		vertexZ[mult]->Write();
   	}

   	for(mult = 0; mult < 8; mult++){

   		ksSpectra_gen[mult]->Write();
   		laSpectra_gen[mult]->Write();
   	}*/




}