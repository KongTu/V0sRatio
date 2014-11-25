#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void etaReweigh(){

	//HIJING:
	TFile* file[2];
	file[0] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/Hijing_0_90_2014.root");
  	file[1] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/Hijing_185plus_2014.root");

	TH3D* ksHist[2];
	TH3D* laHist[2];
	TH3D* xiHist[2];

	TH3D* genksHist[2];
	TH3D* genlaHist[2];

	TH1D* vertexZ[2];

	TH1D* eventNumber[2];

	stringstream vertexName;

	for(int mult = 0; mult < 2; mult++){

		vertexName.str("");
		vertexName << "vertexDistZ_";
		vertexName << mult+1;

		vertexZ[mult] = new TH1D(vertexName.str().c_str(),vertexName.str().c_str(),100,-15,15);

		ksHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_ks_underlying");
		laHist[mult] = (TH3D*)file[mult]->Get("ana/InvMass_la_underlying");
		genksHist[mult] = (TH3D*)file[mult]->Get("ana/genKS_underlying");
		genlaHist[mult] = (TH3D*)file[mult]->Get("ana/genLA_underlying");
		xiHist[mult] = (TH3D*)file[mult]->Get("ana/XiDaughter");
		vertexZ[mult] = (TH1D*)file[mult]->Get("ana/vertexDistZ");

		eventNumber[mult] = (TH1D*)file[mult]->Get("ana/eventNumber");

	}

	TH1D* ks_mass[2][5][28];
	TH1D* la_mass[2][5][20];
	TH1D* genks_mass[2][5][28];
	TH1D* genla_mass[2][5][20];
	TH1D* xiHist_mass[2][5][20];

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

	for(mult = 0; mult < 2; mult++){

		for(int rpy = 0; rpy < 5; rpy++){

    	for (int pt = 0; pt < 28; pt++){

    		ksHistName.str("");

    		ksHistName << "ks_";
    		ksHistName << mult+1;
    		ksHistName << "_";
    		ksHistName << rpy+1;
    		ksHistName << "_";
	        ksHistName << pt+1;

	        ks_mass[mult][rpy][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

	        ksHistName.str("");

         	ksHistName << "genks_";
         	ksHistName << mult+1;
    		ksHistName << "_";
    		ksHistName << rpy+1;
    		ksHistName << "_";
	        ksHistName << pt+1;

	        genks_mass[mult][rpy][pt] = genksHist[mult]->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

    	}

	    for (pt = 0; pt < 20; pt++){

	        laHistName.str("");
	        xiHistName.str("");

	        laHistName << "la_";
	        laHistName << mult+1;
	        laHistName << "_";
	        laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt+1;

	        xiHistName << "xiDau_";
	        xiHistName << mult+1;
	        xiHistName << "_";
	        xiHistName << rpy+1;
	        xiHistName << "_";
	        xiHistName << pt+1;

	        la_mass[mult][rpy][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	        xiHist_mass[mult][rpy][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

	        laHistName.str("");

	        laHistName << "genla_";
	        laHistName << mult+1;
	        laHistName << "_";
	        laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt+1;

	        genla_mass[mult][rpy][pt] = genlaHist[mult]->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	   	
	   	}
	   }
   
   }

   double ks_count[2][5];
   double genks_count[2][5];
   double la_count[2][5];
   double genla_count[2][5];
   double xi_count[2][5];

   for(mult = 0; mult < 2; mult++){

   		for(rpy = 0; rpy < 5; rpy++){

   			ks_count[mult][rpy] = 0.0;
   			la_count[mult][rpy] = 0.0;
   			xi_count[mult][rpy] = 0.0;

   			for(pt = 0; pt < 28; pt++){

   				ks_count[mult][rpy] = ks_count[mult][rpy] + (ks_mass[mult][rpy][pt]->GetEntries());
   				genks_count[mult][rpy] = genks_count[mult][rpy] + (genks_mass[mult][rpy][pt]->GetEntries());

   			}

   			for(pt = 0; pt < 20; pt++){

   				la_count[mult][rpy] = la_count[mult][rpy] + (la_mass[mult][rpy][pt]->GetEntries());
   				genla_count[mult][rpy] = genla_count[mult][rpy] + (genla_mass[mult][rpy][pt]->GetEntries());
   				xi_count[mult][rpy] = xi_count[mult][rpy] + (xiHist_mass[mult][rpy][pt]->GetEntries());

   			}
   		}
   }

   double ks_totcount[2];
   double la_totcount[2];
   double genks_totcount[2];
   double genla_totcount[2];
   double xi_totcount[2];

   for(mult = 0; mult < 2; mult++){

   	ks_totcount[mult] = 0.0;
	la_totcount[mult] = 0.0;
	xi_totcount[mult] = 0.0;

   		for(rpy = 0; rpy < 5; rpy++){

   			for(pt = 0; pt < 28; pt++){

   				ks_totcount[mult] = ks_totcount[mult] + (ks_mass[mult][rpy][pt]->GetEntries());
   				genks_totcount[mult] = genks_totcount[mult] + (genks_mass[mult][rpy][pt]->GetEntries());

   			}

   			for(pt = 0; pt < 20; pt++){

   				la_totcount[mult] = la_totcount[mult] + (la_mass[mult][rpy][pt]->GetEntries());
   				genla_totcount[mult] = genla_totcount[mult] + (genla_mass[mult][rpy][pt]->GetEntries());
   				xi_totcount[mult] = xi_totcount[mult] + (xiHist_mass[mult][rpy][pt]->GetEntries());

   			}
   		}
   }

   double ks_totratio = ks_totcount[1]/ks_totcount[0];
   double la_totratio = la_totcount[1]/la_totcount[0];
   double genks_totratio = genks_totcount[1]/genks_totcount[0];
   double genla_totratio = genla_totcount[1]/genla_totcount[0];
   double xi_totratio = xi_totcount[1]/xi_totcount[0];


   double ks_ratio[5];
   double genks_ratio[5];
   double la_ratio[5];
   double genla_ratio[5];
   double xi_ratio[5];

	for(rpy = 0; rpy < 5; rpy++){

		ks_ratio[rpy] = ks_count[0][rpy]/ks_count[1][rpy];
		la_ratio[rpy] = la_count[0][rpy]/la_count[1][rpy];
		xi_ratio[rpy] = xi_count[0][rpy]/xi_count[1][rpy];

		genks_ratio[rpy] = genks_count[0][rpy]/genks_count[1][rpy];
		genla_ratio[rpy] = genla_count[0][rpy]/genla_count[1][rpy];

	}

	for(rpy = 0; rpy < 5; rpy++){

		cout << " " << ks_ratio[rpy] << endl; 
		cout << " " << la_ratio[rpy] << endl; 
		cout << " " << xi_ratio[rpy] << endl;

		cout << " " << genks_ratio[rpy] << endl;
		cout << " " << genla_ratio[rpy] << endl;

	}
/*
	for(rpy = 0; rpy < 5; rpy++){

		for(pt = 0; pt < 28; pt++){

			ks_mass[1][rpy][pt]->Scale( ks_ratio[rpy] );
			genks_mass[1][rpy][pt]->Scale( genks_ratio[rpy] );
		}



		for(pt = 0; pt < 20; pt++){

			la_mass[1][rpy][pt]->Scale( la_ratio[rpy] );
			genla_mass[1][rpy][pt]->Scale( genla_ratio[rpy] );
			xiHist_mass[1][rpy][pt]->Scale( xi_ratio[rpy] );
		}

	}


	for(rpy = 0; rpy < 5; rpy++){

		for(pt = 0; pt < 28; pt++){

			ks_mass[1][rpy][pt]->Scale( ks_totratio );
			genks_mass[1][rpy][pt]->Scale( genks_totratio );
		}

			

		for(pt = 0; pt < 20; pt++){

			la_mass[1][rpy][pt]->Scale( la_totratio );
			genla_mass[1][rpy][pt]->Scale( genla_totratio );
			xiHist_mass[1][rpy][pt]->Scale( xi_totratio );
		}

	}
*/


	for(mult = 0; mult < 2; mult++){

		for(pt = 0; pt < 28; pt++){

			for(rpy = 1; rpy < 5; rpy++){

				ks_mass[mult][0][pt]->Add(ks_mass[mult][rpy][pt],+1);
				genks_mass[mult][0][pt]->Add(genks_mass[mult][rpy][pt],+1);
			}

		}

		for(pt = 0; pt < 20; pt++){

			for(rpy = 1; rpy < 5; rpy++){

				la_mass[mult][0][pt]->Add(la_mass[mult][rpy][pt],+1);
				genla_mass[mult][0][pt]->Add(genla_mass[mult][rpy][pt],+1);
				xiHist_mass[mult][0][pt]->Add(xiHist_mass[mult][rpy][pt],+1);
			}

		}

	}

	int genksYield[2][28];
   	int genlaYield[2][20];

   	for(mult = 0; mult < 2; mult++){

   		for( pt = 0; pt < 28; pt++){

   			genksYield[mult][pt] = genks_mass[mult][0][pt]->GetEntries();
   		}
	   	
	   	for (int i = 0; i < 20; i++){

	   		genlaYield[mult][i] = genla_mass[mult][0][i]->GetEntries();

	   	}
	   		
	}

	double ksYield[2][28];
   	double laYield[2][20];

   	for(mult = 0; mult < 2; mult++){

   		for (pt = 0; pt < 28; pt++){

   			ksYield[mult][pt] = ks_YieldCal( ks_mass[mult][0][pt] );
   			
   		}

	   	for ( i = 0; i < 20; i++){
	
			la_mass[mult][0][i]->Add( xiHist_mass[mult][0][i],-1 );
	   			laYield[mult][i] = la_YieldCal( la_mass[mult][0][i] );

	   		
	   	}
   	
   	}

  
/**
 * Efficiency in different mults bins:
 */

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


   	TFile f1("HIJING_2multBins_6M_v10_test.root","new");

 	for(mult = 0; mult < 2; mult++){
	   	
   		ks_eff_mult[mult]->Write();
   		la_eff_mult[mult]->Write();
   	}
	   



}