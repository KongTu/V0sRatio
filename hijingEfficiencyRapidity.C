#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiencyRapidity(){

//HIJING:
	//TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_HIJING_3Dhisto/HIJING_TH3D_Oct2_rpyDependent_18M_2014.root");
	//TFile* file = new TFile("/Users/kongkong/2014Research/ROOT_file/newHIJINGsample/HIJING_18M_TH3D_Nov5_rpyDependent_vtxReweight_2014.root");
//EPOS:
	//TFile* file = new TFile("~/2014Research/ROOT_file/V0reco_pPb_EPOS_3Dhisto/EPOS_TH3D_Nov3_rpyDependent_2014.root");
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

	TH1D* ks_mass[5][26];
	TH1D* la_mass[5][20];
	TH1D* genks_mass[5][26];
	TH1D* genla_mass[5][20];
	TH1D* xiHist_mass[5][20];

    double pTbinsBound[] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};

    double ks_pTbinsBound[27] = {0,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[27] = {0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[26] = {0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double rpybins[6] = {6,16,26,35,44,55};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for ( int rpy = 0; rpy < 5; rpy++){

    	for (int pt = 0; pt < 26; pt++){

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

	        la_mass[rpy][pt] = laHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );
	        xiHist_mass[rpy][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );

	  
	        laHistName.str("");

	        laHistName << "genla_";
          	laHistName << rpy+1;
	        laHistName << "_";
	        laHistName << pt;

	        genla_mass[rpy][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );
	   	
	   	}
   }

   	int genksYield[5][26];
   	int genlaYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){


   		for( pt = 0; pt < 26; pt++){

   			genksYield[rpy][pt] = genks_mass[rpy][pt]->GetEntries();
   		}
	   	
	   	for (int i = 0; i < 20; i++){

	   		genlaYield[rpy][i] = genla_mass[rpy][i]->GetEntries();

	   	}
   	}	

   	double ksYield[5][26];
   	double laYield[5][20];

   	for ( rpy = 0; rpy < 5; rpy++){

   		for (pt = 0; pt < 26; pt++){

   			ksYield[rpy][pt] = ks_YieldCal( ks_mass[rpy][pt] );

   		}

	   	for ( i = 0; i < 20; i++){
	
			la_mass[rpy][i]->Add( xiHist_mass[rpy][i],-1 );
	   			laYield[rpy][i] = la_YieldCal( la_mass[rpy][i] );
	   	}
   	}

/**
 * Efficiency in differentrpy bins:
 */

 	TH1D* ks_eff_rpy[5];
 	TH1D* la_eff_rpy[5];
 	
 	stringstream ks_effName;
 	stringstream la_effName;

 	for (rpy = 0; rpy < 5; rpy++){

 		ks_effName.str("");
 		ks_effName << "ks_eff_rpy_vtx_";
 		ks_effName << rpy + 1;
 		ks_eff_rpy[rpy] = new TH1D(ks_effName.str().c_str(), ks_effName.str().c_str(),26,ks_ptbins);

 		la_effName.str("");
 		la_effName << "la_eff_rpy_vtx_";
 		la_effName << rpy + 1;
 		la_eff_rpy[rpy] = new TH1D(la_effName.str().c_str(), la_effName.str().c_str(),20,ptbins);

 		for (pt = 0; pt < 26; pt++){
 			
 			ks_eff_rpy[rpy]->SetBinContent(pt+1, ksYield[rpy][pt]/genksYield[rpy][pt]);
 				ks_eff_rpy[rpy]->SetBinError(pt+1, errorCal( ksYield[rpy][pt], genksYield[rpy][pt] ) );
		}

		for (pt = 0; pt < 20; pt++){

 			la_eff_rpy[rpy]->SetBinContent(pt+1, laYield[rpy][pt]/genlaYield[rpy][pt]);
 				la_eff_rpy[rpy]->SetBinError(pt+1, errorCal( laYield[rpy][pt], genlaYield[rpy][pt] ) );
 		}
 	}

/**
 * 2D efficiency Table:
 */

   	TH2D* ks_eff = new TH2D("ks_eff","ks_eff",5,rpybins,26,ks_ptbins);
   	TH2D* la_eff = new TH2D("la_eff","la_eff",5,rpybins,20,ptbins);

   	for ( rpy = 0; rpy < 5; rpy++){

   		for( pt = 0; pt < 26; pt++){

   			double temp_ks = ksYield[rpy][pt]/genksYield[rpy][pt];

   			ks_eff->SetBinContent(rpy+1,pt+1, temp_ks );
	   			ks_eff->SetBinError(rpy+1,pt+1, errorCal(ksYield[rpy][pt], genksYield[rpy][pt]) );


   		}

	   	for ( i = 0; i < 20; i++){
	   		
	   		double temp_la = laYield[rpy][i]/genlaYield[rpy][i];

	   		la_eff->SetBinContent(rpy+1,i+1, temp_la );
	   			la_eff->SetBinError(rpy+1,i+1, errorCal(laYield[rpy][i], genlaYield[rpy][i]) );


	   	}
   	}

   	TFile f1("EPOSeffTable_withVTXreweight_5M_Nov14_rapidity_v3_26ks_pTbins.root","new");
   	ks_eff->Write();
   	la_eff->Write();
 

   	for ( rpy = 0; rpy < 5; rpy++){

   		ks_eff_rpy[rpy]->Write();
   		la_eff_rpy[rpy]->Write();
   	}



}