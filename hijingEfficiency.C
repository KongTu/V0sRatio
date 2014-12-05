#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void hijingEfficiency(){

	TFile* file = new TFile("/Users/kongkong/2014Research/ROOT_file/3MCsamples_comparison/allEfficiency_v1.root");
	
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

	TH1D* ks_mass[6][26];
	TH1D* la_mass[6][20];
	TH1D* genks_mass[6][26];
	TH1D* genla_mass[6][20];
	TH1D* xiHist_mass[6][20];

 	double ks_pTbinsBound[27] = {1,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[27] = {0.1,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[26] = {0.2,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};
    double ks_ptbincenter[26] = {0.2,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

    double la_pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.65,0.85,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
    double la_binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double etabinsBound[] = {11,19,27,35,43,51,59};
    double etabins[] = {-2.4,-1.6,-0.8,0,0.8,1.6,2.4};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

  
    for ( int eta = 0; eta < 6; eta++){

	    for (int pt = 0; pt < 26; pt++){

	        ksHistName.str("");
	        
	        ksHistName << "ks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        ks_mass[eta][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),etabinsBound[eta]+1,etabinsBound[eta+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );
	        
	        ksHistName.str("");

	        ksHistName << "genks_";
	        ksHistName << eta+1;
	        ksHistName << "_";
	        ksHistName << pt;

	        genks_mass[eta][pt] = genksHist->ProjectionZ( ksHistName.str().c_str(),etabinsBound[eta]+1,etabinsBound[eta+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

		}

       for(pt = 0; pt < 20; pt++){


        	laHistName.str("");
	        xiHistName.str("");

	        laHistName << "la_";
	        laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        xiHistName << "xiDau_";
	        xiHistName << eta+1;
	        xiHistName << "_";
	        xiHistName << pt;

	        la_mass[eta][pt] = laHist->ProjectionZ( laHistName.str().c_str(),etabinsBound[eta]+1,etabinsBound[eta+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
        	xiHist_mass[eta][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),etabinsBound[eta]+1,etabinsBound[eta+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );

  			laHistName.str("");

	        laHistName << "genla_";
          	laHistName << eta+1;
	        laHistName << "_";
	        laHistName << pt;

	        genla_mass[eta][pt] = genlaHist->ProjectionZ( laHistName.str().c_str(),etabinsBound[eta]+1,etabinsBound[eta+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
	   	
        }
	        
   }

   	int genksYield[6][26];
   	int genlaYield[6][20];

   	for ( eta = 0; eta < 6; eta++){
	   	
	   	for (pt = 0; pt < 26; pt++){

	   		genksYield[eta][pt] = genks_mass[eta][pt]->GetEntries();
	   		
	   	}
	   	for(pt = 0; pt < 20; pt++){

	   		genlaYield[eta][pt] = genla_mass[eta][pt]->GetEntries();

	   	}
   	}	

   	double ksYield[6][26];
   	double laYield[6][20];

   	for ( eta = 0; eta < 6; eta++){

	   	for (pt = 0; pt < 26; pt++){

	   		ksYield[eta][pt] = ks_YieldCal( ks_mass[eta][pt] );
	   			
	   	}
	   	for(pt = 0; pt < 20; pt++){

	   		la_mass[eta][pt]->Add( xiHist_mass[eta][pt],-1 );
	   		laYield[eta][pt] = la_YieldCal( la_mass[eta][pt] );
	   	}
   	}


/**
 * Efficiency in different eta bins:
 */

 	TH1D* ks_eff_eta[6];
 	TH1D* la_eff_eta[6];
 	
 	stringstream ks_effName;
 	stringstream la_effName;

 	for (eta = 0; eta < 6; eta++){

 		ks_effName.str("");
 		ks_effName << "ks_eff_all_";
 		ks_effName << eta + 1;
 		ks_eff_eta[eta] = new TH1D(ks_effName.str().c_str(), ks_effName.str().c_str(),26,ks_ptbins);

 		la_effName.str("");
 		la_effName << "la_eff_all_";
 		la_effName << eta + 1;
 		la_eff_eta[eta] = new TH1D(la_effName.str().c_str(), la_effName.str().c_str(),20,la_ptbins);

 		for (pt = 0; pt < 26; pt++){
 			
 			ks_eff_eta[eta]->SetBinContent(pt+1, ksYield[eta][pt]/genksYield[eta][pt]);
 				ks_eff_eta[eta]->SetBinError(pt+1, errorCal( ksYield[eta][pt], genksYield[eta][pt] ) );

 		}

 		for(pt = 0; pt < 20; pt++){

 			la_eff_eta[eta]->SetBinContent(pt+1, laYield[eta][pt]/genlaYield[eta][pt]);
 				la_eff_eta[eta]->SetBinError(pt+1, errorCal( laYield[eta][pt], genlaYield[eta][pt] ) );
 		}
 	}

/**
 * 2D efficiency Table:
 */


   	TH2D* ks_eff = new TH2D("ks_eff","ks_eff",6,etabins,26,ks_ptbins);
   	TH2D* la_eff = new TH2D("la_eff","la_eff",6,etabins,20,la_ptbins);

   	for ( eta = 0; eta < 6; eta++){

	   	for (pt = 3; pt < 26; pt++){

	   		double temp_ks = ksYield[eta][pt]/genksYield[eta][pt];
	   		
	   		ks_eff->SetBinContent(eta+1,pt+1, temp_ks );
   			ks_eff->SetBinError(eta+1,pt+1, errorCal(ksYield[eta][pt], genksYield[eta][pt]) );
	   		
	   	}

	   	for(pt = 0; pt < 20; pt++){

	   		double temp_la = laYield[eta][pt]/genlaYield[eta][pt];
	   		
	   		la_eff->SetBinContent(eta+1,pt+1, temp_la );
   			la_eff->SetBinError(eta+1,pt+1, errorCal(laYield[eta][pt], genlaYield[eta][pt]) );
	   	}
   	}


   	TFile f1("allEfficiency_29M_v1.root","new");
   	ks_eff->Write();
   	la_eff->Write();


   	for ( eta = 0; eta < 6; eta++){

   		ks_eff_eta[eta]->Write();
   		la_eff_eta[eta]->Write();
   	}
   

}