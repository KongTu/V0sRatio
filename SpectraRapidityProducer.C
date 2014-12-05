#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void SpectraRapidityProducer(){

	 gStyle->SetErrorX(0);

/*
Getting the 3D histograms, and store in a 1D 3dimentional histogram:
 */

    TH3D* ksHist;
    TH3D* laHist;
    TH3D* xiHist;

    TH1D* eventNumber;

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_rpyDependence_sample/HM_Nov18_185_300_2014.root");
    ksHist = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist = (TH3D*)file1->Get("ana/XiDaughter");

    eventNumber = (TH1D*)file1->Get("ana/eventNumber");

    double Nev = eventNumber->GetBinContent(2);
    
    TH1D* ks_HM[5][26];
    TH1D* la_HM[5][20];
    TH1D* xi_HM[5][20];

    double ks_pTbinsBound[27] = {0,3,4,5,6,7,8,9,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double ks_ptbins[27] = {0,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double ks_binwidth[26] = {0.3,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};
    double ks_ptbincenter[26] = {0.15,0.35,0.45,0.55,0.65,0.75,0.85,0.95,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};

    double la_pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double la_ptbins[21] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double la_ptbincenter[20] = {0.65,0.85,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.2,3.6,4.0,4.4,4.8,5.3,6.1,7.8};
    double la_binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,2.4};

    double rpybins[6] = {6,16,26,35,44,55};
    double rpybinwidth[5] = {1.0,0.97,0.9,0.93,1.0};

 	stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for (int rpy = 0; rpy < 5; rpy++){

        for (int pt = 0; pt < 26; pt++){
            
            ksHistName.str("");

            ksHistName << "ks1_";
            ksHistName << rpy;
            ksHistName << "_";
            ksHistName << pt;

            ks_HM[rpy][pt] = ksHist->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

        }

        for ( pt = 0; pt < 20; pt++){

            laHistName.str("");
            xiHistName.str("");

            laHistName << "la1_";
            laHistName << rpy;
            laHistName << "_";
            laHistName << pt;

            xiHistName << "xi1_";
            xiHistName << rpy;
            xiHistName << "_";
            xiHistName << pt;

            
            la_HM[rpy][pt] = laHist->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
            xi_HM[rpy][pt] = xiHist->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
        }
    }

/*
EPOS efficiency with vertex reweigh. 
 */

    TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/eposEfficiencyRapidityTable/EPOSeffTable_withVTXreweight_5M_Nov14_rapidity_v3_26ks_pTbins.root");

    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    double ks_eff[5][26];
    double la_eff[5][20];

    double ks_eff_err[5][26];
    double la_eff_err[5][20];

    for (rpy = 0; rpy < 5; rpy++){
        
        for (pt = 0; pt < 26; pt++){
            
            ks_eff[rpy][pt] = hnew1->GetBinContent(rpy+1,pt+1);
                ks_eff_err[rpy][pt] = hnew1->GetBinError(rpy+1,pt+1);
        }

        for (pt = 0; pt < 20; pt++){

            la_eff[rpy][pt] = hnew2->GetBinContent(rpy+1,pt+1);
                la_eff_err[rpy][pt] = hnew2->GetBinError(rpy+1,pt+1);
        }
    }


/*
Fitting all histograms to obtain the raw yield and then apply efficiency, calculate the uncertainties
 */

	double ks_HM_yield[5][26];
    double la_HM_yield[5][20];

    double temp_ks_err[5][26];
    double temp_la_err[5][20];

	for (rpy = 0; rpy < 5; rpy++){

		for(pt = 0; pt < 26; pt++){

			ks_HM_yield[rpy][pt] = ks_YieldCal( ks_HM[rpy][pt] );
				double ksYield_err = sqrt( ks_HM_yield[rpy][pt] );

				temp_ks_err[rpy][pt] = errorCal_num(ks_HM_yield[rpy][pt], ksYield_err, ks_eff[rpy][pt], ks_eff_err[rpy][pt]);

				ks_HM_yield[rpy][pt] = ks_HM_yield[rpy][pt]/ks_eff[rpy][pt];

		}

		for(pt = 0; pt <20; pt++){

			la_HM[rpy][pt]->Add( xi_HM[rpy][pt], -1);

			la_HM_yield[rpy][pt] = la_YieldCal( la_HM[rpy][pt] );
				double laYield_err = sqrt( la_HM_yield[rpy][pt] );

				temp_la_err[rpy][pt] = errorCal_num(la_HM_yield[rpy][pt], laYield_err, la_eff[rpy][pt], la_eff_err[rpy][pt]);

				la_HM_yield[rpy][pt] = la_HM_yield[rpy][pt]/la_eff[rpy][pt];
		}


	}

/*
Filling the yield into the spectra and assign uncertainties:
 */
    
	TH1D* ksSpectra[5];
	TH1D* laSpectra[5];

	for(rpy = 0; rpy < 5; rpy++){

		ksHistName.str("");
        laHistName.str("");

        ksHistName << "ksSpectra_";
        ksHistName << rpy+1;
        
        laHistName << "laSpectra_";
        laHistName << rpy+1;
      
        ksSpectra[rpy] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),26,ks_ptbins);
        laSpectra[rpy] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),20,la_ptbins);



		for (pt = 0; pt < 26; pt++){

			double ks_temp = (ks_HM_yield[rpy][pt]/ks_binwidth[pt])/(2*3.1415926*ks_ptbincenter[pt]*rpybinwidth[rpy]*Nev);

			ksSpectra[rpy]->SetBinContent(pt+1, ks_temp );
			ksSpectra[rpy]->SetBinError(pt+1, temp_ks_err[rpy][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*rpybinwidth[rpy]*Nev ));

		}

		for (pt = 0; pt < 20; pt++){

			double la_temp = (la_HM_yield[rpy][pt]/la_binwidth[pt])/(2*3.1415926*la_ptbincenter[pt]*rpybinwidth[rpy]*Nev);

			laSpectra[rpy]->SetBinContent(pt+1, la_temp );
			laSpectra[rpy]->SetBinError(pt+1, temp_la_err[rpy][pt]/( 2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*rpybinwidth[rpy]*Nev ));

		}
	}

	TFile f1("SpectraRapidityProducer_185_300_vtxReweighEPOS_26ksPTbins_v6.root","new");
    
    for( rpy = 0; rpy < 5; rpy++){

        ksSpectra[rpy]->Write();
        laSpectra[rpy]->Write();
       
    }
    
}