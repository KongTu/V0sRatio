#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void histoSpectraRatio(){

    gStyle->SetErrorX(0);

/*
Getting the 3D histograms, and store in a 1D 3dimentional histogram:
 */

    TH3D* ksHist[9];
    TH3D* laHist[9];

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM1_July22_2014.root");
    ksHist[0] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[0] = (TH3D*)file1->Get("ana/InvMass_la_underlying");

    TFile* file2 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM2_July22_2014.root");
    ksHist[1] = (TH3D*)file2->Get("ana/InvMass_ks_underlying");
    laHist[1] = (TH3D*)file2->Get("ana/InvMass_la_underlying");

    TFile* file3 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM3_July22_2014.root");
    ksHist[2] = (TH3D*)file3->Get("ana/InvMass_ks_underlying");
    laHist[2] = (TH3D*)file3->Get("ana/InvMass_la_underlying");

    TFile* file4 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM4_July22_2014.root");
    ksHist[3] = (TH3D*)file4->Get("ana/InvMass_ks_underlying");
    laHist[3] = (TH3D*)file4->Get("ana/InvMass_la_underlying");

    TFile* file5 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM5_July22_2014.root");
    ksHist[4] = (TH3D*)file5->Get("ana/InvMass_ks_underlying");
    laHist[4] = (TH3D*)file5->Get("ana/InvMass_la_underlying");
 
   
   double ks_norm[5];
   double la_norm[5];

        ks_norm[0] = ksHist[0]->GetEntries();
        ks_norm[1] = ksHist[1]->GetEntries();
        ks_norm[2] = ksHist[2]->GetEntries();
        ks_norm[3] = ksHist[3]->GetEntries();
        ks_norm[4] = ksHist[4]->GetEntries();

        la_norm[0] = laHist[0]->GetEntries();
        la_norm[1] = laHist[1]->GetEntries();
        la_norm[2] = laHist[2]->GetEntries();
        la_norm[3] = laHist[3]->GetEntries();
        la_norm[4] = laHist[4]->GetEntries();


    TH1D* ks_HM[5][6][20];
    TH1D* la_HM[5][6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};

    stringstream ksHistName;
    stringstream laHistName;

    for (int mult = 0; mult < 5; mult++){

        for (int eta = 0; eta < 6; eta++){

            for (int pt = 0; pt < 20; pt++){

                ksHistName.str("");
                laHistName.str("");

                ksHistName << "ks1_";
                ksHistName << mult;
                ksHistName << "_";
                ksHistName << eta;
                ksHistName << "_";
                ksHistName << pt;

                laHistName << "la1_";
                laHistName << mult;
                laHistName << "_";
                laHistName << eta;
                laHistName << "_";
                laHistName << pt;

                ks_HM[mult][eta][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
                la_HM[mult][eta][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound[pt]+1,pTbinsBound[pt+1] );
            }
        }
    }

/*
****************************************
*/

/**
 * Getting efficiency from the table:
 */

    //TFile* t1 = new TFile("~/2014Research/Code/Jet'study/gitV0sRatio/eff_2Dtable.root");
    TFile* t1 = new TFile("~/Desktop/Efficiency2D_V0_all_smooth10.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks2Dnew");
    TH2D* hnew2 = (TH2D*)t1->Get("la2Dnew");

    double ks_eff[6][20];
    double la_eff[6][20];

    for (int i = 0;i < 6;i++){

        for (int r = 0; r < 20; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+2);
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+2);

        }
    }


/*
Start to fit all histograms to obtain the eff_corr yields:
 */

    double ks_HM_yield[5][6][20];
    double la_HM_yield[5][6][20];

    for (mult = 0; mult < 5; mult++){
        
        for (eta = 0; eta < 6; eta++){

            for (pt = 0; pt < 20; pt++){

                ks_HM_yield[mult][eta][pt] = ks_YieldCal( ks_HM[mult][eta][pt] );
                    ks_HM_yield[mult][eta][pt] = ks_HM_yield[mult][eta][pt]/ks_eff[eta][pt];
                la_HM_yield[mult][eta][pt] = la_YieldCal( la_HM[mult][eta][pt] ); 
                    la_HM_yield[mult][eta][pt] = la_HM_yield[mult][eta][pt]/la_eff[eta][pt];


            }
        }
    }
    
/*
*****************************************
 */

    double ks_HM_pTyield[5][20];
    double la_HM_pTyield[5][20];

    for (mult = 0; mult < 5; mult++){

        for (pt = 0; pt < 20; pt++){

            ks_HM_pTyield[mult][pt] = 0;
            la_HM_pTyield[mult][pt] = 0;

            for (eta = 0; eta < 6; eta++){

                ks_HM_pTyield[mult][pt] = ks_HM_yield[mult][eta][pt] + ks_HM_pTyield[mult][pt];
                la_HM_pTyield[mult][pt] = la_HM_yield[mult][eta][pt] + la_HM_pTyield[mult][pt];
            }

            cout << "ks_HM_pTyield: " << ks_HM_pTyield[mult][pt] << endl;
            cout << "la_HM_pTyield: " << la_HM_pTyield[mult][pt] << endl;

        }
    }

    double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};

    TH1D* ksSpectra[5];
    TH1D* laSpectra[5];
    TH1D* ratioHist[5];

    stringstream ratioHistName;

    for (mult = 0; mult < 5; mult++){
       
        ksHistName.str("");
        laHistName.str("");
        ratioHistName.str("");

            ksHistName << "ksSpectra_";
            ksHistName << mult+1;

            laHistName << "laSpectra_";
            laHistName << mult+1;

            ratioHistName << "ratioHist_";
            ratioHistName << mult+1;

        ksSpectra[mult] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),20,ptbins);
        laSpectra[mult] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),20,ptbins);
        ratioHist[mult] = new TH1D(ratioHistName.str().c_str(),ratioHistName.str().c_str(),20,ptbins);
        
        double etarange = 4.8;

        for (pt = 0; pt < 20; pt++){

            double ks_temp = (ks_HM_pTyield[mult][pt]/binwidth[pt])/(2*3.1415926*ptbins[pt+1]*etarange*ks_norm[mult]);
            double la_temp = (la_HM_pTyield[mult][pt]/binwidth[pt])/(2*3.1415926*ptbins[pt+1]*etarange*la_norm[mult]);

            ksSpectra[mult]->SetBinContent(pt+1, ks_temp );
            ksSpectra[mult]->SetBinError(pt+1, sqrt((ks_HM_pTyield[mult][pt]/binwidth[pt]))/(2*3.1415926*ptbins[pt+1]*etarange*ks_norm[mult]));

            laSpectra[mult]->SetBinContent(pt+1, la_temp );
            laSpectra[mult]->SetBinError(pt+1, sqrt((la_HM_pTyield[mult][pt]/binwidth[pt]))/(2*3.1415926*ptbins[pt+1]*etarange*la_norm[mult]));

            ratioHist[mult]->SetBinContent(pt+1, la_HM_pTyield[mult][pt]/(2*ks_HM_pTyield[mult][pt]));
            double err = errorCal( (la_HM_pTyield[mult][pt]/binwidth[pt]), (ks_HM_pTyield[mult][pt]/binwidth[pt]) );
            ratioHist[mult]->SetBinError(pt+1, err );
        }

    }


    TFile f1("HMbins.root","new");
    for( int you = 0; you < 5; you++){

        ksSpectra[you]->Write();
        laSpectra[you]->Write();
        ratioHist[you]->Write();
    }
        


}