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
    //TH3D* xiHist[9];

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM1_July22_2014.root");
    ksHist[0] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[0] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    //xiHist[0] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file2 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM2_July22_2014.root");
    ksHist[1] = (TH3D*)file2->Get("ana/InvMass_ks_underlying");
    laHist[1] = (TH3D*)file2->Get("ana/InvMass_la_underlying");
    //xiHist[1] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file3 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM3_July22_2014.root");
    ksHist[2] = (TH3D*)file3->Get("ana/InvMass_ks_underlying");
    laHist[2] = (TH3D*)file3->Get("ana/InvMass_la_underlying");
    //xiHist[2] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file4 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM4_July22_2014.root");
    ksHist[3] = (TH3D*)file4->Get("ana/InvMass_ks_underlying");
    laHist[3] = (TH3D*)file4->Get("ana/InvMass_la_underlying");
    //xiHist[3] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file5 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_3Dhisto/HMpPb_HM5_July22_2014.root");
    ksHist[4] = (TH3D*)file5->Get("ana/InvMass_ks_underlying");
    laHist[4] = (TH3D*)file5->Get("ana/InvMass_la_underlying");
    //xiHist[4] = (TH3D*)file1->Get("ana/XiDaughter");
 
   
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
    TH1D* xi_HM[5][6][20];

    double pTbinsBound[21] = {6,8,10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};

    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for (int mult = 0; mult < 5; mult++){

        for (int eta = 0; eta < 6; eta++){

            for (int pt = 3; pt < 15; pt++){

                ksHistName.str("");
                laHistName.str("");
                xiHistName.str("");

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

                xiHistName << "xi1_";
                xiHistName << mult;
                xiHistName << "_";
                xiHistName << eta;
                xiHistName << "_";
                xiHistName << pt;

                ks_HM[mult][eta][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
                la_HM[mult][eta][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
                //xi_HM[mult][eta][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),eta+1,eta+1,pTbinsBound_1[pt]+1,pTbinsBound_1[pt+1] );
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
    TFile* t1 = new TFile("~/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyTable/HIJING_withXiRemoval_effKongNew2DTable_10M_Sep21_v2_12pTbins.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    double ks_eff[6][20];
    double la_eff[6][20];

    double ks_eff_err[6][20];
    double la_eff_err[6][20];

    for (int i = 0; i < 6; i++){

        for (int r = 3; r < 15; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+1);
                
                ks_eff_err[i][r] = hnew1->GetBinError(i+1,r+1);
            
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+1);
                
                la_eff_err[i][r] = hnew2->GetBinError(i+1,r+1);

        }
    }


/*
Start to fit all histograms to obtain the eff_corr yields:
 */

    double ks_HM_yield[5][6][20];
    double la_HM_yield[5][6][20];

    double ks_HM_pTyield[5][20];
    double la_HM_pTyield[5][20];

    double temp_ks_err[5][6];
    double temp_la_err[5][6];
    double num_ks_err[5][20];
    double num_la_err[5][20];


    for (mult = 0; mult < 5; mult++){
        
        for (pt = 3; pt < 15; pt++){

            ks_HM_pTyield[mult][pt] = 0;
            la_HM_pTyield[mult][pt] = 0;

            for (eta = 0; eta < 6; eta++){

                ks_HM_yield[mult][eta][pt] = ks_YieldCal( ks_HM[mult][eta][pt] );

                double ksYield_err = sqrt( ks_HM_yield[mult][eta][pt] );
                temp_ks_err[mult][eta] = errorCal_num( ks_HM_yield[mult][eta][pt], ksYield_err, ks_eff[eta][pt], ks_eff_err[eta][pt] );

                    ks_HM_yield[mult][eta][pt] = ks_HM_yield[mult][eta][pt]/ks_eff[eta][pt];

                        ks_HM_pTyield[mult][pt] = ks_HM_yield[mult][eta][pt] + ks_HM_pTyield[mult][pt];

                
                //la_HM[mult][eta][pt]->Add( xi_HM[mult][eta][pt], -1);
                la_HM_yield[mult][eta][pt] = la_YieldCal( la_HM[mult][eta][pt] ); 
                    
                double laYield_err = sqrt( la_HM_yield[mult][eta][pt] );
                temp_la_err[mult][eta] = errorCal_num( la_HM_yield[mult][eta][pt], laYield_err, la_eff[eta][pt], la_eff_err[eta][pt] );

                    la_HM_yield[mult][eta][pt] = la_HM_yield[mult][eta][pt]/la_eff[eta][pt];

                        la_HM_pTyield[mult][pt] = la_HM_yield[mult][eta][pt] + la_HM_pTyield[mult][pt];


            }

            num_ks_err[mult][pt] = errorCal_sum( temp_ks_err[mult][0], temp_ks_err[mult][1], temp_ks_err[mult][2], temp_ks_err[mult][3], temp_ks_err[mult][4], temp_ks_err[mult][5] );
            num_la_err[mult][pt] = errorCal_sum( temp_la_err[mult][0], temp_la_err[mult][1], temp_la_err[mult][2], temp_la_err[mult][3], temp_la_err[mult][4], temp_la_err[mult][5] );
        }
    }


    
/*
*****************************************
 */


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

        ksSpectra[mult] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),15,ptbins_1);
        laSpectra[mult] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),15,ptbins_1);
        ratioHist[mult] = new TH1D(ratioHistName.str().c_str(),ratioHistName.str().c_str(),15,ptbins_1);
        
        double etarange = 4.8;

        for (pt = 3; pt < 15; pt++){

            double ks_temp = (ks_HM_pTyield[mult][pt]/binwidth_1[pt])/(2*3.1415926*ptbins_1[pt+1]*etarange*ks_norm[mult]);
            double la_temp = (la_HM_pTyield[mult][pt]/binwidth_1[pt])/(2*3.1415926*ptbins_1[pt+1]*etarange*la_norm[mult]);

            ksSpectra[mult]->SetBinContent(pt+1, ks_temp );
            ksSpectra[mult]->SetBinError(pt+1, num_ks_err[mult][pt]/( 2*3.1415926*binwidth_1[pt]*ptbins_1[pt+1]*etarange*ks_norm[mult] ));

            laSpectra[mult]->SetBinContent(pt+1, la_temp );
            laSpectra[mult]->SetBinError(pt+1, num_la_err[mult][pt]/( 2*3.1415926*binwidth_1[pt]*ptbins_1[pt+1]*etarange*la_norm[mult] ));

            ratioHist[mult]->SetBinContent(pt+1, la_HM_pTyield[mult][pt]/(2*ks_HM_pTyield[mult][pt]));

                double err = errorCal_lambdakshort( la_HM_pTyield[mult][pt], ks_HM_pTyield[mult][pt], num_la_err[mult][pt], num_ks_err[mult][pt] );
            
                ratioHist[mult]->SetBinError(pt+1, err );
        }

    }


    TFile f1("HMbins_12pTBins_wErr.root","new");
    for( int you = 0; you < 5; you++){

        ksSpectra[you]->Write();
        laSpectra[you]->Write();
        ratioHist[you]->Write();
    }
        


}