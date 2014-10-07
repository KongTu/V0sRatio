#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void histoSpectraRatioRapidity(){

    gStyle->SetErrorX(0);

/*
Getting the 3D histograms, and store in a 1D 3dimentional histogram:
 */

    TH3D* ksHist[5];
    TH3D* laHist[5];
    TH3D* xiHist[5];

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HM_TH3D_Oct3_rpyDependent_120_150_2014.root");
    ksHist[0] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[0] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[0] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HM_TH3D_Oct3_rpyDependent_150_185_2014.root");
    ksHist[1] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[1] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[1] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HM_TH3D_Oct3_rpyDependent_185_220_2014.root");
    ksHist[2] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[2] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[2] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HM_TH3D_Oct3_rpyDependent_220_260_2014.root");
    ksHist[3] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[3] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[3] = (TH3D*)file1->Get("ana/XiDaughter");

    TFile* file1 = new TFile("~/2014Research/ROOT_file/V0reco_pPb_rpyDependent_3Dhisto/HM_TH3D_Oct3_rpyDependent_260plus_2014.root");
    ksHist[4] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[4] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[4] = (TH3D*)file1->Get("ana/XiDaughter");

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
        

    TH1D* ks_HM[5][5][20];
    TH1D* la_HM[5][5][20];
    TH1D* xi_HM[5][5][20];

    double pTbinsBound[19] = {10,12,14,16,18,20,22,24,26,28,30,34,38,42,46,50,56,66,90};
    double pTbinsBound_1[16] = {0,2,4,6,8,10,12,14,16,18,20,26,32,42,60,90};

    double ptbins_1[] = {0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
    double binwidth_1[15] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.6,0.6,1.0,1.8,3.0};

    double rpybins[6] = {6,16,26,35,44,55};
    //double rpybins[6] = {33,91,145,198,250,314};

    stringstream ksHistName;
    stringstream laHistName;
    stringstream xiHistName;

    for (int mult = 0; mult < 5; mult++){

        for (int rpy = 0; rpy < 5; rpy++){

            for (int pt = 0; pt < 18; pt++){

                ksHistName.str("");
                laHistName.str("");
                xiHistName.str("");

                ksHistName << "ks1_";
                ksHistName << mult;
                ksHistName << "_";
                ksHistName << rpy;
                ksHistName << "_";
                ksHistName << pt;

                laHistName << "la1_";
                laHistName << mult;
                laHistName << "_";
                laHistName << rpy;
                laHistName << "_";
                laHistName << pt;

                xiHistName << "xi1_";
                xiHistName << mult;
                xiHistName << "_";
                xiHistName << rpy;
                xiHistName << "_";
                xiHistName << pt;

                ks_HM[mult][rpy][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );
                la_HM[mult][rpy][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );
                xi_HM[mult][rpy][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],pTbinsBound[pt]+1,pTbinsBound[pt+1] );
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
    TFile* t1 = new TFile("~/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Oct2_rapidity_v1_18pTbins.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    TH1D* ks_eff_hist[5];
    TH1D* la_eff_hist[5];

    double ks_eff[5][20];
    double la_eff[5][20];

    double ks_eff_err[5][20];
    double la_eff_err[5][20];

    stringstream ksEffName;
    stringstream laEffName;

    for(rpy = 0; rpy < 5; rpy++){

        ksEffName.str("");
        laEffName.str("");

        ksEffName << "kshist_";
        ksEffName << rpy+1;

        laEffName << "lahist_";
        laEffName << rpy+1;

        ks_eff_hist[rpy] = hnew1->ProjectionY(ksEffName.str().c_str(), rpy+1,rpy+1 );
            ks_eff_hist[rpy]->Smooth(10);
        la_eff_hist[rpy] = hnew2->ProjectionY(laEffName.str().c_str(), rpy+1,rpy+1 );
            la_eff_hist[rpy]->Smooth(10);

    }

/*
smoothing the efficiency:
 */

    for(rpy = 0; rpy < 5; rpy++){

        for(pt = 0; pt < 18; pt++){

            ks_eff[rpy][pt] = ks_eff_hist[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_hist[rpy]->GetBinError( pt+1 );

            la_eff[rpy][pt] = la_eff_hist[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_hist[rpy]->GetBinError( pt+1 );  
        }
    }



/*
without smoothing efficiency:
 */
    /*for (int i = 0; i < 5; i++){

        for (int r = 3; r < 15; r++){

            ks_eff[i][r] = hnew1->GetBinContent(i+1,r+1);
                
                ks_eff_err[i][r] = hnew1->GetBinError(i+1,r+1);
            
            la_eff[i][r] = hnew2->GetBinContent(i+1,r+1);
                
                la_eff_err[i][r] = hnew2->GetBinError(i+1,r+1);

        }
    }*/


/*
Start to fit all histograms to obtain the eff_corr yields:
 */

    double ks_HM_yield[5][5][20];
    double la_HM_yield[5][5][20];

    double temp_ks_err[5][5][20];
    double temp_la_err[5][5][20];


    for (mult = 0; mult < 5; mult++){
        
        for (pt = 0; pt < 18; pt++){

            for (rpy = 0; rpy < 5; rpy++){

                ks_HM_yield[mult][rpy][pt] = ks_YieldCal( ks_HM[mult][rpy][pt] );

                double ksYield_err = sqrt( ks_HM_yield[mult][rpy][pt] );
                
                    /*
                    This is the error on the yield for each mult, each rpy and each pT:
                     */               
                    temp_ks_err[mult][rpy][pt] = errorCal_num( ks_HM_yield[mult][rpy][pt], ksYield_err, ks_eff[rpy][pt], ks_eff_err[rpy][pt] );

                        ks_HM_yield[mult][rpy][pt] = ks_HM_yield[mult][rpy][pt]/ks_eff[rpy][pt];

                
                la_HM[mult][rpy][pt]->Add( xi_HM[mult][rpy][pt], -1);
                    la_HM_yield[mult][rpy][pt] = la_YieldCal( la_HM[mult][rpy][pt] );
                    
                double laYield_err = sqrt( la_HM_yield[mult][rpy][pt] );
                    temp_la_err[mult][rpy][pt] = errorCal_num( la_HM_yield[mult][rpy][pt], laYield_err, la_eff[rpy][pt], la_eff_err[rpy][pt] );

                        la_HM_yield[mult][rpy][pt] = la_HM_yield[mult][rpy][pt]/la_eff[rpy][pt];

            }


        }
    }

/*
*****************************************
 */

    double ptbins[] = {1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double binwidth[18] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};

    double rpybinwidth[5] = {1.0,0.97,0.9,0.93,1.0};

/*
This is transformation from pT to mT:
 */

    double mtm0bins_ks[19];
    double mtm0bins_la[19];

    double mtm0bins_ks_width[18];
    double mtm0bins_la_width[18];

    for(pt = 0; pt < 19; pt++){

        mtm0bins_ks[pt] = sqrt(ptbins[pt] * ptbins[pt] + 0.497 * 0.497);
        mtm0bins_la[pt] = sqrt(ptbins[pt] * ptbins[pt] + 1.115 * 1.115);
    }

    for(pt = 0; pt < 18; pt++){

        mtm0bins_ks_width[pt] = mtm0bins_ks[pt+1] - mtm0bins_ks[pt];
        mtm0bins_la_width[pt] = mtm0bins_la[pt+1] - mtm0bins_la[pt];

    }
/*
**************************************
 */

    TH1D* ksSpectra[5][5];
    TH1D* laSpectra[5][5];
    TH1D* ratioHist[5][5];

    stringstream ratioHistName;

    for (mult = 0; mult < 5; mult++){

        for (rpy = 0; rpy < 5; rpy++){
       
        ksHistName.str("");
        laHistName.str("");
        ratioHistName.str("");

            ksHistName << "ksSpectra_";
            ksHistName << mult+1;
            ksHistName << "_";
            ksHistName << rpy+1;

            laHistName << "laSpectra_";
            laHistName << mult+1;
            laHistName << "_";
            laHistName << rpy+1;

            ratioHistName << "ratioHist_";
            ratioHistName << mult+1;
            ratioHistName << "_";
            ratioHistName << rpy+1;

        ksSpectra[mult][rpy] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),18,mtm0bins_ks);
        laSpectra[mult][rpy] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),18,mtm0bins_la);
        ratioHist[mult][rpy] = new TH1D(ratioHistName.str().c_str(),ratioHistName.str().c_str(),18,ptbins);
    
            for (pt = 0; pt < 18; pt++){

                double ks_temp = (ks_HM_yield[mult][rpy][pt]/mtm0bins_ks_width[pt])/(2*3.1415926*ptbins[pt+1]*rpybinwidth[rpy]*ks_norm[mult]);
                double la_temp = (la_HM_yield[mult][rpy][pt]/mtm0bins_la_width[pt])/(2*3.1415926*ptbins[pt+1]*rpybinwidth[rpy]*la_norm[mult]);

                ksSpectra[mult][rpy]->SetBinContent(pt+1, ks_temp );
                ksSpectra[mult][rpy]->SetBinError(pt+1, temp_ks_err[mult][rpy][pt]/( 2*3.1415926*mtm0bins_ks_width[pt]*ptbins[pt+1]*rpybinwidth[rpy]*ks_norm[mult] ));

                laSpectra[mult][rpy]->SetBinContent(pt+1, la_temp );
                laSpectra[mult][rpy]->SetBinError(pt+1, temp_la_err[mult][rpy][pt]/( 2*3.1415926*mtm0bins_la_width[pt]*ptbins[pt+1]*rpybinwidth[rpy]*la_norm[mult] ));

                ratioHist[mult][rpy]->SetBinContent(pt+1, la_HM_yield[mult][rpy][pt]/(2*ks_HM_yield[mult][rpy][pt]));

                double temp_la_withoutEff_err = sqrt(la_HM_yield[mult][rpy][pt]);
                double temp_ks_withoutEff_err = sqrt(ks_HM_yield[mult][rpy][pt]);

                    double err = errorCal_lambdakshort( la_HM_yield[mult][rpy][pt], ks_HM_yield[mult][rpy][pt], temp_la_err[mult][rpy][pt], temp_ks_err[mult][rpy][pt] );
                
                    //double err = errorCal_lambdakshort( la_HM_yield[mult][rpy][pt], ks_HM_yield[mult][rpy][pt], temp_la_withoutEff_err, temp_ks_withoutEff_err );
                    ratioHist[mult][rpy]->SetBinError(pt+1, err );
            }
        
        }
    }


    TFile f1("HMbins_18pTBins5rpyBins_wSmooth10Eff_FullStats_5Multbins_v1.root","new");
    
    for( mult = 0; mult < 5; mult++){

        for (rpy = 0; rpy < 5; rpy++){

            ksSpectra[mult][rpy]->Write();
            laSpectra[mult][rpy]->Write();
            ratioHist[mult][rpy]->Write();
        }
    }
        


}