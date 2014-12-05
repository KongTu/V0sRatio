#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void histoSpectraRatioRapidity(){

    gStyle->SetErrorX(0);

/*
Getting the 3D histograms, and store in a 1D 3dimentional histogram:
 */

    TH3D* ksHist[8];
    TH3D* laHist[8];
    TH3D* xiHist[8];
    TH1D* eventHist[8];

/**
 * 4 MB datasets:
 */

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov10_0_35.root");
    ksHist[0] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[0] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[0] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[0] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov10_35_60.root");
    ksHist[1] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[1] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[1] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[1] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov10_60_90.root");
    ksHist[2] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[2] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[2] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[2] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov10_90_120.root");
    ksHist[3] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[3] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[3] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[3] = (TH1D*)file1->Get("ana/eventNumber");

/**
 * 4 HM datasets:
 */

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/HM_Nov10_120_150.root");
    ksHist[4] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[4] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[4] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[4] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/HM_Nov10_150_185.root");
    ksHist[5] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[5] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[5] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[5] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/HM_Nov10_185_220.root");
    ksHist[6] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[6] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[6] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[6] = (TH1D*)file1->Get("ana/eventNumber");

    TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/HM_Nov10_220plus.root");
    ksHist[7] = (TH3D*)file1->Get("ana/InvMass_ks_underlying");
    laHist[7] = (TH3D*)file1->Get("ana/InvMass_la_underlying");
    xiHist[7] = (TH3D*)file1->Get("ana/XiDaughter");

    eventHist[7] = (TH1D*)file1->Get("ana/eventNumber");

    double norm[8];

    norm[0] = eventHist[0]->GetEntries();
    norm[1] = eventHist[1]->GetEntries();
    norm[2] = eventHist[2]->GetEntries();
    norm[3] = eventHist[3]->GetEntries();
    norm[4] = eventHist[4]->GetEntries();
    norm[5] = eventHist[5]->GetEntries();
    norm[6] = eventHist[6]->GetEntries();
    norm[7] = eventHist[7]->GetEntries();
        
    TH1D* ks_HM[9][5][26];
    TH1D* la_HM[9][5][20];
    TH1D* xi_HM[9][5][20];

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

    for (int mult = 0; mult < 8; mult++){

        for (int rpy = 0; rpy < 5; rpy++){

            for (int pt = 0; pt < 26; pt++){
                
                ksHistName.str("");

                ksHistName << "ks1_";
                ksHistName << mult;
                ksHistName << "_";
                ksHistName << rpy;
                ksHistName << "_";
                ksHistName << pt;

                ks_HM[mult][rpy][pt] = ksHist[mult]->ProjectionZ( ksHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],ks_pTbinsBound[pt]+1,ks_pTbinsBound[pt+1] );

            }

            for ( pt = 0; pt < 20; pt++){

                laHistName.str("");
                xiHistName.str("");

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

                
                la_HM[mult][rpy][pt] = laHist[mult]->ProjectionZ( laHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
                xi_HM[mult][rpy][pt] = xiHist[mult]->ProjectionZ( xiHistName.str().c_str(),rpybins[rpy]+1,rpybins[rpy+1],la_pTbinsBound[pt]+1,la_pTbinsBound[pt+1] );
            }
        }
    }

/*
****************************************
*/

/**
 * Getting efficiency from the table:
 */


/*  TFile* t1 = new TFile("~/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Oct2_rapidity_v1_28ks_pTbins.root");
    
    TH2D* hnew1 = (TH2D*)t1->Get("ks_eff");
    TH2D* hnew2 = (TH2D*)t1->Get("la_eff");

    TH1D* ks_eff_hist[5];
    TH1D* la_eff_hist[5];

    double ks_eff[5][28];
    double la_eff[5][20];

    double ks_eff_err[5][28];
    double la_eff_err[5][20];
*/
/*
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
/*
    for(rpy = 0; rpy < 5; rpy++){

        for(pt = 0; pt < 28; pt++){

            ks_eff[rpy][pt] = ks_eff_hist[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_hist[rpy]->GetBinError( pt+1 );
        }

        for(pt = 0; pt < 20; pt++){

            la_eff[rpy][pt] = la_eff_hist[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_hist[rpy]->GetBinError( pt+1 );  
        }
    }

*/

/*
without smoothing efficiency:
 */

    //TFile* t1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
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
with partially smoothing efficiency, only smoothing after 2GeV:
 */
/*
    TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/hijingEfficiencyRapidityTable/effKongNew2DTable_18M_Nov3_rapidity_v2_28ks_pTbins.root");
    TFile* file1 = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/EffSmoothProducer/EffSmoothProducer_v2.root");

    TH1D* ks_eff_rpy[5];
    TH1D* la_eff_rpy[5];

    TH1D* ks_eff_rpy_new[5];
    TH1D* la_eff_rpy_new[5];

    double ks_eff[5][28];
    double la_eff[5][20];

    double ks_eff_err[5][28];
    double la_eff_err[5][20];

    stringstream ksName;
    stringstream laName;

    for (int rpy = 0; rpy < 5; rpy++){

        ksName.str("");
        laName.str("");

        ksName << "ks_eff_rpy_";
        ksName << rpy+1;

        laName << "la_eff_rpy_";
        laName << rpy+1;

        ks_eff_rpy[rpy] = (TH1D*) file->Get( ksName.str().c_str() );
        la_eff_rpy[rpy] = (TH1D*) file->Get( laName.str().c_str() );

        ksName.str("");
        laName.str("");

        ksName << "ks_eff_rpy_new_";
        ksName << rpy+1;

        laName << "la_eff_rpy_new_";
        laName << rpy+1;

        ks_eff_rpy_new[rpy] = (TH1D*) file1->Get( ksName.str().c_str() );
        la_eff_rpy_new[rpy] = (TH1D*) file1->Get( laName.str().c_str() );

    }

    for(rpy = 0; rpy < 5; rpy++){

        for(pt = 0; pt < 15; pt++){

            ks_eff[rpy][pt] = ks_eff_rpy[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_rpy[rpy]->GetBinError( pt+1 );
        }

        for(pt = 0; pt < 7; pt++){

            la_eff[rpy][pt] = la_eff_rpy[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_rpy[rpy]->GetBinError( pt+1 );  
        }

      //adding smoothing after 2GeV:
        
        for(pt = 15; pt < 28; pt++){

            ks_eff[rpy][pt] = ks_eff_rpy_new[rpy]->GetBinContent( pt+1 );

                ks_eff_err[rpy][pt] = ks_eff_rpy_new[rpy]->GetBinError( pt+1 );
        }

        for(pt = 7; pt < 20; pt++){

            la_eff[rpy][pt] = la_eff_rpy_new[rpy]->GetBinContent( pt+1 );

                la_eff_err[rpy][pt] = la_eff_rpy_new[rpy]->GetBinError( pt+1 );  
        }


    }

*/
/*
Start to fit all histograms to obtain the eff_corr yields:
 */

    double ks_HM_yield[9][5][26];
    double la_HM_yield[9][5][20];

    double ks_HM_pTyield[9][26];
    double la_HM_pTyield[9][20];

    double temp_ks_err[9][5];
    double temp_la_err[9][5];

    double num_ks_err[9][26];
    double num_la_err[9][20];

    for(mult = 0; mult < 8; mult++){

        for(pt = 0; pt < 26; pt++){

            ks_HM_pTyield[mult][pt] = 0.0;
        }

        for(pt = 0; pt < 20; pt++){

            la_HM_pTyield[mult][pt] = 0.0;
        }
    }

    for (pt = 0; pt < 26; pt++){
    
        for (mult = 0; mult < 8; mult++){

            for (rpy = 0; rpy < 5; rpy++){

                ks_HM_yield[mult][rpy][pt] = ks_YieldCal( ks_HM[mult][rpy][pt] );

                double ksYield_err = sqrt( ks_HM_yield[mult][rpy][pt] );
                
                    /*
                    This is the error on the yield for each mult, each rpy and each pT:
                     */               
                    temp_ks_err[mult][rpy] = errorCal_num( ks_HM_yield[mult][rpy][pt], ksYield_err, ks_eff[rpy][pt], ks_eff_err[rpy][pt] );

                        ks_HM_yield[mult][rpy][pt] = ks_HM_yield[mult][rpy][pt]/ks_eff[rpy][pt];

                            //summing over all the rpy bins:
                                ks_HM_pTyield[mult][pt] = ks_HM_yield[mult][rpy][pt] + ks_HM_pTyield[mult][pt];

            }

            num_ks_err[mult][pt] = errorCal_sum(temp_ks_err[mult][0],temp_ks_err[mult][1],temp_ks_err[mult][2],temp_ks_err[mult][3],temp_ks_err[mult][4],0.0);
            
        }
    }

    for (pt = 0; pt < 20; pt++){
    
        for (mult = 0; mult < 8; mult++){

            for (rpy = 0; rpy < 5; rpy++){

                
                la_HM[mult][rpy][pt]->Add( xi_HM[mult][rpy][pt], -1);
                    la_HM_yield[mult][rpy][pt] = la_YieldCal( la_HM[mult][rpy][pt] );
                    
                double laYield_err = sqrt( la_HM_yield[mult][rpy][pt] );
                    temp_la_err[mult][rpy] = errorCal_num( la_HM_yield[mult][rpy][pt], laYield_err, la_eff[rpy][pt], la_eff_err[rpy][pt] );

                        la_HM_yield[mult][rpy][pt] = la_HM_yield[mult][rpy][pt]/la_eff[rpy][pt];

                        //summing over all the rpy bins:
                            la_HM_pTyield[mult][pt] = la_HM_yield[mult][rpy][pt] + la_HM_pTyield[mult][pt];

            }

            num_la_err[mult][pt] = errorCal_sum(temp_la_err[mult][0],temp_la_err[mult][1],temp_la_err[mult][2],temp_la_err[mult][3],temp_la_err[mult][4],0.0);

        }
    }

    
/*
**************************************
 */

    TH1D* ksSpectra[8];
    TH1D* laSpectra[8];
    TH1D* ratioHist[8];

    stringstream ratioHistName;

    for (mult = 0; mult < 8; mult++){
       
        ksHistName.str("");
        laHistName.str("");
        ratioHistName.str("");

        ksHistName << "ksSpectra_vtx_";
        ksHistName << mult+1;
        
        laHistName << "laSpectra_vtx_";
        laHistName << mult+1;

        ratioHistName << "ratio_vtx_";
        ratioHistName << mult+1;
      
        ksSpectra[mult] = new TH1D(ksHistName.str().c_str(),ksHistName.str().c_str(),26,ks_ptbins);
        laSpectra[mult] = new TH1D(laHistName.str().c_str(),laHistName.str().c_str(),20,la_ptbins);
        ratioHist[mult] = new TH1D(ratioHistName.str().c_str(),ratioHistName.str().c_str(),20,la_ptbins);

        for (pt = 0; pt < 26; pt++){

            double ks_temp = (ks_HM_pTyield[mult][pt]/ks_binwidth[pt])/(2*3.1415926*ks_ptbincenter[pt]*4.8*norm[mult]);

            ksSpectra[mult]->SetBinContent(pt+1, ks_temp );
            ksSpectra[mult]->SetBinError(pt+1, num_ks_err[mult][pt]/(2*3.1415926*ks_binwidth[pt]*ks_ptbincenter[pt]*4.8*norm[mult] ));

        }

        for (pt = 0; pt < 20; pt++){

            double la_temp = (la_HM_pTyield[mult][pt]/la_binwidth[pt])/(2*3.1415926*la_ptbincenter[pt]*4.8*norm[mult]);

            laSpectra[mult]->SetBinContent(pt+1, la_temp );
            laSpectra[mult]->SetBinError(pt+1, num_la_err[mult][pt]/( 2*3.1415926*la_binwidth[pt]*la_ptbincenter[pt]*4.8*norm[mult] ));
            
           
        }

        for(pt = 2; pt < 20; pt++){

            ratioHist[mult]->SetBinContent(pt+1, la_HM_pTyield[mult][pt]/(2*ks_HM_pTyield[mult][pt+6]));

            double err = errorCal_lambdakshort( la_HM_pTyield[mult][pt], ks_HM_pTyield[mult][pt+6], num_la_err[mult][pt], num_ks_err[mult][pt+6] );
            ratioHist[mult]->SetBinError(pt+1, err );

        }
        
    }


    TFile f1("new8Multbins_EPOSvtx_FullStats_v6_26ksbins_pt.root","new");
    
    for( mult = 0; mult < 8; mult++){

        ksSpectra[mult]->Write();
        laSpectra[mult]->Write();
        ratioHist[mult]->Write();
       
    }
        


}