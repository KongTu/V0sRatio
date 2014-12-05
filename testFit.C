#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void testFit(){

	TFile* file = new TFile("/Users/kongkong/2014Research/Code/Jet'study/gitV0sRatio/histoSpectraFolder_rpy/HMbins_18pTBins5rpyBins_wSmooth10Eff_FullStats_5Multbins.root");

	TH1D* ratioHist[5][5];
    TH1D* ksSpectra[5][5];
    TH1D* laSpectra[5][5];

	stringstream ratioHistName;
    stringstream ksSpectraName;
    stringstream laSpectraName;

	for(int mult = 0; mult < 5; mult++){

		for(int rpy = 0; rpy < 5; rpy++){

			ratioHistName.str("");
			ratioHistName << "ratioHist_";
			ratioHistName << mult+1;
			ratioHistName << "_";
			ratioHistName << rpy+1;

            ksSpectraName.str("");
            ksSpectraName << "ksSpectra_";
            ksSpectraName << mult+1;
            ksSpectraName << "_";
            ksSpectraName << rpy+1;

            laSpectraName.str("");
            laSpectraName << "laSpectra_";
            laSpectraName << mult+1;
            laSpectraName << "_";
            laSpectraName << rpy+1;

			ratioHist[mult][rpy] = (TH1D*)file->Get( ratioHistName.str().c_str() );
            ksSpectra[mult][rpy] = (TH1D*)file->Get( ksSpectraName.str().c_str() );
            laSpectra[mult][rpy] = (TH1D*)file->Get( laSpectraName.str().c_str() );


		}
	}

	EasyBlastWaveFit( ksSpectra[2][3] );
	
}