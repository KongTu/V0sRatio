#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void plotEPOSclosureTest(){

	gStyle->SetErrorX(0);

	double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};
    double binwidth[20] = {0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.4,0.6,1.0,1.4};

	TFile* file = new TFile("./EPOShist.root");

	TH1D* ksSpectra = (TH1D*)file->Get("ksSpectra");
	TH1D* genksSpectra = (TH1D*)file->Get("genksSpectra");

	TH1D* laSpectra = (TH1D*)file->Get("laSpectra");
	TH1D* genlaSpectra = (TH1D*)file->Get("genlaSpectra");

	TH1D* closureTest;
	



}