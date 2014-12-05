#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include <stdio.h>
#include <iomanip>
#include <vector>
#include "TGraph.h"
#include "TF1.h"
#include "TH1D.h"
#include "TFile.h"

using namespace std;

void plotBeta_T(){

	TGraph* g = new TGraph(7);
	g->SetPoint(0,0.824,0.152);
	g->SetPoint(1,0.815,0.157);
	g->SetPoint(2,0.810,0.158);
	g->SetPoint(3,0.803,0.164);
	g->SetPoint(4,0.786,0.169);
	g->SetPoint(5,0.772,0.176);
	g->SetPoint(6,0.776,0.165);

	g->SetMarkerStyle(20);
	g->SetMarkerSize(1.3);
	g->SetMarkerColor(kBlue);

	g->Draw("AP");
}