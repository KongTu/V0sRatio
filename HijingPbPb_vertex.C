#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;


void HijingPbPb_vertex(){

	TFile* file[8];

	file[0] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_0_35_2014.root");
	file[1] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_35_60_2014.root");
	file[2] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_60_90_2014.root");
	file[3] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_90_120_2014.root");
	file[4] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_120_150_2014.root");
	file[5] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_150_185_2014.root");
	file[6] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_185_220_2014.root");
	file[7] = new TFile("/Users/kongkong/2014Research/ROOT_file/HijingPbPb/HijingPbPb_TH3D_Nov14_220plus_2014.root");

	TH1D* vertex[8];
	double Nev[8];

	stringstream vertexName;

	for(int mult = 0; mult < 8; mult++){

		vertexName.str("");
		vertexName << "vertexDistZ_";
		vertexName << mult+1;

		vertex[mult] = new TH1D(vertexName.str().c_str(),vertexName.str().c_str(),100,-15,15);

		vertex[mult] = (TH1D*)file[mult]->Get("ana/vertexDistZ");
		Nev[mult] = vertex[mult]->GetEntries();

	}

	TH1D* vertextot = new TH1D("vertextot","vertextot",100,-15,15);

	for(mult = 0; mult < 8; mult++){

		vertextot->Add(vertex[mult],1);
	}

	double weight[8];

	for(mult = 0; mult < 8; mult++){

		double temp = vertextot->GetEntries();
		weight[mult] = temp/Nev[mult];
		vertex[mult]->Scale(weight[mult]);
	}

	TLatex* ratio[9];

	ratio[0] = new TLatex(-5,0.92,"0 < N^{offline}_{trk} < 35");
	ratio[1] = new TLatex(-5,0.92,"35 < N^{offline}_{trk} < 60");
	ratio[2] = new TLatex(-5,0.92,"60 < N^{offline}_{trk} < 90");
	ratio[3] = new TLatex(-5,0.92,"90 < N^{offline}_{trk} < 120");
	ratio[4] = new TLatex(-5,0.92,"120 < N^{offline}_{trk} < 150");
	ratio[5] = new TLatex(-5,0.92,"150 < N^{offline}_{trk} < 185");
	ratio[6] = new TLatex(-5,0.92,"185 < N^{offline}_{trk} < 220");
	ratio[7] = new TLatex(-5,0.92,"220 < N^{offline}_{trk} < #infty");

	TCanvas* c1 = new TCanvas();
	c1->Divide(4,2,0,0);
	for(mult = 0; mult < 8; mult++){

		c1->cd(mult+1);
		vertex[mult]->Divide(vertextot);
		vertex[mult]->GetYaxis()->SetRangeUser(0.9,1.1);
		vertex[mult]->Draw("P");
		ratio[mult]->Draw("same");
	}

	
	for(mult = 0; mult < 8; mult++){
		
		vertexName.str("");
		vertexName << "vertexReweight_HijingPbPb_";
		vertexName << mult+1;
		vertexName << ".root";

		TFile f1(vertexName.str().c_str(),"new");
		vertex[mult]->Write();
	}


	





}