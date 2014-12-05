#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void vertexCheck(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/3MCsamples_comparison/EPOS_pPb_TH3D_Dec3_eta_vertexreweighed_2014.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/MB_Nov10_90_120.root");
	TH1D* vertex_hijing = (TH1D*)file1->Get("ana/vertexDistZ")->Clone("hijing");
	TH1D* vertex_data = (TH1D*)file2->Get("ana/vertexDistZ")->Clone("data");

	double ratio = (vertex_hijing->GetEntries()/vertex_data->GetEntries());
	vertex_data->Scale(ratio);
	vertex_data->Divide( vertex_hijing );
	vertex_data->Draw();

	//TFile f1("vertex_hijing.root","new");
	//vertex_data->Write();

}