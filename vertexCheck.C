#include "fitting.h"

using namespace std;

TSystem->Load("libRooFit");

using namespace RooFit;

void vertexCheck(){

	TFile* file1 = new TFile("/Users/kongkong/2014Research/ROOT_file/newEPOSsample/EPOS_TH3D_Nov21_vtx_2014.root");
	TFile* file2 = new TFile("/Users/kongkong/2014Research/ROOT_file/new_pPb_sample/HM_Nov10_220plus.root");

	TH1D* vertex_EPOS = (TH1D*)file1->Get("ana/vertexDistZ")->Clone("epos");
	TH1D* vertex_data = (TH1D*)file2->Get("ana/vertexDistZ")->Clone("data");

	double ratio = (vertex_data->GetEntries()/vertex_EPOS->GetEntries());
	vertex_EPOS->Scale(ratio);
	vertex_data->Divide( vertex_EPOS );
	vertex_data->Draw();

	TFile f1("vertex_epos.root","new");
	vertex_data->Write();

}