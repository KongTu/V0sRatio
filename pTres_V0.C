#include "fitting.h"

using namespace std;

void pTres_V0()
{
  TFile* f = new TFile("~/Desktop/v0validation_V0EffMatchAna_v14.root");
  TH2D* hKs = (TH2D*)f->Get("v0Validator/ksResolutionPt");
  TH2D* hLam = (TH2D*)f->Get("v0Validator/lamResolutionPt");

  TCanvas* c = new TCanvas("c","",800,400);
  c->Divide(2,1);
  c->cd(1);
  hKs->SetAxisRange(-0.3,0.3,"Y");
  hKs->Draw("colz");
  c->cd(2);
  hLam->SetAxisRange(-0.3,0.3,"Y");
  hLam->Draw("colz");

  double ptbins[] = {0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.6,3.2,4.2,6.0,9.0};
  //double ptbins[] = {0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,2.4,2.6,2.8,3.0,3.4,3.8,4.2,4.6,5.0,5.6,6.6,9.0};

  TH1D* hKsPtRes[15];
  TH1D* hLamPtRes[15];
  TF1* fitKs[15];
  TF1* fitLam[15];

  TH1D* hResVsPt_Ks = new TH1D("resvspt_ks",";p_{T} (GeV);#sigma(p_{T})/p_{T} for K_{s}^{0}",15,&ptbins[0]);
  TH1D* hResVsPt_Lam = new TH1D("resvspt_lam",";p_{T} (GeV);#sigma(p_{T})/p_{T} for #Lambda",15,&ptbins[0]);
  TH1D* hResVsPtAbs_Ks = new TH1D("resvsptabs_ks",";p_{T} (GeV);#sigma(p_{T}) (GeV) for K_{s}^{0}",15,&ptbins[0]);
  TH1D* hResVsPtAbs_Lam = new TH1D("resvsptabs_lam",";p_{T} (GeV);#sigma(p_{T}) (GeV) for #Lambda",15,&ptbins[0]);

  for(int i=0;i<15;i++)
  {
    int nminpt = hKs->GetXaxis()->FindBin(ptbins[i]+0.0001);
    int nmaxpt = hKs->GetXaxis()->FindBin(ptbins[i+1]-0.0001);

    hKsPtRes[i] = hKs->ProjectionY(Form("ks_%d",i),nminpt,nmaxpt,"e");
    hLamPtRes[i] = hLam->ProjectionY(Form("lam_%d",i),nminpt,nmaxpt,"e");

    hKsPtRes[i]->Scale(1.0/hKsPtRes[i]->Integral());
    hLamPtRes[i]->Scale(1.0/hLamPtRes[i]->Integral());
/*
    fitKs[i] = new TF1(Form("fit_ks_%d",i),"[0]*exp(-0.5*x*x/[1]/[1])+[2]*exp(-0.5*x*x/[1]/[1]/9)",-0.3,0.3);
    fitLam[i] = new TF1(Form("fit_lam_%d",i),"[0]*exp(-0.5*x*x/[1]/[1])+[2]*exp(-0.5*x*x/[1]/[1]/9)",-0.3,0.3);
    fitKs[i]->SetParameters(0.1,0.02,0.05,0.04);
    fitLam[i]->SetParameters(0.1,0.02,0.05,0.04);
*/
    fitKs[i] = new TF1(Form("fit_ks_%d",i),"[0]*exp(-0.5*x*x/[1]/[1])",-0.2,0.2);
    fitLam[i] = new TF1(Form("fit_lam_%d",i),"[0]*exp(-0.5*x*x/[1]/[1])",-0.2,0.2);
    fitKs[i]->SetParameters(0.1,0.01);
    fitLam[i]->SetParameters(0.1,0.005);
    hKsPtRes[i]->Fit(Form("fit_ks_%d",i),"RNO");
    hLamPtRes[i]->Fit(Form("fit_lam_%d",i),"RNO");

    hResVsPt_Ks->SetBinContent(i+1,fitKs[i]->GetParameter(1)/hResVsPt_Ks->GetBinCenter(i+1));
    hResVsPt_Lam->SetBinContent(i+1,fitLam[i]->GetParameter(1)/hResVsPt_Lam->GetBinCenter(i+1));
    hResVsPtAbs_Ks->SetBinContent(i+1,fitKs[i]->GetParameter(1));
    hResVsPtAbs_Lam->SetBinContent(i+1,fitLam[i]->GetParameter(1));
  }

  TCanvas* c1 = new TCanvas("c1","",1000,1000);  
  c1->Divide(4,4);
  for(int i=0;i<15;i++)
  {
    c1->cd(i+1);
    hKsPtRes[i]->SetMarkerStyle(15);
    hKsPtRes[i]->SetTitle(Form("%.1f-%.1f",ptbins[i],ptbins[i+1]));   
    hKsPtRes[i]->Draw("PE");
    fitKs[i]->Draw("Lsame");
  }

  TCanvas* c2 = new TCanvas("c2","",1000,1000);
  c2->Divide(4,4);
  for(int i=0;i<15;i++)
  {
    c2->cd(i+1);
    hLamPtRes[i]->SetMarkerStyle(15);
    hLamPtRes[i]->SetTitle(Form("%.1f-%.1f",ptbins[i],ptbins[i+1]));
    hLamPtRes[i]->Draw("PE");
    fitLam[i]->Draw("Lsame");
  }

  TCanvas* c3 = new TCanvas("c3","",800,400);
  c3->Divide(2,1);
  c3->cd(1);
  hResVsPt_Ks->Draw("PE");
  c3->cd(2);
  hResVsPt_Lam->SetAxisRange(0.4,9.0,"X");
  hResVsPt_Lam->Draw("PE");

  TCanvas* c4 = new TCanvas("c4","",800,400);
  c4->Divide(2,1);
  c4->cd(1);
  hResVsPtAbs_Ks->Draw("PE");
  c4->cd(2);
  hResVsPtAbs_Lam->SetAxisRange(0.4,9.0,"X");
  hResVsPtAbs_Lam->Draw("PE");
/*
  TFile file("momentumResSmearHist.root","new");

  for(int i = 0; i < 15; i++){

  hKsPtRes[i]->Write();
  hLamPtRes[i]->Write();
}
*/

}


