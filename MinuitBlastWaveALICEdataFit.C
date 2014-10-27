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

vector<double> x,ex, y,ey,x1,ex1,y1,ey1;

double getChi2;
int getNdf;

const double R = 1.0;
const double dr = 1e-3; // FIXME

double integral(double beta_s, double T, double n, double pt, double mt)
{
  
  double s = 0.;
  for(double r = dr/2; r < R; r += dr)
  {
    double beta_r = beta_s * TMath::Power((r/R),n);
    double rho = TMath::ATanH(beta_r);

    s += r * dr * TMath::BesselK1((mt*cosh(rho))/T) * TMath::BesselI0((pt*sinh(rho))/T);
  }

  return s;
}

void function(int &npar, double *gin, double &f, double *par, int flag)
{
  double aka    = par[0];
  double aka1   = par[4];
  double T      = par[1];
  double beta_s = par[2];
  double n = par[3];

  double chi2 = 0.;

  for(int i = 0; i < int(x.size()); i++)
  {
    double m = 0.497;

    double pt = x[i];
    double mt = sqrt(m*m + pt*pt);

    double a = 0.;
    a = aka;
  
    double theo = a * mt * integral(beta_s, T, n, pt, mt);
    double q = (y[i] - theo)/ey[i];

    chi2 += q*q;
  }

  for(i = 0; i < int(x1.size()); i++){

  	double m1 = 1.115;

  	double pt1 = x1[i];
  	double mt1 = sqrt(m1*m1 + pt1*pt1);

  	double a1 = 0.;
  	a1 = aka1;

  	double theo1 = a1 * mt1 * integral(beta_s,T,n,pt1,mt1);
  	double q1 = (y1[i]-theo1)/ey1[i];

  	chi2 += q1*q1;
  }

  f = chi2;

  getChi2 = f;
  getNdf  = x.size() + x1.size() - 5 - 1;
}

double MyFunc( double *pt, double *p){

  double mass = 0.497;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}

double MyFunc1( double *pt, double *p){

  double mass = 1.115;

  double temp = 0.;
    double mt = sqrt(pt[0]*pt[0]+mass*mass);

    temp = p[2] * mt * integral(p[0],p[1],p[3],pt[0],mt);

    return temp;
  
}


void MinuitBlastWaveALICEdataFit(){

/**
 * manually read all the data points from ALICE's result and put into vectors, x, ex, y, ey;
 */

	double ptbins[] = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.3,3.6,3.9,4.2,4.6,5.0,5.5,6.0,8.0};

	double ptbincenter_ks[] = { 0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 
    0.95, 1.05, 1.15, 1.25, 1.35, 1.45, 1.55, 1.65, 1.75, 1.85, 
    1.95, 2.1, 2.3, 2.5, 2.7, 2.9, 3.15, 3.45, 3.75, 4.05, 
    4.4, 4.8, 5.25, 5.75, 7.0 };

    double ptbincenter_la[] = { 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.3, 1.5, 1.7, 
    1.9, 2.1, 2.3, 2.5, 2.7, 3.0, 3.45, 3.95, 4.6, 5.5, 
    7.0 };
	

	double yields_ks1[] = { 1.532, 1.373, 1.22, 1.026, 0.8654, 0.7002, 0.5571, 0.4442, 0.3476, 
	0.2745, 0.2161, 0.1731, 0.1342, 0.1068, 0.0867, 0.06848, 0.0553, 0.04467, 0.03537, 
	0.02874, 0.02089, 0.01407, 0.009662, 0.006504, 0.004485, 0.002855, 0.001833, 0.00112, 7.169E-4, 
	4.67E-4, 2.456E-4, 1.467E-4, 9.063E-5, 2.868E-5 };

	double yields_ks_error1[] = { 0.161, 0.027, 0.015, 0.011, 0.0076, 0.0053, 0.0039, 0.003, 0.0023, 
    0.0019, 0.0015, 0.0013, 0.0011, 9.0E-4, 7.9E-4, 6.7E-4, 5.8E-4, 5.0E-4, 4.3E-4, 
    3.8E-4, 2.2E-4, 1.7E-4, 1.34E-4, 1.06E-4, 8.5E-5, 5.2E-5, 4.0E-5, 3.0E-5, 2.28E-5, 
    1.57E-5, 1.09E-5, 7.1E-6, 5.37E-6, 1.39E-6 };

	double yields_ks2[] = { 1.282, 1.125, 0.9841, 0.8456, 0.6969, 0.5672, 0.447, 0.3563, 0.2742, 
	0.2176, 0.1696, 0.1349, 0.1055, 0.08433, 0.06682, 0.05295, 0.04317, 0.03504, 0.02797, 
	0.02255, 0.01616, 0.01112, 0.007565, 0.005198, 0.00362, 0.002319, 0.001436, 9.019E-4, 5.672E-4, 
	3.625E-4, 2.299E-4, 1.311E-4, 6.245E-5, 2.493E-5 };

	double yields_ks_error2[] = { 0.145, 0.024, 0.0137, 0.0096, 0.0067, 0.0047, 0.0034, 0.0026, 0.002, 
	0.0016, 0.0013, 0.0011, 9.0E-4, 7.9E-4, 6.7E-4, 5.7E-4, 4.9E-4, 4.3E-4, 3.7E-4, 
	3.2E-4, 1.8E-4, 1.4E-4, 1.14E-4, 9.1E-5, 7.3E-5, 4.5E-5, 3.4E-5, 2.62E-5, 1.97E-5, 
	1.34E-5, 1.02E-5, 6.6E-6, 4.24E-6, 1.26E-6 };

	double yields_ks3[] = { 1.032, 0.9984, 0.8801, 0.7261, 0.5949, 0.4825, 0.3799, 0.2988, 0.2292, 
	0.1805, 0.141, 0.1106, 0.0862, 0.06891, 0.05424, 0.04371, 0.03475, 0.02833, 0.02264, 
	0.01838, 0.01337, 0.009004, 0.006065, 0.004339, 0.003005, 0.001894, 0.001185, 7.542E-4, 4.686E-4, 
	3.037E-4, 1.802E-4, 1.081E-4, 6.501E-5, 2.173E-5 };

	double yields_ks_error3[] = { 0.094, 0.0159, 0.0092, 0.0063, 0.0044, 0.0031, 0.0023, 0.0017, 0.0014, 
	0.0011, 9.0E-4, 8.0E-4, 6.3E-4, 5.4E-4, 4.6E-4, 4.0E-4, 3.4E-4, 3.0E-4, 2.5E-4, 
	2.2E-4, 1.3E-4, 1.0E-4, 7.8E-5, 6.4E-5, 5.1E-5, 3.1E-5, 2.4E-5, 1.84E-5, 1.37E-5, 
	9.4E-6, 7.0E-6, 4.6E-6, 3.42E-6, 9.1E-7 };

	double yields_ks4[] = { 0.8664, 0.8167, 0.7138, 0.5872, 0.4752, 0.3726, 0.2879, 0.2252, 0.172, 
	0.1339, 0.1043, 0.0819, 0.06291, 0.04975, 0.03974, 0.03109, 0.02512, 0.02004, 0.01625, 
	0.01301, 0.009402, 0.006303, 0.004398, 0.003022, 0.002128, 0.001354, 8.566E-4, 5.426E-4, 3.412E-4, 
	2.217E-4, 1.309E-4, 7.942E-5, 4.63E-5, 1.598E-5 };

	double yields_ks_error4[] = { 0.0587, 0.0103, 0.0059, 0.0041, 0.0029, 0.002, 0.0015, 0.0011, 9.0E-4, 
	7.0E-4, 6.0E-4, 5.0E-4, 4.1E-4, 3.5E-4, 3.0E-4, 2.5E-4, 2.2E-4, 1.9E-4, 1.6E-4, 
	1.4E-4, 8.2E-5, 6.4E-5, 5.1E-5, 4.1E-5, 3.3E-5, 2.0E-5, 1.55E-5, 1.2E-5, 9.0E-6, 
	6.1E-6, 4.5E-6, 3.03E-6, 2.22E-6, 6.1E-7 };

	double yields_ks5[] = { 0.6602, 0.6161, 0.529, 0.4355, 0.3441, 0.2617, 0.2011, 0.1554, 0.1174, 
	0.08973, 0.06793, 0.05317, 0.04105, 0.03174, 0.0248, 0.01966, 0.01595, 0.01258, 0.009976, 
	0.008011, 0.005954, 0.003914, 0.002692, 0.001828, 0.001312, 8.455E-4, 5.218E-4, 3.186E-4, 2.097E-4, 
	1.396E-4, 8.241E-5, 4.926E-5, 2.905E-5, 9.493E-6 };

	double yields_ks_error5[] = { 0.0534, 0.0087, 0.005, 0.0034, 0.0024, 0.0016, 0.0012, 9.0E-4, 7.0E-4, 
	5.5E-4, 4.4E-4, 3.7E-4, 3.1E-4, 2.6E-4, 2.2E-4, 1.8E-4, 1.6E-4, 1.4E-4, 1.17E-4, 
	1.01E-4, 5.9E-5, 4.5E-5, 3.6E-5, 2.8E-5, 2.3E-5, 1.45E-5, 1.09E-5, 8.2E-6, 6.3E-6, 
	4.4E-6, 3.22E-6, 2.14E-6, 1.59E-6, 4.19E-7 };

	double yields_ks6[] = { 0.5099, 0.4278, 0.3651, 0.2804, 0.2207, 0.1663, 0.1246, 0.09248, 0.06821, 
	0.0513, 0.03895, 0.02908, 0.02266, 0.01756, 0.01333, 0.01052, 0.008219, 0.006613, 0.005288, 
	0.004173, 0.003059, 0.002031, 0.001332, 9.469E-4, 6.523E-4, 4.37E-4, 2.7E-4, 1.664E-4, 1.1E-4, 
	6.73E-5, 4.216E-5, 2.427E-5, 1.468E-5, 4.802E-6 };

	double yields_ks_error6[] = { 0.0427, 0.007, 0.0039, 0.0026, 0.0018, 0.0012, 9.0E-4, 6.6E-4, 4.9E-4, 
	3.9E-4, 3.1E-4, 2.5E-4, 2.1E-4, 1.7E-4, 1.4E-4, 1.2E-4, 1.02E-4, 8.9E-5, 7.6E-5, 
	6.5E-5, 3.8E-5, 2.9E-5, 2.2E-5, 1.83E-5, 1.47E-5, 9.3E-6, 7.0E-6, 5.3E-6, 4.1E-6, 
	2.71E-6, 2.07E-6, 1.34E-6, 9.9E-7, 2.61E-7 };

	double yields_ks7[] = { 0.2596, 0.2172, 0.1874, 0.1448, 0.1059, 0.07555, 0.05461, 0.03886, 0.0274, 
	0.02033, 0.01438, 0.01057, 0.008109, 0.006127, 0.004504, 0.003565, 0.002764, 0.002139, 0.001738, 
	0.001358, 9.955E-4, 6.173E-4, 4.244E-4, 2.994E-4, 1.98E-4, 1.296E-4, 7.599E-5, 5.017E-5, 2.996E-5, 
	1.953E-5, 1.094E-5, 8.026E-6, 4.519E-6, 1.483E-6 };

	double yields_ks_error7[] = { 0.0312, 0.0051, 0.0029, 0.0019, 0.0013, 8.7E-4, 6.1E-4, 4.4E-4, 3.2E-4, 
	2.5E-4, 1.9E-4, 1.5E-4, 1.22E-4, 1.01E-4, 8.2E-5, 6.9E-5, 5.7E-5, 4.9E-5, 4.2E-5, 
	3.6E-5, 2.07E-5, 1.54E-5, 1.21E-5, 9.8E-6, 7.7E-6, 4.9E-6, 3.46E-6, 2.73E-6, 2.0E-6, 
	1.4E-6, 9.9E-7, 7.13E-7, 5.1E-7, 1.4E-7 };
	

	

	double yields_la1[] = { 0.1833, 0.1867, 0.1737, 0.1574, 0.1356, 0.1232, 0.1006, 0.07763, 0.05522, 
    0.03978, 0.02885, 0.02024, 0.01422, 0.009873, 0.005879, 0.002652, 0.001182, 4.378E-4, 1.2E-4, 
    2.194E-5 };

  	double yields_la_error1[] = { 0.0035, 0.0028, 0.0024, 0.0021, 0.0018, 0.0017, 0.001, 8.2E-4, 6.3E-4, 
    5.0E-4, 4.0E-4, 3.2E-4, 2.6E-4, 2.04E-4, 1.08E-4, 6.0E-5, 3.7E-5, 1.62E-5, 6.5E-6, 
    1.62E-6 };

    double yields_la2[] = { 0.1526, 0.1549, 0.1423, 0.1275, 0.1129, 0.09953, 0.07854, 0.06002, 0.04264, 
    0.03068, 0.02204, 0.01533, 0.01071, 0.007423, 0.004456, 0.002029, 9.097E-4, 3.453E-4, 8.98E-5, 
    1.818E-5 };

    double yields_la_error2[] = { 0.0031, 0.0025, 0.0021, 0.0018, 0.0016, 0.00145, 8.6E-4, 6.9E-4, 5.3E-4, 
    4.1E-4, 3.3E-4, 2.6E-4, 2.1E-4, 1.64E-4, 8.7E-5, 4.9E-5, 3.03E-5, 1.35E-5, 5.25E-6, 
    1.42E-6 };

    double yields_la3[] = { 0.1339, 0.1341, 0.122, 0.111, 0.09506, 0.08272, 0.06756, 0.04923, 0.03472, 
    0.02451, 0.01749, 0.01222, 0.008372, 0.005975, 0.003456, 0.001598, 7.156E-4, 2.673E-4, 7.685E-5, 
    1.365E-5 };

    double yields_la_error3[] = { 0.0021, 0.0017, 0.0015, 0.0013, 0.00113, 0.00102, 6.2E-4, 4.9E-4, 3.8E-4, 
    2.9E-4, 2.3E-4, 1.8E-4, 1.46E-4, 1.18E-4, 6.1E-5, 3.4E-5, 2.13E-5, 9.3E-6, 3.84E-6, 
    9.6E-7 };

    double yields_la4[] = { 0.1125, 0.1103, 0.0972, 0.08599, 0.07332, 0.0637, 0.04987, 0.03599, 0.02512, 
    0.01761, 0.01211, 0.008437, 0.005732, 0.003984, 0.002376, 0.001083, 4.723E-4, 1.81E-4, 5.091E-5, 
    9.54E-6 };

    double yields_la_error4[] = { 0.0014, 0.0012, 0.001, 8.8E-4, 7.7E-4, 6.9E-4, 4.2E-4, 3.3E-4, 2.5E-4, 
    1.9E-4, 1.5E-4, 1.19E-4, 9.4E-5, 7.5E-5, 3.9E-5, 2.2E-5, 1.34E-5, 5.9E-6, 2.43E-6, 
    6.18E-7 };

    double yields_la5[] = { 0.08493, 0.08226, 0.07099, 0.06082, 0.05103, 0.04309, 0.03307, 0.02325, 0.01569, 
    0.01069, 0.007296, 0.004859, 0.003287, 0.002275, 0.001353, 5.792E-4, 2.73E-4, 1.059E-4, 2.818E-5, 
    5.662E-6 };

    double yields_la_error5[] = { 0.00118, 9.6E-4, 8.0E-4, 6.8E-4, 5.9E-4, 5.2E-4, 3.0E-4, 2.3E-4, 1.7E-4, 
    1.3E-4, 1.0E-4, 7.6E-5, 6.0E-5, 4.7E-5, 2.5E-5, 1.33E-5, 8.6E-6, 3.9E-6, 1.55E-6, 
    4.21E-7 };

    double yields_la6[] = { 0.05778, 0.05376, 0.04462, 0.03812, 0.0301, 0.02502, 0.01825, 0.01223, 0.008084, 
    0.005206, 0.003501, 0.002313, 0.001543, 0.001023, 5.947E-4, 2.701E-4, 1.134E-4, 4.36E-5, 1.311E-5, 
    2.355E-6 };

    double yields_la_error6[] = { 9.1E-4, 7.3E-4, 5.8E-4, 4.9E-4, 4.0E-4, 3.5E-4, 2.0E-4, 1.4E-4, 1.05E-4, 
    7.7E-5, 5.8E-5, 4.4E-5, 3.4E-5, 2.6E-5, 1.35E-5, 7.5E-6, 4.5E-6, 2.04E-6, 9.1E-7, 
    2.39E-7 };

    double yields_la7[] = { 0.02595, 0.02352, 0.01842, 0.0153, 0.01131, 0.008895, 0.006325, 0.00397, 0.002447, 
    0.001529, 9.52E-4, 6.399E-4, 4.256E-4, 2.622E-4, 1.625E-4, 6.704E-5, 3.039E-5, 1.031E-5, 3.099E-6, 
    6.841E-7 };

    double yields_la_error7[] = { 6.0E-4, 4.6E-4, 3.6E-4, 2.9E-4, 2.3E-4, 1.92E-4, 1.07E-4, 7.5E-5, 5.2E-5, 
    3.8E-5, 2.74E-5, 2.09E-5, 1.62E-5, 1.2E-5, 6.3E-6, 3.34E-6, 2.12E-6, 9.2E-7, 4.2E-7, 
    1.215E-7 };



	for (int pt = 0; pt < 16; pt++){

		x.push_back( ptbincenter_ks[pt] );
		ex.push_back(0.0);
		y.push_back( yields_ks7[pt] );
		ey.push_back( yields_ks_error7[pt] );

	}

	for (pt = 0; pt < 15; pt++){

		x1.push_back( ptbincenter_la[pt] );
		ex1.push_back(0.0);
		y1.push_back( yields_la7[pt] );
		ey1.push_back( yields_la_error7[pt] );

	}

/**
 * simultaneous fit for K0s and Lambda:
 */

	TMinuit * gMinuit = new TMinuit(5);

	//set the function to minimize with minuit;
	gMinuit->SetFCN(function);

	double arglist[10];
	arglist[0] = 0.001;
	gMinuit->mnexcm("SET ERR", arglist, 1, 0);

	gMinuit->mnparm(0, "aka",    10,   0.1,  1,    10000, 0);
	gMinuit->mnparm(1, "Tkin",   0.15, 0.01, 0.05, 1.0,   0);
	gMinuit->mnparm(2, "beta_s", 0.70, 0.01, 0.15, 1.0,   0);
	gMinuit->mnparm(3, "n",      1.1,  0.01, 1.0,  10.0,  0);
	gMinuit->mnparm(4, "aka1",   10,   0.1,  1,    10000, 0);

	
	gMinuit->mnexcm("MIGRAD",    arglist,  1,   0);
	gMinuit->mnexcm("CALL FCN",  arglist,  1,   0);
	gMinuit->mnexcm("HESSE",     arglist,  1,   0);

	double beta_s,Tkin,aka,aka1,n;
	double ebeta_s,eTkin,eaka,eaka1,en;

	gMinuit->GetParameter(0,aka,eaka);
	gMinuit->GetParameter(1,Tkin,eTkin);
	gMinuit->GetParameter(2,beta_s,ebeta_s);
	gMinuit->GetParameter(3,n,en);
	gMinuit->GetParameter(4,aka1,eaka1);


/**
 * define the fit function by using the parameters from fit:
 */

 	double xmin = 0.0;
 	double xmax = 1.5;

 	double xmin1 = 0.6;
 	double xmax1 = 1.5;

	TF1* f1 = new TF1("f1",MyFunc,xmin,xmax,4);

	f1->SetParameter(0,beta_s);
	f1->SetParameter(1,Tkin);
	f1->SetParameter(2,aka);
	f1->SetParameter(3,n);

 	TF1* f2 = new TF1("f2",MyFunc1,xmin1,xmax1,4);

	f2->SetParameter(0,beta_s);
	f2->SetParameter(1,Tkin);
	f2->SetParameter(2,aka1);
	f2->SetParameter(3,n);


	TGraph* g_ks = new TGraph(34);
	TGraph* g_la = new TGraph(20);

	for (pt = 0; pt < 15; pt++){

		g_ks->SetPoint(pt,x[pt],y[pt]);

	}

	g_ks->SetMarkerStyle(20);
	g_ks->SetMarkerSize(1.3);
	g_ks->SetMarkerColor(kBlue);

	for (pt = 0; pt < 15; pt++){

		g_la->SetPoint(pt,x1[pt],y1[pt]);
		
	}

	g_la->SetMarkerStyle(20);
	g_la->SetMarkerSize(1.3);
	g_la->SetMarkerColor(kRed);


	TCanvas* c1 = new TCanvas();

	c1->SetLogy(1);

	g_ks->Draw("AP");
	//g_ks-GetXaxis()->SetRangeUser(0,3.0);
	
	g_la->Draw("Psame");
	f2->Draw("same");
	f1->Draw("same");






}