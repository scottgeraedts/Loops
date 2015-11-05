#ifndef _RODS_H
#define _RODS_H

#include "lattice.h"
#include <complex>
#include <sstream>

class RODS:public LATTICE
{
public:
	RODS(map<string,double> &params);
	RODS();
	vector<int (RODS::*)(int)>updateFuncs;
	vector<double> updateWeights;
	
	double magnetization();
	double Z2magnetization();
	double rho();
	double QQ();
	double energy();
	void updateCorrelators();
	void printCorrelators(int);

private:
	int updateA(int site);
	int updatePP(int site);
	int updateGamma(int site);
	int updatePhi(int site);
	int compAP(int site);
	int starUpdate(int site);
	
	double cosTerm(int site,int d);
	double angle_step;
	double shiftedphi(int);	
	double rhoHelper(const vector< vector<int> > &in, int,int);
	vector<double> gamma;
	vector<double> phi;
	vector<vector<int> > a;
	vector<vector<int> > pp;
	double t1,t2,theta;
	double lambda1,lambda2;
	int edge;
	vector<vector<vector<complex<double> > > > cors;
	vector< complex<double> > avgphi;
};

#endif
