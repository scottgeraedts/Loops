#ifndef _RODS_H
#define _RODS_H

#include "lattice.h"

class RODS:public LATTICE
{
public:
	RODS(map<string,double> &params);
	RODS();
	double magnetization();
	double rho();
	double QQ();
	double energy();
	vector<int (RODS::*)(int)>updateFuncs;
	vector<double> updateWeights;
private:
	int updateA(int site);
	int updatePP(int site);
	int updateGamma(int site);
	int updatePhi(int site);
	int updateP(int site);
	double cosTerm(int site,int d);
	double angle_step;	
	double rhoHelper(const vector< vector<int> > &in, int,int);
	vector<double> gamma;
	vector<double> phi;
	vector<vector<int> > a;
	vector<vector<int> > pp;
	double t1,t2,theta;
};

#endif
