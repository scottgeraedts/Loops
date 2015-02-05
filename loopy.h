#ifndef _LOOPY_H
#define _LOOPY_H

#include "lattice.h"

class LOOPY:public LATTICE
{
public:
	LOOPY(map<string,double> &params);
	LOOPY();
	double magnetization();
	double rho();
	double energy();
	vector<int (LOOPY::*)(int)>updateFuncs;
	vector<double> updateWeights;
private:
	VTABLE table;
	int updateA(int site);
	int updateGamma(int site);
	int updatePhi(int site);
	double cosTerm(int site,int d);
	double angle_step;	
	double rhoHelper(int,int);
	vector<double> gamma;
	vector<double> phi;
	vector<vector<int> > a;
	double t1,t2,theta;
};

#endif
