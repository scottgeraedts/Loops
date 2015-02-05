#ifndef _LOOPY_H
#define _LOOPY_H

#include "lattice.h"
#include <stack>

class ISING:public LATTICE
{
public:
	ISING(map<string,double> &params);
	ISING();
	~ISING();
	double magnetization();
	double energy();
	vector<int (ISING::*)(int)>updateFuncs;
	vector<double> updateWeights;
	void draw();
private:
	int updateSpin(int site);
	int updateWolff(int site);
	vector<int> s;
	double J,h;
	int drawn;
};

#endif
