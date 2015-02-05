#ifndef	_LATTICE_H
#define _LATTICE_H

#include "MersenneTwister.h"
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <map>

using namespace std;

const double pi=3.141592654;
template <class T>
class UPDATER
{
	public:
		UPDATER(map<string,double> params){
			lat=T(params);
			sum=0.0;
			if (lat.updateWeights.size()!=lat.updateFuncs.size()){
				cout<<"less weights than functions!"<<endl;
				exit(0);
			}
			for (int i=0;i<lat.updateWeights.size();i++){
				sum+=lat.updateWeights[i];
			}
			nTries=vector<int>(lat.updateWeights.size(),0); 
			nAccepts=vector<int>(lat.updateWeights.size(),0); 
		}		
		void random_update(int nsteps){
			double w;
			int site;
			double runningSum=0.0;
			for(int j=0;j<nsteps;j++){
				for(int i=0;i<lat.getN();i++){
					runningSum=0.0;
					w=lat.ran.randExc(sum);
					site=lat.ran.randInt(lat.getD()*lat.getN()-1);
					for(int k=0;k<lat.updateWeights.size();k++){
						runningSum+=lat.updateWeights[k];
						if(w<runningSum){
							call(k,site);
							break;
						}
					}
				}
			}
				
		}
		void rate(){	
			for(int i=0;i<lat.updateWeights.size();i++){
				cout<<(1.0*nAccepts[i])/(1.0*nTries[i])<<" ";
			}cout<<endl;
		}
		
		UPDATER<T>& operator=(const UPDATER<T> &rhs){
			if(this==&rhs)
				return *this;

			return *this;
		}
		double sum;
		T lat;
private:
		vector<int> nTries;
		vector<int> nAccepts;
		void call(int func,int site){
			nTries[func]++;
			nAccepts[func]+=(lat.*lat.updateFuncs[func])(site);
		}
//		int (T::*foo)(int site);
		double weight;
};
class LATTICE{
public:
	LATTICE(map<string,double> &params);
	LATTICE();
	virtual double energy()=0;
	double E();//the running energy
	int getN();int getD();
	MTRand ran;
	double mod2pi(double);
	int pow(int,int);
	double pow(double,int);
	int leviCivita(int a ,int b){
		if((b-a)%3==1) return 1;
		else return -1;
	}
	template <class T>
	double curl(const vector<vector <T> > &in,int i,int d1,int d2){
		return leviCivita(d1,d2)*(in[(this->*p)(i,d1)][d2]-in[i][d2]-in[(this->*p)(i,d2)][d1]+in[i][d1]);
	}
	
protected:
	int L,Lt;//linear size
	int N;//total number of sites
	int DIMS;
	int (LATTICE::*p)(int site,int dir);
	int (LATTICE::*m)(int site,int dir);
	virtual int p_normal(int site,int dir);//connection functions
	virtual int m_normal(int site,int dir);
	double runningE;
	vector<vector<int> >end;//1 if that site is at the end of the lattice in that direction
};	

const double macheps=0.000000000000000001;

class VTABLE{
public:
    void initVTABLE(double t,int func);
    void differentiate(VTABLE&);
    double value(double);
    double value(int);
    double v1(double t,double j);
    double getStep();
    double getXmin();
    double getL();

private:
    int L;
    double step;
    double xmin;
    vector<double> V;
    double (VTABLE::*eval)(double,double);
    double simple(double,double);
    double energy1(double,double);
    double Dsimple(double,double);
    double D2simple(double,double);
};

#endif
