#include "lattice.h"

LATTICE::LATTICE(map<string,double> &params){
	L=params["L"];
	if(params.find("Lt")!=params.end()) Lt=params["Lt"];
	else Lt=L;
	if (params.count("DIMS")) DIMS=params["DIMS"];
	else DIMS=3;
	if (DIMS>3){
		cout<<"dimensions greater than 3 not supported, update p and m"<<endl;
		exit(0);
	}
	N=Lt*pow(L,DIMS-1);
	runningE=0.0;
	ran=MTRand(params["seed"]);
	vector<int>temp3(DIMS,0);
	end=vector<vector<int> >(N,temp3);
	for (int i=0;i<N;i++){
		if (i%Lt==Lt-1) end[i][0]==1;
		for (int d=1;d<DIMS;d++){
			if( (i/(Lt*pow(L,d-1)) )%L==L-1) end[i][d]==1;
		}
	}
	p=&LATTICE::p_normal;
	m=&LATTICE::m_normal;
}
//LATTICE& LATTICE::operator=(const LATTICE &rhs){
//			if(this==&rhs)
//				return *this;

//			return *this;
//}
LATTICE::LATTICE(){}
double LATTICE::E(){ return runningE; }


//default lattice connection functions
inline int LATTICE::p_normal(int i,int dir){
    if(dir==0){
        if(i%Lt==Lt-1) return i-(Lt-1);
        else return i+1;
    }
    else if(dir==1){
        if((i/Lt)%L==L-1) return i-(L-1)*Lt;
        else return i+Lt;
    }else{
        if((i/(Lt*L))%L==L-1) return i-(L-1)*Lt*L;
        else return i+Lt*L;
    }
}
inline int LATTICE::m_normal(int i,int dir){
    if(dir==0){
        if(i%Lt==0) return i+(Lt-1);
        else return i-1;
    }else if(dir==1){
        if((i/Lt)%L==0) return i+(L-1)*Lt;
        else return i-L;
    }else{
        if((i/(L*Lt))%L==0) return i+(L-1)*Lt*L;
        else return i-Lt*L;
    }
}
//template<class T>
//double LATTICE::curl(const vector<vector <T> > &in,int i,int d1,int d2){
////	return leviCivita(d1,d2)*(in[(this->*p)(i,d1)][d2]-in[i][d2]-in[(this->*p)(i,d2)][d1]+in[i][d1]);
//	return leviCivita(d1,d2)*(in[p_normal(i,d1)][d2]-in[i][d2]-in[p_normal(i,d2)][d1]+in[i][d1]);
//}

int LATTICE::pow(int x,int y){
    int i,out=1;
    for(i=0;i<y;i++){
       out=out*x;
    }
    return out;
}
double LATTICE::pow(double x,int y){
    int i;
	double out=1.0;
    for(i=0;i<y;i++){
       out=out*x;
    }
    return out;
}
int LATTICE::getD(){return DIMS;}
int LATTICE::getN(){return N;}
double LATTICE::mod2pi(double x){
    while(x>pi) x-=2.0*pi;
    while(x<(-1.0*pi)) x+=2.0*pi;
    return x;
}
