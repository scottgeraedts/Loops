#include "loopy.h"

LOOPY::LOOPY(map<string,double> &params):LATTICE( params){
	
	if (DIMS!=3){
		cout<<"wrong number of dimensions!"<<endl;
		exit(0);
	}
	t1=params["t1"];
	t2=params["t2"];
	phi=vector<double>(N,0);
	vector<int> temp3(3,0);
	a=vector<vector<int> >(N,temp3);
	gamma=vector<double>(3,0);
	updateFuncs.push_back(&LOOPY::updatePhi);
	updateFuncs.push_back(&LOOPY::updateA);
	updateFuncs.push_back(&LOOPY::updateGamma);
	updateWeights.push_back(1.0*N);
	updateWeights.push_back(3.0*N);
	updateWeights.push_back(3.0);
	angle_step=1.0;
	theta=1.0/3.0;
	table.initVTABLE(t1,1);
	runningE=energy(); //need to call this in every derived member of lattice (because it needs to be run after everything is instantiated
	
}
LOOPY::LOOPY(){};
double LOOPY::energy(){
	double out=0.0;
	for (int i=0;i<N;i++){
		for (int d=0;d<DIMS;d++){
			out+=table.value(phi[i]-phi[(this->*p)(i,d)]-2*pi*theta*a[i][d]+end[i][d]*gamma[d]);
			for(int d2=d+1;d2<3;d2++)
				out+=1.0/(2.0*t2)*pow(curl(a,i,d,d2),2);
		}
	}
	return out;
}

int LOOPY::updatePhi(int i){
	i=i/3;
	double oldE=0.0, newE=0.0;
	double step=ran.rand(2.0*angle_step)-angle_step;

	for(int d=0;d<3;d++){
		oldE+=table.value(cosTerm(i,d))+table.value(cosTerm((this->*m)(i,d),d));
		newE+=table.value(cosTerm(i,d)+step )+table.value(cosTerm((this->*m)(i,d),d)-step) ;
	}
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept){
		phi[i]=mod2pi(phi[i]+step);
		runningE+=newE-oldE;
	}
	return accept;
}
int LOOPY::updateA(int in){
	int site=in/3;
	int d=in%3;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	int step=1;
	if (r<0.5) step=-1;
	oldE+=table.value(cosTerm(site,d));
	for (int d2=1;d2<3;d2++)
		oldE+=1.0/(2.0*t2)*pow(curl(a,site,d,(d+d2)%3),2)+1.0/(2.0*t2)*pow(curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3),2);
	a[site][d]+=step;	
	newE+=table.value(cosTerm(site,d));
	for (int d2=1;d2<3;d2++)
		newE+=1.0/(2.0*t2)*pow(curl(a,site,d,(d+d2)%3),2)+1.0/(2.0*t2)*pow(curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3),2);
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept) runningE+=newE-oldE;
	else a[site][d]-=step;
	return accept;
}

int LOOPY::updateGamma(int in){
	int d=in%3;
	double oldE=0.0, newE=0.0;
	double step=ran.rand(2.0*angle_step)-angle_step;
	for(int i=0;i<N;i++){
		oldE+=table.value(cosTerm(i,d));
		newE+=table.value(cosTerm(i,d)+end[i][d]*step);
	}
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept){
		gamma[d]=mod2pi(gamma[d]+step);
		runningE+=newE-oldE;
	}
	return accept;	
}
double LOOPY::cosTerm(int i,int d){ return phi[i]-phi[(this->*p)(i,d)]-2*pi*theta*a[i][d]+end[i][d]*gamma[d]; }
double LOOPY::magnetization(){
	double x=0.0,y=0.0;
	for(int i=0;i<N;i++){
		x+=cos(phi[i]);
		y+=sin(phi[i]);
	}
	x=x/(1.0*N);
	y=y/(1.0*N);
	return (pow(x,2)+pow(y,2));
}

double LOOPY::rho(){
	double out;
	for(unsigned d1=0;d1<3;d1++){
		for(unsigned d2=d1+1;d2<3;d2++)
			out+=rhoHelper(d1,d2);
	}
	return out/6.0;
}
double LOOPY::rhoHelper(int m,int n){
//calculates curl of A in a given direction m, wrt to a minimum q in dimension n
    double angle;
    double Lu;
    double px=0.0, py=0.0;
	int pos;
    int u1=(m+1)%3;
    int u2=(m+2)%3;
    for(int i=0;i<N;i++){
		pos=(this->*p)(i,m);
		Lu=curl(a,pos,u1,u2);
        angle=(2.0*pi/(1.0*L))*((int)(i/pow(L,n))%L);
        px+=Lu*cos(angle);
	    py+=Lu*sin(angle);
    }
    return (py*py+px*px)/(1.0*N);
}
