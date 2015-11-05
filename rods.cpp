#include "rods.h"

RODS::RODS(map<string,double> &params):LATTICE( params){
	
	if (DIMS!=3){
		cout<<"wrong number of dimensions!"<<endl;
		exit(0);
	}
	t1=params["t1"];
	t2=params["t2"];
	lambda1=params["lambda1"];
	lambda2=params["lambda2"];
	edge=params["edge"];
	phi=vector<double>(N,0);
	vector<int> temp3(3,0);
	a=vector<vector<int> >(N,temp3);
	pp=vector<vector<int> >(N,temp3);
	gamma=vector<double>(3,0);
	
	//add functons and their weights to the function list
	updateFuncs.push_back(&RODS::updatePhi);
	updateFuncs.push_back(&RODS::updateA);
	updateFuncs.push_back(&RODS::updatePP);
	updateFuncs.push_back(&RODS::compAP);
	if(edge) updateFuncs.push_back(&RODS::starUpdate);
	updateFuncs.push_back(&RODS::updateGamma);
	updateWeights.push_back(1.0*N);//phi
	updateWeights.push_back(1.0*N);//a
	updateWeights.push_back(1.0*N);//pp
	updateWeights.push_back(1.0*N);//pp+a
	if(edge) updateWeights.push_back(1.0*N);//star
	updateWeights.push_back(3.0);  //gamma
	
	//zero the correlators. correlators are labelled [start][direction][distance]
	vector< complex<double> > tempc(L,0);
	vector< vector <complex<double> > > temp2(3,tempc);
	cors=vector< vector< vector< complex<double> > > > (L,temp2);
	avgphi=vector< complex<double> >(L,0);
	
	angle_step=1.0;
	theta=1.0;
	runningE=energy(); //need to call this in every derived member of lattice (because it needs to be run after everything is instantiated
	
}
RODS::RODS(){};
double RODS::energy(){
	double out=0.0;
	int JminusQ;
	for (int i=0;i<N;i++){
		for (int d=0;d<DIMS;d++){
			out+=0.5*t1*pow(cosTerm(i,d),2);
			out-=lambda1*( cos( 2*(phi[i]-phi[(this->*p)(i,d)]) ) + cos( 2*(phi[i]-phi[(this->*m)(i,d)]) ) );
			for(int d2=d+1;d2<3;d2++){
				JminusQ=curl(a,i,d,d2)-curl(pp,i,d,d2);
				out+=1.0/(2.0*t2)*pow(JminusQ,2);
				out-=lambda2*cos(pi*JminusQ);
			}
		}
	}
	return out;
}

int RODS::updatePhi(int i){
	i=i/3;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	double step;
//	if (r<0.8) step=r*2.0*angle_step-angle_step;
//	else step=(1-r)*2.0*angle_step-angle_step+pi;
	step=r*2.0*angle_step-angle_step;
	
	double newphi=mod2pi(phi[i]+step);

	for(int d=0;d<3;d++){
		oldE+=0.5*t1*(pow(cosTerm(i,d),2)+pow(cosTerm((this->*m)(i,d),d), 2) );
		oldE-=lambda1*( cos( 2*(phi[i]-phi[(this->*p)(i,d)]) ) + cos( 2*(phi[i]-phi[(this->*m)(i,d)]) ) );
		newE+=0.5*t1*(pow(newphi-phi[(this->*p)(i,d)]-2*pi*pp[i][d]+end[i][d]*gamma[d],2)+pow(phi[(this->*m)(i,d)]-newphi-2*pi*pp[(this->*m)(i,d)][d]+end[(this->*m)(i,d)][d]*gamma[d], 2) );
		newE-=lambda1*( cos( 2*(newphi-phi[(this->*p)(i,d)]) ) + cos( 2*(newphi-phi[(this->*m)(i,d)]) ) );
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
int RODS::updateA(int in){
	int site=in/3;
	int d=in%3;
	if( edge && ( site%L==0 || site%L==L/2) && (d==1 || d==2) ) return 0;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	int step=1,JminusQ;
	if (r<0.25) step=-2;
	else if(r<0.5) step=-1;
	else if(r<0.75) step=1;
	else step=2;
	
	for (int d2=1;d2<3;d2++){
		JminusQ=curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3);
		oldE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
		JminusQ=curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3);
		oldE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
	}	
	a[site][d]+=step;	
	for (int d2=1;d2<3;d2++){
		JminusQ=curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3);
		newE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
		JminusQ=curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3);
		newE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
	}	
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
int RODS::updatePP(int in){
	int site=in/3;
	int d=in%3;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	int step=1,JminusQ;
	if (r<0.5) step=-1;
	oldE+=0.5*t1*pow(cosTerm(site,d),2);
	for (int d2=1;d2<3;d2++){
		JminusQ=curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3);
		oldE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
		JminusQ=curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3);
		oldE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
	}	
	pp[site][d]+=step;	
	newE+=0.5*t1*pow(cosTerm(site,d),2);
	for (int d2=1;d2<3;d2++){
		JminusQ=curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3);
		newE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
		JminusQ=curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3);
		newE+=1.0/(2.0*t2)*pow(JminusQ,2)-lambda2*cos(pi*JminusQ);
	}	
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept) runningE+=newE-oldE;
	else pp[site][d]-=step;
	return accept;
}
int RODS::compAP(int in){
	int site=in/3;
	int d=in%3;
	if( edge && ( site%L==0 || site%L==L/2) && (d==1 || d==2) ) return 0;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	int step=1;
	if (r<0.5) step=-1;
	oldE+=0.5*t1*pow(cosTerm(site,d),2);
//	for (int d2=1;d2<3;d2++)
//		oldE+=1.0/(2.0*t2)*pow(curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3),2)+1.0/(2.0*t2)*pow(curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3),2);
	pp[site][d]+=step;	
	a[site][d]+=step;
	newE+=0.5*t1*pow(cosTerm(site,d),2);
//	for (int d2=1;d2<3;d2++)
//		newE+=1.0/(2.0*t2)*pow(curl(a,site,d,(d+d2)%3)-curl(pp,site,d,(d+d2)%3),2)+1.0/(2.0*t2)*pow(curl(a,(this->*m)(site,(d+d2)%3),d,(d+d2)%3)-curl(pp,(this->*m)(site,(d+d2)%3),d,(d+d2)%3),2);
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept) runningE+=newE-oldE;
	else{
		pp[site][d]-=step;
		a[site][d]-=step;
	}
	return accept;
}
//updates pp and a in a cross pattern, which allows pairs of vortices to cross the boundary
int RODS::starUpdate(int in){
	int i=in/3;
	int d=in%3;
	double oldE=0.0, newE=0.0;
	double r=ran.rand();
	int step;
	if (r<0.25) step=-2;
	else if (r<0.5) step=-1;
	else if (r<0.75) step=1;
	else step=2;
	int site,sign;
	for(int d2=0;d2<2;d2++){
		for(int k=0;k<2;k++){
			if(k==0){
				site=i;
				sign=1;
			}
			else{
				site=(this->*m)(i, (d+d2)%3);
				sign=-1;
			}
			oldE+=0.5*t1*pow(cosTerm(site,(d+d2)%3),2);
			pp[site][(d+d2)%3]+=sign*step;
			a[site][(d+d2)%3]+=sign*step;
			newE+=0.5*t1*pow(cosTerm(site,(d+d2)%3),2);
		}
	}
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept){
		runningE+=newE-oldE;
	}else{
		for(int d2=0;d2<2;d2++){
			for(int k=0;k<2;k++){
				if(k==0){
					site=i;
					sign=1;
				}
				else{
					site=(this->*m)(i, (d+d2)%3);
					sign=-1;
				}
				pp[site][(d+d2)%3]-=sign*step;
				a[site][(d+d2)%3]-=sign*step;
			}	
		}
	}
	return accept;	
}				

int RODS::updateGamma(int in){
	int d=in%3;
	double oldE=0.0, newE=0.0;
	double step=ran.rand(2.0*angle_step)-angle_step;
	double newgamma=mod2pi(gamma[d]+step);
	for(int i=0;i<N;i++){
		oldE+=0.5*t1*pow(cosTerm(i,d),2);
		newE+=0.5*t1*pow(phi[i]-phi[(this->*p)(i,d)]-2*pi*pp[i][d]+end[i][d]*newgamma,2);
	}
	int accept=0;
	if(newE<oldE) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(oldE-newE)) accept=1;
	}
	if(accept){
		gamma[d]=newgamma;
		runningE+=newE-oldE;
	}
	return accept;	
}
double RODS::cosTerm(int i,int d){ return phi[i]-phi[(this->*p)(i,d)]-2*pi*pp[i][d]+end[i][d]*gamma[d]; }
double RODS::magnetization(){
	double x=0.0,y=0.0;
	for(int i=0;i<N;i++){
		x+=cos(phi[i]);
		y+=sin(phi[i]);
	}
	x=x/(1.0*N);
	y=y/(1.0*N);
	return (pow(x,2)+pow(y,2));
}
double RODS::Z2magnetization(){
	complex<double> out(0,0);
	for(int i=0;i<N;i++){
		out+=polar(1.,2*phi[i]);
	}
	out/=(1.0*N);
	return abs(out);
}	
double RODS::rho(){
	double out;
	for(unsigned d1=0;d1<3;d1++){
		for(unsigned d2=d1+1;d2<3;d2++)
			out+=rhoHelper(a, d1,d2);
	}
	return out/6.0;
}
double RODS::Z2rho(){
	double out;
    double angle;
    double Lu;
    complex<double> px=0.0;
	int pos;
    int u1,u2;
	for(unsigned d1=0;d1<3;d1++){
		for(unsigned d2=d1+1;d2<3;d2++){
			out+=rhoHelper(a, d1,d2);
			u1=(d1+1)%3;
			u2=(d1+2)%3;
			for(int i=0;i<N;i++){
				pos=(this->*p)(i,d1);
				Lu=cos(pi*curl(a,pos,u1,u2));
				angle=(2.0*pi/(1.0*L))*((int)(i/pow(L,d2))%L);
				px+=Lu*polar(1.,angle);
			}
			out+=abs(px)/(1.0*N);
		}
	}
	return out/6.0;
}
double RODS::QQ(){
	double out;
	for(unsigned d1=0;d1<3;d1++){
		for(unsigned d2=d1+1;d2<3;d2++)
			out+=rhoHelper(pp, d1,d2);
	}
	return out/6.0;
}
double RODS::rhoHelper(const vector<vector <int> > &in, int m, int n){
//calculates curl of A in a given direction m, wrt to a minimum q in dimension n
    double angle;
    double Lu;
    double px=0.0, py=0.0;
	int pos;
    int u1=(m+1)%3;
    int u2=(m+2)%3;
    for(int i=0;i<N;i++){
		pos=(this->*p)(i,m);
		Lu=curl(in,pos,u1,u2);
        angle=(2.0*pi/(1.0*L))*((int)(i/pow(L,n))%L);
        px+=Lu*cos(angle);
	    py+=Lu*sin(angle);
    }
    return (py*py+px*px)/(1.0*N);
}
void RODS::updateCorrelators(){
	int i;
	for(int start=0;start<L;start++){
		for(int k=0;k<N;k+=L){
			avgphi[start]+=polar(1.,shiftedphi(start+k));
			for(int d=0;d<3;d++){
				i=start+k;
				for(int j=0;j<L;j++){
					cors[start][d][j]+=polar(1.,shiftedphi(i)-shiftedphi(start+k));
					i=(this->*p)(i,d);
				}
			}
		}
	}
}
void RODS::printCorrelators(int NROD){
	ofstream cfout;
	stringstream filename;
	for(int start=0;start<L;start++){
		filename.str("");
		filename<<start<<"cors";
		cfout.open(filename.str().c_str());
		for(int i=0;i<L;i++){
			cfout<<i<<" ";
			for(int d=0;d<3;d++) cfout<<cors[start][d][i].real()/(1.*NROD*L*L)<<" ";//<<cors[start][d][i].imag()/(1.*NROD*L*L)<<" ";
			cfout<<pow(abs(avgphi[start]/(1.*NROD*L*L)),2)<<" ";
			cfout<<endl;
		}
		cfout.close();
	}
}
double RODS::shiftedphi(int i){ 
        int x=i%L;
        int y=(i/L)%L;
        int z=i/(L*L);
        double out=phi[i]+1.0*x*gamma[0]/(1.0*L)+1.0*y*gamma[1]/(1.0*L)+1.0*z*gamma[2]/(1.0*L);
        return out;
}

