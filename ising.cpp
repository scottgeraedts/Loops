#include "ising.h"

ISING::ISING(map<string,double> &params):LATTICE( params){
	
	J=params["J"];
	s=vector<int>(N,1);
	double r;
	for (int i=0;i<N;i++){
		r=ran.rand();
		if (r<0.5) s[i]=-1;
	}
	updateFuncs.push_back(&ISING::updateSpin);
	updateFuncs.push_back(&ISING::updateWolff);
	updateWeights.push_back(1.0*N);
	updateWeights.push_back(1.0*L);
	runningE=energy(); //need to call this in every derived member of lattice (because it needs to be run after everything is instantiated)
	drawn=0;
}
ISING::ISING(){}
double ISING::energy(){
	double out=0.0;
	for (int i=0;i<N;i++){
		out+=h*s[i];
		for (int d=0;d<DIMS;d++)
			out+=-J*s[i]*s[(this->*p)(i,d)];
	}
	return out;
}		
int ISING::updateSpin(int i){
	i=i/DIMS;
	double dE=0.0;

	dE+=-h*s[i];
	for(int d=0;d<DIMS;d++)
		dE+=2*J*s[i]*(s[(this->*p)(i,d)]+s[(this->*m)(i,d)]);

	int accept=0;
	if(dE<0) accept=1;
	else{
		double r=ran.rand();
		if(r<exp(-dE)) accept=1;
	}
	if(accept){
		s[i]=-s[i];
		runningE+=dE;
	}
	return accept;
}
int ISING::updateWolff(int i){
	i=i/DIMS;
	int site,count=0;
	double r;
	vector<int> flipped(N,0);
	stack<int> edge;
	double dE=0.0;

	s[i]=-s[i];
	edge.push(i);
	flipped[i]=1;
	while(!edge.empty()){
		site=edge.top();
		edge.pop();
		for(unsigned d=0;d<DIMS;d++){
			if(s[site]==s[(this->*p)(site,d)]) continue;
			if(flipped[(this->*p)(site,d)]) continue;
			r=ran.rand();
			if(r> (1-exp(-2.*J)) ) continue;
			s[(this->*p)(site,d)]=-s[(this->*p)(site,d)];
			flipped[(this->*p)(site,d)]=1;
			edge.push((this->*p)(site,d));
			count++;
		}	
		for(unsigned d=0;d<DIMS;d++){
			if(s[site]==s[(this->*m)(site,d)]) continue;
			if(flipped[(this->*m)(site,d)]) continue;		
			r=ran.rand();
			if(r> (1-exp(-2.*J)) ) continue;
			s[(this->*m)(site,d)]=-s[(this->*m)(site,d)];
			flipped[(this->*m)(site,d)]=1;
			edge.push((this->*m)(site,d)); 
			count++;
		}
	}
	for(unsigned j=0;j<N;j++){
		if (flipped[j]){
			for(unsigned d=0;d<DIMS;d++){
				if(!flipped[(this->*p)(j,d)]) dE+=-2*J*s[j]*s[(this->*p)(j,d)];
				if(!flipped[(this->*m)(j,d)]) dE+=-2*J*s[j]*s[(this->*m)(j,d)];
			}
		}
	}
	runningE+=dE;
	cout<<dE<<endl;
	return count;
}	

double ISING::magnetization(){
	double out=0.0;
	for(int i=0;i<N;i++)
		out+=s[i];
	return (1.0*out)/(1.0*N);
}
//creates a snapshot of the edge using pstricks
void ISING::draw(){

	double step=1;
	if (DIMS>2){
		cout<<"draw not supported to >2 dimensions!"<<endl;
		exit(0);
	}
	//open file and write header
	ofstream drawout;
	if (!drawn){	
		drawout.open("drawing.tex");
		drawn=1;
		drawout.precision(8);
		drawout<<"\\documentclass[letterpaper,10pt]{article}"<<endl;
		drawout<<"\\usepackage{amsmath}"<<endl;
		drawout<<"\\usepackage{pstricks}"<<endl;
		drawout<<"\\usepackage{pstricks-add}"<<endl;
	//  sketch puts this in the header but im not sure that we need it
		drawout<<"\\oddsidemargin 0in"<<endl;
		drawout<<"\\evensidemargin 0in"<<endl;
		drawout<<"\\topmargin 0in"<<endl;
		drawout<<"\\headheight 0in"<<endl;
		drawout<<"\\headsep 0in"<<endl;
		drawout<<"\\textheight 9in"<<endl;
		drawout<<"\\textwidth 6.5in"<<endl;
		drawout<<"\\begin{document}"<<endl;
		drawout<<"\\pagestyle{empty}"<<endl;
		drawout<<"\\vspace*{\\fill}"<<endl;

		//border
//		drawout<<"\\psline[linecolor=black](0,0,0)("<<-0.5*step<<","<<(L-0.5)*step<<")\\psline[linecolor=black](0,0)("<<L*(step-0.5)<<","<<-0.5*step<<")"<<endl;
//		drawout<<"\\psline[linecolor=black]("<<-0.5*step<<","<<L*(step-0.5)<<")("<<L*(step-0.5)<<","<<L*(step-0.5)<<")"<<endl;
//		drawout<<"\\psline[linecolor=black]("<<L*(step-0.5)<<","<<-0.5*step<<")("<<L*(step-0.5)<<","<<L*(step-0.5)<<")"<<endl;
	}else{
		drawout.open("drawing.tex",ofstream::app);
		drawout<<"\\newpage"<<endl;
	}
			
	drawout<<"\\begin{center}"<<endl;
	drawout<<"\\begin{pspicture}(-.2,-.2)("<<(L+1)*step<<","<<(L+1)*step<<")"<<endl;
	drawout<<"\\pstVerb{1 setlinejoin}"<<endl;
	int x,y;
	string color;
    for(int i=0;i<N;i++){
    	y=(i/L);
    	x=i%L;
		if (s[i]==1) color="green";
		else color="red";
		drawout<<"\\psdot[linecolor="<<color<<",dotsize=10pt]("<<x*step<<","<<y*step<<")"<<endl;    	
    }
	drawout<<"\\end{pspicture}"<<endl;
	drawout<<"\\end{center}"<<endl;
}
ISING::~ISING(){
	if(drawn){
		ofstream drawout;
		drawout.open("drawing.tex",ofstream::app);
		//footer for draw file, if necessary
		drawout<<"\\vspace*{\\fill}"<<endl;
		drawout<<"\\end{document}"<<endl;
		drawout.close();
	}	
}
