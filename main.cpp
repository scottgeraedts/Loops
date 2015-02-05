/* 
 * Created on January 29, 2011, 8:57 PM
 * Used for two species of loops with 1/r^2 interactions
 */

#include "ising.h"

using namespace std;

int main(int argc, char** argv) {
	
    ifstream cfin;
    cfin.open("params");
	map<string,double> params;

	cfin>>params["L"];
	cfin>>params["NWarmUps"];
	cfin>>params["NMeas"];
	cfin>>params["NSteps"];
	cfin>>params["seed"];

	//get model-dependent parameters	
	cfin>>params["J"];
	cfin.close();

	params["DIMS"]=2;
	map<string,double> counter;
	counter["E"]=0;
	counter["E2"]=0;

	//the things you want to print
	counter["M"]=0;
	counter["M2"]=0;

	ofstream hist;
	ostringstream filename;

	//choose how to name your hist file
	filename<<params["J"]<<"_"<<"_"<<params["L"];
	hist.open(filename.str().c_str());

	//create lattice, choose what kind you want
	UPDATER<ISING> sim(params);

	sim.random_update(params["NWarmUps"]);
	double e,m;
	double steps=1.0*params["NMeas"];
	for(int j=0;j<params["NMeas"];j++){
		sim.random_update(params["NSteps"]);
		e=sim.lat.E();
		counter["E"]+=e/steps;
		counter["E2"]+=e*e/steps;
		m=abs(sim.lat.magnetization());
		counter["M"]+=m/steps;
		counter["M2"]+=m*m/steps;
		//make your measurements
//        for(int m=0;m<3;m++){
//           for(int n=1;n<3;n++){
//              	counter["JJ"]+=lat.rho(m,(m+n)%3)/(6.0*steps);
//           }
//        }

		hist<<e<<" "<<sim.lat.energy()<<" "<<m<<endl;//comment this our when you're done checking
	}
	sim.lat.draw();

	ofstream out;
	out.open("out");
    out<<"# Loop lattice model for L="<<params["L"]<<" with "<<params["NWarmUps"]<<" warmup sweeps, "<<params["NMeas"]<<" measurements and "<<params["NSteps"]<<" sweeps between measurements"<<endl;
	out<<"#J";
	for (map<string,double>::iterator it=counter.begin(); it!=counter.end();++it)
		out<<it->first<<" ";
	out<<endl;
	out<<params["J"]<<" ";
	for (map<string,double>::iterator it=counter.begin(); it!=counter.end();++it)
		out<<it->second<<" ";
	out<<endl;
	out.close();
	sim.rate();
	hist.close();

}
