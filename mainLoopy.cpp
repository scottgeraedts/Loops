/* 
 * Created on January 29, 2011, 8:57 PM
 * Used for two species of loops with 1/r^2 interactions
 */

#include "rods.h"

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
	cfin>>params["t1"];
	cfin>>params["t2"];
	cfin>>params["lambda1"];
	cfin>>params["lambda2"];
	cfin.close();

	map<string,double> counter;
	counter["E"]=0;
	counter["E2"]=0;

	//the things you want to print
	counter["JJ"]=0;
	counter["QQ"]=0;
	counter["M"]=0;
	counter["Z2M"]=0;

	ofstream hist;
	ostringstream filename;

	//choose how to name your hist file
	filename<<params["t1"]<<"_"<<params["t2"]<<"_"<<params["L"];
	hist.open(filename.str().c_str());

	//create lattice, choose what kind you want
	UPDATER<RODS> sim(params);

	sim.random_update(params["NWarmUps"]);
	double e;
	double steps=1.0*params["NMeas"];
	for(int j=0;j<params["NMeas"];j++){
		sim.random_update(params["NSteps"]);
		sim.lat.updateCorrelators();
		e=sim.lat.E();
		counter["E"]+=e/steps;
		counter["E2"]+=e*e/steps;
		counter["M"]+=sim.lat.magnetization()/steps;
		counter["JJ"]+=sim.lat.rho()/steps;
		counter["QQ"]+=sim.lat.rho()/steps;
		counter["Z2M"]+=sim.lat.Z2magnetization()/steps;
		//make your measurements
//        for(int m=0;m<3;m++){
//           for(int n=1;n<3;n++){
//              	counter["JJ"]+=lat.rho(m,(m+n)%3)/(6.0*steps);
//           }
//        }

		hist<<e<<" "<<sim.lat.energy()<<endl;//comment this our when you're done checking
	}
	sim.lat.printCorrelators(params["NMeas"]);
	ofstream out;
	out.open("out");
    out<<"# Loop lattice model for L="<<params["L"]<<" with "<<params["NWarmUps"]<<" warmup sweeps, "<<params["NMeas"]<<" measurements and "<<params["NSteps"]<<" sweeps between measurements"<<endl;
	out<<"#t1 t2";
	for (map<string,double>::iterator it=counter.begin(); it!=counter.end();++it)
		out<<it->first<<" ";
	out<<endl;
	out<<params["t1"]<<" "<<params["t2"]<<" ";
	for (map<string,double>::iterator it=counter.begin(); it!=counter.end();++it)
		out<<it->second<<" ";
	out<<endl;
	out.close();
	sim.rate();
	hist.close();

}
