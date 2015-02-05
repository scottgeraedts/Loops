#include <lattice.h>

using namespace std;
//contains the functions for the table of the villian potential

//initializes data table
void VTABLE::initVTABLE(double t,int func){

    //this method runs ~60000 steps from -pi to pi
    L=62833;
    step=0.0001;
    xmin=-3.1415;

    //determine which function will be used to evaluate the potential
    if(func==1) eval=&VTABLE::simple;
    else if(func==2) eval=&VTABLE::D2simple;
    else if(func==3) eval=&VTABLE::Dsimple;
    else eval=&VTABLE::energy1;

    double x=xmin;
    V=vector<double>(L,0.0);

    V[0]=(*this.*eval)(t,-1.0*pi); //fill the array with the appropriate values
    V[L-1]=(*this.*eval)(t,pi);
    for(int k=1;k<L-1;k++){
        V[k]=(*this.*eval)(t,x);
        x+=step;
    }
}

//the function v1 that is part of the Villain potential
//change this if you want a different function
double VTABLE::v1(double t,double j){
    return pow(j,2)/t;
}

//returns data value for a given input
double VTABLE::value(double x){
    //mod the value to get it between -pi and pi
    while(x>pi){x-=2.0*pi;}
    while(x<-1.0*pi){x+=2.0*pi;}

    double dx;
    int i=floor((x-xmin)/step)+1;//find the point immediately to the left of x
    if(i<1) {i=0; dx=x-pi;}
    else{
        if(i>L-2) i=L-1;
        dx=x-((i-1)*step+xmin);
    }
 //   cout<<x<<"  "<<i<<"  "<<dx<<"  "<<(V[i+1]-V[i])/step*dx;
    //interpolate a line between the points on either side of x, then calculate x on that line
    return (V[i]+(V[i+1]-V[i])/step*dx);
}

double VTABLE::value(int i){
    return V[i];
}
//creates this VTABLE is the derivative of a previous VTABLE
void VTABLE::differentiate(VTABLE& in){
    step=in.getStep();
    xmin=in.getXmin();
    L=in.getL();
    V=vector<double>(L,0);
    V[0]=(in.value(1)-in.value(0))/step;
    V[L-1]=(in.value(L-1)-in.value(L-2))/step;
   for(int i=1;i<L-1;i++){
        V[i]=(in.value(i+1)-in.value(i-1))/(2.0*step);
    }
}

double VTABLE::getStep(){return step;}
double VTABLE::getXmin(){return xmin;}
double VTABLE::getL(){return L;}
//this is really the only part that needs to be changed to accomodate a different function
double VTABLE::simple(double t,double x){
    double next,oldout=0.0,out=1.0;
    for(double j=1;j>-1;j++){
        next=2.0*exp(-v1(t,j)/2.0);
        out+=next*cos(j*x);
        if(next<macheps){
        	break;
        }
        oldout=out;
    }
//	if(out<macheps) out=macheps;
    return -1.0*log(out);
}
//double derivative of V?
double VTABLE::D2simple(double t,double x){
    double next,top=0.0;
    for(double j=1;j>-1;j++){
        next=2*exp(-v1(t,j)/2.0)*j*j;
        top+=next*cos(j*x);
        if(next<macheps) break;
    }
    double bottom=1.0;
    for(double j=1;j>-1;j++){
        next=2.0*exp(-v1(t,j)/2.0);
        bottom+=next*cos(j*x);
        if(next<macheps) break;
    }
    return top/bottom+Dsimple(t,x)*Dsimple(t,x);
}
//the derivative of simple
double VTABLE::Dsimple(double t,double x){
    double next,top=0.0;
    for(double j=1;j>-1;j++){
        next=2*exp(-v1(t,j)/2.0)*j;
        top+=next*sin(j*x);
        if(next<macheps) break;
    }
    double bottom=1.0;
    for(double j=1;j>-1;j++){
        next=2.0*exp(-v1(t,j)/2.0);
        bottom+=next*cos(j*x);
        if(next<macheps) break;
    }
//	if(bottom<macheps) bottom=macheps;
    return top/bottom;
}
//the derivative of simple
double VTABLE::energy1(double t,double x){
    double top=0.0;
    double next;
    for(double j=1;j>-1;j++){
        next=exp(-v1(t,j)/2.0)*v1(t,j);
        top+=next*cos(j*x);
        if(next<macheps) break;
    }

    return top;
}

