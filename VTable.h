/*
 * Header for Potential for Villian Model
 */

#ifndef _VTABLE_H
#define	_VTABLE_H

#include <math.h>
#include <iostream>

const double pi=3.141592654;
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
    ~VTABLE();

private:
    int L;
    double step;
    double xmin;
    double * V;
    double (VTABLE::*eval)(double,double);
    double simple(double,double);
    double energy1(double,double);
    double Dsimple(double,double);
    double D2simple(double,double);
};

#endif	/* _VTABLE_H */

