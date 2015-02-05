#ifndef LINEARPROGRAMMING_H
#define LINEARPROGRAMMING_H

#include "MersenneTwister.h"
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <time.h>

using namespace std;

class LinearProgramming{
public:
/*!\class LinearProgramming
 *@brief  uses the simplex algorithm to solve a linear programming problem
 *details of the algorithm can be found in Numerical Recipies, section 10.8
 *this codes essentially comes from there as well
 */
    const static double EPS=1.0e-6;
    static void simplx(double* a,int m,int n,int m1,int m2,int m3,int *icase,int izrov[],int iposv[]);
private:
    static void simp1(const double *a,int mm,const int ll[],int nll,int iabf,int *kp,double *bmax,int n);
    static void simp2(const double *a,int m,int n,int *ip,int kp);
    static void simp3(double *a,int i1,int k1,int ip,int kp);
    static int getA(int,int,int);
    
};

#endif
