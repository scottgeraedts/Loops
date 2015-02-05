/** \file LinearProgramming.cc
 * @brief contains the functions used to solve a linear programming problem 
 * via the simplex method
 * @author Numerical Recipies, sec 10.8
 * @data 2009-07-31
 */
#include "LinearProgramming.h"
#include <iostream>

void LinearProgramming::simplx(double *a,int m,int n,int m1,int m2,int m3,int *icase,int izrov[],int iposv[]){
    int ip,is,k,kh,kp,nl1;
    int *l1=new int[n+1];
    int *l3=new int[m];
    double q1,bmax; 
    nl1=n;
    
    //the three different types of constraints should add up to equal the total number of constraints
    if(m!=(m1+m2+m3)){
        cout<<"Bad input constraint count in simplx"<<endl;
        *icase=2;
        delete [] l1; delete [] l3;
        return;
    }
    for(k=1;k<=n;k++){ l1[k-1]=k; //initialize list of columns available for exchange, and make all variables right-hand
    	izrov[k-1]=k;
    }
    
    //the first column should have all positive numbers
    for(int i=1;i<=m;i++){
        if(a[getA(i+1,1,n)]<0.0){
            cout<<"Bad input tableaux in simplx "<<getA(i+1,1,n)<<endl;
            *icase=2;
            delete [] l1; delete [] l3;
            return;
        }
        iposv[i-1]=n+i;
    }
    //phase one: generate a vector that satisfies the constraints
    if(m2+m3){//if only < constraints, the origin satisfies the constraints and we can skip this step
        for(int i=1;i<=m2;i++) l3[i-1]=1; //list of m2 constraints that havent been exchanged out yet
        
        for(k=1;k<=(n+1);k++){//compute auxiliary objective function
            q1=0.0;
            for(int i=m1+1;i<=m;i++) q1+=a[getA(i+1,k,n)];
            a[getA(m+2,k,n)]=-q1;
        }
        while(1){
            simp1(a,m+1,l1,nl1,0,&kp,&bmax,n);//find max coeff of auxiliary objective function
        /*cout<<"----------------"<<endl;
        for(int j=1;j<=m+2;j++){
            for(int l=1;l<=n+1;l++){
                printf("% 6.4f ",a[getA(j,l,n)]);
            }
            cout<<endl;
        }*/

            if(bmax<=EPS && a[getA(m+2,1,n)]<-EPS){//auxiliary function is negative and cant be improved, no solution exists
                *icase=-1;
                cout<<"aux func neg "<<getA(m+2,1,n)<<" "<<bmax<<endl;
                delete [] l1; delete [] l3;
                return;
            }else if(bmax<=EPS&&a[getA(m+2,1,n)]<=EPS){//the auxiliary objective function is a good starting vector, clean up at 
                //phase one and then continue to phase two
                for(ip=m1+m2+1;ip<=m;ip++){
                    if(iposv[ip-1]==(ip+n)){
                        simp1(a,ip,l1,nl1,1,&kp,&bmax,n);//find pivot column
                        if(bmax<EPS) goto one;
                    }
                }
                for(int i=m1+1;i<=m1+m2;i++)
                    if(l3[i-m1-1]==1){
                        for(k=1;k<=n+1;k++)
                            a[getA(i+1,k,n)]=-a[getA(i+1,k,n)];
                    }
                break;
            }
            simp2(a,m,n,&ip,kp);//find pivot element
            if(ip==0){//solution is unbounded
                *icase=-1;
                cout<<"unbounded in phase 1"<<endl;
                delete [] l1; delete [] l3;
                return;
            }
one:    simp3(a,m+1,n,ip,kp);//exchange the rows and columns of list
            if(iposv[ip-1]>=(n+m1+m2+1)){//deal with m3 constraints
                for(k=1;k<=nl1;k++)
                    if(l1[k-1]==kp) break;
                --nl1;
                for(is=k;is<=nl1;is++) l1[is-1]=l1[is];
            }else{//compensate for < constraint
                kh=iposv[ip-1]-m1-n;
                if(kh>=1 && l3[kh-1]){
                    l3[kh-1]=0;
                    ++a[getA(m+2,kp+1,n)];
                    for(int i=1;i<=m+2;i++)
                        a[getA(i,kp+1,n)]=-a[getA(i,kp+1,n)];
                }
            }
            //cout<<"switch in phase 1 "<<ip<<" "<<kp<<endl;
            is=izrov[kp-1];
            izrov[kp-1]=iposv[ip-1];
            iposv[ip-1]=is;
        }//loop
    }//end of phase one
    //phase two: optimize the solution
    while(1){
        simp1(a,0,l1,nl1,0,&kp,&bmax,n);//is the initial equation all zero?
        if(bmax<=EPS){//good solution found
            *icase=0;
            delete [] l1; delete [] l3;
            return;
        }
        simp2(a,m,n,&ip,kp);
        if(ip==0){//solution is unbounded
            *icase=1;
            cout<<"unbounded"<<endl;
            delete [] l1;delete[]l3;
            return;
        }
            //cout<<"switch in phase 2 "<<ip<<" "<<kp<<endl;
        simp3(a,m,n,ip,kp); //switch the rows
        is=izrov[kp-1];
        izrov[kp-1]=iposv[ip-1];
        iposv[ip-1]=is;
    }//end of phase two
}
void LinearProgramming::simp1(const double *a,int mm,const int ll[],int nll,int iabf,int *kp,double *bmax,int n){
//finds maximum of elements in the row mm, columns ll
    double test;
    
    if(nll<=0)
        *bmax=0.0;
    else{
        *kp=ll[0];//get the first value
        *bmax=a[getA(mm+1,*kp+1,n)];
        for(int k=2;k<=nll;k++){//get the next value, if larger than the previous largest value exchange them
            if(iabf==0)
                test=a[getA(mm+1,ll[k-1]+1,n)]-(*bmax);
            else
                test=fabs(a[getA(mm+1,ll[k-1]+1,n)])-fabs(*bmax);
            if(test>0.0){
                *bmax=a[getA(mm+1,ll[k-1]+1,n)];
                *kp=ll[k-1];
            }
        }
    }
}
void LinearProgramming::simp2(const double *a,int m,int n,int *ip,int kp){
    int i;
    double qp,q0,q,q1;
    
    *ip=0;
    for(i=1;i<=m;i++)
        if(a[getA(i+1,kp+1,n)]<-EPS) break;
    if(i>m) return;
    q1=-a[getA(i+1,1,n)]/a[getA(i+1,kp+1,n)];
    *ip=i;
    for(i=*ip+1;i<=m;i++){
        if(a[getA(i+1,kp+1,n)]<-EPS){
            q=-a[getA(i+1,1,n)]/a[getA(i+1,kp+1,n)];
            if(q<q1){
                *ip=i;
                q1=q;
            }else if(q==q1){
                for(int k=1;k<=n;k++){//we have a degeneracy
                    qp=-a[getA(*ip+1,k+1,n)]/a[getA(*ip+1,kp+1,n)];
                    q0=-a[getA(i+1,k+1,n)]/a[getA(i+1,kp+1,n)];
                    if(q0!=qp) break;
                }
                if(q0<qp) *ip=i;
            }
        }
    }
}
void LinearProgramming::simp3(double *a,int i1,int k1,int ip,int kp){
    int kk,ii;
    double piv;
    
    piv=1.0/a[getA(ip+1,kp+1,k1)];
    for(ii=1;ii<=i1+1;ii++)
        if(ii-1!=ip){
            a[getA(ii,kp+1,k1)]*=piv;
            for(kk=1;kk<=k1+1;kk++)
                if(kk-1!=kp)
                    a[getA(ii,kk,k1)]-=a[getA(ip+1,kk,k1)]*a[getA(ii,kp+1,k1)];
        }
    for(kk=1;kk<=k1+1;kk++)
        if(kk-1!=kp) a[getA(ip+1,kk,k1)]*=-piv;
    a[getA(ip+1,kp+1,k1)]=piv;
}
int LinearProgramming::getA(int y,int x,int n){
    return (y-1)*(n+1)+(x-1);
}
