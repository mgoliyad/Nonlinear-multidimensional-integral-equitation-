/* 
 * File:   newNonComp.h
 * Author: masha
 *
 * Created on March 5, 2019, 7:01 PM
 */

#ifndef NEWNONCOMP_H
#define	NEWNONCOMP_H


#include <complex>
using namespace std;

extern  double c1; 
extern  double c2; 
extern  double nu1_3;
extern  double nu2_3;
extern  double alpha1;
extern  const double alpha2;
extern const double beta1;
extern const double beta2;
extern double *r;
extern double *w;
extern const double delta;
extern int choiceT;
extern int choiceZ;
extern int choiceX;



double phi1(int , int  );

double psi1(int , int  );

double  Omega1 ( int , int ,int);


double phi2(int , int  );

double psi2(int , int  );

double  Omega2 (int , int ,int);

double R(int , int , int  );

double FuncForXN1 (int , int ,int , int,int );

double FuncForXN2 (int , int ,int , int,int );

double FuncForZN1 (int , int , int );

double FuncForZN2 (int , int , int );

double FuncForTU ( int , int );

double FuncForTW ( int , int );

double callMethodXN(double (*f)(int ,int , int ), int , int ,int , int );

double callMethodZN(double (*f)(int ,int , int ), int , int ,int  );

double callMethodTN(double (*f)(int ,int  ), int , int,char  );


#endif	/* NEWNONCOMP_H */

