/* 
 * File:   parameters.h
 * Author: masha
 *
 * Created on March 14, 2019, 7:43 PM
 */

#ifndef PARAMETERS_H
#define	PARAMETERS_H
#include <complex>

using namespace std;

using std::conj;
double c1;
double c2;
double nu1;
double alpha1;
double nu2;
extern int Interval_t;
extern const complex<double> _Complex_I(0.0,1.0); 

extern const double alpha2=0.0;
extern const complex <double> ialpha2_2=0.0;
extern complex<double> ialpha1_2;
extern const double beta1=0.0;
extern const double beta2=0.0;
extern const int num=10;
//const double error_r=0.000001;
//const double error_w=0.00000001;

extern double nu1_3;
extern double nu2_3;

extern complex<double> *A;

extern const double delta=0.05;

extern int max_stepsX=1;
extern int max_stepsZ=1;
extern int max_stepsT=1;

extern int choiceT;
extern int choiceZ;
extern int choiceX;
extern int choiceA;

#endif	/* PARAMETERS_H */

