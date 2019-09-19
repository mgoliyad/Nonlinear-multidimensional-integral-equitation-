/* 
 * File:   Romb.h
 * Author: masha
 *
 * Created on February 26, 2019, 1:29 AM
 */

#ifndef ROMB_H
#define	ROMB_H
#include <complex>
using namespace std;
//using complex;
extern  const double delta; 

extern  double c1; 
extern  double c2; 
extern  double nu1;
extern  double alpha1;
extern  double nu2;

extern const double alpha2;
extern const double beta1;
extern const double beta2;
extern  int max_stepsX; 
extern  int max_stepsZ; 
extern  int max_stepsT; 
extern complex<double> *A;
extern int choiceT;
extern int choiceZ;
extern int choiceX;
 
/*extern const double error_x1;
extern const double error_x2;
extern const double error_z;
 */
 
 /*extern complex<double> ialpha1_2;
 extern complex<double> ialpha2_2;
 extern double nu1_3;
 extern double nu2_3;*/


//using f_X= double _Comple (*) (complex<double> * , double  ,double , double , double );
 
//typedef double _Comple* (*f_Z) (complex<double> * , double  ,double , double  );

//typedef double _Comple (*f_T) (complex<double> * , double  ,double  );

void dump_row(int , complex<double> *);

complex<double> A_value(complex<double> *,double,double );

complex<double> f_valueX (double, double ,double );

complex<double> romberg(complex<double> * (double , double , double) , double , double,double , double);

complex<double> romberg(complex<double>* (double , double  ), double , double,double);

complex<double> romberg(complex<double>* (double ), double , double,int);

complex<double> rombergA(complex<double>*, double,  double , double);

complex<double> f_valueZ (double , double);

complex<double> f_valueT (double );

complex<double> f_valueTnew (double,int );

complex<double> R_integral (int,int  );

complex<double> R_integral_new (int  );
#endif	/* ROMB_H */

