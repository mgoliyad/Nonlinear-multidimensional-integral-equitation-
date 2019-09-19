/* 
 * File:   Comp.h
 * Author: masha
 *
 * Created on February 16, 2019, 8:48 AM
 */

#ifndef COMP_H
#define	COMP_H
#include <complex>
using namespace std;
extern const double delta;
extern  double c1; 
extern  double c2; 

extern double nu1_3;
extern double nu2_3;
extern const complex<double> _Complex_I;
extern complex<double> ialpha1_2;
extern const  complex<double> ialpha2_2;



extern const double beta1;
extern const double beta2;
extern const double error_x1;
extern const double error_x2;
extern const double error_z;

//using namespace std;
//using std::complex;

//extern const   complex<double> _Complex_I; 

std::complex<double>  Mround(std::complex<double>);

std::complex<double> FuncForX (std::complex<double>*, int ,int , int ,int,double );

std::complex<double> X_integral(std::complex<double>*, int , int ,int,double );

std::complex<double> SimpsonX(std::complex<double>*, int , int ,int,double );

std::complex<double> FuncForZ(std::complex<double>*, int , int,double);

std::complex<double> callSimpsonX(std::complex<double>*, int , int ,int,double );

std::complex<double> FuncForZSimp(std::complex<double>*, int , int,double);

std::complex<double> Z_integral(std::complex<double>*, int,double);

std::complex<double> SimpsonZ(std::complex<double>*, int,double );

std::complex<double> callSimpsonZ(std::complex<double>*, int,double);

void simple_integration(std::complex<double>*, int,double );

void Simpson_integrationT(std::complex<double>*, int,double );


           /* #ifdef	__cplusplus
            extern "C" {
            #endif




            #ifdef	__cplusplus
            }
            #endif*/

#endif	/* COMP_H */

