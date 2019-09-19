/* 
 * File:   NonComp.h
 * Author: masha
 *
 * Created on February 16, 2019, 9:04 AM
 */

#ifndef NONCOMP_H
#define	NONCOMP_H

/*extern  double delta;
extern const double c1; 
extern const double c2; 
extern const double nu1;
extern const double alpha1;
extern const double nu2;*/

extern  double c1; 
extern  double c2; 
extern  double nu1;
extern  double alpha1;
extern  double nu2;


extern const double alpha2;
extern const double beta1;
extern const double beta2;
extern const double error_r;
extern const double error_w;


double phi(size_t , size_t ,double ,double  );

double psi(size_t , size_t ,double ,double  );

double  Omega (double  *, size_t , size_t ,size_t,double ,double );

double R(double  *, size_t , size_t , size_t ,double  );

double FuncForX_N (double *, double * ,size_t ,size_t, size_t ,size_t ,double );

double  X_size_tegral_N(double *, double * ,size_t ,size_t, size_t,double  );

double  SimpsonX_N(double *, double * ,size_t ,size_t, size_t,double  );

double callSimpsonX_N(double *, double * ,size_t ,size_t,size_t,double  );

double FuncForZ_N(double *, double * ,size_t ,size_t, size_t,double  );

double Z_integral_N (double *, double * ,size_t ,size_t,double  );

double   SimpsonZ_N(double *, double * ,size_t ,size_t,double  );

double callSimpsonZ_N  (double *, double * ,size_t ,size_t,double );

double FuncForUN(double *, double * ,size_t ,double );

double FuncForWN(double *, double * ,size_t,double  );

double  Simpson_integrationUN(double *, double * ,size_t,double  );

double  Simpson_integrationWN(double *, double * ,size_t,double  );

void Iteration(double *, double * ,size_t ,double );

void simple_diff_eq(double *, double * ,size_t,double );

        /*    #ifdef	__cplusplus
            extern "C" {
            #endif




            #ifdef	__cplusplus
            }
            #endif*/

#endif	/* NONCOMP_H */

