/* 
 * File:   algorithm.h
 * Author: masha
 *
 * Created on March 3, 2019, 9:38 AM
 */
/* 
 * File:   algorithm.h
 * Author: masha
 *
 * Created on March 3, 2019, 9:38 AM
 */

#ifndef ALGORITHMN_H
#define	ALGORITHMN_H

extern  double c1; 
extern  double c2; 
extern  double nu1_3;
extern  double nu2_3;
extern   complex<double> ialpha1_2;
extern const  complex<double> ialpha2_2;
extern const double beta1;
extern const double beta2;

extern complex<double> *A;
extern const double delta;



extern  double c1; 
extern  double c2; 
extern  double nu1_3;
extern  double nu2_3;

extern const double beta1;
extern const double beta2;

extern double *r;
extern double *w;
extern const double delta;



double  Simpson_extendedN(double (*f)(size_t , size_t, size_t ),size_t ,size_t  , size_t );

double  Simpson_extendedTN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t );


double  BooleN(double (*f)(size_t , size_t, size_t ),size_t ,size_t  , size_t );

double  BooleTN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t );

double  Simpson_narrowN(double (*f)(size_t , size_t, size_t ),
                                                        size_t ,size_t  , size_t );

double  Simpson_narrowTN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t );


double  Simpson_3_8N(double (*f)(size_t , size_t, size_t ),size_t ,size_t  , size_t );

double  Simpson_3_8TN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t );

double  Simpson_ruleN(double (*f)(size_t , size_t, size_t ),
                                                        size_t ,size_t  , size_t );
double  Simpson_ruleTN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t );


double  Trapezoid_ruleN(double (*f)(size_t , size_t, size_t ),
                                                        size_t ,size_t  , size_t );
double  Trapezoid_ruleTN(double (*f)(size_t , size_t ),
                                                        size_t ,size_t,char );


#endif	/* ALGORITHMN_H */

