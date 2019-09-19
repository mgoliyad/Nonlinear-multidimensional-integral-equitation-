/* 
 * File:   Comp.h
 * Author: masha
 *
 * Created on February 16, 2019, 8:48 AM
 */

#ifndef ROMBERG_H
#define	ROMBERG_H
/*extern  double delta;
extern const double c1; 
extern const double c2; 
extern const double nu1;
extern const double alpha1;
extern const double nu2;
*/
extern  double c1; 
extern  double c2; 
extern  double nu1;
extern  double alpha1;
extern  double nu2;

extern const double alpha2;
extern const double beta1;
extern const double beta2;
extern const double error_x1;
extern const double error_x2;
extern const double error_z;
 
 
 /*extern double  _Complex ialpha1_2;
 extern double  _Complex ialpha2_2;
 extern double nu1_3;
 extern double nu2_3;*/


//using f_X= double _Comple (*) (double  _Complex * , double  ,double , double , double );
 
//typedef double _Comple* (*f_Z) (double  _Complex * , double  ,double , double  );

//typedef double _Comple (*f_T) (double  _Complex * , double  ,double  );

void dump_row(size_t , double _Complex *);

double _Complex f_valueX (double  _Complex *, double, double ,double , double );

double _Complex romberg(double  _Complex * (double  _Complex * , double  ,double , double , double ) , 
        double , double,double , double,double _Complex*,double );
double _Complex romberg(double  _Complex* (double  _Complex * , double  ,double , double  ), double , double,double ,
        double _Complex* ,double );
double _Complex romberg(double  _Complex* (double  _Complex * , double  ,double  ), double , double,
        double _Complex* ,double );

double _Complex f_valueZ (double  _Complex *,double,  double , double );


double _Complex f_valueT (double  _Complex *, double, double );


double _Complex R_integral (double  _Complex *, size_t , double );


#endif	/* COMP_H */

