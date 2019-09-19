/* 
 * File:   algorithm.h
 * Author: masha
 *
 * Created on March 3, 2019, 9:38 AM
 */

#ifndef ALGORITHM_H
#define	ALGORITHM_H

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



complex<double>  Simpson_extended(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  Simpson_extendedT(complex<double> (*f)(int ,int ), int , int  );

complex<double>  Simpson_extendedR(complex<double> (*f)(int ,int ) , int  );

complex<double>  Boole (complex<double> (*f)(int , int, int ),int ,int  , int );


complex<double>  BooleT(complex<double> (*f)(int ,int ), int , int  );

complex<double>  BooleR(complex<double> (*f)(int ,int ),  int  );

complex<double>  Simpson_narrow(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  Simpson_narrowT(complex<double> (*f)(int ,int ), int , int  );

complex<double>  Simpson_narrowR(complex<double> (*f)(int ,int ),  int  );

complex<double>  Simpson_3_8(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  Simpson_3_8T(complex<double> (*f)(int ,int ), int , int  );

complex<double>  Simpson_3_8R(complex<double> (*f)(int ,int ),  int  );

complex<double>  Simpson_rule(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  Simpson_ruleT(complex<double> (*f)(int ,int ),int , int  );


complex<double>  Simpson_ruleR(complex<double> (*f)(int ,int ),int  );

complex<double>  Trapezoid_rule(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  Trapezoid_ruleT(complex<double> (*f)(int ,int ), int , int  );

complex<double>  Trapezoid_ruleR(complex<double> (*f)(int,int  ), int  );

complex<double>  NoInt(complex<double> (*f)(int , int, int ),int ,int  , int );

complex<double>  NoIntT(complex<double> (*f)(int ,int ), int , int  );

complex<double>  NoIntR(complex<double> (*f)(int,int  ), int  );


complex<double>  Trapezoid_ruleR_new(complex<double> (*f)(int  ), int  );

void checkForSaw(int,int );



#endif	/* ALGORITHM_H */

