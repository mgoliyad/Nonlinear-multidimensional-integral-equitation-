
#ifndef NEWFILE_H
#define	NEWFILE_H

#include <complex>
using namespace std;

extern  double c1; 
extern  double c2; 
extern  double nu1_3;
extern  double nu2_3;
extern   complex<double> ialpha1_2;
extern  const complex<double> ialpha2_2;
extern const double beta1;
extern const double beta2;
extern complex<double> *A;
extern const double delta;
extern int choiceT;
extern int choiceZ;
extern int choiceX;
 


complex<double> FuncForX (int , int ,int , int,int );



complex<double> FuncForZ (int , int , int );



//complex<double> FuncForT (  complex<double> *, int , int );

complex<double> FuncForT (  int , int );

complex<double> callMethodX(complex<double> (*f)(int ,int , int ), int , int ,int , int );

complex<double> callMethodZ(complex<double> (*f)(int ,int , int ), int , int ,int  );

complex<double> callMethodT(complex<double> (*f)(  int ,int  ),  int , int  );
complex<double> callMethodR(complex<double> (*f)(  int ,int  ),  int , int  );
void choseMethod(int , int,complex<double>* );
#endif	/* NEWFILE_H */

