#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>
#include "Romberg.h"

#include <complex>       
using namespace std;

/*using f_X= double _Comple (*) (double  _Complex * , double  ,double , double , double );
using f_Z= double _Comple (*) (double  _Complex * , double  ,double , double  );
using f_T= double _Comple (*) (double  _Complex * , double  ,double   );
*/
const size_t max_steps=3; 
const size_t number_R =(1 << (max_steps));
const double acc=1E-4;


void
dump_row(size_t i, complex<double> *R){
   printf("R[%2zu] = ", i);
   for(size_t j = 0; j <= i; ++j){
      printf("%f; %f ", cabs(R[j]),carg(R[j]));
   }
   printf("\n");
}

complex<double> f_valueX (double  _Complex *A, double x ,double z, double t, double delta)
{
    // printf("Enter f_valueX, x=%f, z=%f, t=%f, delta=%f;nR=%2zu\n",x,z,t,delta,number_R);
     
    double  _Complex A_product, sum;

        complex<double> arg1;
        complex<double> arg2;
        double  arg3;
        double  arg4;
        
        double  _Complex ialpha1_2=_Complex_I*alpha1*alpha1;
 double  _Complex ialpha2_2=_Complex_I*alpha2*alpha2;
 double nu1_3=nu1*nu1*nu1;
 double nu2_3=10.*nu1_3;

            arg1=ialpha1_2*(z*x+z*z);
            arg2=ialpha2_2*(z*x+z*z);
            
            arg3=-nu1_3*(2.*z*z*z/3.+z*z*x)-beta1*(2.*z+x);
            arg4=-nu2_3*(2.*z*z*z/3.+z*z*x)-beta2*(2.*z+x);
            
            int idx1=(t-z-x)/delta;
            int idx2=(t-2*z-x)/delta;
            complex<double> A1=A[idx1];
            complex<double> A2=conj(A[idx2]);
            A_product=A1*A2;
            
            sum=z*z*A_product*(c1*cexp(arg1)*exp(arg3)+c2*cexp(arg2)*exp(arg4));
    //  printf("Exit f_valueX, idx1=%i, idx2=%i,A_p=%f, A1=%f,A2=%f, sum=%f\n",idx1, idx2, cabs(A_product), cabs(A1),
//              cabs(A2),cabs(sum));
     
return sum;}


complex<double>

romberg( complex<double> *f(double  _Complex * , double  ,double , double , double ), double a, double b,double z, double t, complex<double> *A, double delta){
    
   complex<double> R1[max_steps], R2[max_steps]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   Rp[0] = (f(A,a,z,t,delta) + f(A,b,z,t,delta))*h*.5; //first trapezoidal step

   dump_row(0, Rp);

   for(size_t i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(size_t j = 1; j <= ep; ++j){
         c += f(A,a+(2*j-1)*h,z,t,delta);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(size_t j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
      dump_row(i, Rc);

      if(i > 1 && cabs(Rp[i-1]-Rc[i]) < acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_steps-1]; //return our best guess
}

complex<double> romberg( complex<double> *f(double  _Complex * , double  ,double , double ), double a, double b, double t, complex<double> *A, double delta ){
    
   complex<double> R1[max_steps], R2[max_steps]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   Rp[0] = (f(A,a,t,delta) + f(A,b,t,delta))*h*.5; //first trapezoidal step

   dump_row(0, Rp);

   for(size_t i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(size_t j = 1; j <= ep; ++j){
         c += f(A,a+(2*j-1)*h,t,delta);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(size_t j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
      dump_row(i, Rc);

      if(i > 1 && cabs(Rp[i-1]-Rc[i]) < acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_steps-1]; //return our best guess
}

complex<double> romberg( complex<double> *f(double  _Complex * , double  ,double), double a, double b, complex<double> *A , double delta){
    
   complex<double> R1[max_steps], R2[max_steps]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   Rp[0] = (f(A,a,delta) + f(A,b,delta))*h*.5; //first trapezoidal step

   dump_row(0, Rp);

   for(size_t i = 1; i < max_steps; ++i){
      h /= 2.;
      double c = 0;
      size_t ep = 1 << (i-1); //2^(n-1)
      for(size_t j = 1; j <= ep; ++j){
         c += f(A,a+(2*j-1)*h,delta);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(size_t j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
      dump_row(i, Rc);

      if(i > 1 && cabs(Rp[i-1]-Rc[i]) < acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_steps-1]; //return our best guess
}


complex<double> f_valueZ (double  _Complex *A, double t, double z,double delta)
{
   // printf("Enter f_valueZ,  z=%f, t=%f, delta=%f;nR=%2zu\n",z,t,delta,number_R);
    double h=delta/number_R;
    size_t NUM=(t-2.*z)/delta;
    complex<double> sum=0.0;
    
  //  complex<double> *f=( double  _Complex*)malloc((number_R+1)*sizeof(double  _Complex));
    if ((NUM>0)&&(z>0)){
    for(size_t i = 1; i <= NUM; ++i){
        
      /*  for(size_t j = 0; j <= number_R; ++j){
            
            f[j]=f_valueX(A,(i-1)*delta+j*h,z,t,delta);
        }*/
    
      sum+= romberg(f_valueX,(i-1)*delta, i*delta, z,t,A,delta);
    }  
    }
    else{
        return 0.0;
    }
    return sum;
}

complex<double> f_valueT (double  _Complex *A, double t, double delta)
{
   // printf("Enter f_valueT, t=%f, delta=%f, nR=%2zu\n",t,delta, number_R);
    double h=delta/number_R;
    size_t NUM=t/(2.*delta);
    complex<double> sum=0.0;
    
 //  complex<double> *f=( double  _Complex*)malloc((number_R+1)*sizeof(double  _Complex));
    if (NUM>0){
    for(size_t i = 1; i <= NUM; ++i){
        
       /* for(size_t j = 0; j <= number_R; ++j){
            f[j]=f_valueZ(A,t,(i-1)*delta+j*h,delta);
        }*/
    
      sum+= romberg(f_valueZ,(i-1)*delta, i*delta, t,A,delta);
    }  
    }
    else{
        return 0.0;
    }
    return sum;
}

complex<double> R_integral (double  _Complex *A, size_t idx_t, double delta)
{
    double h=delta/number_R;
    complex<double> sum=0.0;
    double tt=idx_t*delta;
    double tt_1=(idx_t-1)*delta;
    
  /* complex<double> *f=( double  _Complex*)malloc((number_R+1)*sizeof(double  _Complex));
    
    for(size_t i = 0; i <= number_R; ++i){
        f[i]=f_valueT(A,(idx_t-1)*delta+i*h,delta);
    }*/
    
   //return (A[t-delta]*(1+delta/2.0)-0.25*delta*(f_valueT(A,t)+f_valueT(A,t-delta)))/(1-delta/2.0);
   // return A[idx_t-1]+romberg(A,tt_1,tt)-0.5*romberg(f,tt_1,tt);
   return (A[idx_t-1]*(1+delta)-0.5*romberg(f_valueT,tt_1, tt,A,delta))/(1-delta/2.0);
}