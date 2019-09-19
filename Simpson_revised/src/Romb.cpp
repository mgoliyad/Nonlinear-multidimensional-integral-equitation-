
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>
#include "Romb.h"
#include <fstream> 


#include <complex>       
using namespace std;

/*using f_X= double _Comple (*) (complex<double> * , double  ,double , double , double );
using f_Z= double _Comple (*) (complex<double> * , double  ,double , double  );
using f_T= double _Comple (*) (complex<double> * , double  ,double   );
*/
//const int max_steps=5; 
//const int number_R =(1 << (max_steps));
const double acc=1E-6;


void
dump_row(int i, complex<double> *R){
  printf("R[%2d] = ", i);
   for(int j = 0; j <= i; ++j){
    printf("%f  ", real(R[j]));
   }
   printf("\n");
}


complex<double> f_valueX (double x ,double z, double t)
{
   // double fa=t-z-z-delta;
        
    if (fabs(t-2.*z)<delta-1E-15){
     //   printf("000, t=%f, z=%f, delta=%f; fa=%f\n",t,z,delta,fa);
        return 0.0;
    }
    else{
        complex<double> arg1,arg2, A_product, sum=0.0;
        complex<double> _Complex_I(0.0,1.0); 
        complex<double> ialpha1_2=_Complex_I*alpha1*alpha1;
        complex<double> ialpha2_2=_Complex_I*alpha2*alpha2;
        double nu1_3=nu1*nu1*nu1;
        double nu2_3=10.*nu1_3;

            
         //   arg1=ialpha1_2*(z*x)-nu1_3*(z*z*x)-beta1*(x);
         //   arg2=ialpha2_2*(z*x)-nu2_3*(z*z*x)-beta2*(x);
        
        arg1=ialpha1_2*z*(z+x)-nu1_3*z*z*(2.*z/3.+x)-beta1*(2.*z+x);
        arg2=ialpha2_2*z*(z+x)-nu2_3*z*z*(2.*z/3.+x)-beta2*(2.*z+x);

        
            int idx1=round ((t-z-x)/delta);
            int idx2=round((t-2.*z-x)/delta);
            
            if ((idx1<0.0)||(idx2<0)){
            printf ("XNEGATIVE, x=%f, z=%f, t=%f\n",x,z,t);
            }
            complex<double> A1=*(A+idx1);
            complex<double> A2=*(A+idx2);
            
            A_product=A1*conj(A2);
            
            sum=A_product*(c1*exp(arg1)+c2*exp(arg2));
          // if (fabs(t-0.9)<1E-5){
         //   printf("R_FX: sum=%f; idx1=%d, idx2=%d, t=%f, Ap=%f\n",real(sum),idx1,idx2,t,real(A_product));
           return sum;
    }
    
}


complex<double>

romberg( complex<double> (*f)(double , double , double  ), double a, double b,double z, double t){
    
    
   
   complex<double> R1[max_stepsX], R2[max_stepsX]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   
   Rp[0] = (f(b,z,t) + f(a,z,t))*h*.5; //first trapezoidal step
  // printf("R_X_Int for a=%f,b=%f,z=%f,t=%f\n",a,b,z,t);
  // dump_row(0, Rp);

   for(int i = 1; i < max_stepsX; ++i){
      h /= 2.;
      complex<double> c = 0;
      int ep = 1 << (i-1); //2^(n-1)
      for(int j = 1; j <= ep; ++j){
             c += f(a+(2*j-1)*h,z,t);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(int j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
   //   dump_row(i, Rc);

      if(i > 1 && norm(Rp[i-1]-Rc[i]) < norm(Rp[i-1])*acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }

   return Rp[max_stepsX-1]; //return our best guess
}

complex<double> romberg( complex<double> (*f)(double , double  ), double a, double b, double t ){
    
   complex<double> R1[max_stepsZ], R2[max_stepsZ]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   
    Rp[0] = (f(b,t) + f(a,t))*h*.5; //first trapezoidal step

         // dump_row(0, Rp);

   for(int i = 1; i < max_stepsZ; ++i){
      h /= 2.;
      complex<double> c = 0;
      int ep = 1 << (i-1); //2^(n-1)
      for(int j = 1; j <= ep; ++j){
            c += f(a+(2*j-1)*h,t);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(int j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
   //   dump_row(i, Rc);

      if(i > 1 && norm(Rp[i-1]-Rc[i]) < norm(Rp[i-1])*acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
   return Rp[max_stepsZ-1]; //return our best guess

}

complex<double> romberg( complex<double> (*f)(double,int ), double a, double b,int idx_t){
    
   complex<double> R1[max_stepsT], R2[max_stepsT]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size

       Rp[0] = (f(b,idx_t) + f(a,idx_t))*h*.5;
   
 //  dump_row(0, Rp);

   for(int i = 1; i < max_stepsT; ++i){
      h /= 2.;
      complex<double> c = 0;
      int ep = 1 << (i-1); //2^(n-1)
      for(int j = 1; j <= ep; ++j){
          c += f(a+(2*j-1)*h,idx_t);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(int j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
   
        //      dump_row(i, Rc);

      if(i > 1 && norm(Rp[i-1]-Rc[i]) < norm(Rp[i-1])*acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
              
   return Rp[max_stepsT-1]; //return our best guess
  

}


complex<double> f_valueZ (double z, double t)
{
    
    if (z<delta) {
        return 0.0;
    }
    else
    {
        complex<double> _Complex_I(0.0,1.0); 
        
        complex<double> ialpha1_2=_Complex_I*alpha1*alpha1;
        complex<double> ialpha2_2=_Complex_I*alpha2*alpha2;
        double nu1_3=nu1*nu1*nu1;
        double nu2_3=10.*nu1_3;
        
        int NUM= (fabs(t-2.*z)<delta) ? 0 : floor((t-2.*z)/delta);
        
        complex<double> sum=0.0;
        // printf("f_valueZ for z=%f, t=%f; NUM=%d\n",z,t,NUM);
        fflush(stdout);
        for(int i = 1; i <= NUM; ++i){
            sum+= romberg(f_valueX,(i-1)*delta, i*delta, z,t);
        
        }
    
        double interval;
        interval=(t-2.*z-NUM*delta);
 
        if (abs(interval)>1E-20){
            
            fflush(stdout);
                sum+=abs(interval)*f_valueX(NUM*delta,z,t);
                complex<double> help=f_valueX(NUM*delta,z,t);
             //   printf("Int2: NUM=%d, z=%f, t=%f, sum=%f, interval=%f, fx=%f\n",NUM,z,t,real(sum), interval, real(help));
                
        }
       
           int idx1=round ((t-z)/delta);

             sum*=A[idx1]*z*z;
        //   printf("f_valueZ:  t=%f; z=%f; sum=%f, idx=%d, A=%f\n",t,z,real(sum), idx1, real(A[idx1]));
           return sum;
    }
    
}

complex<double> f_valueT (double t,int idx_tM)
{ 
    int idx_t=floor(t/delta);
    int NUM=(idx_t)>>1;
    complex<double> sum=0.0;
    
   // printf("Start f_valueT for idx_t=%d, t=%f, NUm=%d\n",idx_t,t, NUM);
    fflush(stdout);

    int i = 1;
        while (i<=NUM){    
            sum+=romberg(f_valueZ,(i-1)*delta, i*delta, t);
          //  printf("Ri=%f, i=%d\n",real(sum),i);
            ++i;
        }  
    
        if (abs(t-2.*NUM*delta)>1E-20){
                sum+=0.5*abs(t-2.*NUM*delta)*f_valueZ(NUM*delta,t);
             //   printf("extra Ri=%f\n",real(sum));
        }
        //  printf("11finish f_valueT for idx_t=%d, t=%f,  NUm=%d, sum=%f\n",idx_t,t,NUM,real(sum));
          fflush(stdout);
        return sum;

}

complex<double> f_valueTnew (double t,int idx_tM)
{ 
    int idx_t=floor(t/delta);
    int NUM=(idx_t)>>1;
    complex<double> sum=0.0;
    
   // printf("Start f_valueT for idx_t=%d, t=%f, NUm=%d\n",idx_t,t, NUM);
    fflush(stdout);
  
    int i = 1;
        while (i<=NUM){    
            sum+=romberg(f_valueZ,(i-1)*delta, i*delta, t);
          //  printf("Ri=%f, i=%d\n",real(sum),i);
            ++i;
        }  
         
        if (abs(t-2.*NUM*delta)>1E-20){
                sum+=0.5*abs(t-2.*NUM*delta)*f_valueZ(NUM*delta,t);
             //   printf("extra Ri=%f\n",real(sum));
        }
        //  printf("finish f_valueT for idx_t=%d, t=%f, NUm=%d\n",idx_t,t, NUM);
          fflush(stdout);
          double t_exp=idx_tM*delta-t;
         // printf("finish f_valueT for idx_t=%d, t=%f, t_exp=%f, NUm=%d, sum=%f\n",idx_t,t, t_exp,NUM,real(sum));
          fflush(stdout);
          sum/=exp(t);
         // sum/=exp(t_exp);
        //  printf(" f_valueT for idx_t=%d, t=%f, t_exp=%f, NUm=%d, sum=%f\n",idx_t,t, t_exp,NUM,real(sum));
          fflush(stdout);
        return sum;
       //  return A[idx_t]-0.5*sum;
}
/*********/

/*complex<double> A_value(complex<double> *A, double v,double delta)
{
    int idx=floor(v/delta);
    
    complex<double> A1=A[idx];
    complex<double> A2=A[idx+1];   
    complex<double> A_new=(A2-A1)*(v/delta-idx-1.)+A2;
    
     printf("Romb: min=%f, max=%f, A_mid=%f\n", real(A1),real(A2),real(A_new));
     
      return      A_new;
            
}

complex<double> rombergA( complex<double> *A, double a, double b, double delta){
    
   complex<double> R1[max_stepsT], R2[max_stepsT]; //buffers
   complex<double> *Rp = &R1[0], *Rc = &R2[0]; //Rp is previous row, Rc is current row
   double h = (b-a); //step size
   
    int idx_a=floor(a/delta);
    int idx_b=floor(b/delta);
  
       Rp[0] = .5*(A[idx_a] + A[idx_b])*h;
   
  // dump_row(0, Rp);

   for(int i = 1; i < max_stepsT; ++i){
      h /= 2.;
      
      complex<double> c = 0;
      int ep = 1 << (i-1); //2^(n-1)
      for(int j = 1; j <= ep; ++j){
         // idx_a=floor((a+(2*j-1)*h)/delta);
        //  c += A[idx_a];
          c +=A_value(A,a+(2*j-1)*h,delta);
      }
      Rc[0] = h*c + .5*Rp[0]; //R(i,0)

      for(int j = 1; j <= i; ++j){
         double n_k = pow(4, j);
         Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
      }

      //Dump ith column of R, R[i,i] is the best estimate so far
   
            //  dump_row(i, Rc);

      if(i > 1 && norm(Rp[i-1]-Rc[i]) < norm(Rp[i-1])*acc){
         return Rc[i-1];
      }

      //swap Rn and Rc as we only need the last row
      complex<double> *rt = Rp;
      Rp = Rc;
      Rc = rt;
   }
              
   return Rp[max_stepsT-1]; //return our best guess
  

}*/
/*********/
complex<double> R_integral (int choice, int idx_t)
{
  //   printf("startRomb, choice=%d, t=%d, delta=%f\n",choice, idx_t, delta);
     fflush(stdout);
     if (choice!=0){
       return romberg(f_valueT,(idx_t-choice)*delta, (idx_t)*delta, idx_t);
     }
     else{
        return f_valueT((idx_t)*delta, idx_t);
     }
//return f_valueT(idx_t*delta);
}

complex<double> R_integral_new (int idx_t)
{
    // printf("startRomb,  t=%d, delta=%f\n", idx_t, delta);
     fflush(stdout);
   
      // return romberg(f_valueTnew,(idx_t-2)*delta, (idx_t)*delta,idx_t);
     double delta_3=delta/3.;

        return delta_3*(  exp(delta-2)*f_valueTnew((idx_t-2)*delta,idx_t) +  4.*exp(delta)*(  f_valueTnew((idx_t-1)*delta,idx_t)  )+
                                    f_valueTnew((idx_t)*delta,idx_t) );
//return f_valueT(idx_t*delta);
}