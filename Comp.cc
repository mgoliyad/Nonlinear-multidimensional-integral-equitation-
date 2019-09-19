#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>

#include "Comp.h"

complex<double> FuncForX (complex<double> *A, int idx_x,int idx_z, int idx_t,int trig,double delta)
//Here we calculate value of first function(x) to be integrated by x for x=idx_x*delta.
{
    complex<double> A_product, sum,A_e;
    double z_x, z2_x;
    double delta2=delta*delta;
    double delta3=delta2*delta;

        //printf("c1=%f,c2=%f\n",c1,c2);
        z_x=static_cast<double>(idx_z*idx_x)*delta2;
        //  z_x=delta2*(double)idx_z*(double)idx_x;
        z2_x=delta3*(double)idx_z*(double)idx_x*(double)idx_z;

        complex<double> arg1;
        complex<double> arg2;

        if (trig==1){
            arg1=ialpha1_2*z_x;
            arg2=-nu1_3*z2_x-beta1*idx_x*delta;
            A_product=(A[idx_t-idx_z-idx_x]*conj(A[idx_t-2*idx_z-idx_x]));
            sum=A_product*exp(arg1)*exp(arg2);
        }
        else{
            arg1=ialpha2_2*z_x;
            arg2=-nu2_3*z2_x-beta2*idx_x*delta;
            A_product=(A[idx_t-idx_z-idx_x]*conj(A[idx_t-2*idx_z-idx_x]));
            sum=A_product*exp(arg1)*exp(arg2);
        }

return sum;
	}

/**************************************************************************************************/

complex<double> X_integral (complex<double> *A,int idx_z, int idx_t,int trig,double delta)
{
    complex<double> sum_x=0.0;
    int idx_x;
    if((idx_t-2*idx_z)>0)
    {
         sum_x=0.5*delta*(FuncForX (A, 0,idx_z,  idx_t,trig,delta) + 
                 FuncForX (A, (idx_t-2*idx_z),idx_z,  idx_t,trig,delta));

		for(idx_x=1;idx_x<(idx_t-2*idx_z);idx_x++)
		{
                    sum_x+=FuncForX (A, idx_x,idx_z,  idx_t,trig,delta)*delta;
		}
    }

return sum_x;
	}

/**************************************************************************************************/

complex<double>  SimpsonX(complex<double> *A, int idx_z, int idx_t, int trig,double delta)

{
    int idx_x,x_Max;
    complex<double> sum_z=0.0,sum_last=0.0;
    double delta_3=delta/3.0;

    if ((idx_t-2*idx_z)%2==0)
    {
        x_Max=idx_t-2*idx_z;
    }
    else{
        x_Max=idx_t-2*idx_z-1;
        
        sum_last=0.5*delta*(FuncForX(A,x_Max+1,idx_z, idx_t, trig,delta)+
                FuncForX(A,x_Max,idx_z, idx_t, trig,delta));
    }

    sum_z= FuncForX(A,0,idx_z, idx_t, trig,delta)+
            FuncForX(A,x_Max,idx_z, idx_t, trig,delta);


    for (idx_x=1;idx_x<x_Max;idx_x++)
    {
        if ((idx_x)%2==1){
            sum_z+= 4.0*FuncForX(A,idx_x,idx_z, idx_t,trig,delta);
        }
        else{
            sum_z+= 2.0*FuncForX(A,idx_x,idx_z, idx_t,trig,delta);
        }
       
    }

    sum_z=sum_z*delta_3+sum_last;


return sum_z;
}

/**************************************************************************************************/
complex<double> FuncForZ (complex<double> *A, int idx_z, int idx_t,double delta)
{
    complex<double>  X_int1, X_int2,sum;
    double  z2, z3;

    double delta2=delta*delta;
    double delta3=delta2*delta;

    z3=delta3*(double)idx_z*(double)idx_z*(double)idx_z;
    z2=delta2*(double)idx_z*(double)idx_z;


        X_int1=X_integral(A,idx_z, idx_t,1,delta);
	X_int2=X_integral(A,idx_z, idx_t,2,delta);

        complex<double> arg1=ialpha1_2*z2-2.*nu1_3*z3/3.-beta1*2.0*idx_z*delta;

        sum=c1*z2*A[idx_t-idx_z]*X_int1*exp(arg1);
        
        arg1=ialpha2_2*z2-2.*nu2_3*z3/3.-beta2*2.0*idx_z*delta;
        
        sum+=c2*z2*A[idx_t-idx_z]*X_int2*exp(arg1);
      
     //   printf("FuncZ %f for z=%zu\n",real(sum),idx_z);
return sum;
}

/**************************************************************************************************/

complex<double> callSimpsonX(complex<double> *A, int idx_z1, int idx_t,int trig,double delta)
{
    if (idx_t-2*idx_z1<2){
        return X_integral(A,idx_z1, idx_t,trig,delta);
    }
    else{
        return SimpsonX(A,idx_z1, idx_t,trig,delta);
    }
}
/**************************************************************************************************/
complex<double> FuncForZSimp (complex<double> *A, int idx_z, int idx_t,double delta)
{
    complex<double>  X_int1, X_int2,sum=0.0;
    double delta2=delta*delta;
    double delta3=delta2*delta;

        double z3=delta3*(double)idx_z*(double)idx_z*(double)idx_z;
        double z2=delta2*(double)idx_z*(double)idx_z;

        X_int1=callSimpsonX(A,idx_z, idx_t,1,delta);
	X_int2=callSimpsonX(A,idx_z, idx_t,2,delta);

        complex<double> arg1=ialpha1_2*z2-2.*nu1_3*z3/3.-beta1*2.0*idx_z*delta;

        sum=z2*A[idx_t-idx_z]*X_int1*c1*exp(arg1);
     
        arg1=ialpha2_2*z2-2.*nu2_3*z3/3.-beta2*2.0*idx_z*delta;
        
        sum+=z2*A[idx_t-idx_z]*X_int2*c2*exp(arg1);

return sum;
}

/**************************************************************************************************/

complex<double> Z_integral (complex<double> *A, int idx_t,double delta)
{
    complex<double> sum_z=0.0;
	int idx_z, idx_t2;

    idx_t2=idx_t>>1;
     if(idx_t2>0){
            for(idx_z=1;idx_z<=idx_t2;idx_z++){
            sum_z+=0.5*delta*( FuncForZ(A,idx_z,idx_t,delta) +FuncForZ(A,idx_z-1,idx_t,delta));
            }
 
            if ((idx_t%2)==1){
                sum_z+=delta*(FuncForZ(A,idx_t2, idx_t,delta))*0.5;
            }
  
        }
      
return sum_z;
}

/**************************************************************************************************/

complex<double>  SimpsonZ(complex<double> *A, int idx_t,double delta)
{
int idx_z, z_Max=0;
    complex<double> sum_z=0.0,sum_last=0.0;
    double delta_3=delta/3.0;

    int idx_t2=idx_t>>1;

    if ((idx_t2)%2==0){
        z_Max=idx_t2;
    }
    else{
        z_Max=idx_t2-1;
        sum_last=0.5*delta*( FuncForZSimp(A,idx_t2, idx_t,delta)+ FuncForZSimp(A,idx_t2-1, idx_t,delta));

    }

    sum_z= (FuncForZSimp(A,z_Max, idx_t,delta));

    for (idx_z=1;idx_z<z_Max;idx_z++)
    {
        if ((idx_z)%2==1){
            sum_z+= 4.*FuncForZSimp(A,idx_z, idx_t,delta);
        }
        else{
            sum_z+= 2.*FuncForZSimp(A,idx_z, idx_t,delta);
        }

    }

    sum_z=sum_z*delta_3 + sum_last;

    if (idx_t%2==1){
            sum_z+=0.5*delta*FuncForZSimp(A,idx_t2, idx_t,delta);
        }

return sum_z;
}

/**************************************************************************************************/

complex<double> callSimpsonZ (complex<double> *A, int idx_t,double delta)
{
   if (idx_t<4){
        return Z_integral(A,idx_t,delta);
    }
    else {
        return SimpsonZ(A,idx_t,delta);
    }
}
/**************************************************************************************************/
void  simple_integration (complex<double> *A, int idx,double delta)
{
	A[idx]=(  A[idx-1]*(1.+0.5*delta) - 0.25*delta*
                ( callSimpsonZ(A,idx,delta) + callSimpsonZ(A,idx-1,delta)) )
                                    /(1.0-0.5*delta);
}
/**************************************************************************************************/
void  Simpson_integrationT(complex<double> *A, int idx,double delta)

{
complex<double> sum;
double delta_3=delta/3.;
double delta_3_8=3.*delta/8.;
double delta_24=delta/24.;
double delta_48=delta/48.;

     if (idx<3){

        sum=-(0.5*delta_3)*(  callSimpsonZ(A,idx,delta) +  4.*(  callSimpsonZ(A,idx-1,delta)  )+
                                    callSimpsonZ(A,idx-2,delta) );
        
      //  printf("sum=%f\n",real(sum));
        sum+= A[idx-2]*(1.0+delta_3);
       //  printf("sum=%f\n",real(sum));
        sum+=(4.* delta_3)*(A[idx-1]);
       //  printf("sum=%f\n",real(sum));
        A[idx]=sum/(1.-delta_3);
       //  printf("sum=%f\n",real(A[idx]));
        

    }
    if ((idx>=3)&& (idx<6))
  
        {
        sum=-(0.5*delta_3_8)*(
                              callSimpsonZ(A,idx,delta) +
                              3.*(  callSimpsonZ(A,idx-1,delta) + callSimpsonZ(A,idx-2,delta) )+
                                    callSimpsonZ(A,idx-3,delta)
                               );
        //printf("Int for simps sum=%f\n",real(sum));
         sum+=A[idx-3]*(1.+delta_3_8);
        // printf("sum=%f\n",real(sum));
         sum+=(3.*delta_3_8)*(A[idx-2]+A[idx-1]);
       //  printf("sum=%f\n",real(sum));
        A[idx]=sum/(1.-delta_3_8);
      //  printf("sum=%f\n",real(A[idx]));

    }

     if ((idx>=6)&& (idx<8))
        {
        
        sum=-(0.5*delta_24)*(  9.*(callSimpsonZ(A,idx, delta) +callSimpsonZ(A,idx-6, delta))+
                               28.*(  callSimpsonZ(A,idx-1, delta) + callSimpsonZ(A,idx-5, delta) )+
                                  23.*(  callSimpsonZ(A,idx-2, delta) + callSimpsonZ(A,idx-4, delta) )+
                             24.*callSimpsonZ(A,idx-3, delta));

      //   printf("Sim: idx=%zu: sum=%f\n",idx, real(sum));
        
        sum+=A[idx-6]*(1.+9.*delta_24);
         sum+=(delta_24)*(28.*A[idx-5]+23.*A[idx-4]+24.*A[idx-3]+23.*A[idx-2]+28.*A[idx-1]);
         
        A[idx]=sum/(1.-9.*delta_24);
       
    }
if (idx>=8)
        {
      //  printf("Smp\n");
        sum=-(0.5*delta_48)*(  17.*(callSimpsonZ(A,idx, delta) +callSimpsonZ(A,idx-8, delta))+
                               59.*(  callSimpsonZ(A,idx-1, delta) + callSimpsonZ(A,idx-7, delta) )+
                               43.*(  callSimpsonZ(A,idx-2, delta) + callSimpsonZ(A,idx-6, delta) )+
                               49.*(  callSimpsonZ(A,idx-3, delta) + callSimpsonZ(A,idx-5, delta) )+
                             48.*callSimpsonZ(A,idx-4, delta));
        
      //   printf("isum=%f\n",real(sum));
        sum+=A[idx-8]*(1.+17.*delta_48);
         sum+=(delta_48)*(59.*A[idx-7]+43.*A[idx-6]+49.*A[idx-5]+48.*A[idx-4]+49.*A[idx-3]+43.*A[idx-2]+59.*A[idx-1]);
       //  printf("isum=%f\n",real(sum));
        A[idx]=sum/(1.-17.*delta_48);
      // printf("isum=%f\n",real(A[idx]));
    }
}
