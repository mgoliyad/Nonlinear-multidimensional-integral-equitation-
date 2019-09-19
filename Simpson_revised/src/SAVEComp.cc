#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
//#include <complex.h>
#include <complex>       
using namespace std;
using std::complex;


#include "Comp.h"
#include "NonComp.h"
#include "Romb.h"



complex<double> Mround(complex<double> var)
{
    // we use array of chars to store number
    // as a string.
    char str1[40];
    char str2[40];
    
    complex<double> _Complex_I(0.0,1.0); 

   // double A1=std::real(var);
     double A1=real(var);
    double A2=imag(var);

    // Prsize_t in string the value of var
    // with two decimal posize_t
    sprintf(str1, "%.2f", A1);
    sprintf(str2, "%.2f", A2);
    // scan string value in var
    sscanf(str1, "%lf", &A1);
    sscanf(str2, "%lf", &A2);
    var=A1+_Complex_I*A2;

    return var;
}
/**************************************************************************************************/

complex<double> FuncForX (complex<double> *A, size_t idx_x,size_t idx_z, size_t idx_t,size_t trig,double delta)
//Here we calculate value of first function(x) to be size_tegrated by x for x=idx_x*delta.
{
    complex<double> _Complex_I(0.0,1.0); 
    complex<double> A_product, sum,A_e;
    double z_x, z2_x;
    double delta2=delta*delta;
    double delta3=delta2*delta;
    double nu1_3=nu1*nu1*nu1;
    double nu2_3=nu2*nu2*nu2;
    complex<double> ialpha1_2=_Complex_I*alpha1*alpha1;
    complex<double> ialpha2_2=_Complex_I*alpha2*alpha2;
    //prsize_tf("c1=%f,c2=%f\n",c1,c2);
        z_x=delta2*(double)idx_z*(double)idx_x;
        z2_x=delta3*(double)idx_z*(double)idx_x*(double)idx_z;

        complex<double> arg;
        
             
        if (trig==1){
            arg=ialpha1_2*z_x-nu1_3*z2_x-beta1*idx_x*delta;
            
            A_product=(A[idx_t-idx_z-idx_x]*conj(A[idx_t-2*idx_z-idx_x]));           
            sum=A_product*exp(arg);
            
        }
        else{
            arg=ialpha2_2*z_x-nu2_3*z2_x-beta2*idx_x*delta;
            A_product=(A[idx_t-idx_z-idx_x]*conj(A[idx_t-2*idx_z-idx_x]));
            sum=A_product*exp(arg);
        }
             
 // if(cabs(sum)<error_x1){
      //  printf("X1:sum=%f, x=%f, z=%f, t=%f; idx1=%zu; idx2=%zu \n",real(sum),idx_x*delta, idx_z*delta, idx_t*delta,idx_t-idx_z-idx_x,idx_t-2*idx_z-idx_x);
  //  sum=0.0;}
             
return sum;
	}

/**************************************************************************************************/

complex<double> X_integral (complex<double> *A,size_t idx_z, size_t idx_t,size_t trig,double delta)
{
    complex<double> sum_x=0.0;
    complex<double> help;
    size_t idx_x;
    if((idx_t-2*idx_z)>0)
    {
        sum_x=0.5*delta*(FuncForX (A, 0,idx_z,  idx_t,trig,delta) + FuncForX (A, (idx_t-2*idx_z),idx_z,  idx_t,trig,delta));

		for(idx_x=1;idx_x<(idx_t-2*idx_z);idx_x++)
		{             
                    sum_x+=FuncForX (A, idx_x,idx_z,  idx_t,trig,delta)*delta;
		}
         
    }

       // if(cabs(sum_x)<error_x1){
            // prsize_tf("X1_I:sum=%f, idx_x=%i, idx_z=%i, idx_t=%i; \n",cabs(sum_x),idx_x, idx_z, idx_t);
         //   sum_x=0.0;}
return sum_x;
	}

/**************************************************************************************************/

complex<double>  SimpsonX(complex<double> *A, size_t idx_z, size_t idx_t, size_t trig,double delta)

//Function to calculate first size_tergral by x  using Simpson  method.
{
    size_t idx_x,x_Max;
    complex<double> sum_z=0.0,sum_last=0.0;
    double delta_3=delta/3.0;

    if ((idx_t-2*idx_z)%2==0)
    {
        x_Max=idx_t-2*idx_z;
    }
    else{
        x_Max=idx_t-2*idx_z-1;
        //Since Simpson method needs even number of nodes we need to check if t-2*z was odd than to apply trapezoid method for the last segment.
        sum_last=0.5*delta*(FuncForX(A,x_Max+1,idx_z, idx_t, trig,delta)+
                FuncForX(A,x_Max,idx_z, idx_t, trig,delta));
    }

    sum_z= FuncForX(A,0,idx_z, idx_t, trig,delta)+
            FuncForX(A,x_Max,idx_z, idx_t, trig,delta);
//This loop implements Simpson method.

    for (idx_x=1;idx_x<x_Max;idx_x++)
    {
        if ((idx_x)%2==1){
            sum_z+= 4.0*FuncForX(A,idx_x,idx_z, idx_t,trig,delta);
        }
        else{
            sum_z+= 2.0*FuncForX(A,idx_x,idx_z, idx_t,trig,delta);
        }
       //  if (idx_x%20==0)
         //   {
              // prsize_tf("A2 beforeR=%lf\n",cimag(sum_z));
            //    sum_z=Mround(sum_z);
             //  prsize_tf("A2 afterR=%lf\n",cimag(sum_z));
          //  }
    }

    sum_z=sum_z*delta_3+sum_last;

   // if(cabs(sum_z)<error_x1){
          //  prsize_tf("X1_sim: sum=%f, idx_x=%i, idx_z=%i, idx_t=%i; x_Max=%i\n",cabs(sum_z),idx_x, idx_z, idx_t,x_Max);
          //  sum_z=0.0;}
 
return sum_z;
}

/**************************************************************************************************/
complex<double> FuncForZ (complex<double> *A, size_t idx_z, size_t idx_t,double delta)

//Same as before; this function will calculate the value of function (z) to be size_tegrated by z for z=delta*idx_z using trapezoid method.
{
    complex<double> _Complex_I(0.0,1.0); 
    complex<double>  X_int1, X_int2,sum;
    double  z2, z3;

    double delta2=delta*delta;
    double delta3=delta2*delta;
    double nu1_3=nu1*nu1*nu1;
    double nu2_3=nu2*nu2*nu2;
    double  nu1_3_2=(nu1_3*2.0)/3.0;
    double  nu2_3_2=(nu2_3*2.0)/3.0;
   // double  nu2_3_20=nu2_3*10.0;
    complex<double> ialpha1_2=_Complex_I*alpha1*alpha1;
    complex<double> ialpha2_2=_Complex_I*alpha2*alpha2;



    z3=delta3*(double)idx_z*(double)idx_z*(double)idx_z;
    z2=delta2*(double)idx_z*(double)idx_z;


        X_int1=X_integral(A,idx_z, idx_t,1,delta);
	X_int2=X_integral(A,idx_z, idx_t,2,delta);

        complex<double> arg=ialpha1_2*z2-nu1_3_2*z3-beta1*2.0*idx_z*delta;
        
        sum=c1*z2*A[idx_t-idx_z]*X_int1*exp(arg);
        
        arg=ialpha2_2*z2-nu2_3_2*z3-beta2*2.0*idx_z*delta;
        
         sum+=c2*z2*A[idx_t-idx_z]*X_int2*exp(arg);
      
      //    sum=z2*A[idx_t-idx_z]*(X_int1*c1*cexp(arg)*exp(-nu3_2*z3));
      //    sum=z2*A[idx_t-idx_z]*(X_int1*c1*cexp(arg)*exp(-nu3_2*z3));

       // if(cabs(sum)<error_z){
              //   prsize_tf("Func_Z: sum=%f, idx_z=%i, idx_t=%i; \n",cabs(sum),idx_z, idx_t);
      //  sum=0.0;}

return sum;
}

/**************************************************************************************************/

complex<double> callSimpsonX(complex<double> *A, size_t idx_z1, size_t idx_t,size_t trig,double delta)

//Simpson method works for n>2 posize_ts. If we have less then 2 posize_ts than we use trapezoid method.
{
    if (idx_t-2*idx_z1<2){     
        return X_integral(A,idx_z1, idx_t,trig,delta);}
    
        else{
        return SimpsonX(A,idx_z1, idx_t,trig,delta);}

}
/**************************************************************************************************/
complex<double> FuncForZSimp (complex<double> *A, size_t idx_z, size_t idx_t,double delta)

//This function will calculate the value of function (z) to be size_tegrated by z for z=delta*idx_z using Simpson method.
{
    complex<double> _Complex_I(0.0,1.0); 
    complex<double>  X_int1, X_int2,sum=0.0;
    double delta2=delta*delta;
    double delta3=delta2*delta;
    double nu1_3=nu1*nu1*nu1;
    double nu2_3=nu2*nu2*nu2;
    double  nu1_3_2=(nu1_3*2.0)/3.0;
    double  nu2_3_2=(nu2_3*2.0)/3.0;
//    double  nu2_3_20=nu2_3_2*10.0;
    complex<double> ialpha1_2=_Complex_I*alpha1*alpha1;
    complex<double> ialpha2_2=_Complex_I*alpha2*alpha2;



        double z3=delta3*(double)idx_z*(double)idx_z*(double)idx_z;
        double z2=delta2*(double)idx_z*(double)idx_z;

        X_int1=callSimpsonX(A,idx_z, idx_t,1,delta);
	X_int2=callSimpsonX(A,idx_z, idx_t,2,delta);

        complex<double> arg;
        arg=ialpha1_2*z2-nu1_3_2*z3-beta1*2.0*idx_z*delta;

        sum=c1*z2*A[idx_t-idx_z]*X_int1*exp(arg);
        
        arg=ialpha2_2*z2-nu2_3_2*z3-beta2*2.0*idx_z*delta;
        
        sum+=c2*z2*A[idx_t-idx_z]*X_int2*exp(arg);
       

   // if(cabs(sum)<error_z){
           //  prsize_tf("Func_Z_sim:sum=%f,idx_z=%i, idx_t=%i; \n",cabs(sum), idx_z, idx_t);
  //  sum=0.0;}

return sum;
}

/**************************************************************************************************/

complex<double> Z_integral (complex<double> *A, size_t idx_t,double delta)
{
//This function will calculate z-size_tegral by trapezoid method.
    complex<double> sum_z=0.0;
	size_t idx_z, idx_t2;

    idx_t2=idx_t>>1;

         if(idx_t2>1E-15){
          
            for(idx_z=1;idx_z<idx_t2;idx_z++){
                
                sum_z+=delta*( FuncForZ(A,idx_z,idx_t,delta) );
                  
            }
            
                sum_z+=0.5*delta*( FuncForZ(A,idx_t2, idx_t,delta) );
             
            if ((idx_t%2)==1){
                 
               sum_z+=delta*(FuncForZ(A,idx_t2, idx_t, delta))*0.5;
            }
        }
      //  if(cabs(sum_z)<error_z){
              //   prsize_tf("Z_I:sum=%f, idx_z=%i, idx_t=%i; \n",cabs(sum_z), idx_z, idx_t);
      //  sum_z=0.0;}
    
return sum_z;
}

/**************************************************************************************************/

complex<double>  SimpsonZ(complex<double> *A, size_t idx_t,double delta)
{
    //This function will calculate z-size_tegral by Simpson method. Same as for x-size_tegral we need to check if we have odd or even number of nodes.
    size_t idx_z, z_Max=0;
    complex<double> sum_z=0.0,sum_last=0.0;
    double delta_3=delta/3.0;

    size_t idx_t2=idx_t>>1;

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

   // if(cabs(sum_z)<error_z){
     //   sum_z=0.0;
      //  prsize_tf("Z_sim: sum=%f, idx_z=%i, idx_t=%i; z_Max=%i\n",cabs(sum_z), idx_z, idx_t,z_Max);
      //  }
return sum_z;
}

/**************************************************************************************************/

complex<double> callSimpsonZ (complex<double> *A, size_t idx_t,double delta)
{
    //same as for X: we need to use trapezoid method if we do not have enough nodes.

//prsize_tf("delta=%f, nu=%f, alpha=%f, delta2=%f, delta3=%f, delta_3=%f, nu3=%f, nu3_2=%f, nu3_20=%f, \n",delta,nu, alpha, delta2,delta3,delta_3,nu3, nu3_2,nu3_20);
    if (idx_t<4){
        return Z_integral(A,idx_t,delta);
    }
    else {
        return SimpsonZ(A,idx_t,delta);
    }
}
/**************************************************************************************************/
void  simple_integration (complex<double> *A, size_t idx,double delta)
{
    //Calculating the main differential equation for t=delta*idx by trapezoid methond on size_terval=(idx-1)*delta - idx*delta.

	A[idx]=(  A[idx-1]*(1.+0.5*delta) - 0.25*delta*
                 ( callSimpsonZ(A,idx,delta) + callSimpsonZ(A,idx-1,delta)) )/(1.0-0.5*delta);
}
/**************************************************************************************************/

void  Simpson_integrationT(complex<double> *A, size_t idx,double delta)
//Calculating the main differential equation for t=delta*idx by Simpson methond on size_terval t=(idx-2)*delta - idx*delta.
 //This function gave better results for C2=1, C1=0 so we decided  to use it in main().
{
complex<double> sum;
double delta_3=delta/3.;
double delta_3_8=3.*delta/8.;
double delta_24=delta/24.;
double delta_48=delta/48.;

     if (idx<3){
        
        sum=-(0.5*delta_3)*(  callSimpsonZ(A,idx,delta) +  4.*(  callSimpsonZ(A,idx-1,delta)  )+
                                    callSimpsonZ(A,idx-2,delta) );
        
       // printf("Sim- idx=%zu: sum=%f\n",idx, real(sum));
        
        sum+= A[idx-2]*(1.0+delta_3);
        sum+=(4.* delta_3)*(A[idx-1]);
        A[idx]=sum/(1.-delta_3);

    }
    if ((idx>=3)&& (idx<6))
  // if (idx>=3)
        {
           
        sum=-(0.5*delta_3_8)*(
                              callSimpsonZ(A,idx,delta) +
                              3.*(  callSimpsonZ(A,idx-1,delta) + callSimpsonZ(A,idx-2,delta) )+
                                    callSimpsonZ(A,idx-3,delta)
                               );
         
         sum+=A[idx-3]*(1.+delta_3_8);
         sum+=(3.*delta_3_8)*(A[idx-2]+A[idx-1]);
        A[idx]=sum/(1.-delta_3_8);

    }

     if ((idx>=6)&& (idx<8))
        {
        
        sum=-(0.5*delta_24)*(  9.*(callSimpsonZ(A,idx, delta) +callSimpsonZ(A,idx-6, delta))+
                               28.*(  callSimpsonZ(A,idx-1, delta) + callSimpsonZ(A,idx-5, delta) )+
                                  23.*(  callSimpsonZ(A,idx-2, delta) + callSimpsonZ(A,idx-4, delta) )+
                             24.*callSimpsonZ(A,idx-3, delta));

         printf("Sim- idx=%zu: sum=%f\n",idx, real(sum));
        
        sum+=A[idx-6]*(1.+9.*delta_24);
         sum+=(delta_24)*(28.*A[idx-5]+23.*A[idx-4]+24.*A[idx-3]+23.*A[idx-2]+28.*A[idx-1]);
         
        A[idx]=sum/(1.-9.*delta_24);
       // prsize_tf("I>3: i=%i, r=%f\n",idx,cabs(A[idx]));
    }
if (idx>=8)
        {
        
        sum=-(0.5*delta_48)*(  17.*(callSimpsonZ(A,idx, delta) +callSimpsonZ(A,idx-8, delta))+
                               59.*(  callSimpsonZ(A,idx-1, delta) + callSimpsonZ(A,idx-7, delta) )+
                               43.*(  callSimpsonZ(A,idx-2, delta) + callSimpsonZ(A,idx-6, delta) )+
                               49.*(  callSimpsonZ(A,idx-3, delta) + callSimpsonZ(A,idx-5, delta) )+
                             48.*callSimpsonZ(A,idx-4, delta));
        
         printf("Sim- idx=%zu: sum=%f\n",idx, real(sum));
        
        sum+=A[idx-8]*(1.+17.*delta_48);
         sum+=(delta_48)*(59.*A[idx-7]+43.*A[idx-6]+49.*A[idx-5]+48.*A[idx-4]+49.*A[idx-3]+43.*A[idx-2]+59.*A[idx-1]);
         
        A[idx]=sum/(1.-17.*delta_48);
       // prsize_tf("I>3: i=%i, r=%f\n",idx,cabs(A[idx]));
    }
}
