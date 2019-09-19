#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
#include <unistd.h>


#include "newComp.h"
#include "algorithm.h"
#include "Romb.h"
#include <complex>       
using namespace std;

complex<double>  Simpson_extended(complex<double> (*f)(int , int, int ),
                                    int idx_t,int idx_z , int idx_x )
{
    double delta_48=delta/48.;
    
        return delta_48*(  17.*(f(idx_t,idx_z, idx_x)  + f(idx_t,idx_z, idx_x-8))+
                               59.*(f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z, idx_x-7) )+
                               43.*(f(idx_t,idx_z,idx_x-2) + f(idx_t,idx_z, idx_x-6) )+
                               49.*(f(idx_t,idx_z,idx_x-3) + f(idx_t,idx_z, idx_x-5) )+
                               48.*f(idx_t,idx_z,idx_x-4));

}
complex<double>  Simpson_extendedT(complex<double> (*f)(int ,int ), 
                                    int choice, int idx_t )

{
    complex<double> sum=0.0;
    double delta_48=delta/48.;
    
        sum=-(0.5*delta_48)*(  17.*(f(choice, idx_t)  + f(choice, idx_t-8))+
                               59.*(f(choice, idx_t-1) + f(choice,idx_t-7) )+
                               43.*(f(choice, idx_t-2) + f(choice,idx_t-6) )+
                               49.*(f(choice, idx_t-3) + f(choice,idx_t-5) )+
                               48.*f(choice,  idx_t-4));
      
        sum+=A[idx_t-8]*(1.+17.*delta_48);
        sum+=(delta_48)*(59.*A[idx_t-7]+43.*A[idx_t-6]+49.*A[idx_t-5]+48.*A[idx_t-4]+49.*A[idx_t-3]+43.*A[idx_t-2]+
                59.*A[idx_t-1]);
        
        sum/=(1.-17.*delta_48);
        return sum;
}
complex<double>  Simpson_extendedR(complex<double> (*f)(int ,int ), int idx_t )
{
    complex<double> sum=0.0;
    double delta_48=delta/48.;
    printf("RR\n");
      sum=-0.5*f(8, idx_t);
      printf("isum=%f\n",abs(sum));
      
        sum+=A[idx_t-8]*(1.+17.*delta_48);
        sum+=(delta_48)*(59.*A[idx_t-7]+43.*A[idx_t-6]+49.*A[idx_t-5]+48.*A[idx_t-4]+49.*A[idx_t-3]+43.*A[idx_t-2]+
                59.*A[idx_t-1]);
      //  printf("isum=%f\n",abs(sum));
        sum/=(1.-17.*delta_48);
    
        return sum;
}
/*******************************************************************************************************************/

complex<double>  Boole(complex<double> (*f)(int ,int , int ), 
            int idx_t,int idx_z, int idx_x)
{

    double delta_2_45=2.*delta/45.;
    
    return delta_2_45*(
                              7.*(  f(idx_t,idx_z,idx_x) +f(idx_t,idx_z,idx_x-4) )+
                
                              32.*(  f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-3) )+
                
                                    12.*f(idx_t,idx_z,idx_x-2)
                               );

}

complex<double>  BooleT(complex<double> (*f)(int ,int ), int choice, int idx_t)
{
    complex<double> sum=0.0;
    double delta_2_45=2.*delta/45.;
    
    
    sum=-(0.5*delta_2_45)*(
                              7.*(f(choice, idx_t) +    f(choice, idx_t-4) )+
                              32.*(f(choice, idx_t-1) + f(choice, idx_t-3) )+
                              12.*f(choice, idx_t-2)
                               );
         
         sum+=A[idx_t-4]*(1.+7.*delta_2_45);
         sum+=(32.*delta_2_45)*(A[idx_t-3]+A[idx_t-1]);
         sum+=(12.*delta_2_45)*(A[idx_t-2]);
        sum/=(1.-7.*delta_2_45);
        
    return sum;
}

complex<double>  BooleR(complex<double> (*f)(int ,int ),  int idx_t)
{
    complex<double> sum=0.0;
    double delta_2_45=2.*delta/45.;
    
    
    sum=-0.5*f(4, idx_t);
         
         sum+=A[idx_t-4]*(1.+7.*delta_2_45);
         sum+=(32.*delta_2_45)*(A[idx_t-3]+A[idx_t-1]);
         sum+=(12.*delta_2_45)*(A[idx_t-2]);
        sum/=(1.-7.*delta_2_45);
        
    return sum;
}

/*******************************************************************************************************************/
complex<double>  Simpson_narrow(complex<double> (*f)(int ,int , int ), 
            int idx_t,int idx_z, int idx_x)
{
    double delta_24=delta/24.;
    
            return delta_24*(  9.*(f(idx_t,idx_z,idx_x) +f(idx_t,idx_z,idx_x-6))+
                               28.*(f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-5) )+
                               23.*(f(idx_t,idx_z,idx_x-2) + f(idx_t,idx_z,idx_x-4) )+
                               24.*f(idx_t,idx_z,idx_x-3));

}

complex<double>  Simpson_narrowT(complex<double> (*f)(int ,int  ), int choice, int idx_t)
{
    complex<double> sum=0.0;
    double delta_24=delta/24.;
    
    sum=-(0.5*delta_24)*(  9.*(f(choice, idx_t) +f(choice, idx_t-6))+
                               28.*(  f(choice, idx_t-1) + f(choice, idx_t-5) )+
                               23.*(  f(choice, idx_t-2) + f(choice, idx_t-4) )+
                               24.*f(choice, idx_t-3));
 
        sum+=A[idx_t-6]*(1.+9.*delta_24);
         sum+=(delta_24)*(28.*A[idx_t-5]+23.*A[idx_t-4]+24.*A[idx_t-3]+23.*A[idx_t-2]+28.*A[idx_t-1]);
         
        sum/=(1.-9.*delta_24);
        return sum;
}

complex<double>  Simpson_narrowR(complex<double> (*f)(int ,int ),int idx_t )
{
    complex<double> sum=0.0;
    double delta_24=delta/24.;
  
        sum=-0.5*f(6, idx_t);
  
        sum+=A[idx_t-6]*(1.+9.*delta_24);
         sum+=(delta_24)*(28.*A[idx_t-5]+23.*A[idx_t-4]+24.*A[idx_t-3]+23.*A[idx_t-2]+28.*A[idx_t-1]);
         
        sum/=(1.-9.*delta_24);
        return sum;
}

/*******************************************************************************************************************/

complex<double>  Simpson_3_8(complex<double> (*f)(int ,int , int ), 
            int idx_t,int idx_z, int idx_x)
{
    double delta_3_8=3.*delta/8.;
    
    return delta_3_8*(
                              f(idx_t,idx_z,idx_x) +
                              3.*(  f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-2) )+
                                    f(idx_t,idx_z,idx_x-3)
                               );

}

complex<double>  Simpson_3_8T(complex<double> (*f)(int ,int  ), int choice, int idx_t)
{
    complex<double> sum=0.0;
    double delta_3_8=3.*delta/8.;
  // printf("3_8T\n");
    sum=-0.5*delta_3_8*(f(choice, idx_t) +
                  3.*(  f(choice, idx_t-1) + 
                        f(choice, idx_t-2) )+
                        f(choice, idx_t-3));

    sum+=A[idx_t-3]*(1.+delta_3_8)+(3.*delta_3_8)*(A[idx_t-2]+A[idx_t-1]);
    sum/=(1.-delta_3_8);
       // printf("sum=%f\n",abs(sum));
        return sum;
}

complex<double>  Simpson_3_8R(complex<double> (*f)(int ,int ),int idx_t )
{
    complex<double> sum=0.0;
    double delta_3_8=3.*delta/8.;
 // printf("3_8R\n");
     sum=-0.5*f(3, idx_t);
     
     sum+=A[idx_t-3]*(1.+delta_3_8)+(3.*delta_3_8)*(A[idx_t-2]+A[idx_t-1]);
     sum/=(1.-delta_3_8);
      //  printf("sum=%f\n",abs(sum));
        return sum;
}
/*******************************************************************************************************************/

complex<double>  Simpson_rule(complex<double> (*f)(int , int, int ), 
                                     int idx_t,int idx_z , int idx_x )
{
double delta_3=delta/3.;

        return delta_3*(  f(idx_t,idx_z,idx_x) +  4.*(  f(idx_t,idx_z,idx_x-1)  )+
                                    f(idx_t,idx_z,idx_x-2) );
        
}

complex<double>  Simpson_ruleT(complex<double> (*f)(int ,int  ), int choice, int idx_t)
{
    complex<double> sum=0.0;
    double delta_3=delta/3.;
   // printf("Simpson:\n");
        sum=-(0.5*delta_3)*(  f(choice,idx_t) +  4.*(  f(choice,idx_t-1)  )+
                                    f(choice,idx_t-2) );
        sum+= A[idx_t-2]*(1.0+delta_3)+(4.* delta_3)*(A[idx_t-1]);

        sum/=(1.-delta_3);
        return sum;
}
complex<double>  Simpson_ruleR(complex<double> (*f)(int ,int ),int idx_t )
{
    complex<double> sum=0.0;
    double delta_3=delta/3.;
    
         sum=-0.5*f(2, idx_t);
         sum+= A[idx_t-2]*(1.0+delta_3)+(4.* delta_3)*(A[idx_t-1]);

        sum/=(1.-delta_3);        
        return sum;


}
/*******************************************************************************************************************/

complex<double>  Trapezoid_rule(complex<double> (*f)(int , int, int ), int idx_t,int idx_z , int idx_x )
{
    complex<double> sum=0.0;
    double delta_2=delta/2.;

        return delta_2*(  f(idx_t,idx_z, idx_x) +
                          f(idx_t,idx_z,idx_x-1) );

}

complex<double>  NoInt(complex<double> (*f)(int , int, int ), int idx_t,int idx_z , int idx_x )
{
    complex<double> sum=0.0;
    
        return f(idx_t,idx_z, idx_x);

}

complex<double>  Trapezoid_ruleT(complex<double> (*f)(int ,int ), 
                                    int choice, int idx_t )
{
    complex<double> sum=0;
    double delta_2=delta/2.;

    complex<double> f_f=f(choice,idx_t);
    complex<double> f_f1=f(choice,idx_t-1);    
    sum=(A[idx_t-1]*(1.0+delta_2)-(0.5*delta_2)*(  f_f +f_f1))/(1.-delta_2);
    
    return sum;

}
complex<double>  NoIntT(complex<double> (*f)(int ,int ), 
                                    int choice, int idx_t )
{
    complex<double> sum=0;
    complex<double> f_f=f(choice,idx_t) / 2.;
    sum=A[idx_t-1] / delta;
    complex<double> sum1=sum;
    sum -= f_f;
    sum/=(1./delta-1.);
    return sum;

}

complex<double>  Trapezoid_ruleR(complex<double> (*f)(int,int  ),int idx_t)
{
    complex<double> sum=0.0;
    double delta_2=delta/2.;

        sum=-0.5*f(1,idx_t);
        fflush(stdout);
        sum+=A[idx_t-1]*(1.0+delta_2);
        sum/=(1-0.5*delta);
    return sum;
   
}

complex<double>  NoIntR(complex<double> (*f)(int,int  ),int idx_t)
{
    complex<double> sum=0.0;
    double delta_2=delta/2.;
    sum=(A[idx_t-1]/delta-0.5*f(0,idx_t))/(1.0/delta-1.0);
    return sum;
}

complex<double>  Trapezoid_ruleR_new(complex<double> (*f)(int  ),int idx_t)
{
    complex<double> sum=0.0;
    sum=-0.5*f(idx_t);
    sum+=A[idx_t-2]*exp(delta+delta);
    return sum;
}

