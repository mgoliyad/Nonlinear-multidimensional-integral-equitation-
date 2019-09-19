
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
#include <unistd.h>


#include "newNonComp.h"
#include "algorithm.h"
#include <complex>       
using namespace std;


double  Simpson_extendedN(double (*f)(size_t , size_t, size_t ),
                                    size_t idx_t,size_t idx_z , size_t idx_x )
{
    double delta_48=delta/48.;
    
        return delta_48*(  17.*(f(idx_t,idx_z, idx_x)  + f(idx_t,idx_z, idx_x-8))+
                               59.*(f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z, idx_x-7) )+
                               43.*(f(idx_t,idx_z,idx_x-2) + f(idx_t,idx_z, idx_x-6) )+
                               49.*(f(idx_t,idx_z,idx_x-3) + f(idx_t,idx_z, idx_x-5) )+
                               48.*f(idx_t,idx_z,idx_x-4));

}

double  Simpson_extendedTN(double (*f)(size_t , size_t ),
                                    size_t idx_t,size_t idx_z  )
{
    double delta_48=delta/48.;
    
        return delta_48*(  17.*(f(idx_t,idx_z)  + f(idx_t,idx_z-8))+
                               59.*(f(idx_t,idx_z-1) + f(idx_t,idx_z-7) )+
                               43.*(f(idx_t,idx_z-2) + f(idx_t,idx_z-6) )+
                               49.*(f(idx_t,idx_z-3) + f(idx_t,idx_z-5) )+
                               48.*f(idx_t,idx_z-4));

}


double  BooleN(double (*f)(size_t ,size_t , size_t ), 
            size_t idx_t,size_t idx_z, size_t idx_x)
{

    double delta_2_45=2.*delta/45.;
    
    return delta_2_45*(
                              7.*(  f(idx_t,idx_z,idx_x) +f(idx_t,idx_z,idx_x-4) )+
                
                              32.*(  f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-3) )+
                
                                    12.*f(idx_t,idx_z,idx_x-2)
                               );

}

double  BooleTN(double (*f)(size_t ,size_t  ), 
            size_t idx_t,size_t idx_z)
{

    double delta_2_45=2.*delta/45.;
    
    return delta_2_45*(
                              7.*(  f(idx_t,idx_z) +f(idx_t,idx_z-4) )+
                
                              32.*(  f(idx_t,idx_z-1) + f(idx_t,idx_z-3) )+
                
                                    12.*f(idx_t,idx_z-2)
                               );

}

double  Simpson_narrowN(double (*f)(size_t ,size_t , size_t ), 
            size_t idx_t,size_t idx_z, size_t idx_x)
{
    double delta_24=delta/24.;
    
            return delta_24*(  9.*(f(idx_t,idx_z,idx_x) +f(idx_t,idx_z,idx_x-6))+
                               28.*(f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-5) )+
                               23.*(f(idx_t,idx_z,idx_x-2) + f(idx_t,idx_z,idx_x-4) )+
                               24.*f(idx_t,idx_z,idx_x-3));

}

double  Simpson_narrowTN(double (*f)(size_t ,size_t  ), 
            size_t idx_t,size_t idx_z)
{
    double delta_24=delta/24.;
    
            return delta_24*(  9.*(f(idx_t,idx_z) +f(idx_t,idx_z-6))+
                               28.*(f(idx_t,idx_z-1) + f(idx_t,idx_z-5) )+
                               23.*(f(idx_t,idx_z-2) + f(idx_t,idx_z-4) )+
                               24.*f(idx_t,idx_z-3));

}



double  Simpson_3_8N(double (*f)(size_t ,size_t , size_t ), 
            size_t idx_t,size_t idx_z, size_t idx_x)
{
    double delta_3_8=3.*delta/8.;
    
    return delta_3_8*(
                              f(idx_t,idx_z,idx_x) +
                              3.*(  f(idx_t,idx_z,idx_x-1) + f(idx_t,idx_z,idx_x-2) )+
                                    f(idx_t,idx_z,idx_x-3)
                               );

}
double  Simpson_3_8TN(double (*f)(size_t ,size_t  ), 
            size_t idx_t,size_t idx_z)
{
    double delta_3_8=3.*delta/8.;
    
    return delta_3_8*(
                              f(idx_t,idx_z) +
                              3.*(  f(idx_t,idx_z-1) + f(idx_t,idx_z-2) )+
                                    f(idx_t,idx_z-3)
                               );

}



double  Simpson_ruleN(double (*f)(size_t , size_t, size_t ), 
                                     size_t idx_t,size_t idx_z , size_t idx_x )
{
double delta_3=delta/3.;

        return delta_3*(  f(idx_t,idx_z,idx_x) +  4.*(  f(idx_t,idx_z,idx_x-1)  )+
                                    f(idx_t,idx_z,idx_x-2) );
        
}
 
double  Simpson_ruleTN(double (*f)(size_t , size_t ), 
                                     size_t idx_t,size_t idx_z  )
{
double delta_3=delta/3.;
double s1, s2, s3;
s1=f(idx_t,idx_z);
s2=f(idx_t,idx_z-1);
s3=f(idx_t,idx_z-2);
double res=delta_3*(s1+4.*s2+s3);
//printf("Sim: %f, %f, %f, d=%f, res=%f\n",s1,s2,s3,delta_3,res);
fflush(stdout);
return res;
       // return delta_3*(  f(idx_t,idx_z) +  4.*(  f(idx_t,idx_z-1)  )+
         //                           f(idx_t,idx_z-2) );
        
}

double  Trapezoid_ruleN(double (*f)(size_t , size_t, size_t ), 
                                     size_t idx_t,size_t idx_z , size_t idx_x )
{
    double sum=0.0;
    double delta_2=delta/2.;
    
        return delta_2*(  f(idx_t,idx_z, idx_x) +
                          f(idx_t,idx_z,idx_x-1) );

}

double  Trapezoid_ruleTN(double (*f)(size_t , size_t ), 
                                     size_t choice, size_t idx_t,char flag )
{
    double sum=0.0;
    double delta_2=delta/2.;
    
      /*  return delta_2*(  f(idx_t,idx_z) +
                          f(idx_t,idx_z-1) );*/
        
    if (flag=='r'){   
        
        double  help=-0.5*(  f(choice,idx_t) +
                              f(choice,idx_t-1));
    printf("doing r\n");
        printf ("N A_1=%f, help=%f\n", real(r[idx_t-1]), help);
        
    /*   sum=-(0.5*delta_2)*(  f(choice,idx_t) +
                              f(choice,idx_t-1));*/
        sum=delta_2*help;
        printf ("N sum=%f\n",(sum));
        sum+= r[idx_t-1]*(1.0+delta_2);
        printf ("N sum=%f\n",(sum));
        sum/=(1.-delta_2);
printf ("N sum=%f\n",(sum));
    }
    else{
        printf("doing w\n");
       // double help1=f(choice,idx_t);
        double help2=f(choice,idx_t-1);
        printf ("N  help2=%f\n", help2);
        sum=help2;
        printf ("N sum=%f\n",(sum));
        sum*=-0.5*delta;
        printf ("N sum=%f\n",(sum));
      /*  sum=delta_2*(  f(choice,idx_t) +
                              f(choice,idx_t-1));*/
        
    }
return sum;
}
