#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>
#include "Comp.h"
#include "newNonComp.h"
#include <complex>   
#include "algorithmN.h"
using namespace std;

double phi1(size_t idx_z, size_t idx_x)
{
    double result;
    double z=(double)idx_z*delta;
    double z2=z*z;
    double x=(double)idx_x*delta;
//    double nu1_3=nu*nu*nu;
    result=(2.0*z/3.0+x);
    result*=z2*nu1_3;
   return result;
}
double phi2(size_t idx_z, size_t idx_x)
{
    double result;
    double z=(double)idx_z*delta;
    double z2=z*z;
    double x=(double)idx_x*delta;
//    double nu1_3=nu*nu*nu;
    result=(2.0*z/3.0+x);
    result*=z2*nu2_3;
   return result;
}
double psi1(size_t idx_z, size_t idx_x)
{
    double  result;
    double alpha1_2=alpha1*alpha1;
    double z=(double)idx_z*delta;
    double x=(double)idx_x*delta;
    result=alpha1_2*z*(x+z);
  //  printf("psi1: %f\n",result);
return result;
}
double psi2(size_t idx_z, size_t idx_x)
{
  /* double  result;
    double alpha2_2=alpha2*alpha2;
    double z=(double)idx_z*delta;
    double x=(double)idx_x*delta;
    result=alpha2_2*z*(x+z);
  //  printf("psi2: %f\n",result);
return result;*/
    return 0.0;
}

double  Omega1 ( size_t idx_t, size_t idx_z,size_t idx_x)
{
   double result;
   size_t i1=idx_t-idx_z;
    size_t i2=idx_t-idx_z-idx_x;
    size_t i3=idx_t-idx_z-idx_z-idx_x;
   result= w[i2]+w[ i2]-w[ i3];
   result+=psi1(idx_z,idx_x);
  // if (idx_t==51) {printf("Om1: %f\n",result);}
   return result;
}
double  Omega2 (size_t idx_t, size_t idx_z,size_t idx_x)
{
   double result;
   size_t i1=idx_t-idx_z;
    size_t i2=idx_t-idx_z-idx_x;
    size_t i3=idx_t-idx_z-idx_z-idx_x;
   result= w[i2]+w[ i2]-w[ i3];
   result+=psi2(idx_z,idx_x);
 // if (idx_t==51) {printf("Om2: %f\n",result);}
   return result;
}

double R(size_t idx_t, size_t idx_z, size_t idx_x)
{double  result;
size_t i1=idx_t-idx_z;
size_t i2=idx_t-idx_z-idx_x;
size_t i3=idx_t-idx_z-idx_z-idx_x;

double r1=r[i1];
double r2=r[i2];
double r3=r[i3];
result=r1*r2*r3;
   // result=r[idx_t-idx_z]*r[idx_t-idx_z-idx_x]*r[idx_t-2*idx_z-idx_x];
return result;
}
/************************************************/
double FuncForXN1 (size_t idx_t,size_t idx_z, size_t idx_x)
{
    double  result, result1, result2;
    double z2=delta*delta*(double)idx_z*(double)idx_z;
  
     //   result1=exp(-phi1(idx_z,idx_x)-beta1*(double)(2*idx_z+idx_x)*delta);
    result1=exp(-phi1(idx_z,idx_x));
       //  if (idx_t==51){ printf("XN1 s1=%f, \n",result1);}
        result1*=R(idx_t,idx_z,idx_x);
       //  if (idx_t==51){ printf("XN1 s1=%f, \n",result1);}
        result1*=cos(Omega1(idx_t,idx_z,idx_x));
    // if (idx_t==51){ printf("XN1 s1=%f, om=%f\n",result1,Omega1(idx_t,idx_z,idx_x) );}
        
      //  result2=exp(-phi2(idx_z,idx_x)-beta2*(double)(2*idx_z+idx_x)*delta);
        result2=exp(-phi2(idx_z,idx_x));
      //   if (idx_t==51){ printf("XN1-2 s1=%f, \n",result);}
        result2*=R(idx_t,idx_z,idx_x);
      //  if (idx_t==51){ printf("XN1-2 s1=%f, \n",result);}
        result2*=cos(Omega2(idx_t,idx_z,idx_x));
      //  if (idx_t==51){ printf("XN1-2 s1=%f,om=%f\n",result1,Omega2(idx_t,idx_z,idx_x) );}
    
    result=(c1*result1+c2*result2);
  // printf("XN1 s1=%f, s2=%f, s=%f, t=%zu, z=%zu, x=%zu\n",result, result1, result2, idx_t, idx_z,  idx_x);
return result;
}
double FuncForXN2 (size_t idx_t,size_t idx_z, size_t idx_x)
{
    double  result, result1, result2;
    
   
      //  result1=exp(-phi1(idx_z,idx_x)-beta1*(double)(2*idx_z+idx_x)*delta);
    result1=exp(-phi1(idx_z,idx_x));
      //  if (idx_t==51){ printf("XN2-1 s1=%f, \n",result);}
        result1*=R(idx_t,idx_z,idx_x);
       // if (idx_t==51){ printf("XN2-1 s1=%f, \n",result);}
        result1*=sin(Omega1(idx_t,idx_z,idx_x));
       // if (idx_t==51){ printf("XN2-1 s1=%f, om=%f\n",result, Omega1(idx_t,idx_z,idx_x));}
    
       // result2=exp(-phi2(idx_z,idx_x)-beta2*(double)(2*idx_z+idx_x)*delta);
        result2=exp(-phi2(idx_z,idx_x));
      //  if (idx_t==51){ printf("XN2-2 s1=%f, \n",result);}
        result2*=R(idx_t,idx_z,idx_x);
      //  if (idx_t==51){ printf("XN2-2 s1=%f, \n",result);}
        result2*=sin(Omega2(idx_t,idx_z,idx_x));
       // if (idx_t==51){ printf("XN2-2 s1=%f, om2=%f\n",result, Omega2(idx_t,idx_z,idx_x));}
    
    
     result=(c1*result1+c2*result2);
return result;
}

double FuncForZN1 (size_t choiceX, size_t idx_t,size_t idx_z)
{
  //  printf("ZN1: chX=%zu, t=%zu, z=%zu\n",choiceX,idx_t, idx_z);
    if (idx_z==0){return 0.0;}
    double sum=0.0;
 double z2=delta*delta*(double)idx_z*(double)idx_z;
    size_t idx=1;
        
    while (idx*choiceX<=idx_t-idx_z-idx_z) {
        sum+=callMethodXN(FuncForXN1,idx_t, idx_z,idx*choiceX,choiceX);      
        idx++;       
        
    }
    
    size_t new_choiceX=idx_t-idx_z-idx_z-(idx-1)*choiceX;
    if (new_choiceX!=0){ 
        sum+=callMethodXN(FuncForXN1,idx_t, idx_z, idx_t-idx_z-idx_z, new_choiceX);
        
    }
 // if (idx_t==51){  printf("zN1s=%f, s2=%f, s=%f\n",sum, sum, sum);}
    return z2*sum;
}
double FuncForZN2 (size_t choiceX, size_t idx_t,size_t idx_z)
{
   // printf("ZN2: chX=%zu, t=%zu, z=%zu\n",choiceX,idx_t, idx_z);
    if (idx_z==0){return 0.0;}
    double sum=0.0;
 double z2=delta*delta*(double)idx_z*(double)idx_z;
    size_t idx=1;
        
    while (idx*choiceX<=idx_t-idx_z-idx_z) {
        sum+=callMethodXN(FuncForXN2,idx_t, idx_z,idx*choiceX,choiceX);      
        idx++;   
        
    }
    
    size_t new_choiceX=idx_t-idx_z-idx_z-(idx-1)*choiceX;
    if (new_choiceX!=0){ 
        sum+=callMethodXN(FuncForXN2,idx_t, idx_z, idx_t-idx_z-idx_z, new_choiceX);
        
    }
  // if (idx_t==20){  printf("zN2s1=%f, s2=%f, s=%f\n",sum, sum, sum);}
return z2*sum;
}


double FuncForTU ( size_t choiceZ, size_t idx_t)
{
   // printf("TU: chX=%zu, t=%zu, \n",choiceZ,idx_t);
    double sum, sum1=0, sum2=0;
    double delta2=delta*delta;
    double t=(double)idx_t*delta;
    
     size_t idx_t2=floor(idx_t/2);

    size_t idx=1;
     
    while (idx*choiceZ<=idx_t2) {
        
        sum1+=callMethodZN(FuncForZN1,choiceZ, idx_t, idx*choiceZ);
        sum2+=callMethodZN(FuncForZN2,choiceZ, idx_t, idx*choiceZ);
     
               idx++;   
    }
   
   size_t new_choiceZ=idx_t2-(idx-1)*choiceZ;
    if (new_choiceZ!=0){ 
        sum1+=callMethodZN(FuncForZN1,new_choiceZ, idx_t,idx_t2);
        sum2+=callMethodZN(FuncForZN2,new_choiceZ, idx_t,idx_t2);
        
    }
  
        
        if (idx_t%2==1){
            sum2+=0.5*delta*(FuncForZN2(choiceZ,idx_t,idx_t2));
            sum1+=0.5*delta*(FuncForZN1(choiceZ,idx_t,idx_t2));
          
        }
      sum1*=cos(w[idx]);
      sum2*=sin(w[idx]);
      //  sum=-0.5*exp(-t)*(sum1+sum2);
      //  sum=-0.5*(sum1+sum2)+r[idx_t];
        sum=(sum1+sum2);
return sum;
}

double FuncForTW ( size_t choiceZ, size_t idx_t)
{
    double sum, sum1=0, sum2=0;
    double delta2=delta*delta;
    double t=(double)idx_t*delta;
    
     size_t idx_t2=floor(idx_t/2);

    size_t idx=1;
     
    while (idx*choiceZ<=idx_t2) {
        
        sum1+=callMethodZN(FuncForZN1,choiceZ, idx_t, idx*choiceZ);
        sum2+=callMethodZN(FuncForZN2,choiceZ, idx_t, idx*choiceZ);
     
               ++idx;   
    }
   
   size_t new_choiceZ=idx_t2-(idx-1)*choiceZ;
    if (new_choiceZ!=0){ 
        sum1+=callMethodZN(FuncForZN1,new_choiceZ, idx_t,idx_t2);
        sum2+=callMethodZN(FuncForZN2,new_choiceZ, idx_t,idx_t2);
        
    }
  
        
        if (idx_t%2==1){
            sum2+=0.5*delta*(FuncForZN2(choiceZ,idx_t,idx_t2));
            sum1+=0.5*delta*(FuncForZN1(choiceZ,idx_t,idx_t2));
          
        }
      sum1*=-sin(w[idx]);
      sum2*=cos(w[idx]);
      //  sum=-0.5*exp(-t)*(sum1+sum2);
        sum=(sum1+sum2)/r[idx];
        
return sum;
}


double callMethodXN(double (*f)( size_t ,size_t , size_t ), 
                                  size_t idx_t,size_t idx_z, size_t idx_x,size_t choiceX)
{
   if ((choiceX)==0 ){
    return 0.0;   
   }
   
   if ((choiceX)==1 ){
     return  Trapezoid_ruleN(f, idx_t, idx_z,  idx_x);
      
   }
   
   if ((choiceX)==2) {
     return  Simpson_ruleN(f, idx_t, idx_z,  idx_x);
      
   }
   
   if ((choiceX)==3) {
     return  Simpson_3_8N(f, idx_t, idx_z,  idx_x);
      
   }
   
    if (choiceX==4){
     return  BooleN(f, idx_t, idx_z,  idx_x);
      
   }
   if (choiceX==5){
     return  0.5*delta*(f(idx_t,idx_z,idx_x-1)+f(idx_t,idx_z,idx_x))+ 
             
             BooleN(f, idx_t-1, idx_z,  idx_x-1);
      
   }
   if (choiceX==6) {
     return  Simpson_narrowN(f, idx_t, idx_z,  idx_x);
      
   }
   if (choiceX==7){
     return  0.5*delta*(f(idx_t,idx_z,idx_x-1)+f(idx_t,idx_z,idx_x))+ 
             
             Simpson_narrowN(f, idx_t-1, idx_z,  idx_x-1);
      
   }
   if ((choiceT==8)||(choiceT==9)) {
     return  Simpson_extendedN(f, idx_t, idx_z,  idx_x);
      
   }
   return 0;
}

double callMethodZN(double (*f)(size_t ,size_t , size_t ), 
                            size_t choiceZ, size_t idx_t,size_t idx_z)
{
   if ((choiceZ)==0 ){
    return 0.0;   
   }
   
   if ((choiceZ)==1 ){
     return  Trapezoid_ruleN(f ,choiceX,idx_t, idx_z);
      
   }
   
   if ((choiceZ)==2) {
     return  Simpson_ruleN(f ,choiceX, idx_t, idx_z);
      
   }
   
   if ((choiceZ)==3) {
     return  Simpson_3_8N(f ,choiceX, idx_t, idx_z);
   }
    if (choiceZ==4){
     return  BooleN(f ,choiceX, idx_t, idx_z);
      
   }
   if (choiceZ==5){
     return  0.5*delta*(f(choiceX,idx_t,idx_z-1)+f(choiceX,idx_t,idx_z))+ 
            
             BooleN(f ,choiceX, idx_t, idx_z-1);
      
   }
   if (choiceZ==6) {
     return  Simpson_narrowN (f ,choiceX, idx_t, idx_z);
      
   }
   if (choiceZ==7){
     return  0.5*delta*(f(choiceX,idx_t,idx_z-1)+f(choiceX,idx_t,idx_z))+ 
             
             Simpson_narrowN(f ,choiceX, idx_t, idx_z-1);
      
   }
   if ((choiceT==8)||(choiceT==9)) {
     return  Simpson_extendedN(f ,choiceX, idx_t, idx_z);
      
   }
   return 0;
}

double callMethodTN(double (*f)(size_t, size_t ), 
                            size_t choiceT, size_t idx_t,char flag)
{
   if ((choiceT)==0 ){
    return 0.0;   
   }
   
   if ((choiceT)==1 ){
     return  Trapezoid_ruleTN(f,choiceZ, idx_t,flag);
      
   }
   
   if ((choiceT)==2) {
       
     return  Simpson_ruleTN(f, choiceZ, idx_t);
     
      
   }
   
   if ((choiceT)==3) {
     return  Simpson_3_8TN(f,choiceZ, idx_t);
      
   }
    if (choiceT==4){
     return  BooleTN(f,choiceZ, idx_t);
      
   }
   if (choiceT==5){
     return  BooleTN(f,choiceZ, idx_t-2)+
             delta*(f(choiceZ,idx_t-2)+f(choiceZ,idx_t));
      
   }
   
   if (choiceT==6) {
     return  Simpson_narrowTN(f,choiceZ, idx_t);
      
   }
   if (choiceT==7){
     return  Simpson_narrowTN(f,choiceZ, idx_t-2)+
             delta*(f(choiceZ,idx_t-2)+f(choiceZ,idx_t));
      
   }
   if ((choiceT==8)||(choiceT==9)) {
     return  Simpson_extendedTN(f,choiceZ, idx_t);
      
   }
   
   
   return 0;
}