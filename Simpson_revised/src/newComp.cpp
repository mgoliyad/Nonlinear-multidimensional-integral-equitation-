#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <fcntl.h>
//#include <unistd.h>


#include "newComp.h"
#include "Romb.h"
#include "algorithm.h"
#include <complex>     

using namespace std;
 
    
complex<double> FuncForX (int idx_t,int idx_z, int idx_x)
{

    double delta2=delta*delta;
    double delta3=delta2*delta;
    
        double z=delta*(double)idx_z;
        double x=delta*(double)idx_x;
        double z2=delta2*(double)idx_z*(double)idx_z;

        complex<double> arg1=ialpha1_2*z*(z+x)-nu1_3*z2*(2.*z/3.+x)-beta1*(2.*z+x);
        complex<double> arg2=ialpha2_2*z*(z+x)-nu2_3*z2*(2.*z/3.+x)-beta2*(2.*z+x);
        complex<double> A_product=(A[idx_t-idx_z-idx_x]*conj(A[idx_t-2*idx_z-idx_x]));           
        complex<double> sum=c1*exp(arg1)+c2*exp(arg2);
        sum*=A_product;
        
     //   printf("FinalFuncForX sum=%f, z=%int\n, idx=%int, A_product=%f\n",real(sum),idx_z, idx_t-idx_z, real(A_product));
return sum;
}

complex<double> FuncForZ (int choiceX, int idx_t,int idx_z )
{
  //  printf("FuncForZ for z=%int, t=%int\n",idx_z,idx_t);
    complex<double> sum;
    double delta2=delta*delta;
    double z2=delta2*(double)idx_z*(double)idx_z;
 
    int idx=1;
        
    while (idx*choiceX<=idx_t-idx_z-idx_z) {
        sum+=callMethodX(FuncForX,idx_t, idx_z,idx*choiceX,choiceX); 
       // printf("FuncForZ sum=%f, i=%int\n",real(sum),idx);
        ++idx;       
    }
    
    int new_choiceX=idx_t-idx_z-idx_z-(idx-1)*choiceX;
    if (new_choiceX!=0){ 
       // printf("ExFuncForZ sum=%f, i=%int\n",real(sum),idx);
        sum+=callMethodX(FuncForX,idx_t, idx_z, idx_t-idx_z-idx_z, new_choiceX);
    }

        sum*=z2*A[idx_t-idx_z];
//printf("FinalFuncForZ sum=%f, z=%int\n, t=%int, idx=%int, A=%f\n",real(sum),idx_z,idx_t, idx_t-idx_z, real(A[idx_t-idx_z]));
return sum;
}

complex<double> FuncForT (  int choiceZ, int idx_t)
{
    complex<double> sum=0.0;
    double delta2=delta*delta;
   
     int idx_t2=floor(idx_t/2);
     int idx=1;
    
    while (idx*choiceZ<=idx_t2) {
           sum+=callMethodZ(FuncForZ,choiceZ, idx_t, idx*choiceZ);
           //    printf("FuncForT sum=%f, i=%int, choiceZ=%int, idx_t2=%int\n",real(sum),idx,choiceZ,idx_t2);
               ++idx; 
    }
   
   int new_choiceZ=idx_t2-(idx-1)*choiceZ;
    if (new_choiceZ!=0){ 
        sum+=callMethodZ(FuncForZ,new_choiceZ, idx_t,idx_t2);
       //  printf("Ex1FuncForT sum=%f, i=%int, new_choiceZ=%int\n",real(sum),idx, new_choiceZ);
    }
  
        
        if (idx_t%2==1){
            sum+=0.5*delta*(FuncForZ(choiceZ,idx_t,idx_t2));
          //   printf("Ex2FuncForT sum=%f, i=%int\n",real(sum),idx);
         
        }
  // printf("F_T:t=%int %f\n",idx_t,real(sum));
   return sum;
}


complex<double> callMethodX(complex<double> (*f)(int ,int , int ), 
                                  int idx_t,int idx_z, int idx_x,int choiceX)
{
   if ((choiceX)==0 ){
    return 0.0;   
   }
   
   if ((choiceX)==1 ){
     return  Trapezoid_rule(f, idx_t, idx_z,  idx_x);
      
   }
   
   if ((choiceX)==2) {
     return  Simpson_rule(f, idx_t, idx_z,  idx_x);
      
   }
   
   if ((choiceX)==3) {
     return  Simpson_3_8(f, idx_t, idx_z,  idx_x);
      
   }
   
    if (choiceX==4){
     return  Boole(f, idx_t, idx_z,  idx_x);
      
   }
   if (choiceX==5){
     return  0.5*delta*(f(idx_t,idx_z,idx_x-1)+f(idx_t,idx_z,idx_x))+ 
             
             Boole(f, idx_t-1, idx_z,  idx_x-1);
      
   }
   if (choiceX==6) {
     return  Simpson_narrow(f, idx_t, idx_z,  idx_x);
      
   }
   if (choiceX==7){
     return  0.5*delta*(f(idx_t,idx_z,idx_x-1)+f(idx_t,idx_z,idx_x))+ 
             
             Simpson_narrow(f, idx_t-1, idx_z,  idx_x-1);
      
   }
   if (choiceX==8) {
     return  Simpson_extended(f, idx_t, idx_z,  idx_x);
   }
   /*if (choiceX==9) {
     return  romberg(f,(idx_x-1)*delta, (idx_x)*delta, idx_z*delta,idx_t*delta);
   }*/
   return 0;
}

complex<double> callMethodZ(complex<double> (*f)(int ,int , int ), 
                            int choiceZ, int idx_t,int idx_z)
{
   if ((choiceZ)==0 ){
     return  0.0;
    // return NoInt(f ,choiceX,idx_t, idx_z);    
   }
   
   if ((choiceZ)==1 ){
       complex<double> sum=Trapezoid_rule(f ,choiceX,idx_t, idx_z);
    //   printf("t:=%int; z:=%int;\n",idx_t,idx_z );
    //   printf("IntX_r:=%f; IntX_w:=%f\n",abs(sum), arg(sum));
     return  Trapezoid_rule(f ,choiceX,idx_t, idx_z);
      
   }
   
   if ((choiceZ)==2) {
     return  Simpson_rule(f ,choiceX, idx_t, idx_z);
      
   }
   
   if ((choiceZ)==3) {
     return  Simpson_3_8(f ,choiceX, idx_t, idx_z);
   }
    if (choiceZ==4){
     return  Boole(f ,choiceX, idx_t, idx_z);
      
   }
   if (choiceZ==5){
     return  0.5*delta*(f(choiceX,idx_t,idx_z-1)+f(choiceX,idx_t,idx_z))+ 
            
             Boole(f ,choiceX, idx_t, idx_z-1);
      
   }
   if (choiceZ==6) {
     return  Simpson_narrow (f ,choiceX, idx_t, idx_z);
      
   }
   if (choiceZ==7){
     return  0.5*delta*(f(choiceX,idx_t,idx_z-1)+f(choiceX,idx_t,idx_z))+ 
             
             Simpson_narrow(f ,choiceX, idx_t, idx_z-1);
      
   }
   if (choiceZ==8) {
     return  Simpson_extended(f ,choiceX, idx_t, idx_z);
      
   }
   /*if (choiceZ==9) {
     return  
    sum+=romberg(f,(idx_z-1)*delta, idx_z*delta, idx_t*delta);
     romberg( complex<double> (*f)(double , double  ), double a, double b, double t ){
   }*/
   return 0;
}

//complex<double> callMethodT(complex<double> (*f)(complex<double> *, int, int ), complex<double> *A,int choiceT, int idx_t)
complex<double> callMethodT(complex<double> (*f)(int, int ), int choiceT, int idx_t)
{
  
   if ((choiceT)==0 ){
     return NoIntT(f,choiceZ, idx_t);
   }
   
   if ((choiceT)==1 ){
     return  Trapezoid_ruleT(f,choiceZ, idx_t);
      
   }
   
   if ((choiceT)==2) {
     return  Simpson_ruleT(f,choiceZ, idx_t);
       
   }
   
   if ((choiceT)==3) {
     return  Simpson_3_8T(f,choiceZ, idx_t);
      
   }
    if (choiceT==4){
     return  BooleT(f,choiceZ, idx_t);
      
   }
   
   if (choiceT==6) {
     return  Simpson_narrowT(f,choiceZ, idx_t);
      
   }
  
   if (choiceT==8) {
     return  Simpson_extendedT(f,choiceZ, idx_t);
      
   }
   return 0;
}

complex<double> callMethodR(complex<double> (*f)(int,int ), int choiceT, int idx_t)
{
  
   if ((choiceT)==0 ){
    return NoIntR(f, idx_t);   
   }
   
   if ((choiceT)==1 ){
     return  Trapezoid_ruleR(f, idx_t);
      
   }
   
   if ((choiceT)==2) {
      
     return  Simpson_ruleR(f, idx_t);
      
   }
   
   if ((choiceT)==3) {
     return  Simpson_3_8R(f, idx_t);
      
   }
    if (choiceT==4){
     return  BooleR(f, idx_t);
      
   }
 
   if (choiceT==6) {
     return  Simpson_narrowR(f, idx_t);
      
   }
  
   if (choiceT==8) {
     return  Simpson_extendedR(f, idx_t);
      
   }
   
   return 0;
}


void choseMethod(int idx,int myChoice, complex<double> * A_new)
{
    
    if (choiceT!=9){
        
                        A_new[idx]=callMethodT(FuncForT, myChoice, idx);
                        
    
                        
                    }else{
        
                        A_new[idx]=callMethodR(R_integral, myChoice, idx);
                       // printf("Multiple choice method:r[%int]=%f,w[%int]=%f\n",idx, abs(A[idx]),idx, arg(A[idx]));
                       // fflush(stdout);
                       // A[idx]=Trapezoid_ruleR_new(R_integral_new, idx);
                           
                    }
    
}