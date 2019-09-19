#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <complex.h>
#include "Comp.h"
#include "NonComp.h"

#include <complex>       
using namespace std;

double phi(size_t idx_x, size_t idx_z, double nu,double delta)
{
    double result;
    double nu3=nu*nu*nu;
    result=(2.0*(double)idx_z/3.0+(double)idx_x)*delta;
    result*=(double)idx_z*(double)idx_z*nu3*delta*delta;
   return result;
}

double psi(size_t idx_x, size_t idx_z,double alpha,double delta)
{double  result;
double alpha2=alpha*alpha;
    result=alpha2*(double)idx_z*delta*((double)idx_x+(double)idx_z)*delta;
return result;
}

double  Omega (double  *w, size_t idx_x, size_t idx_z,size_t idx_t,double alpha,double delta)
{
   double result;
   result= w[idx_t-idx_z]+w[ idx_t-idx_z-idx_x]-w[ idx_t-2*idx_z-idx_x];
   result+=psi(idx_x,idx_z,alpha,delta);
   return result;
}

double R(double  *r, size_t idx_x, size_t idx_z, size_t idx_t,double delta)
{double  result;
size_t i1=idx_t-idx_z;
size_t i2=idx_t-idx_z-idx_x;
size_t i3=idx_t-2*idx_z-idx_x;

double r1=r[i1];
double r2=r[i2];
double r3=r[i3];
if (idx_t<10){
//printf("r1=%f, r2=%f,r3=%f,i1=%zu,i2=%zu,i3=%zu,idx_x=%zu,idx_z=%zu,idx_t=%zu\n",r1,r2,r3,i1,i2,i3,idx_x,idx_z,idx_t);
}
    result=r[idx_t-idx_z]*r[idx_t-idx_z-idx_x]*r[idx_t-2*idx_z-idx_x];
return result;
}
/************************************************/
double FuncForX_N (double  *r, double  *w, size_t idx_x,size_t idx_z, size_t idx_t,size_t trig,double delta)
{
    double  result=0.0;
    
    
    if (trig==1){
        result=exp(-phi(idx_x,idx_z,nu1,delta)-beta1*(double)(2*idx_z+idx_x)*delta);
        result*=R(r, idx_x,idx_z,idx_t,delta);
        result*=cos(Omega(w,idx_x,idx_z,idx_t,alpha1,delta));
    }
    if (trig==2){
        result=exp(-phi(idx_x,idx_z,nu2,delta)-beta2*(double)(2*idx_z+idx_x)*delta);
        result*=R(r, idx_x,idx_z,idx_t,delta);
        result*=cos(Omega(w,idx_x,idx_z,idx_t,alpha2,delta));
    }
    if (trig==3){
        result=exp(-phi(idx_x,idx_z,nu1,delta)-beta1*(double)(2*idx_z+idx_x)*delta);
        result*=R(r, idx_x,idx_z,idx_t,delta);
        result*=sin(Omega(w,idx_x,idx_z,idx_t,alpha1,delta));
    }
    if (trig==4){
        result=exp(-phi(idx_x,idx_z,nu2,delta)-beta2*(double)(2*idx_z+idx_x)*delta);
        result*=R(r, idx_x,idx_z,idx_t,delta);
        result*=sin(Omega(w,idx_x,idx_z,idx_t,alpha2,delta));
    }
return result;
}

/************************************************/
double  X_integral_N (double *r, double *w, size_t idx_z, size_t idx_t,size_t trig,double delta)
{
    double sum_x=0.0;
    size_t idx_x;
    if((idx_t-2*idx_z)!=0)
    {
         sum_x=0.5*(FuncForX_N (r,w,0,idx_z,idx_t,trig,delta) + 
                 FuncForX_N(r,w,(idx_t-2*idx_z),idx_z,idx_t,trig,delta));

		for(idx_x=1;idx_x<(idx_t-2*idx_z);idx_x++)
		{
                    sum_x+=FuncForX_N (r,w,idx_x,idx_z,idx_t,trig,delta);
		}
    }
return sum_x*delta;
}
/************************************************/
double  SimpsonX_N(double *r, double *w, size_t idx_z, size_t idx_t,size_t trig,double delta)

{
    size_t idx_x,x_Max;
    double sum_z=0.0,sum_last=0.0;
    double delta_3=delta/3.0;

    if ((idx_t-2*idx_z)%2==0)
    {
        x_Max=idx_t-2*idx_z;
    }
    else{
        x_Max=idx_t-2*idx_z-1;
        sum_last=0.5*(FuncForX_N(r,w,x_Max+1,idx_z,idx_t,trig,delta)+
                            FuncForX_N(r,w,x_Max,idx_z, idx_t,trig,delta));
    }

    sum_z= FuncForX_N(r,w,0,idx_z,idx_t,trig,delta)+
            FuncForX_N(r,w,x_Max, idx_z,idx_t,trig,delta);

    for (idx_x=1;idx_x<x_Max;idx_x++)
    {
        if ((idx_x)%2==1){
            sum_z+= 4.0*FuncForX_N(r,w,idx_x,idx_z, idx_t,trig,delta);
        }
        else{
            sum_z+= 2.0*FuncForX_N(r,w,idx_x,idx_z, idx_t,trig,delta);
        }
    }

    sum_z=delta*(0.3*sum_z+sum_last);
return sum_z;
}

double  SimpsonX_N38 (double *r, double *w, size_t idx_z, size_t idx_t,size_t trig,double delta)

{
    size_t idx_x,x_Max;
    double sum_z=0.0,sum_last=0.0;
    double delta_3_8=3.*delta/8.;

    
    x_Max=idx_t-2*idx_z;
    

    sum_z= FuncForX_N(r,w,0,idx_z,idx_t,trig,delta)+
            FuncForX_N(r,w,x_Max, idx_z,idx_t,trig,delta);

    for (idx_x=1;idx_x<x_Max;idx_x++)
    {
        if (idx_x%3!=0){
            sum_z+= 3.0*FuncForX_N(r,w,idx_x,idx_z, idx_t,trig,delta);
        }
        else{
            
            sum_z+= 2.0*FuncForX_N(r,w,idx_x,idx_z, idx_t,trig,delta);
        }
    }
    sum_z=sum_z*delta_3_8;
return sum_z;
}
/************************************************/

double callSimpsonX_N(double *r, double *w, size_t idx_z, size_t idx_t,size_t trig,double delta)
{
    if (idx_t-2*idx_z<2){
        return X_integral_N(r,w,idx_z, idx_t, trig,delta);
    }
    else{
       // if ((idx_t-2*idx_z)%3!=0){
        return SimpsonX_N(r,w,idx_z, idx_t,trig,delta);
        //}
       // else{
       //     return SimpsonX_N38(r,w,idx_z, idx_t,trig,delta);
       // }
    }
}
/************************************************/
double FuncForZ_N (double  *r, double *w, size_t idx_z, size_t idx_t,size_t trig,double delta)
{
    double  z2, result;
    
    z2=delta*delta*(double)idx_z*(double)idx_z;
    result=z2*callSimpsonX_N(r,w, idx_z, idx_t,trig,delta);
return result;
}
/************************************************/
double Z_integral_N (double *r, double *w, size_t idx_t,size_t trig,double delta)
{
    double sum_z=0.0;
    size_t idx_z, idx_t2;

    idx_t2=idx_t>>1;
     if(idx_t2!=0){
            sum_z=0.5*( FuncForZ_N(r,w,idx_t2, idx_t,trig,delta) );

            for(idx_z=1;idx_z<idx_t2;idx_z++){
                sum_z+=( FuncForZ_N(r,w,idx_z,idx_t,trig,delta) );
            }

        }
    if (idx_t%2==1){
        sum_z+=0.5*FuncForZ_N(r,w,idx_t2,idx_t,trig,delta);
    }
    sum_z*=delta;
return sum_z;
}

/**************************************************************************************************/
double   SimpsonZ_N(double   *r, double *w, size_t idx_t,size_t trig,double delta)
{
    size_t idx_z, z_Max=0;
    double sum_z=0.0,sum_last=0.0;
   // double delta_3=delta/3.0;

    size_t idx_t2=floor(idx_t/2);

    if ((idx_t2)%2==0){
        z_Max=idx_t2;
    }
    else{
        z_Max=idx_t2-1;
        sum_last=0.5*delta*( FuncForZ_N(r,w,idx_t2, idx_t,trig,delta)+ 
                FuncForZ_N(r,w,idx_t2-1, idx_t,trig,delta));
    }
    sum_z= (FuncForZ_N(r,w,z_Max,idx_t,trig,delta));

    for (idx_z=1;idx_z<z_Max;idx_z++)
    {
        if ((idx_z)%2==1){
            sum_z+= 4.*FuncForZ_N(r,w,idx_z, idx_t,trig,delta);
        }
        else{
            sum_z+= 2.*FuncForZ_N(r,w,idx_z, idx_t,trig,delta);
        }
    }

    sum_z=(0.3*sum_z + sum_last)*delta;

return sum_z;
}

double   SimpsonZ_N38(double   *r, double *w, size_t idx_t,size_t trig,double delta)
{
    size_t idx_z, z_Max=0;
    double sum_z=0.0,sum_last=0.0;
    double delta_3_8=3.*delta/8.;
    size_t idx_t2=idx_t>>1;
    
    z_Max=idx_t2;
    
    sum_z= FuncForZ_N(r,w,0,idx_t,trig,delta)+
            FuncForZ_N(r,w,z_Max,idx_t,trig,delta);

    for (idx_z=1;idx_z<z_Max;idx_z++)
    {
        if (idx_z%3!=0){
            sum_z+= 3.0*FuncForZ_N(r,w,idx_z, idx_t,trig,delta);
        }
        else{
            
            sum_z+= 2.0*FuncForZ_N(r,w,idx_z, idx_t,trig,delta);
        }
    }
    
    sum_z=sum_z*delta_3_8;
return sum_z;
}
/**************************************************************************************************/


double callSimpsonZ_N (double *r, double *w, size_t idx_t,size_t trig,double delta)
{
    if (idx_t<4){
        return Z_integral_N(r,w,idx_t,trig,delta);
    }
    else{
       // if (idx_t%3!=0)
        return SimpsonZ_N(r,w, idx_t,trig,delta);
       // else{
       //     return SimpsonZ_N38(r,w, idx_t,trig,delta);
      //  }
    }
}
/**************************************************************************************************/
double FuncForUN(double *r, double *w, size_t idx,double delta)
{
double q1, q2,result;

    q1=-0.5*(c1*callSimpsonZ_N (r,w,idx,1,delta)-
        c2*callSimpsonZ_N (r,w,idx,2,delta));
    q2=-0.5*(c1*callSimpsonZ_N (r,w,idx,3,delta)-
        c2*callSimpsonZ_N (r,w,idx,4,delta));
    result=(q1*cos(w[idx])+q2*sin(w[idx]))*exp(-(double)idx*delta);
    //printf("FUN3: result=%f\n",result); 
    return result;
}
double FuncForWN(double *r, double *w, size_t idx,double delta)
{
    double q1, q2,result;
    q1=-0.5*(c1*callSimpsonZ_N (r,w,idx,1,delta)-
        c2*callSimpsonZ_N (r,w,idx,2,delta));
    q2=-0.5*(c1*callSimpsonZ_N (r,w,idx,3,delta)-
        c2*callSimpsonZ_N (r,w,idx,4,delta));
    result=-q1*sin(w[idx])+q2*cos(w[idx]);
    result/=r[idx];
    return result;
}
double  Simpson_integrationUN(double *r, double *w, size_t idx,double delta)
{
double  sum;
double delta_3=delta/3.;
double delta_3_8=3.*delta/8.;
double help1, help2, help3;
//double delta_24=delta/24.;
//printf("UN1: idx=%zu\n",idx);
//if (idx<3){
        sum= r[idx-2]*exp(-(double)(idx-2)*delta);
       /* printf("UN:sum=%f\n",sum);
        help1=FuncForUN(r,w,idx,delta);
        help2=FuncForUN(r,w,idx-1,delta);
        help3=FuncForUN(r,w,idx-2,delta);
        printf("WU:help1=%f;help2=%f;help3=%f\n",help1,help2,help3);*/
//printf("UN2: sum=%f\n",sum); 
        sum+=(delta_3)*(  FuncForUN(r,w,idx,delta) +  
                    4.*(  FuncForUN(r,w,idx-1,delta  ))+
                            FuncForUN(r,w,idx-2,delta )); 
        // printf("UN:sum=%f\n",sum);
  /*}
  //  if ((idx>=3)||(idx<6))
 if (idx>=3)
        {
            sum= r[idx-3]*exp(-(double)(idx-3)*delta);
            sum+=(delta_3_8)*(  
                    FuncForUN(r,w,idx,delta) +  
                    3.*(  
                    FuncForUN(r,w,idx-1,delta) +FuncForUN(r,w,idx-2,delta ))+
                            FuncForUN(r,w,idx-3,delta )); 
    }*/
        return sum;
}

double  Simpson_integrationWN(double *r, double *w, size_t idx,double delta)
{
double sum;
double delta_3=delta/3.;
double help1, help2, help3;
                    //sum= w[idx-1]+delta*FuncForWN(r,w,idx,delta);
double delta_3_8=3.*delta/8.;
                    //double delta_24=delta/24.;
// if (idx<3){
        sum= w[idx-2];
      /*  printf("WN:sum=%f\n",sum);
        help1=FuncForWN(r,w,idx,delta);
        help2=FuncForWN(r,w,idx-1,delta);
        help3=FuncForWN(r,w,idx-2,delta);
        printf("WN:help1=%f;help2=%f;help3=%f\n",help1,help2,help3);*/
        sum+=(delta_3)*(  FuncForWN(r,w,idx,delta) +  
                        4.*(  FuncForWN(r,w,idx-1 ,delta ))+
                            FuncForWN(r,w,idx-2,delta ));
      //  printf("WN:sum=%f\n",sum);
 /* }
                           //  if ((idx>=3)||(idx<6))
 if (idx>=3)
 {
        sum= w[idx-3];
        sum+=(delta_3_8)*(  FuncForWN(r,w,idx,delta) +  
         3.*(  FuncForWN(r,w,idx-1 ,delta )+ FuncForWN(r,w,idx-2 ,delta ))+
                   FuncForWN(r,w,idx-3,delta ));
  }*/

return sum;
}
/**************************************************************************************************/
void Iteration(double *r, double *w ,size_t idx,double delta)
{
    
    double  w_temp,r_temp,diff_w,diff_r;
    diff_w=1.0;
    diff_r=1.0;
    size_t count=0;
    size_t noMoreR=0;
    size_t noMoreW=0;
            
    double u;
    
  //  printf("r_comp[%zu]:=%f\n",idx,r[idx]);
  //   printf("w_comp[%zu]:=%f\n",idx, w[idx]);
     
    do{
        count++;
        
        if (noMoreR==0)
        {
            
            u=Simpson_integrationUN(r,w,idx,delta);
            //r[idx]=u*exp((double)idx*delta);
            r_temp=u*exp((double)idx*delta);
          //  printf("r[%zu]:=%f\n",idx, r_temp);
        }
    
        if (noMoreW==0)
        {
           
            //w[idx]=w[idx-1]+delta*(-q1*sin(w[idx])+q2*cos(w[idx]))/u;
            //w[idx]=Simpson_size_tegrationWN(r,w,idx,delta);
            w_temp=Simpson_integrationWN(r,w,idx,delta);
          //  printf("w[%zu]=%f\n",idx, w_temp);
        }
   
    
            diff_w=fabs(w_temp-w[idx]);
            diff_r=fabs(r_temp-r[idx]);
    
           // w_temp=w[idx];
          //  r_temp=r[idx];
            
            if ((diff_r<=error_r*fabs(r_temp))||(count>100))
            {
                noMoreR=1;
              //  printf("No more R, diff_r=%f, tempr=%f, curR=%f\n", diff_r, r_temp, r[idx]);
                r[idx]=r_temp;
                
            }
            else{
                r[idx]=r_temp;
            }
            if ((diff_w<=error_w*fabs(w_temp))||(count>100))
            {
                noMoreW=1;
              //  printf("No more w, diff_r=%f, tempw=%f, curw=%f\n", diff_r, w_temp, w[idx]);
                w[idx]=w_temp;
                
            }
            else{
                w[idx]=w_temp;
            }
          //  printf("r_comp[%zu]:=%f\n",idx,r[idx]);
  //   printf("w_comp[%zu]:=%f\n",idx, w[idx]);
    } while ((noMoreW==0)||(noMoreR==0));
    
  //  printf("exit: count=%zu, idx=%zu;r=%f,w=%f, diff_w=%f, diff_r=%f, error_r=%f, error_w=%f\n",count,idx, r[idx],w[idx],diff_w,diff_r, error_r, error_w);
    
    
    /*
    size_t num_PI=(size_t)floor(fabs(w[idx])/M_PI);
    
    if (num_PI>0){
        if (w[idx]>0)
        {
        w[idx]-=num_PI*2.*M_PI;
        }
        else
        {
        w[idx]+=num_PI*2.*M_PI;
        }
    }*/
    
}
/**************************************************************************************************/
void simple_diff_eq(double *r, double *w ,size_t idx,double delta)
{
    //printf("diff idx=%zu\n",idx);
    double q1,q2,q1_1,q2_1, w_temp,r_temp,diff_w,diff_r;
    diff_w=1.0;
    diff_r=1.0;
    size_t count=0;
    size_t noMoreR=0;
    size_t noMoreW=0;
                /* if (idx>1){
                 w[idx]=2.*w[idx-1]-w[idx-2];
                 printf("w=%f; w1=%f,w2=%f\n",w[idx],w[idx-1],w[idx-2]);
                 }
                 else{
                     w[idx]=w[idx-1];
                 }*/
    q1=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,1,delta)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,2,delta)/2.0;
    q2=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,3,delta)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,4,delta)/2.0;
    
    q1_1=-c1*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,1,delta)/2.0-
        c2*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,2,delta)/2.0;
    q2_1=-c1*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,3,delta)/2.0-
        c2*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,4,delta)/2.0;
    
    double u_1=r[idx-1]*exp(-(double)(idx-1)*delta);
    double Trap,Trap1;
    //Trap=0.5*delta*(q1*cos(w[idx])+q2*sin(w[idx])+q1_1*cos(w[idx-1])+q2_1*sin(w[idx-1]));
    
    Trap=delta*(q1*cos(w[idx])+q2*sin(w[idx]));
    
    double u=u_1+Trap;
    printf("old r=%f\n",r[idx]);
     r[idx]=u*exp((double)idx*delta);
     printf("new r=%f\n",r[idx]);
     
    // printf("exiteq1: count=%zu, idx=%zu; r=%f,w=%f, u=%f, u_1=%f, Trap=%f\n",count,idx, r[idx],w[idx],u, u_1, Trap);
    
        w_temp=w[idx];
        r_temp=r[idx];
    
    do{
        count++;
        
        q1=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,1,delta)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,2,delta)/2.0;
    q2=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,3,delta)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,4,delta)/2.0;
    
    q1_1=-c1*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,1,delta)/2.0-
        c2*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,2,delta)/2.0;
    q2_1=-c1*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,3,delta)/2.0-
        c2*exp(-(double)(idx-1)*delta)*callSimpsonZ_N (r,w,idx-1,4,delta)/2.0;
    
       // Trap=0.5*delta*(q1*cos(w[idx])+q2*sin(w[idx])+q1_1*cos(w[idx-1])+q2_1*sin(w[idx-1]));
      //  Trap1=0.5*delta*(-q1*sin(w[idx])+q2*cos(w[idx])-q1_1*sin(w[idx-1])+q2_1*cos(w[idx-1]));
        
        Trap=delta*(q1*cos(w[idx])+q2*sin(w[idx]));
        Trap1=delta*(-q1*sin(w[idx])+q2*cos(w[idx]));
        if (noMoreR==0)
        {
            u=u_1+Trap;
            r[idx]=u*exp((double)idx*delta);
            printf("newnew r=%f\n",r[idx]);
           // printf("exiteq2: count=%zu, idx=%zu; r=%f,w=%f, diff_w=%f, diff_r=%f\n",count,idx, r[idx],w[idx],diff_w,diff_r);
        }
    
        if (noMoreW==0)
        {
            w[idx]=w[idx-1]+Trap1/u;
        }
   
    
            diff_w=fabs(w_temp-w[idx]);
            diff_r=fabs(r_temp-r[idx]);
    
            w_temp=w[idx];
            r_temp=r[idx];
            
            if ((diff_r<=error_r*fabs(r_temp))||(count>100))
            {
                noMoreR=1;
            }
            if ((diff_w<=error_w*fabs(w_temp))||(count>100))
            {
                noMoreW=1;
            }
    } while ((noMoreW==0)&&(noMoreW==0));
    
    
  /*  size_t num_PI=(size_t)floor(fabs(w[idx])/M_PI);
    
    if (num_PI>0){
        if (w[idx]>0)
        {
        w[idx]-=num_PI*2.*M_PI;
        }
        else
        {
        w[idx]+=num_PI*2.*M_PI;
        }
    }
    q1=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,1)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,2)/2.0;
    q2=-c1*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,3)/2.0-
        c2*exp(-(double)idx*delta)*callSimpsonZ_N (r,w,idx,4)/2.0;
    
    double smallDelta,smallDiff;
    smallDiff=fabs(u-u_1);
    smallDelta=delta*sqrt(q1*q1+q2*q2);
    if(smallDiff>smallDelta)
    {
  //  printf("idx=%zu: smallDelta=%f; smallDiff=%f\n",idx , smallDelta,smallDiff);
    }*/
   // printf("exiteq: count=%zu, idx=%zu; r=%f,w=%f, diff_w=%f, diff_r=%f\n",count,idx, r[idx],w[idx],diff_w,diff_r);
}