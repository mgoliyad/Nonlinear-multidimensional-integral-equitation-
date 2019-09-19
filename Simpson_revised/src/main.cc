#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <fcntl.h>
#include <unistd.h>
#include <fstream> 
//#include <complex.h>
#include "newComp.h"
#include "newNonComp.h"
#include "NonComp.h"
#include "Romb.h"
#include "algorithm.h"
#include "Comp.h"


#include <complex>       
using namespace std;

using std::conj;

using namespace std;

using std::conj;
double c1;
double c2;
double nu1;
double alpha1;
double nu2;
int Interval_t;
int Interval_t2;

extern const complex<double> _Complex_I(0.0,1.0); 
extern const double alpha2=0.0;
extern const complex <double> ialpha2_2=0.0;
complex<double> ialpha1_2;
extern const double beta1=0.0;
extern const double beta2=0.0;
extern const int num=10;
extern const double delta=0.05;
//const double error_r=0.000001;
//const double error_w=0.00000001;

 double nu1_3;
 double nu2_3;

 complex<double> *A;

 int max_stepsX=1;
 int max_stepsZ=1;
 int max_stepsT=1;

 int choiceT;
 int choiceZ;
 int choiceX;
 int choiceA;


int main(int argc, char** argv)
{
    printf ("Enter your preferred calculation method for X_Integral: 1 - Trapezoid, 2 - Simpson, 3 - Simpson_3_8; \n 4 - Boole, "
                                                                    "6 - Simpson_narrow, 8 - Simpson_extended; 9 - Romberg;\n");
    fflush(stdout);
    scanf ("%d", &choiceX);
        
    if (choiceX==9){
        printf ("Enter number of X-Iteration for Romberg method\n");
        fflush(stdout);
        scanf ("%d", &max_stepsX);
        }
        
    printf ("Enter your preferred method for Z-Integral: 1 - Trapezoid, 2 - Simpson, 3 - Simpson_3_8; 4 - Boole, 6 - Simpson_narrow, 8 - Simpson_extended; 9-Romberg\n");
    fflush(stdout);
    scanf ("%d", &choiceZ);
        
    if (choiceZ==9){
        printf ("Enter number of Z-Iteration for Romberg method\n");
        fflush(stdout);
        scanf ("%d", &max_stepsZ);
        }
        
    printf ("Enter your preferred method for T-Integral: 0 - No Integration;  1 - Trapezoid, 2 - Simpson, 3 - Simpson_3_8; 4 - Boole, 6 - Simpson_narrow, 8 - Simpson_extended; 9-Romberg\n");
    fflush(stdout);
    scanf ("%d", &choiceT);
        
    if (choiceT==9){
        printf ("Enter number of T-Iteration for Romberg method\n");
        fflush(stdout);
        scanf ("%d", &max_stepsT);
        }

    if (choiceT==9){
        printf ("Enter your preferred method for A-Integral: 0 - No Integration; 1 - Trapezoid, 2 - Simpson, 3 - Simpson_3_8; 4 - Boole, 6 - Simpson_narrow, 8 - Simpson_extended;\n");
        fflush(stdout);
        scanf ("%d", &choiceA);
        }
    
    printf ("Enter c1:\n");
    fflush(stdout);
    scanf ("%lf", &c1);
    printf ("Enter nu1:\n");
    fflush(stdout);
    scanf ("%lf", &nu1);
    printf ("Enter alpha1:\n");
    fflush(stdout);
    scanf ("%lf", &alpha1);
    printf ("Enter Interval:\n");
    fflush(stdout);
    scanf ("%d", &Interval_t);
        
    
    ialpha1_2=_Complex_I*alpha1*alpha1;
    nu1_3=nu1*nu1*nu1;
    nu2_3=10.*nu1_3;
    c2=1.0-c1;
         
    int idx=0, idx_delta=0;

    double  tt=0.0;
    double deltaS=delta/num;
    int num_t=Interval_t/delta;
    
    A=new complex<double> [num_t+1];
    complex<double> * A_D=new complex<double> [num+1];
        
    A[0]=1.0;
    A_D[0]=1.0;
     
    simple_integration (A_D,1,deltaS);
  
    for(int idx_delta=2;idx_delta<=num;idx_delta++){
        Simpson_integrationT (A_D,idx_delta,deltaS);
        }

    A[1]=A_D[num];
       

        FILE *OUTFILE0;
	char Filename[sizeof "/Users/masha/Results/ComplexOut_c1_0.00_c2_0.00_nu1_0.0_al1_0.0.txt"];
	sprintf(Filename,"/Users/masha/Results/ComplexOut_c1_%.1f_c2_%.1f_nu1_%.1f_al1_%.1f.txt",c1,c2,nu1,alpha1);
	OUTFILE0=fopen(Filename, "w");

        FILE* gnuplot_pipeA0 = popen ("gnuplot -persistent", "w");
        fprintf(gnuplot_pipeA0, "set title '%s for C1=%f; Nu1=%f; Alpha1=%f;'\n", "A0 modul plot",c1,nu1,alpha1);
        fprintf(gnuplot_pipeA0, "plot '-'\n");
        
        fprintf(gnuplot_pipeA0, "%f %f\n", tt, abs(A[idx]));
        fprintf(OUTFILE0,"t[%2d]:=%f; A1[%2d]:=%f; A2[%2d]:=%f;\n",idx, tt, idx,abs(A[idx]), idx,arg(A[idx]));
        
        
        FILE* gnuplot_pipeARG = popen ("gnuplot -persistent", "w");
        fprintf(gnuplot_pipeARG, "set title '%s for C1=%f; Nu1=%f; Alpha1=%f;'\n", "A2 argument plot",c1,nu1,alpha1);
        fprintf(gnuplot_pipeARG, "plot '-'\n");
        fprintf(gnuplot_pipeARG, "%f %f\n", tt, arg(A[idx]));
        
        
        
        idx=1;
        tt+=delta;

        printf("Multiple choice method: r[%d]=%f,w[%d]=%f\n",idx, abs(A[idx]),idx, arg(A[idx]));
        fflush(stdout);
        
        int myChoice;
        
        ++idx;
        bool checkAgain=false;
        complex<double> temp=0;
        while (tt<Interval_t){
            tt=idx*delta;
            printf("*************  t=%f; idx=%d ********************\n",tt,idx);
            fflush(stdout);
            if (choiceT!=9){
                myChoice= (choiceT < idx) ? choiceT : idx;
                }else{
                    myChoice= (choiceA < idx) ? choiceA : idx;
                    }
                
            if ((myChoice==5)|| (myChoice==7)){
                --myChoice;
                }
            
            choseMethod(idx, myChoice, A);
                        
            fprintf(gnuplot_pipeA0, "%f %f\n", tt, abs(A[idx]));
            fprintf(gnuplot_pipeARG, "%f %f\n", tt, arg(A[idx]));
            fprintf(OUTFILE0,"t[%2d]:=%f; A1[%2d]:=%8f; A2[%2d]:=%8f;\n",idx, tt, idx,abs(A[idx]), idx,arg(A[idx]));
            
            printf("r[%d]=%f,w[%d]=%f\n",idx, abs(A[idx]),idx, arg(A[idx]));
            fflush(stdout);
            
            if ((isnan(abs(A[idx])))==1 || isnan(arg(A[idx]))==1 ){
		fclose(OUTFILE0);
		printf("Unstable\n");
		exit(1);
            }
                        
            idx++;
        }
        
        fprintf(gnuplot_pipeA0, "e\n");
        fprintf(gnuplot_pipeA0, "refresh\n");
        pclose (gnuplot_pipeA0);
        
        fprintf(gnuplot_pipeARG, "e\n");
        fprintf(gnuplot_pipeARG, "refresh\n");
        pclose (gnuplot_pipeARG);
        
        fclose(OUTFILE0);
        free(A);
        free(A_D);
	exit(2);
}
