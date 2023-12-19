#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#include <pthread.h>
#include <BoundaryShape.h>
#include <Variable.h>
#include <Array.h>
#include <Energy.h>
#include <InOut.h>
#include <Elastic_Field.h>
#include <Modify_NL.h>
#include <OPGL.h>

extern int N_Energy_Terms;
pthread_mutex_t **Signal_Synchronize_NL_Mutex, *M_Synchronize_NL_Mutex;
pthread_cond_t **Signal_Synchronize_NL_Cond;
extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
extern pthread_mutex_t OPGL_Mutex;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
extern OPGL_Parameters_GUI OPGL_GUI;
extern GUI_Parameter GPM;
DiscreatPoint_Run_NL_Parameter DRNP;
using namespace std;


double Minimal_Polynome(double *Coefficient, int Degree, double *Interval, double Precision)
{
    double a, b, lambda, mu, f0=0, f1=0;
    int i;
    a     =Interval[0];     b =Interval[1];
    lambda=a+0.382*(b-a);   mu=a+0.618*(b-a);

    do
    {
        for(i=0;i<=Degree;i++)
        {
            f0=f0+Coefficient[i]*pow(lambda,i);
            f1=f1+Coefficient[i]*pow(mu,i);
        }

        if(f0>f1)
        {
            a     =lambda;       b =b;
            lambda=mu;           mu=a+0.618*(b-a);

        }else
        {
            a     =a;             b =mu;
            lambda=a+0.382*(b-a); mu=lambda;
        }
    }while((b-a)>Precision);

    if(f0>f1)
    {
        return lambda;
    }else
    {
        return mu;
    }
}

double Minimal_G_L_Fun(Modify_NL_Parameters *MNP, double *Gradient_Local, double *Interval, double Precision)
{
    double a, b, lambda, mu, f0=0, f1=0, *M0;
    double ***M_Local;

    a = Interval[0];                b = Interval[1];
    lambda = a + 0.382 * (b - a);   mu = a + 0.618 * (b - a);

    M_Local=MNP->EV->M_Local;    M0 =MNP->EV->M0;

    do
    {
        M0[0] = *M_Local[1][1] + lambda * Gradient_Local[0];
        M0[1] = *M_Local[2][1] + lambda * Gradient_Local[1];
        M0[2] = *M_Local[3][1] + lambda * Gradient_Local[2];
        f0 = Energy_Local_M0(MNP->EV);

        M0[0] = *M_Local[1][1] + lambda * Gradient_Local[0];
        M0[1] = *M_Local[2][1] + lambda * Gradient_Local[1];
        M0[2] = *M_Local[3][1] + lambda * Gradient_Local[2];
        f1 = Energy_Local_M0(MNP->EV);

        //printf("%g\t%g\t%g\t%g\t%g\t%g\n",a,b,lambda,mu,f0,f1);


        if(f0>f1)
        {
            a      = lambda;       b  = b;
            lambda = mu;           mu = a+0.618*(b-a);

        }else
        {
            a     = a;             b      = mu;
            mu    = lambda;        lambda = a+0.382*(b-a);
        }
    }while((b-a)>Precision);

    if(f0>f1)
    {
        return lambda;
    }else
    {
        return mu;
    }
}

int Modify_Nonlinear_OnePoint(Modify_NL_Parameters *MNP)
{
    double Interval[2], Precision=0.01, P;
    double *Gradient_Local0, *Gradient_Local1, *s0, *s1;
    int r, s;
    //i=MNP->EV->i; j=MNP->EV->j; k=MNP->EV->k;
    Gradient_Local0 = Make1DArray(3);
    Gradient_Local1 = Make1DArray(3);
    s0              = Make1DArray(3);
    s1              = Make1DArray(3);

    //printf("%d\t%d\t%d\n",i,j,k);
    double *M0;
    double ***M_Local, W, W1, W0, M_Module;

    M_Local =MNP->EV->M_Local;    M0 =MNP->EV->M0;
    //printf("%lf\n",Mx[i+1][j][k]);

    M0[0] = *M_Local[1][1];
    M0[1] = *M_Local[2][1];
    M0[2] = *M_Local[3][1];
    W0 = Energy_Local_M0(MNP->EV);
    MNP->EV->W_Local_BC = W0;

    double alpha=0, beta=0, numerator=0, denominator=0, fraction;

    Interval[0]=-0.1; Interval[1]=0.1;

    Effective_Field_LocalCell_M0(MNP->EV,Gradient_Local0);

    for(r=0;r<3;r++)
    {
        Gradient_Local0[r] = -Gradient_Local0[r];
        s0[r] = Gradient_Local0[r];
    }

    alpha=Minimal_G_L_Fun(MNP,Gradient_Local0,Interval,Precision);
    //printf("%g\n",alpha);

    do
    {
        M0[0] = *M_Local[1][1] + alpha * Gradient_Local0[0];
        M0[1] = *M_Local[2][1] + alpha * Gradient_Local0[1];
        M0[1] = *M_Local[3][1] + alpha * Gradient_Local0[2];
        M_Module = sqrt(pow(M0[0], 2) + pow(M0[1], 2) + pow(M0[2], 2));
        if (M_Module > 4)
        {
            alpha=0.5*alpha;
        }

        if (fabs(alpha) < 1.e-5)
        {
            alpha=0;
            M0[0] = *M_Local[1][1];
            M0[1] = *M_Local[2][1];
            M0[2] = *M_Local[3][1];
        }
    } while (M_Module > 2 && alpha != 0.);
    s=0;
    do
    {
        W = Energy_Local_M0(MNP->EV);
        Effective_Field_LocalCell_M0(MNP->EV,Gradient_Local0);

        //printf("%g\t%g\t%g\n",Gradient1[0][0],Gradient1[1][0],Gradient1[2][0]);

        numerator=0;
        denominator=0;

        for(r=0;r<3;r++)
        {
            Gradient_Local0[r] = -Gradient_Local0[r];
            numerator = numerator + Gradient_Local1[r] * (Gradient_Local1[r] - Gradient_Local0[r]);
            denominator = denominator + Gradient_Local0[r] * Gradient_Local0[r];
        }

        fraction=numerator/denominator;

        if (fraction > 0 && fraction < 0.1)
        {
            beta = fraction;
        }else
        {
            beta = 0;
        }

        //printf("%g\t%g\t%g\t%g\n",W,numerator,denominator,beta);
        for(r=0;r<3;r++)
        {
            s1[r] = Gradient_Local1[r] + beta * s0[r];
        }

        alpha = Minimal_G_L_Fun(MNP, s1, Interval, Precision);

        do
        {
            M0[0] = *M_Local[1][1] + alpha * s1[0];
            M0[1] = *M_Local[2][1] + alpha * s1[1];
            M0[2] = *M_Local[3][1] + alpha * s1[2];
            M_Module = sqrt(pow(M0[0], 2) + pow(M0[1], 2) + pow(M0[2], 2));
            if(M_Module>4)
            {
                alpha=0.5*alpha;
            }

            if(fabs(alpha)<1.e-5)
            {
                alpha=0;
                M0[0] = *M_Local[1][1];
                M0[1] = *M_Local[2][1];
                M0[2] = *M_Local[3][1];
            }
        } while (M_Module > 2 && alpha != 0.);

        for (r = 0; r < 3; r++)
        {
            Gradient_Local0[r] = Gradient_Local1[r];
            s0[r]              = s1[r];
        }

        W1 = Energy_Local_M0(MNP->EV);
        //printf("%g\t%g\n",W1,alpha);

        P = fabs(W - W1) / fabs(W);

        s = s + 1;
        if (s > 10 || P <= 1.e-8)
        {
            if (W1 <= W0)
            {
                P = 0.;
            }else
            {
                P=0.;
                M0[0] = *M_Local[1][1];
                M0[1] = *M_Local[2][1];
                M0[2] = *M_Local[3][1];
            }
        }else
        {
/*
            if(W1>=W)
            {
                P=0.;
                *M_Local[1][1]=M0[0];
                *M_Local[2][1]=M0[1];
                *M_Local[3][1]=M0[2];
            }else
            {
                P=fabs(W0-W1)/fabs(W0);
            }
            P=fabs(W0-W1)/fabs(W0);
*/
        }

    }while(fabs(P)>1.e-8);

    *M_Local[1][1] = M0[0];
    *M_Local[2][1] = M0[1];
    *M_Local[3][1] = M0[2];
    //while(fabs(denominator)>1.e-6);

    //printf("%lf\t%lf\n",W0,W1);
    //printf("%lf\t%lf\t%lf\n",M_Unchanged[0],M_Unchanged[1],M_Unchanged[2]);
    //printf("%lf\t%lf\t%lf\n",*M_Local[1][1],*M_Local[2][1],*M_Local[3][1]);
    //printf("%lf\t%lf\t%lf\n",*M_Local[1][1]-M_Unchanged[0],*M_Local[2][1]-M_Unchanged[1],*M_Local[3][1]-M_Unchanged[2]);
    free(Gradient_Local0);
    free(Gradient_Local1);
    free(s0);
    free(s1);
    MNP->EV->W_Local_AC = Energy_Local_M0(MNP->EV);

    return 1;
}

int Modify_Pthread_NL(Modify_NL_Pthreads *MNPt)
{
    int Pthread_Index, NCycles, *Signal;
    int Signal_Synchronize;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int i, N, Lower, Signal_RunAndStop;

    Pthread_Index  = MNPt->Pthread_Index;
    //Pthread_Number = MNPt->Pthread_Number;
    NCycles        = MNPt->NCycles;
    Signal         = MNPt->Signal;

    CalculatePointsThreads_Lower = MNPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MNPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;

    if(CalculatePointsThreads_Lower==0)
    {
        Lower=1;
    }else
    {
        Lower=0;
    }

    pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Run;
    pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);//1 is Pthread State.
////////////////////////////////////////////////////////////////////////////////////////////////////
/*
    Energy_Variables_Distribution EVD;
    EVD.CST = MNPt->CST;    EVD.dxy = MNPt->dxy;
    EVD.h   = MNPt->h;      EVD.t   = MNPt->t;
    EVD.m   = MNPt->m;      EVD.n   = MNPt->n;
    EVD.l   = MNPt->l;
    EVD.Mx  = MNPt->Mx_Global_Pthread;
    EVD.My  = MNPt->My_Global_Pthread;
    EVD.Mz  = MNPt->Mz_Global_Pthread;
    EVD.MxV = MNPt->MxV_Global_Pthread;
    EVD.MyV = MNPt->MyV_Global_Pthread;
    EVD.MzV = MNPt->MzV_Global_Pthread;
    EVD.w   = Make3DArray(EVD.m,EVD.n,EVD.l);
    EVD.Local_Points_1D=MNPt->Local_Points_1D;
    EVD.Position_3DTo1D=MNPt->Position_3DTo1D;
*/
////////////////////////////////////////////////////////////////////////////////////////////////////
    for(N=0;N<NCycles;++N)
    {
        //printf("%d\t%d\n",Pthread_Index,N);
        pthread_mutex_lock(&M_Synchronize_NL_Mutex[Pthread_Index]);
        MNPt->PTPS->Energy = 0.;  MNPt->PTPS->SKN = 0.;

        for(i=Lower;i<CalculatePointsThreads;++i)
        {
            Modify_Nonlinear_OnePoint(&(MNPt->MNP[i]));
            if(MNPt->Continuity_Modify_Mode==1)
            {
                Continue_Modify(MNPt->MNP[i].EV);
            }
            MNPt->PTPS->Energy = MNPt->PTPS->Energy + Energy_Single(MNPt->MNP[i].EV);
            MNPt->PTPS->SKN    = MNPt->PTPS->SKN + Skyrmion_Number_Single(MNPt->MNP[i].EV);
        }

        //Modify_Nonlinear_Part(MNPt);
        pthread_mutex_unlock(&M_Synchronize_NL_Mutex[Pthread_Index]);

        pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
        Signal_Synchronize = Signal[3];
        pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

        if (Signal_Synchronize == SignalSyn_Stop)
        {
            pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            if(Signal_Synchronize == SignalSyn_PrepareToSynchronize)
            {
                pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
                Signal[1] = StatePthread_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            }

            pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_Run && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_NL_Cond[Pthread_Index][3], &Signal_Synchronize_NL_Mutex[Pthread_Index][3]);
            }
            pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_Run;
            pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            
        }
        
        /////////////////////////// Suspend thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        while (GPM.Thread_State[Pthread_Index + 1][1] == Signal_GUI_Break)//Signal_ContinueAndBreak
        {
            pthread_cond_wait(&GPM.Pthread_GUI_Cond[Pthread_Index + 1], &GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        }
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        
        /////////////////////////// Stop thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        Signal_RunAndStop = GPM.Thread_State[Pthread_Index + 1][0];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        if (Signal_RunAndStop == Signal_GUI_Stop)
        {
            break;
        }
    }

    pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Stop;
    pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[Pthread_Index][1]);//1 is Pthread State.

    //free3DArray(EVD.w,EVD.m,EVD.n);
    printf("Thread Fin\t%d\n",Pthread_Index);
    return 1;


}


void* Modify_Pthread_NL_Transfer(void *arg)
{
    Modify_NL_Pthreads *MNPt;
    MNPt=(Modify_NL_Pthreads*)arg;
    Modify_Pthread_NL(MNPt);
    return NULL;
}

int Modify_Pthread_NL_Synchronize(Modify_NL_Synchronize *MNS)
{
    int i, j, k, v, p, q, Point_Index, Boundary_Type, Synchronize_Times, N_Total;
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points;
    int Pthread_Number, ***Signal, Signal_RunAndStop, Signal_ContinueAndBreak;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST, SKN, Energy_Average;

    clock_t start_time, end_time, start_time_P, end_time_P;
    double duration;

    m=MNS->m;   n=MNS->n;   l=MNS->l;   Virtual_Points=MNS->Virtual_Points;
    Pthread_Number               = MNS->Pthread_Number;
    Bound_CalculatePoint_Threads = MNS->Bound_CalculatePoint_Threads;
    Local_Points_1D              = MNS->Local_Points_1D;
    Position_3DTo1D              = MNS->Position_3DTo1D;
    Mx_Global_Pthread            = MNS->Mx_Global_Pthread;
    My_Global_Pthread            = MNS->My_Global_Pthread;
    Mz_Global_Pthread            = MNS->Mz_Global_Pthread;
    Mx_Global_Synchronize        = MNS->Mx_Global_Synchronize;
    My_Global_Synchronize        = MNS->My_Global_Synchronize;
    Mz_Global_Synchronize        = MNS->Mz_Global_Synchronize;
    MxV_Global_Pthread           = MNS->MxV_Global_Pthread;
    MyV_Global_Pthread           = MNS->MyV_Global_Pthread;
    MzV_Global_Pthread           = MNS->MzV_Global_Pthread;
    MxV_Global_Synchronize       = MNS->MxV_Global_Synchronize;
    MyV_Global_Synchronize       = MNS->MyV_Global_Synchronize;
    MzV_Global_Synchronize       = MNS->MzV_Global_Synchronize;
    Signal                       = MNS->Signal;
    CST                          = MNS->CST;
    N_Total                      = MNS->Inner_Points + MNS->Boundary_Points;
///////////////////////////////////////////////////////////////////////////
    int zc;
    if(l==0)
    {
        zc=0;
    }else
    {
        zc=(l-1)/2;
    }
    Energy_Variables_Distribution EVD;
    EVD.CST = CST;    EVD.dxy = MNS->dxy;
    EVD.h   = MNS->h; EVD.t   = MNS->t;
    EVD.m   = m;      EVD.n   = n;
    EVD.l   = l;
    EVD.Mx    = Mx_Global_Synchronize;
    EVD.My    = My_Global_Synchronize;
    EVD.Mz    = Mz_Global_Synchronize;
    EVD.MxV   = MxV_Global_Synchronize;
    EVD.MyV   = MyV_Global_Synchronize;
    EVD.MzV   = MzV_Global_Synchronize;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D=Local_Points_1D;
    EVD.Position_3DTo1D=Position_3DTo1D;
    EVD.Energy_Coefficient = MNS->Energy_Coefficient;

    FILE *fp_Mx,*fp_My, *fp_Mz, *fp_MxV,*fp_MyV, *fp_MzV, *fp_W, *fp_SKN_D, *fp_Process;
    char filenameMx[500], filenameMy[500], filenameMz[500],filenamew[500],filenameSKN_D[500],
         filenameMxV[500],filenameMyV[500],filenameMzV[500],filenameP[500];

    sprintf(filenameMxV,"Version %d\\Process\\MxV.dat",MNS->Version);
    sprintf(filenameMyV,"Version %d\\Process\\MyV.dat",MNS->Version);
    sprintf(filenameMzV,"Version %d\\Process\\MzV.dat",MNS->Version);
    sprintf(filenameP, "Version %d\\Process.dat",MNS->Version);
    fp_MxV     = fopen(filenameMxV,"w+");
    fp_MyV     = fopen(filenameMyV,"w+");
    fp_MzV     = fopen(filenameMzV,"w+");
    fp_Process = fopen(filenameP, "w+");
    fprintf(fp_Process,"Times\t\t W\t\t SKN\t\t SKN_P\t\t SKN_A\t\t Time\n");
////////////////////////////////////////////////////////////////////////////
    int Pthread_State, Synchronize_State, Pthread_cond_State;
    Synchronize_Times=0;
    do
    {
        start_time_P = clock();
        Synchronize_State = StateSyn_Synchronize;
        for(p=0;p<Pthread_Number;++p)
        {
            pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][1]);  //1 is Pthread State.
            Pthread_State = *Signal[p + 1][1];
            pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][1]);//1 is Pthread State.
            if (Pthread_State != StatePthread_Stop)
            {
                start_time = clock();
                pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p+1][3] = SignalSyn_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][3]);//3 is Synchronize Signal.
                do{
                    pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][1]);//1 is Pthread State.

                    end_time = clock();
                    duration = (double)(end_time-start_time)/CLOCKS_PER_SEC;
                    /*
                    if(duration>5)
                    {
                        break;
                    }
                    */

                    /////////////////////////// Suspend thread by GUI /////////////////////////////
                    pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
                    while (GPM.Thread_State[0][1] == Signal_GUI_Break)//Signal_ContinueAndBreak
                    {
                        pthread_cond_wait(&GPM.Pthread_GUI_Cond[0], &GPM.Pthread_GUI_Mutex[0]);
                    }
                    pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);

                    /////////////////////////// Stop thread by GUI /////////////////////////////
                    pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
                    Signal_RunAndStop = GPM.Thread_State[0][0];
                    pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);
                    //Stop Signal
                    if (Signal_RunAndStop == Signal_GUI_Stop)
                    {
                        Pthread_State = StatePthread_Stop;
                        break;
                    }

                }while(Pthread_State!=StatePthread_PrepareToSynchronize && Pthread_State!=StatePthread_Stop);

                if (Pthread_State == StatePthread_PrepareToSynchronize && Pthread_State != StatePthread_Stop)
                {
                    pthread_mutex_lock(&M_Synchronize_NL_Mutex[p]);
                    for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
                    {
                        Point_Index   = MNS->Sort_Calculate_Points[q];
                        Boundary_Type = Local_Points_1D[Point_Index][0][0];
                        i             = Local_Points_1D[Point_Index][0][1];
                        j             = Local_Points_1D[Point_Index][0][2];
                        k             = Local_Points_1D[Point_Index][0][3];
                        v             = Local_Points_1D[Point_Index][0][4];
                        if(Boundary_Type!=0)
                        {
                            Mx_Global_Synchronize[i][j][k] = Mx_Global_Pthread[p][i][j][k];
                            My_Global_Synchronize[i][j][k] = My_Global_Pthread[p][i][j][k];
                            Mz_Global_Synchronize[i][j][k] = Mz_Global_Pthread[p][i][j][k];
                        }else
                        {
                            MxV_Global_Synchronize[v]      = MxV_Global_Pthread[p][v];
                            MyV_Global_Synchronize[v]      = MyV_Global_Pthread[p][v];
                            MzV_Global_Synchronize[v]      = MzV_Global_Pthread[p][v];
                        }
                    }

                    for(i=0;i<m;++i)
                    {
                        for(j=0;j<n;++j)
                        {
                            for(k=0;k<l;++k)
                            {
                                Mx_Global_Pthread[p][i][j][k] = Mx_Global_Synchronize[i][j][k];
                                My_Global_Pthread[p][i][j][k] = My_Global_Synchronize[i][j][k];
                                Mz_Global_Pthread[p][i][j][k] = Mz_Global_Synchronize[i][j][k];
                            }
                        }
                    }

                    for(v=0;v<Virtual_Points;++v)
                    {
                        MxV_Global_Pthread[p][v] = MxV_Global_Synchronize[v];
                        MyV_Global_Pthread[p][v] = MyV_Global_Synchronize[v];
                        MzV_Global_Pthread[p][v] = MzV_Global_Synchronize[v];
                    }
                    pthread_mutex_unlock(&M_Synchronize_NL_Mutex[p]);

                    ////////////////////// Calculate Energy //////////////////////////////////////////////////////
                    Energy_Average = 0.;    SKN = 0.;
                    for (p = 0; p < Pthread_Number; ++p)
                    {
                        Energy_Average   = Energy_Average + MNS->PTPS[p].Energy / N_Total;
                        SKN = SKN + MNS->PTPS[p].SKN;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][3]);  //3 is Synchronize Signal.
                    *Signal[p + 1][3] = SignalSyn_Run;
                    pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][3]);//3 is Synchronize Signal.
                    do
                    {
                        pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][1]);  //1 is Pthread State.
                        Pthread_State = *Signal[p + 1][1];
                        pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][1]);//1 is Pthread State.
                        if (Pthread_State != StatePthread_PrepareToSynchronize)
                        {
                            break;
                        }

                        pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][3]); //3 is Synchronize Signal.
                        Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_NL_Cond[p][3]);
                        pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][3]);//3 is Synchronize Signal.
                    } while (Pthread_cond_State != 0);
                }
                
            }else
            {
                Synchronize_State = StateSyn_Stop;
            }

            /////////////////////////// Stop thread by GUI /////////////////////////////
            pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
            Signal_RunAndStop = GPM.Thread_State[0][0];
            pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);
            //Stop Signal
            if (Signal_RunAndStop == Signal_GUI_Stop)
            {
                Synchronize_State = StateSyn_Stop;
                break;
            }
        }

        if(Synchronize_State==StateSyn_Stop)
        {
            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Stop;
                pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][3]);//3 is Synchronize Signal.
            }
        }
////////////////////////////////////////////////////////////////////////////////////////
        //Energy_Distribution(&EVD);
        ////////////////////// Calculate Energy //////////////////////////////////////////////////////
        Energy_Average = 0.;    SKN = 0.;
        for (p = 0; p < Pthread_Number; ++p)
        {
            Energy_Average = Energy_Average + MNS->PTPS[p].Energy / N_Total;
            SKN = SKN + MNS->PTPS[p].SKN;
        }
        pthread_mutex_lock(&OPGL_Mutex);
        for (i = 0; i < m; ++i)
        {
            for(j=0;j<n;++j)
            {
                for(k=0;k<l;++k)
                {
                    Point_Index = Position_3DTo1D[i][j][k];
                    Boundary_Type = Local_Points_1D[Point_Index][0][0];
                    if(Boundary_Type!=0)
                    {
                        Mx_OPGL[i][j][k] = Mx_Global_Synchronize[i][j][k];
                        My_OPGL[i][j][k] = My_Global_Synchronize[i][j][k];
                        Mz_OPGL[i][j][k] = Mz_Global_Synchronize[i][j][k];
                    }else
                    {
                        Mx_OPGL[i][j][k] = 0.;
                        My_OPGL[i][j][k] = 0.;
                        Mz_OPGL[i][j][k] = 0.;
                    }
                    
                }
            }
        }

        for(v=0;v<Virtual_Points;++v)
        {
            MxV_OPGL[v] = MxV_Global_Synchronize[v];
            MyV_OPGL[v] = MyV_Global_Synchronize[v];
            MzV_OPGL[v] = MzV_Global_Synchronize[v];
        }
        //Energy_Distribution(&EVD);
        OPGL_GUI.Energy = Energy_Average;
        OPGL_GUI.SKN = SKN;
        pthread_mutex_unlock(&OPGL_Mutex);
/////////////////////////////////////////////////////////////////////////////////////////
        if(Synchronize_Times%50==0)
        {
            for(k=0;k<l;k++)
            {
                sprintf(filenameMx,"Version %d\\Process\\z=%d\\Mx.dat",MNS->Version,k-zc);
                sprintf(filenameMy,"Version %d\\Process\\z=%d\\My.dat",MNS->Version,k-zc);
                sprintf(filenameMz,"Version %d\\Process\\z=%d\\Mz.dat",MNS->Version,k-zc);
                sprintf(filenamew, "Version %d\\Process\\z=%d\\W.dat",MNS->Version,k-zc);
                sprintf(filenameSKN_D, "Version %d\\Process\\z=%d\\SKN_D.dat",MNS->Version,k-zc);

                fp_Mx=fopen(filenameMx,"w+");
                fp_My=fopen(filenameMy,"w+");
                fp_Mz=fopen(filenameMz,"w+");
                fp_W =fopen(filenamew, "w+");
                fp_SKN_D=fopen(filenameSKN_D, "w+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fprintf(fp_Mx,"%0.8f\t",Mx_Global_Synchronize[i][j][k]);
                        fprintf(fp_My,"%0.8f\t",My_Global_Synchronize[i][j][k]);
                        fprintf(fp_Mz,"%0.8f\t",Mz_Global_Synchronize[i][j][k]);
                        fprintf(fp_W, "%0.8f\t",EVD.w[i][j][k]);
                        fprintf(fp_SKN_D, "%0.8f\t",EVD.SKN_D[i][j][k]);
                    }
                    fprintf(fp_Mx,"\n");
                    fprintf(fp_My,"\n");
                    fprintf(fp_Mz,"\n");
                    fprintf(fp_W, "\n");
                    fprintf(fp_SKN_D, "\n");
                }

                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
                fclose(fp_W);
                fclose(fp_SKN_D);
            }

            for(v=0;v<Virtual_Points;++v)
            {
                fprintf(fp_MxV,"%0.8f\t",MxV_Global_Synchronize[v]);
                fprintf(fp_MyV,"%0.8f\t",MyV_Global_Synchronize[v]);
                fprintf(fp_MzV,"%0.8f\t",MzV_Global_Synchronize[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);
        }
        end_time_P = clock();
        duration = (double)(end_time_P-start_time_P)/CLOCKS_PER_SEC;
        fprintf(fp_Process,"%d\t\t %0.8f\t %0.8f\t %0.8f\n",
                Synchronize_Times,Energy_Average,SKN,duration);
        printf("%d\t\t %0.8f\t %0.8f\t %0.8f\n",
               Synchronize_Times,Energy_Average,SKN,duration);

        Synchronize_Times=Synchronize_Times+1;
        
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
        Signal_ContinueAndBreak = GPM.Thread_State[0][1];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);

        /////////////////////////// Suspend thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
        while (GPM.Thread_State[0][1] == Signal_GUI_Break)//Signal_ContinueAndBreak
        {
            pthread_cond_wait(&GPM.Pthread_GUI_Cond[0], &GPM.Pthread_GUI_Mutex[0]);
        }
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);

        /////////////////////////// Stop thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
        Signal_RunAndStop = GPM.Thread_State[0][0];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);
        //Stop Signal
        if (Signal_RunAndStop == Signal_GUI_Stop)
        {
            break;
        }
    } while (Synchronize_State != StateSyn_Stop);

    for(p=0;p<Pthread_Number;++p)
    {
        pthread_mutex_lock(&Signal_Synchronize_NL_Mutex[p][3]);  //3 is Synchronize Signal.
        *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
        pthread_cond_signal(&Signal_Synchronize_NL_Cond[p][3]);
        pthread_mutex_unlock(&Signal_Synchronize_NL_Mutex[p][3]);//3 is Synchronize Signal.
    }
    fclose(fp_Process);
///////////////////////////////////////////////////////////////////////////////////////////
    for(p=0;p<Pthread_Number;++p)
    {
        for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
        {
            Point_Index   = MNS->Sort_Calculate_Points[q];
            Boundary_Type = Local_Points_1D[Point_Index][0][0];
            i             = Local_Points_1D[Point_Index][0][1];
            j             = Local_Points_1D[Point_Index][0][2];
            k             = Local_Points_1D[Point_Index][0][3];
            v             = Local_Points_1D[Point_Index][0][4];
            if(Boundary_Type!=0)
            {
                Mx_Global_Synchronize[i][j][k] = Mx_Global_Pthread[p][i][j][k];
                My_Global_Synchronize[i][j][k] = My_Global_Pthread[p][i][j][k];
                Mz_Global_Synchronize[i][j][k] = Mz_Global_Pthread[p][i][j][k];
            }else
            {
                MxV_Global_Synchronize[v]      = MxV_Global_Pthread[p][v];
                MyV_Global_Synchronize[v]      = MyV_Global_Pthread[p][v];
                MzV_Global_Synchronize[v]      = MzV_Global_Pthread[p][v];
            }
        }
    }

    //////////////////////// Stop OPGL //////////////////////////////////////////////////////////
    pthread_mutex_lock(&OPGL_Mutex);
    OPGL_GUI.State_Syn = StateSyn_Stop;
    pthread_mutex_unlock(&OPGL_Mutex);

    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    return 1;
}

void* Modify_Pthread_NL_Synchronize_Transfer(void *arg)
{
    Modify_NL_Synchronize *MNS;
    MNS=(Modify_NL_Synchronize*)arg;
    Modify_Pthread_NL_Synchronize(MNS);
    return NULL;
}

int OPGL_NL_Pthread(OPGL_NL_Pthread_Parameters *ONPP)
{
    OPGL_Parameters OPM;
    OPM.IPN = ONPP->IPN;
    OPM.Mx_OPGL = ONPP->Mx_OPGL;
    OPM.My_OPGL = ONPP->My_OPGL;
    OPM.Mz_OPGL = ONPP->Mz_OPGL;
    OPM.MxV_OPGL = ONPP->MxV_OPGL;
    OPM.MyV_OPGL = ONPP->MyV_OPGL;
    OPM.MzV_OPGL = ONPP->MzV_OPGL;
    OPM.CST = ONPP->CST;
    OPGL_NMD(&OPM);
    return 1;
}

void* OPGL_NL_Pthread_Transfer(void *arg)
{
    OPGL_NL_Pthread_Parameters *ONPP;
    ONPP=(OPGL_NL_Pthread_Parameters*)arg;
    OPGL_NL_Pthread(ONPP);
    return NULL;
}

int DiscreatPoint_Run_NL_Part_Run(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, Version, Threads_Number, Calculate_Points, CalculatePointsThreads;
    double h, t, *dxy;
    double ***M_Local, ****CST;
    int i_G, j_G, k_G, v, r, s, i, j, k, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    int ****Local_Points_1D_LocalEnergy, **Neighbour_Points_3Dto1D_Array, Point_Index;
    int Boundary_Type;

    clock_t start_time;
    start_time = clock();
    DRNP.start_time = start_time;
    Parameter_Transfer_Pthread_Synchronize *PTPS;

////////////////////////// Assignment ///////////////////
    xn  = IPN->xn;  yn  = IPN->yn;
    zn  = IPN->zn;  dxy = IPN->dxy;
    h   = IPN->h;   t   = IPN->t;
    Version    =IPN->Version;
    //NCycles    =IPN->NCycles;
    Threads_Number   = IPN->Threads_Number;
    PTPS = (Parameter_Transfer_Pthread_Synchronize *)malloc(sizeof(Parameter_Transfer_Pthread_Synchronize) * (Threads_Number));

    m = 2*xn+3; n = 2*yn+3;
    if(zn==0)
    {
        l=1;
    }else
    {
        l = 2*zn+3;
    }

///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape *BS;
    Calculate_Points_Parameters *CPP;

    BS=&DRNP.BS;

    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    Calculate_Points=BS->Calculate_Points;

    CPP = &DRNP.CPP;
    CPP->Position_3DTo1D  = Make3DArrayinteger(m,n,l);
    CPP->Local_Points_1D  = Make3DArrayinteger(BS->Calculate_Points,7,7);
    CPP->Virtual_Position = Make1DArrayinteger(BS->Virtual_Points);
    CPP->Local_Points_1D_LocalEnergy       = Make4DArrayinteger(BS->Calculate_Points,7,7,7);
    CPP->Neighbour_Points_3Dto1D_Array     = Make2DArrayinteger(14,2);
    CPP->Sort_Calculate_Points             = Make1DArrayinteger(BS->Calculate_Points);
    CalculPointsParameters_Initial(CPP,BS);
    Neighbour_Points_3Dto1D(CPP);
    Local_Points_1D_LocalEnergy   = CPP->Local_Points_1D_LocalEnergy;
    Neighbour_Points_3Dto1D_Array = CPP->Neighbour_Points_3Dto1D_Array;
    //ReSort_Calcualte_Points(CPP,BS);
    //ReSort_Calcualte_Points_ByThreads(CPP,BS,Threads_Number);
    //ReSort_Calcualte_Points_ByOrders(CPP,BS);
    //ReSort_Calcualte_Points_BySymmetry_Xaxis(CPP,BS);
    //ReSort_Calcualte_Points_BySymmetry_Yaxis(CPP,BS);
    ReSort_Calcualte_Points_ByRow(CPP,BS);

    Energy_Variables **EV;
    Modify_NL_Pthreads *MNPt;

    MNPt=(Modify_NL_Pthreads*)malloc(sizeof(Modify_NL_Pthreads)*(Threads_Number));

    M_Local       = Make3DArray(Threads_Number,4,14);
    CST           = Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);
    pthread_mutex_init(&OPGL_Mutex, 0);
    Signal_Synchronize_NL_Mutex = (pthread_mutex_t**)malloc(sizeof(pthread_mutex_t*)*(Threads_Number));
    Signal_Synchronize_NL_Cond  = (pthread_cond_t**)malloc(sizeof(pthread_cond_t*)*(Threads_Number));
    M_Synchronize_NL_Mutex      = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(Threads_Number));

    int ***Signal;
    Signal=Make2DArrayintegerPointer(Threads_Number+1,4);
    EV=(Energy_Variables**)malloc(sizeof(Energy_Variables*)*(Threads_Number));
    int p, q, o, I, J;
///////////////////////// Pthreads Parameters Initial /////////////////////////
    for(p=0;p<Threads_Number;++p)
    {
        //MNPt[p].CST           =Make4DArray(7,m,n,l);
        //Elastic_Field_Initial(MNPt[p].CST,IPN);
        MNPt[p].Pthread_Index =p;
        MNPt[p].NCycles       =IPN->NCycles;
        MNPt[p].Pthread_Number=Threads_Number;

        MNPt[p].Mx_Global_Pthread  = Make3DArray(m,n,l);
        MNPt[p].My_Global_Pthread  = Make3DArray(m,n,l);
        MNPt[p].Mz_Global_Pthread  = Make3DArray(m,n,l);

        MNPt[p].Sort_Calculate_Points = Make1DArrayinteger(BS->Calculate_Points);
        for(i=0;i<BS->Calculate_Points;++i)
        {
            MNPt[p].Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
        }
        MNPt[p].Local_Points_1D = CPP->Local_Points_1D;
        MNPt[p].Position_3DTo1D = CPP->Position_3DTo1D;
        MNPt[p].h   = h;    MNPt[p].t   = t;
        MNPt[p].m   = m;    MNPt[p].n   = n;
        MNPt[p].l   = l;
        MNPt[p].dxy = dxy;  MNPt[p].CST = CST;
        MNPt[p].Continuity_Modify_Mode = IPN->Continuity_Modify_Mode;
        MNPt[p].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(k=0;k<l;++k)
                {
                    MNPt[p].Mx_Global_Pthread[i][j][k]=Mx_Global[i][j][k];
                    MNPt[p].My_Global_Pthread[i][j][k]=My_Global[i][j][k];
                    MNPt[p].Mz_Global_Pthread[i][j][k]=Mz_Global[i][j][k];
                }
            }
        }

        MNPt[p].MxV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MNPt[p].MyV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MNPt[p].MzV_Global_Pthread = Make1DArray(BS->Virtual_Points);

        for(i=0;i<BS->Virtual_Points;++i)
        {
            MNPt[p].MxV_Global_Pthread[i]=MxV_Global[i];
            MNPt[p].MyV_Global_Pthread[i]=MyV_Global[i];
            MNPt[p].MzV_Global_Pthread[i]=MzV_Global[i];
        }

        CalculatePointsThreads_Lower=floor(Calculate_Points/Threads_Number)*p;
        if(p==Threads_Number-1)
        {
            CalculatePointsThreads       = Calculate_Points-floor(Calculate_Points/Threads_Number)*p;
            CalculatePointsThreads_Upper = Calculate_Points;
        }else
        {
            CalculatePointsThreads       = floor(Calculate_Points/Threads_Number);
            CalculatePointsThreads_Upper = floor(Calculate_Points/Threads_Number)*(p+1);
        }
        MNPt[p].CalculatePointsThreads_Lower=CalculatePointsThreads_Lower;
        MNPt[p].CalculatePointsThreads_Upper=CalculatePointsThreads_Upper;

        /////////////////////////// Synchronize Parameter Initial ////////////////////////////////
        MNPt[p].Signal          = Make1DArrayinteger(4);
        MNPt[p].Signal[1]       = 1;//0 is Pthread State.
        MNPt[p].Signal[2]       = 0;//1 is Pthread Signal.
        MNPt[p].Signal[3]       = 0;//2 is Synchronize Signal.
        Signal[p+1][1]          = &MNPt[p].Signal[1];
        Signal[p+1][2]          = &MNPt[p].Signal[2];
        Signal[p+1][3]          = &MNPt[p].Signal[3];
        //printf("%p\t%p\t%p\n",&MNPt[p].Signal[1],&MNPt[p].Signal[2],&MNPt[p].Signal[3]);
        //printf("%p\t%p\t%p\n",Signal[p+1][1],Signal[p+1][2],Signal[p+1][3]);
        //printf("%d\t%d\t%d\n",MNPt[p].Signal[1],MNPt[p].Signal[2],MNPt[p].Signal[3]);
        //printf("%d\t%d\t%d\n",*Signal[p+1][1],*Signal[p+1][2],*Signal[p+1][3]);
        Signal_Synchronize_NL_Mutex[p] = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(4));
        Signal_Synchronize_NL_Cond[p]  = (pthread_cond_t*)malloc(sizeof(pthread_cond_t)*(4));
        pthread_mutex_init(&Signal_Synchronize_NL_Mutex[p][1],0);//1 is Pthread State.
        pthread_mutex_init(&Signal_Synchronize_NL_Mutex[p][2],0);//2 is Pthread Signal.
        pthread_mutex_init(&Signal_Synchronize_NL_Mutex[p][3],0);//3 is Synchronize Signal.
        pthread_cond_init(&Signal_Synchronize_NL_Cond[p][1],NULL);//1 is Pthread State.
        pthread_cond_init(&Signal_Synchronize_NL_Cond[p][2],NULL);//2 is Pthread Signal.
        pthread_cond_init(&Signal_Synchronize_NL_Cond[p][3],NULL);//3 is Synchronize Signal.
        pthread_mutex_init(&M_Synchronize_NL_Mutex[p],0);

        /////////////////////////// Discrete Point Parameter Initial //////////////////////////
        MNPt[p].MNP    = (Modify_NL_Parameters*)malloc(sizeof(Modify_NL_Parameters)*(CalculatePointsThreads));
        EV[p]          = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));

        for(r=0;r<14;++r)
        {
            M_Local[p][1][r]=0.;//Mx
            M_Local[p][2][r]=0.;//My
            M_Local[p][3][r]=0.;//Mz
        }

        for(q=0;q<CalculatePointsThreads;++q)
        {
            MNPt[p].MNP[q].EV=&EV[p][q];
            EV[p][q].dxy=Make1DArray(3);
            EV[p][q].dxy[0] = dxy[0];    EV[p][q].dxy[1] = dxy[1];
            EV[p][q].dxy[2] = dxy[2];    EV[p][q].i      = 3;
            EV[p][q].j      = 3;         EV[p][q].k      = 3;
            EV[p][q].m      = m;         EV[p][q].n      = n;
            EV[p][q].l      = l;
            EV[p][q].h      = h;         EV[p][q].t      = t;
            EV[p][q].M0          = Make1DArray(3);
            EV[p][q].LocalEnergy = Make1DArrayinteger(7);
            EV[p][q].CST         = Make2DArray(7,7);
            EV[p][q].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
            EV[p][q].M_Local     = Make2DArrayPointer(4,15);
            EV[p][q].Energy_Coefficient = Make1DArray(N_Energy_Terms);

            for (r = 0; r < N_Energy_Terms; ++r)
            {
                EV[p][q].Energy_Coefficient[r] = IPN->Energy_Coefficient[r];
            }

            for(r=0;r<5;++r)
            {
                EV[p][q].M_Local[1][r] = &M_Local[p][1][r];
                EV[p][q].M_Local[2][r] = &M_Local[p][2][r];
                EV[p][q].M_Local[3][r] = &M_Local[p][3][r];
            }

            Point_Index=q+CalculatePointsThreads_Lower;
            Point_Index=CPP->Sort_Calculate_Points[Point_Index];
            //printf("%d\n",Point_Index);
/////////////////////////////////////////////////////////////////////////////////////////////////
            for(I=0;I<7;++I)
            {
                EV[p][q].LocalEnergy[I] = CPP->Local_Points_1D[Point_Index][I][5];

                Boundary_Type           = CPP->Local_Points_1D[Point_Index][I][0];
                i_G                     = CPP->Local_Points_1D[Point_Index][I][1];
                j_G                     = CPP->Local_Points_1D[Point_Index][I][2];
                k_G                     = CPP->Local_Points_1D[Point_Index][I][3];
                //printf("%d\t%d\t%d\t%d\t%d\n",Boundary_Type,i_G,j_G,k_G,v);
                for(J=0;J<7;++J)
                {
                    if(Boundary_Type!=0)
                    {
                        EV[p][q].CST[J][I]  = CST[J][i_G][j_G][k_G];
                    }else
                    {
                        if(CPP->Local_Points_1D[Point_Index][I][5]==1)
                        {
                            EV[p][q].CST[J][I] = CST[J][i_G][j_G][k_G];
                        }else
                        {
                            EV[p][q].CST[J][I] = 0.;
                        }
                    }
                }
            }

            for(o=1;o<14;++o)
            {
                r = Neighbour_Points_3Dto1D_Array[o][0];
                s = Neighbour_Points_3Dto1D_Array[o][1];

                Boundary_Type = Local_Points_1D_LocalEnergy[Point_Index][r][s][0];
                i_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][1];
                j_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][2];
                k_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][3];
                v             = Local_Points_1D_LocalEnergy[Point_Index][r][s][4];

                if(Boundary_Type!=0)
                {
                    EV[p][q].M_Local[1][o] = &MNPt[p].Mx_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[2][o] = &MNPt[p].My_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[3][o] = &MNPt[p].Mz_Global_Pthread[i_G][j_G][k_G];

                }else
                {
                    EV[p][q].M_Local[1][o] = &MNPt[p].MxV_Global_Pthread[v];
                    EV[p][q].M_Local[2][o] = &MNPt[p].MyV_Global_Pthread[v];
                    EV[p][q].M_Local[3][o] = &MNPt[p].MzV_Global_Pthread[v];
                }
            }
        }
    }

///////////////////////// Free Part /////////////////////////
    free3DArrayinteger(BS->Boundary,m,n);
    free4DArrayinteger(BS->LocalEnergy,m,n,l);
    free(CPP->Virtual_Position);
    free4DArrayinteger(Local_Points_1D_LocalEnergy,BS->Calculate_Points,7,7);
    free2DArrayinteger(CPP->Neighbour_Points_3Dto1D_Array,14);

///////////////////////// Synchronize Parameters Initial /////////////////////////
    Modify_NL_Synchronize *MNS;
    MNS = &DRNP.MNS;
    MNS->m=m;    MNS->n=n;    MNS->l=l;    MNS->Virtual_Points=BS->Virtual_Points;
    MNS->Pthread_Number=Threads_Number;
    MNS->Bound_CalculatePoint_Threads = Make2DArrayinteger(Threads_Number,3);
    MNS->Mx_Global_Synchronize        = Make3DArray(m,n,l);
    MNS->My_Global_Synchronize        = Make3DArray(m,n,l);
    MNS->Mz_Global_Synchronize        = Make3DArray(m,n,l);
    MNS->Signal                       = Signal;

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                MNS->Mx_Global_Synchronize[i][j][k]=Mx_Global[i][j][k];
                MNS->My_Global_Synchronize[i][j][k]=My_Global[i][j][k];
                MNS->Mz_Global_Synchronize[i][j][k]=Mz_Global[i][j][k];
            }
        }
    }
    MNS->h   = h;      MNS->t = t;
    MNS->dxy = dxy;
    MNS->Energy_Coefficient = IPN->Energy_Coefficient;
    MNS->MxV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MNS->MyV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MNS->MzV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MNS->Sort_Calculate_Points        = Make1DArrayinteger(BS->Calculate_Points);
    MNS->Virtual_Points               = BS->Virtual_Points;
    MNS->Version                      = Version;
    MNS->CST                          = CST;
    MNS->Local_Points_1D              = CPP->Local_Points_1D;
    MNS->Position_3DTo1D              = CPP->Position_3DTo1D;
    MNS->Mx_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MNS->My_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MNS->Mz_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MNS->MxV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MNS->MyV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MNS->MzV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);

    for(p=0;p<BS->Calculate_Points;++p)
    {
        MNS->Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
    }

    for(p=0;p<Threads_Number;++p)
    {
        MNS->Bound_CalculatePoint_Threads[p][0] = MNPt[p].CalculatePointsThreads_Upper-MNPt[p].CalculatePointsThreads_Lower;
        MNS->Bound_CalculatePoint_Threads[p][1] = MNPt[p].CalculatePointsThreads_Lower;
        MNS->Bound_CalculatePoint_Threads[p][2] = MNPt[p].CalculatePointsThreads_Upper;

        MNS->Mx_Global_Pthread[p]               = MNPt[p].Mx_Global_Pthread;
        MNS->My_Global_Pthread[p]               = MNPt[p].My_Global_Pthread;
        MNS->Mz_Global_Pthread[p]               = MNPt[p].Mz_Global_Pthread;
        MNS->MxV_Global_Pthread[p]              = MNPt[p].MxV_Global_Pthread;
        MNS->MyV_Global_Pthread[p]              = MNPt[p].MyV_Global_Pthread;
        MNS->MzV_Global_Pthread[p]              = MNPt[p].MzV_Global_Pthread;
    }

    ///////////////////////// Run Parameter /////////////////////////
    DRNP.CST = CST;
    DRNP.EV = EV;
    DRNP.M_Local = M_Local;
    DRNP.MNPt = MNPt;
    DRNP.Signal = Signal;
    DRNP.PTPS = PTPS;

    ///////////////////////// Pthreads Part /////////////////////////
    int err;
    pthread_t *tid;
    tid = (pthread_t *)malloc(sizeof(pthread_t) * (Threads_Number + 2));
    DRNP.tid = tid;
    for(p=0;p<Threads_Number;++p)
    {
        err=pthread_create(&tid[p],NULL,Modify_Pthread_NL_Transfer,&MNPt[p]);
        if(err != 0)
        {
            printf("create thread error\n");
            Sleep(1000);
        }
    }

    err=pthread_create(&tid[Threads_Number],NULL,Modify_Pthread_NL_Synchronize_Transfer,MNS);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }

    
    ///////////////////////// OPGL Part /////////////////////////
    OPGL_GUI.State_Syn = StateSyn_Run;
    OPGL_NL_Pthread_Parameters *ONPP;
    ONPP = &DRNP.ONPP;
    ONPP->Mx_OPGL = Mx_OPGL; 
    ONPP->My_OPGL = My_OPGL;
    ONPP->Mz_OPGL = Mz_OPGL;
    ONPP->MxV_OPGL = MxV_OPGL; 
    ONPP->MyV_OPGL = MyV_OPGL;
    ONPP->MzV_OPGL = MzV_OPGL;
    ONPP->IPN = IPN;
    ONPP->CST = Make4DArray(7,m,n,l);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                for (r = 0; r < 7; ++r)
                {
                    ONPP->CST[r][i][j][k] = CST[r][i][j][k];
                }
            }
        }
    }

    OPGL_Mode_GUI OMG;
    OMG = (OPGL_Mode_GUI)IPN->OPGL_Mode;
    switch (OMG)
    {
        case OPGL_Mode_On:
        {
            err = pthread_create(&tid[Threads_Number + 1], NULL, OPGL_NL_Pthread_Transfer, ONPP);
            if(err != 0)
            {
                printf("create thread error\n");
                Sleep(1000);
            }
        }
            
        default:
        break;
    }
    ////////////////////// Free pthread //////////////////////////////////////////
    pthread_t tid_Part_Free;
    err = pthread_create(&tid_Part_Free, NULL, DiscreatPoint_Run_NL_Part_Free_Transfer, IPN);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    
    return 1;
}

int DiscreatPoint_Run_NL_Part_Free(Input_Parameter_NMD *IPN)
{
    int m, n, l, i, j, k, v, p, q, CalculatePointsThreads, N_Thread;
    OPGL_Mode_GUI OMG;
    OMG = (OPGL_Mode_GUI)IPN->OPGL_Mode;
    switch (OMG)
    {
        case OPGL_Mode_On:
        {
            N_Thread = IPN->Threads_Number + 2;
            break;
        }

        case OPGL_Mode_Off:
        {
            N_Thread = IPN->Threads_Number + 1;
            break;
        }
            
        default:
        break;
    }
    printf("Free Thread\n");
    for (p = 0; p < N_Thread; ++p)
    {
        pthread_join(DRNP.tid[p], NULL);
    }
    double SKN, W, end_time, duration;
    m = 2 * IPN->xn + 3;    n = 2 * IPN->yn + 3;
    if (IPN->zn == 0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }
    
    //////////////////////////////// Assignment to the main Array ///////////////////////////////////
    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                Mx_Global[i][j][k] = DRNP.MNS.Mx_Global_Synchronize[i][j][k];
                My_Global[i][j][k] = DRNP.MNS.My_Global_Synchronize[i][j][k];
                Mz_Global[i][j][k] = DRNP.MNS.Mz_Global_Synchronize[i][j][k];
            }
        }
    }

    for (v = 0; v < DRNP.BS.Virtual_Points; ++v)
    {
        MxV_Global[v] = DRNP.MNS.MxV_Global_Synchronize[v];
        MyV_Global[v] = DRNP.MNS.MyV_Global_Synchronize[v];
        MzV_Global[v] = DRNP.MNS.MzV_Global_Synchronize[v];
    }
    ///////////////////////// Result Part /////////////////////////
    W=0.; SKN=0.;
    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        W   = W + DRNP.MNPt[p].W;
        SKN = SKN + DRNP.MNPt[p].SKN;
    }
    W = W / (DRNP.BS.Boundary_Points + DRNP.BS.Inner_Points);
    IPN->W = W; IPN->SKN = SKN;

    Energy_Variables_Distribution EVD;
    EVD.CST = DRNP.CST;    EVD.dxy = IPN->dxy;
    EVD.h   = IPN->h;      EVD.t   = IPN->t;
    EVD.m   = m;      EVD.n   = n;
    EVD.l     = l;
    EVD.Mx    = Mx_Global;
    EVD.My    = My_Global;
    EVD.Mz    = Mz_Global;
    EVD.MxV   = MxV_Global;
    EVD.MyV   = MyV_Global;
    EVD.MzV   = MzV_Global;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D = DRNP.CPP.Local_Points_1D;
    EVD.Position_3DTo1D = DRNP.CPP.Position_3DTo1D;
    Energy_Distribution(&EVD);
    Write_NMD_Energy(IPN,EVD.w);
    Write_NMD_SKN_D(IPN, EVD.SKN_D);
    IPN->W   = EVD.W_Average;
    IPN->SKN = EVD.SKN;
    IPN->SKN_Positif = EVD.SKN_Positif;
    IPN->SKN_Area = EVD.SKN_Area;
    printf("%lf\t%lf\t%lf\t%lf\n",IPN->W,IPN->SKN,IPN->SKN_Positif,IPN->SKN_Area);

    ///////////////////////// Free Part /////////////////////////

    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        free3DArray(DRNP.MNPt[p].Mx_Global_Pthread,m,n);
        free3DArray(DRNP.MNPt[p].My_Global_Pthread,m,n);
        free3DArray(DRNP.MNPt[p].Mz_Global_Pthread,m,n);
        free(DRNP.MNPt[p].MxV_Global_Pthread);
        free(DRNP.MNPt[p].MyV_Global_Pthread);
        free(DRNP.MNPt[p].MzV_Global_Pthread);
        free(DRNP.MNPt[p].Signal);

        CalculatePointsThreads = DRNP.MNPt[p].CalculatePointsThreads_Upper - DRNP.MNPt[p].CalculatePointsThreads_Lower;
        for(q=0;q<CalculatePointsThreads;++q)
        {
            free(DRNP.EV[p][q].dxy);
            free(DRNP.EV[p][q].M0);
            free(DRNP.EV[p][q].LocalEnergy);
            free2DArray(DRNP.EV[p][q].CST,7);
            free2DArrayPointer(DRNP.EV[p][q].M_Local,4);
            free(DRNP.EV[p][q].Energy_Coefficient);
        }

        for (q = 1; q < 4; ++q)
        {
            pthread_mutex_destroy(&Signal_Synchronize_NL_Mutex[p][q]);
            pthread_cond_destroy(&Signal_Synchronize_NL_Cond[p][q]);
        }

        free(DRNP.MNPt[p].MNP);
        free(DRNP.EV[p]);
        free(Signal_Synchronize_NL_Mutex[p]);
        free(DRNP.MNPt[p].Sort_Calculate_Points);
    }
    free4DArray(DRNP.CST,7,m,n);
    free2DArrayinteger(DRNP.MNS.Bound_CalculatePoint_Threads,IPN->Threads_Number);
    free3DArray(DRNP.MNS.Mx_Global_Synchronize,m,n);
    free3DArray(DRNP.MNS.My_Global_Synchronize,m,n);
    free3DArray(DRNP.MNS.Mz_Global_Synchronize,m,n);
    free(DRNP.MNS.MxV_Global_Synchronize);
    free(DRNP.MNS.MyV_Global_Synchronize);
    free(DRNP.MNS.MzV_Global_Synchronize);
    free(DRNP.MNS.Mx_Global_Pthread);
    free(DRNP.MNS.My_Global_Pthread);
    free(DRNP.MNS.Mz_Global_Pthread);
    free(DRNP.MNS.MxV_Global_Pthread);
    free(DRNP.MNS.MyV_Global_Pthread);
    free(DRNP.MNS.MzV_Global_Pthread);
    free(DRNP.MNS.Sort_Calculate_Points);
    free3DArray(DRNP.M_Local,IPN->Threads_Number,4);
    free3DArrayinteger(DRNP.CPP.Local_Points_1D,DRNP.BS.Calculate_Points,7);
    free3DArrayinteger(DRNP.CPP.Position_3DTo1D,m,n);
    free(DRNP.CPP.Sort_Calculate_Points);
    free(DRNP.MNPt);
    free(DRNP.EV);
    free2DArrayintegerPointer(DRNP.Signal,IPN->Threads_Number+1);
    free(Signal_Synchronize_NL_Mutex);
    free(M_Synchronize_NL_Mutex);
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free4DArray(DRNP.ONPP.CST,7,m,n);
    free(DRNP.tid);
    free(DRNP.PTPS);

    ///////////////////////// Result Part /////////////////////////
    end_time  = clock();
    DRNP.end_time = end_time;
    duration  = (double)(end_time-DRNP.start_time)/CLOCKS_PER_SEC;
    IPN->Time = duration;

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Finish;
    pthread_cond_signal(&GPM.Run_GUI_Cond);
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);

    printf("Free Done\n");
    return 1;
}

void* DiscreatPoint_Run_NL_Part_Free_Transfer(void *arg)
{
    Input_Parameter_NMD *IPN;
    IPN=(Input_Parameter_NMD*)arg;
    DiscreatPoint_Run_NL_Part_Free(IPN);
    return NULL;
}

int DiscreatPoint_Run_NL(Input_Parameter_NMD *IPN)
{
    int State_GUI_Run;
    Thread_State_GUI TSG;
    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    State_GUI_Run = GPM.State_GUI_Run;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    TSG = (Thread_State_GUI)State_GUI_Run;
    switch (TSG)
    {
        case TSG_Prepare:
        {
            pthread_mutex_lock(&GPM.Run_GUI_Mutex);
            GPM.State_GUI_Run = TSG_Run;
            pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
            DiscreatPoint_Run_NL_Part_Run(IPN);
            break;
        }

        case TSG_Run:
        {
            break;
        }

        case TSG_Finish:
        {
            break;
        }
        
        default:
        break;
    }
    return 1;
}

