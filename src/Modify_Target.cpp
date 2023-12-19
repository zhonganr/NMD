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
#include <Modify_Target.h>
#include <malloc.h>
#include <OPGL.h>
#include <Initial_State.h>

extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
pthread_mutex_t **Signal_Synchronize_Target_Mutex, *M_Synchronize_Target_Mutex;
extern pthread_mutex_t OPGL_Mutex;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;


using namespace std;


int Modify_Target_OneStep(Modify_Target_Parameters *MTP)
{
    //Assignment
    double ***M_Local, *M0, *M_Target, alpha=1, P;
    int End_Index, N_Cell=0, i;

    M0 = MTP->EV->M0;
    M_Target = MTP->EV->M_Target;
    M_Local = MTP->EV->M_Local;

    double Mx0, My0, Mz0, dMx, dMy, dMz;
    //double Mx1, My1, Mz1;
    //double M_Module0, M_Module1, theta0, theta1, phi0, phi1, M_Module;
    double W0, W1;

    M0[0] = *M_Local[1][1];
    M0[1] = *M_Local[2][1];
    M0[2] = *M_Local[3][1];
    Mx0   = *M_Local[1][1];
    My0   = *M_Local[2][1];
    Mz0   = *M_Local[3][1];
    dMx = M_Target[0] - Mx0;
    dMy = M_Target[1] - My0;
    dMz = M_Target[2] - Mz0;

    for (i = 0; i < 7; ++i)
    {
        if(MTP->EV->LocalEnergy[i]==1)
        {
            N_Cell = N_Cell + 1;
        }
    }
    /*
    M_Module0 = sqrt(Mx0 * Mx0 + My0 * My0 + Mz0 * Mz0);
    if(M_Module0==0)
    {
        M_Module0 = 1;
    }
    theta0 = acos(Mz0 / M_Module0);

    M_Module0 = sqrt(Mx0 * Mx0 + My0 * My0 + Mz0 * Mz0);
    phi0 = atan2(My0, Mx0);
    */
    W0=Energy_Local(MTP->EV)/N_Cell;
    W1 = W0;
    MTP->EV->W_Local_BC=W0;

    alpha = 1;
    P = 1.e-5;
/*    
    do{
        *M_Local[1][1] = Mx0 + alpha * dMx;
        *M_Local[2][1] = My0 + alpha * dMy;
        *M_Local[3][1] = Mz0 + alpha * dMz;

        Mx1 = Mx0 + alpha * dMx;
        My1 = My0 + alpha * dMy;
        Mz1 = Mz0 + alpha * dMz;
        M_Module1 = sqrt(Mx1 * Mx1 + My1 * My1 + Mz1 * Mz1);
        if(M_Module1==0)
        {
            M_Module1 = 1;
        }
        theta1 = acos(Mz1 / M_Module1);

        M_Module1 = sqrt(Mx1 * Mx1 + My1 * My1 + Mz1 * Mz1);
        phi1 = atan2(My1, Mx1);

        if(fabs(M_Module1-M_Module0)<0.1 && fabs(theta1-theta0)<0.1 && fabs(phi1-phi0)<0.2)
        {
            W1 = Energy_Local(MTP->EV);
            if(W1-W0<P)
            {
                break;
            }else
            {
                W1 = W0;
                alpha = 0.5 * alpha;
                *M_Local[1][1] = Mx0;
                *M_Local[2][1] = My0;
                *M_Local[3][1] = Mz0;
            }
        }else
        {
            W1 = W0;
            alpha = 0.5 * alpha;
            *M_Local[1][1] = Mx0;
            *M_Local[2][1] = My0;
            *M_Local[3][1] = Mz0;
        }
    }while (alpha >= 0.001);
*/
    do{
        *M_Local[1][1] = Mx0 + alpha * dMx;
        *M_Local[2][1] = My0 + alpha * dMy;
        *M_Local[3][1] = Mz0 + alpha * dMz;


        W1 = Energy_Local(MTP->EV)/N_Cell;
        if(W1-W0<P)
        {
                break;
        }else
        {
            W1 = W0;
            alpha = 0.5 * alpha;
            *M_Local[1][1] = Mx0;
            *M_Local[2][1] = My0;
            *M_Local[3][1] = Mz0;
        }

    }while (alpha >= 0.0001);

    //MTP->W_Change = MTP->W_Change + W1 - W0;
    MTP->W_Change = W1 - W0;
    M0[0] = *M_Local[1][1];
    M0[1] = *M_Local[2][1];
    M0[2] = *M_Local[3][1];
    if(W1-W0<0)
    {
        //printf("%0.8f\t%0.8f\t%0.8f\t%0.6f\n", W0, W1, W1-W0,alpha);
    }
    //printf("%0.8f\t%0.8f\t%0.8f\t%0.6f\n", W0, W1, W1-W0,alpha);
    return 1;
}

int Select_Points_Random(Modify_Target_Pthreads *MTPt)
{
    int N, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int Lower;
    CalculatePointsThreads_Lower = MTPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MTPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;

    if(CalculatePointsThreads_Lower==0)
    {
        Lower=1;
    }else
    {
        Lower=0;
    }

    static bool first = true;
    if(first)
    {
        srand((unsigned)time(NULL));
        first=false;
    }

    N = (rand()%(CalculatePointsThreads-Lower))+Lower;

    return N;
}

int Modify_Pthread_Target(Modify_Target_Pthreads *MTPt)
{
    int Pthread_Index, NCycles, *Signal;
    int Signal_Synchronize;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int i, j, N, Lower;
    double W1=0, W0=0;

    Pthread_Index  = MTPt->Pthread_Index;
    //Pthread_Number = MTPt->Pthread_Number;
    NCycles        = MTPt->NCycles;
    Signal         = MTPt->Signal;
    CalculatePointsThreads_Lower = MTPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MTPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;
    pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1]=1;
    pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][1]);//1 is Pthread State.

    if(CalculatePointsThreads_Lower==0)
    {
        Lower=1;
    }else
    {
        Lower=0;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
    for(j=Lower;j<CalculatePointsThreads;++j)
    {
        //MTPt->MTP[j].W_Change = 0;
        W0 = W0 + Energy_Single(MTPt->MTP[j].EV);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0;i<NCycles;++i)
    {
        pthread_mutex_lock(&M_Synchronize_Target_Mutex[Pthread_Index]);
        for(j=0;j<CalculatePointsThreads;++j)
        {
            N=Select_Points_Random(MTPt);
            //N = j;
            Modify_Target_OneStep(&MTPt->MTP[N]);
            //MTPt->Energy_Change[1] = MTPt->Energy_Change[1] + MTPt->MTP[N].W_Change;
        }
        W1 = 0;
        for(j=Lower;j<CalculatePointsThreads;++j)
        {
            //MTPt->MTP[j].W_Change = 0;
            W1 = W1 + Energy_Single(MTPt->MTP[j].EV);
        }

        MTPt->Energy_Change[1] = W1 - W0 + MTPt->Energy_Change[1];
        W0 = W1;
        if(MTPt->Energy_Change[1]>=MTPt->Energy_Change[2])//1 Actual 2 Max
        {
            MTPt->Energy_Change[2] = MTPt->Energy_Change[1];
        }
        
        pthread_mutex_unlock(&M_Synchronize_Target_Mutex[Pthread_Index]);

        pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
        Signal_Synchronize=Signal[3];
        pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.
        if(Signal_Synchronize==1)
        {
            pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2]=1;
            pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            do{
                pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
                Signal_Synchronize=Signal[3];
                pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.
            }while(Signal_Synchronize==1);
            pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2]=0;
            pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
        }else if(Signal_Synchronize==2)
        {
            pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2]=0;
            pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }
    }

    pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1]=0;
    pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[Pthread_Index][1]);//1 is Pthread State.

    printf("fin\t%d\n",Pthread_Index);
    return 1;
}

void* Modify_Pthread_Target_Transfer(void *arg)
{
    Modify_Target_Pthreads *MTPt;
    MTPt=(Modify_Target_Pthreads*)arg;
    Modify_Pthread_Target(MTPt);
    return NULL;
}

int Modify_Pthread_Target_Synchronize(Modify_Target_Synchronize *MTS)
{
    int i, j, k, v, p, q, Point_Index, Boundary_Type, Synchronize_Times;
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points;
    int Pthread_Number, ***Signal;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST;
    double *Energy_Change_Actual, *Energy_Change_Max;
    double W_Actual, W_Max, W0, W1, W_Original, Energy_Ratio;
    int N_Total=0;

    clock_t start_time, end_time, start_time_P, end_time_P;
    double duration;

    m=MTS->m;   n=MTS->n;   l=MTS->l;   Virtual_Points=MTS->Virtual_Points;
    Pthread_Number               = MTS->Pthread_Number;
    Bound_CalculatePoint_Threads = MTS->Bound_CalculatePoint_Threads;
    Local_Points_1D              = MTS->Local_Points_1D;
    Position_3DTo1D              = MTS->Position_3DTo1D;
    Mx_Global_Pthread            = MTS->Mx_Global_Pthread;
    My_Global_Pthread            = MTS->My_Global_Pthread;
    Mz_Global_Pthread            = MTS->Mz_Global_Pthread;
    Mx_Global_Synchronize        = MTS->Mx_Global_Synchronize;
    My_Global_Synchronize        = MTS->My_Global_Synchronize;
    Mz_Global_Synchronize        = MTS->Mz_Global_Synchronize;
    MxV_Global_Pthread           = MTS->MxV_Global_Pthread;
    MyV_Global_Pthread           = MTS->MyV_Global_Pthread;
    MzV_Global_Pthread           = MTS->MzV_Global_Pthread;
    MxV_Global_Synchronize       = MTS->MxV_Global_Synchronize;
    MyV_Global_Synchronize       = MTS->MyV_Global_Synchronize;
    MzV_Global_Synchronize       = MTS->MzV_Global_Synchronize;
    Signal                       = MTS->Signal;
    CST                          = MTS->CST;

    Energy_Change_Actual = Make1DArray(Pthread_Number);
    Energy_Change_Max    = Make1DArray(Pthread_Number);
    for (p = 0; p < Pthread_Number; ++p)
    {
        Energy_Change_Actual[p] = 0;
        Energy_Change_Max[p]    = 0;
    }
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
    EVD.CST = CST;    EVD.dxy = MTS->dxy;
    EVD.h   = MTS->h; EVD.t   = MTS->t;
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
    EVD.Local_Points_1D = Local_Points_1D;
    EVD.Position_3DTo1D = Position_3DTo1D;

    FILE *fp_Mx,*fp_My, *fp_Mz, *fp_MxV,*fp_MyV, *fp_MzV, *fp_W, *fp_SKN_D, *fp_Process;
    char filenameMx[500], filenameMy[500], filenameMz[500],filenamew[500],filenameSKN_D[500],
         filenameMxV[500],filenameMyV[500],filenameMzV[500],filenameP[500];

    sprintf(filenameMxV,"Version %d\\Process\\MxV.dat",MTS->Version);
    sprintf(filenameMyV,"Version %d\\Process\\MyV.dat",MTS->Version);
    sprintf(filenameMzV,"Version %d\\Process\\MzV.dat",MTS->Version);
    sprintf(filenameP, "Version %d\\Process.dat",MTS->Version);
    fp_MxV     = fopen(filenameMxV,"w+");
    fp_MyV     = fopen(filenameMyV,"w+");
    fp_MzV     = fopen(filenameMzV,"w+");
    fp_Process = fopen(filenameP, "w+");
    fprintf(fp_Process,"Times\t\t W\t\t SKN\t\t SKN_P\t\t SKN_A\t\t Time\n");

    Energy_Distribution(&EVD);
    W0 = EVD.W_Average;
    W_Original = W0;
    ////////////////////////////////////////////////////////////////////////////
    int Pthread_State, Synchronize_State, Signal_Pthread;
    Synchronize_Times=0;
    do
    {
        start_time_P = clock();
        Synchronize_State=1;
        for(p=0;p<Pthread_Number;++p)
        {
            pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][1]);  //1 is Pthread State.
            Pthread_State=*Signal[p+1][1];
            pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][1]);//1 is Pthread State.
            if(Pthread_State==1)
            {
                pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p+1][3] = 1;
                pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][3]);//3 is Synchronize Signal.

            }else
            {
                Synchronize_State=0;
            }
        }

        for(p=0;p<Pthread_Number;++p)
        {
            if(Synchronize_State==0)
            {
                pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p+1][3] = 2;                                    //Stop all threads
                pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][3]);//3 is Synchronize Signal.
            }

            start_time = clock();
            do{
                pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][2]);  //2 is Pthread Signal.
                Signal_Pthread = *Signal[p+1][2];
                pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][2]);//2 is Pthread Signal.
                pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][1]);  //1 is Pthread State.
                Pthread_State=*Signal[p+1][1];
                pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][1]);//1 is Pthread State.

                end_time = clock();
                duration = (double)(end_time-start_time)/CLOCKS_PER_SEC;
                if(duration>4 || Pthread_State==0)
                {
                    Signal_Pthread=1;
                }
            }while(Signal_Pthread==0);
        }

        N_Total = 0;
        for(p=0;p<Pthread_Number;++p)
        {
            pthread_mutex_lock(&M_Synchronize_Target_Mutex[p]);
            for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
            {
                Point_Index=MTS->Sort_Calculate_Points[q];
                Boundary_Type = Local_Points_1D[Point_Index][0][0];
                i             = Local_Points_1D[Point_Index][0][1];
                j             = Local_Points_1D[Point_Index][0][2];
                k             = Local_Points_1D[Point_Index][0][3];
                v             = Local_Points_1D[Point_Index][0][4];
                if(Local_Points_1D[Point_Index][0][5]==1)
                {
                    N_Total = N_Total + 1;
                }
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

            Energy_Change_Actual[p] = *MTS->Energy_change[p + 1][1];
            Energy_Change_Max[p]    = *MTS->Energy_change[p + 1][2];
            pthread_mutex_unlock(&M_Synchronize_Target_Mutex[p]);
        }

        for(p=0;p<Pthread_Number;++p)
        {
            pthread_mutex_lock(&M_Synchronize_Target_Mutex[p]);
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
            pthread_mutex_unlock(&M_Synchronize_Target_Mutex[p]);

            pthread_mutex_lock(&Signal_Synchronize_Target_Mutex[p][3]);  //3 is Synchronize Signal.
            *Signal[p+1][3] = 0;
            pthread_mutex_unlock(&Signal_Synchronize_Target_Mutex[p][3]);//3 is Synchronize Signal.
        }
////////////////////////////////////////////////////////////////////////////////////////
        pthread_mutex_lock(&OPGL_Mutex);
        for(i=0;i<m;++i)
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

        Energy_Distribution(&EVD);
        pthread_mutex_unlock(&OPGL_Mutex);
        W_Actual = 0.;  W_Max = 0.;
        //W1 = EVD.W_Average;
        //W_Actual = W1 - W0;
        if(W_Actual>W_Max)
        {
            //W_Max = W_Actual;
        }
        W0 = W1;
                
        for (p = 0; p < Pthread_Number; ++p)
        {
            W_Actual = W_Actual + Energy_Change_Actual[p]/N_Total;
            W_Max    = W_Max + Energy_Change_Max[p]/N_Total;
        }

        if(W_Max==0)
        {
            Energy_Ratio = 0;
        }else
        {
            Energy_Ratio = (EVD.W_Average - W_Original) / W_Max;
        }
        /////////////////////////////////////////////////////////////////////////////////////////
        //Energy_Distribution(&EVD);
        if (Synchronize_Times % 50 == 0)
        {
            for (k = 0; k < l; k++)
            {
                sprintf(filenameMx, "Version %d\\Process\\z=%d\\Mx.dat", MTS->Version, k - zc);
                sprintf(filenameMy, "Version %d\\Process\\z=%d\\My.dat", MTS->Version, k - zc);
                sprintf(filenameMz, "Version %d\\Process\\z=%d\\Mz.dat", MTS->Version, k - zc);
                sprintf(filenamew,  "Version %d\\Process\\z=%d\\W.dat",  MTS->Version, k - zc);
                sprintf(filenameSKN_D, "Version %d\\Process\\z=%d\\SKN_D.dat", MTS->Version, k - zc);

                fp_Mx = fopen(filenameMx, "w+");
                fp_My = fopen(filenameMy, "w+");
                fp_Mz = fopen(filenameMz, "w+");
                fp_W = fopen(filenamew, "w+");
                fp_SKN_D = fopen(filenameSKN_D, "w+");

                for (i = 0; i < m; i++)
                {
                    for (j = 0; j < n; j++)
                    {
                        fprintf(fp_Mx, "%0.8f\t", Mx_Global_Synchronize[i][j][k]);
                        fprintf(fp_My, "%0.8f\t", My_Global_Synchronize[i][j][k]);
                        fprintf(fp_Mz, "%0.8f\t", Mz_Global_Synchronize[i][j][k]);
                        fprintf(fp_W, "%0.8f\t", EVD.w[i][j][k]);
                        fprintf(fp_SKN_D, "%0.8f\t", EVD.SKN_D[i][j][k]);
                    }
                    fprintf(fp_Mx, "\n");
                    fprintf(fp_My, "\n");
                    fprintf(fp_Mz, "\n");
                    fprintf(fp_W, "\n");
                    fprintf(fp_SKN_D, "\n");
                }

                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
                fclose(fp_W);
                fclose(fp_SKN_D);
            }

            for (v = 0; v < Virtual_Points; ++v)
            {
                fprintf(fp_MxV, "%0.8f\t", MxV_Global_Synchronize[v]);
                fprintf(fp_MyV, "%0.8f\t", MyV_Global_Synchronize[v]);
                fprintf(fp_MzV, "%0.8f\t", MzV_Global_Synchronize[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);
        }
        end_time_P = clock();
        duration = (double)(end_time_P-start_time_P)/CLOCKS_PER_SEC;
        fprintf(fp_Process,"%d\t\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\n",
                Synchronize_Times,EVD.W_Average,EVD.SKN,EVD.SKN_Positif,W_Max,EVD.W_Average-W_Original,Energy_Ratio,duration);
        printf("%d\t\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t% 0.8f\t %0.8f\n",
               Synchronize_Times,EVD.W_Average,EVD.SKN,EVD.SKN_Positif,EVD.SKN_Area,W_Actual,W_Max,W_Original,Energy_Ratio,duration);

        Synchronize_Times=Synchronize_Times+1;
    }while(Synchronize_State==1);
    fclose(fp_Process);
///////////////////////////////////////////////////////////////////////////////////////////
    for(p=0;p<Pthread_Number;++p)
    {
        for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
        {
            Point_Index   = MTS->Sort_Calculate_Points[q];
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
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free(Energy_Change_Actual);
    free(Energy_Change_Max);
    return 1;
}

void* Modify_Pthread_Target_Synchronize_Transfer(void *arg)
{
    Modify_Target_Synchronize *MTS;
    MTS=(Modify_Target_Synchronize*)arg;
    Modify_Pthread_Target_Synchronize(MTS);
    return NULL;
}

int OPGL_Target_Pthread(OPGL_Target_Pthread_Parameters *OTPP)
{
    OPGL_Parameters OPT;
    OPT.IPN = OTPP->IPN;
    OPT.Mx_OPGL = OTPP->Mx_OPGL;
    OPT.My_OPGL = OTPP->My_OPGL;
    OPT.Mz_OPGL = OTPP->Mz_OPGL;
    OPT.MxV_OPGL = OTPP->MxV_OPGL;
    OPT.MyV_OPGL = OTPP->MyV_OPGL;
    OPT.MzV_OPGL = OTPP->MzV_OPGL;
    OPT.CST = OTPP->CST;
    OPGL_NMD(&OPT);
    return 1;
}

void* OPGL_Target_Pthread_Transfer(void *arg)
{
    OPGL_Target_Pthread_Parameters *OTPP;
    OTPP=(OPGL_Target_Pthread_Parameters*)arg;
    OPGL_Target_Pthread(OTPP);
    return NULL;
}

int DiscreatPoint_Run_Target(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, Version, Threads_Number, Calculate_Points, CalculatePointsThreads;
    double h, t, *dxy, W, SKN;
    double ***M_Local, ****CST;
    int i_G, j_G, k_G, v, r, s, i, j, k, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    int ****Local_Points_1D_LocalEnergy, **Neighbour_Points_3Dto1D_Array, Point_Index;
    int Boundary_Type;
    double ***Mx_Target, ***My_Target, ***Mz_Target, *MxV_Target, *MyV_Target, *MzV_Target;

    clock_t start_time, end_time;
    double duration;
    start_time = clock();

////////////////////////// Assignment ///////////////////
    xn  = IPN->xn;  yn  = IPN->yn;
    zn  = IPN->zn;  dxy = IPN->dxy;
    h   = IPN->h;   t   = IPN->t;
    Version    =IPN->Version;
    //NCycles    =IPN->NCycles;
    Threads_Number   = IPN->Threads_Number;

    m = 2*xn+3; n = 2*yn+3;
    if(zn==0)
    {
        l=1;
    }else
    {
        l = 2*zn+3;
    }

    Mx_Target = Make3DArray(m, n, l);
    My_Target = Make3DArray(m, n, l);
    Mz_Target = Make3DArray(m, n, l);
    ///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape BS0, *BS;
    Calculate_Points_Parameters CPP0, *CPP;

    BS=&BS0;
    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    Calculate_Points=BS->Calculate_Points;
    //printf("%d\n",Calculate_Points);

    CPP=&CPP0;
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
    ReSort_Calcualte_Points_BySymmetry_Xaxis(CPP,BS);

    MxV_Target = Make1DArray(BS->Virtual_Points);
    MyV_Target = Make1DArray(BS->Virtual_Points);
    MzV_Target = Make1DArray(BS->Virtual_Points);

    InitialM_Target(Mx_Target, My_Target, Mz_Target, MxV_Target, MyV_Target, MzV_Target, BS);

    Energy_Variables **EV;
    Modify_Target_Pthreads *MTPt;

    MTPt=(Modify_Target_Pthreads*)malloc(sizeof(Modify_Target_Pthreads)*(Threads_Number));

    M_Local       = Make3DArray(Threads_Number,4,14);
    CST           = Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);
    Signal_Synchronize_Target_Mutex = (pthread_mutex_t**)malloc(sizeof(pthread_mutex_t*)*(Threads_Number));
    M_Synchronize_Target_Mutex      = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(Threads_Number));

    int ***Signal;
    double ***Energy_Change;
    Signal=Make2DArrayintegerPointer(Threads_Number+1,4);
    Energy_Change = Make2DArrayPointer(Threads_Number + 1, 3);
    EV=(Energy_Variables**)malloc(sizeof(Energy_Variables*)*(Threads_Number));
    int p, q, o, I, J;
///////////////////////// Pthreads Parameters Initial /////////////////////////
    for(p=0;p<Threads_Number;++p)
    {
        //MTPt[p].CST           =Make4DArray(7,m,n,l);
        //Elastic_Field_Initial(MTPt[p].CST,IPN);
        MTPt[p].Pthread_Index  = p;
        MTPt[p].NCycles        = IPN->NCycles;
        MTPt[p].Pthread_Number = Threads_Number;

        MTPt[p].Mx_Global_Pthread  = Make3DArray(m,n,l);
        MTPt[p].My_Global_Pthread  = Make3DArray(m,n,l);
        MTPt[p].Mz_Global_Pthread  = Make3DArray(m,n,l);
        MTPt[p].Continuity_Modify_Mode = IPN->Continuity_Modify_Mode;
        MTPt[p].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;

        MTPt[p].Sort_Calculate_Points = Make1DArrayinteger(BS->Calculate_Points);
        for(i=0;i<BS->Calculate_Points;++i)
        {
            MTPt[p].Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
        }
        MTPt[p].Local_Points_1D = CPP->Local_Points_1D;
        MTPt[p].Position_3DTo1D = CPP->Position_3DTo1D;
        MTPt[p].h   = h;    MTPt[p].t   = t;
        MTPt[p].m   = m;    MTPt[p].n   = n;
        MTPt[p].l   = l;
        MTPt[p].dxy = dxy;  MTPt[p].CST = CST;

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(k=0;k<l;++k)
                {
                    MTPt[p].Mx_Global_Pthread[i][j][k]=Mx_Global[i][j][k];
                    MTPt[p].My_Global_Pthread[i][j][k]=My_Global[i][j][k];
                    MTPt[p].Mz_Global_Pthread[i][j][k]=Mz_Global[i][j][k];
                }
            }
        }

        MTPt[p].MxV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MTPt[p].MyV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MTPt[p].MzV_Global_Pthread = Make1DArray(BS->Virtual_Points);

        for(i=0;i<BS->Virtual_Points;++i)
        {
            MTPt[p].MxV_Global_Pthread[i]=MxV_Global[i];
            MTPt[p].MyV_Global_Pthread[i]=MyV_Global[i];
            MTPt[p].MzV_Global_Pthread[i]=MzV_Global[i];
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
        MTPt[p].CalculatePointsThreads_Lower=CalculatePointsThreads_Lower;
        MTPt[p].CalculatePointsThreads_Upper=CalculatePointsThreads_Upper;

        /////////////////////////// Synchronize Parameter Initial ////////////////////////////////
        MTPt[p].Signal           = Make1DArrayinteger(4);
        MTPt[p].Signal[1]        = 1;//0 is Pthread State.
        MTPt[p].Signal[2]        = 0;//1 is Pthread Signal.
        MTPt[p].Signal[3]        = 0;//2 is Synchronize Signal.
        Signal[p+1][1]           = &MTPt[p].Signal[1];
        Signal[p+1][2]           = &MTPt[p].Signal[2];
        Signal[p+1][3]           = &MTPt[p].Signal[3];

        MTPt[p].Energy_Change    = Make1DArray(3);
        MTPt[p].Energy_Change[1] = 0.;//Actual Change
        MTPt[p].Energy_Change[2] = 0.;//Max Change
        Energy_Change[p + 1][1]  = &MTPt[p].Energy_Change[1];
        Energy_Change[p + 1][2]  = &MTPt[p].Energy_Change[2];

        Signal_Synchronize_Target_Mutex[p] = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(4));
        pthread_mutex_init(&Signal_Synchronize_Target_Mutex[p][1],0);//1 is Pthread State.
        pthread_mutex_init(&Signal_Synchronize_Target_Mutex[p][2],0);//2 is Pthread Signal.
        pthread_mutex_init(&Signal_Synchronize_Target_Mutex[p][3],0);//3 is Synchronize Signal.
        pthread_mutex_init(&M_Synchronize_Target_Mutex[p],0);

        /////////////////////////// Discrete Point Parameter Initial //////////////////////////
        MTPt[p].MTP    = (Modify_Target_Parameters*)malloc(sizeof(Modify_Target_Parameters)*(CalculatePointsThreads));
        EV[p]          = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));

        for(r=0;r<14;++r)
        {
            M_Local[p][1][r]=0.;//Mx
            M_Local[p][2][r]=0.;//My
            M_Local[p][3][r]=0.;//Mz
        }

        for(q=0;q<CalculatePointsThreads;++q)
        {
            MTPt[p].MTP[q].EV = &EV[p][q];
            EV[p][q].dxy    = Make1DArray(3);
            EV[p][q].dxy[0] = dxy[0];    EV[p][q].dxy[1] = dxy[1];
            EV[p][q].dxy[2] = dxy[2];    EV[p][q].i      = 3;
            EV[p][q].j      = 3;         EV[p][q].k      = 3;
            EV[p][q].m      = m;         EV[p][q].n      = n;
            EV[p][q].l      = l;
            EV[p][q].h      = h;         EV[p][q].t      = t;
            EV[p][q].M0          = Make1DArray(3);
            EV[p][q].M_Target    = Make1DArray(3);
            EV[p][q].LocalEnergy = Make1DArrayinteger(7);
            EV[p][q].CST         = Make2DArray(7,7);
            EV[p][q].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
            EV[p][q].M_Local     = Make2DArrayPointer(4,15);

            for(r=0;r<5;++r)
            {
                EV[p][q].M_Local[1][r] = &M_Local[p][1][r];
                EV[p][q].M_Local[2][r] = &M_Local[p][2][r];
                EV[p][q].M_Local[3][r] = &M_Local[p][3][r];
            }

            Point_Index=q+CalculatePointsThreads_Lower;
            Point_Index=CPP->Sort_Calculate_Points[Point_Index];
            //printf("%d\n",Point_Index);
            i_G                  = CPP->Local_Points_1D[Point_Index][0][1];
            j_G                  = CPP->Local_Points_1D[Point_Index][0][2];
            k_G                  = CPP->Local_Points_1D[Point_Index][0][3];
            EV[p][q].M_Target[0] = Mx_Target[i_G][j_G][k_G];
            EV[p][q].M_Target[1] = My_Target[i_G][j_G][k_G];
            EV[p][q].M_Target[2] = Mz_Target[i_G][j_G][k_G];
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
                //printf("%d\t%d\t%d\n",o,r,s);

                Boundary_Type = Local_Points_1D_LocalEnergy[Point_Index][r][s][0];
                i_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][1];
                j_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][2];
                k_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][3];
                v             = Local_Points_1D_LocalEnergy[Point_Index][r][s][4];

                if(Boundary_Type!=0)
                {
                    EV[p][q].M_Local[1][o] = &MTPt[p].Mx_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[2][o] = &MTPt[p].My_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[3][o] = &MTPt[p].Mz_Global_Pthread[i_G][j_G][k_G];
                }else
                {
                    EV[p][q].M_Local[1][o] = &MTPt[p].MxV_Global_Pthread[v];
                    EV[p][q].M_Local[2][o] = &MTPt[p].MyV_Global_Pthread[v];
                    EV[p][q].M_Local[3][o] = &MTPt[p].MzV_Global_Pthread[v];
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
    Modify_Target_Synchronize MTS;
    MTS.m=m;    MTS.n=n;    MTS.l=l;    MTS.Virtual_Points=BS->Virtual_Points;
    MTS.Pthread_Number=Threads_Number;
    MTS.Bound_CalculatePoint_Threads = Make2DArrayinteger(Threads_Number,3);
    MTS.Mx_Global_Synchronize        = Make3DArray(m,n,l);
    MTS.My_Global_Synchronize        = Make3DArray(m,n,l);
    MTS.Mz_Global_Synchronize        = Make3DArray(m,n,l);
    MTS.Signal                       = Signal;
    MTS.Energy_change                = Energy_Change;
    /*    MTS.Signal=*Signal;
    for(p=0;p<Threads_Number;++p)
    {
        MTS.Signal[p]=*Signal[p];
    }*/

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                MTS.Mx_Global_Synchronize[i][j][k]=Mx_Global[i][j][k];
                MTS.My_Global_Synchronize[i][j][k]=My_Global[i][j][k];
                MTS.Mz_Global_Synchronize[i][j][k]=Mz_Global[i][j][k];
            }
        }
    }

    MTS.h = h;
    MTS.t = t;
    MTS.dxy = dxy;
    MTS.MxV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MTS.MyV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MTS.MzV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MTS.Sort_Calculate_Points        = Make1DArrayinteger(BS->Calculate_Points);
    MTS.Virtual_Points               = BS->Virtual_Points;
    MTS.Version                      = Version;
    MTS.CST                          = CST;
    MTS.Local_Points_1D              = CPP->Local_Points_1D;
    MTS.Position_3DTo1D              = CPP->Position_3DTo1D;
    MTS.Mx_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MTS.My_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MTS.Mz_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MTS.MxV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MTS.MyV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MTS.MzV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);

    for (v = 0; v < BS->Virtual_Points; ++v)
    {
        MTS.MxV_Global_Synchronize[v] = MxV_Global[v];
        MTS.MyV_Global_Synchronize[v] = MyV_Global[v];
        MTS.MzV_Global_Synchronize[v] = MzV_Global[v];
    }
    
    for(p=0;p<BS->Calculate_Points;++p)
    {
        MTS.Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
    }

    for(p=0;p<Threads_Number;++p)
    {
        MTS.Bound_CalculatePoint_Threads[p][0] = MTPt[p].CalculatePointsThreads_Upper-MTPt[p].CalculatePointsThreads_Lower;
        MTS.Bound_CalculatePoint_Threads[p][1] = MTPt[p].CalculatePointsThreads_Lower;
        MTS.Bound_CalculatePoint_Threads[p][2] = MTPt[p].CalculatePointsThreads_Upper;

        MTS.Mx_Global_Pthread[p]               = MTPt[p].Mx_Global_Pthread;
        MTS.My_Global_Pthread[p]               = MTPt[p].My_Global_Pthread;
        MTS.Mz_Global_Pthread[p]               = MTPt[p].Mz_Global_Pthread;
        MTS.MxV_Global_Pthread[p]              = MTPt[p].MxV_Global_Pthread;
        MTS.MyV_Global_Pthread[p]              = MTPt[p].MyV_Global_Pthread;
        MTS.MzV_Global_Pthread[p]              = MTPt[p].MzV_Global_Pthread;
    }

    ///////////////////////// Pthreads Part /////////////////////////
    int err;
    pthread_t tid[Threads_Number+2];
    for(p=0;p<Threads_Number;++p)
    {
        err=pthread_create(&tid[p],NULL,Modify_Pthread_Target_Transfer,&MTPt[p]);
        if(err != 0)
        {
            printf("create thread error\n");
            Sleep(1000);
        }
    }

    err=pthread_create(&tid[Threads_Number],NULL,Modify_Pthread_Target_Synchronize_Transfer,&MTS);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }

    ///////////////////////// OPGL Part /////////////////////////

    OPGL_Target_Pthread_Parameters OTPP;
    OTPP.Mx_OPGL = Mx_OPGL; 
    OTPP.My_OPGL = My_OPGL;
    OTPP.Mz_OPGL = Mz_OPGL;
    OTPP.MxV_OPGL = MxV_OPGL; 
    OTPP.MyV_OPGL = MyV_OPGL;
    OTPP.MzV_OPGL = MzV_OPGL;
    OTPP.IPN = IPN;
    OTPP.CST = Make4DArray(7,m,n,l);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                for (r = 0; r < 7; ++r)
                {
                    OTPP.CST[r][i][j][k] = CST[r][i][j][k];
                }
            }
        }
    }
    err=pthread_create(&tid[Threads_Number+1],NULL,OPGL_Target_Pthread_Transfer,&OTPP);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    /////////////////////////////////////////////////////////////////////
    for(p=0;p<Threads_Number+2;++p)
    {
        pthread_join(tid[p],NULL);
    }
    printf("1\n");

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                Mx_Global[i][j][k]=MTS.Mx_Global_Synchronize[i][j][k];
                My_Global[i][j][k]=MTS.My_Global_Synchronize[i][j][k];
                Mz_Global[i][j][k]=MTS.Mz_Global_Synchronize[i][j][k];
                //Mx_Global[i][j][k]=MTPt[0].Mx_Global_Pthread[i][j][k];
                //My_Global[i][j][k]=MTPt[0].My_Global_Pthread[i][j][k];
                //Mz_Global[i][j][k]=MTPt[0].Mz_Global_Pthread[i][j][k];
            }
        }
    }

    for(v=0;v<BS->Virtual_Points;++v)
    {
        MxV_Global[v]=MTS.MxV_Global_Synchronize[v];
        MyV_Global[v]=MTS.MyV_Global_Synchronize[v];
        MzV_Global[v]=MTS.MzV_Global_Synchronize[v];
    }
    ///////////////////////// Result Part /////////////////////////
    W=0.; SKN=0.;
    for(p=0;p<Threads_Number;++p)
    {
        W   = W+MTPt[p].W;
        SKN = SKN+MTPt[p].SKN;
    }
    W   = W/(BS->Boundary_Points+BS->Inner_Points);
    IPN->W = W; IPN->SKN = SKN;

    Energy_Variables_Distribution EVD;
    EVD.CST = CST;    EVD.dxy = dxy;
    EVD.h   = h;      EVD.t   = t;
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
    EVD.Local_Points_1D = CPP->Local_Points_1D;
    EVD.Position_3DTo1D = CPP->Position_3DTo1D;
    Energy_Distribution(&EVD);
    Write_NMD_Energy(IPN,EVD.w);
    Write_NMD_SKN_D(IPN,EVD.SKN_D);
    IPN->W   = EVD.W_Average;
    IPN->SKN = EVD.SKN;
    IPN->SKN_Positif = EVD.SKN_Positif;
    IPN->SKN_Area = EVD.SKN_Area;
    printf("%lf\t%lf\t%lf\t%lf\n",IPN->W,IPN->SKN,IPN->SKN_Positif,IPN->SKN_Area);

    ///////////////////////// Free Part /////////////////////////

    for(p=0;p<Threads_Number;++p)
    {
        free3DArray(MTPt[p].Mx_Global_Pthread,m,n);
        free3DArray(MTPt[p].My_Global_Pthread,m,n);
        free3DArray(MTPt[p].Mz_Global_Pthread,m,n);
        free(MTPt[p].MxV_Global_Pthread);
        free(MTPt[p].MyV_Global_Pthread);
        free(MTPt[p].MzV_Global_Pthread);
        free(MTPt[p].Signal);

        CalculatePointsThreads=MTPt[p].CalculatePointsThreads_Upper-MTPt[p].CalculatePointsThreads_Lower;
        for(q=0;q<CalculatePointsThreads;++q)
        {
            free(EV[p][q].dxy);
            free(EV[p][q].M0);
            free(EV[p][q].M_Target);
            free(EV[p][q].LocalEnergy);
            free2DArray(EV[p][q].CST,7);
            free2DArrayPointer(EV[p][q].M_Local,4);
        }
        free(MTPt[p].MTP);
        free(EV[p]);
        free(Signal_Synchronize_Target_Mutex[p]);
        free(MTPt[p].Sort_Calculate_Points);
    }
    free4DArray(CST,7,m,n);
    free2DArrayinteger(MTS.Bound_CalculatePoint_Threads,Threads_Number);
    free3DArray(MTS.Mx_Global_Synchronize,m,n);
    free3DArray(MTS.My_Global_Synchronize,m,n);
    free3DArray(MTS.Mz_Global_Synchronize,m,n);
    free(MTS.MxV_Global_Synchronize);
    free(MTS.MyV_Global_Synchronize);
    free(MTS.MzV_Global_Synchronize);
    free(MTS.Mx_Global_Pthread);
    free(MTS.My_Global_Pthread);
    free(MTS.Mz_Global_Pthread);
    free(MTS.MxV_Global_Pthread);
    free(MTS.MyV_Global_Pthread);
    free(MTS.MzV_Global_Pthread);
    free(MTS.Sort_Calculate_Points);
    free3DArray(M_Local,Threads_Number,4);
    free3DArrayinteger(CPP->Local_Points_1D,BS->Calculate_Points,7);
    free3DArrayinteger(CPP->Position_3DTo1D,m,n);
    free(CPP->Sort_Calculate_Points);
    free(MTPt);
    free(EV);
    free2DArrayintegerPointer(Signal,Threads_Number+1);
    free(Signal_Synchronize_Target_Mutex);
    free(M_Synchronize_Target_Mutex);
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free4DArray(OTPP.CST,7,m,n);
    free3DArray(Mx_Target,m,n);
    free3DArray(My_Target,m,n);
    free3DArray(Mz_Target,m,n);
    free(MxV_Target);
    free(MyV_Target);
    free(MzV_Target);
    ///////////////////////// Result Part /////////////////////////
    end_time  = clock();
    duration  = (double)(end_time-start_time)/CLOCKS_PER_SEC;
    IPN->Time = duration;

    return 1;
}

