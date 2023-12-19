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
#include <Modify_MC.h>
#include <OPGL.h>

extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
extern int N_Energy_Terms;
pthread_mutex_t **Signal_Synchronize_MC_Mutex, *M_Synchronize_MC_Mutex;
pthread_cond_t **Signal_Synchronize_MC_Cond;
extern pthread_mutex_t OPGL_Mutex;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
extern OPGL_Parameters_GUI OPGL_GUI;
extern GUI_Parameter GPM;
DiscreatPoint_Run_MC_Parameter DRMP;
using namespace std;

int Minimal(double *W, int n)
{
    int i, k;
    k=0;
    double W0=W[0];
    for(i=0;i<n;i++)
    {
        if(W0>=W[i])
        {
            W0=W[i];
            k=i;
        }
    }
    return k;
}

int Selection(int n_Mx, int n_My, int n_Mz, int n_dx, int n_dy, double *Cv, double *dM)
{
    if(n_Mx==1)
    {
        dM[0]=Cv[2];
    }else if(n_Mx==0)
    {
        dM[0]=0.;
    }else if(n_Mx==2)
    {
        dM[0]=-Cv[2];
    }

    if(n_My==1)
    {
        dM[1]=Cv[4];
    }else if(n_My==0)
    {
        dM[1]=0.;
    }else if(n_My==2)
    {
        dM[1]=-Cv[4];
    }

    if(n_Mz==1)
    {
        dM[2]=Cv[6];
    }else if(n_Mz==0)
    {
        dM[2]=0.;
    }else if(n_Mz==2)
    {
        dM[2]=-Cv[6];
    }

    if(n_dx==1)
    {
        dM[3]=Cv[8];
    }else if(n_dx==0)
    {
        dM[3]=0.;
    }else if(n_dx==2)
    {
        dM[3]=-Cv[8];
    }

    if(n_dy==1)
    {
        dM[4]=Cv[10];
    }else if(n_dy==0)
    {
        dM[4]=0.;
    }else if(n_dy==2)
    {
        dM[4]=-Cv[10];
    }

    return 1;
}

int Modify_MC_OneStep(Modify_MC_Parameters *MMP)
{
    //Assignment
    double ***M_Local, *M0, *Cv, M_Module;

    Cv = MMP->Cv;   M0 = MMP->EV->M0;
    M_Local = MMP->EV->M_Local;

    double  Mx0, My0, Mz0, dM[5];
    double *W_Array, W0;
    int r=0, n_Mx, n_My, n_Mz, n_total, k1=0;
    W_Array=Make1DArray(26);
    Cv[0]=Cv[0]+1;

    M0[0] = *M_Local[1][1];
    M0[1] = *M_Local[2][1];
    M0[2] = *M_Local[3][1];
    Mx0   = *M_Local[1][1];
    My0   = *M_Local[2][1];
    Mz0   = *M_Local[3][1];
/*
    M_Module = sqrt(Mx0 * Mx0 + My0 * My0 + Mz0 * Mz0);

    *M_Local[1][1] = (Mx0) / M_Module;
    *M_Local[2][1] = (My0) / M_Module;
    *M_Local[3][1] = (Mz0) / M_Module;

    M0[0] = *M_Local[1][1];
    M0[1] = *M_Local[2][1];
    M0[2] = *M_Local[3][1];
*/
    W0=Energy_Local(MMP->EV);
    MMP->EV->W_Local_BC=W0;

    for(n_total=0;n_total<27;n_total++)
    {
        n_Mx  =n_total%3;
        n_My  =floor(n_total/3);
        n_My  =n_My%3;
        n_Mz  =floor(floor(n_total/3)/3);
        n_Mz  =n_Mz%3;

        if(n_Mx!=0 || n_My!=0 || n_Mz!=0)
        {
            Selection(n_Mx,n_My,n_Mz,0,0,Cv,dM);
            //*M_Local[1][1]=Mx0+dM[0];
            //*M_Local[2][1]=My0+dM[1];
            //*M_Local[3][1]=Mz0+dM[2];

            M_Module = sqrt((Mx0 + dM[0]) * (Mx0 + dM[0]) + (My0 + dM[1]) * (My0 + dM[1]) + (Mz0 + dM[2]) * (Mz0 + dM[2]));
            if(fabs(M_Module)<=0.001)
            {
                M_Module = 1;
            }
            M_Module = 1;

            *M_Local[1][1]=(Mx0+dM[0])/M_Module;
            *M_Local[2][1]=(My0+dM[1])/M_Module;
            *M_Local[3][1]=(Mz0+dM[2])/M_Module;

            W_Array[r]=Energy_Local(MMP->EV);
            r=r+1;
        }
    }
    r=Minimal(W_Array,26);

    MMP->EV->W_Local_AC=W_Array[r];

    for(n_total=0;n_total<27;n_total++)
    {
        n_Mx  =n_total%3;
        n_My  =floor(n_total/3);
        n_My  =n_My%3;
        n_Mz  =floor(floor(n_total/3)/3);
        n_Mz  =n_Mz%3;

        if(n_Mx!=0 || n_My!=0 || n_Mz!=0)
        {
            if(k1==r)
            {
                Selection(n_Mx,n_My,n_Mz,0,0,Cv,dM);
                //*M_Local[1][1]=Mx0+dM[0];
                //*M_Local[2][1]=My0+dM[1];
                //*M_Local[3][1]=Mz0+dM[2];

                M_Module = sqrt((Mx0 + dM[0]) * (Mx0 + dM[0]) + (My0 + dM[1]) * (My0 + dM[1]) + (Mz0 + dM[2]) * (Mz0 + dM[2]));
                if(fabs(M_Module)<=0.001)
                {
                    M_Module = 1;
                }
                M_Module = 1;

                *M_Local[1][1] = (Mx0 + dM[0]) / M_Module;
                *M_Local[2][1] = (My0 + dM[1]) / M_Module;
                *M_Local[3][1] = (Mz0 + dM[2]) / M_Module;

                if(W_Array[r]<=W0)
                {
                    if(n_Mx!=0)
                    {
                        Cv[1]=Cv[1]+1;
                    }

                    if(n_My!=0)
                    {
                        Cv[3]=Cv[3]+1;
                    }

                    if(n_Mz!=0)
                    {
                        Cv[5]=Cv[5]+1;
                    }
                }
            }
            k1=k1+1;
        }
    }


    if(Cv[0]>=1000)
    {
        if(Cv[1]>=800)
        {
            Cv[2]=2*Cv[2];

        }else if(Cv[1]<=20)
        {
            if(Cv[2]>=4*1.e-10)
            {
                Cv[2]=0.5*Cv[2];
            }
        }

        if(Cv[3]>=800)
        {
            if(Cv[4]<=5)
            {
                Cv[4]=2*Cv[4];
            }

        }else if(Cv[3]<=20)
        {
            if(Cv[4]>=4*1.e-10)
            {
                Cv[4]=0.5*Cv[4];
            }
        }

        if(Cv[5]>=800)
        {
            if(Cv[6]<=5)
            {
                Cv[6]=2*Cv[6];
            }

        }else if(Cv[5]<=20)
        {
            if(Cv[6]>=4*1.e-10)
            {
                Cv[6]=0.5*Cv[6];
            }
        }


        Cv[0]=0; Cv[1]=0; Cv[3]=0; Cv[5]=0;
    }

    free(W_Array);
    //printf("%0.8f\t%0.8f\n",W0,MMP->W_AfterModify);
    return 1;
}

int Rotation(Modify_MC_Parameters *MMP)
{
    //Assignment
    double ***M_Local, *T, *M0;
    T  = MMP->T;
    M0 = MMP->EV->M0;
    M_Local = MMP->EV->M_Local;


    double P0, R=0.5, W0, W1, kb=0.5*1.e-7;
    T[1]=0;

    MMP->Energy_Variation=0;
    Modify_MC_OneStep(MMP);
    W0 = MMP->EV->W_Local_BC;
    //W   =W0;
    W1 = MMP->EV->W_Local_AC;
    //P   =exp(-(W1-W)/(T[0]*kb));
    P0 = exp(-(W1-W0)/(T[0]*kb));

    if(W1<W0)
    {
        //W0     =W1;
        //W      =W1;
        T[2]   =0;//Ttotal
        T[3]   =0;//Tn1
        R      =0.5;
        MMP->Energy_Variation=1;
    }

    if(P0>=R)
    {
        //W  =W1;
        MMP->EV->W_Local_BC = W1;
        if(R<0.9)
        {
            T[3]=T[3]+1;
        }
    }else
    {
        MMP->EV->W_Local_AC = W0;
        *M_Local[1][1]      = M0[0];
        *M_Local[2][1]      = M0[1];
        *M_Local[3][1]      = M0[2];
    }

    if(P0<0.3 && P0>=0.1)
    {
        R=0.6;
    }else if(P0<0.1 && P0>=0.01)
    {
        R=0.8;
    }else if(P0<0.01 && P0>=0.0001)
    {
        R=0.95;
    }else if(P0<0.0001)
    {
        R=0.99;
    }
    else
    {
        R=0.5;
    }
    T[1]=R;

    T[2]=T[2]+1;
    if(T[2]>=100)
    {
        if(T[3]<5)
        {
            T[0]=T[0]+1.;
            T[2]=0;
        }else if(T[3]>30 && T[0]>=2.)
        {
            T[0]=T[0]-1.;
            T[2]=0;
        }
    }

    return 1;
}

int Converge_Paramter_Initial(double *Cv, double *T)
{
    T[0]=14.; T[2]=0; T[3]=0;
    Cv[0]=0; Cv[1]=0; Cv[3]=0; Cv[5]=0; Cv[7]=0; Cv[9]=0;
    Cv[2]=0.001; Cv[4]=0.001; Cv[6]=0.001; Cv[8]=0.0001; Cv[10]=0.0001;

    return 1;
}

int Rotation_EnergyDecrease(Modify_MC_Parameters *MMP)
{
    //Assignment
    double ***M_Local, *M0;
    int  N;
    M_Local = MMP->EV->M_Local;


    double W0, W1;
    M0 = MMP->EV->M0;

    W0 = MMP->EV->W_Local_BC;

    for(N=0;N<100;++N)
    {
        Modify_MC_OneStep(MMP);
        W0 = MMP->EV->W_Local_BC;
        W1 = MMP->EV->W_Local_AC;
        if(W1<W0)
        {
            MMP->Energy_Variation=1;
        }else
        {
            if(N==0)
            {
                MMP->Energy_Variation=0;
            }
            *M_Local[1][1] = M0[0];
            *M_Local[2][1] = M0[1];
            *M_Local[3][1] = M0[2];
            break;
        }
    }
    //printf("%f\t %f\n", W0, W);
    return 1;
}

int Select_Points_Random(Modify_MC_Pthreads *MMPt)
{
    int N, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int Lower;
    CalculatePointsThreads_Lower = MMPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MMPt->CalculatePointsThreads_Upper;
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

int Modify_Pthread_MC(Modify_MC_Pthreads *MMPt)
{
    int Pthread_Index, NCycles, *Signal;
    int Signal_Synchronize, Signal_RunAndStop;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int i, j, N;

    Pthread_Index  = MMPt->Pthread_Index;
    //Pthread_Number = MMPt->Pthread_Number;
    NCycles        = MMPt->NCycles;
    Signal         = MMPt->Signal;
    CalculatePointsThreads_Lower = MMPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MMPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;
    pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Run;
    pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);//1 is Pthread State.
////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=0;i<NCycles;++i)
    {
        pthread_mutex_lock(&M_Synchronize_MC_Mutex[Pthread_Index]);
        for(j=0;j<CalculatePointsThreads;++j)
        {
            N=Select_Points_Random(MMPt);
            Rotation(&MMPt->MMP[N]);
            if(MMPt->MMP->Energy_Variation==1)
            {
                Rotation_EnergyDecrease(&MMPt->MMP[N]);
                //Neighbour_Points_Pool(N,MMPt);
            }

            if(MMPt->Continuity_Modify_Mode==1)
            {
                Continue_Modify(MMPt->MMP[N].EV);
            }
            
        }
        pthread_mutex_unlock(&M_Synchronize_MC_Mutex[Pthread_Index]);

        pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
        Signal_Synchronize = Signal[3];
        pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

        if (Signal_Synchronize == SignalSyn_Stop)
        {
            pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            if(Signal_Synchronize == SignalSyn_PrepareToSynchronize)
            {
                pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
                Signal[1] = StatePthread_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            }

            pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_Run && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_MC_Cond[Pthread_Index][3], &Signal_Synchronize_MC_Mutex[Pthread_Index][3]);
            }
            pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_Run;
            pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            
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
    printf("Thread\n");

    pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Stop;
    pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[Pthread_Index][1]);//1 is Pthread State.

    //free3DArray(EVD.w,EVD.m,EVD.n);
    printf("fin\t%d\n",Pthread_Index);
    return 1;
}

void* Modify_Pthread_MC_Transfer(void *arg)
{
    Modify_MC_Pthreads *MMPt;
    MMPt=(Modify_MC_Pthreads*)arg;
    Modify_Pthread_MC(MMPt);
    return NULL;
}

int Modify_Pthread_MC_Synchronize(Modify_MC_Synchronize *MMS)
{
    int i, j, k, v, p, q, Point_Index, Boundary_Type, Synchronize_Times;
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points;
    int Pthread_Number, ***Signal, Signal_RunAndStop, Signal_ContinueAndBreak;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST;

    clock_t start_time, end_time, start_time_P, end_time_P;
    double duration;

    m=MMS->m;   n=MMS->n;   l=MMS->l;   Virtual_Points=MMS->Virtual_Points;
    Pthread_Number               = MMS->Pthread_Number;
    Bound_CalculatePoint_Threads = MMS->Bound_CalculatePoint_Threads;
    Local_Points_1D              = MMS->Local_Points_1D;
    Position_3DTo1D              = MMS->Position_3DTo1D;
    Mx_Global_Pthread            = MMS->Mx_Global_Pthread;
    My_Global_Pthread            = MMS->My_Global_Pthread;
    Mz_Global_Pthread            = MMS->Mz_Global_Pthread;
    Mx_Global_Synchronize        = MMS->Mx_Global_Synchronize;
    My_Global_Synchronize        = MMS->My_Global_Synchronize;
    Mz_Global_Synchronize        = MMS->Mz_Global_Synchronize;
    MxV_Global_Pthread           = MMS->MxV_Global_Pthread;
    MyV_Global_Pthread           = MMS->MyV_Global_Pthread;
    MzV_Global_Pthread           = MMS->MzV_Global_Pthread;
    MxV_Global_Synchronize       = MMS->MxV_Global_Synchronize;
    MyV_Global_Synchronize       = MMS->MyV_Global_Synchronize;
    MzV_Global_Synchronize       = MMS->MzV_Global_Synchronize;
    Signal                       = MMS->Signal;
    CST                          = MMS->CST;
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
    EVD.CST = CST;    EVD.dxy = MMS->dxy;
    EVD.h   = MMS->h; EVD.t   = MMS->t;
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
    EVD.Energy_Coefficient = MMS->Energy_Coefficient;

    FILE *fp_Mx,*fp_My, *fp_Mz, *fp_MxV,*fp_MyV, *fp_MzV, *fp_W, *fp_SKN_D, *fp_Process;
    char filenameMx[500], filenameMy[500], filenameMz[500],filenamew[500],filenameSKN_D[500],
         filenameMxV[500],filenameMyV[500],filenameMzV[500],filenameP[500];

    sprintf(filenameMxV,"Version %d\\Process\\MxV.dat",MMS->Version);
    sprintf(filenameMyV,"Version %d\\Process\\MyV.dat",MMS->Version);
    sprintf(filenameMzV,"Version %d\\Process\\MzV.dat",MMS->Version);
    sprintf(filenameP, "Version %d\\Process.dat",MMS->Version);
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
            pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][1]);  //1 is Pthread State.
            Pthread_State = *Signal[p + 1][1];
            pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][1]);//1 is Pthread State.
            if (Pthread_State != StatePthread_Stop)
            {
                start_time = clock();
                pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p+1][3] = SignalSyn_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][3]);//3 is Synchronize Signal.
                do{
                    pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][1]);//1 is Pthread State.

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
                    pthread_mutex_lock(&M_Synchronize_MC_Mutex[p]);
                    for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
                    {
                        Point_Index   = MMS->Sort_Calculate_Points[q];
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
                    pthread_mutex_unlock(&M_Synchronize_MC_Mutex[p]);

                    pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][3]);  //3 is Synchronize Signal.
                    *Signal[p + 1][3] = SignalSyn_Run;
                    pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][3]);//3 is Synchronize Signal.
                    do
                    {
                        pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][1]);  //1 is Pthread State.
                        Pthread_State = *Signal[p + 1][1];
                        pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][1]);//1 is Pthread State.
                        if (Pthread_State != StatePthread_PrepareToSynchronize)
                        {
                            break;
                        }

                        pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][3]); //3 is Synchronize Signal.
                        Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_MC_Cond[p][3]);
                        pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][3]);//3 is Synchronize Signal.
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
                pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Stop;
                pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][3]);//3 is Synchronize Signal.
            }
        }
////////////////////////////////////////////////////////////////////////////////////////
        Energy_Distribution(&EVD);
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
        //Energy_Distribution(&EVD);
        OPGL_GUI.Energy = EVD.W_Average;
        OPGL_GUI.SKN = EVD.SKN;
        pthread_mutex_unlock(&OPGL_Mutex);
/////////////////////////////////////////////////////////////////////////////////////////
        if(Synchronize_Times%50==0)
        {
            for(k=0;k<l;k++)
            {
                sprintf(filenameMx,"Version %d\\Process\\z=%d\\Mx.dat",MMS->Version,k-zc);
                sprintf(filenameMy,"Version %d\\Process\\z=%d\\My.dat",MMS->Version,k-zc);
                sprintf(filenameMz,"Version %d\\Process\\z=%d\\Mz.dat",MMS->Version,k-zc);
                sprintf(filenamew, "Version %d\\Process\\z=%d\\W.dat",MMS->Version,k-zc);
                sprintf(filenameSKN_D, "Version %d\\Process\\z=%d\\SKN_D.dat",MMS->Version,k-zc);

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
        fprintf(fp_Process,"%d\t\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\n",
                Synchronize_Times,EVD.W_Average,EVD.SKN,EVD.SKN_Positif,EVD.SKN_Area,duration);
        printf("%d\t\t %0.8f\t %0.8f\t %0.8f\t %0.8f\t %0.8f\n",
               Synchronize_Times,EVD.W_Average,EVD.SKN,EVD.SKN_Positif,EVD.SKN_Area,duration);

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
        pthread_mutex_lock(&Signal_Synchronize_MC_Mutex[p][3]);  //3 is Synchronize Signal.
        *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
        pthread_cond_signal(&Signal_Synchronize_MC_Cond[p][3]);
        pthread_mutex_unlock(&Signal_Synchronize_MC_Mutex[p][3]);//3 is Synchronize Signal.
    }
    fclose(fp_Process);
///////////////////////////////////////////////////////////////////////////////////////////
    for(p=0;p<Pthread_Number;++p)
    {
        for(q=Bound_CalculatePoint_Threads[p][1];q<Bound_CalculatePoint_Threads[p][2];++q)
        {
            Point_Index   = MMS->Sort_Calculate_Points[q];
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

void* Modify_Pthread_MC_Synchronize_Transfer(void *arg)
{
    Modify_MC_Synchronize *MMS;
    MMS=(Modify_MC_Synchronize*)arg;
    Modify_Pthread_MC_Synchronize(MMS);
    return NULL;
}

int OPGL_MC_Pthread(OPGL_MC_Pthread_Parameters *OMPP)
{
    OPGL_Parameters OPM;
    OPM.IPN = OMPP->IPN;
    OPM.Mx_OPGL = OMPP->Mx_OPGL;
    OPM.My_OPGL = OMPP->My_OPGL;
    OPM.Mz_OPGL = OMPP->Mz_OPGL;
    OPM.MxV_OPGL = OMPP->MxV_OPGL;
    OPM.MyV_OPGL = OMPP->MyV_OPGL;
    OPM.MzV_OPGL = OMPP->MzV_OPGL;
    OPM.CST = OMPP->CST;
    OPGL_NMD(&OPM);
    return 1;
}

void* OPGL_MC_Pthread_Transfer(void *arg)
{
    OPGL_MC_Pthread_Parameters *OMPP;
    OMPP=(OPGL_MC_Pthread_Parameters*)arg;
    OPGL_MC_Pthread(OMPP);
    return NULL;
}

int DiscreatPoint_Run_MC_Part_Run(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, Version, Threads_Number, Calculate_Points, CalculatePointsThreads;
    double h, t, *dxy;
    double ***M_Local, ****CST;
    int i_G, j_G, k_G, v, r, s, i, j, k, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    int ****Local_Points_1D_LocalEnergy, **Neighbour_Points_3Dto1D_Array, Point_Index;
    int Boundary_Type;

    clock_t start_time;
    start_time = clock();
    DRMP.start_time = start_time;

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

///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape *BS;
    Calculate_Points_Parameters *CPP;

    BS=&DRMP.BS;
    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    Calculate_Points=BS->Calculate_Points;
    //printf("%d\n",Calculate_Points);

    CPP=&DRMP.CPP;
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
    ReSort_Calcualte_Points(CPP,BS);
    //ReSort_Calcualte_Points_ByThreads(CPP,BS,Threads_Number);
    //ReSort_Calcualte_Points_ByOrders(CPP,BS);
    //ReSort_Calcualte_Points_BySymmetry_Xaxis(CPP,BS);

    Energy_Variables **EV;
    Modify_MC_Pthreads *MMPt;

    MMPt=(Modify_MC_Pthreads*)malloc(sizeof(Modify_MC_Pthreads)*(Threads_Number));

    M_Local       = Make3DArray(Threads_Number,4,14);
    CST           = Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);
    pthread_mutex_init(&OPGL_Mutex, 0);
    Signal_Synchronize_MC_Mutex = (pthread_mutex_t**)malloc(sizeof(pthread_mutex_t*)*(Threads_Number));
    Signal_Synchronize_MC_Cond  = (pthread_cond_t**)malloc(sizeof(pthread_cond_t*)*(Threads_Number));
    M_Synchronize_MC_Mutex      = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(Threads_Number));

    int ***Signal;
    Signal=Make2DArrayintegerPointer(Threads_Number+1,4);
    EV=(Energy_Variables**)malloc(sizeof(Energy_Variables*)*(Threads_Number));
    int p, q, o, I, J;
///////////////////////// Pthreads Parameters Initial /////////////////////////
    for(p=0;p<Threads_Number;++p)
    {
        //MMPt[p].CST           =Make4DArray(7,m,n,l);
        //Elastic_Field_Initial(MMPt[p].CST,IPN);
        MMPt[p].Pthread_Index  = p;
        MMPt[p].NCycles        = IPN->NCycles;
        MMPt[p].Pthread_Number = Threads_Number;

        MMPt[p].Mx_Global_Pthread  = Make3DArray(m,n,l);
        MMPt[p].My_Global_Pthread  = Make3DArray(m,n,l);
        MMPt[p].Mz_Global_Pthread  = Make3DArray(m,n,l);
        MMPt[p].Cv                 = Make1DArray(11);
        MMPt[p].T                  = Make1DArray(5);
        MMPt[p].Continuity_Modify_Mode = IPN->Continuity_Modify_Mode;
        MMPt[p].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;

        Converge_Paramter_Initial(MMPt[p].Cv,MMPt[p].T);

        MMPt[p].Sort_Calculate_Points = Make1DArrayinteger(BS->Calculate_Points);
        for(i=0;i<BS->Calculate_Points;++i)
        {
            MMPt[p].Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
        }
        MMPt[p].Local_Points_1D = CPP->Local_Points_1D;
        MMPt[p].Position_3DTo1D = CPP->Position_3DTo1D;
        MMPt[p].h   = h;    MMPt[p].t   = t;
        MMPt[p].m   = m;    MMPt[p].n   = n;
        MMPt[p].l   = l;
        MMPt[p].dxy = dxy;  MMPt[p].CST = CST;

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(k=0;k<l;++k)
                {
                    MMPt[p].Mx_Global_Pthread[i][j][k]=Mx_Global[i][j][k];
                    MMPt[p].My_Global_Pthread[i][j][k]=My_Global[i][j][k];
                    MMPt[p].Mz_Global_Pthread[i][j][k]=Mz_Global[i][j][k];
                }
            }
        }

        MMPt[p].MxV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MMPt[p].MyV_Global_Pthread = Make1DArray(BS->Virtual_Points);
        MMPt[p].MzV_Global_Pthread = Make1DArray(BS->Virtual_Points);

        for(i=0;i<BS->Virtual_Points;++i)
        {
            MMPt[p].MxV_Global_Pthread[i]=MxV_Global[i];
            MMPt[p].MyV_Global_Pthread[i]=MyV_Global[i];
            MMPt[p].MzV_Global_Pthread[i]=MzV_Global[i];
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
        MMPt[p].CalculatePointsThreads_Lower=CalculatePointsThreads_Lower;
        MMPt[p].CalculatePointsThreads_Upper=CalculatePointsThreads_Upper;

        /////////////////////////// Synchronize Parameter Initial ////////////////////////////////
        MMPt[p].Signal          = Make1DArrayinteger(4);
        MMPt[p].Signal[1]       = 1;//0 is Pthread State.
        MMPt[p].Signal[2]       = 0;//1 is Pthread Signal.
        MMPt[p].Signal[3]       = 0;//2 is Synchronize Signal.
        Signal[p+1][1]          = &MMPt[p].Signal[1];
        Signal[p+1][2]          = &MMPt[p].Signal[2];
        Signal[p+1][3]          = &MMPt[p].Signal[3];
        //printf("%p\t%p\t%p\n",&MMPt[p].Signal[1],&MMPt[p].Signal[2],&MMPt[p].Signal[3]);
        //printf("%p\t%p\t%p\n",Signal[p+1][1],Signal[p+1][2],Signal[p+1][3]);
        //printf("%d\t%d\t%d\n",MMPt[p].Signal[1],MMPt[p].Signal[2],MMPt[p].Signal[3]);
        //printf("%d\t%d\t%d\n",*Signal[p+1][1],*Signal[p+1][2],*Signal[p+1][3]);
        Signal_Synchronize_MC_Mutex[p] = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(4));
        Signal_Synchronize_MC_Cond[p]  = (pthread_cond_t*)malloc(sizeof(pthread_cond_t)*(4));
        pthread_mutex_init(&Signal_Synchronize_MC_Mutex[p][1],0);//1 is Pthread State.
        pthread_mutex_init(&Signal_Synchronize_MC_Mutex[p][2],0);//2 is Pthread Signal.
        pthread_mutex_init(&Signal_Synchronize_MC_Mutex[p][3],0);//3 is Synchronize Signal.
        pthread_cond_init(&Signal_Synchronize_MC_Cond[p][1],NULL);//1 is Pthread State.
        pthread_cond_init(&Signal_Synchronize_MC_Cond[p][2],NULL);//2 is Pthread Signal.
        pthread_cond_init(&Signal_Synchronize_MC_Cond[p][3],NULL);//3 is Synchronize Signal.
        pthread_mutex_init(&M_Synchronize_MC_Mutex[p],0);

        /////////////////////////// Discrete Point Parameter Initial //////////////////////////
        MMPt[p].MMP    = (Modify_MC_Parameters*)malloc(sizeof(Modify_MC_Parameters)*(CalculatePointsThreads));
        EV[p]          = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));

        for(r=0;r<14;++r)
        {
            M_Local[p][1][r]=0.;//Mx
            M_Local[p][2][r]=0.;//My
            M_Local[p][3][r]=0.;//Mz
        }

        for(q=0;q<CalculatePointsThreads;++q)
        {
            MMPt[p].MMP[q].Cv = MMPt[p].Cv;
            MMPt[p].MMP[q].T  = MMPt[p].T;
            MMPt[p].MMP[q].EV = &EV[p][q];
            EV[p][q].dxy    = Make1DArray(3);
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
                //printf("%d\t%d\t%d\n",o,r,s);

                Boundary_Type = Local_Points_1D_LocalEnergy[Point_Index][r][s][0];
                i_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][1];
                j_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][2];
                k_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][3];
                v             = Local_Points_1D_LocalEnergy[Point_Index][r][s][4];

                if(Boundary_Type!=0)
                {
                    EV[p][q].M_Local[1][o] = &MMPt[p].Mx_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[2][o] = &MMPt[p].My_Global_Pthread[i_G][j_G][k_G];
                    EV[p][q].M_Local[3][o] = &MMPt[p].Mz_Global_Pthread[i_G][j_G][k_G];

                }else
                {
                    EV[p][q].M_Local[1][o] = &MMPt[p].MxV_Global_Pthread[v];
                    EV[p][q].M_Local[2][o] = &MMPt[p].MyV_Global_Pthread[v];
                    EV[p][q].M_Local[3][o] = &MMPt[p].MzV_Global_Pthread[v];
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
    Modify_MC_Synchronize *MMS;
    MMS = &DRMP.MMS;
    MMS->m=m;    MMS->n=n;    MMS->l=l;    MMS->Virtual_Points=BS->Virtual_Points;
    MMS->Pthread_Number=Threads_Number;
    MMS->Bound_CalculatePoint_Threads = Make2DArrayinteger(Threads_Number,3);
    MMS->Mx_Global_Synchronize        = Make3DArray(m,n,l);
    MMS->My_Global_Synchronize        = Make3DArray(m,n,l);
    MMS->Mz_Global_Synchronize        = Make3DArray(m,n,l);
    MMS->Signal                       = Signal;
/*    MMS.Signal=*Signal;
    for(p=0;p<Threads_Number;++p)
    {
        MMS.Signal[p]=*Signal[p];
    }*/

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                MMS->Mx_Global_Synchronize[i][j][k]=Mx_Global[i][j][k];
                MMS->My_Global_Synchronize[i][j][k]=My_Global[i][j][k];
                MMS->Mz_Global_Synchronize[i][j][k]=Mz_Global[i][j][k];
            }
        }
    }


    MMS->h = h;
    MMS->t = t;
    MMS->dxy = dxy;
    MMS->Energy_Coefficient = IPN->Energy_Coefficient;
    MMS->MxV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MMS->MyV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MMS->MzV_Global_Synchronize       = Make1DArray(BS->Virtual_Points);
    MMS->Sort_Calculate_Points        = Make1DArrayinteger(BS->Calculate_Points);
    MMS->Virtual_Points               = BS->Virtual_Points;
    MMS->Version                      = Version;
    MMS->CST                          = CST;
    MMS->Local_Points_1D              = CPP->Local_Points_1D;
    MMS->Position_3DTo1D              = CPP->Position_3DTo1D;
    MMS->Mx_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MMS->My_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MMS->Mz_Global_Pthread            = (double****)malloc(sizeof(double***)*Threads_Number);
    MMS->MxV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MMS->MyV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);
    MMS->MzV_Global_Pthread           = (double**)malloc(sizeof(double*)*Threads_Number);

    for (v = 0; v < BS->Virtual_Points; ++v)
    {
        MMS->MxV_Global_Synchronize[v]=MxV_Global[v];
        MMS->MyV_Global_Synchronize[v]=MyV_Global[v];
        MMS->MzV_Global_Synchronize[v]=MzV_Global[v];
    }

    for(p=0;p<BS->Calculate_Points;++p)
    {
        MMS->Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
    }

    for(p=0;p<Threads_Number;++p)
    {
        MMS->Bound_CalculatePoint_Threads[p][0] = MMPt[p].CalculatePointsThreads_Upper-MMPt[p].CalculatePointsThreads_Lower;
        MMS->Bound_CalculatePoint_Threads[p][1] = MMPt[p].CalculatePointsThreads_Lower;
        MMS->Bound_CalculatePoint_Threads[p][2] = MMPt[p].CalculatePointsThreads_Upper;

        MMS->Mx_Global_Pthread[p]               = MMPt[p].Mx_Global_Pthread;
        MMS->My_Global_Pthread[p]               = MMPt[p].My_Global_Pthread;
        MMS->Mz_Global_Pthread[p]               = MMPt[p].Mz_Global_Pthread;
        MMS->MxV_Global_Pthread[p]              = MMPt[p].MxV_Global_Pthread;
        MMS->MyV_Global_Pthread[p]              = MMPt[p].MyV_Global_Pthread;
        MMS->MzV_Global_Pthread[p]              = MMPt[p].MzV_Global_Pthread;

    }
    ///////////////////////// Run Parameter /////////////////////////
    DRMP.CST = CST;
    DRMP.EV = EV;
    DRMP.M_Local = M_Local;
    DRMP.MMPt = MMPt;
    DRMP.Signal = Signal;
    
    ///////////////////////// Pthreads Part /////////////////////////
    int err;
    pthread_t *tid;
    tid = (pthread_t *)malloc(sizeof(pthread_t) * (Threads_Number + 2));
    DRMP.tid = tid;
    for(p=0;p<Threads_Number;++p)
    {
        err=pthread_create(&tid[p],NULL,Modify_Pthread_MC_Transfer,&MMPt[p]);
        if(err != 0)
        {
            printf("create thread error\n");
            Sleep(1000);
        }
    }
    err=pthread_create(&tid[Threads_Number],NULL,Modify_Pthread_MC_Synchronize_Transfer,MMS);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    ///////////////////////// OPGL Part /////////////////////////
    OPGL_GUI.State_Syn = StateSyn_Run;
    OPGL_MC_Pthread_Parameters *OMPP;
    OMPP = &DRMP.OMPP;
    OMPP->Mx_OPGL = Mx_OPGL; 
    OMPP->My_OPGL = My_OPGL;
    OMPP->Mz_OPGL = Mz_OPGL;
    OMPP->MxV_OPGL = MxV_OPGL; 
    OMPP->MyV_OPGL = MyV_OPGL;
    OMPP->MzV_OPGL = MzV_OPGL;
    OMPP->IPN = IPN;
    OMPP->CST = Make4DArray(7,m,n,l);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                for (r = 0; r < 7; ++r)
                {
                    OMPP->CST[r][i][j][k] = CST[r][i][j][k];
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
            err = pthread_create(&tid[Threads_Number + 1], NULL, OPGL_MC_Pthread_Transfer, OMPP);
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
    err = pthread_create(&tid_Part_Free, NULL, DiscreatPoint_Run_MC_Part_Free_Transfer, IPN);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    
    return 1;
}

int DiscreatPoint_Run_MC_Part_Free(Input_Parameter_NMD *IPN)
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
        pthread_join(DRMP.tid[p], NULL);
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
                Mx_Global[i][j][k]=DRMP.MMS.Mx_Global_Synchronize[i][j][k];
                My_Global[i][j][k]=DRMP.MMS.My_Global_Synchronize[i][j][k];
                Mz_Global[i][j][k]=DRMP.MMS.Mz_Global_Synchronize[i][j][k];
            }
        }
    }

    for(v=0;v<DRMP.BS.Virtual_Points;++v)
    {
        MxV_Global[v] = DRMP.MMS.MxV_Global_Synchronize[v];
        MyV_Global[v] = DRMP.MMS.MyV_Global_Synchronize[v];
        MzV_Global[v] = DRMP.MMS.MzV_Global_Synchronize[v];
    }
    ///////////////////////// Result Part /////////////////////////
    W=0.; SKN=0.;
    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        W   = W + DRMP.MMPt[p].W;
        SKN = SKN + DRMP.MMPt[p].SKN;
    }
    W = W / (DRMP.BS.Boundary_Points + DRMP.BS.Inner_Points);
    IPN->W = W; IPN->SKN = SKN;

    Energy_Variables_Distribution EVD;
    EVD.CST = DRMP.CST;    EVD.dxy = IPN->dxy;
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
    EVD.Local_Points_1D = DRMP.CPP.Local_Points_1D;
    EVD.Position_3DTo1D = DRMP.CPP.Position_3DTo1D;
    EVD.Energy_Coefficient = IPN->Energy_Coefficient;
    Energy_Distribution(&EVD);
    Write_NMD_Energy(IPN,EVD.w);
    Write_NMD_SKN_D(IPN,EVD.SKN_D);
    IPN->W   = EVD.W_Average;
    IPN->SKN = EVD.SKN;
    IPN->SKN_Positif = EVD.SKN_Positif;
    IPN->SKN_Area = EVD.SKN_Area;
    printf("%lf\t%lf\t%lf\t%lf\n",IPN->W,IPN->SKN,IPN->SKN_Positif,IPN->SKN_Area);

    ///////////////////////// Free Part /////////////////////////

    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        free3DArray(DRMP.MMPt[p].Mx_Global_Pthread,m,n);
        free3DArray(DRMP.MMPt[p].My_Global_Pthread,m,n);
        free3DArray(DRMP.MMPt[p].Mz_Global_Pthread,m,n);
        free(DRMP.MMPt[p].MxV_Global_Pthread);
        free(DRMP.MMPt[p].MyV_Global_Pthread);
        free(DRMP.MMPt[p].MzV_Global_Pthread);
        free(DRMP.MMPt[p].Signal);
        free(DRMP.MMPt[p].Cv);
        free(DRMP.MMPt[p].T);

        CalculatePointsThreads = DRMP.MMPt[p].CalculatePointsThreads_Upper - DRMP.MMPt[p].CalculatePointsThreads_Lower;
        for(q=0;q<CalculatePointsThreads;++q)
        {
            free(DRMP.EV[p][q].dxy);
            free(DRMP.EV[p][q].M0);
            free(DRMP.EV[p][q].LocalEnergy);
            free2DArray(DRMP.EV[p][q].CST,7);
            free2DArrayPointer(DRMP.EV[p][q].M_Local,4);
            free(DRMP.EV[p][q].Energy_Coefficient);
        }

        for (q = 1; q < 4; ++q)
        {
            pthread_mutex_destroy(&Signal_Synchronize_MC_Mutex[p][q]);
            pthread_cond_destroy(&Signal_Synchronize_MC_Cond[p][q]);
        }

        free(DRMP.MMPt[p].MMP);
        free(DRMP.EV[p]);
        free(Signal_Synchronize_MC_Mutex[p]);
        free(DRMP.MMPt[p].Sort_Calculate_Points);
    }
    free4DArray(DRMP.CST,7,m,n);
    free2DArrayinteger(DRMP.MMS.Bound_CalculatePoint_Threads,IPN->Threads_Number);
    free3DArray(DRMP.MMS.Mx_Global_Synchronize,m,n);
    free3DArray(DRMP.MMS.My_Global_Synchronize,m,n);
    free3DArray(DRMP.MMS.Mz_Global_Synchronize,m,n);
    free(DRMP.MMS.MxV_Global_Synchronize);
    free(DRMP.MMS.MyV_Global_Synchronize);
    free(DRMP.MMS.MzV_Global_Synchronize);
    free(DRMP.MMS.Mx_Global_Pthread);
    free(DRMP.MMS.My_Global_Pthread);
    free(DRMP.MMS.Mz_Global_Pthread);
    free(DRMP.MMS.MxV_Global_Pthread);
    free(DRMP.MMS.MyV_Global_Pthread);
    free(DRMP.MMS.MzV_Global_Pthread);
    free(DRMP.MMS.Sort_Calculate_Points);
    free3DArray(DRMP.M_Local,IPN->Threads_Number,4);
    free3DArrayinteger(DRMP.CPP.Local_Points_1D,DRMP.BS.Calculate_Points,7);
    free3DArrayinteger(DRMP.CPP.Position_3DTo1D,m,n);
    free(DRMP.CPP.Sort_Calculate_Points);
    free(DRMP.MMPt);
    free(DRMP.EV);
    free2DArrayintegerPointer(DRMP.Signal,IPN->Threads_Number+1);
    free(Signal_Synchronize_MC_Mutex);
    free(M_Synchronize_MC_Mutex);
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free4DArray(DRMP.OMPP.CST,7,m,n);
    free(DRMP.tid);

    ///////////////////////// Result Part /////////////////////////
    end_time  = clock();
    DRMP.end_time = end_time;
    duration  = (double)(end_time-DRMP.start_time)/CLOCKS_PER_SEC;
    IPN->Time = duration;

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Finish;
    pthread_cond_signal(&GPM.Run_GUI_Cond);
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);

    printf("Free Done\n");
    return 1;
}

void* DiscreatPoint_Run_MC_Part_Free_Transfer(void *arg)
{
    Input_Parameter_NMD *IPN;
    IPN=(Input_Parameter_NMD*)arg;
    DiscreatPoint_Run_MC_Part_Free(IPN);
    return NULL;
}

int DiscreatPoint_Run_MC(Input_Parameter_NMD *IPN)
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
            DiscreatPoint_Run_MC_Part_Run(IPN);
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