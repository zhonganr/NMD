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
#include <GNEB.h>
#include <malloc.h>
#include <OPGL.h>
#include <Initial_State.h>

extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
extern int N_Energy_Terms;
pthread_mutex_t **Signal_Synchronize_GNEB_Mutex, *M_Synchronize_GNEB_Mutex;
pthread_cond_t **Signal_Synchronize_GNEB_Cond;

extern pthread_mutex_t OPGL_Mutex;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
double ***Mx_Initial, ***My_Initial, ***Mz_Initial, *MxV_Initial, *MyV_Initial, *MzV_Initial, 
       ***Mx_Final, ***My_Final, ***Mz_Final, *MxV_Final, *MyV_Final, *MzV_Final,
       ****Mx_Image, ****My_Image, ****Mz_Image, **MxV_Image, **MyV_Image, **MzV_Image;
double Energy_Variation_Image;
extern GUI_Parameter GPM;
DiscreatPoint_Run_GNEB_Parameter DRGP;

using namespace std;

extern double l1_Strain, l2_Strain, l3_Strain, q_Strain, lo1_Strain, lo2_Strain, lo3_Strain, d_epsilon,
              l1_Stress, l2_Stress, l3_Stress, q_Stress, lo1_Stress, lo2_Stress, lo3_Stress, k_eff, k_0;

int Case_Energy(double E1, double E2, double E3)
{
    if(E1>=E2 && E2>=E3)
    {
        return 1;
    }else if(E1<=E2 && E2<=E3)
    {
        return 2;
    }else
    {
        return 3;
    }
    
}

int Cross_Product(double *V1, double *V2, double *V_Result)
{
    //V1 x V2 = V_Result
    V_Result[0] = V1[1] * V2[2] - V1[2] * V2[1];
    V_Result[1] = V1[2] * V2[0] - V1[0] * V2[2];
    V_Result[2] = V1[0] * V2[1] - V1[1] * V2[0];

    return 1;
}


int Modify_GNEB_OneStep(Modify_GNEB_Parameters *MGP)
{
    //Assignment
    double **t_Image, *W_Image, **F_Image, k_Spring, **Heff_Image, **D_Image, *Heff_Local, *Heff_Initial, *Heff_Final, dT;
    double *M_Initial, *M_Final, **M_Image, W_Initial, W_Final, Module, *Temp1_Array, *Temp2_Array, *Temp3_Array;
    int N_Image, i_Image, i, Case_W;

    dT = 0.0002;   k_Spring = 200;
    N_Image = MGP->N_Image;

    Heff_Local   = Make1DArray(3);
    Heff_Initial = Make1DArray(3);
    Heff_Final   = Make1DArray(3);
    M_Initial    = Make1DArray(3);
    M_Final      = Make1DArray(3);
    Temp1_Array  = Make1DArray(3);
    Temp2_Array  = Make1DArray(3);
    Temp3_Array  = Make1DArray(3);
    M_Image      = Make2DArray(N_Image,3);
    Heff_Image   = Make2DArray(N_Image, 3);
    F_Image      = Make2DArray(N_Image, 3);
    t_Image      = Make2DArray(N_Image, 3);
    D_Image      = Make2DArray(N_Image, 3);
    W_Image      = Make1DArray(N_Image);

    MGP->EV_Initial->M0[0] = *MGP->EV_Initial->M_Local[1][1];
    MGP->EV_Initial->M0[1] = *MGP->EV_Initial->M_Local[2][1];
    MGP->EV_Initial->M0[2] = *MGP->EV_Initial->M_Local[3][1];
    M_Initial[0]           = *MGP->EV_Initial->M_Local[1][1];
    M_Initial[1]           = *MGP->EV_Initial->M_Local[2][1];
    M_Initial[2]           = *MGP->EV_Initial->M_Local[3][1];

    MGP->EV_Final->M0[0] = *MGP->EV_Final->M_Local[1][1];
    MGP->EV_Final->M0[1] = *MGP->EV_Final->M_Local[2][1];
    MGP->EV_Final->M0[2] = *MGP->EV_Final->M_Local[3][1];
    M_Final[0]           = *MGP->EV_Final->M_Local[1][1];
    M_Final[1]           = *MGP->EV_Final->M_Local[2][1];
    M_Final[2]           = *MGP->EV_Final->M_Local[3][1];

    //printf("%0.8f\t%0.8f\t%0.8f\n", M_Initial[0], M_Initial[1], M_Initial[2]);

    Effective_Field_LocalCell_M0(MGP->EV_Initial, Temp3_Array);
    for (i = 0; i < 3; ++i)
    {
        Heff_Initial[i] = Temp3_Array[i];
    }

    Effective_Field_LocalCell_M0(MGP->EV_Final, Temp3_Array);
    for (i = 0; i < 3; ++i)
    {
        Heff_Final[i] = Temp3_Array[i];
    }
    W_Initial = Energy_Local_M0(MGP->EV_Initial);
    W_Final   = Energy_Local_M0(MGP->EV_Final);
    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        MGP->EV_Image[i_Image].M0[0] = *MGP->EV_Image[i_Image].M_Local[1][1];
        MGP->EV_Image[i_Image].M0[1] = *MGP->EV_Image[i_Image].M_Local[2][1];
        MGP->EV_Image[i_Image].M0[2] = *MGP->EV_Image[i_Image].M_Local[3][1];
        M_Image[i_Image][0]          = *MGP->EV_Image[i_Image].M_Local[1][1];
        M_Image[i_Image][1]          = *MGP->EV_Image[i_Image].M_Local[2][1];
        M_Image[i_Image][2]          = *MGP->EV_Image[i_Image].M_Local[3][1];

        //printf("%d\t%0.8f\t%0.8f\t%0.8f\n", i_Image,M_Initial[0], M_Initial[1], M_Initial[2]);
        //printf("%d\t%0.8f\t%0.8f\t%0.8f\n", i_Image,M_Image[i_Image][0], M_Image[i_Image][1], M_Image[i_Image][2]);
        //printf("%d\t%0.8f\t%0.8f\t%0.8f\n", i_Image,M_Final[0], M_Final[1], M_Final[2]);

        Effective_Field_LocalCell_M0(&(MGP->EV_Image[i_Image]), Heff_Local);
        W_Image[i_Image] = Energy_Local_M0(&(MGP->EV_Image[i_Image]));
        for (i = 0; i < 3; ++i)
        {
            Heff_Image[i_Image][i] = Heff_Local[i];
        }
    }

    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        if (i_Image == 0)
        {
            //////////////////////////// t /////////////////////////////////////////
            Case_W = Case_Energy(W_Initial, W_Image[i_Image], W_Image[i_Image + 1]);
            if(Case_W==1)
            {
                Module = sqrt(pow(M_Image[i_Image][0] - M_Initial[0], 2) + pow(M_Image[i_Image][1] - M_Initial[1], 2) + pow(M_Image[i_Image][2] - M_Initial[2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image][0] - M_Initial[0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image][1] - M_Initial[1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image][2] - M_Initial[2])/Module;
            }else if(Case_W==2)
            {
                Module = sqrt(pow(M_Image[i_Image+1][0] - M_Image[i_Image][0], 2) + pow(M_Image[i_Image+1][1] - M_Image[i_Image][1], 2) + pow(M_Image[i_Image+1][2] - M_Image[i_Image][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image+1][0] - M_Image[i_Image][0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image+1][1] - M_Image[i_Image][1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image+1][2] - M_Image[i_Image][2])/Module;
            }else
            {
                Module = sqrt(pow(M_Image[i_Image+1][0] - M_Initial[0], 2) + pow(M_Image[i_Image+1][1] - M_Initial[1], 2) + pow(M_Image[i_Image+1][2] - M_Initial[2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image+1][0] - M_Initial[0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image+1][1] - M_Initial[1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image+1][2] - M_Initial[2])/Module;
            }
            /////////////////////////////////// F //////////////////////////////////////////////////
            Module = sqrt(pow(M_Image[i_Image+1][0] - M_Image[i_Image][0], 2) + pow(M_Image[i_Image+1][1] - M_Image[i_Image][1], 2) + pow(M_Image[i_Image+1][2] - M_Image[i_Image][2], 2))-
                     sqrt(pow(M_Image[i_Image][0] - M_Initial[0], 2) + pow(M_Image[i_Image][1] - M_Initial[1], 2) + pow(M_Image[i_Image][2] - M_Initial[2], 2));

            F_Image[i_Image][0] = k_Spring * Module * t_Image[i_Image][0];
            F_Image[i_Image][1] = k_Spring * Module * t_Image[i_Image][1];
            F_Image[i_Image][2] = k_Spring * Module * t_Image[i_Image][2];
        }else if(i_Image==N_Image-1)
        {
            //////////////////////////// t /////////////////////////////////////////
            Case_W = Case_Energy(W_Image[i_Image-1], W_Image[i_Image], W_Final);
            if(Case_W==1)
            {
                Module = sqrt(pow(M_Image[i_Image][0] - M_Image[i_Image-1][0], 2) + pow(M_Image[i_Image][1] - M_Image[i_Image-1][1], 2) + pow(M_Image[i_Image][2] - M_Image[i_Image-1][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image][0] - M_Image[i_Image-1][0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image][1] - M_Image[i_Image-1][1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image][2] - M_Image[i_Image-1][2])/Module;
            }else if(Case_W==2)
            {
                Module = sqrt(pow(M_Final[0] - M_Image[i_Image][0], 2) + pow(M_Final[1] - M_Image[i_Image][1], 2) + pow(M_Final[2] - M_Image[i_Image][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Final[0] - M_Image[i_Image][0])/Module;
                t_Image[i_Image][1] = (M_Final[1] - M_Image[i_Image][1])/Module;
                t_Image[i_Image][2] = (M_Final[2] - M_Image[i_Image][2])/Module;
            }else
            {
                Module = sqrt(pow(M_Final[0] - M_Image[i_Image-1][0], 2) + pow(M_Final[1] - M_Image[i_Image-1][1], 2) + pow(M_Final[2] - M_Image[i_Image-1][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Final[0] - M_Image[i_Image-1][0])/Module;
                t_Image[i_Image][1] = (M_Final[1] - M_Image[i_Image-1][1])/Module;
                t_Image[i_Image][2] = (M_Final[2] - M_Image[i_Image-1][2])/Module;
            }
            /////////////////////////////////// F //////////////////////////////////////////////////
            Module = sqrt(pow(M_Final[0] - M_Image[i_Image][0], 2) + pow(M_Final[1] - M_Image[i_Image][1], 2) + pow(M_Final[2] - M_Image[i_Image][2], 2))-
                     sqrt(pow(M_Image[i_Image][0] - M_Image[i_Image-1][0], 2) + pow(M_Image[i_Image][1] - M_Image[i_Image-1][1], 2) + pow(M_Image[i_Image][2] - M_Image[i_Image-1][2], 2));

            F_Image[i_Image][0] = k_Spring * Module * t_Image[i_Image][0];
            F_Image[i_Image][1] = k_Spring * Module * t_Image[i_Image][1];
            F_Image[i_Image][2] = k_Spring * Module * t_Image[i_Image][2];
        }else
        {
            //////////////////////////// t /////////////////////////////////////////
            Case_W = Case_Energy(W_Image[i_Image-1], W_Image[i_Image], W_Image[i_Image+1]);
            if(Case_W==1)
            {
                Module = sqrt(pow(M_Image[i_Image][0] - M_Image[i_Image-1][0], 2) + pow(M_Image[i_Image][1] - M_Image[i_Image-1][1], 2) + pow(M_Image[i_Image][2] - M_Image[i_Image-1][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image][0] - M_Image[i_Image-1][0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image][1] - M_Image[i_Image-1][1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image][2] - M_Image[i_Image-1][2])/Module;
            }else if(Case_W==2)
            {
                Module = sqrt(pow(M_Image[i_Image+1][0] - M_Image[i_Image][0], 2) + pow(M_Image[i_Image+1][1] - M_Image[i_Image][1], 2) + pow(M_Image[i_Image+1][2] - M_Image[i_Image][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image+1][0] - M_Image[i_Image][0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image+1][1] - M_Image[i_Image][1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image+1][2] - M_Image[i_Image][2])/Module;
            }else
            {
                Module = sqrt(pow(M_Image[i_Image+1][0] - M_Image[i_Image-1][0], 2) + pow(M_Image[i_Image+1][1] - M_Image[i_Image-1][1], 2) + pow(M_Image[i_Image+1][2] - M_Image[i_Image-1][2], 2));
                if(Module==0)
                {
                    Module = 1;
                }
                t_Image[i_Image][0] = (M_Image[i_Image+1][0] - M_Image[i_Image-1][0])/Module;
                t_Image[i_Image][1] = (M_Image[i_Image+1][1] - M_Image[i_Image-1][1])/Module;
                t_Image[i_Image][2] = (M_Image[i_Image+1][2] - M_Image[i_Image-1][2])/Module;
            }
            /////////////////////////////////// F //////////////////////////////////////////////////
            Module = sqrt(pow(M_Image[i_Image+1][0] - M_Image[i_Image][0], 2) + pow(M_Image[i_Image+1][1] - M_Image[i_Image][1], 2) + pow(M_Image[i_Image+1][2] - M_Image[i_Image][2], 2))-
                     sqrt(pow(M_Image[i_Image][0] - M_Image[i_Image-1][0], 2) + pow(M_Image[i_Image][1] - M_Image[i_Image-1][1], 2) + pow(M_Image[i_Image][2] - M_Image[i_Image-1][2], 2));

            F_Image[i_Image][0] = k_Spring * Module * t_Image[i_Image][0];
            F_Image[i_Image][1] = k_Spring * Module * t_Image[i_Image][1];
            F_Image[i_Image][2] = k_Spring * Module * t_Image[i_Image][2];
        }
        
        /////////////////////////////////// D //////////////////////////////////////////////////
        Module = Heff_Image[i_Image][0] * t_Image[i_Image][0] + Heff_Image[i_Image][1] * t_Image[i_Image][1] + Heff_Image[i_Image][2] * t_Image[i_Image][2];
        D_Image[i_Image][0] = Heff_Image[i_Image][0] - Module * t_Image[i_Image][0] + F_Image[i_Image][0];
        D_Image[i_Image][1] = Heff_Image[i_Image][1] - Module * t_Image[i_Image][1] + F_Image[i_Image][1];
        D_Image[i_Image][2] = Heff_Image[i_Image][2] - Module * t_Image[i_Image][2] + F_Image[i_Image][2];

        //printf("%0.8f\t%0.8f\t%0.8f\n", D_Image[i_Image][0], D_Image[i_Image][1], D_Image[i_Image][2]);

        ///////////////////////////////// M //////////////////////////////////////////////////
        Cross_Product(MGP->EV_Image[i_Image].M0, D_Image[i_Image], Temp1_Array);
        Cross_Product(MGP->EV_Image[i_Image].M0, Temp1_Array, Temp2_Array);
        //printf("%0.8f\t%0.8f\t%0.8f\n", D_Image[i_Image][0], D_Image[i_Image][1], D_Image[i_Image][2]);
        if(fabs(Temp2_Array[0])>1000 || fabs(Temp2_Array[1])>1000 || fabs(Temp2_Array[2])>1000)
        {
            printf("%d\n", i_Image);
            printf("%0.8f\t%0.8f\t%0.8f\n", Temp2_Array[0], Temp2_Array[1], Temp2_Array[2]);
            printf("%0.8f\t%0.8f\t%0.8f\n", D_Image[i_Image][0], D_Image[i_Image][1], D_Image[i_Image][2]);
            printf("%0.8f\t%0.8f\t%0.8f\n", F_Image[i_Image][0], F_Image[i_Image][1], F_Image[i_Image][2]);
            printf("%0.8f\t%0.8f\t%0.8f\n", t_Image[i_Image][0], t_Image[i_Image][1], t_Image[i_Image][2]);
            printf("%0.8f\t%0.8f\t%0.8f\n", Heff_Image[i_Image][0], Heff_Image[i_Image][1], Heff_Image[i_Image][2]);
        }

        Module = sqrt(pow(M_Image[i_Image][0], 2) + pow(M_Image[i_Image][1], 2) + pow(M_Image[i_Image][2], 2));
        MGP->EV_Image[i_Image].M0[0] = MGP->EV_Image[i_Image].M0[0] - dT * Temp2_Array[0] / Module;
        MGP->EV_Image[i_Image].M0[1] = MGP->EV_Image[i_Image].M0[1] - dT * Temp2_Array[1] / Module;
        MGP->EV_Image[i_Image].M0[2] = MGP->EV_Image[i_Image].M0[2] - dT * Temp2_Array[2] / Module;
        //MGP->EV_Image[i_Image].M0[0] = MGP->EV_Image[i_Image].M0[0] + dT * D_Image[i_Image][0] / Module;
        //MGP->EV_Image[i_Image].M0[1] = MGP->EV_Image[i_Image].M0[1] + dT * D_Image[i_Image][1] / Module;
        //MGP->EV_Image[i_Image].M0[2] = MGP->EV_Image[i_Image].M0[2] + dT * D_Image[i_Image][2] / Module;

    }

    ////////////////////////////////////// free ////////////////////////////////////////////
    free(Heff_Local);
    free(Heff_Initial);
    free(Heff_Final);
    free(M_Initial);
    free(M_Final);
    free(Temp1_Array);
    free(Temp2_Array);
    free(Temp3_Array);
    free2DArray(M_Image, N_Image);
    free2DArray(Heff_Image, N_Image);
    free2DArray(F_Image, N_Image);
    free2DArray(t_Image, N_Image);
    free2DArray(D_Image, N_Image);
    free(W_Image);
    return 1;
}


int Modify_Pthread_GNEB(Modify_GNEB_Pthreads *MGPt)
{
    int Pthread_Index, NCycles, *Signal;
    int Signal_Synchronize, N_Image;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int i, j, Lower, i_Image;
    int Signal_RunAndStop;

    Pthread_Index  = MGPt->Pthread_Index;
    //Pthread_Number = MGPt->Pthread_Number;
    NCycles        = MGPt->NCycles;
    Signal         = MGPt->Signal;
    N_Image        = MGPt->N_Image;

    CalculatePointsThreads_Lower = MGPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MGPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;
    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1]=1;
    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);//1 is Pthread State.

    if(CalculatePointsThreads_Lower==0)
    {
        Lower=1;
    }else
    {
        Lower=0;
    }
////////////////////////////////////////////////////////////////////////////////////////////////////
    for (i = 0; i < NCycles + 1; ++i)
    {
        for (i_Image = 0; i_Image < N_Image;++i_Image)
        {
            MGPt->PTPS->Energy_Image[i_Image] = 0;
            MGPt->PTPS->SKN_Image[i_Image] = 0;
            //printf("1\n");
        }
        pthread_mutex_lock(&M_Synchronize_GNEB_Mutex[Pthread_Index]);
        for(j=Lower;j<CalculatePointsThreads;++j)
        {
            Modify_GNEB_OneStep(&MGPt->MGP[j]);
            for (i_Image = 0; i_Image < N_Image;++i_Image)
            {
                MGPt->PTPS->Energy_Image[i_Image] = MGPt->PTPS->Energy_Image[i_Image] + Energy_Single(&MGPt->MGP[j].EV_Image[i_Image]);
                MGPt->PTPS->SKN_Image[i_Image]    = MGPt->PTPS->SKN_Image[i_Image] + Skyrmion_Number_Single(&MGPt->MGP[j].EV_Image[i_Image]);
            }
        }
        
        pthread_mutex_unlock(&M_Synchronize_GNEB_Mutex[Pthread_Index]);

        pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
        Signal_Synchronize = Signal[3];
        pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.
        if (Signal_Synchronize != SignalSyn_Stop)
        {
            ////////////////////////// Prepare to Synchronize ///////////////////////////
            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]); //1 is Pthread State.
            Signal[1] = StatePthread_PrepareToSynchronize;
            //printf("%d\tStatePthread_PrepareToSynchronize\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);//1 is Pthread State.

            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_Synchronize && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_GNEB_Cond[Pthread_Index][3], &Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);
            }
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            //printf("%d\tSignalPthread_Stop\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }
        
        if (Signal_Synchronize != SignalSyn_Stop)
        {
            ////////////////////////// Synchronize ///////////////////////////
            pthread_mutex_lock(&M_Synchronize_GNEB_Mutex[Pthread_Index]);
            for(j=Lower;j<CalculatePointsThreads;++j)
            {
                for (i_Image = 0; i_Image < N_Image; ++i_Image)
                {
                    *MGPt->MGP[j].EV_Image[i_Image].M_Local[1][1] = MGPt->MGP[j].EV_Image[i_Image].M0[0];
                    *MGPt->MGP[j].EV_Image[i_Image].M_Local[2][1] = MGPt->MGP[j].EV_Image[i_Image].M0[1];
                    *MGPt->MGP[j].EV_Image[i_Image].M_Local[3][1] = MGPt->MGP[j].EV_Image[i_Image].M0[2];
                }
                /*
                for (i_Image = 0; i_Image < N_Image;++i_Image)
                {
                    MGPt->PTPS->Energy_Image[i_Image] = MGPt->PTPS->Energy_Image[i_Image] + Energy_Single(&MGPt->MGP[j].EV_Image[i_Image]);
                    MGPt->PTPS->SKN_Image[i_Image]    = MGPt->PTPS->SKN_Image[i_Image] + Skyrmion_Number_Single(&MGPt->MGP[j].EV_Image[i_Image]);
                }
                */
            }
            pthread_mutex_unlock(&M_Synchronize_GNEB_Mutex[Pthread_Index]);

            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_PrepareToRun;
            //printf("%d\tStatePthread_PrepareToRun\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);//1 is Pthread State.

            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_PrepareToSynchronize && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_GNEB_Cond[Pthread_Index][3], &Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);
            }
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_Run;
            //printf("%d\tStatePthread_Run\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);//1 is Pthread State.
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            //printf("%d\tSignalPthread_Stop\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
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

    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Stop;
    //printf("%d\tStatePthread_Stop\n", Pthread_Index);
    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[Pthread_Index][1]);//1 is Pthread State.

    printf("fin\t%d\n",Pthread_Index);
    return 1;
}

void* Modify_Pthread_GNEB_Transfer(void *arg)
{
    Modify_GNEB_Pthreads *MGPt;
    MGPt=(Modify_GNEB_Pthreads*)arg;
    Modify_Pthread_GNEB(MGPt);
    return NULL;
}

int Modify_Pthread_GNEB_Synchronize(Modify_GNEB_Synchronize *MGS)
{
    int p, Synchronize_Times;
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l;
    int Pthread_Number, ***Signal;
    int **Bound_CalculatePoint_Threads;
    double ****CST;
    double *W_Image, *SKN_Image;
    int N_Total, N_Image, i_Image, Pthread_cond_State;
    int Signal_RunAndStop;

    clock_t start_time_P, end_time_P;
    double duration;

    m=MGS->m;   n=MGS->n;   l=MGS->l;
    Pthread_Number               = MGS->Pthread_Number;
    Bound_CalculatePoint_Threads = MGS->Bound_CalculatePoint_Threads;
    Local_Points_1D              = MGS->Local_Points_1D;
    Position_3DTo1D              = MGS->Position_3DTo1D;
    Signal                       = MGS->Signal;
    CST                          = MGS->CST;
    N_Image                      = MGS->N_Image;
    N_Total                      = MGS->Inner_Points + MGS->Boundary_Points;

    ///////////////////////////////////////////////////////////////////////////


    W_Image   = Make1DArray(N_Image);
    SKN_Image = Make1DArray(N_Image);

    Energy_Variables_Distribution EVD;
    EVD.CST = CST;    EVD.dxy = MGS->dxy;
    EVD.h   = MGS->h; EVD.t   = MGS->t;
    EVD.m   = m;      EVD.n   = n;
    EVD.l   = l;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D = Local_Points_1D;
    EVD.Position_3DTo1D = Position_3DTo1D;
    EVD.Energy_Coefficient = MGS->Energy_Coefficient;

    //FILE *fp_Mx, *fp_My, *fp_Mz, *fp_W, *fp_SKN_D;
    FILE *fp_Process;
    //char filenameMx[500], filenameMy[500], filenameMz[500], filenamew[500], filenameSKN_D[500];
    char filenameP[500];

    sprintf(filenameP, "Version %d\\Process.dat",MGS->Version);

    fp_Process = fopen(filenameP, "w+");
    fprintf(fp_Process,"Times\t\t Image\t\t W\t\t SKN\t\t Time\n");
    
    ////////////////////////////////////////////////////////////////////////////
    int Pthread_State, Synchronize_State;
    Synchronize_Times=0;

    for(p=0;p<Pthread_Number;++p)
    {
        pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
        *Signal[p + 1][3] = SignalSyn_PrepareToSynchronize;                        
        pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
    }
    
    do
    {
        start_time_P = clock();
        ///////////////////////// Make sure all child threads are waiting ///////////////////////////////
        do
        {
            Synchronize_State = StateSyn_Synchronize;
            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][1]);  //1 is Pthread State.
                Pthread_State=*Signal[p+1][1];
                pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][1]);//1 is Pthread State.
                if (Pthread_State == StatePthread_Stop)
                {
                    Synchronize_State = StateSyn_Stop;
                    break;
                }else if(Pthread_State != StatePthread_PrepareToSynchronize)
                {
                    Synchronize_State = StateSyn_PrepareToSynchronize;
                }
            }
        } while (Synchronize_State == StateSyn_PrepareToSynchronize);

        ///////////////////// Calculate Energy ////////////////////////////////////
        for(i_Image=0;i_Image<N_Image;++i_Image)
        {
            W_Image[i_Image] = 0;
            SKN_Image[i_Image]=0;

            for (p = 0; p < Pthread_Number; ++p)
            {
                W_Image[i_Image]   = W_Image[i_Image] + MGS->PTPS[p].Energy_Image[i_Image] / N_Total;
                SKN_Image[i_Image] = SKN_Image[i_Image] + MGS->PTPS[p].SKN_Image[i_Image];
            }
        }

        ///////////////////////// Stop child threads ///////////////////////////////
        if (Synchronize_State == StateSyn_Stop)
        {
            for(p=0;p<Pthread_Number;++p)
            {
                Sleep(100);
                pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
                pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State=*Signal[p+1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun && Pthread_State != StatePthread_PrepareToSynchronize)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_GNEB_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }
        }
        /////////////////////// Synchronize /////////////////////////////////////////
        if (Synchronize_State == StateSyn_Synchronize)
        {
            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Synchronize;
                pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToSynchronize)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]); //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_GNEB_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }

            do
            {
                Synchronize_State = StateSyn_PrepareToRun;
                for(p=0;p<Pthread_Number;++p)
                {
                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State=*Signal[p+1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun)
                    {
                        Synchronize_State = StateSyn_Synchronize;
                    }
                }
            } while (Synchronize_State == StateSyn_Synchronize);

            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_GNEB_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }
        }

        ///////////////////// Calculate Energy ////////////////////////////////////
        /*
        for(i_Image=0;i_Image<N_Image;++i_Image)
        {
            W_Image[i_Image] = 0;
            SKN_Image[i_Image]=0;

            for (p = 0; p < Pthread_Number; ++p)
            {
                W_Image[i_Image]   = W_Image[i_Image] + MGS->PTPS[p].Energy_Image[i_Image] / N_Total;
                SKN_Image[i_Image] = SKN_Image[i_Image] + MGS->PTPS[p].SKN_Image[i_Image];
            }
        }
        */

        end_time_P = clock();
        duration = (double)(end_time_P-start_time_P)/CLOCKS_PER_SEC;
        for (i_Image = 0; i_Image < N_Image;++i_Image)
        {
            fprintf(fp_Process, "%d\t\t %d\t\t %0.8f\t %0.8f\t %0.8f\n",
                    Synchronize_Times, i_Image, W_Image[i_Image], SKN_Image[i_Image], duration);
            printf("%d\t\t %d\t\t %0.8f\t %0.8f\t %0.8f\n",
                   Synchronize_Times, i_Image, W_Image[i_Image], SKN_Image[i_Image], duration);
        }
        Synchronize_Times=Synchronize_Times+1;

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
        pthread_mutex_lock(&Signal_Synchronize_GNEB_Mutex[p][3]);  //3 is Synchronize Signal.
        *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
        pthread_cond_signal(&Signal_Synchronize_GNEB_Cond[p][3]);
        pthread_mutex_unlock(&Signal_Synchronize_GNEB_Mutex[p][3]);//3 is Synchronize Signal.
    }

    fclose(fp_Process);

    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free(W_Image);
    free(SKN_Image);

    return 1;
}

void* Modify_Pthread_GNEB_Synchronize_Transfer(void *arg)
{
    Modify_GNEB_Synchronize *MGS;
    MGS=(Modify_GNEB_Synchronize*)arg;
    Modify_Pthread_GNEB_Synchronize(MGS);
    return NULL;
}

int DiscreatPoint_Run_GNEB_Part_Run(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, Version, Threads_Number, Calculate_Points, CalculatePointsThreads;
    double h, t, *dxy;
    double ***M_Local, ****CST, *Energy_Image, *SKN_Image;
    int i_G, j_G, k_G, v, r, s, i, j, k, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, i_Image;
    int ****Local_Points_1D_LocalEnergy, **Neighbour_Points_3Dto1D_Array, Point_Index;
    int Boundary_Type, ****LocalEnergy;

    int N_Image = IPN->N_Image;

    clock_t start_time;
    start_time = clock();
    DRGP.start_time = start_time;
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

    Energy_Image = Make1DArray(N_Image + 2);     SKN_Image    = Make1DArray(N_Image + 2);
    
    Mx_Initial = Make3DArray(m, n, l);            My_Initial = Make3DArray(m, n, l);            Mz_Initial = Make3DArray(m, n, l);
    Mx_Final   = Make3DArray(m, n, l);            My_Final   = Make3DArray(m, n, l);            Mz_Final   = Make3DArray(m, n, l);
    Mx_Image   = Make4DArray(N_Image, m, n, l);   My_Image   = Make4DArray(N_Image, m, n, l);   Mz_Image   = Make4DArray(N_Image, m, n, l);

    ///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape *BS;
    Calculate_Points_Parameters *CPP;

    BS = &DRGP.BS;
    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    Calculate_Points=BS->Calculate_Points;
    LocalEnergy     =BS->LocalEnergy;
    //printf("%d\n",Calculate_Points);

    CPP = &DRGP.CPP;
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
////////////////////////////////////////////////////////////////////////////
    MxV_Initial = Make1DArray(BS->Virtual_Points);
    MyV_Initial = Make1DArray(BS->Virtual_Points);
    MzV_Initial = Make1DArray(BS->Virtual_Points);

    MxV_Final = Make1DArray(BS->Virtual_Points);
    MyV_Final = Make1DArray(BS->Virtual_Points);
    MzV_Final = Make1DArray(BS->Virtual_Points);

    MxV_Image = Make2DArray(N_Image, BS->Virtual_Points);
    MyV_Image = Make2DArray(N_Image, BS->Virtual_Points);
    MzV_Image = Make2DArray(N_Image, BS->Virtual_Points);
    
    ////////////////////////////////////////////////////////////////////////////
    InitialM_Target(Mx_Final, My_Final, Mz_Final, MxV_Final, MyV_Final, MzV_Final, BS);
    //InitialM_3D(IPN, Mx_Initial, My_Initial, Mz_Initial, MxV_Initial, MyV_Initial, MzV_Initial, BS);
    //Magnetization_Transfer_By_Rotation(IPN, Mx_Initial, My_Initial, Mz_Initial, 0.1);
    
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                Mx_Initial[i][j][k] = Mx_Global[i][j][k];
                My_Initial[i][j][k] = My_Global[i][j][k];
                Mz_Initial[i][j][k] = Mz_Global[i][j][k];
                //printf("%0.8f\t%0.8f\t%0.8f\n", Mx_Initial[i][j][k], My_Initial[i][j][k], Mz_Initial[i][j][k]);
                //printf("%0.8f\t%0.8f\t%0.8f\n", Mx_Final[i][j][k], My_Final[i][j][k], Mz_Final[i][j][k]);
            }
        }
    }

    for (v = 0; v < BS->Virtual_Points; ++v)
    {
        MxV_Initial[v] = MxV_Global[v];
        MyV_Initial[v] = MyV_Global[v];
        MzV_Initial[v] = MzV_Global[v];
    }


    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        for (i = 0; i < m; ++i)
        {
            for (j = 0; j < n; ++j)
            {
                for (k = 0; k < l; ++k)
                {
                    Mx_Image[i_Image][i][j][k] = (N_Image - i_Image) / (N_Image + 1.) * Mx_Initial[i][j][k] + (i_Image + 1.) / (N_Image + 1.) * Mx_Final[i][j][k];
                    My_Image[i_Image][i][j][k] = (N_Image - i_Image) / (N_Image + 1.) * My_Initial[i][j][k] + (i_Image + 1.) / (N_Image + 1.) * My_Final[i][j][k];
                    Mz_Image[i_Image][i][j][k] = (N_Image - i_Image) / (N_Image + 1.) * Mz_Initial[i][j][k] + (i_Image + 1.) / (N_Image + 1.) * Mz_Final[i][j][k];

                }
            }
        }

        for (v = 0; v < BS->Virtual_Points; ++v)
        {
            MxV_Image[i_Image][v] = (N_Image - i_Image) / (N_Image + 1.) * MxV_Initial[v] + (i_Image + 1.) / (N_Image + 1.) * MxV_Final[v];
            MyV_Image[i_Image][v] = (N_Image - i_Image) / (N_Image + 1.) * MyV_Initial[v] + (i_Image + 1.) / (N_Image + 1.) * MyV_Final[v];
            MzV_Image[i_Image][v] = (N_Image - i_Image) / (N_Image + 1.) * MzV_Initial[v] + (i_Image + 1.) / (N_Image + 1.) * MzV_Final[v];
        }
    }

    //Read_Image(IPN, BS, Mx_Image, My_Image, Mz_Image, MxV_Image, MyV_Image, MzV_Image);
    Energy_Variables **EV_Initial, **EV_Final;

    Modify_GNEB_Pthreads *MGPt;

    MGPt=(Modify_GNEB_Pthreads*)malloc(sizeof(Modify_GNEB_Pthreads)*(Threads_Number));

    M_Local       = Make3DArray(Threads_Number,4,14);
    CST           = Make4DArray(7,m,n,l);

    Elastic_Field_Initial(CST,IPN);
    Signal_Synchronize_GNEB_Mutex = (pthread_mutex_t**)malloc(sizeof(pthread_mutex_t*)*(Threads_Number));
    M_Synchronize_GNEB_Mutex      = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(Threads_Number));
    Signal_Synchronize_GNEB_Cond  = (pthread_cond_t**)malloc(sizeof(pthread_cond_t*)*(Threads_Number));

    int ***Signal;
    double ***Energy_Change;

    Signal=Make2DArrayintegerPointer(Threads_Number+1,4);
    Energy_Change = Make2DArrayPointer(Threads_Number + 1, 3);
    EV_Initial = (Energy_Variables **)malloc(sizeof(Energy_Variables *) * (Threads_Number));
    EV_Final   = (Energy_Variables **)malloc(sizeof(Energy_Variables *) * (Threads_Number));
    int p, q, o, I, J;
///////////////////////// Pthreads Parameters Initial /////////////////////////
    for(p=0;p<Threads_Number;++p)
    {
        //MGPt[p].CST           =Make4DArray(7,m,n,l);
        //Elastic_Field_Initial(MGPt[p].CST,IPN);
        MGPt[p].Pthread_Index  = p;
        MGPt[p].NCycles        = IPN->NCycles;
        MGPt[p].Pthread_Number = Threads_Number;
        MGPt[p].N_Image        = N_Image;

        MGPt[p].Continuity_Modify_Mode = IPN->Continuity_Modify_Mode;
        MGPt[p].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;

        MGPt[p].Sort_Calculate_Points = Make1DArrayinteger(BS->Calculate_Points);
        for(i=0;i<BS->Calculate_Points;++i)
        {
            MGPt[p].Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
        }
        MGPt[p].Local_Points_1D = CPP->Local_Points_1D;
        MGPt[p].Position_3DTo1D = CPP->Position_3DTo1D;
        MGPt[p].h   = h;    MGPt[p].t   = t;
        MGPt[p].m   = m;    MGPt[p].n   = n;
        MGPt[p].l   = l;
        MGPt[p].dxy = dxy;  MGPt[p].CST = CST;


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
        MGPt[p].CalculatePointsThreads_Lower=CalculatePointsThreads_Lower;
        MGPt[p].CalculatePointsThreads_Upper=CalculatePointsThreads_Upper;

        MGPt[p].PTPS = &PTPS[p];
        PTPS[p].Energy_Image = Make1DArray(N_Image);
        PTPS[p].SKN_Image = Make1DArray(N_Image);

        /////////////////////////// Synchronize Parameter Initial ////////////////////////////////
        MGPt[p].Signal           = Make1DArrayinteger(4);
        MGPt[p].Signal[1]        = StatePthread_Run;//0 is Pthread State.
        MGPt[p].Signal[2]        = SignalPthread_Run;//1 is Pthread Signal.
        MGPt[p].Signal[3]        = SignalSyn_PrepareToSynchronize;//2 is Synchronize Signal.
        Signal[p+1][1]           = &MGPt[p].Signal[1];
        Signal[p+1][2]           = &MGPt[p].Signal[2];
        Signal[p+1][3]           = &MGPt[p].Signal[3];

        Signal_Synchronize_GNEB_Mutex[p] = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(4));
        Signal_Synchronize_GNEB_Cond[p]  = (pthread_cond_t*)malloc(sizeof(pthread_cond_t)*(4));
        pthread_mutex_init(&Signal_Synchronize_GNEB_Mutex[p][1],0);//1 is Pthread State.
        pthread_mutex_init(&Signal_Synchronize_GNEB_Mutex[p][2],0);//2 is Pthread Signal.
        pthread_mutex_init(&Signal_Synchronize_GNEB_Mutex[p][3],0);//3 is Synchronize Signal.
        pthread_cond_init(&Signal_Synchronize_GNEB_Cond[p][1],NULL);//1 is Pthread State.
        pthread_cond_init(&Signal_Synchronize_GNEB_Cond[p][2],NULL);//2 is Pthread Signal.
        pthread_cond_init(&Signal_Synchronize_GNEB_Cond[p][3],NULL);//3 is Synchronize Signal.
        pthread_mutex_init(&M_Synchronize_GNEB_Mutex[p],0);

        /////////////////////////// Discrete Point Parameter Initial //////////////////////////
        MGPt[p].MGP    = (Modify_GNEB_Parameters*)malloc(sizeof(Modify_GNEB_Parameters)*(CalculatePointsThreads));
        EV_Initial[p]  = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));
        EV_Final[p]    = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));

        for(r=0;r<14;++r)
        {
            M_Local[p][1][r]=0.;//Mx
            M_Local[p][2][r]=0.;//My
            M_Local[p][3][r]=0.;//Mz
        }

        for(q=0;q<CalculatePointsThreads;++q)
        {
            MGPt[p].MGP[q].N_Image    = N_Image;
            MGPt[p].MGP[q].EV_Initial = &EV_Initial[p][q];
            MGPt[p].MGP[q].EV_Final   = &EV_Final[p][q];
            MGPt[p].MGP[q].EV_Image = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(N_Image));
            for (i_Image = 0; i_Image < N_Image; ++i_Image)
            {
                MGPt[p].MGP[q].EV_Image[i_Image].dxy = Make1DArray(3);
                MGPt[p].MGP[q].EV_Image[i_Image].dxy[0] = dxy[0];    MGPt[p].MGP[q].EV_Image[i_Image].dxy[1] = dxy[1];
                MGPt[p].MGP[q].EV_Image[i_Image].dxy[2] = dxy[2];    MGPt[p].MGP[q].EV_Image[i_Image].i      = 3;
                MGPt[p].MGP[q].EV_Image[i_Image].j      = 3;         MGPt[p].MGP[q].EV_Image[i_Image].k      = 3;
                MGPt[p].MGP[q].EV_Image[i_Image].m      = m;         MGPt[p].MGP[q].EV_Image[i_Image].n      = n;
                MGPt[p].MGP[q].EV_Image[i_Image].l      = l;
                MGPt[p].MGP[q].EV_Image[i_Image].h      = h;         MGPt[p].MGP[q].EV_Image[i_Image].t      = t;
                MGPt[p].MGP[q].EV_Image[i_Image].M0          = Make1DArray(3);
                MGPt[p].MGP[q].EV_Image[i_Image].LocalEnergy = Make1DArrayinteger(7);
                MGPt[p].MGP[q].EV_Image[i_Image].CST         = Make2DArray(7,7);
                MGPt[p].MGP[q].EV_Image[i_Image].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
                MGPt[p].MGP[q].EV_Image[i_Image].M_Local     = Make2DArrayPointer(4,15);
                MGPt[p].MGP[q].EV_Image[i_Image].Energy_Coefficient = Make1DArray(N_Energy_Terms);

                for(r=0;r<5;++r)
                {
                    MGPt[p].MGP[q].EV_Image[i_Image].M_Local[1][r] = &M_Local[p][1][r];
                    MGPt[p].MGP[q].EV_Image[i_Image].M_Local[2][r] = &M_Local[p][2][r];
                    MGPt[p].MGP[q].EV_Image[i_Image].M_Local[3][r] = &M_Local[p][3][r];
                }

                for (r = 0; r < N_Energy_Terms; ++r)
                {
                    MGPt[p].MGP[q].EV_Image[i_Image].Energy_Coefficient[r] = IPN->Energy_Coefficient[r];
                }

                Point_Index = q + CalculatePointsThreads_Lower;
                Point_Index=CPP->Sort_Calculate_Points[Point_Index];
                //printf("%d\n",Point_Index);
                i_G                  = CPP->Local_Points_1D[Point_Index][0][1];
                j_G                  = CPP->Local_Points_1D[Point_Index][0][2];
                k_G                  = CPP->Local_Points_1D[Point_Index][0][3];
                /////////////////////////////////////////////////////////////////////////////////////////////////
                for(I=0;I<7;++I)
                {
                    MGPt[p].MGP[q].EV_Image[i_Image].LocalEnergy[I] = CPP->Local_Points_1D[Point_Index][I][5];

                    Boundary_Type           = CPP->Local_Points_1D[Point_Index][I][0];
                    i_G                     = CPP->Local_Points_1D[Point_Index][I][1];
                    j_G                     = CPP->Local_Points_1D[Point_Index][I][2];
                    k_G                     = CPP->Local_Points_1D[Point_Index][I][3];
                    //printf("%d\t%d\t%d\t%d\t%d\n",Boundary_Type,i_G,j_G,k_G,v);
                    for(J=0;J<7;++J)
                    {
                        if(Boundary_Type!=0)
                        {
                            MGPt[p].MGP[q].EV_Image[i_Image].CST[J][I]  = CST[J][i_G][j_G][k_G];
                        }else
                        {
                            if(CPP->Local_Points_1D[Point_Index][I][5]==1)
                            {
                                MGPt[p].MGP[q].EV_Image[i_Image].CST[J][I] = CST[J][i_G][j_G][k_G];
                            }else
                            {
                                MGPt[p].MGP[q].EV_Image[i_Image].CST[J][I] = 0.;
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
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[1][o] = &Mx_Image[i_Image][i_G][j_G][k_G];
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[2][o] = &My_Image[i_Image][i_G][j_G][k_G];
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[3][o] = &Mz_Image[i_Image][i_G][j_G][k_G];
                    }else
                    {
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[1][o] = &MxV_Image[i_Image][v];
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[2][o] = &MyV_Image[i_Image][v];
                        MGPt[p].MGP[q].EV_Image[i_Image].M_Local[3][o] = &MzV_Image[i_Image][v];
                    }
                }

                MGPt[p].MGP[q].EV_Image[i_Image].M0[0] = *MGPt[p].MGP[q].EV_Image[i_Image].M_Local[1][1];
                MGPt[p].MGP[q].EV_Image[i_Image].M0[1] = *MGPt[p].MGP[q].EV_Image[i_Image].M_Local[2][1];
                MGPt[p].MGP[q].EV_Image[i_Image].M0[2] = *MGPt[p].MGP[q].EV_Image[i_Image].M_Local[3][1];
            }

            ///////////////////////////////////////////////////////////////////////////////////////
            EV_Initial[p][q].dxy = Make1DArray(3);
            EV_Initial[p][q].dxy[0] = dxy[0];    EV_Initial[p][q].dxy[1] = dxy[1];
            EV_Initial[p][q].dxy[2] = dxy[2];    EV_Initial[p][q].i      = 3;
            EV_Initial[p][q].j      = 3;         EV_Initial[p][q].k      = 3;
            EV_Initial[p][q].m      = m;         EV_Initial[p][q].n      = n;
            EV_Initial[p][q].l      = l;
            EV_Initial[p][q].h      = h;         EV_Initial[p][q].t      = t;
            EV_Initial[p][q].M0          = Make1DArray(3);
            EV_Initial[p][q].LocalEnergy = Make1DArrayinteger(7);
            EV_Initial[p][q].CST         = Make2DArray(7,7);
            EV_Initial[p][q].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
            EV_Initial[p][q].M_Local     = Make2DArrayPointer(4,15);
            EV_Initial[p][q].Energy_Coefficient = Make1DArray(N_Energy_Terms);

            EV_Final[p][q].dxy = Make1DArray(3);
            EV_Final[p][q].dxy[0] = dxy[0];    EV_Final[p][q].dxy[1] = dxy[1];
            EV_Final[p][q].dxy[2] = dxy[2];    EV_Final[p][q].i      = 3;
            EV_Final[p][q].j      = 3;         EV_Final[p][q].k      = 3;
            EV_Final[p][q].m      = m;         EV_Final[p][q].n      = n;
            EV_Final[p][q].l      = l;
            EV_Final[p][q].h      = h;         EV_Final[p][q].t      = t;
            EV_Final[p][q].M0          = Make1DArray(3);
            EV_Final[p][q].LocalEnergy = Make1DArrayinteger(7);
            EV_Final[p][q].CST         = Make2DArray(7,7);
            EV_Final[p][q].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
            EV_Final[p][q].M_Local     = Make2DArrayPointer(4,15);
            EV_Final[p][q].Energy_Coefficient = Make1DArray(N_Energy_Terms);


            for(r=0;r<5;++r)
            {
                EV_Initial[p][q].M_Local[1][r] = &M_Local[p][1][r];
                EV_Initial[p][q].M_Local[2][r] = &M_Local[p][2][r];
                EV_Initial[p][q].M_Local[3][r] = &M_Local[p][3][r];

                EV_Final[p][q].M_Local[1][r] = &M_Local[p][1][r];
                EV_Final[p][q].M_Local[2][r] = &M_Local[p][2][r];
                EV_Final[p][q].M_Local[3][r] = &M_Local[p][3][r];
            }

            for (r = 0; r < N_Energy_Terms; ++r)
            {
                EV_Initial[p][q].Energy_Coefficient[r] = IPN->Energy_Coefficient[r];
                EV_Final[p][q].Energy_Coefficient[r]   = IPN->Energy_Coefficient[r];
            }

            Point_Index=q+CalculatePointsThreads_Lower;
            Point_Index=CPP->Sort_Calculate_Points[Point_Index];
            //printf("%d\n",Point_Index);
            i_G                  = CPP->Local_Points_1D[Point_Index][0][1];
            j_G                  = CPP->Local_Points_1D[Point_Index][0][2];
            k_G                  = CPP->Local_Points_1D[Point_Index][0][3];

            /////////////////////////////////////////////////////////////////////////////////////////////////
            for(I=0;I<7;++I)
            {
                EV_Initial[p][q].LocalEnergy[I] = CPP->Local_Points_1D[Point_Index][I][5];
                EV_Final[p][q].LocalEnergy[I] = CPP->Local_Points_1D[Point_Index][I][5];

                Boundary_Type           = CPP->Local_Points_1D[Point_Index][I][0];
                i_G                     = CPP->Local_Points_1D[Point_Index][I][1];
                j_G                     = CPP->Local_Points_1D[Point_Index][I][2];
                k_G                     = CPP->Local_Points_1D[Point_Index][I][3];
                //printf("%d\t%d\t%d\t%d\t%d\n",Boundary_Type,i_G,j_G,k_G,v);
                for(J=0;J<7;++J)
                {
                    if(Boundary_Type!=0)
                    {
                        EV_Initial[p][q].CST[J][I]  = CST[J][i_G][j_G][k_G];
                        EV_Final[p][q].CST[J][I]  = CST[J][i_G][j_G][k_G];
                    }else
                    {
                        if(CPP->Local_Points_1D[Point_Index][I][5]==1)
                        {
                            EV_Initial[p][q].CST[J][I] = CST[J][i_G][j_G][k_G];
                            EV_Final[p][q].CST[J][I] = CST[J][i_G][j_G][k_G];
                        }else
                        {
                            EV_Initial[p][q].CST[J][I] = 0.;
                            EV_Final[p][q].CST[J][I] = 0.;
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
                    EV_Initial[p][q].M_Local[1][o] = &Mx_Initial[i_G][j_G][k_G];
                    EV_Initial[p][q].M_Local[2][o] = &My_Initial[i_G][j_G][k_G];
                    EV_Initial[p][q].M_Local[3][o] = &Mz_Initial[i_G][j_G][k_G];

                    EV_Final[p][q].M_Local[1][o]   = &Mx_Final[i_G][j_G][k_G];
                    EV_Final[p][q].M_Local[2][o]   = &My_Final[i_G][j_G][k_G];
                    EV_Final[p][q].M_Local[3][o]   = &Mz_Final[i_G][j_G][k_G];
                }else
                {
                    EV_Initial[p][q].M_Local[1][o] = &MxV_Initial[v];
                    EV_Initial[p][q].M_Local[2][o] = &MyV_Initial[v];
                    EV_Initial[p][q].M_Local[3][o] = &MzV_Initial[v];

                    EV_Final[p][q].M_Local[1][o]   = &MxV_Final[v];
                    EV_Final[p][q].M_Local[2][o]   = &MyV_Final[v];
                    EV_Final[p][q].M_Local[3][o]   = &MzV_Final[v];
                }
            }
        }
    }

///////////////////////// Free Part /////////////////////////
    free4DArrayinteger(BS->LocalEnergy,m,n,l);
    free(CPP->Virtual_Position);
    free4DArrayinteger(Local_Points_1D_LocalEnergy,BS->Calculate_Points,7,7);
    free2DArrayinteger(CPP->Neighbour_Points_3Dto1D_Array,14);

///////////////////////// Synchronize Parameters Initial /////////////////////////
    Modify_GNEB_Synchronize *MGS;
    MGS = &DRGP.MGS;
    MGS->m=m;    MGS->n=n;    MGS->l=l;    MGS->Virtual_Points  = BS->Virtual_Points;
    MGS->Inner_Points=BS->Inner_Points;  MGS->Boundary_Points = BS->Boundary_Points;
    MGS->Pthread_Number=Threads_Number;
    MGS->N_Image       =N_Image;
    MGS->Bound_CalculatePoint_Threads = Make2DArrayinteger(Threads_Number,3);
    MGS->Signal                       = Signal;
    MGS->PTPS                         = PTPS;

    MGS->h = h;
    MGS->t = t;
    MGS->dxy = dxy;

    MGS->Sort_Calculate_Points        = Make1DArrayinteger(BS->Calculate_Points);
    MGS->Virtual_Points               = BS->Virtual_Points;
    MGS->Version                      = Version;
    MGS->CST                          = CST;
    MGS->Local_Points_1D              = CPP->Local_Points_1D;
    MGS->Position_3DTo1D              = CPP->Position_3DTo1D;
    MGS->Energy_Coefficient           = IPN->Energy_Coefficient;

    for(p=0;p<BS->Calculate_Points;++p)
    {
        MGS->Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
    }

    for(p=0;p<Threads_Number;++p)
    {
        MGS->Bound_CalculatePoint_Threads[p][0] = MGPt[p].CalculatePointsThreads_Upper-MGPt[p].CalculatePointsThreads_Lower;
        MGS->Bound_CalculatePoint_Threads[p][1] = MGPt[p].CalculatePointsThreads_Lower;
        MGS->Bound_CalculatePoint_Threads[p][2] = MGPt[p].CalculatePointsThreads_Upper;
    }
    
    DRGP.PTPS = PTPS;
    DRGP.EV_Initial = EV_Initial;
    DRGP.EV_Final   = EV_Final;
    DRGP.Energy_Change = Energy_Change;
    DRGP.Signal = Signal;
    DRGP.M_Local = M_Local;
    DRGP.CST = CST;
    DRGP.MGPt = MGPt;
    DRGP.Energy_Image = Energy_Image;
    DRGP.SKN_Image    = SKN_Image;
    ///////////////////////// Pthreads Part /////////////////////////
    int err;
    pthread_t *tid;
    tid = (pthread_t *)malloc(sizeof(pthread_t) * (Threads_Number + 1));
    DRGP.tid = tid;
    for(p=0;p<Threads_Number;++p)
    {
        err=pthread_create(&tid[p],NULL,Modify_Pthread_GNEB_Transfer,&MGPt[p]);
        if(err != 0)
        {
            printf("create thread error\n");
            Sleep(1000);
        }
    }

    err=pthread_create(&tid[Threads_Number],NULL,Modify_Pthread_GNEB_Synchronize_Transfer,MGS);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }

    ////////////////////// Free pthread //////////////////////////////////////////
    pthread_t tid_Part_Free;
    err = pthread_create(&tid_Part_Free, NULL, DiscreatPoint_Run_GNEB_Part_Free_Transfer, IPN);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }

    return 1;
}

int DiscreatPoint_Run_GNEB_Part_Free(Input_Parameter_NMD *IPN)
{
    int m, n, l, p, q, i_Image, CalculatePointsThreads;
    for (p = 0; p < IPN->Threads_Number + 1; ++p)
    {
        pthread_join(DRGP.tid[p],NULL);
    }
    printf("free thread\n");
    m = 2 * IPN->xn + 3;    n = 2 * IPN->yn + 3;
    if (IPN->zn == 0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }


    ///////////////////////// Result Part /////////////////////////
    double W, SKN;
    W=0.; SKN=0.;
    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        W   = W + DRGP.MGPt[p].W;
        SKN = SKN + DRGP.MGPt[p].SKN;
    }
    W = W / (DRGP.BS.Boundary_Points + DRGP.BS.Inner_Points);
    IPN->W = W; IPN->SKN = SKN;
    printf("1\n");
    Energy_Variables_Distribution EVD;
    EVD.CST = DRGP.CST;    EVD.dxy = IPN->dxy;
    EVD.h   = IPN->h;      EVD.t   = IPN->t;
    EVD.m   = m;      EVD.n   = n;
    EVD.l     = l;
    EVD.Mx    = Mx_Initial;
    EVD.My    = My_Initial;
    EVD.Mz    = Mz_Initial;
    EVD.MxV   = MxV_Initial;
    EVD.MyV   = MyV_Initial;
    EVD.MzV   = MzV_Initial;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D = DRGP.CPP.Local_Points_1D;
    EVD.Position_3DTo1D = DRGP.CPP.Position_3DTo1D;
    EVD.Energy_Coefficient = IPN->Energy_Coefficient;
    Energy_Distribution(&EVD);
    DRGP.Energy_Image[0] = EVD.W_Average;    DRGP.SKN_Image[0] = EVD.SKN;
    EVD.Mx    = Mx_Final;
    EVD.My    = My_Final;
    EVD.Mz    = Mz_Final;
    EVD.MxV   = MxV_Final;
    EVD.MyV   = MyV_Final;
    EVD.MzV   = MzV_Final;
    Energy_Distribution(&EVD);
    DRGP.Energy_Image[IPN->N_Image+1] = EVD.W_Average;    DRGP.SKN_Image[IPN->N_Image+1] = EVD.SKN;
    Energy_Variation_Image = 0.;
    for (i_Image = 1; i_Image <= IPN->N_Image; ++i_Image)
    {
        EVD.Mx  = Mx_Image[i_Image-1];
        EVD.My  = My_Image[i_Image-1];
        EVD.Mz  = Mz_Image[i_Image-1];
        EVD.MxV = MxV_Image[i_Image-1];
        EVD.MyV = MyV_Image[i_Image-1];
        EVD.MzV = MzV_Image[i_Image-1];
        Energy_Distribution(&EVD);
        DRGP.Energy_Image[i_Image] = EVD.W_Average;    DRGP.SKN_Image[i_Image] = EVD.SKN;
        Energy_Variation_Image = Energy_Variation_Image + (DRGP.Energy_Image[i_Image] - DRGP.Energy_Image[IPN->N_Image + 1]) / IPN->N_Image;
    }
    Write_Image(IPN, &DRGP.BS, Mx_Image, My_Image, Mz_Image, MxV_Image, MyV_Image, MzV_Image);
    Write_Image_Result(DRGP.Energy_Image, DRGP.SKN_Image, IPN->N_Image);
    ///////////////////////// Free Part /////////////////////////
    for(p=0;p<IPN->Threads_Number;++p)
    {
        free(DRGP.MGPt[p].Signal);
        free(DRGP.PTPS[p].Energy_Image);
        free(DRGP.PTPS[p].SKN_Image);

        CalculatePointsThreads = DRGP.MGPt[p].CalculatePointsThreads_Upper - DRGP.MGPt[p].CalculatePointsThreads_Lower;
        for(q=0;q<CalculatePointsThreads;++q)
        {
            for (i_Image = 0; i_Image < IPN->N_Image; ++i_Image)
            {
                free(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].dxy);
                free(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].M0);
                free(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].LocalEnergy);
                free2DArray(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].CST,7);
                free2DArrayPointer(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].M_Local,4);
                free(DRGP.MGPt[p].MGP[q].EV_Image[i_Image].Energy_Coefficient);
            }
            free(DRGP.MGPt[p].MGP[q].EV_Image);
            free(DRGP.EV_Initial[p][q].dxy);
            free(DRGP.EV_Initial[p][q].M0);
            free(DRGP.EV_Initial[p][q].LocalEnergy);
            free2DArray(DRGP.EV_Initial[p][q].CST,7);
            free2DArrayPointer(DRGP.EV_Initial[p][q].M_Local,4);
            free(DRGP.EV_Initial[p][q].Energy_Coefficient);

            free(DRGP.EV_Final[p][q].dxy);
            free(DRGP.EV_Final[p][q].M0);
            free(DRGP.EV_Final[p][q].LocalEnergy);
            free2DArray(DRGP.EV_Final[p][q].CST,7);
            free2DArrayPointer(DRGP.EV_Final[p][q].M_Local,4);
            free(DRGP.EV_Final[p][q].Energy_Coefficient);
        }
        free(DRGP.MGPt[p].MGP);
        free(DRGP.EV_Initial[p]);
        free(DRGP.EV_Final[p]);

        for (q = 0; q < 3; ++q)
        {
            pthread_mutex_destroy(&Signal_Synchronize_GNEB_Mutex[p][q]);
            pthread_cond_destroy(&Signal_Synchronize_GNEB_Cond[p][q]);
        }
        free(Signal_Synchronize_GNEB_Mutex[p]);
        free(Signal_Synchronize_GNEB_Cond[p]);
        free(DRGP.MGPt[p].Sort_Calculate_Points);
        pthread_mutex_destroy(&M_Synchronize_GNEB_Mutex[p]);
    }
    free4DArray(DRGP.CST,7,m,n);
    free2DArrayinteger(DRGP.MGS.Bound_CalculatePoint_Threads,IPN->Threads_Number);
    free(DRGP.MGS.Sort_Calculate_Points);
    free3DArray(DRGP.M_Local,IPN->Threads_Number,4);
    free3DArrayinteger(DRGP.CPP.Local_Points_1D,DRGP.BS.Calculate_Points,7);
    free3DArrayinteger(DRGP.CPP.Position_3DTo1D,m,n);
    free(DRGP.CPP.Sort_Calculate_Points);
    free(DRGP.MGPt);
    free(DRGP.EV_Initial);
    free(DRGP.EV_Final);
    free2DArrayintegerPointer(DRGP.Signal,IPN->Threads_Number+1);
    free(Signal_Synchronize_GNEB_Mutex);
    free(Signal_Synchronize_GNEB_Cond);
    free(M_Synchronize_GNEB_Mutex);
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);

    free4DArray(Mx_Image, IPN->N_Image, m, n);
    free4DArray(My_Image, IPN->N_Image, m, n);
    free4DArray(Mz_Image, IPN->N_Image, m, n);
    free2DArray(MxV_Image, IPN->N_Image);
    free2DArray(MyV_Image, IPN->N_Image);
    free2DArray(MzV_Image, IPN->N_Image);

    free3DArrayinteger(DRGP.BS.Boundary,m,n);
    free(DRGP.Energy_Image);
    free(DRGP.SKN_Image);
    free(DRGP.PTPS);
    free(DRGP.tid);
    ///////////////////////// Result Part /////////////////////////
    double end_time, duration;
    end_time  = clock();
    duration  = (double)(end_time-DRGP.start_time)/CLOCKS_PER_SEC;
    IPN->Time = duration;

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Finish;
    pthread_cond_signal(&GPM.Run_GUI_Cond);
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    printf("free done\n");
    return 1;
}

void* DiscreatPoint_Run_GNEB_Part_Free_Transfer(void *arg)
{
    Input_Parameter_NMD *IPN;
    IPN=(Input_Parameter_NMD*)arg;
    DiscreatPoint_Run_GNEB_Part_Free(IPN);
    return NULL;
}

int DiscreatPoint_Run_GNEB(Input_Parameter_NMD *IPN)
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
            DiscreatPoint_Run_GNEB_Part_Run(IPN);
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

