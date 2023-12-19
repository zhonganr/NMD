#ifndef VARIABLE_H_INCLUDED
#define VARIABLE_H_INCLUDED

#include <iostream>
#include "pthread.h"

//Enumuration
enum Initial_Mode_NMD{Initial_Helical=-1, Read_M_3D, Read_M_2D, Initial_M_FE, Read_M_Process, Initial_LLG, IMN_1Circle, IMN_2Circles, IMN_3Circles, 
                      IMN_4Circles, IMN_5Circles, IMN_6Circles, IMN_7Circles, IMN_8Circles, IMN_9Circles, IMN_10Circles, IMN_11Circles};
enum Energy_Terms{Exchange, DM_Bloch, DM_Neel, Zeeman, Landau, Magnetoelastic, Axial_Anisotropy};
enum Elastic_Field_Mode{EFM_Zerofield, EFM_UniaxialStrain, EFM_UniaxialStress, EFM_EdgeDislocation, EFM_ScrewDislocation, EFM_CrackI, 
                        EFM_ConcentratedForce, EFM_BeamUniformForce, EFM_BeamConcentratedForce, EFM_BeamBending, EFM_ReadAnsys, EFM_HollowTorsion, 
                        EFM_EllipseTorsion};

enum Boundary_Initial_Type_NMD{BITN_Rectangle, BITN_Circle, BITN_Circle_Hole, BITN_Ellipse, BITN_Triangle};
enum Calculate_Mode_GUI{CMG_MC, CMG_NL, CMG_LLG_SOT, CMG_LLG_STT, CMG_GNEB};
enum Thread_State_GUI{TSG_Prepare, TSG_Run, TSG_Finish};
//Signal of Synchronize thread to the Child thread
enum Signal_Synchronize{SignalSyn_Run, SignalSyn_PrepareToSynchronize, SignalSyn_Synchronize, SignalSyn_Stop, SignalSyn_NotActive};
//Signal of Child thread to the Synchronize thread
enum Signal_Pthread{SignalPthread_Stop, SignalPthread_Run};
enum State_Pthread{StatePthread_PrepareToRun, StatePthread_Run, StatePthread_Wait, StatePthread_PrepareToSynchronize, StatePthread_Synchronize, StatePthread_Stop};
enum State_Synchronize{StateSyn_PrepareToRun, StateSyn_Run, StateSyn_PrepareToSynchronize, StateSyn_Synchronize, StateSyn_Wait, StateSyn_Stop};
enum Signal_GUI{Signal_GUI_Stop, Signal_GUI_Run, Signal_GUI_Break, Signal_GUI_Continue};
enum OPGL_Mode_GUI{OPGL_Mode_Off, OPGL_Mode_On};
enum Version_Mode_GUI{VM_Plus, VM_Minus};
enum State_Pthread_Strain_GUI{SPS_Prepare, SPS_Run, SPS_Finish};

struct Boundary_Shape
{
    int xn, yn, zn, ***Boundary, ****LocalEnergy;
    int Boundary_Type_x, Boundary_Type_y, Boundary_Type_z;
    //Consider the a point (i,j,k), Boundary[i][j][k] represents the Boundary Type on this point.
    //Boundary[i][j][k]=1 is inner side, 2 is on the boundary, 0 is outside.
    int Boundary_Points, Inner_Points, Virtual_Points, Calculate_Points;
};

struct Calculate_Points_Parameters
{
    int ***Local_Points_1D, ***Position_3DTo1D, *Virtual_Position, ****Local_Points_1D_LocalEnergy;
    int ***Neighbour_Points_Local_Coordinate, **Neighbour_Points_3Dto1D_Array;
    int *Sort_Calculate_Points;
};

struct Input_Parameter_NMD
{
    int xn, yn, zn, Version, VersionMode, InitialMode, NCycles, Threads_Number, ElasticFieldMode;
    int MagneticFieldMode, MagneticFieldPeriod, TemperatureFieldMode, TemperatureFieldPeriod;
    int Boundary_Type_x, Boundary_Type_y, Boundary_Type_z, Boundary_Shape_Type, CalculateMode;
    int Continuity_Modify_Mode, Inplane_Density, N_Image, OPGL_Mode;
    double *dxy, *CST_Uniaxial, h, t, Continuity_Modify_Coefficient;
    double b, a, FC, FU_Beam, FC_Beam, B_M, ElasticFieldScaleFactor;
    double W, SKN, Time, SKN_Area, SKN_Positif;
    double dT, MpX, MpY, MpZ, CurrentDensity, DerivativeM;
    double *Energy_Coefficient;
};

struct GUI_Parameter
{
    //Run: Thread_State[i][0]=1, Stop: Thread_State[i][0]=0, Continue: Thread_State[i][1]=1, Break: Thread_State[i][1]=0
    int **Thread_State, State_GUI_Run;
    pthread_cond_t *Pthread_GUI_Cond, Run_GUI_Cond;
    pthread_mutex_t *Pthread_GUI_Mutex, Run_GUI_Mutex;
    
};

struct Pthread_State_Run_Parameter
{
    Boundary_Shape BS;
    Calculate_Points_Parameters CPP;
    pthread_t *tid_Part_Run, *tid_Part_Free;
    char *file_SKN;
    FILE *fp_Note;
    int Version, *State_Pthread_Strain;
    pthread_mutex_t State_Pthread_Strain_Mutex;
};

struct Pthread_State_Run_Strain_Parameter
{
    int Index_Pthread, N_Pthread;
    char *dir_SKN;
};

struct Energy_Variables
{
    double h, t, *dxy, W_Single, W_Local_BC, W_Local_AC, dT, CurrentDensity, MpX, MpY, MpZ;
    double ***M_Local, *M0, **CST, Continuity_Modify_Coefficient, *M_Target;
    int i, j, k, m, n, l, *LocalEnergy;
    double *Energy_Coefficient;
};

struct Parameter_Transfer_Pthread_Synchronize
{
    double *Energy_Image, *SKN_Image;
    double Energy, SKN;
};


struct Modify_MC_Parameters
{
    Energy_Variables *EV;
    double W_BeforeModify, W_AfterModify, *Cv, *T, *M0;
    int Energy_Variation;
    //Energy_Variation=0 is not decreased and 1 is decreased.
};

struct Modify_MC_Pthreads
{
    Boundary_Shape *BS;
    Modify_MC_Parameters *MMP;
    int Pthread_Index, Pthread_Number, NCycles, *Signal;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    double ***Mx_Global_Pthread, ***My_Global_Pthread, ***Mz_Global_Pthread;
    double *MxV_Global_Pthread, *MyV_Global_Pthread, *MzV_Global_Pthread;
    double W, SKN, *Cv, *T;
    int ***Local_Points_1D, ***Position_3DTo1D, *Sort_Calculate_Points;
    double ****CST, *dxy, h, t;
    int m, n, l;
    int Continuity_Modify_Mode;
    double Continuity_Modify_Coefficient;
};

struct Modify_MC_Synchronize
{
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points, Boundary_Points, Inner_Points;
    int Pthread_Number, ***Signal, Version;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST, *dxy, h, t;
    int *Sort_Calculate_Points;
    double W, SKN, *Energy_Coefficient;
    Input_Parameter_NMD IPN;
};

struct OPGL_MC_Pthread_Parameters
{
    Input_Parameter_NMD *IPN;
    double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
    double ****CST;
};

struct DiscreatPoint_Run_MC_Parameter
{
    double ****CST, ***M_Local, start_time, end_time;
    int ***Signal;
    Boundary_Shape BS;
    Calculate_Points_Parameters CPP;
    Energy_Variables **EV;
    Modify_MC_Pthreads *MMPt;
    Modify_MC_Synchronize MMS;
    OPGL_MC_Pthread_Parameters OMPP;
    pthread_t *tid;
};

struct Modify_GNEB_Parameters
{
    Energy_Variables *EV_Initial, *EV_Final, *EV_Image;
    double *M0;
    int N_Image;
    //Energy_Variation=0 is not decreased and 1 is decreased.
};

struct Modify_GNEB_Pthreads
{
    Boundary_Shape *BS;
    Modify_GNEB_Parameters *MGP;
    int Pthread_Index, Pthread_Number, NCycles, *Signal;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    double W, SKN, *Energy_Image, *SKN_Image;
    int ***Local_Points_1D, ***Position_3DTo1D, *Sort_Calculate_Points;
    double ****CST, *dxy, h, t;
    int m, n, l;
    int Continuity_Modify_Mode, N_Image;
    double Continuity_Modify_Coefficient;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct Modify_GNEB_Synchronize
{
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points, Boundary_Points, Inner_Points;
    int Pthread_Number, ***Signal, Version;
    int **Bound_CalculatePoint_Threads;
    double ****CST, *dxy, h, t;
    int *Sort_Calculate_Points, N_Image;
    double W, SKN, *Energy_Coefficient;
    Input_Parameter_NMD IPN;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct DiscreatPoint_Run_GNEB_Parameter
{
    double ****CST, ***M_Local, start_time, end_time, ***Energy_Change, *Energy_Image, *SKN_Image;
    int ***Signal;
    Boundary_Shape BS;
    Calculate_Points_Parameters CPP;
    Energy_Variables **EV_Initial, **EV_Final;
    Modify_GNEB_Pthreads *MGPt;
    Modify_GNEB_Synchronize MGS;
    pthread_t *tid;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct Modify_Target_Parameters
{
    Energy_Variables *EV;
    double W_BeforeModify, W_AfterModify, *M0;
    int Energy_Variation;
    double W_Change;
    //Energy_Variation=0 is not decreased and 1 is decreased.
};

struct Modify_Target_Pthreads
{
    Boundary_Shape *BS;
    Modify_Target_Parameters *MTP;
    int Pthread_Index, Pthread_Number, NCycles, *Signal;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    double ***Mx_Global_Pthread, ***My_Global_Pthread, ***Mz_Global_Pthread;
    double *MxV_Global_Pthread, *MyV_Global_Pthread, *MzV_Global_Pthread;
    double W, SKN;
    int ***Local_Points_1D, ***Position_3DTo1D, *Sort_Calculate_Points;
    double ****CST, *dxy, h, t;
    int m, n, l;
    int Continuity_Modify_Mode;
    double Continuity_Modify_Coefficient;
    double *Energy_Change;
};

struct Modify_Target_Synchronize
{
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points, Boundary_Points, Inner_Points;
    int Pthread_Number, ***Signal, Version;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST, *dxy, h, t;
    int *Sort_Calculate_Points;
    double W, SKN;
    double ***Energy_change;
    Input_Parameter_NMD IPN;
};

struct OPGL_Target_Pthread_Parameters
{
    Input_Parameter_NMD *IPN;
    double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
    double ****CST;
};

struct OPGL_NL_Pthread_Parameters
{
    Input_Parameter_NMD *IPN;
    double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
    double ****CST;
};

struct OPGL_LLG_Pthread_Parameters
{
    Input_Parameter_NMD *IPN;
    double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
    double ****CST;
};

struct OPGL_Parameters
{
    Input_Parameter_NMD *IPN;
    double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
    double ****CST;
};

struct OPGL_Parameters_GUI
{
    double Energy, SKN;
    int State_Syn;
};


struct Modify_NL_Parameters
{
    Energy_Variables *EV;
    double W_BeforeModify, W_AfterModify, *Gradient, *M0;
    int Energy_Variation;
    //Energy_Variation=0 is not decreased and 1 is decreased.
};

struct Modify_NL_Pthreads
{
    Boundary_Shape *BS;
    //Calculate_Points_Parameters *CPP;
    //Energy_Variables *EV;
    Modify_NL_Parameters *MNP;
    int Pthread_Index, Pthread_Number, NCycles, *Signal;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    double ***Mx_Global_Pthread, ***My_Global_Pthread, ***Mz_Global_Pthread;
    double *MxV_Global_Pthread, *MyV_Global_Pthread, *MzV_Global_Pthread;
    double W, SKN;
    int ***Local_Points_1D, ***Position_3DTo1D, *Sort_Calculate_Points;
    double ****CST, *dxy, h, t;
    int m, n, l;
    int Continuity_Modify_Mode;
    double Continuity_Modify_Coefficient;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct Modify_NL_Synchronize
{
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points, Boundary_Points, Inner_Points;
    int Pthread_Number, ***Signal, Version;
    int **Bound_CalculatePoint_Threads;
    double ****Mx_Global_Pthread, ****My_Global_Pthread, ****Mz_Global_Pthread;
    double **MxV_Global_Pthread, **MyV_Global_Pthread, **MzV_Global_Pthread;
    double ***Mx_Global_Synchronize, ***My_Global_Synchronize, ***Mz_Global_Synchronize;
    double *MxV_Global_Synchronize, *MyV_Global_Synchronize, *MzV_Global_Synchronize;
    double ****CST, *dxy, h, t, *Energy_Coefficient;
    int *Sort_Calculate_Points;
    Input_Parameter_NMD IPN;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct DiscreatPoint_Run_NL_Parameter
{
    double ****CST, ***M_Local, start_time, end_time;
    int ***Signal;
    Boundary_Shape BS;
    Calculate_Points_Parameters CPP;
    Energy_Variables **EV;
    Modify_NL_Pthreads *MNPt;
    Modify_NL_Synchronize MNS;
    OPGL_NL_Pthread_Parameters ONPP;
    pthread_t *tid;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct Modify_LLG_Parameters
{
    Energy_Variables *EV;
    double W_BeforeModify, W_AfterModify, *Gradient, *M0;
    //double dT, MpX, MpY, MpZ, CurrentDensity;
    int Energy_Variation;
    //Energy_Variation=0 is not decreased and 1 is decreased.
};

struct Modify_LLG_Pthreads
{
    Boundary_Shape *BS;
    //Calculate_Points_Parameters *CPP;
    //Energy_Variables *EV;
    Modify_LLG_Parameters *MLP;
    int Pthread_Index, Pthread_Number, NCycles, *Signal, MagneticFieldMode, MagneticFieldPeriod;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    double W, SKN, Energy_Pthread, SKN_Pthread;
    int ***Local_Points_1D, ***Position_3DTo1D, *Sort_Calculate_Points;
    double ****CST, *dxy, h, t;
    int m, n, l, CalculateMode;
    int Continuity_Modify_Mode;
    double Continuity_Modify_Coefficient;
    double DerivativeM;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
    //double dT, CurrentDensity, MpX, MpY, MpZ;
};

struct Modify_LLG_Synchronize
{
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points, Boundary_Points, Inner_Points;
    int Pthread_Number, ***Signal, Version;
    int **Bound_CalculatePoint_Threads;
    double ****CST, *dxy, h, t;
    int *Sort_Calculate_Points;
    double W, SKN, *Energy_Coefficient;
    Input_Parameter_NMD IPN;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
};

struct DiscreatPoint_Run_LLG_Parameter
{
    double ****CST, ***M_Local, start_time, end_time;
    int ***Signal;
    Boundary_Shape BS;
    Calculate_Points_Parameters CPP;
    Energy_Variables **EV;
    Modify_LLG_Pthreads *MLPt;
    Modify_LLG_Synchronize MLS;
    pthread_t *tid;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
    OPGL_LLG_Pthread_Parameters OLPP;
};

struct NMD_Result_Output_Parameter
{
    double W, SKN, h, t, *dxy;
    int xn, yn, zn, ***Boundary;
    double ***Mx, ***My, ***Mz, *MxV, *MyV, *MzV, ****CST;
};

struct Energy_Variables_Distribution
{
    double h, t, *dxy, W_Average, SKN, SKN_Area, SKN_Positif;
    double ***Mx, ***My, ***Mz, *MxV, *MyV, *MzV, ***w, ***SKN_D, ****CST;
    int m, n, l, ***Local_Points_1D, ***Position_3DTo1D;
    double *Energy_Coefficient;
};

#endif // VARIABLE_H_INCLUDED
