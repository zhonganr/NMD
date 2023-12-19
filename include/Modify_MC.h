#ifndef MODIFY_MC_H_INCLUDED
#define MODIFY_MC_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

int Minimal(double *W, int n);

int Selection(int n_Mx, int n_My, int n_Mz, int n_dx, int n_dy, double *Cv, double *dM);

int Modify_MC_OneStep(Modify_MC_Parameters *MMP);

int Rotation(Modify_MC_Parameters *MMP);

int Converge_Paramter_Initial(double *Cv, double *T);

int Rotation_EnergyDecrease(Modify_MC_Parameters *MMP);

int Select_Points_Random(Modify_MC_Pthreads *MMPt);

int Modify_Pthread_MC(Modify_MC_Pthreads *MMPt);

void* Modify_Pthread_MC_Transfer(void *arg);


int Modify_Pthread_MC_Synchronize(Modify_MC_Synchronize *MMS);

void* Modify_Pthread_MC_Synchronize_Transfer(void *arg);

int OPGL_MC_Pthread(OPGL_MC_Pthread_Parameters *OMPP);

void *OPGL_MC_Pthread_Transfer(void *arg);

int DiscreatPoint_Run_MC(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_MC_Part_Run(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_MC_Part_Free(Input_Parameter_NMD *IPN);

void *DiscreatPoint_Run_MC_Part_Free_Transfer(void *arg);

#endif // MODIFY_MC_H_INCLUDED
