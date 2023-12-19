#ifndef GNEB_H_INCLUDED
#define GNEB_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

int Case_Energy(double E1, double E2, double E3);

int Cross_Product(double *V1, double *V2, double *V_Result);

int Effective_Field_LocalCell_M0(Energy_Variables *EV, double *H_eff_Local);

int Modify_GNEB_OneStep(Modify_GNEB_Parameters *MGP);

int Modify_Pthread_GNEB(Modify_GNEB_Pthreads *MGPt);

void *Modify_Pthread_GNEB_Transfer(void *arg);

int Modify_Pthread_GNEB_Synchronize(Modify_GNEB_Synchronize *MGS);

void *Modify_Pthread_GNEB_Synchronize_Transfer(void *arg);

int DiscreatPoint_Run_GNEB_Part_Run(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_GNEB_Part_Free(Input_Parameter_NMD *IPN);

void *DiscreatPoint_Run_GNEB_Part_Free_Transfer(void *arg);

int DiscreatPoint_Run_GNEB(Input_Parameter_NMD *IPN);

#endif // GNEB_H_INCLUDED
