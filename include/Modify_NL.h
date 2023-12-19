#ifndef MODIFY_NL_H_INCLUDED
#define MODIFY_NL_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

double Minimal_Polynome(double *Coefficient, int Degree, double *Interval, double Precision);

double Minimal_G_L_Fun(Modify_NL_Parameters *MNP, double *Gradient_Local, double *Interval, double Precision);

int Modify_Nonlinear_OnePoint(Modify_NL_Parameters *MNP);

int Modify_Pthread_NL(Modify_NL_Pthreads *MNPt);

void* Modify_Pthread_NL_Transfer(void *arg);

int Modify_Pthread_NL_Synchronize(Modify_NL_Synchronize *MNS);

void* Modify_Pthread_NL_Synchronize_Transfer(void *arg);

int OPGL_NL_Pthread(OPGL_NL_Pthread_Parameters *ONPP);

void *OPGL_NL_Pthread_Transfer(void *arg);

int DiscreatPoint_Run_NL(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_NL_Part_Run(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_NL_Part_Free(Input_Parameter_NMD *IPN);

void *DiscreatPoint_Run_NL_Part_Free_Transfer(void *arg);

#endif // MODIFY_NL_H_INCLUDED
