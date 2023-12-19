#ifndef LLG_H_INCLUDED
#define LLG_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

int Gradient_M_LocalCell(Energy_Variables *EV, double *Gradient_M_Local);

int StochasticLLG_SOT(Modify_LLG_Parameters *MLP);

int StochasticLLG_STT(Modify_LLG_Parameters *MLP);

int Modify_Pthread_LLG(Modify_LLG_Pthreads *MLPt);

void *Modify_Pthread_LLG_Transfer(void *arg);

int Modify_Pthread_LLG_Synchronize(Modify_LLG_Synchronize *MLS);

void *Modify_Pthread_LLG_Synchronize_Transfer(void *arg);

int OPGL_LLG_Pthread(OPGL_LLG_Pthread_Parameters *OLPP);

void *OPGL_LLG_Pthread_Transfer(void *arg);

int DiscreatPoint_Run_LLG_Part_Run(Input_Parameter_NMD *IPN);

int DiscreatPoint_Run_LLG_Part_Free(Input_Parameter_NMD *IPN);

void *DiscreatPoint_Run_LLG_Part_Free_Transfer(void *arg);

int DiscreatPoint_Run_LLG(Input_Parameter_NMD *IPN);



#endif // LLG_H_INCLUDED