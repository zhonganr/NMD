#ifndef MODIFY_TARGET_H_INCLUDED
#define MODIFY_TARGET_H_INCLUDED

#include <iostream>
#include "Variable.h"

int Modify_Target_OneStep(Modify_Target_Parameters *MTP);

int Select_Points_Random(Modify_Target_Pthreads *MTPt);

int Modify_Pthread_Target(Modify_Target_Pthreads *MTPt);

void *Modify_Pthread_Target_Transfer(void *arg);

int Modify_Pthread_Target_Synchronize(Modify_Target_Synchronize *MTS);

void *Modify_Pthread_Target_Synchronize_Transfer(void *arg);

int OPGL_Target_Pthread(OPGL_Target_Pthread_Parameters *OTPP);

void *OPGL_Target_Pthread_Transfer(void *arg);

int DiscreatPoint_Run_Target(Input_Parameter_NMD *IPN);



#endif // MODIFY_TARGET_H_INCLUDED
