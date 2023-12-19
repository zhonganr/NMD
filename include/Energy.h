#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED
#include <iostream>
#include "Variable.h"

using namespace std;

int Effective_Field_LocalCell_M0(Energy_Variables *EV, double *H_eff_Local);

double Energy_Single(Energy_Variables *EV);

double Energy_Local(Energy_Variables *EV);

double Energy_Local_M0(Energy_Variables *EV);

int Continue_Modify(Energy_Variables *EV);

double Skyrmion_Number_Single(Energy_Variables *EV);

int Energy_Distribution(Energy_Variables_Distribution *EVD);

double Skyrmion_Number_Incircle(Input_Parameter_NMD *IPN);

#endif // ENERGY_H_INCLUDED
