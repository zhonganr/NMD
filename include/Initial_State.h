#ifndef INITIAL_STATE_H_INCLUDED
#define INITIAL_STATE_H_INCLUDED

#include <iostream>
#include "Variable.h"
int Initial_Magnetization_FE(int m, int n, int l, double ***Mx, double ***My, double ***Mz);

int Initial_Magnetization_Center_Circles(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz);

int Initial_Magnetization_Multiple_Circles(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz, int Surround_number);

int Initial_Magnetization_Multiple_Circles_Hollow(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz, int Surround_number);

int InitialM_3D(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                double *MxV, double *MyV, double *MzV, Boundary_Shape *BS);

int InitialM_3D_SKN(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                    double *MxV, double *MyV, double *MzV, Boundary_Shape *BS, char *dir_SKN);

int Initial_Magnetization_LLG(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz);

int Trim_Boundary_M3D_VirtualToBoundary(Boundary_Shape *BS, Calculate_Points_Parameters *CPP, Input_Parameter_NMD *IPN,
                                        double ***Mx, double ***My, double ***Mz, double *MxV, double *MyV, double *MzV);

int Trim_Boundary_M3D_BoundaryToVirtual(Boundary_Shape *BS, Calculate_Points_Parameters *CPP, Input_Parameter_NMD *IPN,
                                        double ***Mx, double ***My, double ***Mz, double *MxV, double *MyV, double *MzV);

int Initial_Magnetization_Helical(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz);

int Initial_Magnetization_10_Circles(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz);

int Initial_Magnetization_11_Circles(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz);

int Initial_Magnetization_LLG_BC(int m, int n, int l, double ***Mx, double ***My, double ***Mz);

int InitialM_Target(double ***Mx, double ***My, double ***Mz,
                    double *MxV, double *MyV, double *MzV, Boundary_Shape *BS);

#endif // INITIAL_STATE_H_INCLUDED
