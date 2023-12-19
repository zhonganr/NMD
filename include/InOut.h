#ifndef INOUT_H_INCLUDED
#define INOUT_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

int WriteInput_NMD(Input_Parameter_NMD *IPN);

int ReadInput_NMD(Input_Parameter_NMD *IPN);

int WriteVectorM_3D_Monolayer(int m, int n, double ***Mx, double ***My, double ***Mz);

int Make_NMD_Result_Dir(Input_Parameter_NMD *IPN);

int Make_NMD_Result_Dir_SKN(Input_Parameter_NMD *IPN, char *dir_SKN);

int Write_NMD_ElastiField(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS);

int Write_NMD_Magnetization(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                            double *MxV, double *MyV, double *MzV, Boundary_Shape *BS);

int Write_NMD_Magnetization_SKN(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                                double *MxV, double *MyV, double *MzV, Boundary_Shape *BS, char *dir_SKN);

int Write_NMD_Energy(Input_Parameter_NMD *IPN, double ***w);

int Write_NMD_SKN_D(Input_Parameter_NMD *IPN, double ***SKN_D);

int Write_NMD_Result(Input_Parameter_NMD *IPN);

int Simplify_Matrix_3D(Input_Parameter_NMD *IPN, Boundary_Shape *BS);

int Write_Inplane_M(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ***Mx, double ***My);

int Write_Inplane_M_SKN(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ***Mx, double ***My, char *dir_SKN);

int Write_Judge_Phase_Curve(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****CST);

int Read_NMD_SKN_D(Input_Parameter_NMD *IPN, double **SKN_D);

int SKN_Boundary_Found(Input_Parameter_NMD *IPN);

int Stress_Analysis(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS);

int Write_Principal_Stress(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS);

int Write_Principal_Strain(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS);

double Determinant_4D(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4);

int Calculate_Interpolation(double **Interpolate_Coordinate, double *Interpolate_Values);

int FE_Neighbour_Points_Find(double **Grid, double **Coordinate, int ***Node_4Points, int NodeMax, int i, int j);

int ReadAnsys(Input_Parameter_NMD *IPN, double ****CST);

int Write_Image(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****Mx_Image,
                double ****My_Image, double ****Mz_Image, double **MxV_Image, double **MyV_Image, double **MzV_Image);

int Read_Image(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****Mx_Image,
               double ****My_Image, double ****Mz_Image, double **MxV_Image, double **MyV_Image, double **MzV_Image);

int Write_Image_Result(double *Energy_Image, double *SKN_Image, int N_Image);

int Interpolation_NMD_2D(int N_Interpolate_Points, double **Interpolate_Points_Parameter, double *Interpolate_Coefficient);

#endif // IO_H_INCLUDED
