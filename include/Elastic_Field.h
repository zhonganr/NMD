#ifndef ELASTIC_FIELD_H_INCLUDED
#define ELASTIC_FIELD_H_INCLUDED

#include <iostream>
#include "Variable.h"

using namespace std;

int StressToStrain(double *Stress, double *Strain);

int StrainToStress(double *Stress, double *Strain);

int BeamCircleHole_Torsion(int m, int n, int l, double ****CST, double *dxy, double B_M);

int BeamEllipse_Torsion(int m, int n, int l, double ****CST, double *dxy, double B_M);

int Screw_Dislocation(int m, int n, int l, double ****CST, double *dxy, double b);

int Edge_Dislocation(int m, int n, int l, double ****CST, double *dxy, double b);

int Uniaxial_Strain(int m, int n, int l, double ****CST, double *CST_Uniaxial);

int Uniaxial_Stress(int m, int n, int l, double ****CST, double *CST_Uniaxial);

int ZeroField(int m, int n, int l, double ****CST);

int CrackI_Stress_CC(int m, int n, int l, double ****CST, double *dxy, double a);

int Concentrated_Force(int m, int n, int l, double ****CST, double *dxy, double FC);

int Elastic_Field_Initial(double ****CST, Input_Parameter_NMD *IPN);

int StressToPS(double *Stress, double *PrincipalStress);

int Read_Elastic_Field_ANSYS(Input_Parameter_NMD *IPN, double ****CST);

#endif // ELASTIC_FILELD_H_INCLUDED
