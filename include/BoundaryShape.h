#ifndef BOUNDARYSHAPE_H_INCLUDED
#define BOUNDARYSHAPE_H_INCLUDED
#include <iostream>
#include "Variable.h"

using namespace std;

int Boundary_Ellipse(Boundary_Shape *BS);

int Boundary_Circle(int rn, Boundary_Shape *BS);

int Boundary_Circle_Hole(Boundary_Shape *BS);

int Regular_Polygon_Plane_Shape(int rn, int N_Polygon, double **Line_Polygon);

int Boundary_Polygon(Boundary_Shape *BS, int N_Polygon);

int Boundary_Rectangle(Boundary_Shape *BS);

int Boundary_Rectangle_20Grid_Crack(Boundary_Shape *BS);

int Boundary_Rectangle_Rectangle_Hole(Boundary_Shape *BS);

int Boundary_Rectangle_4Circle_Holes(Boundary_Shape *BS);

int Boundary_Rectangle_Circle(Boundary_Shape *BS);

int Boundary_Rectangle_Ellipse(Boundary_Shape *BS);

int Shape_Boundary_Assignment(Boundary_Shape *BS);

int LocalEnergy_Initial(Boundary_Shape *BS);

int CalculPoints_Numbers(Boundary_Shape *BS);

int CalculPointsParameters_Initial(Calculate_Points_Parameters *CPP, Boundary_Shape *BS);

int Neighbour_Points_Local_Initial(Calculate_Points_Parameters *CPP);

int Neighbour_Points_3Dto1D(Calculate_Points_Parameters *CPP);

int Boundary_Initial_Type(Input_Parameter_NMD *IPN, Boundary_Shape *BS);

int ReSort_Calcualte_Points(Calculate_Points_Parameters *CPP,Boundary_Shape *BS);

int ReSort_Calcualte_Points_ByThreads(Calculate_Points_Parameters *CPP,Boundary_Shape *BS, int Threads_Number);

int ReSort_Calcualte_Points_BySymmetry_Xaxis(Calculate_Points_Parameters *CPP,Boundary_Shape *BS);

int ReSort_Calcualte_Points_BySymmetry_Yaxis(Calculate_Points_Parameters *CPP,Boundary_Shape *BS);

int ReSort_Calcualte_Points_ByRow(Calculate_Points_Parameters *CPP,Boundary_Shape *BS);

int ReSort_Calcualte_Points_ByOrders(Calculate_Points_Parameters *CPP,Boundary_Shape *BS);

#endif // BOUNDARYSHAPE_H_INCLUDED
