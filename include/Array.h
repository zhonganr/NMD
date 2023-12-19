#ifndef ARRAY_H_INCLUDED
#define ARRAY_H_INCLUDED

#include <malloc.h>
///////////////////////////////////////////// 1 D ///////////////////////////////////////
double* Make1DArray(int m);
float *Make1DArray_float(int m);
double** Make1DArrayPointer(int m);
int*  Make1DArrayinteger(int m);
char*  Make1DString(int m);
///////////////////////////////////////////// 2 D ///////////////////////////////////////
double** Make2DArray(int m,int n);
double*** Make2DArrayPointer(int m,int n);
int double_2DArray_Pointer_Transfer_to_2DArray(int m, double ***a_Pointer, double **a);

int** Make2DArrayinteger(int m, int n);
int*** Make2DArrayintegerPointer(int m, int n);
int **int_2DArrayInteger_Pointer_Transfer_to_2DArrayInteger(int m, int ***a_Pointer);

///////////////////////////////////////////// 3 D ///////////////////////////////////////
double*** Make3DArray(int m, int n, int l);
double**** Make3DArrayPointer(int m, int n, int l);
double ***double_3DArray_Pointer_Transfer_to_3DArray(int m, int n, double ****a_Pointer);

int*** Make3DArrayinteger(int m, int n, int l);

char*** Make3DString(int m, int n, int l);

///////////////////////////////////////////// 4 D ///////////////////////////////////////
double**** Make4DArray(int m, int n, int l, int d);
int**** Make4DArrayinteger(int m, int n, int l, int d);

///////////////////////////////////////////// 5 D ///////////////////////////////////////
double****** Make5DArrayPointer(int m, int n, int l, int p, int q);

///////////////////////////////////////////// Free Part ///////////////////////////////////////
void free2DArray(double **a, int m);
void free2DArrayPointer(double ***a, int m);
void free3DArray(double ***a, int m, int n);
void free3DArrayPointer(double ****a, int m, int n);
void free4DArray(double ****a, int m, int n, int l);

void free2DArrayinteger(int **a, int m);
void free2DArrayintegerPointer(int ***a, int m);
void free3DArrayinteger(int ***a, int m, int n);
void free4DArrayinteger(int ****a, int m, int n, int l);
void free3DArrayString(char ***a, int m, int n);

#endif // ARRAY_H_INCLUDED
