#include <malloc.h>
#include <Array.h>

///////////////////////////////////////////// 1 D ///////////////////////////////////////
double* Make1DArray(int m)
{
	double *a;
	a=(double*)malloc(sizeof(double)*m);
	return a;
}

float* Make1DArray_float(int m)
{
	float *a;
	a=(float*)malloc(sizeof(float)*m);
	return a;
}

double** Make1DArrayPointer(int m)
{
	double **a;
	a=(double**)malloc(sizeof(double*)*m);
	return a;
}


int*  Make1DArrayinteger(int m)
{
	int *a;
	a=(int*)malloc(sizeof(int)*m);
	return a;
}

char*  Make1DString(int m)
{
	char *a;
	a=(char*)malloc(sizeof(char)*m);
	return a;
}

///////////////////////////////////////////// 2 D ///////////////////////////////////////
double** Make2DArray(int m,int n)
{
	int i;
	double **a;
	a= (double**)malloc(sizeof(double*)*m);
	for(i=0;i<m;i++)
		a[i]= (double*)malloc(sizeof(double)*n);
	return a;
}

double*** Make2DArrayPointer(int m,int n)
{
	int i;
	double ***a;
	a= (double***)malloc(sizeof(double**)*m);
	for(i=0;i<m;i++)
		a[i]= (double**)malloc(sizeof(double*)*n);
	return a;
}

int double_2DArray_Pointer_Transfer_to_2DArray(int m, double ***a_Pointer, double **a)
{
    int i;
    a=*a_Pointer;
    for(i=0;i<m;++i)
    {
        a[i]=*a_Pointer[i];
    }

    return 1;
}

int** Make2DArrayinteger(int m, int n)
{
    int i;
	int **a;
	a= (int**)malloc(sizeof(int*)*m);
	for(i=0;i<m;i++)
		a[i]= (int*)malloc(sizeof(int)*n);
	return a;
}

int*** Make2DArrayintegerPointer(int m, int n)
{
    int i;
	int ***a;
	a= (int***)malloc(sizeof(int**)*m);
	for(i=0;i<m;i++)
		a[i]= (int**)malloc(sizeof(int*)*n);
	return a;
}

int **int_2DArrayInteger_Pointer_Transfer_to_2DArrayInteger(int m, int ***a_Pointer)
{
    int **a, i;
    a=*a_Pointer;
    for(i=0;i<m;++i)
    {
        a[i]=*a_Pointer[i];
    }

    return a;
}
///////////////////////////////////////////// 3 D ///////////////////////////////////////
double*** Make3DArray(int m, int n, int l)
{
    int i, j;
	double ***a;
	a= (double***)malloc(sizeof(double**)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (double**)malloc(sizeof(double*)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (double*)malloc(sizeof(double)*l);
        }
    }

	return a;
}

double**** Make3DArrayPointer(int m, int n, int l)
{
    int i, j;
	double ****a;
	a= (double****)malloc(sizeof(double***)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (double***)malloc(sizeof(double**)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (double**)malloc(sizeof(double*)*l);
        }
    }

	return a;
}

double ***double_3DArray_Pointer_Transfer_to_3DArray(int m, int n, double ****a_Pointer)
{
    double ***a;
    int i, j;
    a=*a_Pointer;
    for(i=0;i<m;++i)
    {
        a[i]=*a_Pointer[i];
        for(j=0;j<n;++j)
        {
            a[i][j]=*a_Pointer[i][j];
        }
    }

    return a;
}

int*** Make3DArrayinteger(int m, int n, int l)
{
    int i, j;
	int ***a;
	a= (int***)malloc(sizeof(int**)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (int**)malloc(sizeof(int*)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (int*)malloc(sizeof(int)*l);
        }
    }

	return a;
}

char*** Make3DString(int m, int n, int l)
{
    int i, j;
	char ***a;
	a= (char***)malloc(sizeof(char**)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (char**)malloc(sizeof(char*)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (char*)malloc(sizeof(char)*l);
        }
    }

	return a;
}

///////////////////////////////////////////// 4 D ///////////////////////////////////////
double**** Make4DArray(int m, int n, int l, int d)
{
    int i, j, k;
	double ****a;
	a= (double****)malloc(sizeof(double***)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (double***)malloc(sizeof(double**)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (double**)malloc(sizeof(double*)*l);
            for(k=0;k<l;k++)
            {
                a[i][j][k]= (double*)malloc(sizeof(double)*d);
            }
        }
    }

	return a;
}


int**** Make4DArrayinteger(int m, int n, int l, int d)
{
    int i, j, k;
	int ****a;
	a= (int****)malloc(sizeof(int***)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (int***)malloc(sizeof(int**)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (int**)malloc(sizeof(int*)*l);
            for(k=0;k<l;k++)
            {
                a[i][j][k]= (int*)malloc(sizeof(int)*d);
            }
        }
    }

	return a;
}

///////////////////////////////////////////// 5 D ///////////////////////////////////////
double****** Make5DArrayPointer(int m, int n, int l, int p, int q)
{
    int i, j, k, r;
	double ******a;
	a= (double******)malloc(sizeof(double*****)*m);
	for(i=0;i<m;i++)
    {
        a[i]= (double*****)malloc(sizeof(double****)*n);
        for(j=0;j<n;j++)
        {
            a[i][j]= (double****)malloc(sizeof(double***)*l);
            for(k=0;k<l;k++)
            {
                a[i][j][k]= (double***)malloc(sizeof(double**)*p);
                for(r=0;r<p;++p)
                {
                    a[i][j][k][r]= (double**)malloc(sizeof(double*)*q);
                }
            }
        }
    }

	return a;
}
///////////////////////////////////////////// Free Part ///////////////////////////////////////

void free2DArray(double **a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

void free2DArrayPointer(double ***a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

void free2DArrayinteger(int **a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

void free2DArrayintegerPointer(int ***a, int m)
{
	int i;
	for( i=0; i<m; i++ )
		free(a[i]);
	free(a);
}

void free3DArray(double ***a, int m, int n)
{
	int i, j;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}

void free3DArrayPointer(double ****a, int m, int n)
{
	int i, j;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}

void free4DArray(double ****a, int m, int n, int l)
{
	int i, j, k;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                free(a[i][j][k]);
            }
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}

void free3DArrayinteger(int ***a, int m, int n)
{
	int i, j;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}

void free4DArrayinteger(int ****a, int m, int n, int l)
{
	int i, j, k;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                free(a[i][j][k]);
            }
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}


void free3DArrayString(char ***a, int m, int n)
{
	int i, j;
	for( i=0; i<m; i++ )
    {
        for(j=0;j<n;j++)
        {
            free(a[i][j]);
        }
        free(a[i]);
    }

	free(a);
}
