#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Array.h>
#include <Initial_State.h>
#include <Variable.h>

#define pi   3.14159265358

using namespace std;

int Initial_Magnetization_FE(int m, int n, int l, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }
    return 1;
}

int Initial_Magnetization_Center_Circles(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k, xc, yc, rn, xn, yn;
    double D, Theta, Phi;
    xc = (m-1)/2; yc = (n-1)/2;
    xn = xc-1;    yn = yc-1;

    if(yn<xn)
    {
        rn = yn;
    }else
    {
        rn = xn;
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }

    if (IPN->Energy_Coefficient[DM_Bloch] > 0) //Bloch
    {
        Theta = pi / rn;
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    D = sqrt(pow((i-xc),2)+pow((j-yc),2));
                    if(D<=rn)
                    {
                        Phi = atan2((j - yc), (i - xc));
                        Mz[i][j][k] = -cos(D * Theta);
                        Mx[i][j][k] = -sin(D * Theta)*sin(Phi);
                        My[i][j][k] = sin(D * Theta)*cos(Phi);
                    }
                }
            }
        }
    }
    else if(IPN->Energy_Coefficient[DM_Neel] > 0)//Neel
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    D = sqrt(pow((i-xc),2)+pow((j-yc),2));
                    if(D<=rn)
                    {
                        Theta = pi * (1 - D / rn);
                        Phi = atan2(j - yc, i - xc);
                        Mz[i][j][k] = cos(Theta);
                        Mx[i][j][k] = -sin(Theta) * cos(Phi);
                        My[i][j][k] = -sin(Theta) * sin(Phi);
                    }
                }
            }
        }
    }
    return 1;
}

// int Initial_Magnetization_Multiple_Circles(int m, int n, int l, double ***Mx, double ***My, double ***Mz)
// {
//     int i, j, k, xc, yc, rn, xn, yn, Layer_number, Surround_number, Center_number, p, q, Cn;
//     int **Circle_Center;

//     double D, Theta, r, Phi;

//     Layer_number = 3; Surround_number = 5;
//     xc = (m-1)/2; yc = (n-1)/2;
//     xn = xc-1;    yn = yc-1;

//     Center_number = 1+Layer_number*Surround_number;
//     Circle_Center = Make2DArrayinteger(Center_number,2);

//     Circle_Center[0][0] = xn;   Circle_Center[0][1] = yn;

//     if(yn<xn)
//     {
//         rn = yn;
//     }else
//     {
//         rn = xn;
//     }



//     Cn = 1;
//     for(p=1;p<Layer_number;++p)
//     {
//         for(q=0;q<Surround_number;++q)
//         {
//             r = rn/(2*Layer_number-1)*(2*p);
//             Theta = 2*pi/Surround_number*q;

//             Cn = Cn+1;
//             Circle_Center[Cn][0] = xn+r*cos(Theta);     Circle_Center[Cn][1] = yn+r*sin(Theta);
//         }
//     }

//     for(i=0;i<m;i++)
//     {
//         for(j=0;j<n;j++)
//         {
//             for(k=0;k<l;k++)
//             {
//                 Mx[i][j][k]=0;
//                 My[i][j][k]=0;
//                 Mz[i][j][k]=1;
//             }
//         }
//     }

//     r = rn/(2*Layer_number-1);

//     for(k=0;k<l;++k)
//     {
//         for(i=0;i<m;++i)
//         {
//             for(j=0;j<n;++j)
//             {
//                 for(p=0;p<Center_number;++p)
//                 {
//                     D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
//                     if(D<=r)
//                     {
//                         //Mz[i][j][k] = -(1-D/r)+D/r;
//                         Phi = (1 - D / r) * pi;
//                         Mz[i][j][k] = cos(Phi);
//                         Mx[i][j][k] = sin(Phi)/sqrt(2);
//                         My[i][j][k] = sin(Phi)/sqrt(2);
//                     }
//                 }
//             }
//         }
//     }
//     return 1;
// }

int Initial_Magnetization_Multiple_Circles(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz, int Surround_number)
{
    int i, j, k, xc, yc, rn, xn, yn, Layer_number, Center_number, p, q, Cn;
    int **Circle_Center;

    double D, Theta, Theta1, r, Phi;

    Layer_number = 2; //Surround_number = 5;
    xc = (m-1)/2; yc = (n-1)/2;
    xn = xc-1;    yn = yc-1;

    Center_number = 1+Layer_number*Surround_number;
    Circle_Center = Make2DArrayinteger(Center_number,2);

    Circle_Center[0][0] = xn;   Circle_Center[0][1] = yn;

    if(yn<xn)
    {
        rn = yn;
    }else
    {
        rn = xn;
    }



    Cn = 1;
    for(p=1;p<Layer_number;++p)
    {
        for(q=0;q<Surround_number;++q)
        {
            r = rn/(2*Layer_number-1)*(2*p);
            Theta = 2*pi/Surround_number*q;

            Cn = Cn+1;
            Circle_Center[Cn][0] = xn+r*cos(Theta);     Circle_Center[Cn][1] = yn+r*sin(Theta);
        }
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }

    r = rn/(2*Layer_number-1);
    if (IPN->Energy_Coefficient[DM_Bloch] > 0) //Bloch
    {
        Theta = pi / r;
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    for(p=0;p<Center_number;++p)
                    {
                        D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
                        if(D<=r)
                        {
                            Phi = atan2((j-Circle_Center[p][1]), (i-Circle_Center[p][0]));
                            Mz[i][j][k] = -cos(D * Theta);
                            Mx[i][j][k] = -sin(D * Theta)*sin(Phi);
                            My[i][j][k] = sin(D * Theta)*cos(Phi);
                        }
                    }
                }
            }
        }
    }else if (IPN->Energy_Coefficient[DM_Neel] > 0) //Neel
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    for(p=0;p<Center_number;++p)
                    {
                        D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
                        if(D<=r)
                        {
                            //Mz[i][j][k] = -(1-D/r)+D/r;
                            Theta1 = pi * (1 - D / r);
                            Phi = atan2(j - Circle_Center[p][1], i - Circle_Center[p][0]);
                            Mz[i][j][k] = cos(Theta1);
                            Mx[i][j][k] = -sin(Theta1) * cos(Phi);
                            My[i][j][k] = -sin(Theta1) * sin(Phi);
                        }
                    }
                }
            }
        }
    }
    free2DArrayinteger(Circle_Center, Center_number);
    return 1;
}

int Initial_Magnetization_Multiple_Circles_Hollow(int m, int n, int l, double ***Mx, double ***My, double ***Mz, int Surround_number)
{
    int i, j, k, xc, yc, rn, xn, yn, Center_number, p, q, Cn;
    int **Circle_Center;

    double D, Theta, Theta1, r, Phi;

    //Surround_number = 5;
    xc = (m-1)/2; yc = (n-1)/2;
    xn = xc-1;    yn = yc-1;

    Center_number = Surround_number;
    Circle_Center = Make2DArrayinteger(Center_number,2);

    //Circle_Center[0][0] = xn;   Circle_Center[0][1] = yn;

    if(yn<xn)
    {
        rn = yn;
    }else
    {
        rn = xn;
    }



    Cn = 0;
    // for(p=1;p<Layer_number;++p)
    // {
    for(q=0;q<Surround_number;++q)
    {
        r = rn*(1-1/(1/sin(pi/Surround_number)+1));
        Theta = 2*pi/Surround_number*q;

        
        Circle_Center[Cn][0] = xn+r*cos(Theta);     Circle_Center[Cn][1] = yn+r*sin(Theta);
        Cn = Cn+1;
    }
    // }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }

    r = rn/(1/sin(pi/Surround_number)+1);
//    r = rn*sin(pi/Surround_number)/2;


    for(k=0;k<l;++k)
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(p=0;p<Center_number;++p)
                {
                    D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
                    if(D<=r)
                    {
                        //Mz[i][j][k] = -(1-D/r)+D/r;
                        Theta1 = pi * (1 - D / r);
                        Phi = atan2(j - Circle_Center[p][1], i - Circle_Center[p][0]);
                        Mz[i][j][k] = cos(Theta1);
                        Mx[i][j][k] = -sin(Theta1) * cos(Phi);
                        My[i][j][k] = -sin(Theta1) * sin(Phi);
                    }
                }
            }
        }
    }
    free2DArrayinteger(Circle_Center, Center_number);
    return 1;
}

/*
int InitialM_3D(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                double *MxV, double *MyV, double *MzV, Boundary_Shape *BS)
{
    int i, j, j_sym, k, v, m, n, l, zc, Version;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV;
    enum Initial_Mode_NMD{Initial_M_FE=-2, Initial_Helical, Read_M_3D, M_SKN1, M_SKN2, M_SKN3, M_SKN4, M_SKN5, M_SKN5_C, M_SKN6_C, M_SKN7_C, M_SKN8_C, M_SKN9_C, M_SKN10, M_SKN11};
    Initial_Mode_NMD IMN;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);
    Version = IPN->Version;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }

    IMN = (Initial_Mode_NMD) IPN->InitialMode;

    switch(IMN)
    {
        case Read_M_3D :
        {
            sprintf(filename_MxV,"Version %d\\M3D\\MxV.dat",Version);
            sprintf(filename_MyV,"Version %d\\M3D\\MyV.dat",Version);
            sprintf(filename_MzV,"Version %d\\M3D\\MzV.dat",Version);
            fp_MxV = fopen(filename_MxV,"r+");
            fp_MyV = fopen(filename_MyV,"r+");
            fp_MzV = fopen(filename_MzV,"r+");

            for(v=0;v<BS->Virtual_Points;++v)
            {
                fscanf(fp_MxV,"%lf",&MxV[v]);
                fscanf(fp_MyV,"%lf",&MyV[v]);
                fscanf(fp_MzV,"%lf",&MzV[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);

            for(k=0;k<l;k++)
            {
                sprintf(filename_Mx,"Version %d\\M3D\\z=%d\\Mx.dat",Version,k-zc);
                sprintf(filename_My,"Version %d\\M3D\\z=%d\\My.dat",Version,k-zc);
                sprintf(filename_Mz,"Version %d\\M3D\\z=%d\\Mz.dat",Version,k-zc);
                fp_Mx = fopen(filename_Mx,"r+");
                fp_My = fopen(filename_My,"r+");
                fp_Mz = fopen(filename_Mz,"r+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fscanf(fp_Mx,"%lf",&Mx[i][j][k]);
                        fscanf(fp_My,"%lf",&My[i][j][k]);
                        fscanf(fp_Mz,"%lf",&Mz[i][j][k]);
                    }
                }
                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
            }
            break;
        }

        case M_SKN1:
        {
            Initial_Magnetization_Center_Circles(m,n,l,Mx,My,Mz);
            break;
        }

        case M_SKN2:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,2);
            break;
        }

        case M_SKN3:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,3);
            break;
        }

        case M_SKN4:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,4);
            break;
        }

        case M_SKN5:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,5);
            break;
        }

        case M_SKN5_C:
        {
            Initial_Magnetization_Multiple_Circles(m,n,l,Mx,My,Mz,5);
            break;
        }

        case M_SKN6_C:
        {
            Initial_Magnetization_Multiple_Circles(m,n,l,Mx,My,Mz,6);
            break;
        }

        case M_SKN7_C:
        {
            Initial_Magnetization_Multiple_Circles(m,n,l,Mx,My,Mz,7);
            break;
        }

        case M_SKN8_C:
        {
            Initial_Magnetization_Multiple_Circles(m,n,l,Mx,My,Mz,8);
            break;
        }

        case M_SKN9_C:
        {
            Initial_Magnetization_Multiple_Circles(m,n,l,Mx,My,Mz,9);
            break;
        }

        case M_SKN10:
        {
            Initial_Magnetization_10_Circles(IPN, Mx, My, Mz);
            break;
        }

        case M_SKN11:
        {
            Initial_Magnetization_11_Circles(IPN, Mx, My, Mz);
            break;
        }
        
        case Initial_M_FE:
        {
            Initial_Magnetization_FE(m,n,l,Mx,My,Mz);
            break;
        }

        case Initial_Helical:
        {
            Initial_Magnetization_Helical(IPN,Mx,My,Mz);
            break;
        }
        default :
            break;
    }

    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    return 1;
}
*/

int InitialM_3D(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                double *MxV, double *MyV, double *MzV, Boundary_Shape *BS)
{
    int i, j, k, v, m, n, l, zc, Version;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV;
    Initial_Mode_NMD IMN;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);
    Version = IPN->Version;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }

    IMN = (Initial_Mode_NMD) IPN->InitialMode;

    switch(IMN)
    {
        case Read_M_3D :
        {
            sprintf(filename_MxV,"Version %d\\M3D\\MxV.dat",Version);
            sprintf(filename_MyV,"Version %d\\M3D\\MyV.dat",Version);
            sprintf(filename_MzV,"Version %d\\M3D\\MzV.dat",Version);
            fp_MxV = fopen(filename_MxV,"r+");
            fp_MyV = fopen(filename_MyV,"r+");
            fp_MzV = fopen(filename_MzV,"r+");

            for(v=0;v<BS->Virtual_Points;++v)
            {
                fscanf(fp_MxV,"%lf",&MxV[v]);
                fscanf(fp_MyV,"%lf",&MyV[v]);
                fscanf(fp_MzV,"%lf",&MzV[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);

            for(k=0;k<l;k++)
            {
                sprintf(filename_Mx,"Version %d\\M3D\\z=%d\\Mx.dat",Version,k-zc);
                sprintf(filename_My,"Version %d\\M3D\\z=%d\\My.dat",Version,k-zc);
                sprintf(filename_Mz,"Version %d\\M3D\\z=%d\\Mz.dat",Version,k-zc);
                fp_Mx = fopen(filename_Mx,"r+");
                fp_My = fopen(filename_My,"r+");
                fp_Mz = fopen(filename_Mz,"r+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fscanf(fp_Mx,"%lf",&Mx[i][j][k]);
                        fscanf(fp_My,"%lf",&My[i][j][k]);
                        fscanf(fp_Mz,"%lf",&Mz[i][j][k]);
                    }
                }
                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
            }
            break;
        }

        case Read_M_2D:
        {
            double ***MxT, ***MyT, ***MzT;
            //Transposition
            MxT=Make3DArray(m,n,l);
            MyT=Make3DArray(m,n,l);
            MzT=Make3DArray(m,n,l);

            for(k=0;k<l;k++)
            {
                sprintf(filename_Mx,"M2D\\Mx.dat");
                sprintf(filename_My,"M2D\\My.dat");
                sprintf(filename_Mz,"M2D\\Mz.dat");
                fp_Mx = fopen(filename_Mx,"r+");
                fp_My = fopen(filename_My,"r+");
                fp_Mz = fopen(filename_Mz,"r+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fscanf(fp_Mx,"%lf",&MxT[i][j][k]);
                        fscanf(fp_My,"%lf",&MyT[i][j][k]);
                        fscanf(fp_Mz,"%lf",&MzT[i][j][k]);
                    }
                }
                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);


                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        Mx[i][j][k]=MxT[j][i][k];
                        My[i][j][k]=MyT[j][i][k];
                        Mz[i][j][k]=MzT[j][i][k];
                    }
                }
            }
            break;
        }

        case Initial_M_FE:
        {
            Initial_Magnetization_FE(m,n,l,Mx,My,Mz);
            break;
        }


        case Read_M_Process:
        {
            sprintf(filename_MxV,"Version %d\\Process\\MxV.dat",Version);
            sprintf(filename_MyV,"Version %d\\Process\\MyV.dat",Version);
            sprintf(filename_MzV,"Version %d\\Process\\MzV.dat",Version);
            fp_MxV = fopen(filename_MxV,"r+");
            fp_MyV = fopen(filename_MyV,"r+");
            fp_MzV = fopen(filename_MzV,"r+");

            for(v=0;v<BS->Virtual_Points;++v)
            {
                fscanf(fp_MxV,"%lf",&MxV[v]);
                fscanf(fp_MyV,"%lf",&MyV[v]);
                fscanf(fp_MzV,"%lf",&MzV[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);

            for(k=0;k<l;k++)
            {
            //sprintf(filename,"Version %d\\M3D",IPP->Version);
                sprintf(filename_Mx,"Version %d\\Process\\z=%d\\Mx.dat",Version,k-zc);
                sprintf(filename_My,"Version %d\\Process\\z=%d\\My.dat",Version,k-zc);
                sprintf(filename_Mz,"Version %d\\Process\\z=%d\\Mz.dat",Version,k-zc);
                fp_Mx=fopen(filename_Mx,"r+");
                fp_My=fopen(filename_My,"r+");
                fp_Mz=fopen(filename_Mz,"r+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fscanf(fp_Mx,"%lf",&Mx[i][j][k]);
                        fscanf(fp_My,"%lf",&My[i][j][k]);
                        fscanf(fp_Mz,"%lf",&Mz[i][j][k]);
                    }
                }
                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
            }
            break;
        }

        case Initial_LLG:
        {
            Initial_Magnetization_LLG(IPN,m,n,l,Mx,My,Mz);
            break;
        }
        case Initial_Helical:
        {
            Initial_Magnetization_Helical(IPN,Mx,My,Mz);
            break;
        }

        case IMN_1Circle:
        {
            Initial_Magnetization_Center_Circles(IPN,m,n,l,Mx,My,Mz);
            break;
        }

        case IMN_2Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,2);
            break;
        }

        case IMN_3Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,3);
            break;
        }

        case IMN_4Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,4);
            break;
        }

        case IMN_5Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,5);
            break;
        }

        case IMN_6Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,5);
            break;
        }

        case IMN_7Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,6);
            break;
        }

        case IMN_8Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,7);
            break;
        }

        case IMN_9Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,8);
            break;
        }

        case IMN_10Circles:
        {
            Initial_Magnetization_10_Circles(IPN, Mx, My, Mz);
            break;
        }

        case IMN_11Circles:
        {
            Initial_Magnetization_11_Circles(IPN, Mx, My, Mz);
            break;
        }
        

        default :
            break;
    }


    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    return 1;
}

int InitialM_3D_SKN(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                    double *MxV, double *MyV, double *MzV, Boundary_Shape *BS, char *dir_SKN)
{
    int i, j, k, v, m, n, l, zc, Version;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV, *String_Temp;
    Initial_Mode_NMD IMN;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);
    String_Temp  = Make1DString(500);
    Version = IPN->Version;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }

    IMN = (Initial_Mode_NMD) IPN->InitialMode;

    switch(IMN)
    {
        case Read_M_3D :
        {
            sprintf(String_Temp,"Version %d\\M3D\\MxV.dat",Version);
            sprintf(filename_MxV, "%s\\%s", dir_SKN, String_Temp);

            sprintf(String_Temp,"Version %d\\M3D\\MyV.dat",Version);
            sprintf(filename_MyV, "%s\\%s", dir_SKN, String_Temp);

            sprintf(String_Temp,"Version %d\\M3D\\MzV.dat",Version);
            sprintf(filename_MzV, "%s\\%s", dir_SKN, String_Temp);

            fp_MxV = fopen(filename_MxV,"r+");
            fp_MyV = fopen(filename_MyV,"r+");
            fp_MzV = fopen(filename_MzV,"r+");

            for(v=0;v<BS->Virtual_Points;++v)
            {
                fscanf(fp_MxV,"%lf",&MxV[v]);
                fscanf(fp_MyV,"%lf",&MyV[v]);
                fscanf(fp_MzV,"%lf",&MzV[v]);
            }
            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);

            for(k=0;k<l;k++)
            {
                sprintf(String_Temp,"Version %d\\M3D\\z=%d\\Mx.dat",Version,k-zc);
                sprintf(filename_Mx, "%s\\%s", dir_SKN, String_Temp);

                sprintf(String_Temp,"Version %d\\M3D\\z=%d\\My.dat",Version,k-zc);
                sprintf(filename_My, "%s\\%s", dir_SKN, String_Temp);

                sprintf(String_Temp,"Version %d\\M3D\\z=%d\\Mz.dat",Version,k-zc);
                sprintf(filename_Mz, "%s\\%s", dir_SKN, String_Temp);

                fp_Mx = fopen(filename_Mx,"r+");
                fp_My = fopen(filename_My,"r+");
                fp_Mz = fopen(filename_Mz,"r+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fscanf(fp_Mx,"%lf",&Mx[i][j][k]);
                        fscanf(fp_My,"%lf",&My[i][j][k]);
                        fscanf(fp_Mz,"%lf",&Mz[i][j][k]);
                    }
                }
                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
            }
            break;
        }

        case IMN_1Circle:
        {
            Initial_Magnetization_Center_Circles(IPN,m,n,l,Mx,My,Mz);
            break;
        }

        case IMN_2Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,2);
            break;
        }

        case IMN_3Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,3);
            break;
        }

        case IMN_4Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,4);
            break;
        }

        case IMN_5Circles:
        {
            Initial_Magnetization_Multiple_Circles_Hollow(m,n,l,Mx,My,Mz,5);
            break;
        }

        case IMN_6Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,5);
            break;
        }

        case IMN_7Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,6);
            break;
        }

        case IMN_8Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,7);
            break;
        }

        case IMN_9Circles:
        {
            Initial_Magnetization_Multiple_Circles(IPN,m,n,l,Mx,My,Mz,8);
            break;
        }

        case IMN_10Circles:
        {
            Initial_Magnetization_10_Circles(IPN, Mx, My, Mz);
            break;
        }

        case IMN_11Circles:
        {
            Initial_Magnetization_11_Circles(IPN, Mx, My, Mz);
            break;
        }
        
        case Initial_M_FE:
        {
            Initial_Magnetization_FE(m,n,l,Mx,My,Mz);
            break;
        }

        case Initial_Helical:
        {
            Initial_Magnetization_Helical(IPN,Mx,My,Mz);
            break;
        }
        default :
            break;
    }

    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    free(String_Temp);
    return 1;
}

int Trim_Boundary_M3D_VirtualToBoundary(Boundary_Shape *BS, Calculate_Points_Parameters *CPP, Input_Parameter_NMD *IPN,
                                        double ***Mx, double ***My, double ***Mz, double *MxV, double *MyV, double *MzV)
{
    int i, j, k, m, n, l, v, N, p;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;
    m = 2*BS->xn+3; n = 2*BS->yn+3;
    if(BS->zn==0)
    {
        l=1;
    }else
    {
        l=2*BS->zn+3;
    }

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(BS->Boundary[i][j][k]==0)
                {
                    Mx[i][j][k]=0;
                    My[i][j][k]=0;
                    Mz[i][j][k]=0;
                }
            }
        }
    }

    if(IPN->Boundary_Shape_Type==0)
    {
        for(v=1;v<BS->Virtual_Points;++v)
        {
            N = CPP->Virtual_Position[v];
            for(p=0;p<7;++p)
            {
                LPP=(Local_Point_Postion) p;
                switch(LPP)
                {
                    case North:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            i = i+1;
                            Mx[i][j][k] = MxV[v];
                            My[i][j][k] = MyV[v];
                            Mz[i][j][k] = MzV[v];
                        }
                        break;
                    }

                    case South:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            i = i-1;
                            Mx[i][j][k] = MxV[v];
                            My[i][j][k] = MyV[v];
                            Mz[i][j][k] = MzV[v];
                        }
                        break;
                    }

                    case West:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            j = j+1;
                            Mx[i][j][k] = MxV[v];
                            My[i][j][k] = MyV[v];
                            Mz[i][j][k] = MzV[v];
                        }
                        break;
                    }

                    case East:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            j = j-1;
                            Mx[i][j][k] = MxV[v];
                            My[i][j][k] = MyV[v];
                            Mz[i][j][k] = MzV[v];
                        }
                        break;
                    }

                    default:
                        break;
                }
            }
        }
    }

    return 1;
}

int Trim_Boundary_M3D_BoundaryToVirtual(Boundary_Shape *BS, Calculate_Points_Parameters *CPP, Input_Parameter_NMD *IPN,
                                        double ***Mx, double ***My, double ***Mz, double *MxV, double *MyV, double *MzV)
{
    int i, j, k, v, N, p;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;

    if(IPN->Boundary_Shape_Type==0)
    {
        for(v=1;v<BS->Virtual_Points;++v)
        {
            N = CPP->Virtual_Position[v];
            for(p=0;p<7;++p)
            {
                LPP=(Local_Point_Postion) p;
                switch(LPP)
                {
                    case North:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            i = i+1;
                            MxV[v] = Mx[i][j][k];
                            MyV[v] = My[i][j][k];
                            MzV[v] = Mz[i][j][k];
                        }
                        break;
                    }

                    case South:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            i = i-1;
                            MxV[v] = Mx[i][j][k];
                            MyV[v] = My[i][j][k];
                            MzV[v] = Mz[i][j][k];
                        }
                        break;
                    }

                    case West:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            j = j+1;
                            MxV[v] = Mx[i][j][k];
                            MyV[v] = My[i][j][k];
                            MzV[v] = Mz[i][j][k];
                        }
                        break;
                    }

                    case East:
                    {
                        if(CPP->Local_Points_1D[N][p][5]==1)
                        {
                            i = CPP->Local_Points_1D[N][p][1];
                            j = CPP->Local_Points_1D[N][p][2];
                            k = CPP->Local_Points_1D[N][p][3];

                            j = j-1;
                            MxV[v] = Mx[i][j][k];
                            MyV[v] = My[i][j][k];
                            MzV[v] = Mz[i][j][k];
                        }
                        break;
                    }

                    default:
                        break;
                }
            }
        }
    }

    return 1;
}

int Initial_Magnetization_LLG(Input_Parameter_NMD *IPN, int m, int n, int l, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k, xc, yc, rn;
    double D, Theta, Mmax, Phi;
    Mmax = 0.65;
    xc = (m-1)/2; yc = (n-1)/2;   
    //xc = (m-1)/2-200; yc = (n-1)/2;
    //xc = yc;

    //rn = 25;
    rn = n*1/8;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=Mmax;
            }
        }
    }
    
    if (IPN->Energy_Coefficient[DM_Bloch] > 0) //Bloch
    {
        Theta = pi / rn;
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    D = sqrt(pow((i-xc),2)+pow((j-yc),2));
                    if(D<=rn)
                    {
                        Phi = atan2((j - yc), (i - xc));
                        Mz[i][j][k] = -Mmax * cos(D * Theta);
                        Mx[i][j][k] = -Mmax * sin(D * Theta)*sin(Phi);
                        My[i][j][k] = Mmax * sin(D * Theta)*cos(Phi);
                    }
                }
            }
        }
    }else if (IPN->Energy_Coefficient[DM_Neel] > 0) //Neel
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    D = sqrt(pow((i-xc),2)+pow((j-yc),2));
                    if(D<=rn)
                    {
                        Theta = pi * (1 - D / rn);
                        Phi = atan2(j - yc, i - xc);
                        Mz[i][j][k] = Mmax * cos(Theta);
                        Mx[i][j][k] = -Mmax * sin(Theta) * cos(Phi);
                        My[i][j][k] = -Mmax * sin(Theta) * sin(Phi);
                    }
                }
            }
        }
    }
    
    return 1;
}

int Initial_Magnetization_Helical(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz)
{
    int m, n, l;
    m = 2 * IPN->xn + 3;
    n = 2 * IPN->yn + 3;
    if(IPN->zn==0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }

    int i, j, k, xc, yc;
    double x, y;
    //double D, Theta, MInplan, Mmax, Phi;
    xc = (m-1)/2;   yc = (n-1)/2;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                x = (i - xc) * IPN->dxy[0];
                y = (j - yc) * IPN->dxy[1];
                //Mx[i][j][k]=sin(2*pi*x);
                My[i][j][k]=0;
                //Mz[i][j][k]=cos(2*pi*x);
                Mx[i][j][k]=sin(x);
                My[i][j][k]=0;
                Mz[i][j][k]=cos(x);
            }
        }
    }
    return 1;
}

int Initial_Magnetization_10_Circles(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k, xn, yn, rn, Center_number, p;
    //int xc, yc;
    xn = IPN->xn;   yn=IPN->yn;
    int m, n, l;
    m = 2 * xn + 3;
    n = 2 * yn + 3;
    if(IPN->zn==0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }

    
    double **Circle_Center;
    double D, Theta, Theta1, r, rr, Phi;

    //Layer_number = 2; //Surround_number = 5;
    //xc = (m-1)/2; yc = (n-1)/2;
    //xn = xc-1;    yn = yc-1;

    Center_number = 10;
    Circle_Center = Make2DArray(Center_number,2);

    //Circle_Center[0][0] = xn;   Circle_Center[0][1] = yn;

    if(yn<xn)
    {
        rn = yn;
    }else
    {
        rn = xn;
    }
    rr = rn * 0.26225;
    
    for(p=0;p<8;++p)
    {
            r = rn-rr;
            Theta = 3*pi/2-3.5 * 2 * asin(rr / r) + p*2*asin(rr/r);
            //printf("Theta=%lf\n", Theta);
            //Cn = Cn+1;
            Circle_Center[p][0] = xn+r*cos(Theta);     Circle_Center[p][1] = yn+r*sin(Theta);
    }
    //Theta = pi/2;
    //r = sqrt((rn - rr) * (rn - rr) - rr * rr) - rr * sqrt(3);
    //Circle_Center[8][0] = xn+r*cos(Theta);     Circle_Center[8][1] = yn+r*sin(Theta);

    Circle_Center[8][0] = (Circle_Center[3][0]+Circle_Center[4][0])/2;     
    Circle_Center[8][1] = Circle_Center[3][1]+sqrt(3)*rr;

    Circle_Center[9][0] = (Circle_Center[3][0]+Circle_Center[4][0])/2;     
    Circle_Center[9][1] = Circle_Center[3][1]+sqrt(3)*rr+2*rr;


    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }

    r = rn * 0.26225;

    for(k=0;k<l;++k)
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(p=0;p<Center_number;++p)
                {
                    D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
                    if(D<=r)
                    {
                        //Mz[i][j][k] = -(1-D/r)+D/r;
                        Theta1 = pi * (1 - D / r);
                        Phi = atan2(j - Circle_Center[p][1], i - Circle_Center[p][0]);
                        Mz[i][j][k] = cos(Theta1);
                        Mx[i][j][k] = -sin(Theta1) * cos(Phi);
                        My[i][j][k] = -sin(Theta1) * sin(Phi);
                    }
                }
            }
        }
    }
    free2DArray(Circle_Center, 10);
    return 1;
}

int Initial_Magnetization_11_Circles(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k, xn, yn, rn, Center_number, p;
    xn = IPN->xn;   yn=IPN->yn;
    int m, n, l;
    m = 2 * xn + 3;
    n = 2 * yn + 3;
    if(IPN->zn==0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }

    
    double **Circle_Center;
    double D, Theta, Theta1, r, rr, Phi;

    //Layer_number = 2; //Surround_number = 5;
    //xc = (m-1)/2; yc = (n-1)/2;
    //xn = xc-1;    yn = yc-1;

    Center_number = 11;
    Circle_Center = Make2DArray(Center_number,2);

    //Circle_Center[0][0] = xn;   Circle_Center[0][1] = yn;

    if(yn<xn)
    {
        rn = yn;
    }else
    {
        rn = xn;
    }
    rr = rn * 0.25485;
    
    for(p=0;p<8;++p)
    {
            r = rn-rr;
            Theta = 3*pi/2-3.5 * 2 * asin(rr / r) + p*2*asin(rr/r);
            //printf("Theta=%lf\n", Theta);
            //Cn = Cn+1;
            Circle_Center[p][0] = xn+r*cos(Theta);     Circle_Center[p][1] = yn+r*sin(Theta);
    }
    //Theta = pi/2;
    //r = sqrt((rn - rr) * (rn - rr) - rr * rr) - rr * sqrt(3);
    //Circle_Center[8][0] = xn+r*cos(Theta);     Circle_Center[8][1] = yn+r*sin(Theta);

    Circle_Center[8][0] = (Circle_Center[3][0]+Circle_Center[4][0])/2;
    Circle_Center[8][1] = Circle_Center[0][1] - sqrt(4 * rr * rr - (Circle_Center[0][0] - Circle_Center[7][0]) * (Circle_Center[0][0] - Circle_Center[7][0]) / 4);

    Circle_Center[9][0] = Circle_Center[3][0];     
    Circle_Center[9][1] = Circle_Center[8][1]-sqrt(3)*rr;

    Circle_Center[10][0] = Circle_Center[4][0];     
    Circle_Center[10][1] = Circle_Center[8][1]-sqrt(3)*rr;


    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                Mx[i][j][k]=0;
                My[i][j][k]=0;
                Mz[i][j][k]=1;
            }
        }
    }

    r = rn * 0.25485;

    for(k=0;k<l;++k)
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(p=0;p<Center_number;++p)
                {
                    D = sqrt(pow((i-Circle_Center[p][0]),2)+pow((j-Circle_Center[p][1]),2));
                    if(D<=r)
                    {
                        //Mz[i][j][k] = -(1-D/r)+D/r;
                        Theta1 = pi * (1 - D / r);
                        Phi = atan2(j - Circle_Center[p][1], i - Circle_Center[p][0]);
                        Mz[i][j][k] = cos(Theta1);
                        Mx[i][j][k] = -sin(Theta1) * cos(Phi);
                        My[i][j][k] = -sin(Theta1) * sin(Phi);
                    }
                }
            }
        }
    }
    free2DArray(Circle_Center, 10);
    return 1;
}

int Initial_Magnetization_LLG_BC(int m, int n, int l, double ***Mx, double ***My, double ***Mz)
{
    int i, j, k, xc, yc;
    double D, Mmax, Mmin, rn;
    int XB = 7 * m / 8, XE = m, YB = 3*n / 4, YE = n;

    xc = (m-1)/2; yc = (n-1)/2;
    xc = yc;

    if (m/16<n/6)
    {
        rn = n / 6;
    }
    else
    {
        rn = m / 16;
    }
    
    //rn = m*1/8;
    Mmin = Mz[0][0][0];
    Mmax = Mz[0][0][0];
    for(i=XB;i<XE;i++)
    {
        for(j=YB;j<YE;j++)
        {
            for(k=0;k<l;k++)
            {
                if (Mz[i][j][k]<Mmin)
                {
                    Mmin = Mz[i][j][k];
                    xc = i;
                    yc = j;
                }
            }
        }
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                if (Mz[i][j][k]>Mmax)
                {
                    Mmax = Mz[i][j][k];
                }
            }
        }
    }


    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                D = sqrt(pow((i-xc),2)+pow((j-yc),2));
                //if (D>rn||Mz[i][j][k]>0.45)
                if (D>rn)
                {
                    Mx[i][j][k]=0;
                    My[i][j][k]=0;
                    Mz[i][j][k]=Mmax;
                }
            }
        }
    }
    return 1;
}

int InitialM_Target(double ***Mx, double ***My, double ***Mz,
                    double *MxV, double *MyV, double *MzV, Boundary_Shape *BS)
{
    int i, j, v, m, n;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);

    m=2*BS->xn+3;   n=2*BS->yn+3;


    sprintf(filename_Mx, "M_Target\\Mx.dat");
    sprintf(filename_My, "M_Target\\My.dat");
    sprintf(filename_Mz, "M_Target\\Mz.dat");
    sprintf(filename_MxV, "M_Target\\MxV.dat");
    sprintf(filename_MyV, "M_Target\\MyV.dat");
    sprintf(filename_MzV, "M_Target\\MzV.dat");

    fp_MxV = fopen(filename_MxV,"r+");
    fp_MyV = fopen(filename_MyV,"r+");
    fp_MzV = fopen(filename_MzV,"r+");

    for(v=0;v<BS->Virtual_Points;++v)
    {
        fscanf(fp_MxV,"%lf",&MxV[v]);
        fscanf(fp_MyV,"%lf",&MyV[v]);
        fscanf(fp_MzV,"%lf",&MzV[v]);
    }
    fclose(fp_MxV);
    fclose(fp_MyV);
    fclose(fp_MzV);

    fp_Mx = fopen(filename_Mx,"r+");
    fp_My = fopen(filename_My,"r+");
    fp_Mz = fopen(filename_Mz,"r+");

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            fscanf(fp_Mx,"%lf",&Mx[i][j][0]);
            fscanf(fp_My,"%lf",&My[i][j][0]);
            fscanf(fp_Mz,"%lf",&Mz[i][j][0]);
        }
    }
    fclose(fp_Mx);
    fclose(fp_My);
    fclose(fp_Mz);

    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    return 1;
}