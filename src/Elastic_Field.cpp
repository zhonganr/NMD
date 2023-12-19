#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex>
#include <Elastic_Field.h>
#include <Variable.h>
#include <Array.h>

using namespace std;

#define pi   3.14159265358
/*
CST[0]=epsilonx;
CST[1]=epsilony;
CST[2]=epsilonz;
CST[3]=gammaxy;
CST[4]=gammaxz;
CST[5]=gammayz;
*/
//For Mnsi
// static double A  =1.270*1.e-23,
//               D  =1.140*1.e-14,
//               C11=2.830*1.e11,
//               C12=0.641*1.e11,
//               C44=1.179*1.e11,
//               SigmaS=1.e9;
//For Co
static double A = 3.45 * 1.e-23,
              D = 2.00 * 1.e-15,
              C11 = 2.42 * 1.e11,
              C12 = 1.60 * 1.e11,
              C44 = 1.28 * 1.e11,
              SigmaS=1.e9;

int StressToStrain(double *Stress, double *Strain)
{
    Strain[0]= Stress[0]*(C11+C12)/(C11-C12)/(C11+2*C12)+
              (Stress[1]+Stress[2])*(-C12/(C11-C12)/(C11+2*C12));
    Strain[1]= Stress[1]*(C11+C12)/(C11-C12)/(C11+2*C12)+
              (Stress[0]+Stress[2])*(-C12/(C11-C12)/(C11+2*C12));
    Strain[2]= Stress[2]*(C11+C12)/(C11-C12)/(C11+2*C12)+
              (Stress[0]+Stress[1])*(-C12/(C11-C12)/(C11+2*C12));
    Strain[3]= Stress[3]/C44;
    Strain[4]= Stress[4]/C44;
    Strain[5]= Stress[5]/C44;

    return 1;
}

int StrainToStress(double *Stress, double *Strain)
{
    Stress[0]= Strain[0]*C11+Strain[1]*C12+Strain[2]*C12;
    Stress[1]= Strain[0]*C12+Strain[1]*C11+Strain[2]*C12;
    Stress[2]= Strain[0]*C12+Strain[1]*C12+Strain[2]*C11;
    Stress[3]= Strain[3]*C44;
    Stress[4]= Strain[4]*C44;
    Stress[5]= Strain[5]*C44;

    return 1;
}

int BeamEllipse_Torsion(int m, int n, int l, double ****CST, double *dxy, double B_M)//B_M:Torsion moment
{
    int i, j, k;
    double x, y, Stress[6], an, bn, yc, xc;
    double LD=2*A/D;

    bn=(n-3)*dxy[1]*LD/2;
    an=(m-3)*dxy[0]*LD/2;
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x=(i-xc)*dxy[0]*LD;
                y=(j-yc)*dxy[1]*LD;

                if(x==0. && y==0.)
                {
                    y=0.1*dxy[1]*LD;
                    x=0.1*dxy[0]*LD;
                }


                Stress[0]=0;

                Stress[1]=0;

                Stress[2]=0;

                Stress[3]=0;

                Stress[4]=-2*B_M*y/pi/an/bn/bn/bn;

                Stress[5]=2*B_M*x/pi/an/an/an/bn;

                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= Stress[2];
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= Stress[4];
                CST[5][i][j][k]= Stress[5];
                CST[6][i][j][k]= 1.;

            }
        }
    }
    return 1;
}

int BeamCircleHole_Torsion(int m, int n, int l, double ****CST, double *dxy, double B_M)//B_M:Torsion moment
{
    int i, j, k, s;
    double x, y, Stress[6], yc, xc, an, bn;
    double LD=2*A/D;
    // printf("Torsion_Mode_to_be_set");
    an=(n-3)*dxy[1]*LD/2;
    s=(n-3)*0.5/2;
    bn=s*dxy[1]*LD;
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x=(i-xc)*dxy[0]*LD;
                y=(j-yc)*dxy[1]*LD;

                if(x==0. && y==0.)
                {
                    y=0.1*dxy[1]*LD;
                    x=0.1*dxy[0]*LD;
                }


                Stress[0]=0;

                Stress[1]=0;

                Stress[2]=0;

                Stress[3]=0;

                Stress[4]=-2*B_M*y/pi/(an*an*an*an-bn*bn*bn*bn);

                Stress[5]=2*B_M*x/pi/(an*an*an*an-bn*bn*bn*bn);

                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= Stress[2];
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= Stress[4];
                CST[5][i][j][k]= Stress[5];
                CST[6][i][j][k]= 1.;

            }
        }
    }
    return 1;
}

int Screw_Dislocation(int m, int n, int l, double ****CST, double *dxy, double b)
{
    int i, j, k, xc, yc;
    double x1, x2;
    double LD=2*A/D;
    // v is the poisson ratio of MnSi
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x1=(i-xc)*dxy[0]*LD;
                x2=(j-yc)*dxy[1]*LD;

                if(x1==0. && x2==0.)
                {
                    x1=0.1*dxy[0]*LD;
                    x2=0.1*dxy[1]*LD;
                }

                CST[0][i][j][k]= 0;
                CST[1][i][j][k]= 0;
                CST[2][i][j][k]= 0;
                CST[3][i][j][k]= 0;
                CST[4][i][j][k]= -b/(2*pi)*x2/(x1*x1+x2*x2);
                CST[5][i][j][k]= b/(2*pi)*x1/(x1*x1+x2*x2);
                CST[6][i][j][k]= 0;
            }
        }
    }

    return 1;
}


int Edge_Dislocation(int m, int n, int l, double ****CST, double *dxy, double b)
{
    int i, j, k, xc, yc;
    double x1, x2;
    double LD=2*A/D;
    double v=C12/(C11+C12);
    // v is the poisson ratio of MnSi
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x1=(i-xc)*dxy[0]*LD;
                x2=(j-yc)*dxy[1]*LD;

                if(x1==0. && x2==0.)
                {
                    x1=0.1*dxy[0]*LD;
                    x2=0.1*dxy[1]*LD;
                }

                CST[0][i][j][k]= b/(2*pi)*(-x2)/(x1*x1+x2*x2)
                                +b/(4*pi*(1-v))*x2*(x2*x2-x1*x1)/pow((x1*x1+x2*x2),2);
                CST[1][i][j][k]= b*(2*v-1)/(8*pi*(1-v))*(2*x2/(x1*x1+x2*x2))
                                +b/(4*pi*(1-v))*(2*x2*x1*x1)/pow((x1*x1+x2*x2),2);
                CST[2][i][j][k]=0;
                CST[3][i][j][k]= b/(2*pi)*(x1/(x1*x1+x2*x2))
                                +b/(4*pi*(1-v))*x1*(x1*x1-x2*x2)/pow((x1*x1+x2*x2),2)
                                +b*(2*v-1)/(8*pi*(1-v))*(2*x1/(x1*x1+x2*x2))
                                +b/(4*pi*(1-v))*(-2*x1*x2*x2)/pow((x1*x1+x2*x2),2);
                CST[4][i][j][k]=0;
                CST[5][i][j][k]=0;
                CST[6][i][j][k]=0;
            }
        }
    }

    return 1;
}


int Uniaxial_Strain(int m, int n, int l, double ****CST, double *CST_Uniaxial)
{
    int i, j, k, r;
    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                for(r=0;r<6;r++)
                {
                    CST[r][i][j][k]=CST_Uniaxial[r];
                }
                CST[6][i][j][k]=0;
            }
        }
    }

    return 1;
}

int Uniaxial_Stress(int m, int n, int l, double ****CST, double *CST_Uniaxial)
{
    int i, j, k, r;
    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                for(r=0;r<6;r++)
                {
                    CST[r][i][j][k]=CST_Uniaxial[r];
                }
                CST[6][i][j][k]=1;
            }
        }
    }

    return 1;
}


int ZeroField(int m, int n, int l, double ****CST)
{
    int i, j, k, r;
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l;k++)
            {
                for(r=0;r<7;r++)
                {
                    CST[r][i][j][k]=0.;
                }
            }
        }
    }
    return 1;
}

int CrackI_Stress_CC(int m, int n, int l, double ****CST, double *dxy, double a)
{
    int i, j, k, xc, yc;
    double x1, x2, theta, phi, phi1, chi, eta, xi, r, Stress[6], KI;
    double LD=2*A/D;

    chi=(C11*C11-C12*C12-2*C12*C44)/(2*C11*C44);
    //chi<1
    xi =sqrt((1-chi)/2);    eta=sqrt((1+chi)/2);


    xc=(m-1)/2; yc=(n-1)/2;
    KI=sqrt(2*pi*a)*SigmaS*(C11+C12)/(C11);

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x1   =(i-xc)*dxy[0]*LD;
                x2   =(j-yc)*dxy[1]*LD;
                r    =sqrt(x1*x1+x2*x2);
                theta=atan2(x2,x1);
                //theta=atan2(y,x)
                //x1=r*cos(theta), x2=r*sin(theta)
                //phi  =atan2(eta*sin(theta),cos(theta)+xi*sin(theta));
                phi  =atan2(eta*x2,x1+xi*x2);
                phi1 =atan2(eta*x2,x1-xi*x2);

                if(r==0.)
                {
                    r=0.05*sqrt(dxy[0]*dxy[0]+dxy[1]*dxy[1]);
                }


                Stress[0]=KI/sqrt(2*pi*r)*(pow(pow(cos(theta)+xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                          (cos(phi/2)+eta/xi*sin(phi/2))+
                                           pow(pow(cos(theta)-xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                          (cos(phi1/2)-eta/xi*sin(phi1/2)));

                Stress[1]=KI/sqrt(2*pi*r)*(pow(pow(cos(theta)+xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                          (cos(phi/2)-eta/xi*sin(phi/2))+
                                           pow(pow(cos(theta)-xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                          (cos(phi1/2)+eta/xi*sin(phi1/2)));

                Stress[3]=-KI/(sqrt(2*pi*r)*xi)*(pow(pow(cos(theta)+xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                                 cos(phi/2)-
                                                 pow(pow(cos(theta)-xi*sin(theta),2)+pow(eta*sin(theta),2),-0.25)*
                                                 cos(phi1/2));


                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= 0.;
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= 0.;
                CST[5][i][j][k]= 0.;
                CST[6][i][j][k]= 1.;

            }
        }
    }

    return 1;
}

int Concentrated_Force(int m, int n, int l, double ****CST, double *dxy, double FC)
{
    int i, j, k, xc, yc;
    double x, y, z, r, Stress[6], Amplitude;
    double LD=2*A/D;
    double v=C12/(C11+C12);

    xc=(n-1)/2; yc=(m-1)/2;

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<l-1;k++)
            {
                x   =(i-xc)*dxy[0]*LD;
                y   =(j-yc)*dxy[1]*LD;
                z   =(k-l+2)*dxy[2]*LD;

                if(x==0. && y==0. && z==0.)
                {
                    x=0.1*dxy[0]*LD;
                    y=0.1*dxy[1]*LD;
                    z=0.1*dxy[2]*LD;
                }

                r        =sqrt(x*x+y*y+z*z);
                Amplitude=-FC/(2*pi*r*r);
                z        =-z;


                Stress[0]=Amplitude*(3*x*x*z/(r*r*r)-(1-2*v)*(z/r-r/(r+z)+x*x*(2*r+z)/(r*(r+z)*(r+z))));

                Stress[1]=Amplitude*(3*y*y*z/(r*r*r)-(1-2*v)*(z/r-r/(r+z)+y*y*(2*r+z)/(r*(r+z)*(r+z))));

                Stress[2]=Amplitude*(3*z*z*z/(r*r*r));

                Stress[3]=Amplitude*(3*x*y*z/(r*r*r)-(1-2*v)*x*y*(2*r+z)/(r*(r+z)*(r+z)));

                Stress[4]=Amplitude*(3*x*z*z/(r*r*r));

                Stress[5]=Amplitude*(3*y*z*z/(r*r*r));

                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= Stress[2];
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= Stress[4];
                CST[5][i][j][k]= Stress[5];
                CST[6][i][j][k]= 1.;

            }

            CST[0][i][j][l-1]= 0;
            CST[1][i][j][l-1]= 0;
            CST[2][i][j][l-1]= 0;
            CST[3][i][j][l-1]= 0;
            CST[4][i][j][l-1]= 0;
            CST[5][i][j][l-1]= 0;
            CST[6][i][j][l-1]= 1.;
        }
    }

    return 1;
}

int Beam_Concentrated_Force(int m, int n, int l, double ****CST, double *dxy, double FC_Beam)
{
    int i, j, k;
    double x, y, Stress[6], h, yc, xc;
    double LD=2*A/D;

    h=(n-1)*dxy[1]*LD;
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x=(i-xc)*dxy[0]*LD;
                y=(j-yc)*dxy[1]*LD;

                if(x==0. && y==0.)
                {
                    x=0.1*dxy[0]*LD;
                    y=0.1*dxy[1]*LD;
                }


                x=-x+(m-1)*dxy[0]*LD/2;

                Stress[0]=-12*FC_Beam*x*y/(h*h*h);

                Stress[1]=0;

                Stress[2]=0;

                Stress[3]=-3*FC_Beam/(2*h)+6*FC_Beam*y*y/(h*h*h);

                Stress[4]=0;

                Stress[5]=0;

                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= Stress[2];
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= Stress[4];
                CST[5][i][j][k]= Stress[5];
                CST[6][i][j][k]= 1.;

            }
        }
    }
    return 1;
}



int Beam_Uniform_Force(int m, int n, int l, double ****CST, double *dxy, double FU_Beam)
{
    int i, j, k;
    double x, y, Stress[6], h, xc, yc, J;
    double LD=2*A/D, Q;
    //double v=C12/(C11+C12);
    Q = FU_Beam;
    h=(n-1)*dxy[1]*LD;
    xc=(m-1)/2; yc=(n-1)/2;
    J=h*h*h/12;
    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x=(i-xc)*dxy[0]*LD;
                y=(j-yc)*dxy[1]*LD;

                if(x==0. && y==0.)
                {
                    x=0.1*dxy[0]*LD;
                    y=0.1*dxy[1]*LD;
                }


                x=-x+(m-1)*dxy[0]*LD/2;

                Stress[0]=-Q*x*x*y/(2*J)+(Q/(2*J))*(2*y*y*y/3-h*h*y/10);

                Stress[1]=-Q/2*(1-3*y/h+4*y*y*y/(h*h*h));

                Stress[2]=0;

                Stress[3]=(Q/(2*J))*(y*y-h*h/4)*x;

                Stress[4]=0;

                Stress[5]=0;

                CST[0][i][j][k]= Stress[0];
                CST[1][i][j][k]= Stress[1];
                CST[2][i][j][k]= Stress[2];
                CST[3][i][j][k]= Stress[3];
                CST[4][i][j][k]= Stress[4];
                CST[5][i][j][k]= Stress[5];
                CST[6][i][j][k]= 1.;

            }
        }
    }
    return 1;
}

int Beam_Bending(int m, int n, int l, double ****CST, double *dxy, double B_M)//B_M: Bending moment
{
    int i, j, k;
    double x, y, Stress[6], Strain[6], h, yc, xc;
    double LD=2*A/D;

    h=(n-1)*dxy[1]*LD;
    xc=(m-1)/2; yc=(n-1)/2;

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                x=(i-xc)*dxy[0]*LD;
                y=(j-yc)*dxy[1]*LD;

                if(x==0. && y==0.)
                {
                    x=0.1*dxy[0]*LD;
                    y=0.1*dxy[1]*LD;
                }


                Stress[0]=12*B_M*y/(h*h*h);

                Stress[1]=0;

                Stress[2]=0;

                Stress[3]=0;

                Stress[4]=0;

                Stress[5]=0;

                // CST[0][i][j][k]= Stress[0];
                // CST[1][i][j][k]= Stress[1];
                // CST[2][i][j][k]= Stress[2];
                // CST[3][i][j][k]= Stress[3];
                // CST[4][i][j][k]= Stress[4];
                // CST[5][i][j][k]= Stress[5];
                // CST[6][i][j][k]= 1.;


                Strain[0]= Stress[0]*(C11+C12)/(C11-C12)/(C11+2*C12)+
                        (Stress[1]+Stress[2])*(-C12/(C11-C12)/(C11+2*C12));
                Strain[1]= Stress[1]*(C11+C12)/(C11-C12)/(C11+2*C12)+
                        (Stress[0]+Stress[2])*(-C12/(C11-C12)/(C11+2*C12));
                Strain[2]= Stress[2]*(C11+C12)/(C11-C12)/(C11+2*C12)+
                        (Stress[0]+Stress[1])*(-C12/(C11-C12)/(C11+2*C12));
                Strain[3]= Stress[3]/C44;
                Strain[4]= Stress[4]/C44;
                Strain[5]= Stress[5]/C44;


                CST[0][i][j][k]= Strain[0];
                CST[1][i][j][k]= Strain[1];
                CST[2][i][j][k]= Strain[2];
                CST[3][i][j][k]= Strain[3];
                CST[4][i][j][k]= Strain[4];
                CST[5][i][j][k]= Strain[5];
                CST[6][i][j][k]= 0;

            }
        }
    }
    return 1;
}

int Elastic_Field_Initial(double ****CST, Input_Parameter_NMD *IPN)
{
    int m, n, l;
    m=2*IPN->xn+3;  n=2*IPN->yn+3;
    Elastic_Field_Mode EFM;
    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l=2*IPN->zn+3;
    }

    EFM = (Elastic_Field_Mode)IPN->ElasticFieldMode;
    switch (EFM)
    {
        case EFM_Zerofield:
        {
            ZeroField(m,n,l,CST);
            break;
        }

        case EFM_UniaxialStrain:
        {
            Uniaxial_Strain(m,n,l,CST,IPN->CST_Uniaxial);
            break;
        }

        case EFM_UniaxialStress:
        {
            Uniaxial_Stress(m,n,l,CST,IPN->CST_Uniaxial);
            break;
        }

        case EFM_EdgeDislocation:
        {
            Edge_Dislocation(m,n,l,CST,IPN->dxy,IPN->b);
            break;
        }

        case EFM_ScrewDislocation:
        {
            Screw_Dislocation(m,n,l,CST,IPN->dxy,IPN->b);
            break;
        }

        case EFM_CrackI:
        {
            CrackI_Stress_CC(m,n,l,CST,IPN->dxy,IPN->a);
            break;
        }

        case EFM_ConcentratedForce:
        {
            Concentrated_Force(m,n,l,CST,IPN->dxy,IPN->FC);
            break;
        }

        case EFM_BeamUniformForce:
        {
            Beam_Uniform_Force(m,n,l,CST,IPN->dxy,IPN->FU_Beam);
            break;
        }

        case EFM_BeamConcentratedForce:
        {
            Beam_Concentrated_Force(m,n,l,CST,IPN->dxy,IPN->FC_Beam);
            break;
        }

        case EFM_BeamBending:
        {
            Beam_Bending(m,n,l,CST,IPN->dxy,IPN->B_M);
            break;
        }

        case EFM_ReadAnsys:
        {
            Read_Elastic_Field_ANSYS(IPN,CST);
            break;
        }

        case EFM_HollowTorsion:
        {
            BeamCircleHole_Torsion(m,n,l,CST,IPN->dxy,IPN->B_M);
            break;
        }
        
        case EFM_EllipseTorsion:
        {
            BeamEllipse_Torsion(m,n,l,CST,IPN->dxy,IPN->B_M);
            break;
        }
    
        default:
        {
            ZeroField(m,n,l,CST);
        }
        break;
    }

    return 1;
}

/*
int StressToPS(double *Stress, double *PrincipalStress)
{
    double a,b,c,d,u,x1,x2,x3,q;
    complex<double> v,M,N,vvv;
    a=-1.;
    b=Stress[0]+Stress[1]+Stress[2];
    c=-(Stress[0]*Stress[2]+Stress[1]*Stress[2]+Stress[0]*Stress[1]+Stress[4]*Stress[4]+Stress[3]*Stress[3]+Stress[5]*Stress[5]);
    d=Stress[0]*Stress[1]*Stress[2]+2*Stress[3]*Stress[4]*Stress[5]-Stress[4]*Stress[4]*Stress[1]-Stress[5]*Stress[5]*Stress[0]-Stress[3]*Stress[3]*Stress[2];
    u=(9*a*b*c-27*a*a*d-2*b*b*b)/(54*a*a*a);
    q=(3*(4*a*c*c*c-b*b*c*c-18*a*b*c*d+27*a*a*d*d+4*b*b*b*d));
    if (q>=0)
    {
        q=pow(q,0.5);
        complex<double> vv(q,0.0);
        v=vv/(18*a*a);
    }
    else
    {
        q=pow(-q,0.5);
        complex<double> vv(0.0,q);
        v=vv/(18*a*a);
    }

    //v=std::pow(q,0.5);
    if (std::abs(u+v) >= std::abs(u-v))
    {
        //vvv=u+v;
        M=std::pow(u+v,0.333333333333333);
        //std::cout<<std::pow(u+v,0.3333333333333);
    }
    else
    {
        M=std::pow(u-v,0.333333333333333);
    }

    if (real(M) != 0 ||imag(M)  != 0)
    {
        N=(b*b-3*a*c)/(9*a*M);
    }
    else
    {
        N=0;
    }
    complex<double> omiga1(-0.5,sqrt(3)/2);
    complex<double> omiga2(-0.5,-sqrt(3)/2);

    std::cout<<M+N-b/(3*a);
    x1=real(M+N-b/(3*a));
    x2=real(omiga1*M+omiga2*N-b/(3*a));
    x3=real(omiga2*M+omiga1*N-b/(3*a));
    if (x1>=x2)
    {
        PrincipalStress[0]=x1;
        PrincipalStress[1]=x2;
    }
    else
    {
        PrincipalStress[0]=x2;
        PrincipalStress[1]=x1;
    }
    if  (PrincipalStress[0]<x3)
    {
        PrincipalStress[2]=PrincipalStress[1];
        PrincipalStress[1]=PrincipalStress[0];
        PrincipalStress[0]=x3;
    }
    else
    {
        if (PrincipalStress[1]>=x3)
        {
            PrincipalStress[2]=x3;
        }
        else
        {
            PrincipalStress[2]=PrincipalStress[1];
            PrincipalStress[1]=x3;
        }
    }
    return 1;
}
*/
int Read_Elastic_Field_ANSYS(Input_Parameter_NMD *IPN, double ****CST)
{
    //printf("111");
    int m, n, i, j, k;
    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    FILE *fp_S11, *fp_S22, *fp_S33, *fp_S12, *fp_S13, *fp_S23;
    char *filename, *filename_S11, *filename_S22, *filename_S33, *filename_S12, *filename_S13, *filename_S23;
    filename_S11 = Make1DString(500); filename_S22 = Make1DString(500);
    filename_S33 = Make1DString(500); filename_S12 = Make1DString(500);
    filename_S13 = Make1DString(500); filename_S23 = Make1DString(500);
    //filename     = Make1DString(500);

    //sprintf(filename,"Elastic Field of ANSYS");
    //mkdir(filename);
    sprintf(filename_S11,"Elastic Field of ANSYS\\S11.dat");
    sprintf(filename_S22,"Elastic Field of ANSYS\\S22.dat");
    sprintf(filename_S33,"Elastic Field of ANSYS\\S33.dat");
    sprintf(filename_S12,"Elastic Field of ANSYS\\S12.dat");
    sprintf(filename_S13,"Elastic Field of ANSYS\\S13.dat");
    sprintf(filename_S23,"Elastic Field of ANSYS\\S23.dat");
    fp_S11=fopen(filename_S11,"r+"); fp_S22=fopen(filename_S22,"r+");
    fp_S33=fopen(filename_S33,"r+"); fp_S12=fopen(filename_S12,"r+");
    fp_S13=fopen(filename_S13,"r+"); fp_S23=fopen(filename_S23,"r+");
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {

            fscanf(fp_S11,"%lf",&CST[0][i][j][0]);
            fscanf(fp_S22,"%lf",&CST[1][i][j][0]);
            fscanf(fp_S33,"%lf",&CST[2][i][j][0]);
            fscanf(fp_S12,"%lf",&CST[3][i][j][0]);
            fscanf(fp_S13,"%lf",&CST[4][i][j][0]);
            fscanf(fp_S23,"%lf",&CST[5][i][j][0]);
/*
            fscanf(fp_S11,"%lf",&CST[0][j][i][0]);
            fscanf(fp_S22,"%lf",&CST[1][j][i][0]);
            fscanf(fp_S33,"%lf",&CST[2][j][i][0]);
            fscanf(fp_S12,"%lf",&CST[3][j][i][0]);
            fscanf(fp_S13,"%lf",&CST[4][j][i][0]);
            fscanf(fp_S23,"%lf",&CST[5][j][i][0]);
*/
        }
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            for(k=0;k<6;++k)
            {
                CST[k][i][j][0]=IPN->ElasticFieldScaleFactor*CST[k][i][j][0];
            }
            CST[6][i][j][0]=1;
        }
    }

    fclose(fp_S11);
    fclose(fp_S22);
    fclose(fp_S33);
    fclose(fp_S12);
    fclose(fp_S13);
    fclose(fp_S23);

    //free(filename);
    free(filename_S11);
    free(filename_S22);
    free(filename_S33);
    free(filename_S12);
    free(filename_S13);
    free(filename_S23);

    return 1;
}

int StressToPS(double *Stress, double *PrincipalStress)
{
    double a,b,c,d,u,v,x1,x2,x3,M,N;
    a=-1.;
    b=Stress[0]+Stress[1]+Stress[2];
    c=-(Stress[0]*Stress[2]+Stress[1]*Stress[2]+Stress[0]*Stress[1]-Stress[4]*Stress[4]-Stress[3]*Stress[3]-Stress[5]*Stress[5]);
    d=Stress[0]*Stress[1]*Stress[2]+2*Stress[3]*Stress[4]*Stress[5]-Stress[4]*Stress[4]*Stress[1]-Stress[5]*Stress[5]*Stress[0]-Stress[3]*Stress[3]*Stress[2];
    u=-(9*a*b*c-27*a*a*d-2*b*b*b)/(54*a*a*a);
    v=(3*a*c-b*b)/(9*a*a)*(3*a*c-b*b)/(9*a*a)*(3*a*c-b*b)/(9*a*a);
    if (d!=0)
    {
        x1=-b/3/a+pow((-u+pow(u*u+v,0.5)),1/3)+pow((-u-pow(u*u+v,0.5)),1/3);
        M=-b/a-x1;
        N=-d/a/x1;
        x2=M/2+sqrt(M*M-4*N)/2;
        x3=M/2-sqrt(M*M-4*N)/2;
        if (x1>=x2)
        {
            PrincipalStress[0]=x1;
            PrincipalStress[1]=x2;
        }
        else
        {
            PrincipalStress[0]=x2;
            PrincipalStress[1]=x1;
        }
        if  (PrincipalStress[0]<x3)
        {
            PrincipalStress[2]=PrincipalStress[1];
            PrincipalStress[1]=PrincipalStress[0];
            PrincipalStress[0]=x3;
        }
        else
        {
            if (PrincipalStress[1]>=x3)
            {
                PrincipalStress[2]=x3;
            }
            else
            {
                PrincipalStress[2]=PrincipalStress[1];
                PrincipalStress[1]=x3;
            }
        }
    }
    else
    {
        PrincipalStress[0]=(-b-sqrt(b*b-4*a*c))/(2*a);
        PrincipalStress[1]=PrincipalStress[0];
        PrincipalStress[2]=(-b+sqrt(b*b-4*a*c))/(2*a);
    }
    //std::cout<<a*PrincipalStress[0]*PrincipalStress[0]+b*PrincipalStress[0]+c;
    return 1;
}
