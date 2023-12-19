#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Array.h>
#include <Energy.h>
#include <Variable.h>
#include <BoundaryShape.h>
#include <Initial_State.h>

#define pi   3.14159265358

using namespace std;
//For MnSi
/*
static double L1         =-0.70*1.e6,
              L2         = 0.60*1.e6,
              L3         = 1.65*1.e6,
              //G          = 2.56*1.e-6,
              A          = 1.27*1.e-23,
              Alpha      = 6.44*1.e-7,
              Beta       = 3.53*1.e-16,
              Keff       = 0,
              Q          =-2.00*1.e7,
              Ms         = 1.63*1.e5,
              D          = 1.14*1.e-14,
              Lo1        =-0.57*1.e-4,
              Lo2        = 1.15*1.e-4,
              Lo3        =-0.57*1.e-4,
              C11        = 2.830*1.e11,
              C12        = 0.641*1.e11,
              C44        = 1.179*1.e11,
              mu0        = 1.25664*1.e-6,
              gamma      = 1.76*1.e11;
*/
/*
//For Co
static double L1 = -3.1 * 1.e6,
              L2 = 0.,
              L3 = 2.453 * 1.e6,
              //G          = 2.89*1.e-8,
              A = 3.45 * 1.e-23,
              Alpha = 6.44 * 1.e-7,
              Beta = 5.255 * 1.e-16,
              Keff = -2.02 * 1.e-7,
              //Keff = 0,
              Q = 0.,
              Ms = 7.00 * 1.e5,
              D = 2.00 * 1.e-15,
              Lo1 = 0.,
              Lo2 = 0.,
              Lo3 = 0.,
              C11 = 2.42 * 1.e11,
              C12 = 1.60 * 1.e11,
              C44 = 1.28 * 1.e11;
*/
/*
//For Co_1
static double L1 = 0,
              L2 = 0.,
              L3 = 0,
              //G          = 2.89*1.e-8,
              A = 2.06 * 1.e-23,
              Alpha = 6.44 * 1.e-7,
              Beta = 5.255 * 1.e-16,
              Keff = -1.74 * 1.e-7,
              //Keff = 0,
              Q = 0.,
              Ms = 6.97 * 1.e5,
              D = 2.68 * 1.e-15,
              Lo1 = 0.,
              Lo2 = 0.,
              Lo3 = 0.,
              C11 = 2.42 * 1.e11,
              C12 = 1.60 * 1.e11,
              C44 = 1.28 * 1.e11;
*/
//Ld=34.5nm
//Ld=15.38nm

//For Co_2

static double L1 = 0,
              L2 = 0.,
              L3 = 0,
              //G          = 2.89*1.e-8,
              A = 3.45 * 1.e-23,
              Alpha = 6.44 * 1.e-7,
              Beta = 5.255 * 1.e-16,
              Keff = -2.02 * 1.e-7,
              //Keff = 0,
              Q = 0.,
              Ms = 6.97 * 1.e5,
              D = 2.00 * 1.e-15,
              Lo1 = 0.,
              Lo2 = 0.,
              Lo3 = 0.,
              C11 = 2.42 * 1.e11,
              C12 = 1.60 * 1.e11,
              C44 = 1.28 * 1.e11,
              mu0        = 1.25664*1.e-6,
              gamma      = 1.76*1.e11;

double        G          = D*D/(4*A),
              //k_eff      = Keff/(G),
              k_eff      = -0.3,
              l1_Strain  = L1/(G*Ms*Ms),
              l2_Strain  = L2/(G*Ms*Ms),
              l3_Strain  = L3/(G*Ms*Ms),
              q_Strain   = Q/(G*Ms*Ms),
              lo1_Strain = 2*Lo1/(D*Ms*Ms),
              lo2_Strain = 2*Lo2/(D*Ms*Ms),
              lo3_Strain = 2*Lo3/(D*Ms*Ms),

              l1_Stress  = l1_Strain/(C11-C12),
              l2_Stress  = l2_Strain/(C11-C12),
              l3_Stress  = l3_Strain/C44,
              q_Stress   = (q_Strain*(C11-C12)-C12*(l1_Strain+l2_Strain))/(C11-C12)/(C11+2*C12),
              lo1_Stress = (C11*lo1_Strain+C12*(lo1_Strain-lo2_Strain-lo3_Strain))/(C11-C12)/(C11+2*C12),
              lo2_Stress = (C11*lo2_Strain-C12*(lo1_Strain-lo2_Strain+lo3_Strain))/(C11-C12)/(C11+2*C12),
              lo3_Stress = (C11*lo3_Strain-C12*(lo1_Strain+lo2_Strain-lo3_Strain))/(C11-C12)/(C11+2*C12),
              Ms_LLG     = Ms,
              mu0_LLG    = mu0,
              gamma_LLG  = gamma,

              d_epsilon  = -39/2,
              //d_epsilon  = 0,
              k_0        = -50/2*(-0.3)/k_eff,
              //k_0        = -12.5,
              //k_0        = 0,
              A_e        = 0;
            //H_D      = 0.0589;
            //d_epsilon  = 0.;
//enum Energy_Terms{Exchange, DM_Bloch, DM_Neel, Zeeman, Landau, Magnetoelastic, Axial_Anisotropy};
int N_Energy_Terms = 7;

int Effective_Field_LocalCell_M0(Energy_Variables *EV, double *H_eff_Local)
{
    double ***M_Local, *dxy;
    double H_Ex = 0, H_DM_Bloch = 0, H_DM_Neel = 0, H_ZM = 0, H_Ld = 0, H_ME = 0, H_AA = 0;
    double *Energy_Coefficient, C_Ex, C_DM_Bloch, C_DM_Neel, C_ZM, C_LD, C_ME, C_AA;
    
    double h, t, **CST, **H_eff;
    int i, r, s, *LocalEnergy, p, m, n, l;
    double l1, l2, l3, q, lo1, lo2, lo3;

    m   = EV->m;     n   = EV->n;
    l   = EV->l;     dxy = EV->dxy;
    t   = EV->t;     h   = EV->h;
    CST = EV->CST;
    M_Local     = EV->M_Local;
    LocalEnergy = EV->LocalEnergy;
    Energy_Coefficient = EV->Energy_Coefficient;


    double dx=dxy[0], dy=dxy[1], dz=dxy[2];
    double S11, S12, S13, S22, S23, S33, CST_Type;
    double MxO, MxN, MxS, MxW, MxE, MxL, MxU, MxNN, MxSS, MxWW, MxEE, MxLL, MxUU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU, MyNN, MySS, MyWW, MyEE, MyLL, MyUU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU, MzNN, MzSS, MzWW, MzEE, MzLL, MzUU;

    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    enum Local_Point_Postion_L{Null, O, N, S, W, E, L, U, NN, SS, WW, EE, LL, UU};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    //NN(i-2,j,k),SS(i+2,j,k), WW(i,j-2,k), EE(i,j+2,k), LL(i,j,k-2), UU(i,j,k+2);
    Local_Point_Postion LPP;
    Local_Point_Postion_L LPPL;
    Energy_Terms ET;

    for(i=0;i<14;++i)
    {
        LPPL=(Local_Point_Postion_L) i;
        switch(LPPL)
        {
            case O :
            {
                MxO  = EV->M0[0];    MyO  = EV->M0[1];    MzO  = EV->M0[2];
                break;
            }

            case N :
            {
                MxN  = *M_Local[1][i];    MyN  = *M_Local[2][i];    MzN  = *M_Local[3][i];
                break;
            }

            case S :
            {
                MxS  = *M_Local[1][i];    MyS  = *M_Local[2][i];    MzS  = *M_Local[3][i];
                break;
            }

            case W :
            {
                MxW  = *M_Local[1][i];    MyW  = *M_Local[2][i];    MzW  = *M_Local[3][i];
                break;
            }

            case E :
            {
                MxE  = *M_Local[1][i];    MyE  = *M_Local[2][i];    MzE  = *M_Local[3][i];
                break;
            }

            case L :
            {
                MxL  = *M_Local[1][i];    MyL  = *M_Local[2][i];    MzL  = *M_Local[3][i];
                break;
            }


            case U :
            {
                MxU  = *M_Local[1][i];    MyU  = *M_Local[2][i];    MzU  = *M_Local[3][i];
                break;
            }

            case NN :
            {
                MxNN = *M_Local[1][i];    MyNN = *M_Local[2][i];    MzNN = *M_Local[3][i];
                break;
            }

            case SS :
            {
                MxSS = *M_Local[1][i];    MySS = *M_Local[2][i];    MzSS = *M_Local[3][i];
                break;
            }

            case WW :
            {
                MxWW = *M_Local[1][i];    MyWW = *M_Local[2][i];    MzWW = *M_Local[3][i];
                break;
            }

            case EE :
            {
                MxEE = *M_Local[1][i];    MyEE = *M_Local[2][i];    MzEE = *M_Local[3][i];
                break;
            }

            case LL :
            {
                MxLL = *M_Local[1][i];    MyLL = *M_Local[2][i];    MzLL = *M_Local[3][i];
                break;
            }

            case UU :
            {
                MxUU = *M_Local[1][i];    MyUU = *M_Local[2][i];    MzUU = *M_Local[3][i];
                break;
            }

            default:
                break;
        }
    }

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch(ET)
        {
            case Exchange:
            {
                C_Ex = Energy_Coefficient[i];
                break;
            }

            case DM_Bloch:
            {
                C_DM_Bloch = Energy_Coefficient[i];
                break;
            }

            case DM_Neel:
            {
                C_DM_Neel = Energy_Coefficient[i];
                break;
            }

            case Zeeman:
            {
                C_ZM = Energy_Coefficient[i];
                break;
            }

            case Landau:
            {
                C_LD = Energy_Coefficient[i];
                break;
            }

            case Magnetoelastic:
            {
                C_ME = Energy_Coefficient[i];
                break;
            }

            case Axial_Anisotropy:
            {
                C_AA = Energy_Coefficient[i];
                break;
            }

            default:
            {
                break;
            }
        }
    }

    H_eff = Make2DArray(3, 7);

    ///////////////////////////////////////////////////////////
    if(l==1)//2D Case
    {
        for(p=0;p<7;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
///////////////////////////////////////////////////////////
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        /////////////////////////////// x direction ////////////////////////////
                        H_Ex        = 0.;
                        H_DM_Bloch  = 2 * (MzE - MzW) / (2 * dy);
                        H_DM_Neel   = 2 * (1 + (S11 + S22) * d_epsilon) * ((MzS - MzN) / (2 * dx));
                        H_ZM        = 0.;
                        H_Ld        = 2 * t * MxO + 4 * MxO * (pow(MxO, 2) + pow(MyO, 2) + pow(MzO, 2));
                        H_ME        = l3 * (MyO * S12 + MzO * S13) - (lo3 * S11 + lo2 * S22 + lo1 * S33) * (MzE - MzW) / (2 * dy) +
                                      2 * MxO * (l1 * S11 + l2 * S22 + q * (S11 + S22 + S33));
                        H_AA        = 0.;

                        H_eff[0][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// y direction ////////////////////////////
                        H_Ex        = 0.;
                        H_DM_Bloch  = -2 * (MzS - MzN) / (2 * dx);
                        H_DM_Neel   = 2 * (1 + (S11 + S22) * d_epsilon) * ((MzE - MzW) / (2 * dy));
                        H_ZM        = 0.;
                        H_Ld        = 2 * t * MyO + 4 * MyO * (pow(MxO, 2) + pow(MyO, 2) + pow(MzO, 2));
                        H_ME        = l3 * (MxO * S12 + MzO * S23) + (lo2 * S11 + lo3 * S22 + lo1 * S33) * (MzS - MzN) / (2 * dx) +
                                      2 * MyO * (l1 * S22 + l2 * S33 + q * (S11 + S22 + S33));
                        H_AA        = 0.;

                        H_eff[1][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// z direction ////////////////////////////
                        H_Ex        = 0.;
                        H_DM_Bloch  = 2 * (MyS - MyN) / (2 * dx) - 2 * (MxE - MxW) / (2 * dy);
                        H_DM_Neel   = -2 * (1 + (S11 + S22) * d_epsilon) * ((MxS - MxN) / (2 * dx) + (MyE - MyW) / (2 * dy));
                        H_ZM        = -2 * h;
                        H_Ld        = 2 * t * MzO + 4 * MzO * (pow(MxO, 2) + pow(MyO, 2) + pow(MzO, 2));
                        H_ME        = l3 * (MxO * S13 + MyO * S23) - (lo2 * S11 + lo1 * S22 + lo3 * S33) * (MyS - MyN) / (2 * dx) +
                                      (lo1 * S11 + lo2 * S22 + lo3 * S33) * (MxE - MxW) / (2 * dy) +
                                      2 * MzO * (l2 * S11 + l1 * S33 + q * (S11 + S22 + S33));
                        H_AA        = 2 * k_eff * (1 + (S11 + S22) * k_0) * MzO;

                        H_eff[2][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
///////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        /////////////////////////////// x direction ////////////////////////////
                        H_Ex        = (MxO - MxNN) / (2 * dx * dx);
                        H_DM_Bloch  = 0.;
                        H_DM_Neel   = -2 * MzN * (1 + (S11 + S22) * d_epsilon) / (2 * dx);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = 0.;
                        H_AA        = 0.;

                        H_eff[0][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// y direction ////////////////////////////
                        H_Ex        = (MyO - MyNN) / (2 * dx * dx);
                        H_DM_Bloch  = 2 * MzN / (2 * dx);
                        H_DM_Neel   = 0.;
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = -MzN * (lo2 * S11 + lo1 * S22 + lo3 * S33) / (2 * dx);
                        H_AA        = 0.;

                        H_eff[1][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// z direction ////////////////////////////
                        H_Ex        = (MzO - MzNN) / (2 * dx * dx);
                        H_DM_Bloch  = -2 * MyN / (2 * dx);
                        H_DM_Neel   = 2 * MxN * (1 + (S11 + S22) * d_epsilon) / (2 * dx);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = MyN * (lo2 * S11 + lo3 * S22 + lo1 * S33) / (2 * dx);
                        H_AA        = 0.;

                        H_eff[2][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        /////////////////////////////// x direction ////////////////////////////
                        H_Ex        = -(MxSS - MxO) / (2 * dx * dx);
                        H_DM_Bloch  = 0.;
                        H_DM_Neel   = 2 * MzS * (1 + (S11 + S22) * d_epsilon) / (2 * dx);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = 0.;
                        H_AA        = 0.;

                        H_eff[0][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// y direction ////////////////////////////
                        H_Ex        = -(MySS - MyO) / (2 * dx * dx);
                        H_DM_Bloch  = -2 * MzS / (2 * dx);
                        H_DM_Neel   = 0.;
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = MzS * (lo2 * S11 + lo1 * S22 + lo3 * S33) / (2 * dx);
                        H_AA        = 0.;

                        H_eff[1][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// z direction ////////////////////////////
                        H_Ex        = -(MzSS - MzO) / (2 * dx * dx);
                        H_DM_Bloch  = 2 * MyS / (2 * dx);
                        H_DM_Neel   = -2 * MxS * (1 + (S11 + S22) * d_epsilon) / (2 * dx);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = -MyS * (lo2 * S11 + lo3 * S22 + lo1 * S33) / (2 * dx);
                        H_AA        = 0.;

                        H_eff[2][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        /////////////////////////////// x direction ////////////////////////////
                        H_Ex        = (MxO - MxWW) / (2 * dy * dy);
                        H_DM_Bloch  = -2 * MzW / (2 * dy);
                        H_DM_Neel   = 0.;
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = MzW * (lo1 * S11 + lo2 * S22 + lo3 * S33) / (2 * dy);
                        H_AA        = 0.;

                        H_eff[0][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// y direction ////////////////////////////
                        H_Ex        = (MyO - MyWW) / (2 * dy * dy);
                        H_DM_Bloch  = 0.;
                        H_DM_Neel   = -2 * MzW * (1 + (S11 + S22) * d_epsilon);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = 0.;
                        H_AA        = 0.;

                        H_eff[1][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// z direction ////////////////////////////
                        H_Ex        = (MzO - MzWW) / (2 * dy * dy);
                        H_DM_Bloch  = 2 * MxW / (2 * dy);
                        H_DM_Neel   = 2 * MyW * (1 + (S11 + S22) * d_epsilon);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = -MxW * (lo3 * S11 + lo2 * S22 + lo1 * S33) / (2 * dy);
                        H_AA        = 0.;

                        H_eff[2][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        /////////////////////////////// x direction ////////////////////////////
                        H_Ex        = -(MxEE - MxO) / (2 * dy * dy);
                        H_DM_Bloch  = 2 * MzE / (2 * dy);
                        H_DM_Neel   = 0.;
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = -MzE * (lo1 * S11 + lo2 * S22 + lo3 * S33) / (2 * dy);
                        H_AA        = 0.;

                        H_eff[0][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// y direction ////////////////////////////
                        H_Ex        = -(MyEE - MyO) / (2 * dy * dy);
                        H_DM_Bloch  = 0.;
                        H_DM_Neel   = 2 * MzE * (1 + (S11 + S22) * d_epsilon);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = 0.;
                        H_AA        = 0.;

                        H_eff[1][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                        /////////////////////////////// z direction ////////////////////////////
                        H_Ex        = -(MzEE - MzO) / (2 * dy * dy);
                        H_DM_Bloch  = -2 * MxE / (2 * dy);
                        H_DM_Neel   = -2 * MyE * (1 + (S11 + S22) * d_epsilon);
                        H_ZM        = 0.;
                        H_Ld        = 0.;
                        H_ME        = MxE * (lo3 * S11 + lo2 * S22 + lo1 * S33) / (2 * dy);
                        H_AA        = 0.;

                        H_eff[2][p] = C_Ex * H_Ex + C_DM_Bloch * H_DM_Bloch + C_DM_Neel * H_DM_Neel + C_ZM * H_ZM + C_LD * H_Ld + C_ME * H_ME + C_AA * H_AA;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
///////////////////////////////////////////////////////////
                default:
                {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    break;
                }
            }
        }
    }else//3D Case////////////////////////////////////////
    {
        for(p=0;p<7;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
///////////////////////////////////////////////////////////
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
///////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }

                    break;
                }
///////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
///////////////////////////////////////////////////////////
                case Lower:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;

                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
///////////////////////////////////////////////////////////
                case Upper:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];

                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];

                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }

                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }else
                    {
                        H_eff[0][p] = 0;

                        H_eff[1][p] = 0;

                        H_eff[2][p] = 0;
                    }
                    break;
                }
            }
        }
    }

    for(s=0;s<3;s++)
    {
        H_eff_Local[s]=0;
    }

    for(r=0;r<7;r++)
    {
        for(s=0;s<3;s++)
        {
            H_eff_Local[s] = H_eff_Local[s] - H_eff[s][r];
        }
    }

    free2DArray(H_eff, 3);
    return 1;
}

double Energy_Single(Energy_Variables *EV)
{
    double WEx = 0, WDM_Bloch = 0, WDM_Neel = 0, WZM = 0, WLd = 0, WAe = 0, Wme0 = 0, Wme1 = 0, W_Keff = 0;
    double dx, dy, dz, fo1, fo2, fo3;
    int m, n, l;
    double S11, S12, S13, S22, S23, S33, CST_Type;
    //CST Type=1: Sij is the Stress tensor, CST Type=0: Sij is the Strain tensor.
    double l1, l2, l3, q, lo1, lo2, lo3;

    double h, t, *dxy, ***M_Local, **CST, *Energy_Coefficient, C_Ex, C_DM_Bloch, C_DM_Neel, C_ZM, C_LD, C_ME, C_AA;
    enum Local_Point_Postion{Null, Origin, North, South, West, East, Lower, Upper};

    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;
    Energy_Terms ET;

    int i, *LocalEnergy;

    m   = EV->m;   l   = EV->l;
    n   = EV->n;   dxy = EV->dxy;
    h   = EV->h;   t   = EV->t;
    CST = EV->CST;

    LocalEnergy        = EV->LocalEnergy;
    M_Local            = EV->M_Local;
    Energy_Coefficient = EV->Energy_Coefficient;

    dx=dxy[0];  dy=dxy[1];  dz=dxy[2];

    double MxO, MxN, MxS, MxW, MxE, MxL, MxU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU;
    double W_Total;
    for(i=0;i<8;++i)
    {
        LPP=(Local_Point_Postion) i;
        switch(LPP)
        {
            case Origin :
            {
                MxO  = *M_Local[1][i];    MyO  = *M_Local[2][i];    MzO  = *M_Local[3][i];
                break;
            }

            case North :
            {
                MxN  = *M_Local[1][i];    MyN  = *M_Local[2][i];    MzN  = *M_Local[3][i];
                break;
            }

            case South :
            {
                MxS  = *M_Local[1][i];    MyS  = *M_Local[2][i];    MzS  = *M_Local[3][i];
                break;
            }

            case West :
            {
                MxW  = *M_Local[1][i];    MyW  = *M_Local[2][i];    MzW  = *M_Local[3][i];
                break;
            }

            case East :
            {
                MxE  = *M_Local[1][i];    MyE  = *M_Local[2][i];    MzE  = *M_Local[3][i];
                break;
            }

            case Lower :
            {
                MxL  = *M_Local[1][i];    MyL  = *M_Local[2][i];    MzL  = *M_Local[3][i];
                break;
            }


            case Upper :
            {
                MxU  = *M_Local[1][i];    MyU  = *M_Local[2][i];    MzU  = *M_Local[3][i];
                break;
            }

            default:
                break;
        }
    }

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch(ET)
        {
            case Exchange:
            {
                C_Ex = Energy_Coefficient[i];
                break;
            }

            case DM_Bloch:
            {
                C_DM_Bloch = Energy_Coefficient[i];
                break;
            }

            case DM_Neel:
            {
                C_DM_Neel = Energy_Coefficient[i];
                break;
            }

            case Zeeman:
            {
                C_ZM = Energy_Coefficient[i];
                break;
            }

            case Landau:
            {
                C_LD = Energy_Coefficient[i];
                break;
            }

            case Magnetoelastic:
            {
                C_ME = Energy_Coefficient[i];
                break;
            }

            case Axial_Anisotropy:
            {
                C_AA = Energy_Coefficient[i];
                break;
            }

            default:
            {
                break;
            }
        }
    }

    if(l==1)//2D Case
    {
        if(LocalEnergy[0]==1)
        {
            CST_Type=CST[6][0];
            if(CST_Type!=0)
            {
                S11 = CST[0][0];    S22 = CST[1][0];
                S33 = CST[2][0];    S12 = CST[3][0];
                S13 = CST[4][0];    S23 = CST[5][0];

                q   = q_Stress;     l1  = l1_Stress;
                l2  = l2_Stress;    l3  = l3_Stress;
                lo1 = lo1_Stress;   lo2 = lo2_Stress;
                lo3 = lo3_Stress;
            }else
            {
                S11 = CST[0][0];    S22 = CST[1][0];
                S33 = CST[2][0];    S12 = 2*CST[3][0];
                S13 = 2*CST[4][0];  S23 = 2*CST[5][0];

                q   = q_Strain;     l1  = l1_Strain;
                l2  = l2_Strain;    l3  = l3_Strain;
                lo1 = lo1_Strain;   lo2 = lo2_Strain;
                lo3 = lo3_Strain;
            }

            WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2);

            WDM_Bloch = (MzE - MzW) / (2 * dy) * MxO -
                        (MzS - MzN) / (2 * dx) * MyO +
                        ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;

            WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzS - MzN) / (2 * dx) * MxO + (MzE - MzW) / (2 * dy) * MyO -
                                                        ((MxS - MxN) / (2 * dx) + (MyE - MyW) / (2 * dy)) * MzO);

            WZM=MzO*2*h;

            W_Keff=k_eff*(1+k_0*(S11+S22))*MzO*MzO;

            WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);

            WAe = A_e*(pow((MxS - MxN) / (2 * dx), 2) + pow((MxE - MxW) / (2 * dy), 2) + pow((MyE - MyW) / (2 * dy), 2) + pow((MyS - MyN) / (2 * dx), 2));

            //WLd=0;

            Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                 (l1*S22+l2*S33)*pow(MyO,2)+
                 (l1*S33+l2*S11)*pow(MzO,2)+
                 l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                 q*WLd*(S11+S22+S33);

            fo1=S11*(MxE-MxW)/(2*dy)*MzO-
                S22*(MyS-MyN)/(2*dx)*MzO+
                S33*((MzS-MzN)/(2*dx)*MyO-(MzE-MzW)/(2*dy)*MxO);

            fo2=S11*((MzS-MzN)/(2*dx)*MyO-(MyS-MyN)/(2*dx)*MzO)+
                S22*((MxE-MxW)/(2*dy)*MzO-(MzE-MzW)/(2*dy)*MxO);

            fo3=S11*(-(MzE-MzW))/(2*dy)*MxO+
                S22*(MzS-MzN)/(2*dx)*MyO+
                S33*((MxE-MxW)/(2*dy)*MzO-(MyS-MyN)/(2*dx))*MzO;

            Wme1=lo1*fo1+lo2*fo2+lo3*fo3;

            W_Total = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                      C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
            //W_Total=2*WEx+2*WDM;
        }else
        {
            W_Total=0;
        }
    }else//3D Case
    {
        if(LocalEnergy[0]==1)
        {
            CST_Type=CST[6][0];
            if(CST_Type!=0)
            {
                S11 = CST[0][0];    S22 = CST[1][0];
                S33 = CST[2][0];    S12 = CST[3][0];
                S13 = CST[4][0];    S23 = CST[5][0];

                q   = q_Stress;     l1  = l1_Stress;
                l2  = l2_Stress;    l3  = l3_Stress;
                lo1 = lo1_Stress;   lo2 = lo2_Stress;
                lo3 = lo3_Stress;
            }else
            {
                S11 = CST[0][0];    S22 = CST[1][0];
                S33 = CST[2][0];    S12 = 2*CST[3][0];
                S13 = 2*CST[4][0];  S23 = 2*CST[5][0];

                q   = q_Strain;     l1  = l1_Strain;
                l2  = l2_Strain;    l3  = l3_Strain;
                lo1 = lo1_Strain;   lo2 = lo2_Strain;
                lo3 = lo3_Strain;
            }

            WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2)+
                pow((MxU-MxL)/(2*dz),2)+pow((MyU-MyL)/(2*dz),2)+
                pow((MzU-MzL)/(2*dz),2);

            WDM_Bloch = ((MzE - MzW) / (2 * dy) - (MyU - MyL) / (2 * dz)) * MxO +
                        ((MxU - MxL) / (2 * dz) - (MzS - MzN) / (2 * dx)) * MyO +
                        ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;

            WZM=MzO*2*h;

            W_Keff=k_eff*MzO*MzO;

            WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);


            Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                 (l1*S22+l2*S33)*pow(MyO,2)+
                 (l1*S33+l2*S11)*pow(MzO,2)+
                 l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                 q*WLd*(S11+S22+S33);

            fo1=S11*((MxE-MxW)/(2*dy)*MzO-
                    (MxU-MxL)/(2*dz)*MyO)+
                S22*((MyU-MyL)/(2*dz)*MxO-
                    (MyS-MyN)/(2*dx)*MzO)+
                S33*((MzS-MzN)/(2*dx)*MyO-
                    (MzE-MzW)/(2*dy)*MxO);

            fo2=S11*((MzS-MzN)/(2*dx)*MyO-
                    (MyS-MyN)/(2*dx)*MzO)+
                S22*((MxE-MxW)/(2*dy)*MzO-
                    (MzE-MzW)/(2*dy)*MxO)+
                S33*((MyU-MyL)/(2*dz)*MxO-
                    (MxU-MxL)/(2*dz)*MyO);

            fo3=S11*((MyU-MyL)/(2*dz)-
                    (MzE-MzW)/(2*dy))*MxO+
                S22*((MzS-MzN)/(2*dx)-
                    (MxU-MxL)/(2*dz))*MyO+
                S33*((MxE-MxW)/(2*dy)-
                    (MyS-MyN)/(2*dx))*MzO;

            Wme1=lo1*fo1+lo2*fo2+lo3*fo3;

            W_Total = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                      C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
        }else
        {
            W_Total=0;
        }
    }
    if(fabs(W_Total)<2)
    {
        //printf("1\n");
    }else
    {
        //printf("%0.1f\n", W_Total);
    }
    return W_Total;
}

double Energy_Local(Energy_Variables *EV)
{
    double WEx = 0, WDM_Bloch = 0, WDM_Neel = 0, WZM = 0, WLd = 0, WAe = 0, Wme0 = 0, Wme1 = 0, W_Keff = 0;
    double dx, dy, dz, fo1, fo2, fo3;
    int m, n, l;
    double S11, S12, S13, S22, S23, S33, CST_Type;
    //CST Type=1: Sij is the Stress tensor, CST Type=0: Sij is the Strain tensor.
    double l1, l2, l3, q, lo1, lo2, lo3;

    double h, t, *dxy, ***M_Local, **CST, *Energy_Coefficient, C_Ex, C_DM_Bloch, C_DM_Neel, C_ZM, C_LD, C_ME, C_AA;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    enum Local_Point_Postion_L{Null, O, N, S, W, E, L, U, NN, SS, WW, EE, LL, UU};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    //NN(i-2,j,k),SS(i+2,j,k), WW(i,j-2,k), EE(i,j+2,k), LL(i,j,k-2), UU(i,j,k+2);
    Local_Point_Postion LPP;
    Local_Point_Postion_L LPPL;
    Energy_Terms ET;
    int i, p, N_Local=0, *LocalEnergy;

    m   = EV->m;   l   = EV->l;
    n   = EV->n;   dxy = EV->dxy;
    h   = EV->h;   t   = EV->t;
    CST = EV->CST;
    LocalEnergy        = EV->LocalEnergy;
    M_Local            = EV->M_Local;
    Energy_Coefficient = EV->Energy_Coefficient;

    double MxO, MxN, MxS, MxW, MxE, MxL, MxU, MxNN, MxSS, MxWW, MxEE, MxLL, MxUU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU, MyNN, MySS, MyWW, MyEE, MyLL, MyUU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU, MzNN, MzSS, MzWW, MzEE, MzLL, MzUU;
    double W_Total, W_Origin=0, W_North=0, W_South=0, W_West=0, W_East=0, W_Lower=0, W_Upper=0;
    for(i=0;i<14;++i)
    {
        LPPL=(Local_Point_Postion_L) i;
        switch(LPPL)
        {
            case O :
            {
                MxO  = *M_Local[1][i];    MyO  = *M_Local[2][i];    MzO  = *M_Local[3][i];
                break;
            }

            case N :
            {
                MxN  = *M_Local[1][i];    MyN  = *M_Local[2][i];    MzN  = *M_Local[3][i];
                break;
            }

            case S :
            {
                MxS  = *M_Local[1][i];    MyS  = *M_Local[2][i];    MzS  = *M_Local[3][i];
                break;
            }

            case W :
            {
                MxW  = *M_Local[1][i];    MyW  = *M_Local[2][i];    MzW  = *M_Local[3][i];
                break;
            }

            case E :
            {
                MxE  = *M_Local[1][i];    MyE  = *M_Local[2][i];    MzE  = *M_Local[3][i];
                break;
            }

            case L :
            {
                MxL  = *M_Local[1][i];    MyL  = *M_Local[2][i];    MzL  = *M_Local[3][i];
                break;
            }


            case U :
            {
                MxU  = *M_Local[1][i];    MyU  = *M_Local[2][i];    MzU  = *M_Local[3][i];
                break;
            }

            case NN :
            {
                MxNN = *M_Local[1][i];    MyNN = *M_Local[2][i];    MzNN = *M_Local[3][i];
                break;
            }

            case SS :
            {
                MxSS = *M_Local[1][i];    MySS = *M_Local[2][i];    MzSS = *M_Local[3][i];
                break;
            }

            case WW :
            {
                MxWW = *M_Local[1][i];    MyWW = *M_Local[2][i];    MzWW = *M_Local[3][i];
                break;
            }

            case EE :
            {
                MxEE = *M_Local[1][i];    MyEE = *M_Local[2][i];    MzEE = *M_Local[3][i];
                break;
            }

            case LL :
            {
                MxLL = *M_Local[1][i];    MyLL = *M_Local[2][i];    MzLL = *M_Local[3][i];
                break;
            }

            case UU :
            {
                MxUU = *M_Local[1][i];    MyUU = *M_Local[2][i];    MzUU = *M_Local[3][i];
                break;
            }

            default:
                break;
        }
    }

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch(ET)
        {
            case Exchange:
            {
                C_Ex = Energy_Coefficient[i];
                break;
            }

            case DM_Bloch:
            {
                C_DM_Bloch = Energy_Coefficient[i];
                break;
            }

            case DM_Neel:
            {
                C_DM_Neel = Energy_Coefficient[i];
                break;
            }

            case Zeeman:
            {
                C_ZM = Energy_Coefficient[i];
                break;
            }

            case Landau:
            {
                C_LD = Energy_Coefficient[i];
                break;
            }

            case Magnetoelastic:
            {
                C_ME = Energy_Coefficient[i];
                break;
            }

            case Axial_Anisotropy:
            {
                C_AA = Energy_Coefficient[i];
                break;
            }

            default:
            {
                break;
            }
        }
    }

    dx=dxy[0];  dy=dxy[1];  dz=dxy[2];

    if(l==1)//2D Case
    {
        for(p=0;p<5;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                            pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                            pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2);

                        WDM_Bloch = (MzE - MzW) / (2 * dy) * MxO -
                                    (MzS - MzN) / (2 * dx) * MyO +
                                    ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzS - MzN) / (2 * dx) * MxO + (MzE - MzW) / (2 * dy) * MyO -
                                                                    ((MxS - MxN) / (2 * dx) + (MyE - MyW) / (2 * dy)) * MzO);

                        WZM=MzO*2*h;
                        W_Keff=k_eff*(1+k_0*(S11+S22))*MzO*MzO;
                        WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);
                        WAe = A_e*(pow((MxS - MxN) / (2 * dx), 2) + pow((MxE - MxW) / (2 * dy), 2) + pow((MyE - MyW) / (2 * dy), 2) + pow((MyS - MyN) / (2 * dx), 2));
                        //WLd=0;
                        Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                            (l1*S22+l2*S33)*pow(MyO,2)+
                            (l1*S33+l2*S11)*pow(MzO,2)+
                            l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                            q*WLd*(S11+S22+S33);
                        fo1=S11*(MxE-MxW)/(2*dy)*MzO-
                            S22*(MyS-MyN)/(2*dx)*MzO+
                            S33*((MzS-MzN)/(2*dx)*MyO-(MzE-MzW)/(2*dy)*MxO);
                        fo2=S11*((MzS-MzN)/(2*dx)*MyO-(MyS-MyN)/(2*dx)*MzO)+
                            S22*((MxE-MxW)/(2*dy)*MzO-(MzE-MzW)/(2*dy)*MxO);
                        fo3=S11*(-(MzE-MzW))/(2*dy)*MxO+
                            S22*(MzS-MzN)/(2*dx)*MyO+
                            S33*((MxE-MxW)/(2*dy)*MzO-(MyS-MyN)/(2*dx))*MzO;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Origin = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                   C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Origin=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxNN)/(2*dx),2)+pow((MyO-MyNN)/(2*dx),2)+
                            pow((MzO-MzNN)/(2*dx),2);

                        WDM_Bloch = -(MzO - MzNN) / (2 * dx) * MyN +
                                    (MyO - MyNN) / (2 * dx) * MzN;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzO - MzNN) / (2 * dx) * MxN -
                                                                    (MxO - MxNN) / (2 * dx) * MzN);

                        //WZM=MzN*2*h;
                        WZM=0.;
                        //WLd=pow(MxN,2)+pow(MyN,2)+pow(MzN,2);
                        WLd=0.;
                        WAe = A_e*(pow((MxO - MxNN) / (2 * dx), 2) + pow((MyO - MyNN) / (2 * dx), 2));
                        W_Keff=0;
/*
                        Wme0=(l1*S11+l2*S22)*pow(MxN,2)+
                            (l1*S22+l2*S33)*pow(MyN,2)+
                            (l1*S33+l2*S11)*pow(MzN,2)+
                            l3*(S12*MxN*MyN+S13*MxN*MzN+S23*MyN*MzN)+
                            q*WLd*(S11+S22+S33);
*/
                        Wme0=0.;
                        fo1=-S22*(MyO-MyNN)/(2*dx)*MzN+
                            S33*(MzO-MzNN)/(2*dx)*MyN;
                        fo2=S11*((MzO-MzNN)/(2*dx)*MyN-(MyO-MyNN)/(2*dx)*MzN);
                        fo3=S22*(MzO-MzNN)/(2*dx)*MyN+
                            S33*(-(MyO-MyNN)/(2*dx))*MzN;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_North = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_North=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxSS-MxO)/(2*dx),2)+pow((MySS-MyO)/(2*dx),2)+
                            pow((MzSS-MzO)/(2*dx),2);

                        WDM_Bloch = -(MzSS - MzO) / (2 * dx) * MyS +
                                    (MySS - MyO) / (2 * dx) * MzS;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzSS - MzO) / (2 * dx) * MxS -
                                                                    (MxSS - MxO) / (2 * dx) * MzS);

                        WZM=0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxSS - MxO) / (2 * dx), 2) + pow((MySS - MyO) / (2 * dx), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=-S22*(MySS-MyO)/(2*dx)*MzS+
                            S33*(MzSS-MzO)/(2*dx)*MyS;
                        fo2=S11*((MzSS-MzO)/(2*dx)*MyS-(MySS-MyO)/(2*dx)*MzS);
                        fo3=S22*(MzSS-MzO)/(2*dx)*MyS+
                            S33*(-(MySS-MyO)/(2*dx))*MzS;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_South = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_South=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxWW)/(2*dy),2)+pow((MyO-MyWW)/(2*dy),2)+
                            pow((MzO-MzWW)/(2*dy),2);

                        WDM_Bloch = (MzO - MzWW) / (2 * dy) * MxW +
                                    (-(MxO - MxWW) / (2 * dy)) * MzW;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzO - MzWW) / (2 * dy) * MyW -
                                                                    (MyO - MyWW) / (2 * dy) * MzW);

                        WZM = 0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxO - MxWW) / (2 * dy), 2) + pow((MyO - MyWW) / (2 * dy), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxO-MxWW)/(2*dy)*MzW+
                            S33*(-(MzO-MzWW)/(2*dy)*MxW);
                        fo2=S22*((MxO-MxWW)/(2*dy)*MzW-(MzO-MzWW)/(2*dy)*MxW);
                        fo3=S11*(-(MzO-MzWW))/(2*dy)*MxW+
                            S33*((MxO-MxWW)/(2*dy)*MzW);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_West = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_West=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxEE-MxO)/(2*dy),2)+pow((MyEE-MyO)/(2*dy),2)+
                            pow((MzEE-MzO)/(2*dy),2);

                        WDM_Bloch = (MzEE - MzO) / (2 * dy) * MxE +
                                    (-(MxEE - MxO) / (2 * dy)) * MzE;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzEE - MzO) / (2 * dy) * MyE -
                                                                    (MyEE - MyO) / (2 * dy) * MzE);
                        WZM=0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxEE - MxO) / (2 * dy), 2) + pow((MyEE - MyO) / (2 * dy), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxEE-MxO)/(2*dy)*MzE+
                            S33*(-(MzEE-MzO)/(2*dy)*MxE);
                        fo2=S22*((MxEE-MxO)/(2*dy)*MzE-(MzEE-MzO)/(2*dy)*MxE);
                        fo3=S11*(-(MzEE-MzO))/(2*dy)*MxE+
                            S33*((MxEE-MxO)/(2*dy)*MzE);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_East = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_East=0;
                    }
                    break;
                }
                default:
                    break;
            }
            N_Local=N_Local+LocalEnergy[p];
        }
    }else//3D Case
    {
        for(p=0;p<7;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
//////////////////////////////////////////////////////////////////////////////////////
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                            pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                            pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2)+
                            pow((MxU-MxL)/(2*dz),2)+pow((MyU-MyL)/(2*dz),2)+
                            pow((MzU-MzL)/(2*dz),2);
                        WDM_Bloch = ((MzE - MzW) / (2 * dy) - (MyU - MyL) / (2 * dz)) * MxO +
                                    ((MxU - MxL) / (2 * dz) - (MzS - MzN) / (2 * dx)) * MyO +
                                    ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;
                        WZM=MzO*2*h;
                        W_Keff=k_eff*MzO*MzO;
                        WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);
                        Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                            (l1*S22+l2*S33)*pow(MyO,2)+
                            (l1*S33+l2*S11)*pow(MzO,2)+
                            l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                            q*WLd*(S11+S22+S33);
                        fo1=S11*((MxE-MxW)/(2*dy)*MzO-
                                (MxU-MxL)/(2*dz)*MyO)+
                            S22*((MyU-MyL)/(2*dz)*MxO-
                                (MyS-MyN)/(2*dx)*MzO)+
                            S33*((MzS-MzN)/(2*dx)*MyO-
                                (MzE-MzW)/(2*dy)*MxO);
                        fo2=S11*((MzS-MzN)/(2*dx)*MyO-
                                (MyS-MyN)/(2*dx)*MzO)+
                            S22*((MxE-MxW)/(2*dy)*MzO-
                                (MzE-MzW)/(2*dy)*MxO)+
                            S33*((MyU-MyL)/(2*dz)*MxO-
                                (MxU-MxL)/(2*dz)*MyO);
                        fo3=S11*((MyU-MyL)/(2*dz)-
                                (MzE-MzW)/(2*dy))*MxO+
                            S22*((MzS-MzN)/(2*dx)-
                                (MxU-MxL)/(2*dz))*MyO+
                            S33*((MxE-MxW)/(2*dy)-
                                (MyS-MyN)/(2*dx))*MzO;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Origin = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                   C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Origin=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxNN)/(2*dx),2)+pow((MyO-MyNN)/(2*dx),2)+
                            pow((MzO-MzNN)/(2*dx),2);
                        WDM_Bloch = -(MzO - MzNN) / (2 * dx) * MyN +
                                    (MyO - MyNN) / (2 * dx) * MzN;
                        //WZM=MzN*2*h;
                        WZM=0.;
                        //WLd=pow(MxN,2)+pow(MyN,2)+pow(MzN,2);
                        WLd=0.;
                        W_Keff=0;
/*
                        Wme0=(l1*S11+l2*S22)*pow(MxN,2)+
                            (l1*S22+l2*S33)*pow(MyN,2)+
                            (l1*S33+l2*S11)*pow(MzN,2)+
                            l3*(S12*MxN*MyN+S13*MxN*MzN+S23*MyN*MzN)+
                            q*WLd*(S11+S22+S33);
*/
                        Wme0=0.;
                        fo1=-S22*(MyO-MyNN)/(2*dx)*MzN+
                            S33*(MzO-MzNN)/(2*dx)*MyN;
                        fo2=S11*((MzO-MzNN)/(2*dx)*MyN-(MyO-MyNN)/(2*dx)*MzN);
                        fo3=S22*(MzO-MzNN)/(2*dx)*MyN+
                            S33*(-(MyO-MyNN)/(2*dx))*MzN;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_North = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_North=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxSS-MxO)/(2*dx),2)+pow((MySS-MyO)/(2*dx),2)+
                            pow((MzSS-MzO)/(2*dx),2);
                        WDM_Bloch = -(MzSS - MzO) / (2 * dx) * MyS +
                                    (MySS - MyO) / (2 * dx) * MzS;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=-S22*(MySS-MyO)/(2*dx)*MzS+
                            S33*(MzSS-MzO)/(2*dx)*MyS;
                        fo2=S11*((MzSS-MzO)/(2*dx)*MyS-(MySS-MyO)/(2*dx)*MzS);
                        fo3=S22*(MzSS-MzO)/(2*dx)*MyS+
                            S33*(-(MySS-MyO)/(2*dx))*MzS;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_South = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_South=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxWW)/(2*dy),2)+pow((MyO-MyWW)/(2*dy),2)+
                            pow((MzO-MzWW)/(2*dy),2);
                        WDM_Bloch = (MzO - MzWW) / (2 * dy) * MxW +
                                    (-(MxO - MxWW) / (2 * dy)) * MzW;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxO-MxWW)/(2*dy)*MzW+
                            S33*(-(MzO-MzWW)/(2*dy)*MxW);
                        fo2=S22*((MxO-MxWW)/(2*dy)*MzW-(MzO-MzWW)/(2*dy)*MxW);
                        fo3=S11*(-(MzO-MzWW))/(2*dy)*MxW+
                            S33*((MxO-MxWW)/(2*dy)*MzW);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_West = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_West=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxEE-MxO)/(2*dy),2)+pow((MyEE-MyO)/(2*dy),2)+
                            pow((MzEE-MzO)/(2*dy),2);
                        WDM_Bloch = (MzEE - MzO) / (2 * dy) * MxE +
                                    (-(MxEE - MxO) / (2 * dy)) * MzE;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxEE-MxO)/(2*dy)*MzE+
                            S33*(-(MzEE-MzO)/(2*dy)*MxE);
                        fo2=S22*((MxEE-MxO)/(2*dy)*MzE-(MzEE-MzO)/(2*dy)*MxE);
                        fo3=S11*(-(MzEE-MzO))/(2*dy)*MxE+
                            S33*((MxEE-MxO)/(2*dy)*MzE);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_East = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_East=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case Lower:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxLL)/(2*dz),2)+pow((MyO-MyLL)/(2*dz),2)+
                            pow((MzO-MzLL)/(2*dz),2);
                        WDM_Bloch = (-(MyO - MyLL) / (2 * dz)) * MxL +
                                    ((MxO - MxLL) / (2 * dz)) * MyL;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(-(MxO-MxLL)/(2*dz)*MyL)+
                            S22*((MyO-MyLL)/(2*dz)*MxL);
                        fo2=S33*((MyO-MyLL)/(2*dz)*MxL-
                                (MxO-MxLL)/(2*dz)*MyL);
                        fo3=S11*((MyO-MyLL)/(2*dz))*MxL+
                            S22*(-(MxO-MxLL)/(2*dz))*MyL;;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Lower = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Lower=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case Upper:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxUU-MxO)/(2*dz),2)+pow((MyUU-MyO)/(2*dz),2)+
                            pow((MzUU-MzO)/(2*dz),2);
                        WDM_Bloch = (-(MyUU - MyO) / (2 * dz)) * MxU +
                                    ((MxUU - MxO) / (2 * dz)) * MyU;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(-(MxUU-MxO)/(2*dz)*MyU)+
                            S22*((MyUU-MyO)/(2*dz)*MxU);
                        fo2=S33*((MyUU-MyO)/(2*dz)*MxU-
                                (MxUU-MxO)/(2*dz)*MyU);
                        fo3=S11*((MyUU-MyO)/(2*dz))*MxU+
                            S22*(-(MxUU-MxO)/(2*dz))*MyU;;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Upper = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - C_AA * W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Upper=0;
                    }
                    break;
                }
                default:
                    break;
            }
            N_Local=N_Local+LocalEnergy[p];
        }
    }

    W_Total=W_Origin+W_North+W_South+W_West+W_East+W_Lower+W_Upper;

    if(N!=0)
    {
        //W_Total=W_Total/N_Local;
        W_Total=W_Total;
    }else
    {
        W_Total=0;
    }
    //printf("%lf\n",W_Total);

    return W_Total;
}

double Energy_Local_M0(Energy_Variables *EV)
{
    double WEx = 0, WDM_Bloch = 0, WDM_Neel = 0, WZM = 0, WLd = 0, WAe = 0, Wme0 = 0, Wme1 = 0, W_Keff = 0;
    double dx, dy, dz, fo1, fo2, fo3;
    int m, n, l;
    double S11, S12, S13, S22, S23, S33, CST_Type;
    //CST Type=1: Sij is the Stress tensor, CST Type=0: Sij is the Strain tensor.
    double l1, l2, l3, q, lo1, lo2, lo3;

    double h, t, *dxy, ***M_Local, **CST, *Energy_Coefficient, C_Ex, C_DM_Bloch, C_DM_Neel, C_ZM, C_LD, C_ME;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    enum Local_Point_Postion_L{Null, O, N, S, W, E, L, U, NN, SS, WW, EE, LL, UU};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    //NN(i-2,j,k),SS(i+2,j,k), WW(i,j-2,k), EE(i,j+2,k), LL(i,j,k-2), UU(i,j,k+2);
    Local_Point_Postion LPP;
    Local_Point_Postion_L LPPL;
    Energy_Terms ET;
    int i, p, N_Local=0, *LocalEnergy;

    m   = EV->m;   l   = EV->l;
    n   = EV->n;   dxy = EV->dxy;
    h   = EV->h;   t   = EV->t;
    CST = EV->CST;
    LocalEnergy        = EV->LocalEnergy;
    M_Local            = EV->M_Local;
    Energy_Coefficient = EV->Energy_Coefficient;

    double MxO, MxN, MxS, MxW, MxE, MxL, MxU, MxNN, MxSS, MxWW, MxEE, MxLL, MxUU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU, MyNN, MySS, MyWW, MyEE, MyLL, MyUU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU, MzNN, MzSS, MzWW, MzEE, MzLL, MzUU;
    double W_Total, W_Origin=0, W_North=0, W_South=0, W_West=0, W_East=0, W_Lower=0, W_Upper=0;
    for(i=0;i<14;++i)
    {
        LPPL=(Local_Point_Postion_L) i;
        switch(LPPL)
        {
            case O :
            {
                MxO  = EV->M0[0];   MyO  = EV->M0[1];    MzO  = EV->M0[2];
                break;
            }

            case N :
            {
                MxN  = *M_Local[1][i];    MyN  = *M_Local[2][i];    MzN  = *M_Local[3][i];
                break;
            }

            case S :
            {
                MxS  = *M_Local[1][i];    MyS  = *M_Local[2][i];    MzS  = *M_Local[3][i];
                break;
            }

            case W :
            {
                MxW  = *M_Local[1][i];    MyW  = *M_Local[2][i];    MzW  = *M_Local[3][i];
                break;
            }

            case E :
            {
                MxE  = *M_Local[1][i];    MyE  = *M_Local[2][i];    MzE  = *M_Local[3][i];
                break;
            }

            case L :
            {
                MxL  = *M_Local[1][i];    MyL  = *M_Local[2][i];    MzL  = *M_Local[3][i];
                break;
            }


            case U :
            {
                MxU  = *M_Local[1][i];    MyU  = *M_Local[2][i];    MzU  = *M_Local[3][i];
                break;
            }

            case NN :
            {
                MxNN = *M_Local[1][i];    MyNN = *M_Local[2][i];    MzNN = *M_Local[3][i];
                break;
            }

            case SS :
            {
                MxSS = *M_Local[1][i];    MySS = *M_Local[2][i];    MzSS = *M_Local[3][i];
                break;
            }

            case WW :
            {
                MxWW = *M_Local[1][i];    MyWW = *M_Local[2][i];    MzWW = *M_Local[3][i];
                break;
            }

            case EE :
            {
                MxEE = *M_Local[1][i];    MyEE = *M_Local[2][i];    MzEE = *M_Local[3][i];
                break;
            }

            case LL :
            {
                MxLL = *M_Local[1][i];    MyLL = *M_Local[2][i];    MzLL = *M_Local[3][i];
                break;
            }

            case UU :
            {
                MxUU = *M_Local[1][i];    MyUU = *M_Local[2][i];    MzUU = *M_Local[3][i];
                break;
            }

            default:
                break;
        }
    }

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch(ET)
        {
            case Exchange:
            {
                C_Ex = Energy_Coefficient[i];
                break;
            }

            case DM_Bloch:
            {
                C_DM_Bloch = Energy_Coefficient[i];
                break;
            }

            case DM_Neel:
            {
                C_DM_Neel = Energy_Coefficient[i];
                break;
            }

            case Zeeman:
            {
                C_ZM = Energy_Coefficient[i];
                break;
            }

            case Landau:
            {
                C_LD = Energy_Coefficient[i];
                break;
            }

            case Magnetoelastic:
            {
                C_ME = Energy_Coefficient[i];
                break;
            }

            default:
            {
                break;
            }
        }
    }

    dx=dxy[0];  dy=dxy[1];  dz=dxy[2];

    if(l==1)//2D Case
    {
        for(p=0;p<5;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                            pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                            pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2);

                        WDM_Bloch = (MzE - MzW) / (2 * dy) * MxO -
                                    (MzS - MzN) / (2 * dx) * MyO +
                                    ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzS - MzN) / (2 * dx) * MxO + (MzE - MzW) / (2 * dy) * MyO -
                                                                    ((MxS - MxN) / (2 * dx) + (MyE - MyW) / (2 * dy)) * MzO);

                        WZM=MzO*2*h;
                        W_Keff=k_eff*(1+k_0*(S11+S22))*MzO*MzO;
                        WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);
                        WAe = A_e*(pow((MxS - MxN) / (2 * dx), 2) + pow((MxE - MxW) / (2 * dy), 2) + pow((MyE - MyW) / (2 * dy), 2) + pow((MyS - MyN) / (2 * dx), 2));
                        //WLd=0;
                        Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                            (l1*S22+l2*S33)*pow(MyO,2)+
                            (l1*S33+l2*S11)*pow(MzO,2)+
                            l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                            q*WLd*(S11+S22+S33);
                        fo1=S11*(MxE-MxW)/(2*dy)*MzO-
                            S22*(MyS-MyN)/(2*dx)*MzO+
                            S33*((MzS-MzN)/(2*dx)*MyO-(MzE-MzW)/(2*dy)*MxO);
                        fo2=S11*((MzS-MzN)/(2*dx)*MyO-(MyS-MyN)/(2*dx)*MzO)+
                            S22*((MxE-MxW)/(2*dy)*MzO-(MzE-MzW)/(2*dy)*MxO);
                        fo3=S11*(-(MzE-MzW))/(2*dy)*MxO+
                            S22*(MzS-MzN)/(2*dx)*MyO+
                            S33*((MxE-MxW)/(2*dy)*MzO-(MyS-MyN)/(2*dx))*MzO;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Origin = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                   C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Origin=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxNN)/(2*dx),2)+pow((MyO-MyNN)/(2*dx),2)+
                            pow((MzO-MzNN)/(2*dx),2);

                        WDM_Bloch = -(MzO - MzNN) / (2 * dx) * MyN +
                                    (MyO - MyNN) / (2 * dx) * MzN;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzO - MzNN) / (2 * dx) * MxN -
                                                                    (MxO - MxNN) / (2 * dx) * MzN);

                        //WZM=MzN*2*h;
                        WZM=0.;
                        //WLd=pow(MxN,2)+pow(MyN,2)+pow(MzN,2);
                        WLd=0.;
                        WAe = A_e*(pow((MxO - MxNN) / (2 * dx), 2) + pow((MyO - MyNN) / (2 * dx), 2));
                        W_Keff=0;
/*
                        Wme0=(l1*S11+l2*S22)*pow(MxN,2)+
                            (l1*S22+l2*S33)*pow(MyN,2)+
                            (l1*S33+l2*S11)*pow(MzN,2)+
                            l3*(S12*MxN*MyN+S13*MxN*MzN+S23*MyN*MzN)+
                            q*WLd*(S11+S22+S33);
*/
                        Wme0=0.;
                        fo1=-S22*(MyO-MyNN)/(2*dx)*MzN+
                            S33*(MzO-MzNN)/(2*dx)*MyN;
                        fo2=S11*((MzO-MzNN)/(2*dx)*MyN-(MyO-MyNN)/(2*dx)*MzN);
                        fo3=S22*(MzO-MzNN)/(2*dx)*MyN+
                            S33*(-(MyO-MyNN)/(2*dx))*MzN;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_North = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_North=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxSS-MxO)/(2*dx),2)+pow((MySS-MyO)/(2*dx),2)+
                            pow((MzSS-MzO)/(2*dx),2);

                        WDM_Bloch = -(MzSS - MzO) / (2 * dx) * MyS +
                                    (MySS - MyO) / (2 * dx) * MzS;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzSS - MzO) / (2 * dx) * MxS -
                                                                    (MxSS - MxO) / (2 * dx) * MzS);

                        WZM=0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxSS - MxO) / (2 * dx), 2) + pow((MySS - MyO) / (2 * dx), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=-S22*(MySS-MyO)/(2*dx)*MzS+
                            S33*(MzSS-MzO)/(2*dx)*MyS;
                        fo2=S11*((MzSS-MzO)/(2*dx)*MyS-(MySS-MyO)/(2*dx)*MzS);
                        fo3=S22*(MzSS-MzO)/(2*dx)*MyS+
                            S33*(-(MySS-MyO)/(2*dx))*MzS;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_South = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_South=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxWW)/(2*dy),2)+pow((MyO-MyWW)/(2*dy),2)+
                            pow((MzO-MzWW)/(2*dy),2);

                        WDM_Bloch = (MzO - MzWW) / (2 * dy) * MxW +
                                    (-(MxO - MxWW) / (2 * dy)) * MzW;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzO - MzWW) / (2 * dy) * MyW -
                                                                    (MyO - MyWW) / (2 * dy) * MzW);

                        WZM = 0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxO - MxWW) / (2 * dy), 2) + pow((MyO - MyWW) / (2 * dy), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxO-MxWW)/(2*dy)*MzW+
                            S33*(-(MzO-MzWW)/(2*dy)*MxW);
                        fo2=S22*((MxO-MxWW)/(2*dy)*MzW-(MzO-MzWW)/(2*dy)*MxW);
                        fo3=S11*(-(MzO-MzWW))/(2*dy)*MxW+
                            S33*((MxO-MxWW)/(2*dy)*MzW);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_West = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_West=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxEE-MxO)/(2*dy),2)+pow((MyEE-MyO)/(2*dy),2)+
                            pow((MzEE-MzO)/(2*dy),2);

                        WDM_Bloch = (MzEE - MzO) / (2 * dy) * MxE +
                                    (-(MxEE - MxO) / (2 * dy)) * MzE;

                        WDM_Neel = (1 + (S11 + S22) * d_epsilon) * ((MzEE - MzO) / (2 * dy) * MyE -
                                                                    (MyEE - MyO) / (2 * dy) * MzE);
                        WZM=0.;
                        WLd=0.;
                        WAe = A_e*(pow((MxEE - MxO) / (2 * dy), 2) + pow((MyEE - MyO) / (2 * dy), 2));
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxEE-MxO)/(2*dy)*MzE+
                            S33*(-(MzEE-MzO)/(2*dy)*MxE);
                        fo2=S22*((MxEE-MxO)/(2*dy)*MzE-(MzEE-MzO)/(2*dy)*MxE);
                        fo3=S11*(-(MzEE-MzO))/(2*dy)*MxE+
                            S33*((MxEE-MxO)/(2*dy)*MzE);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_East = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_East=0;
                    }
                    break;
                }
                default:
                    break;
            }
            N_Local=N_Local+LocalEnergy[p];
        }
    }else//3D Case
    {
        for(p=0;p<7;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
//////////////////////////////////////////////////////////////////////////////////////
                case Origin:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxS-MxN)/(2*dx),2)+pow((MyS-MyN)/(2*dx),2)+
                            pow((MzS-MzN)/(2*dx),2)+pow((MxE-MxW)/(2*dy),2)+
                            pow((MyE-MyW)/(2*dy),2)+pow((MzE-MzW)/(2*dy),2)+
                            pow((MxU-MxL)/(2*dz),2)+pow((MyU-MyL)/(2*dz),2)+
                            pow((MzU-MzL)/(2*dz),2);
                        WDM_Bloch = ((MzE - MzW) / (2 * dy) - (MyU - MyL) / (2 * dz)) * MxO +
                                    ((MxU - MxL) / (2 * dz) - (MzS - MzN) / (2 * dx)) * MyO +
                                    ((MyS - MyN) / (2 * dx) - (MxE - MxW) / (2 * dy)) * MzO;
                        WZM=MzO*2*h;
                        W_Keff=k_eff*MzO*MzO;
                        WLd=pow(MxO,2)+pow(MyO,2)+pow(MzO,2);
                        Wme0=(l1*S11+l2*S22)*pow(MxO,2)+
                            (l1*S22+l2*S33)*pow(MyO,2)+
                            (l1*S33+l2*S11)*pow(MzO,2)+
                            l3*(S12*MxO*MyO+S13*MxO*MzO+S23*MyO*MzO)+
                            q*WLd*(S11+S22+S33);
                        fo1=S11*((MxE-MxW)/(2*dy)*MzO-
                                (MxU-MxL)/(2*dz)*MyO)+
                            S22*((MyU-MyL)/(2*dz)*MxO-
                                (MyS-MyN)/(2*dx)*MzO)+
                            S33*((MzS-MzN)/(2*dx)*MyO-
                                (MzE-MzW)/(2*dy)*MxO);
                        fo2=S11*((MzS-MzN)/(2*dx)*MyO-
                                (MyS-MyN)/(2*dx)*MzO)+
                            S22*((MxE-MxW)/(2*dy)*MzO-
                                (MzE-MzW)/(2*dy)*MxO)+
                            S33*((MyU-MyL)/(2*dz)*MxO-
                                (MxU-MxL)/(2*dz)*MyO);
                        fo3=S11*((MyU-MyL)/(2*dz)-
                                (MzE-MzW)/(2*dy))*MxO+
                            S22*((MzS-MzN)/(2*dx)-
                                (MxU-MxL)/(2*dz))*MyO+
                            S33*((MxE-MxW)/(2*dy)-
                                (MyS-MyN)/(2*dx))*MzO;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Origin = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                   C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Origin=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case North:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxNN)/(2*dx),2)+pow((MyO-MyNN)/(2*dx),2)+
                            pow((MzO-MzNN)/(2*dx),2);
                        WDM_Bloch = -(MzO - MzNN) / (2 * dx) * MyN +
                                    (MyO - MyNN) / (2 * dx) * MzN;
                        //WZM=MzN*2*h;
                        WZM=0.;
                        //WLd=pow(MxN,2)+pow(MyN,2)+pow(MzN,2);
                        WLd=0.;
                        W_Keff=0;
/*
                        Wme0=(l1*S11+l2*S22)*pow(MxN,2)+
                            (l1*S22+l2*S33)*pow(MyN,2)+
                            (l1*S33+l2*S11)*pow(MzN,2)+
                            l3*(S12*MxN*MyN+S13*MxN*MzN+S23*MyN*MzN)+
                            q*WLd*(S11+S22+S33);
*/
                        Wme0=0.;
                        fo1=-S22*(MyO-MyNN)/(2*dx)*MzN+
                            S33*(MzO-MzNN)/(2*dx)*MyN;
                        fo2=S11*((MzO-MzNN)/(2*dx)*MyN-(MyO-MyNN)/(2*dx)*MzN);
                        fo3=S22*(MzO-MzNN)/(2*dx)*MyN+
                            S33*(-(MyO-MyNN)/(2*dx))*MzN;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_North = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_North=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case South:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxSS-MxO)/(2*dx),2)+pow((MySS-MyO)/(2*dx),2)+
                            pow((MzSS-MzO)/(2*dx),2);
                        WDM_Bloch = -(MzSS - MzO) / (2 * dx) * MyS +
                                    (MySS - MyO) / (2 * dx) * MzS;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=-S22*(MySS-MyO)/(2*dx)*MzS+
                            S33*(MzSS-MzO)/(2*dx)*MyS;
                        fo2=S11*((MzSS-MzO)/(2*dx)*MyS-(MySS-MyO)/(2*dx)*MzS);
                        fo3=S22*(MzSS-MzO)/(2*dx)*MyS+
                            S33*(-(MySS-MyO)/(2*dx))*MzS;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_South = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_South=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case West:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxWW)/(2*dy),2)+pow((MyO-MyWW)/(2*dy),2)+
                            pow((MzO-MzWW)/(2*dy),2);
                        WDM_Bloch = (MzO - MzWW) / (2 * dy) * MxW +
                                    (-(MxO - MxWW) / (2 * dy)) * MzW;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxO-MxWW)/(2*dy)*MzW+
                            S33*(-(MzO-MzWW)/(2*dy)*MxW);
                        fo2=S22*((MxO-MxWW)/(2*dy)*MzW-(MzO-MzWW)/(2*dy)*MxW);
                        fo3=S11*(-(MzO-MzWW))/(2*dy)*MxW+
                            S33*((MxO-MxWW)/(2*dy)*MzW);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_West = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_West=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case East:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxEE-MxO)/(2*dy),2)+pow((MyEE-MyO)/(2*dy),2)+
                            pow((MzEE-MzO)/(2*dy),2);
                        WDM_Bloch = (MzEE - MzO) / (2 * dy) * MxE +
                                    (-(MxEE - MxO) / (2 * dy)) * MzE;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(MxEE-MxO)/(2*dy)*MzE+
                            S33*(-(MzEE-MzO)/(2*dy)*MxE);
                        fo2=S22*((MxEE-MxO)/(2*dy)*MzE-(MzEE-MzO)/(2*dy)*MxE);
                        fo3=S11*(-(MzEE-MzO))/(2*dy)*MxE+
                            S33*((MxEE-MxO)/(2*dy)*MzE);
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_East = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                 C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_East=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case Lower:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxO-MxLL)/(2*dz),2)+pow((MyO-MyLL)/(2*dz),2)+
                            pow((MzO-MzLL)/(2*dz),2);
                        WDM_Bloch = (-(MyO - MyLL) / (2 * dz)) * MxL +
                                    ((MxO - MxLL) / (2 * dz)) * MyL;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(-(MxO-MxLL)/(2*dz)*MyL)+
                            S22*((MyO-MyLL)/(2*dz)*MxL);
                        fo2=S33*((MyO-MyLL)/(2*dz)*MxL-
                                (MxO-MxLL)/(2*dz)*MyL);
                        fo3=S11*((MyO-MyLL)/(2*dz))*MxL+
                            S22*(-(MxO-MxLL)/(2*dz))*MyL;;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Lower = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Lower=0;
                    }
                    break;
                }
//////////////////////////////////////////////////////////////////////////////////////
                case Upper:
                {
                    if(LocalEnergy[p]==1)
                    {
                        CST_Type=CST[6][p];
                        if(CST_Type!=0)
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = CST[3][p];
                            S13 = CST[4][p];    S23 = CST[5][p];
                            q   = q_Stress;     l1  = l1_Stress;
                            l2  = l2_Stress;    l3  = l3_Stress;
                            lo1 = lo1_Stress;   lo2 = lo2_Stress;
                            lo3 = lo3_Stress;
                        }else
                        {
                            S11 = CST[0][p];    S22 = CST[1][p];
                            S33 = CST[2][p];    S12 = 2*CST[3][p];
                            S13 = 2*CST[4][p];  S23 = 2*CST[5][p];
                            q   = q_Strain;     l1  = l1_Strain;
                            l2  = l2_Strain;    l3  = l3_Strain;
                            lo1 = lo1_Strain;   lo2 = lo2_Strain;
                            lo3 = lo3_Strain;
                        }
                        WEx=pow((MxUU-MxO)/(2*dz),2)+pow((MyUU-MyO)/(2*dz),2)+
                            pow((MzUU-MzO)/(2*dz),2);
                        WDM_Bloch = (-(MyUU - MyO) / (2 * dz)) * MxU +
                                    ((MxUU - MxO) / (2 * dz)) * MyU;
                        WZM=0.;
                        WLd=0.;
                        W_Keff=0;
                        Wme0=0.;
                        fo1=S11*(-(MxUU-MxO)/(2*dz)*MyU)+
                            S22*((MyUU-MyO)/(2*dz)*MxU);
                        fo2=S33*((MyUU-MyO)/(2*dz)*MxU-
                                (MxUU-MxO)/(2*dz)*MyU);
                        fo3=S11*((MyUU-MyO)/(2*dz))*MxU+
                            S22*(-(MxUU-MxO)/(2*dz))*MyU;;
                        Wme1=lo1*fo1+lo2*fo2+lo3*fo3;
                        W_Upper = C_Ex * WEx + C_DM_Bloch * 2 * WDM_Bloch + C_DM_Neel * 2 * WDM_Neel - C_ZM * WZM - W_Keff +
                                  C_LD * (t * WLd + WLd * WLd) + C_ME * (Wme0 + Wme1) + WAe;
                    }else
                    {
                        W_Upper=0;
                    }
                    break;
                }
                default:
                    break;
            }
            N_Local=N_Local+LocalEnergy[p];
        }
    }

    W_Total=W_Origin+W_North+W_South+W_West+W_East+W_Lower+W_Upper;

    if(N!=0)
    {
        //W_Total=W_Total/N_Local;
        W_Total=W_Total;
    }else
    {
        W_Total=0;
    }
    //printf("%lf\n",W);
    return W_Total;
}

int Continue_Modify(Energy_Variables *EV)
{
    int l;
    double ***M_Local;
    enum Local_Point_Postion{Null, Origin, North, South, West, East, Lower, Upper};

    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;

    int i, *LocalEnergy;
    l   = EV->l;
    LocalEnergy = EV->LocalEnergy;
    M_Local     = EV->M_Local;

    double MxO, MxN, MxS, MxW, MxE, MxL, MxU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU;

    double Mx_Average, My_Average, Mz_Average, Error_Range, Mx_Difference, My_Difference, Mz_Difference;
    for(i=0;i<8;++i)
    {
        LPP=(Local_Point_Postion) i;
        switch(LPP)
        {
            case Origin :
            {
                MxO  = *M_Local[1][i];    MyO  = *M_Local[2][i];    MzO  = *M_Local[3][i];
                break;
            }

            case North :
            {
                MxN  = *M_Local[1][i];    MyN  = *M_Local[2][i];    MzN  = *M_Local[3][i];
                break;
            }

            case South :
            {
                MxS  = *M_Local[1][i];    MyS  = *M_Local[2][i];    MzS  = *M_Local[3][i];
                break;
            }

            case West :
            {
                MxW  = *M_Local[1][i];    MyW  = *M_Local[2][i];    MzW  = *M_Local[3][i];
                break;
            }

            case East :
            {
                MxE  = *M_Local[1][i];    MyE  = *M_Local[2][i];    MzE  = *M_Local[3][i];
                break;
            }

            case Lower :
            {
                MxL  = *M_Local[1][i];    MyL  = *M_Local[2][i];    MzL  = *M_Local[3][i];
                break;
            }


            case Upper :
            {
                MxU  = *M_Local[1][i];    MyU  = *M_Local[2][i];    MzU  = *M_Local[3][i];
                break;
            }

            default:
                break;
        }
    }

    if(l==1)
    {
        Mx_Average = (MxN + MxS + MxW + MxE)/4;
        My_Average = (MyN + MyS + MyW + MyE)/4;
        Mz_Average = (MzN + MzS + MzW + MzE)/4;
    }else
    {
        Mx_Average = (MxN + MxS + MxW + MxE + MxL + MxU)/6;
        My_Average = (MyN + MyS + MyW + MyE + MyL + MyU)/6;
        Mz_Average = (MzN + MzS + MzW + MzE + MzL + MzU)/6;
    }

    Error_Range = EV->Continuity_Modify_Coefficient*EV->dxy[0];
/////////////////////////////////////////////////////////////
    do{
        Mx_Difference = fabs(Mx_Average-MxO)/(fabs(Mx_Average)+fabs(MxO));
        if(Mx_Difference>=Error_Range)
        {
            MxO = (MxO + Mx_Average)/2;
        }
    }while(Mx_Difference>=Error_Range);
//////////////////////////////////////////////////////////////
    do{
        My_Difference = fabs(My_Average-MyO)/(fabs(My_Average)+fabs(MyO));
        if(My_Difference>=Error_Range)
        {
            MyO = (MyO + My_Average)/2;
        }
    }while(My_Difference>=Error_Range);
///////////////////////////////////////////////////////////
    do{
        Mz_Difference = fabs(Mz_Average-MzO)/(fabs(Mz_Average)+fabs(MzO));
        if(Mz_Difference>=Error_Range)
        {
            MzO = (MzO + Mz_Average)/2;
        }
    }while(Mz_Difference>=Error_Range);

    if(LocalEnergy[0]==1)
    {
        *M_Local[1][1] = MxO;
        *M_Local[2][1] = MyO;
        *M_Local[3][1] = MzO;
    }

    return 1;
}

double Skyrmion_Number_Single(Energy_Variables *EV)
{
    double *dxy, ***M_Local, M,  N, dx, dy;
    double nxO, nyO, nzO,
           nxN, nyN, nzN,
           nxS, nyS, nzS,
           nxW, nyW, nzW,
           nxE, nyE, nzE;
    int p, *LocalEnergy;
    enum Local_Point_Postion{Null, Origin, North, South, West, East};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;
    M_Local     = EV->M_Local;
    dxy         = EV->dxy;
    LocalEnergy = EV->LocalEnergy;
    dx=dxy[0];  dy=dxy[1];

    if(LocalEnergy[0]==1)
    {
        for(p=0;p<6;++p)
        {
            LPP=(Local_Point_Postion) p;
            switch(LPP)
            {
                case Origin:
                {
                    M   = sqrt(pow(*M_Local[1][p],2)+pow(*M_Local[2][p],2)+pow(*M_Local[3][p],2));
                    if(M==0)
                    {
                        M = 1;
                    }
                    nxO = *M_Local[1][p]/M;
                    nyO = *M_Local[2][p]/M;
                    nzO = *M_Local[3][p]/M;
                    break;
                }

                case North:
                {
                    M   = sqrt(pow(*M_Local[1][p],2)+pow(*M_Local[2][p],2)+pow(*M_Local[3][p],2));
                    if(M==0)
                    {
                        M = 1;
                    }
                    nxN = *M_Local[1][p]/M;
                    nyN = *M_Local[2][p]/M;
                    nzN = *M_Local[3][p]/M;
                    break;
                }

                case South:
                {
                    M   = sqrt(pow(*M_Local[1][p],2)+pow(*M_Local[2][p],2)+pow(*M_Local[3][p],2));
                    if(M==0)
                    {
                        M = 1;
                    }
                    nxS = *M_Local[1][p]/M;
                    nyS = *M_Local[2][p]/M;
                    nzS = *M_Local[3][p]/M;
                    break;
                }

                case West:
                {
                    M   = sqrt(pow(*M_Local[1][p],2)+pow(*M_Local[2][p],2)+pow(*M_Local[3][p],2));
                    if(M==0)
                    {
                        M = 1;
                    }
                    nxW = *M_Local[1][p]/M;
                    nyW = *M_Local[2][p]/M;
                    nzW = *M_Local[3][p]/M;
                    break;
                }

                case East:
                {
                    M   = sqrt(pow(*M_Local[1][p],2)+pow(*M_Local[2][p],2)+pow(*M_Local[3][p],2));
                    if(M==0)
                    {
                        M = 1;
                    }
                    nxE = *M_Local[1][p]/M;
                    nyE = *M_Local[2][p]/M;
                    nzE = *M_Local[3][p]/M;
                    break;
                }

                default:
                    break;
            }
        }
        N  = 1/(4*pi)*dx*dy*(nxO*((nyS-nyN)/(2*dx)*(nzE-nzW)/(2*dy)-(nzS-nzN)/(2*dx)*(nyE-nyW)/(2*dy))+
                             nyO*((nzS-nzN)/(2*dx)*(nxE-nxW)/(2*dy)-(nxS-nxN)/(2*dx)*(nzE-nzW)/(2*dy))+
                             nzO*((nxS-nxN)/(2*dx)*(nyE-nyW)/(2*dy)-(nyS-nyN)/(2*dx)*(nxE-nxW)/(2*dy)));
        N=-N;
    }else
    {
        N=0;
    }
    if(fabs(N)<2)
    {
        //printf("1\n");
    }else
    {
        //printf("%0.1f\n", N);
    }
    

    return N;

}

int Energy_Distribution(Energy_Variables_Distribution *EVD)
{
    int i, j, k, p, N, ***Local_Points_1D, ***Position_3DTo1D;
    int I, J, K, v, Boundary_Type, N_Total=0;
    double W_Average=0, SKN=0, SKN_Area=0, SKN_Positif=0;
    Energy_Variables EV;
    EV.h = EVD->h;    EV.t = EVD->t;
    EV.m = EVD->m;    EV.n = EVD->n;
    EV.l = EVD->l;
    EV.M_Local     = Make2DArrayPointer(4,8);
    EV.CST         = Make2DArray(7,1);
    EV.LocalEnergy = Make1DArrayinteger(1);
    EV.dxy         = Make1DArray(3);
    EV.dxy[0]      = EVD->dxy[0];
    EV.dxy[1]      = EVD->dxy[1];
    EV.dxy[2]      = EVD->dxy[2];
    EV.Energy_Coefficient = Make1DArray(N_Energy_Terms);
    for (i = 0; i < N_Energy_Terms; ++i)
    {
        EV.Energy_Coefficient[i] = EVD->Energy_Coefficient[i];
    }

    double ***Mx, ***My, ***Mz, *MxV, *MyV, *MzV, ****CST;
    Mx  = EVD->Mx;   My  = EVD->My;
    Mz  = EVD->Mz;   MxV = EVD->MxV;
    MyV = EVD->MyV;  MzV = EVD->MzV;
    CST = EVD->CST;
    Local_Points_1D = EVD->Local_Points_1D;
    Position_3DTo1D = EVD->Position_3DTo1D;

    //enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    //Local_Point_Postion LPP;
    for(i=0;i<EVD->m;++i)
    {
        for(j=0;j<EVD->n;++j)
        {
            for(k=0;k<EVD->l;++k)
            {
                N = Position_3DTo1D[i][j][k];
                EV.LocalEnergy[0] = Local_Points_1D[N][0][5];
                if(EV.LocalEnergy[0]==1)
                {
                    N_Total=N_Total+1;
                }
                //printf("%d\t%d\t%d\t%d\n",i,j,k,EV.LocalEnergy[0]);
                for(p=0;p<7;++p)
                {
                    Boundary_Type = Local_Points_1D[N][p][0];
                    I             = Local_Points_1D[N][p][1];
                    J             = Local_Points_1D[N][p][2];
                    K             = Local_Points_1D[N][p][3];
                    v             = Local_Points_1D[N][p][4];
                    if(Boundary_Type!=0)
                    {
                        EV.M_Local[1][p+1]=&Mx[I][J][K];
                        EV.M_Local[2][p+1]=&My[I][J][K];
                        EV.M_Local[3][p+1]=&Mz[I][J][K];
                    }else
                    {
                        EV.M_Local[1][p+1]=&MxV[v];
                        EV.M_Local[2][p+1]=&MyV[v];
                        EV.M_Local[3][p+1]=&MzV[v];
                    }
                    EV.CST[p][0]=CST[p][i][j][k];
                    //printf("%d\t%lf\t%lf\t%lf\n",p,*EV.M_Local[1][p+1],*EV.M_Local[2][p+1],*EV.M_Local[3][p+1]);
                    //printf("%d\t%p\t%p\t%p\n",p,EV.M_Local[1][p+1],EV.M_Local[2][p+1],EV.M_Local[3][p+1]);
                }

                EVD->w[i][j][k]     = Energy_Single(&EV);
                EVD->SKN_D[i][j][k] = Skyrmion_Number_Single(&EV);
                W_Average           = W_Average+EVD->w[i][j][k];
                if(Skyrmion_Number_Single(&EV)>0)
                {
                    SKN_Area        = SKN_Area+EVD->dxy[0]*EVD->dxy[0];
                    SKN_Positif     = SKN_Positif+Skyrmion_Number_Single(&EV);
                }
                SKN                 = SKN+Skyrmion_Number_Single(&EV);
                //if(fabs(Skyrmion_Number_Single(&EV))>2)
                {
                    //printf("%d\t%d\t%d\t%lf\n",i,j,k,Skyrmion_Number_Single(&EV));
                    //printf("\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\n");
                }
            }
        }
    }
    EVD->W_Average=W_Average/N_Total;
    //EVD->W_Average=W_Average;
    EVD->SKN=SKN;
    EVD->SKN_Area=SKN_Area;
    EVD->SKN_Positif=SKN_Positif;

    free2DArrayPointer(EV.M_Local,4);
    free2DArray(EV.CST,7);
    free(EV.LocalEnergy);
    free(EV.dxy);
    free(EV.Energy_Coefficient);
    return 1;
}

double Skyrmion_Number_Incircle(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, p, i, j, k, v, i_L, j_L, k_L;
    int xc, yc, zc, rn, R_In, Edges, N, Reduction_Radius;
    double *dxy, ***Mx, ***My, ***Mz, *MxV, *MyV, *MzV, SKN;
    xn  = IPN->xn;  yn  = IPN->yn;
    zn  = IPN->zn;  dxy = IPN->dxy;

    m  = 2*xn+3; n  = 2*yn+3;
    xc = xn+1;   yc = yn+1;
    if(zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*zn+3;
        zc = zn+1;
    }

    if(xn>yn)
    {
        rn=yn;
    }else
    {
        rn=xn;
    }


    if(IPN->Boundary_Shape_Type==0)
    {
        R_In = rn;
    }else if(IPN->Boundary_Shape_Type==1)
    {
        R_In = rn;
        printf("Reduction Radius=");
        scanf("%d",&Reduction_Radius);
        R_In = rn-Reduction_Radius;
    }else if(IPN->Boundary_Shape_Type==2)
    {
        R_In = rn;
    }else
    {
        printf("Reduction Radius=");
        scanf("%d",&Reduction_Radius);
        Edges = IPN->Boundary_Shape_Type;
        R_In  = rn*cos(pi/Edges)-Reduction_Radius;
    }
    printf("%d\n",R_In);

    Boundary_Shape BS0, *BS;
    Calculate_Points_Parameters CPP0, *CPP;

    BS=&BS0;
    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);

    CPP=&CPP0;
    CPP->Position_3DTo1D  = Make3DArrayinteger(m,n,l);
    CPP->Local_Points_1D  = Make3DArrayinteger(BS->Calculate_Points,7,7);
    CPP->Virtual_Position = Make1DArrayinteger(BS->Virtual_Points);
    CPP->Local_Points_1D_LocalEnergy       = Make4DArrayinteger(BS->Calculate_Points,7,7,7);
    CalculPointsParameters_Initial(CPP,BS);

    Mx  = Make3DArray(m,n,l);
    My  = Make3DArray(m,n,l);
    Mz  = Make3DArray(m,n,l);
    MxV = Make1DArray(BS->Virtual_Points);
    MyV = Make1DArray(BS->Virtual_Points);
    MzV = Make1DArray(BS->Virtual_Points);
    InitialM_3D(IPN,Mx,My,Mz,MxV,MyV,MzV,BS);

    Energy_Variables EV;
    EV.dxy=Make1DArray(3);
    EV.dxy[0] = dxy[0];    EV.dxy[1] = dxy[1];
    EV.dxy[2] = dxy[2];
    EV.LocalEnergy = Make1DArrayinteger(1);
    EV.M_Local     = Make2DArrayPointer(4,6);

    SKN = 0.;
    for(k=0;k<l;++k)
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(sqrt((i-xc)*(i-xc)+(j-yc)*(j-yc))<=R_In)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        N = CPP->Position_3DTo1D[i][j][k];
                        EV.LocalEnergy[0] = CPP->Local_Points_1D[N][0][5];
                        for(p=0;p<5;++p)
                        {
                            if(CPP->Local_Points_1D[N][p][0]!=0)
                            {
                                i_L = CPP->Local_Points_1D[N][p][1];
                                j_L = CPP->Local_Points_1D[N][p][2];
                                k_L = CPP->Local_Points_1D[N][p][3];

                                EV.M_Local[1][p+1] = &Mx[i_L][j_L][k_L];
                                EV.M_Local[2][p+1] = &My[i_L][j_L][k_L];
                                EV.M_Local[3][p+1] = &Mz[i_L][j_L][k_L];
                            }else
                            {
                                v = CPP->Local_Points_1D[N][p][4];

                                EV.M_Local[1][p+1] = &MxV[v];
                                EV.M_Local[2][p+1] = &MyV[v];
                                EV.M_Local[3][p+1] = &MzV[v];
                            }
                        }
                        SKN = SKN+Skyrmion_Number_Single(&EV);
                    }
                }
            }
        }
    }
    free3DArray(Mx,m,n);
    free3DArray(My,m,n);
    free3DArray(Mz,m,n);
    free(MxV);
    free(MyV);
    free(MzV);
    free3DArrayinteger(BS->Boundary,m,n);
    free4DArrayinteger(BS->LocalEnergy,m,n,l);
    free(CPP->Virtual_Position);
    free4DArrayinteger(CPP->Local_Points_1D_LocalEnergy,BS->Calculate_Points,7,7);
    free3DArrayinteger(CPP->Local_Points_1D,BS->Calculate_Points,7);
    free3DArrayinteger(CPP->Position_3DTo1D,m,n);
    free(EV.dxy);
    free(EV.LocalEnergy);
    free2DArrayPointer(EV.M_Local,4);

    return SKN;
}
