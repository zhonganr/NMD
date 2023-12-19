#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <windows.h>
#include <pthread.h>
#include <BoundaryShape.h>
#include <Variable.h>
#include <Array.h>
#include <Energy.h>
#include <InOut.h>
#include <Elastic_Field.h>
#include <OPGL.h>
#include <LLG.h>

#define pi   3.14159265358

extern int N_Energy_Terms;
pthread_mutex_t **Signal_Synchronize_LLG_Mutex, *M_Synchronize_LLG_Mutex;
pthread_cond_t **Signal_Synchronize_LLG_Cond;

extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
extern pthread_mutex_t OPGL_Mutex;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
extern double Ms_LLG, mu0_LLG, gamma_LLG;
double ***Mx_LLG, ***My_LLG, ***Mz_LLG, *MxV_LLG, *MyV_LLG, *MzV_LLG;
extern OPGL_Parameters_GUI OPGL_GUI;
extern GUI_Parameter GPM;
DiscreatPoint_Run_LLG_Parameter DRLP;


using namespace std;

extern double l1_Strain, l2_Strain, l3_Strain, q_Strain, lo1_Strain, lo2_Strain, lo3_Strain, d_epsilon,
              l1_Stress, l2_Stress, l3_Stress, q_Stress, lo1_Stress, lo2_Stress, lo3_Stress, k_eff;


int Gradient_M_LocalCell(Energy_Variables *EV, double *Gradient_M_Local)
{
    double ***M_Local, *dxy;
    int i, *LocalEnergy, l, s;

    l   = EV->l;     dxy = EV->dxy;
    //CST = MLP->EV->CST;
    M_Local     = EV->M_Local;
    LocalEnergy = EV->LocalEnergy;

    double dx = dxy[0];
    //double dy = dxy[1], dz = dxy[2];
    double MxN, MxS, MyN, MyS, MzN, MzS;
    /*
    double MxO, MxN, MxS, MxW, MxE, MxL, MxU, MxNN, MxSS, MxWW, MxEE, MxLL, MxUU,
           MyO, MyN, MyS, MyW, MyE, MyL, MyU, MyNN, MySS, MyWW, MyEE, MyLL, MyUU,
           MzO, MzN, MzS, MzW, MzE, MzL, MzU, MzNN, MzSS, MzWW, MzEE, MzLL, MzUU;
    */
    enum Local_Point_Postion_L{Null, O, N, S, W, E, L, U, NN, SS, WW, EE, LL, UU};

    Local_Point_Postion_L LPPL;

    for(s=0;s<3;s++)
    {
        Gradient_M_Local[s]=0;
    }

    for(i=0;i<14;++i)
    {
        LPPL=(Local_Point_Postion_L) i;
        switch(LPPL)
        {
            /*
            case O :
            {
                MxO  = *M_Local[1][i];    MyO  = *M_Local[2][i];    MzO  = *M_Local[3][i];
                break;
            }
            */
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
            /*
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
            */
            default:
                break;
        }
    }


    ///////////////////////////////////////////////////////////
    if(l==1)//2D Case
    {
        
        if(LocalEnergy[0]==1)
        {
            //Bloch
            Gradient_M_Local[0] = (MxS - MxN) / (2 * dx);

            Gradient_M_Local[1] = (MyS - MyN) / (2 * dx);

            Gradient_M_Local[2] = (MzS - MzN) / (2 * dx);
        
        }
    }
    else//3D
    {
        printf("LLG with STT is not yet for 3D !");
    }

    return 1;
}


int StochasticLLG_SOT(Modify_LLG_Parameters *MLP)
{
	//double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	double Rx, Ry, Rz;// random variables - thermal fluctuations
	//double ax, ay, az;// temp
	double Ax, Ay, Az;// total matrix
	double detMi;	 // detMi = 1/detM
    int Precession = 1;
    double alpha=1;
    double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);
    double temperature = 0;
    double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
    double CurrentDensity, MpX, MpY, MpZ;
    CurrentDensity = MLP->EV->CurrentDensity;
    MpX = MLP->EV->MpX;
    MpY = MLP->EV->MpY;
    MpZ = MLP->EV->MpZ;
    // CurrentDensity = 0.3;
    // MpX = 1.;
    // MpY = 0.;
    // MpZ = 0.;

    double Mx0, My0, Mz0, Mx1, My1, Mz1, ***M_Local, *M0, M_Module;
    M_Local = MLP->EV->M_Local; M0 = MLP->EV->M0;

    double h_SIB=0.001, rh_SIB;
    h_SIB = MLP->EV->dT;    
    rh_SIB = sqrt(h_SIB);
    double *H_eff_Local;
    H_eff_Local = Make1DArray(3);

    Mx0 = *M_Local[1][1];
    My0 = *M_Local[2][1];
    Mz0 = *M_Local[3][1];

    double MxInitial, MyInitial, MzInitial; 
    MxInitial = Mx0;
    MyInitial = My0;
    MzInitial = Mz0;

    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);
    
    Mx1 = Mx0;  My1 = My0;  Mz1 = Mz0;
    /*
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
*/
    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    //Hx = 0;
    //Hy = 0;
    //Hz = 0;
    //printf("Hx Hy Hz %0.6f\t%0.6f\t%0.6f\n", Hx, Hy, Hz);
    //prediction step of midpoint solver:

	// deterministic terms of Landau–Lifshitz equation:
    // Ax = Alpha_p * Hx + Alpha_d * (My0 * Hz - Mz0 * Hy);
    // Ay = Alpha_p * Hy + Alpha_d * (Mz0 * Hx - Mx0 * Hz);
    // Az = Alpha_p * Hz + Alpha_d * (Mx0 * Hy - My0 * Hx);

    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);


    // Ax = Alpha_p * Hx + Alpha_d * (My0 * Hz - Mz0 * Hy) - CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    // Ay = Alpha_p * Hy + Alpha_d * (Mz0 * Hx - Mx0 * Hz) - CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    // Az = Alpha_p * Hz + Alpha_d * (Mx0 * Hy - My0 * Hx) - CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k1 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k1/6 in global temp array
    Mx1 = Cx / 6.0;
    My1 = Cy / 6.0;
    Mz1 = Cz / 6.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k1/2 will be used on the next step
    Mx0 = Mx0 + Cx * 0.5;
    My0 = My0 + Cy * 0.5;
    Mz0 = Mz0 + Cz * 0.5;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

    //k2
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;

    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k2 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k2/3 in global temp array
    Mx1 = Mx1 + Cx / 3.0;
    My1 = My1 + Cy / 3.0;
    Mz1 = Mz1 + Cz / 3.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k2/2 will be used on the next step
    Mx0 = Mx0 + Cx * 0.5;
    My0 = My0 + Cy * 0.5;
    Mz0 = Mz0 + Cz * 0.5;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

//k3
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;

    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k3 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k2/3 in global temp array
    Mx1 = Mx1 + Cx / 3.0;
    My1 = My1 + Cy / 3.0;
    Mz1 = Mz1 + Cz / 3.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k3 will be used on the next step
    Mx0 = Mx0 + Cx;
    My0 = My0 + Cy;
    Mz0 = Mz0 + Cz;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

    //k4
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;

    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k4 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k4/6 in global temp array
    Mx1 = Mx1 + Cx / 6.0;
    My1 = My1 + Cy / 6.0;
    Mz1 = Mz1 + Cz / 6.0;

    //y_{n+1}=y_n+k1/6+k2/3+k3/3+k4/6 - final step:
    Mx0 = Mx0 + Mx1;
    My0 = My0 + My1;
    Mz0 = Mz0 + Mz1;

    //printf("%0.6f\t%0.6f\t%0.6f\n",Mx1, My1, Mz1);
    //normalize spin
    detMi = 1.0;
    //detMi = 1.0 / sqrt(Mx0 * Mx0  + My0  * My0  + Mz0  * Mz0);


    M0[0] = Mx0 * detMi;
    M0[1] = My0 * detMi;
    M0[2] = Mz0 * detMi;

    double *Gradient;
    Gradient = MLP->Gradient;
    Gradient[0] = M0[0] - MxInitial;
    Gradient[1] = M0[1] - MyInitial;
    Gradient[2] = M0[2] - MzInitial;

    int i, j;
    for (i = 1; i < 14;++i)
    {
        for (j = 1; j < 4;++j)
        {
            //printf("%0.6f\t", *M_Local[j][i]);
        }
        //printf("\n");
    }
    //printf("\n");
    //printf("%0.10f\t%0.10f\t%0.10f\n", M0[0], M0[1], M0[2]);

    //*M_Local[1][1] = M0[0];
    //*M_Local[2][1] = M0[1];
    //*M_Local[3][1] = M0[2];
    //printf("%0.6f\t%0.6f\t%0.6f\n", M0[0], M0[1], M0[2]);
    //printf("\n");
    free(H_eff_Local);
    return 1;
}

int StochasticLLG_STT(Modify_LLG_Parameters *MLP)
{
	//double nx, ny, nz;// components of the unit vector
	double Hx, Hy, Hz;// components of the effective field
	float Cx, Cy, Cz;// spin-torque term
	double Rx, Ry, Rz;// random variables - thermal fluctuations
	//double ax, ay, az;// temp
	double Ax, Ay, Az;// total matrix
	double detMi;	 // detMi = 1/detM
    int Precession = 1;
    double alpha=1, beta=1;
    double Alpha_d = alpha / (1.0 + alpha * alpha * Precession);
	double Alpha_p = 1.0f / (1.0 + alpha * alpha);
    double temperature = 0;
    double D = sqrt(2.0 * alpha / (1.0 + alpha * alpha) * temperature);
    double CurrentDensity, MpX, MpY, MpZ;
    CurrentDensity = MLP->EV->CurrentDensity;


    double Mx0, My0, Mz0, Mx1, My1, Mz1, ***M_Local, *M0, M_Module;
    M_Local = MLP->EV->M_Local; M0 = MLP->EV->M0;

    double h_SIB=0.001, rh_SIB;
    h_SIB = MLP->EV->dT;    
    rh_SIB = sqrt(h_SIB);
    double *H_eff_Local, *Gradient_M;
    H_eff_Local = Make1DArray(3);
    Gradient_M = Make1DArray(3);

    Mx0 = *M_Local[1][1];
    My0 = *M_Local[2][1];
    Mz0 = *M_Local[3][1];

    double MxInitial, MyInitial, MzInitial; 
    MxInitial = Mx0;
    MyInitial = My0;
    MzInitial = Mz0;

    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;


    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);
    
    Mx1 = Mx0;  My1 = My0;  Mz1 = Mz0;
    /*
	//electric DC current vector (VCu) and density (Cu)
	Cx = VCu[0] * Cu;
	Cy = VCu[1] * Cu;
	Cz = VCu[2] * Cu;
*/
    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    Gradient_M_LocalCell(MLP->EV, Gradient_M);
    
    MpX=Gradient_M[0];
    MpY=Gradient_M[1];
    MpZ=Gradient_M[2];
    

    Hx=Hx+beta*CurrentDensity*Gradient_M[0];
    Hy=Hy+beta*CurrentDensity*Gradient_M[1];
    Hz=Hz+beta*CurrentDensity*Gradient_M[2];


    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);


    // Ax = Alpha_p * Hx + Alpha_d * (My0 * Hz - Mz0 * Hy) - CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    // Ay = Alpha_p * Hy + Alpha_d * (Mz0 * Hx - Mx0 * Hz) - CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    // Az = Alpha_p * Hz + Alpha_d * (Mx0 * Hy - My0 * Hx) - CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k1 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k1/6 in global temp array
    Mx1 = Cx / 6.0;
    My1 = Cy / 6.0;
    Mz1 = Cz / 6.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k1/2 will be used on the next step
    Mx0 = Mx0 + Cx * 0.5;
    My0 = My0 + Cy * 0.5;
    Mz0 = Mz0 + Cz * 0.5;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

    //k2
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;
/*
    if (fabs(Mx0) < 1 && fabs(My0) < 1 && fabs(Mz0) < 1)
    {

    }else
    {
        printf("k1\t%0.6f\t%0.6f\t%0.6f\n", Mx0, My0, Mz0);
        printf("k1\tH\t%0.6f\t%0.6f\t%0.6f\n", Hx, Hy, Hz);
        printf("k1\tGradient\t%0.6f\t%0.6f\t%0.6f\n", MpX, MpY, MpZ);
        printf("k1\tA\t%0.6f\t%0.6f\t%0.6f\n", Ax, Ay, Az);
        printf("k1\tC\t%0.6f\t%0.6f\t%0.6f\n", Cx, Cy, Cz);
        system("pause");
    }
*/
    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    Gradient_M_LocalCell(MLP->EV, Gradient_M);
    MpX=Gradient_M[0];
    MpY=Gradient_M[1];
    MpZ=Gradient_M[2];

    Hx=Hx+beta*CurrentDensity*Gradient_M[0];
    Hy=Hy+beta*CurrentDensity*Gradient_M[1];
    Hz=Hz+beta*CurrentDensity*Gradient_M[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k2 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k2/3 in global temp array
    Mx1 = Mx1 + Cx / 3.0;
    My1 = My1 + Cy / 3.0;
    Mz1 = Mz1 + Cz / 3.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k2/2 will be used on the next step
    Mx0 = Mx0 + Cx * 0.5;
    My0 = My0 + Cy * 0.5;
    Mz0 = Mz0 + Cz * 0.5;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

//k3
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;


    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    Gradient_M_LocalCell(MLP->EV, Gradient_M);
    MpX=Gradient_M[0];
    MpY=Gradient_M[1];
    MpZ=Gradient_M[2];

    Hx=Hx+beta*CurrentDensity*Gradient_M[0];
    Hy=Hy+beta*CurrentDensity*Gradient_M[1];
    Hz=Hz+beta*CurrentDensity*Gradient_M[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k3 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k2/3 in global temp array
    Mx1 = Mx1 + Cx / 3.0;
    My1 = My1 + Cy / 3.0;
    Mz1 = Mz1 + Cz / 3.0;
    //printf("%0.6f\n", Mx1);

    //y_n+k3 will be used on the next step
    Mx0 = Mx0 + Cx;
    My0 = My0 + Cy;
    Mz0 = Mz0 + Cz;
    M_Module=Mx0*Mx0 + My0*My0 + Mz0*Mz0;//M^2 
    Alpha_d = alpha / (1.0 + alpha * alpha * Precession * M_Module);
	Alpha_p = 1.0f / (1.0 + alpha * alpha * M_Module);

    //k4
    M0[0] = Mx0;
    M0[1] = My0;
    M0[2] = Mz0;

    Effective_Field_LocalCell_M0(MLP->EV, H_eff_Local);
    Hx = H_eff_Local[0];
    Hy = H_eff_Local[1];
    Hz = H_eff_Local[2];

    Gradient_M_LocalCell(MLP->EV, Gradient_M);
    MpX=Gradient_M[0];
    MpY=Gradient_M[1];
    MpZ=Gradient_M[2];

    Hx=Hx+beta*CurrentDensity*Gradient_M[0];
    Hy=Hy+beta*CurrentDensity*Gradient_M[1];
    Hz=Hz+beta*CurrentDensity*Gradient_M[2];

    // deterministic terms of Landau–Lifshitz equation:
    Ax = Alpha_p * (Hx - alpha * M_Module * CurrentDensity * MpX) + Alpha_d * (My0 * Hz - Mz0 * Hy) + Alpha_p * CurrentDensity * (My0 * MpZ - Mz0 * MpY);
    Ay = Alpha_p * (Hy - alpha * M_Module * CurrentDensity * MpY) + Alpha_d * (Mz0 * Hx - Mx0 * Hz) + Alpha_p * CurrentDensity * (Mz0 * MpX - Mx0 * MpZ);
    Az = Alpha_p * (Hz - alpha * M_Module * CurrentDensity * MpZ) + Alpha_d * (Mx0 * Hy - My0 * Hx) + Alpha_p * CurrentDensity * (Mx0 * MpY - My0 * MpX);
/*
    // Spin-torque term
    Ax = Ax + 0.5f * h * (-alpha * Cx + (ny * Cz - nz * Cy)); //pay attention to the signe and factors
    Ay = Ay + 0.5f * h * (-alpha * Cy + (nz * Cx - nx * Cz));
    Az = Az + 0.5f * h * (-alpha * Cz + (nx * Cy - ny * Cx));
*/

    Rx = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Ry = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    Rz = 2.0 * (0.5 - rand() / (float)RAND_MAX);
    // stochastic terms of Landau–Lifshitz equation:
    Ax = Ax + rh_SIB * D * (Alpha_p * Rx + Alpha_d * (My0 * Rz - Mz0 * Ry));
    Ay = Ay + rh_SIB * D * (Alpha_p * Ry + Alpha_d * (Mz0 * Rx - Mx0 * Rz));
    Az = Az + rh_SIB * D * (Alpha_p * Rz + Alpha_d * (Mx0 * Ry - My0 * Rx));

    //k4 (Cx,Cy,Cz here are used as temp variables)
    Cx = -h_SIB * (My0 * Az - Mz0 * Ay);
    Cy = -h_SIB * (Mz0 * Ax - Mx0 * Az);
    Cz = -h_SIB * (Mx0 * Ay - My0 * Ax);

    //save k4/6 in global temp array
    Mx1 = Mx1 + Cx / 6.0;
    My1 = My1 + Cy / 6.0;
    Mz1 = Mz1 + Cz / 6.0;

    //y_{n+1}=y_n+k1/6+k2/3+k3/3+k4/6 - final step:
    Mx0 = Mx0 + Mx1;
    My0 = My0 + My1;
    Mz0 = Mz0 + Mz1;

    //printf("%0.6f\t%0.6f\t%0.6f\n",Mx1, My1, Mz1);
    //normalize spin
    detMi = 1.0;
    //detMi = 1.0 / sqrt(Mx0 * Mx0  + My0  * My0  + Mz0  * Mz0);


    M0[0] = Mx0 * detMi;
    M0[1] = My0 * detMi;
    M0[2] = Mz0 * detMi;

    double *Gradient;
    Gradient = MLP->Gradient;
    Gradient[0] = M0[0] - MxInitial;
    Gradient[1] = M0[1] - MyInitial;
    Gradient[2] = M0[2] - MzInitial;

    int i, j;
    for (i = 1; i < 14;++i)
    {
        for (j = 1; j < 4;++j)
        {
            //printf("%0.6f\t", *M_Local[j][i]);
        }
        //printf("\n");
    }
    //printf("\n");
    //printf("%0.10f\t%0.10f\t%0.10f\n", M0[0], M0[1], M0[2]);

    //*M_Local[1][1] = M0[0];
    //*M_Local[2][1] = M0[1];
    //*M_Local[3][1] = M0[2];
    //printf("%0.6f\t%0.6f\t%0.6f\n", M0[0], M0[1], M0[2]);
    free(H_eff_Local);
    return 1;
}


int Modify_Pthread_LLG(Modify_LLG_Pthreads *MLPt)
{
    int Pthread_Index, NCycles, *Signal, MagneticFieldMode, MagneticFieldPeriod;
    int Signal_Synchronize;
    int CalculatePointsThreads_Lower, CalculatePointsThreads_Upper, CalculatePointsThreads;
    int i, N, Lower, Signal_RunAndStop;
    double h;
    double DerivativeM=0.;
    Calculate_Mode_GUI CMG;
    CMG = (Calculate_Mode_GUI)MLPt->CalculateMode;

    Pthread_Index  = MLPt->Pthread_Index;
    //Pthread_Number = MLPt->Pthread_Number;
    NCycles        = MLPt->NCycles;
    Signal         = MLPt->Signal;

    CalculatePointsThreads_Lower = MLPt->CalculatePointsThreads_Lower;
    CalculatePointsThreads_Upper = MLPt->CalculatePointsThreads_Upper;
    CalculatePointsThreads       = CalculatePointsThreads_Upper-CalculatePointsThreads_Lower;

    if(CalculatePointsThreads_Lower==0)
    {
        Lower=1;
    }else
    {
        Lower=0;
    }

    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = SignalPthread_Run;
    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);//1 is Pthread State.
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    MagneticFieldMode = MLPt->MagneticFieldMode;
    MagneticFieldPeriod = MLPt->MagneticFieldPeriod;
    h = MLPt->h;
    ////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("Pthread Run\t%d\n",Pthread_Index);
    for(N=0;N<NCycles;++N)
    {
        //printf("%d\t%d\n",Pthread_Index,N);
        pthread_mutex_lock(&M_Synchronize_LLG_Mutex[Pthread_Index]);
        MLPt->PTPS->Energy = 0.;  MLPt->PTPS->SKN = 0.;
        if (MagneticFieldMode==1)
        {
            h = h * sin(2 * pi * N / MagneticFieldPeriod);
        }
        for(i=Lower;i<CalculatePointsThreads;++i)
        {
            //printf("%d\t%d\n",Pthread_Index,i);
            MLPt->MLP[i].EV->h = h;
            MLPt->MLP[i].W_BeforeModify = Energy_Single(MLPt->MLP[i].EV);
            MLPt->PTPS->Energy = MLPt->PTPS->Energy + Energy_Single(MLPt->MLP[i].EV);
            MLPt->PTPS->SKN    = MLPt->PTPS->SKN + Skyrmion_Number_Single(MLPt->MLP[i].EV);
            switch (CMG)
            {
                case CMG_LLG_SOT:
                {
                    StochasticLLG_SOT(&(MLPt->MLP[i]));
                }
                break;
            
                case CMG_LLG_STT:
                {
                    StochasticLLG_STT(&(MLPt->MLP[i]));
                }
                break;

                default:
                    break;
            }
            DerivativeM = DerivativeM + sqrt(pow(MLPt->MLP[i].Gradient[0], 2)+pow(MLPt->MLP[i].Gradient[1], 2)+pow(MLPt->MLP[i].Gradient[2], 2));
            if(MLPt->Continuity_Modify_Mode==1)
            {
                Continue_Modify(MLPt->MLP[i].EV);
            }
        }
        pthread_mutex_unlock(&M_Synchronize_LLG_Mutex[Pthread_Index]);

        if (Signal_Synchronize != SignalSyn_Stop)
        {
            ////////////////////////// Prepare to Synchronize ///////////////////////////
            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]); //1 is Pthread State.
            Signal[1] = StatePthread_PrepareToSynchronize;
            //printf("%d\tStatePthread_PrepareToSynchronize\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);//1 is Pthread State.

            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_Synchronize && Signal[3]!=SignalSyn_NotActive && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_LLG_Cond[Pthread_Index][3], &Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);
            }
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            //printf("%d\tSignalPthread_Stop\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }

        if (Signal_Synchronize != SignalSyn_Stop)
        {
            ////////////////////////// Synchronize ///////////////////////////
            pthread_mutex_lock(&M_Synchronize_LLG_Mutex[Pthread_Index]);
            //MLPt->PTPS->Energy = 0.;  MLPt->PTPS->SKN = 0.;
            for (i = Lower; i < CalculatePointsThreads; ++i)
            {
                *MLPt->MLP[i].EV->M_Local[1][1] = MLPt->MLP[i].EV->M0[0];
                *MLPt->MLP[i].EV->M_Local[2][1] = MLPt->MLP[i].EV->M0[1];
                *MLPt->MLP[i].EV->M_Local[3][1] = MLPt->MLP[i].EV->M0[2];
            }
            pthread_mutex_unlock(&M_Synchronize_LLG_Mutex[Pthread_Index]);

            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_PrepareToRun;
            //printf("%d\tStatePthread_PrepareToRun\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);//1 is Pthread State.

            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);  //3 is Synchronize Signal.
            while(Signal[3]!=SignalSyn_PrepareToSynchronize && Signal[3]!=SignalSyn_NotActive && Signal[3]!=SignalSyn_Stop)//3 is Synchronize Signal.
            {
                pthread_cond_wait(&Signal_Synchronize_LLG_Cond[Pthread_Index][3], &Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);
            }
            Signal_Synchronize = Signal[3];
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][3]);//3 is Synchronize Signal.

            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);  //1 is Pthread State.
            Signal[1] = StatePthread_Run;
            //printf("%d\tStatePthread_Run\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);//1 is Pthread State.
        }else
        {
            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][2]);  //2 is Pthread Signal.
            Signal[2] = SignalPthread_Stop;
            //printf("%d\tSignalPthread_Stop\n", Pthread_Index);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][2]);//2 is Pthread Signal.
            break;
        }
        /////////////////////////// Suspend thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        while (GPM.Thread_State[Pthread_Index + 1][1] == Signal_GUI_Break && GPM.Thread_State[Pthread_Index + 1][0] != Signal_GUI_Stop) //Signal_ContinueAndBreak
        {
            pthread_cond_wait(&GPM.Pthread_GUI_Cond[Pthread_Index + 1], &GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        }
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        
        /////////////////////////// Stop thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        Signal_RunAndStop = GPM.Thread_State[Pthread_Index + 1][0];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Pthread_Index + 1]);
        if (Signal_RunAndStop == Signal_GUI_Stop)
        {
            break;
        }
    }

    MLPt->DerivativeM = DerivativeM / NCycles;

    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);  //1 is Pthread State.
    Signal[1] = StatePthread_Stop;
    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[Pthread_Index][1]);//1 is Pthread State.

    printf("fin\t%d\n",Pthread_Index);
    return 1;

}

void* Modify_Pthread_LLG_Transfer(void *arg)
{
    Modify_LLG_Pthreads *MLPt;
    MLPt=(Modify_LLG_Pthreads*)arg;
    Modify_Pthread_LLG(MLPt);
    return NULL;
}


int Modify_Pthread_LLG_Synchronize(Modify_LLG_Synchronize *MLS)
{
    int i, j, k, v, p, Point_Index, Boundary_Type, Synchronize_Times, MinPositionX, MinPositionY;;
    int ***Local_Points_1D, ***Position_3DTo1D, m, n, l, Virtual_Points;
    int Pthread_Number, ***Signal, Pthread_cond_State, Signal_RunAndStop;
    int **Bound_CalculatePoint_Threads, N_Total;
    double ****CST, Energy_Average, SKN;
    double MzMin, RealTime, dT;

    clock_t start_time_P, end_time_P;
    double duration;

    m=MLS->m;   n=MLS->n;   l=MLS->l;   Virtual_Points=MLS->Virtual_Points;
    Pthread_Number               = MLS->Pthread_Number;
    Bound_CalculatePoint_Threads = MLS->Bound_CalculatePoint_Threads;
    Local_Points_1D              = MLS->Local_Points_1D;
    Position_3DTo1D              = MLS->Position_3DTo1D;
    Signal                       = MLS->Signal;
    CST                          = MLS->CST;
    dT                           = MLS->IPN.dT;
    N_Total                      = MLS->Inner_Points + MLS->Boundary_Points;
///////////////////////////////////////////////////////////////////////////
    int zc;
    if(l==0)
    {
        zc=0;
    }else
    {
        zc=(l-1)/2;
    }

    FILE *fp_Mx,*fp_My, *fp_Mz, *fp_MxV,*fp_MyV, *fp_MzV, *fp_Process, *fp_Trace;
    char filenameMx[500], filenameMy[500], filenameMz[500],
         filenameMxV[500],filenameMyV[500],filenameMzV[500],filenameP[500],filenameTr[500];

    sprintf(filenameMxV,"Version %d\\Process\\MxV.dat",MLS->Version);
    sprintf(filenameMyV,"Version %d\\Process\\MyV.dat",MLS->Version);
    sprintf(filenameMzV,"Version %d\\Process\\MzV.dat",MLS->Version);
    sprintf(filenameP, "Version %d\\Process.dat",MLS->Version);
    sprintf(filenameTr, "Version %d\\Trace.dat",MLS->Version);
    fp_MxV     = fopen(filenameMxV,"w+");
    fp_MyV     = fopen(filenameMyV,"w+");
    fp_MzV     = fopen(filenameMzV,"w+");
    fp_Process = fopen(filenameP, "w+");
    fp_Trace   = fopen(filenameTr, "w+");
    fprintf(fp_Process, "Times\t\t W\t\t SKN\t\t Time\n");
    fprintf(fp_Trace,"Times\t\t X\t\t Y\n");
    ////////////////////////////////////////////////////////////////////////////
    int Pthread_State, Synchronize_State;
    Synchronize_Times=0;

    do
    {
        start_time_P = clock();
        ///////////////////////// Make sure all child threads are waiting ///////////////////////////////
        do
        {
            Synchronize_State = StateSyn_Synchronize;
            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
                Pthread_State=*Signal[p+1][1];
                pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
                if (Pthread_State == StatePthread_Stop)
                {
                    Synchronize_State = StateSyn_Stop;
                    break;
                }else if(Pthread_State != StatePthread_PrepareToSynchronize)
                {
                    Synchronize_State = StateSyn_PrepareToSynchronize;
                }
            }
        } while (Synchronize_State == StateSyn_PrepareToSynchronize);

        ////////////////////// Calculate Energy //////////////////////////////////////////////////////
        Energy_Average = 0.;    SKN = 0.;
        for (p = 0; p < Pthread_Number; ++p)
        {
            Energy_Average   = Energy_Average + MLS->PTPS[p].Energy / N_Total;
            SKN = SKN + MLS->PTPS[p].SKN;
        }

        ///////////////////////// Stop child threads ///////////////////////////////
        if (Synchronize_State == StateSyn_Stop)
        {
            for(p=0;p<Pthread_Number;++p)
            {
                Sleep(100);
                pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
                pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State=*Signal[p+1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun && Pthread_State != StatePthread_PrepareToSynchronize)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_LLG_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }
        }
        /////////////////////// Synchronize /////////////////////////////////////////
        if (Synchronize_State == StateSyn_Synchronize)
        {
            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_Synchronize;
                pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToSynchronize)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]); //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_LLG_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }

            do
            {
                Synchronize_State = StateSyn_PrepareToRun;
                for(p=0;p<Pthread_Number;++p)
                {
                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State=*Signal[p+1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun)
                    {
                        Synchronize_State = StateSyn_Synchronize;
                    }
                }
            } while (Synchronize_State == StateSyn_Synchronize);

            for(p=0;p<Pthread_Number;++p)
            {
                pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
                *Signal[p + 1][3] = SignalSyn_PrepareToSynchronize;
                pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                do
                {
                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
                    Pthread_State = *Signal[p + 1][1];
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
                    if (Pthread_State != StatePthread_PrepareToRun)
                    {
                        break;
                    }

                    pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
                    Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_LLG_Cond[p][3]);
                    pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
                } while (Pthread_cond_State != 0);
            }
        }

////////////////////////////////////////////////////////////////////////////////////////

        pthread_mutex_lock(&OPGL_Mutex);
        OPGL_GUI.Energy = Energy_Average;
        OPGL_GUI.SKN = SKN;
        pthread_mutex_unlock(&OPGL_Mutex);
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(k=0;k<l;++k)
                {
                    Point_Index = Position_3DTo1D[i][j][k];
                    Boundary_Type = Local_Points_1D[Point_Index][0][0];
                    if(Boundary_Type!=0)
                    {
                        pthread_mutex_lock(&OPGL_Mutex);
                        Mx_OPGL[i][j][k] = Mx_LLG[i][j][k];
                        My_OPGL[i][j][k] = My_LLG[i][j][k];
                        Mz_OPGL[i][j][k] = Mz_LLG[i][j][k];
                        pthread_mutex_unlock(&OPGL_Mutex);
                    }else
                    {
                        pthread_mutex_lock(&OPGL_Mutex);
                        Mx_OPGL[i][j][k] = 0.;
                        My_OPGL[i][j][k] = 0.;
                        Mz_OPGL[i][j][k] = 0.;
                        pthread_mutex_unlock(&OPGL_Mutex);
                    }
                    
                }
            }
        }

        for(v=0;v<Virtual_Points;++v)
        {
            pthread_mutex_lock(&OPGL_Mutex);
            MxV_OPGL[v] = MxV_LLG[v];
            MyV_OPGL[v] = MyV_LLG[v];
            MzV_OPGL[v] = MzV_LLG[v];
            pthread_mutex_unlock(&OPGL_Mutex);
        }

        ////////////////////// Calculate Energy //////////////////////////////////////////////////////
        /*
        Energy_Average = 0.;    SKN = 0.;
        for (p = 0; p < Pthread_Number; ++p)
        {
            Energy_Average   = Energy_Average + MLS->PTPS[p].Energy / N_Total;
            SKN = SKN + MLS->PTPS[p].SKN;
        }
        */
        //Energy_Distribution(&EVD);
        if(Synchronize_Times%100==0)
        {
            MzMin = Mz_LLG[0][0][0];
            for(k=0;k<l;k++)
            {
                if(Synchronize_Times%2000==0)
                {
                    sprintf(filenameMx,"Version %d\\Process\\z=%d\\Mx%d.dat",MLS->Version,k-zc,Synchronize_Times);
                    sprintf(filenameMy,"Version %d\\Process\\z=%d\\My%d.dat",MLS->Version,k-zc,Synchronize_Times);
                    sprintf(filenameMz,"Version %d\\Process\\z=%d\\Mz%d.dat",MLS->Version,k-zc,Synchronize_Times);
                    fp_Mx=fopen(filenameMx,"w+");
                    fp_My=fopen(filenameMy,"w+");
                    fp_Mz=fopen(filenameMz,"w+");
                    for(i=0;i<m;i++)
                    {
                        for(j=0;j<n;j++)
                        {
                            fprintf(fp_Mx,"%0.8f\t",Mx_LLG[i][j][k]);
                            fprintf(fp_My,"%0.8f\t",My_LLG[i][j][k]);
                            fprintf(fp_Mz,"%0.8f\t",Mz_LLG[i][j][k]);
                        }
                        fprintf(fp_Mx,"\n");
                        fprintf(fp_My,"\n");
                        fprintf(fp_Mz,"\n");
                    }
                    fclose(fp_Mx);
                    fclose(fp_My);
                    fclose(fp_Mz);
                }
                sprintf(filenameMx,"Version %d\\Process\\z=%d\\Mx.dat",MLS->Version,k-zc);
                sprintf(filenameMy,"Version %d\\Process\\z=%d\\My.dat",MLS->Version,k-zc);
                sprintf(filenameMz,"Version %d\\Process\\z=%d\\Mz.dat",MLS->Version,k-zc);
                fp_Mx=fopen(filenameMx,"w+");
                fp_My=fopen(filenameMy,"w+");
                fp_Mz=fopen(filenameMz,"w+");

                for(i=0;i<m;i++)
                {
                    for(j=0;j<n;j++)
                    {
                        fprintf(fp_Mx,"%0.8f\t",Mx_LLG[i][j][k]);
                        fprintf(fp_My,"%0.8f\t",My_LLG[i][j][k]);
                        fprintf(fp_Mz,"%0.8f\t",Mz_LLG[i][j][k]);
                        if (Mz_LLG[i][j][k]<MzMin)
                        {
                            MzMin = Mz_LLG[i][j][k];
                            MinPositionX = i;
                            MinPositionY = j;
                        }
                    }
                    fprintf(fp_Mx,"\n");
                    fprintf(fp_My,"\n");
                    fprintf(fp_Mz,"\n");
                }

                fclose(fp_Mx);
                fclose(fp_My);
                fclose(fp_Mz);
            }

            for(v=0;v<Virtual_Points;++v)
            {
                fprintf(fp_MxV,"%0.8f\t",MxV_LLG[v]);
                fprintf(fp_MyV,"%0.8f\t",MyV_LLG[v]);
                fprintf(fp_MzV,"%0.8f\t",MzV_LLG[v]);
            }

            fclose(fp_MxV);
            fclose(fp_MyV);
            fclose(fp_MzV);
        }
        end_time_P = clock();
        duration = (double)(end_time_P-start_time_P)/CLOCKS_PER_SEC;
        RealTime = Synchronize_Times * dT / (gamma_LLG * mu0_LLG * Ms_LLG) * 1.e12;
        //printf("%lf\n", RealTime);
        fprintf(fp_Trace,"%0.8f\t\t %d\t\t %d\n",
                RealTime,MinPositionX,MinPositionY);
        fprintf(fp_Process, "%d\t\t %0.8f\t %0.8f\t %0.8f\n",
                Synchronize_Times, Energy_Average, SKN, duration);
        printf("%d\t\t %0.8f\t %0.8f\t %0.8f\n",
               Synchronize_Times, Energy_Average, SKN, duration);

        Synchronize_Times=Synchronize_Times+1;
        /////////////////////////// Suspend thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
        while (GPM.Thread_State[0][1] == Signal_GUI_Break)//Signal_ContinueAndBreak
        {
            pthread_cond_wait(&GPM.Pthread_GUI_Cond[0], &GPM.Pthread_GUI_Mutex[0]);
        }
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);

        /////////////////////////// Stop thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[0]);
        Signal_RunAndStop = GPM.Thread_State[0][0];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[0]);
        //Stop Signal
        if (Signal_RunAndStop == Signal_GUI_Stop)
        {
            break;
        }
    }while(Synchronize_State != StateSyn_Stop);

    for(p=0;p<Pthread_Number;++p)
    {
        pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]);  //3 is Synchronize Signal.
        *Signal[p + 1][3] = SignalSyn_Stop;                        //Stop all threads
        pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
        do
        {
            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][1]);  //1 is Pthread State.
            Pthread_State = *Signal[p + 1][1];
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][1]);//1 is Pthread State.
            if (Pthread_State == StatePthread_Stop)
            {
                break;
            }

            pthread_mutex_lock(&Signal_Synchronize_LLG_Mutex[p][3]); //3 is Synchronize Signal.
            Pthread_cond_State = pthread_cond_signal(&Signal_Synchronize_LLG_Cond[p][3]);
            pthread_mutex_unlock(&Signal_Synchronize_LLG_Mutex[p][3]);//3 is Synchronize Signal.
        } while (Pthread_cond_State != 0);
        
    }
    fclose(fp_Process);

//////////////////////// Stop OPGL //////////////////////////////////////////////////////////
    pthread_mutex_lock(&OPGL_Mutex);
    OPGL_GUI.State_Syn = StateSyn_Stop;
    pthread_mutex_unlock(&OPGL_Mutex);


    printf("Syn fin\n");
    return 1;
}

void* Modify_Pthread_LLG_Synchronize_Transfer(void *arg)
{
    Modify_LLG_Synchronize *MLS;
    MLS=(Modify_LLG_Synchronize*)arg;
    Modify_Pthread_LLG_Synchronize(MLS);
    return NULL;
}

int OPGL_LLG_Pthread(OPGL_LLG_Pthread_Parameters *OLPP)
{
    OPGL_Parameters OPM;
    OPM.IPN = OLPP->IPN;
    OPM.Mx_OPGL = OLPP->Mx_OPGL;
    OPM.My_OPGL = OLPP->My_OPGL;
    OPM.Mz_OPGL = OLPP->Mz_OPGL;
    OPM.MxV_OPGL = OLPP->MxV_OPGL;
    OPM.MyV_OPGL = OLPP->MyV_OPGL;
    OPM.MzV_OPGL = OLPP->MzV_OPGL;
    OPM.CST = OLPP->CST;
    OPGL_NMD(&OPM);
    return 1;
}

void* OPGL_LLG_Pthread_Transfer(void *arg)
{
    OPGL_LLG_Pthread_Parameters *OLPP;
    OLPP=(OPGL_LLG_Pthread_Parameters*)arg;
    OPGL_LLG_Pthread(OLPP);
    return NULL;
}

int DiscreatPoint_Run_LLG_Part_Run(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, m, n, l, Version, Threads_Number, Calculate_Points, CalculatePointsThreads, MagneticFieldMode, MagneticFieldPeriod;
    double h, t, *dxy, dT=0.005, CurrentDensity, MpX, MpY, MpZ;
    double ***M_Local, ****CST;
    int i_G, j_G, k_G, v, r, s, i, j, k, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    int ****Local_Points_1D_LocalEnergy, **Neighbour_Points_3Dto1D_Array, Point_Index;
    int Boundary_Type;

    clock_t start_time;
    start_time = clock();
    DRLP.start_time = start_time;
    Parameter_Transfer_Pthread_Synchronize *PTPS;
////////////////////////// Assignment ///////////////////
    xn  = IPN->xn;  yn  = IPN->yn;
    zn  = IPN->zn;  dxy = IPN->dxy;
    h   = IPN->h;   t   = IPN->t;
    Version    =IPN->Version;
    dT   = IPN->dT;
    CurrentDensity=IPN->CurrentDensity;
    MpX = IPN->MpX;
    MpY = IPN->MpY;
    MpZ = IPN->MpZ;
    MagneticFieldMode = IPN->MagneticFieldMode;
    MagneticFieldPeriod = IPN->MagneticFieldPeriod;
    //NCycles    =IPN->NCycles;
    Threads_Number   = IPN->Threads_Number;
    PTPS = (Parameter_Transfer_Pthread_Synchronize *)malloc(sizeof(Parameter_Transfer_Pthread_Synchronize) * (Threads_Number));

    m = 2*xn+3; n = 2*yn+3;
    if(zn==0)
    {
        l=1;
    }else
    {
        l = 2*zn+3;
    }

    Mx_LLG = Make3DArray(m, n, l);
    My_LLG = Make3DArray(m, n, l);
    Mz_LLG = Make3DArray(m, n, l);

    ///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape *BS;
    BS = &DRLP.BS;
    Calculate_Points_Parameters *CPP;
    CPP = &DRLP.CPP;

    BS->xn=xn;  BS->yn=yn;
    BS->zn=zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    Calculate_Points=BS->Calculate_Points;

    CPP->Position_3DTo1D  = Make3DArrayinteger(m,n,l);
    CPP->Local_Points_1D  = Make3DArrayinteger(BS->Calculate_Points,7,7);
    CPP->Virtual_Position = Make1DArrayinteger(BS->Virtual_Points);
    CPP->Local_Points_1D_LocalEnergy       = Make4DArrayinteger(BS->Calculate_Points,7,7,7);
    CPP->Neighbour_Points_3Dto1D_Array     = Make2DArrayinteger(14,2);
    CPP->Sort_Calculate_Points             = Make1DArrayinteger(BS->Calculate_Points);
    CalculPointsParameters_Initial(CPP,BS);
    Neighbour_Points_3Dto1D(CPP);
    Local_Points_1D_LocalEnergy   = CPP->Local_Points_1D_LocalEnergy;
    Neighbour_Points_3Dto1D_Array = CPP->Neighbour_Points_3Dto1D_Array;
    ReSort_Calcualte_Points(CPP,BS);
    //ReSort_Calcualte_Points_ByThreads(CPP,BS,Threads_Number);
    //ReSort_Calcualte_Points_ByOrders(CPP,BS);
    //ReSort_Calcualte_Points_BySymmetry_Xaxis(CPP,BS);
    //ReSort_Calcualte_Points_BySymmetry_Yaxis(CPP,BS);
    //ReSort_Calcualte_Points_ByRow(CPP,BS);
    ////////////////////////////////////////////////////////////////////////////
    MxV_LLG = Make1DArray(BS->Virtual_Points);
    MyV_LLG = Make1DArray(BS->Virtual_Points);
    MzV_LLG = Make1DArray(BS->Virtual_Points);
    ////////////////////////////////////////////////////////////////////////////
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                Mx_LLG[i][j][k] = Mx_Global[i][j][k];
                My_LLG[i][j][k] = My_Global[i][j][k];
                Mz_LLG[i][j][k] = Mz_Global[i][j][k];
            }
        }
    }

    for (v = 0; v < BS->Virtual_Points; ++v)
    {
        MxV_LLG[v] = MxV_Global[v];
        MyV_LLG[v] = MyV_Global[v];
        MzV_LLG[v] = MzV_Global[v];
    }

    Energy_Variables **EV;
    Modify_LLG_Pthreads *MLPt;

    MLPt=(Modify_LLG_Pthreads*)malloc(sizeof(Modify_LLG_Pthreads)*(Threads_Number));

    M_Local       = Make3DArray(Threads_Number,4,14);
    CST           = Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);
    pthread_mutex_init(&OPGL_Mutex, 0);
    Signal_Synchronize_LLG_Mutex = (pthread_mutex_t**)malloc(sizeof(pthread_mutex_t*)*(Threads_Number));
    Signal_Synchronize_LLG_Cond  = (pthread_cond_t**)malloc(sizeof(pthread_cond_t*)*(Threads_Number));
    M_Synchronize_LLG_Mutex      = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(Threads_Number));

    int ***Signal;
    Signal=Make2DArrayintegerPointer(Threads_Number+1,4);
    EV=(Energy_Variables**)malloc(sizeof(Energy_Variables*)*(Threads_Number));
    int p, q, o, I, J;
///////////////////////// Pthreads Parameters Initial /////////////////////////
    for(p=0;p<Threads_Number;++p)
    {
        //MLPt[p].CST           =Make4DArray(7,m,n,l);
        //Elastic_Field_Initial(MLPt[p].CST,IPN);
        MLPt[p].Pthread_Index  = p;
        MLPt[p].CalculateMode  = IPN->CalculateMode;
        MLPt[p].NCycles        = IPN->NCycles;
        MLPt[p].Pthread_Number = Threads_Number;

        MLPt[p].Sort_Calculate_Points = Make1DArrayinteger(BS->Calculate_Points);
        for(i=0;i<BS->Calculate_Points;++i)
        {
            MLPt[p].Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
        }
        MLPt[p].Local_Points_1D = CPP->Local_Points_1D;
        MLPt[p].Position_3DTo1D = CPP->Position_3DTo1D;
        MLPt[p].h   = h;    MLPt[p].t   = t;
        MLPt[p].m   = m;    MLPt[p].n   = n;
        MLPt[p].l   = l;
        MLPt[p].MagneticFieldPeriod = MagneticFieldPeriod;
        MLPt[p].MagneticFieldMode = MagneticFieldMode;
        MLPt[p].dxy = dxy;  MLPt[p].CST = CST;
        MLPt[p].Continuity_Modify_Mode = IPN->Continuity_Modify_Mode;
        MLPt[p].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
        


        CalculatePointsThreads_Lower=floor(Calculate_Points/Threads_Number)*p;
        if(p==Threads_Number-1)
        {
            CalculatePointsThreads       = Calculate_Points-floor(Calculate_Points/Threads_Number)*p;
            CalculatePointsThreads_Upper = Calculate_Points;
        }else
        {
            CalculatePointsThreads       = floor(Calculate_Points/Threads_Number);
            CalculatePointsThreads_Upper = floor(Calculate_Points/Threads_Number)*(p+1);
        }
        MLPt[p].CalculatePointsThreads_Lower=CalculatePointsThreads_Lower;
        MLPt[p].CalculatePointsThreads_Upper=CalculatePointsThreads_Upper;

        MLPt[p].PTPS = &PTPS[p];

        /////////////////////////// Synchronize Parameter Initial ////////////////////////////////
        MLPt[p].Signal          = Make1DArrayinteger(4);
        MLPt[p].Signal[1]       = StatePthread_Run;//0 is Pthread State.
        MLPt[p].Signal[2]       = SignalPthread_Run;//1 is Pthread Signal.
        MLPt[p].Signal[3]       = SignalSyn_PrepareToSynchronize;//2 is Synchronize Signal.
        //MLPt[p].Signal[3]       = SignalSyn_NotActive;//2 is Synchronize Signal.
        Signal[p+1][1]          = &MLPt[p].Signal[1];
        Signal[p+1][2]          = &MLPt[p].Signal[2];
        Signal[p+1][3]          = &MLPt[p].Signal[3];

        Signal_Synchronize_LLG_Mutex[p] = (pthread_mutex_t*)malloc(sizeof(pthread_mutex_t)*(4));
        Signal_Synchronize_LLG_Cond[p]  = (pthread_cond_t*)malloc(sizeof(pthread_cond_t)*(4));
        pthread_mutex_init(&Signal_Synchronize_LLG_Mutex[p][1],0);//1 is Pthread State.
        pthread_mutex_init(&Signal_Synchronize_LLG_Mutex[p][2],0);//2 is Pthread Signal.
        pthread_mutex_init(&Signal_Synchronize_LLG_Mutex[p][3],0);//3 is Synchronize Signal.
        pthread_cond_init(&Signal_Synchronize_LLG_Cond[p][1],NULL);//1 is Pthread State.
        pthread_cond_init(&Signal_Synchronize_LLG_Cond[p][2],NULL);//2 is Pthread Signal.
        pthread_cond_init(&Signal_Synchronize_LLG_Cond[p][3],NULL);//3 is Synchronize Signal.
        pthread_mutex_init(&M_Synchronize_LLG_Mutex[p],0);

        /////////////////////////// Discrete Point Parameter Initial //////////////////////////
        MLPt[p].MLP    = (Modify_LLG_Parameters*)malloc(sizeof(Modify_LLG_Parameters)*(CalculatePointsThreads));
        EV[p]          = (Energy_Variables*)malloc(sizeof(Energy_Variables)*(CalculatePointsThreads));

        for(r=0;r<14;++r)
        {
            M_Local[p][1][r]=0.;//Mx
            M_Local[p][2][r]=0.;//My
            M_Local[p][3][r]=0.;//Mz
        }

        for(q=0;q<CalculatePointsThreads;++q)
        {
            MLPt[p].MLP[q].EV=&EV[p][q];
            MLPt[p].MLP[q].Gradient=Make1DArray(3);
            MLPt[p].MLP[q].Gradient[0] = 0; MLPt[p].MLP[q].Gradient[1] = 0; MLPt[p].MLP[q].Gradient[2] = 0;
            EV[p][q].dxy=Make1DArray(3);
            EV[p][q].dxy[0] = dxy[0];    EV[p][q].dxy[1] = dxy[1];
            EV[p][q].dxy[2] = dxy[2];    EV[p][q].i      = 3;
            EV[p][q].j      = 3;         EV[p][q].k      = 3;
            EV[p][q].m      = m;         EV[p][q].n      = n;
            EV[p][q].l      = l;
            EV[p][q].h      = h;         EV[p][q].t      = t;
            EV[p][q].dT     = dT;
            EV[p][q].M0          = Make1DArray(3);
            EV[p][q].LocalEnergy = Make1DArrayinteger(7);
            EV[p][q].CST         = Make2DArray(7,7);
            EV[p][q].Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
            EV[p][q].dT = dT;
            EV[p][q].CurrentDensity = CurrentDensity;
            EV[p][q].MpX = MpX;
            EV[p][q].MpY = MpY;
            EV[p][q].MpZ = MpZ;

            EV[p][q].M_Local     = Make2DArrayPointer(4,15);
            EV[p][q].Energy_Coefficient = Make1DArray(N_Energy_Terms);

            for(r=0;r<5;++r)
            {
                EV[p][q].M_Local[1][r] = &M_Local[p][1][r];
                EV[p][q].M_Local[2][r] = &M_Local[p][2][r];
                EV[p][q].M_Local[3][r] = &M_Local[p][3][r];
            }

            for (r = 0; r < N_Energy_Terms; ++r)
            {
                EV[p][q].Energy_Coefficient[r] = IPN->Energy_Coefficient[r];
            }

            Point_Index=q+CalculatePointsThreads_Lower;
            Point_Index=CPP->Sort_Calculate_Points[Point_Index];
            //printf("%d\n",Point_Index);
/////////////////////////////////////////////////////////////////////////////////////////////////
            for(I=0;I<7;++I)
            {
                EV[p][q].LocalEnergy[I] = CPP->Local_Points_1D[Point_Index][I][5];

                Boundary_Type           = CPP->Local_Points_1D[Point_Index][I][0];
                i_G                     = CPP->Local_Points_1D[Point_Index][I][1];
                j_G                     = CPP->Local_Points_1D[Point_Index][I][2];
                k_G                     = CPP->Local_Points_1D[Point_Index][I][3];
                //printf("%d\t%d\t%d\t%d\t%d\n",Boundary_Type,i_G,j_G,k_G,v);
                for(J=0;J<7;++J)
                {
                    if(Boundary_Type!=0)
                    {
                        EV[p][q].CST[J][I]  = CST[J][i_G][j_G][k_G];
                    }else
                    {
                        if(CPP->Local_Points_1D[Point_Index][I][5]==1)
                        {
                            EV[p][q].CST[J][I] = CST[J][i_G][j_G][k_G];
                        }else
                        {
                            EV[p][q].CST[J][I] = 0.;
                        }
                    }
                }
            }

            for(o=1;o<14;++o)
            {
                r = Neighbour_Points_3Dto1D_Array[o][0];
                s = Neighbour_Points_3Dto1D_Array[o][1];

                Boundary_Type = Local_Points_1D_LocalEnergy[Point_Index][r][s][0];
                i_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][1];
                j_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][2];
                k_G           = Local_Points_1D_LocalEnergy[Point_Index][r][s][3];
                v             = Local_Points_1D_LocalEnergy[Point_Index][r][s][4];

                if(Boundary_Type!=0)
                {
                    EV[p][q].M_Local[1][o] = &Mx_LLG[i_G][j_G][k_G];
                    EV[p][q].M_Local[2][o] = &My_LLG[i_G][j_G][k_G];
                    EV[p][q].M_Local[3][o] = &Mz_LLG[i_G][j_G][k_G];

                }else
                {
                    EV[p][q].M_Local[1][o] = &MxV_LLG[v];
                    EV[p][q].M_Local[2][o] = &MyV_LLG[v];
                    EV[p][q].M_Local[3][o] = &MzV_LLG[v];
                }
            }

            EV[p][q].M0[0] = *EV[p][q].M_Local[1][1];
            EV[p][q].M0[1] = *EV[p][q].M_Local[2][1];
            EV[p][q].M0[2] = *EV[p][q].M_Local[3][1];
        }
    }
///////////////////////// Free Part /////////////////////////
    free3DArrayinteger(BS->Boundary,m,n);
    free4DArrayinteger(BS->LocalEnergy,m,n,l);
    free(CPP->Virtual_Position);
    free4DArrayinteger(Local_Points_1D_LocalEnergy,BS->Calculate_Points,7,7);
    free2DArrayinteger(CPP->Neighbour_Points_3Dto1D_Array,14);
///////////////////////// Synchronize Parameters Initial /////////////////////////
    Modify_LLG_Synchronize *MLS;
    MLS = &DRLP.MLS;
    MLS->m=m;    MLS->n=n;    MLS->l=l;    MLS->Virtual_Points  = BS->Virtual_Points;
    MLS->Inner_Points=BS->Inner_Points;  MLS->Boundary_Points = BS->Boundary_Points;
    MLS->Pthread_Number=Threads_Number;
    MLS->Bound_CalculatePoint_Threads = Make2DArrayinteger(Threads_Number,3);
    MLS->Signal                       = Signal;
    MLS->PTPS                         = PTPS;
    MLS->IPN = *IPN;

    MLS->h   = h;      MLS->t = t;
    MLS->dxy = dxy;
    MLS->Sort_Calculate_Points        = Make1DArrayinteger(BS->Calculate_Points);
    MLS->Virtual_Points               = BS->Virtual_Points;
    MLS->Version                      = Version;
    MLS->CST                          = CST;
    MLS->Local_Points_1D              = CPP->Local_Points_1D;
    MLS->Position_3DTo1D              = CPP->Position_3DTo1D;
    MLS->Energy_Coefficient           = IPN->Energy_Coefficient;

    for(p=0;p<BS->Calculate_Points;++p)
    {
        MLS->Sort_Calculate_Points[p]=CPP->Sort_Calculate_Points[p];
    }

    for(p=0;p<Threads_Number;++p)
    {
        MLS->Bound_CalculatePoint_Threads[p][0] = MLPt[p].CalculatePointsThreads_Upper-MLPt[p].CalculatePointsThreads_Lower;
        MLS->Bound_CalculatePoint_Threads[p][1] = MLPt[p].CalculatePointsThreads_Lower;
        MLS->Bound_CalculatePoint_Threads[p][2] = MLPt[p].CalculatePointsThreads_Upper;
    }

    DRLP.PTPS = PTPS;
    DRLP.EV = EV;
    DRLP.Signal = Signal;
    DRLP.M_Local = M_Local;
    DRLP.CST = CST;
    DRLP.MLPt = MLPt;

    printf("Run\n");
    ///////////////////////// Pthreads Part /////////////////////////
    int err;
    pthread_t *tid;
    tid = (pthread_t *)malloc(sizeof(pthread_t) * (Threads_Number + 2));
    DRLP.tid = tid;
    for(p=0;p<Threads_Number;++p)
    {
        err=pthread_create(&tid[p],NULL,Modify_Pthread_LLG_Transfer,&MLPt[p]);
        if(err != 0)
        {
            printf("create thread error\n");
            Sleep(1000);
        }
    }

    err=pthread_create(&tid[Threads_Number],NULL,Modify_Pthread_LLG_Synchronize_Transfer,MLS);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    ///////////////////////// OPGL Part /////////////////////////
    OPGL_GUI.State_Syn = StateSyn_Run;
    OPGL_LLG_Pthread_Parameters *OLPP;
    OLPP = &DRLP.OLPP;
    OLPP->Mx_OPGL = Mx_OPGL; 
    OLPP->My_OPGL = My_OPGL;
    OLPP->Mz_OPGL = Mz_OPGL;
    OLPP->MxV_OPGL = MxV_OPGL; 
    OLPP->MyV_OPGL = MyV_OPGL;
    OLPP->MzV_OPGL = MzV_OPGL;
    OLPP->IPN = IPN;
    OLPP->CST = Make4DArray(7,m,n,l);
    for (i = 0; i < m; ++i)
    {
        for (j = 0; j < n; ++j)
        {
            for (k = 0; k < l; ++k)
            {
                for (r = 0; r < 7; ++r)
                {
                    OLPP->CST[r][i][j][k] = CST[r][i][j][k];
                }
            }
        }
    }

    OPGL_Mode_GUI OMG;
    OMG = (OPGL_Mode_GUI)IPN->OPGL_Mode;
    switch (OMG)
    {
        case OPGL_Mode_On:
        {
            err = pthread_create(&tid[Threads_Number + 1], NULL, OPGL_LLG_Pthread_Transfer, OLPP);
            if(err != 0)
            {
                printf("create thread error\n");
                Sleep(1000);
            }
        }
            
        default:
        break;
    }

    ////////////////////// Free pthread //////////////////////////////////////////
    pthread_t tid_Part_Free;
    err = pthread_create(&tid_Part_Free, NULL, DiscreatPoint_Run_LLG_Part_Free_Transfer, IPN);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
    return 1;
}

int DiscreatPoint_Run_LLG_Part_Free(Input_Parameter_NMD *IPN)
{
    int i, j, k, v, m, n, l, p, q, CalculatePointsThreads, N_Thread;
    OPGL_Mode_GUI OMG;
    OMG = (OPGL_Mode_GUI)IPN->OPGL_Mode;
    switch (OMG)
    {
        case OPGL_Mode_On:
        {
            N_Thread = IPN->Threads_Number + 2;
            break;
        }

        case OPGL_Mode_Off:
        {
            N_Thread = IPN->Threads_Number + 1;
            break;
        }
            
        default:
        break;
    }

    for (p = 0; p < N_Thread; ++p)
    {
        pthread_join(DRLP.tid[p],NULL);
    }
    printf("free thread\n");
    m = 2 * IPN->xn + 3;    n = 2 * IPN->yn + 3;
    if (IPN->zn == 0)
    {
        l = 1;
    }else
    {
        l = 2 * IPN->zn + 3;
    }

    ///////////////////////// Result Part /////////////////////////

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                Mx_Global[i][j][k]=Mx_LLG[i][j][k];
                My_Global[i][j][k]=My_LLG[i][j][k];
                Mz_Global[i][j][k]=Mz_LLG[i][j][k];
            }
        }
    }

    for (v = 0; v < DRLP.BS.Virtual_Points; ++v)
    {
        MxV_Global[v] = MxV_LLG[v];
        MyV_Global[v] = MyV_LLG[v];
        MzV_Global[v] = MzV_LLG[v];
    }
    ///////////////////////// Result Part /////////////////////////
    double W, SKN, DerivativeM;
    W=0.; SKN=0.;
    DerivativeM = 0.;
    for (p = 0; p < IPN->Threads_Number; ++p)
    {
        W   = W + DRLP.MLPt[p].W;
        SKN = SKN + DRLP.MLPt[p].SKN;
        DerivativeM = DerivativeM + DRLP.MLPt[p].DerivativeM;
    }
    W   = W/(DRLP.BS.Boundary_Points+DRLP.BS.Inner_Points);
    IPN->W = W; IPN->SKN = SKN;
    DerivativeM = DerivativeM / (DRLP.BS.Boundary_Points + DRLP.BS.Inner_Points);

    Energy_Variables_Distribution EVD;
    EVD.CST   = DRLP.CST;    EVD.dxy = IPN->dxy;
    EVD.h     = IPN->h;      EVD.t   = IPN->t;
    EVD.m     = m;           EVD.n   = n;
    EVD.l     = l;
    EVD.Mx    = Mx_Global;
    EVD.My    = My_Global;
    EVD.Mz    = Mz_Global;
    EVD.MxV   = MxV_Global;
    EVD.MyV   = MyV_Global;
    EVD.MzV   = MzV_Global;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D = DRLP.CPP.Local_Points_1D;
    EVD.Position_3DTo1D = DRLP.CPP.Position_3DTo1D;
    EVD.Energy_Coefficient = IPN->Energy_Coefficient;
    Energy_Distribution(&EVD);
    Write_NMD_Energy(IPN,EVD.w);
    Write_NMD_SKN_D(IPN, EVD.SKN_D);
    IPN->W   = EVD.W_Average;
    IPN->SKN = EVD.SKN;
    IPN->SKN_Positif = EVD.SKN_Positif;
    IPN->SKN_Area = EVD.SKN_Area;
    IPN->DerivativeM = DerivativeM;
    printf("%lf\t%lf\t%lf\t%lf\n",IPN->W,IPN->SKN,IPN->SKN_Positif,IPN->SKN_Area);

    ///////////////////////// Free Part /////////////////////////
    pthread_mutex_destroy(&OPGL_Mutex);
    for(p=0;p<IPN->Threads_Number;++p)
    {
        free(DRLP.MLPt[p].Signal);

        CalculatePointsThreads = DRLP.MLPt[p].CalculatePointsThreads_Upper - DRLP.MLPt[p].CalculatePointsThreads_Lower;
        for(q=0;q<CalculatePointsThreads;++q)
        {
            free(DRLP.EV[p][q].dxy);
            free(DRLP.EV[p][q].M0);
            free(DRLP.EV[p][q].LocalEnergy);
            free2DArray(DRLP.EV[p][q].CST,7);
            free2DArrayPointer(DRLP.EV[p][q].M_Local,4);
            free(DRLP.EV[p][q].Energy_Coefficient);
        }

        for (q = 1; q < 4; ++q)
        {
            pthread_mutex_destroy(&Signal_Synchronize_LLG_Mutex[p][q]);
            pthread_cond_destroy(&Signal_Synchronize_LLG_Cond[p][q]);
        }
        free(DRLP.MLPt[p].MLP);
        free(DRLP.EV[p]);
        free(Signal_Synchronize_LLG_Mutex[p]);
        free(Signal_Synchronize_LLG_Cond[p]);
        free(DRLP.MLPt[p].Sort_Calculate_Points);

        pthread_mutex_destroy(&M_Synchronize_LLG_Mutex[p]);
    }
    free4DArray(DRLP.CST,7,m,n);
    free2DArrayinteger(DRLP.MLS.Bound_CalculatePoint_Threads,IPN->Threads_Number);
    free3DArray(Mx_LLG,m,n);
    free3DArray(My_LLG,m,n);
    free3DArray(Mz_LLG,m,n);
    free(MxV_LLG);
    free(MyV_LLG);
    free(MzV_LLG);
    free(DRLP.MLS.Sort_Calculate_Points);
    free3DArray(DRLP.M_Local,IPN->Threads_Number,4);
    free3DArrayinteger(DRLP.CPP.Local_Points_1D,DRLP.BS.Calculate_Points,7);
    free3DArrayinteger(DRLP.CPP.Position_3DTo1D,m,n);
    free(DRLP.CPP.Sort_Calculate_Points);
    free(DRLP.MLPt);
    free(DRLP.EV);
    free2DArrayintegerPointer(DRLP.Signal, IPN->Threads_Number + 1);
    free(Signal_Synchronize_LLG_Mutex);
    free(Signal_Synchronize_LLG_Cond);
    free(M_Synchronize_LLG_Mutex);
    free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
    free4DArray(DRLP.OLPP.CST,7,m,n);
    free(DRLP.PTPS);
    ///////////////////////// Result Part /////////////////////////
    double end_time, duration;
    end_time  = clock();
    duration  = (double)(end_time-DRLP.start_time)/CLOCKS_PER_SEC;
    IPN->Time = duration;

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Finish;
    pthread_cond_signal(&GPM.Run_GUI_Cond);
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    printf("free done\n");

    return 1;
}

void* DiscreatPoint_Run_LLG_Part_Free_Transfer(void *arg)
{
    Input_Parameter_NMD *IPN;
    IPN=(Input_Parameter_NMD*)arg;
    DiscreatPoint_Run_LLG_Part_Free(IPN);
    return NULL;
}

int DiscreatPoint_Run_LLG(Input_Parameter_NMD *IPN)
{
    int State_GUI_Run;
    Thread_State_GUI TSG;
    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    State_GUI_Run = GPM.State_GUI_Run;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    TSG = (Thread_State_GUI)State_GUI_Run;
    switch (TSG)
    {
        case TSG_Prepare:
        {
            pthread_mutex_lock(&GPM.Run_GUI_Mutex);
            GPM.State_GUI_Run = TSG_Run;
            pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
            DiscreatPoint_Run_LLG_Part_Run(IPN);
            break;
        }

        case TSG_Run:
        {
            break;
        }

        case TSG_Finish:
        {
            break;
        }
        
        default:
        break;
    }
    return 1;
}