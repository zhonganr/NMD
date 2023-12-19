#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <windows.h>
#include <omp.h>
#include <dir.h>
#include <InOut.h>
#include <Array.h>
#include <BoundaryShape.h>
#include <Variable.h>
#include <Modify_MC.h>
#include <Modify_NL.h>
#include <Modify_Target.h>
#include <malloc.h>
#include <LLG.h>
#include <GNEB.h>
#include <Initial_State.h>
#include <Elastic_Field.h>
#include <Energy.h>
#include <OPGL.h>
#include <gui.h>

#include <QApplication>

double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
extern double Energy_Variation_Image;
int i_Angle_Rotation;

extern "C" __declspec(dllexport) int AmdPowerXpressRequestHighPerformance = 0x00000001;
extern "C" __declspec(dllexport) DWORD NvOptimusEnablement = 0x00000001;
using namespace std;
Input_Parameter_NMD *IPN;
extern GUI_Parameter GPM;
extern int N_Energy_Terms;
int main(int argc, char *argv[])
{
    int xn, yn, zn, m, n, l, i, j, k, p;
    int rn=20;

    xn=1; yn=1; zn=1;
    m=2*xn+3;   n=2*yn+3;   l=2*zn+3;
    int Version, InitialMode, NCycles, Threads_Number, ElasticFieldMode;
    int Boundary_Shape_Type, CalculateMode, Boundary_Type_x, Boundary_Type_y, Boundary_Type_z;
    int Continuity_Modify_Mode, InPlane_Density, N_Image;
    int MagneticFieldMode, MagneticFieldPeriod, TemperatureFieldMode, TemperatureFieldPeriod;
    double *dxy, *CST_Uniaxial, h, t, ****CST, Continuity_Modify_Coefficient, dT, CurrentDensity, MpX, MpY, MpZ;
    double b, a, FC, FU_Beam, FC_Beam, B_M, SKN, ElasticFieldScaleFactor;
    char *dir_SKN, *file_SKN;
    dir_SKN = Make1DString(500);
    file_SKN = Make1DString(500);

    //Input_Parameter_NMD IPN0, *IPN;
    Input_Parameter_NMD IPN0;
    IPN=&IPN0;
    IPN->dxy=Make1DArray(3);
    IPN->CST_Uniaxial=Make1DArray(6);
    IPN->Energy_Coefficient=Make1DArray(N_Energy_Terms);

    //ReadInput_NMD(IPN);
    //WriteInput_NMD(IPN);
    ReadInput_NMD(IPN);
    //printf("%d\n", IPN->NCycles);

    m=2*IPN->xn+3;   n=2*IPN->yn+3;
    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l=2*IPN->zn+3;
    }

/*
    do
    {
        SKN = Skyrmion_Number_Incircle(IPN);
        printf("%0.8f\n",SKN);
    }while(1);
*/
    //SKN_Boundary_Found(IPN);

	
/*
    FILE *fp_Note;
    Version = IPN->Version;
    for (j = 1; j < 2; ++j)
    {
        //h = 2.7 - j * 0.05;
        //sprintf(file_SKN, "h=%0.2f\\Note.dat", h);
        sprintf(file_SKN, "Note.dat");
        fp_Note = fopen(file_SKN, "w+");
        fprintf(fp_Note,"Number\t Energy\t SKN\n");
        for (i = 1; i < 12; ++i)
        {
            //IPN->h = h;
            if(i==0)
            {
                //sprintf(dir_SKN, "h=%0.2f\\FE", h);
            }else
            {
                //sprintf(dir_SKN, "h=%0.2f\\%d", h, i);
                sprintf(dir_SKN, "%d", i);
            }
            IPN->Version=Version;
            InitialM_3D_SKN(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, BS, dir_SKN);
            Trim_Boundary_M3D_BoundaryToVirtual(BS,CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);
            IPN->Version=Version+1;
            Make_NMD_Result_Dir_SKN(IPN,dir_SKN);
            if(IPN->CalculateMode==0)
            {
                DiscreatPoint_Run_NL(IPN);
            }else if(IPN->CalculateMode==1)
            {
                DiscreatPoint_Run_MC(IPN);
            }else if(IPN->CalculateMode==2)
            {
                DiscreatPoint_Run_LLG(IPN);
            }
            printf("%d\t %0.8f\t %0.8f\n",i,IPN->W,IPN->SKN);
            fprintf(fp_Note,"%d\t %0.8f\t %0.8f\n",i,IPN->W,IPN->SKN);
            Sleep(100);
            Trim_Boundary_M3D_VirtualToBoundary(BS,CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);
            Write_NMD_Magnetization_SKN(IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global,BS,dir_SKN);
            Write_Inplane_M_SKN(IPN, BS, Mx_Global, My_Global, dir_SKN);
        }
        fclose(fp_Note);       
    }
    WriteInput_NMD(IPN);
*/
    //OPGL_GNEB(IPN);
/////////////////////////////////////////////// Rotation ///////////////////////////////////////////////////
/*
    InitialM_3D(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, BS);
    Version = IPN->Version;
    IPN->Version=IPN->Version+1;
    Make_NMD_Result_Dir(IPN);
    FILE *fp_Rotation_Result;
    fp_Rotation_Result = fopen("Rotation_Result.dat", "w+");
    int N_Angle = 15;
    double d_Angle = 0.1, Angle_Rotation;
    fprintf(fp_Rotation_Result, "Angle\t dE\n");
    for (i = 0; i < N_Angle; ++i)
    {
        IPN->Version = Version;
        InitialM_3D(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, BS);
        i_Angle_Rotation = i;
        IPN->Version = Version + 1;
        Angle_Rotation = d_Angle * i;
        
        Magnetization_Transfer_By_Rotation(IPN, Mx_Global, My_Global, Mz_Global, Angle_Rotation);
        Trim_Boundary_M3D_BoundaryToVirtual(BS, CPP, IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global);
        NCycles = IPN->NCycles;
        IPN->NCycles = 5000;
        IPN->Continuity_Modify_Mode = 1;
        DiscreatPoint_Run_MC(IPN);
        IPN->NCycles = NCycles;
        IPN->Continuity_Modify_Mode = 0;
        DiscreatPoint_Run_GNEB(IPN);
        fprintf(fp_Rotation_Result, "%0.2f\t %0.6f\n", Angle_Rotation, Energy_Variation_Image);
    }
    fclose(fp_Rotation_Result);
*/
/*    
    InitialM_3D(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, BS);
    Trim_Boundary_M3D_BoundaryToVirtual(BS,CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);

    IPN->Version=IPN->Version+1;
    Make_NMD_Result_Dir(IPN);
    //OPGL_NMD(Mx_Global, My_Global, Mz_Global, m, n, InPlane_Density);
    printf("%d\n", IPN->InitialMode);

    if(IPN->CalculateMode==0)
    {
        DiscreatPoint_Run_NL(IPN);
    }else if(IPN->CalculateMode==1)
    {
        DiscreatPoint_Run_MC(IPN);
    }else if(IPN->CalculateMode==2)
    {
        DiscreatPoint_Run_LLG_SOT(IPN);
    }else if(IPN->CalculateMode==3)
    {
        DiscreatPoint_Run_LLG_STT(IPN);
    }else if(IPN->CalculateMode==4)
    {
        DiscreatPoint_Run_Target(IPN);
    }else if(IPN->CalculateMode==5)
    {
        DiscreatPoint_Run_GNEB(IPN);
    }

    //
    Trim_Boundary_M3D_VirtualToBoundary(BS,CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);
    Write_NMD_Magnetization(IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global,BS);
    CST=Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);

    //ReadAnsys(IPN,CST);

    Write_NMD_ElastiField(IPN,CST,BS);
    Write_Principal_Stress(IPN,CST,BS);
    Write_Principal_Strain(IPN,CST,BS);


    Write_NMD_Result(IPN);
    Write_Inplane_M(IPN,BS,Mx_Global,My_Global);
    Write_Judge_Phase_Curve(IPN,BS,CST);
    free4DArray(CST,7,m,n);
    WriteInput_NMD(IPN);
*/
    //WriteVectorM_3D_Monolayer(m,n,Mx_Global,My_Global,Mz_Global);
/*
    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                printf("%d\t%d\t%d\t%d\n",i,j,k,BS->Boundary[i][j][k]);
            }
        }
    }
*/
    QApplication a_GUI(argc, argv);
    GUI w;
    w.show();
    //free(IPN->dxy);
    //free(IPN->CST_Uniaxial);
    //free(IPN->Energy_Coefficient);
    free(dir_SKN);

    //cout << "Hello world!" << endl;
    
    return a_GUI.exec();
}
