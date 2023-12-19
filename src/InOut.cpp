#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <direct.h>
#include <Array.h>
#include <Elastic_Field.h>
#include <Variable.h>
#include <Initial_State.h>
#include <InOut.h>
#include <BoundaryShape.h>
#define PI      3.14159265359
extern double q_Stress;
extern int i_Angle_Rotation;

using namespace std;
int WriteInput_NMD(Input_Parameter_NMD *IPN)
{
    int xn, yn, zn, Version, VersionMode, InitialMode, NCycles, Threads_Number, ElasticFieldMode;
    int Boundary_Shape_Type, Boundary_Type_x, Boundary_Type_y, Boundary_Type_z, CalculateMode;
    int MagneticFieldMode, MagneticFieldPeriod, TemperatureFieldMode, TemperatureFieldPeriod;
    int Continuity_Modify_Mode, OPGL_Mode;
    double h, t, dx, dy, dz, Continuity_Modify_Coefficient;
    double b, a, FC, FU_Beam, FC_Beam, B_M, ElasticFieldScaleFactor;
    double dT, MpX, MpY, MpZ, CurrentDensity;
    double *Energy_Coefficient;

    xn = IPN->xn;  yn = IPN->yn;
    zn = IPN->zn;  b  = IPN->b;
    a  = IPN->a;   FC = IPN->FC;
    h  = IPN->h;   t  = IPN->t;
    FU_Beam     = IPN->FU_Beam;
    FC_Beam     = IPN->FC_Beam;
    B_M         = IPN->B_M;
    Version     = IPN->Version;
    NCycles     = IPN->NCycles;
    InitialMode = IPN->InitialMode;
    dx          = IPN->dxy[0];
    dy          = IPN->dxy[1];
    dz          = IPN->dxy[2];
    VersionMode            = IPN->VersionMode;
    Threads_Number         = IPN->Threads_Number;
    CalculateMode          = IPN->CalculateMode;
    ElasticFieldMode       = IPN->ElasticFieldMode;
    MagneticFieldMode      = IPN->MagneticFieldMode;
    MagneticFieldPeriod    = IPN->MagneticFieldPeriod;
    TemperatureFieldMode   = IPN->TemperatureFieldMode;
    TemperatureFieldPeriod = IPN->TemperatureFieldPeriod;
    Boundary_Shape_Type    = IPN->Boundary_Shape_Type;
    Boundary_Type_x        = IPN->Boundary_Type_x;
    Boundary_Type_y        = IPN->Boundary_Type_y;
    Boundary_Type_z        = IPN->Boundary_Type_z;
    Continuity_Modify_Mode        = IPN->Continuity_Modify_Mode;
    Continuity_Modify_Coefficient = IPN->Continuity_Modify_Coefficient;
    ElasticFieldScaleFactor       = IPN->ElasticFieldScaleFactor;
    Energy_Coefficient            = IPN->Energy_Coefficient;
    dT          = IPN->dT;
    MpX         = IPN->MpX;
    MpY         = IPN->MpY;
    MpZ         = IPN->MpZ;
    CurrentDensity = IPN->CurrentDensity;
    OPGL_Mode = IPN->OPGL_Mode;

    FILE *fp;

    fp=fopen("Input_NMD.dat", "w+");
    fprintf(fp,"hex\t\t t\n");
    fprintf(fp,"%g\t\t %g\n",h,t);
    fprintf(fp,"MagneticFieldMode\t MagneticFieldPeriod\t TemperatureFieldMode\t TemperatureFieldPeriod\n");
    fprintf(fp, "%d\t\t\t %d\t\t\t %d\t\t\t %d\n", MagneticFieldMode, MagneticFieldPeriod, TemperatureFieldMode, TemperatureFieldPeriod);

    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"Exchange\t\t DM_Bloch\t DM_Neel\t Zeeman\t\t Landau\t\t ME\t Axial_Anisotropy\n");
    fprintf(fp, "%0.2f\t\t %0.2f\t\t %0.2f\t\t %0.2f\t\t %0.2f\t\t %0.2f\t %0.2f\t\n", Energy_Coefficient[0], Energy_Coefficient[1],
            Energy_Coefficient[2], Energy_Coefficient[3], Energy_Coefficient[4], Energy_Coefficient[5], Energy_Coefficient[6]);

    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"xn\t\t yn\t\t zn\n");
    fprintf(fp,"%d\t\t %d\t\t %d\n",xn,yn,zn);

    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"Boundary Shape Type\n");
    fprintf(fp,"%d\n",Boundary_Shape_Type);
    fprintf(fp,"Boundary Type\n");
    fprintf(fp,"x\t\t y\t\t z\t\n");
    fprintf(fp,"%d\t\t %d\t\t %d\n",Boundary_Type_x,Boundary_Type_y,Boundary_Type_z);

    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"Continuity Modify Mode\n");
    fprintf(fp,"%d\n",Continuity_Modify_Mode);
    fprintf(fp,"Continuity Modify Coefficient\n");
    fprintf(fp,"%0.8f\n",Continuity_Modify_Coefficient);

    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"e11\t\t %lf\n",IPN->CST_Uniaxial[0]);
    fprintf(fp,"e22\t\t %lf\n",IPN->CST_Uniaxial[1]);
    fprintf(fp,"e33\t\t %lf\n",IPN->CST_Uniaxial[2]);
    fprintf(fp,"e12\t\t %lf\n",IPN->CST_Uniaxial[3]);
    fprintf(fp,"e13\t\t %lf\n",IPN->CST_Uniaxial[4]);
    fprintf(fp,"e23\t\t %lf\n",IPN->CST_Uniaxial[5]);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"Version=\t%d\n",Version);
    fprintf(fp,"VersionMode=\t%d\n",VersionMode);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"CalculateMode=\t%d\n",CalculateMode);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"ThreadsNumber=\t%d\n",Threads_Number);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"dT=\t%g\n",dT);
    fprintf(fp,"Current Density=\t%g\n",CurrentDensity);
    fprintf(fp,"MpX\t\t MpY\t\t MpZ\n");
    fprintf(fp,"%g\t\t %g\t\t %g\n",MpX,MpY,MpZ);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"NCycles=\t%d\n",NCycles);
    fprintf(fp,"InitialMode=\t%d\n",InitialMode);
    fprintf(fp,"ElasticFieldMode=\t%d\n",ElasticFieldMode);
    fprintf(fp,"ElasticFieldScaleFactor=\t%0.3f\n",ElasticFieldScaleFactor);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"b\t\t a\t\t FC\t\t FU_Beam\t\t FC_Beam\t\t B_M\n");
    fprintf(fp,"%g\t\t %g\t\t %g\t\t %g\t\t %g\t\t %g\n",b,a,FC,FU_Beam,FC_Beam,B_M);
    fprintf(fp,"----------------------------------\n");

    fprintf(fp,"dx\t\t dy\t\t dz\n");
    fprintf(fp,"%0.8f\t %0.8f\t %0.8f\n",dx,dy,dz);

    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"In Plane Components Density=\t%d\n",IPN->Inplane_Density);
    fprintf(fp,"N_Image=\t%d\n",IPN->N_Image);
    fprintf(fp,"OPGLMode=\t%d\n",OPGL_Mode);

    fclose(fp);

    return 1;

}

int ReadInput_NMD(Input_Parameter_NMD *IPN)
{
    char s[1024];
    int xn, yn, zn, Version, VersionMode, InitialMode, NCycles, Threads_Number, ElasticFieldMode;
    int Boundary_Shape_Type, Boundary_Type_x, Boundary_Type_y, Boundary_Type_z, CalculateMode;
    int MagneticFieldMode, MagneticFieldPeriod, TemperatureFieldMode, TemperatureFieldPeriod;
    int Continuity_Modify_Mode, Inplane_Density, N_Image, OPGL_Mode;
    double *dxy, *CST_Uniaxial, h, t, Continuity_Modify_Coefficient;
    double b, a, FC, FU_Beam, FC_Beam, B_M, ElasticFieldScaleFactor;
    double dT, MpX, MpY, MpZ, CurrentDensity;
    double *Energy_Coefficient;

    dxy = IPN->dxy; CST_Uniaxial = IPN->CST_Uniaxial;
    Energy_Coefficient = IPN->Energy_Coefficient;

    FILE *fp;

    fp=fopen("Input_NMD.dat", "r+");
    fscanf(fp,"%s\t %s",s,s);
    fscanf(fp,"%lf\t %lf",&h,&t);
    fscanf(fp, "%s %s %s %s", s, s, s, s);
    fscanf(fp, "%d %d %d %d", &MagneticFieldMode, &MagneticFieldPeriod, &TemperatureFieldMode, &TemperatureFieldPeriod);

    fscanf(fp,"%s",s);
    fscanf(fp, "%s %s %s %s %s %s %s", s, s, s, s, s, s, s);
    fscanf(fp, "%lf\t %lf\t %lf\t %lf\t %lf\t %lf\t %lf", &Energy_Coefficient[0], &Energy_Coefficient[1], &Energy_Coefficient[2],
           &Energy_Coefficient[3], &Energy_Coefficient[4], &Energy_Coefficient[5], &Energy_Coefficient[6]);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t %s\t %s",s,s,s);
    fscanf(fp,"%d\t %d\t %d",&xn,&yn,&zn);
//////////////////////////////////////////////////////
    fscanf(fp,"%s",s);
    fscanf(fp,"%s %s %s",s,s,s);
    fscanf(fp,"%d",&Boundary_Shape_Type);
    fscanf(fp,"%s %s",s,s);
    fscanf(fp,"%s %s %s",s,s,s);
    fscanf(fp,"%d\t %d\t %d",&Boundary_Type_x,&Boundary_Type_y,&Boundary_Type_z);
//////////////////////////////////////////////////////////
    fscanf(fp,"%s",s);
    fscanf(fp,"%s %s %s",s,s,s);
    fscanf(fp,"%d",&Continuity_Modify_Mode);
    fscanf(fp,"%s %s %s",s,s,s);
    fscanf(fp,"%lf",&Continuity_Modify_Coefficient);
/////////////////////////////////////////////////////////////////////////
    fscanf(fp,"%s",s);

    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[0]);
    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[1]);
    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[2]);
    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[3]);
    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[4]);
    fscanf(fp,"%s\t%lf",s,&CST_Uniaxial[5]);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t%d",s,&Version);
    fscanf(fp,"%s\t%d",s,&VersionMode);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t%d",s,&CalculateMode);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t%d",s,&Threads_Number);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t%lf",s,&dT);
    fscanf(fp,"%s %s\t%lf",s,s,&CurrentDensity);
    fscanf(fp,"%s\t %s\t %s",s,s,s);
    fscanf(fp,"%lf\t %lf\t %lf",&MpX,&MpY,&MpZ);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t%d",s,&NCycles);
    fscanf(fp,"%s\t%d",s,&InitialMode);
    fscanf(fp,"%s\t%d",s,&ElasticFieldMode);
    fscanf(fp,"%s\t%lf",s,&ElasticFieldScaleFactor);

    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t %s\t %s\t %s\t %s\t %s",s,s,s,s,s,s);
    fscanf(fp,"%lf\t %lf\t %lf\t %lf\t %lf\t %lf",&b,&a,&FC,&FU_Beam,&FC_Beam,&B_M);
    //printf("%g\t%g\t%g\n",b,a,FC);
    fscanf(fp,"%s",s);
    fscanf(fp,"%s\t %s\t %s\n",s,s,s);
    fscanf(fp,"%lf\t %lf\t %lf\n",&dxy[0],&dxy[1],&dxy[2]);
    fscanf(fp,"%s",s);
    fscanf(fp,"%s %s %s %s %d",s,s,s,s,&Inplane_Density);
    fscanf(fp, "%s %d", s, &N_Image);
    fscanf(fp, "%s %d", s, &OPGL_Mode);
    /*
    printf("%lf\t%lf\n",h,t);
    printf("%d\t%d\t%d\n",xn,yn,zn);
    printf("%d\n",Boundary_Type);
    printf("%d\n",Version);
    printf("%d\n",CalculateMode);
    printf("%d\n",Threads_Number);
    printf("%d\n",NCycles);
    printf("%d\n",InitialMode);
    printf("%d\n",ElasticFieldMode);
*/
    IPN->xn=xn;  IPN->yn=yn;
    IPN->zn=zn;
    IPN->h =h;   IPN->t =t;

    IPN->b            = b;
    IPN->a            = a;
    IPN->FC           = FC;
    IPN->FU_Beam      = FU_Beam;
    IPN->FC_Beam      = FC_Beam;
    IPN->B_M          = B_M;
    IPN->Version      = Version;
    IPN->VersionMode  = VersionMode;
    IPN->NCycles      = NCycles;
    IPN->InitialMode  = InitialMode;
    IPN->ElasticFieldMode       = ElasticFieldMode;
    IPN->MagneticFieldMode      = MagneticFieldMode;
    IPN->MagneticFieldPeriod    = MagneticFieldPeriod;
    IPN->TemperatureFieldMode   = TemperatureFieldMode;
    IPN->TemperatureFieldPeriod = TemperatureFieldPeriod;
    IPN->Threads_Number         = Threads_Number;
    IPN->Boundary_Shape_Type    = Boundary_Shape_Type;
    IPN->Boundary_Type_x        = Boundary_Type_x;
    IPN->Boundary_Type_y        = Boundary_Type_y;
    IPN->Boundary_Type_z        = Boundary_Type_z;
    IPN->CalculateMode          = CalculateMode;
    IPN->Continuity_Modify_Mode = Continuity_Modify_Mode;
    IPN->Continuity_Modify_Coefficient = Continuity_Modify_Coefficient;
    IPN->ElasticFieldScaleFactor       = ElasticFieldScaleFactor;
    IPN->Inplane_Density = Inplane_Density;
    IPN->dT            = dT;
    IPN->MpX           = MpX ;
    IPN->MpY           = MpY ;
    IPN->MpZ           = MpZ ;
    IPN->CurrentDensity=CurrentDensity;
    IPN->N_Image=N_Image;
    IPN->OPGL_Mode = OPGL_Mode;

    fclose(fp);

    return 1;
}

int Make_NMD_Result_Dir(Input_Parameter_NMD *IPN)
{
    char filename[500];
    int k, zn, zc, l;
    zn=IPN->zn;
    if(zn==0)
    {
        zc=0;
        l=1;
    }else
    {
        zc = zn+1;
        l  = 2*zn+3;
    }

    sprintf(filename,"Version %d",IPN->Version);
    mkdir(filename);
    sprintf(filename,"Version %d\\M3D",IPN->Version);
    mkdir(filename);
    sprintf(filename,"Version %d\\Elastic Field",IPN->Version);
    mkdir(filename);
    sprintf(filename,"Version %d\\Elastic Field\\Stress",IPN->Version);
    mkdir(filename);
    sprintf(filename,"Version %d\\Elastic Field\\Strain",IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\M3D\\z=%d",IPN->Version,k-zc);
        mkdir(filename);
        sprintf(filename,"Version %d\\Elastic Field\\Stress\\z=%d",IPN->Version,k-zc);
        mkdir(filename);

        sprintf(filename,"Version %d\\Elastic Field\\Strain\\z=%d",IPN->Version,k-zc);
        mkdir(filename);

    }
    sprintf(filename,"Version %d\\Process",IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\Process\\z=%d",IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"Version %d\\In-Plane Components",IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\In-Plane Components\\z=%d",IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"Version %d\\For Printing",IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\For Printing\\z=%d",IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"Version %d\\For Printing\\Elastic Field",IPN->Version);
    mkdir(filename);

    sprintf(filename,"Version %d\\For Printing\\Elastic Field\\Stress",IPN->Version);
    mkdir(filename);
    sprintf(filename,"Version %d\\For Printing\\Elastic Field\\Strain",IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d",IPN->Version,k-zc);
        mkdir(filename);

        sprintf(filename,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d",IPN->Version,k-zc);
        mkdir(filename);

    }

    return 1;
}

int Make_NMD_Result_Dir_SKN(Input_Parameter_NMD *IPN, char *dir_SKN)
{
    char filename[500];
    int k, zn, zc, l;
    zn=IPN->zn;
    if(zn==0)
    {
        zc=0;
        l=1;
    }else
    {
        zc = zn+1;
        l  = 2*zn+3;
    }

    sprintf(filename,"%s",dir_SKN);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d",dir_SKN,IPN->Version);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d\\M3D",dir_SKN,IPN->Version);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d\\Elastic Field",dir_SKN,IPN->Version);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d\\Elastic Field\\Stress",dir_SKN,IPN->Version);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d\\Elastic Field\\Strain",dir_SKN,IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"%s\\Version %d\\M3D\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);
        sprintf(filename,"%s\\Version %d\\Elastic Field\\Stress\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);

        sprintf(filename,"%s\\Version %d\\Elastic Field\\Strain\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);

    }
    sprintf(filename,"%s\\Version %d\\Process",dir_SKN,IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"%s\\Version %d\\Process\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"%s\\Version %d\\In-Plane Components",dir_SKN,IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"%s\\Version %d\\In-Plane Components\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"%s\\Version %d\\For Printing",dir_SKN,IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"%s\\Version %d\\For Printing\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);
    }

    sprintf(filename,"%s\\Version %d\\For Printing\\Elastic Field",dir_SKN,IPN->Version);
    mkdir(filename);

    sprintf(filename,"%s\\Version %d\\For Printing\\Elastic Field\\Stress",dir_SKN,IPN->Version);
    mkdir(filename);
    sprintf(filename,"%s\\Version %d\\For Printing\\Elastic Field\\Strain",dir_SKN,IPN->Version);
    mkdir(filename);
    for(k=0;k<l;k++)
    {
        sprintf(filename,"%s\\Version %d\\For Printing\\Elastic Field\\Stress\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);

        sprintf(filename,"%s\\Version %d\\For Printing\\Elastic Field\\Strain\\z=%d",dir_SKN,IPN->Version,k-zc);
        mkdir(filename);

    }

    return 1;
}

int Write_NMD_ElastiField(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, zc, r, Version;
    double *Strain, *Stress;
    char *filename_Strain_S11, *filename_Strain_S22, *filename_Strain_S33,
         *filename_Strain_S12, *filename_Strain_S13, *filename_Strain_S23,
         *filename_Stress_S11, *filename_Stress_S22, *filename_Stress_S33,
         *filename_Stress_S12, *filename_Stress_S13, *filename_Stress_S23;
    FILE *fp_Strain_S11, *fp_Strain_S22, *fp_Strain_S33, *fp_Strain_S12, *fp_Strain_S13, *fp_Strain_S23,
         *fp_Stress_S11, *fp_Stress_S22, *fp_Stress_S33, *fp_Stress_S12, *fp_Stress_S13, *fp_Stress_S23;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    Strain = Make1DArray(6);
    Stress = Make1DArray(6);
    filename_Strain_S11 = Make1DString(500);
    filename_Strain_S22 = Make1DString(500);
    filename_Strain_S33 = Make1DString(500);
    filename_Strain_S12 = Make1DString(500);
    filename_Strain_S13 = Make1DString(500);
    filename_Strain_S23 = Make1DString(500);

    filename_Stress_S11 = Make1DString(500);
    filename_Stress_S22 = Make1DString(500);
    filename_Stress_S33 = Make1DString(500);
    filename_Stress_S12 = Make1DString(500);
    filename_Stress_S13 = Make1DString(500);
    filename_Stress_S23 = Make1DString(500);

    Version = IPN->Version;
    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }

    for(k=0;k<l;++k)
    {
        sprintf(filename_Strain_S11,"Version %d\\Elastic Field\\Strain\\z=%d\\S11.dat",Version,k-zc);
        sprintf(filename_Strain_S22,"Version %d\\Elastic Field\\Strain\\z=%d\\S22.dat",Version,k-zc);
        sprintf(filename_Strain_S33,"Version %d\\Elastic Field\\Strain\\z=%d\\S33.dat",Version,k-zc);
        sprintf(filename_Strain_S12,"Version %d\\Elastic Field\\Strain\\z=%d\\S12.dat",Version,k-zc);
        sprintf(filename_Strain_S13,"Version %d\\Elastic Field\\Strain\\z=%d\\S13.dat",Version,k-zc);
        sprintf(filename_Strain_S23,"Version %d\\Elastic Field\\Strain\\z=%d\\S23.dat",Version,k-zc);

        sprintf(filename_Stress_S11,"Version %d\\Elastic Field\\Stress\\z=%d\\S11.dat",Version,k-zc);
        sprintf(filename_Stress_S22,"Version %d\\Elastic Field\\Stress\\z=%d\\S22.dat",Version,k-zc);
        sprintf(filename_Stress_S33,"Version %d\\Elastic Field\\Stress\\z=%d\\S33.dat",Version,k-zc);
        sprintf(filename_Stress_S12,"Version %d\\Elastic Field\\Stress\\z=%d\\S12.dat",Version,k-zc);
        sprintf(filename_Stress_S13,"Version %d\\Elastic Field\\Stress\\z=%d\\S13.dat",Version,k-zc);
        sprintf(filename_Stress_S23,"Version %d\\Elastic Field\\Stress\\z=%d\\S23.dat",Version,k-zc);

        fp_Strain_S11=fopen(filename_Strain_S11,"w+");
        fp_Strain_S22=fopen(filename_Strain_S22,"w+");
        fp_Strain_S33=fopen(filename_Strain_S33,"w+");
        fp_Strain_S12=fopen(filename_Strain_S12,"w+");
        fp_Strain_S13=fopen(filename_Strain_S13,"w+");
        fp_Strain_S23=fopen(filename_Strain_S23,"w+");

        fp_Stress_S11=fopen(filename_Stress_S11,"w+");
        fp_Stress_S22=fopen(filename_Stress_S22,"w+");
        fp_Stress_S33=fopen(filename_Stress_S33,"w+");
        fp_Stress_S12=fopen(filename_Stress_S12,"w+");
        fp_Stress_S13=fopen(filename_Stress_S13,"w+");
        fp_Stress_S23=fopen(filename_Stress_S23,"w+");

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);

                    if(j<n-1)
                    {
                        fprintf(fp_Strain_S11,"%0.8f\t",Strain[0]);
                        fprintf(fp_Strain_S22,"%0.8f\t",Strain[1]);
                        fprintf(fp_Strain_S33,"%0.8f\t",Strain[2]);
                        fprintf(fp_Strain_S12,"%0.8f\t",Strain[3]);
                        fprintf(fp_Strain_S13,"%0.8f\t",Strain[4]);
                        fprintf(fp_Strain_S23,"%0.8f\t",Strain[5]);

                        fprintf(fp_Stress_S11,"%0.8f\t",Stress[0]);
                        fprintf(fp_Stress_S22,"%0.8f\t",Stress[1]);
                        fprintf(fp_Stress_S33,"%0.8f\t",Stress[2]);
                        fprintf(fp_Stress_S12,"%0.8f\t",Stress[3]);
                        fprintf(fp_Stress_S13,"%0.8f\t",Stress[4]);
                        fprintf(fp_Stress_S23,"%0.8f\t",Stress[5]);
                    }else
                    {
                        fprintf(fp_Strain_S11,"%0.8f",Strain[0]);
                        fprintf(fp_Strain_S22,"%0.8f",Strain[1]);
                        fprintf(fp_Strain_S33,"%0.8f",Strain[2]);
                        fprintf(fp_Strain_S12,"%0.8f",Strain[3]);
                        fprintf(fp_Strain_S13,"%0.8f",Strain[4]);
                        fprintf(fp_Strain_S23,"%0.8f",Strain[5]);

                        fprintf(fp_Stress_S11,"%0.8f",Stress[0]);
                        fprintf(fp_Stress_S22,"%0.8f",Stress[1]);
                        fprintf(fp_Stress_S33,"%0.8f",Stress[2]);
                        fprintf(fp_Stress_S12,"%0.8f",Stress[3]);
                        fprintf(fp_Stress_S13,"%0.8f",Stress[4]);
                        fprintf(fp_Stress_S23,"%0.8f",Stress[5]);
                    }
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                    StressToStrain(Stress,Strain);

                    if(j<n-1)
                    {
                        fprintf(fp_Strain_S11,"%0.8f\t",Strain[0]);
                        fprintf(fp_Strain_S22,"%0.8f\t",Strain[1]);
                        fprintf(fp_Strain_S33,"%0.8f\t",Strain[2]);
                        fprintf(fp_Strain_S12,"%0.8f\t",Strain[3]);
                        fprintf(fp_Strain_S13,"%0.8f\t",Strain[4]);
                        fprintf(fp_Strain_S23,"%0.8f\t",Strain[5]);

                        fprintf(fp_Stress_S11,"%0.8f\t",Stress[0]);
                        fprintf(fp_Stress_S22,"%0.8f\t",Stress[1]);
                        fprintf(fp_Stress_S33,"%0.8f\t",Stress[2]);
                        fprintf(fp_Stress_S12,"%0.8f\t",Stress[3]);
                        fprintf(fp_Stress_S13,"%0.8f\t",Stress[4]);
                        fprintf(fp_Stress_S23,"%0.8f\t",Stress[5]);
                    }else
                    {
                        fprintf(fp_Strain_S11,"%0.8f",Strain[0]);
                        fprintf(fp_Strain_S22,"%0.8f",Strain[1]);
                        fprintf(fp_Strain_S33,"%0.8f",Strain[2]);
                        fprintf(fp_Strain_S12,"%0.8f",Strain[3]);
                        fprintf(fp_Strain_S13,"%0.8f",Strain[4]);
                        fprintf(fp_Strain_S23,"%0.8f",Strain[5]);

                        fprintf(fp_Stress_S11,"%0.8f",Stress[0]);
                        fprintf(fp_Stress_S22,"%0.8f",Stress[1]);
                        fprintf(fp_Stress_S33,"%0.8f",Stress[2]);
                        fprintf(fp_Stress_S12,"%0.8f",Stress[3]);
                        fprintf(fp_Stress_S13,"%0.8f",Stress[4]);
                        fprintf(fp_Stress_S23,"%0.8f",Stress[5]);
                    }
                }
            }
            fprintf(fp_Strain_S11,"\n");
            fprintf(fp_Strain_S22,"\n");
            fprintf(fp_Strain_S33,"\n");
            fprintf(fp_Strain_S12,"\n");
            fprintf(fp_Strain_S13,"\n");
            fprintf(fp_Strain_S23,"\n");

            fprintf(fp_Stress_S11,"\n");
            fprintf(fp_Stress_S22,"\n");
            fprintf(fp_Stress_S33,"\n");
            fprintf(fp_Stress_S12,"\n");
            fprintf(fp_Stress_S13,"\n");
            fprintf(fp_Stress_S23,"\n");
        }
        fclose(fp_Strain_S11);
        fclose(fp_Strain_S22);
        fclose(fp_Strain_S33);
        fclose(fp_Strain_S12);
        fclose(fp_Strain_S13);
        fclose(fp_Strain_S23);

        fclose(fp_Stress_S11);
        fclose(fp_Stress_S22);
        fclose(fp_Stress_S33);
        fclose(fp_Stress_S12);
        fclose(fp_Stress_S13);
        fclose(fp_Stress_S23);
    }
//////////////////////////////////////////////////////////////////////////////////////
    for(k=0;k<l;++k)
    {
        sprintf(filename_Strain_S11,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S11.dat",Version,k-zc);
        sprintf(filename_Strain_S22,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S22.dat",Version,k-zc);
        sprintf(filename_Strain_S33,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S33.dat",Version,k-zc);
        sprintf(filename_Strain_S12,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S12.dat",Version,k-zc);
        sprintf(filename_Strain_S13,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S13.dat",Version,k-zc);
        sprintf(filename_Strain_S23,"Version %d\\For Printing\\Elastic Field\\Strain\\z=%d\\S23.dat",Version,k-zc);

        sprintf(filename_Stress_S11,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S11.dat",Version,k-zc);
        sprintf(filename_Stress_S22,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S22.dat",Version,k-zc);
        sprintf(filename_Stress_S33,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S33.dat",Version,k-zc);
        sprintf(filename_Stress_S12,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S12.dat",Version,k-zc);
        sprintf(filename_Stress_S13,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S13.dat",Version,k-zc);
        sprintf(filename_Stress_S23,"Version %d\\For Printing\\Elastic Field\\Stress\\z=%d\\S23.dat",Version,k-zc);

        fp_Strain_S11=fopen(filename_Strain_S11,"w+");
        fp_Strain_S22=fopen(filename_Strain_S22,"w+");
        fp_Strain_S33=fopen(filename_Strain_S33,"w+");
        fp_Strain_S12=fopen(filename_Strain_S12,"w+");
        fp_Strain_S13=fopen(filename_Strain_S13,"w+");
        fp_Strain_S23=fopen(filename_Strain_S23,"w+");

        fp_Stress_S11=fopen(filename_Stress_S11,"w+");
        fp_Stress_S22=fopen(filename_Stress_S22,"w+");
        fp_Stress_S33=fopen(filename_Stress_S33,"w+");
        fp_Stress_S12=fopen(filename_Stress_S12,"w+");
        fp_Stress_S13=fopen(filename_Stress_S13,"w+");
        fp_Stress_S23=fopen(filename_Stress_S23,"w+");

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);

                    if(j<n-1)
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Strain_S11,"%0.8f\t",Strain[0]);
                            fprintf(fp_Strain_S22,"%0.8f\t",Strain[1]);
                            fprintf(fp_Strain_S33,"%0.8f\t",Strain[2]);
                            fprintf(fp_Strain_S12,"%0.8f\t",Strain[3]);
                            fprintf(fp_Strain_S13,"%0.8f\t",Strain[4]);
                            fprintf(fp_Strain_S23,"%0.8f\t",Strain[5]);

                            fprintf(fp_Stress_S11,"%0.8f\t",Stress[0]);
                            fprintf(fp_Stress_S22,"%0.8f\t",Stress[1]);
                            fprintf(fp_Stress_S33,"%0.8f\t",Stress[2]);
                            fprintf(fp_Stress_S12,"%0.8f\t",Stress[3]);
                            fprintf(fp_Stress_S13,"%0.8f\t",Stress[4]);
                            fprintf(fp_Stress_S23,"%0.8f\t",Stress[5]);
                        }else
                        {
                            fprintf(fp_Strain_S11,"\t");
                            fprintf(fp_Strain_S22,"\t");
                            fprintf(fp_Strain_S33,"\t");
                            fprintf(fp_Strain_S12,"\t");
                            fprintf(fp_Strain_S13,"\t");
                            fprintf(fp_Strain_S23,"\t");

                            fprintf(fp_Stress_S11,"\t");
                            fprintf(fp_Stress_S22,"\t");
                            fprintf(fp_Stress_S33,"\t");
                            fprintf(fp_Stress_S12,"\t");
                            fprintf(fp_Stress_S13,"\t");
                            fprintf(fp_Stress_S23,"\t");
                        }
                    }else
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Strain_S11,"%0.8f",Strain[0]);
                            fprintf(fp_Strain_S22,"%0.8f",Strain[1]);
                            fprintf(fp_Strain_S33,"%0.8f",Strain[2]);
                            fprintf(fp_Strain_S12,"%0.8f",Strain[3]);
                            fprintf(fp_Strain_S13,"%0.8f",Strain[4]);
                            fprintf(fp_Strain_S23,"%0.8f",Strain[5]);

                            fprintf(fp_Stress_S11,"%0.8f",Stress[0]);
                            fprintf(fp_Stress_S22,"%0.8f",Stress[1]);
                            fprintf(fp_Stress_S33,"%0.8f",Stress[2]);
                            fprintf(fp_Stress_S12,"%0.8f",Stress[3]);
                            fprintf(fp_Stress_S13,"%0.8f",Stress[4]);
                            fprintf(fp_Stress_S23,"%0.8f",Stress[5]);
                        }
                    }
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                    StressToStrain(Stress,Strain);

                    if(j<n-1)
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Strain_S11,"%0.8f\t",Strain[0]);
                            fprintf(fp_Strain_S22,"%0.8f\t",Strain[1]);
                            fprintf(fp_Strain_S33,"%0.8f\t",Strain[2]);
                            fprintf(fp_Strain_S12,"%0.8f\t",Strain[3]);
                            fprintf(fp_Strain_S13,"%0.8f\t",Strain[4]);
                            fprintf(fp_Strain_S23,"%0.8f\t",Strain[5]);

                            fprintf(fp_Stress_S11,"%0.8f\t",Stress[0]);
                            fprintf(fp_Stress_S22,"%0.8f\t",Stress[1]);
                            fprintf(fp_Stress_S33,"%0.8f\t",Stress[2]);
                            fprintf(fp_Stress_S12,"%0.8f\t",Stress[3]);
                            fprintf(fp_Stress_S13,"%0.8f\t",Stress[4]);
                            fprintf(fp_Stress_S23,"%0.8f\t",Stress[5]);
                        }else
                        {
                            fprintf(fp_Strain_S11,"\t");
                            fprintf(fp_Strain_S22,"\t");
                            fprintf(fp_Strain_S33,"\t");
                            fprintf(fp_Strain_S12,"\t");
                            fprintf(fp_Strain_S13,"\t");
                            fprintf(fp_Strain_S23,"\t");

                            fprintf(fp_Stress_S11,"\t");
                            fprintf(fp_Stress_S22,"\t");
                            fprintf(fp_Stress_S33,"\t");
                            fprintf(fp_Stress_S12,"\t");
                            fprintf(fp_Stress_S13,"\t");
                            fprintf(fp_Stress_S23,"\t");
                        }
                    }else
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Strain_S11,"%0.8f",Strain[0]);
                            fprintf(fp_Strain_S22,"%0.8f",Strain[1]);
                            fprintf(fp_Strain_S33,"%0.8f",Strain[2]);
                            fprintf(fp_Strain_S12,"%0.8f",Strain[3]);
                            fprintf(fp_Strain_S13,"%0.8f",Strain[4]);
                            fprintf(fp_Strain_S23,"%0.8f",Strain[5]);

                            fprintf(fp_Stress_S11,"%0.8f",Stress[0]);
                            fprintf(fp_Stress_S22,"%0.8f",Stress[1]);
                            fprintf(fp_Stress_S33,"%0.8f",Stress[2]);
                            fprintf(fp_Stress_S12,"%0.8f",Stress[3]);
                            fprintf(fp_Stress_S13,"%0.8f",Stress[4]);
                            fprintf(fp_Stress_S23,"%0.8f",Stress[5]);
                        }
                    }
                }
            }
            fprintf(fp_Strain_S11,"\n");
            fprintf(fp_Strain_S22,"\n");
            fprintf(fp_Strain_S33,"\n");
            fprintf(fp_Strain_S12,"\n");
            fprintf(fp_Strain_S13,"\n");
            fprintf(fp_Strain_S23,"\n");

            fprintf(fp_Stress_S11,"\n");
            fprintf(fp_Stress_S22,"\n");
            fprintf(fp_Stress_S33,"\n");
            fprintf(fp_Stress_S12,"\n");
            fprintf(fp_Stress_S13,"\n");
            fprintf(fp_Stress_S23,"\n");
        }
        fclose(fp_Strain_S11);
        fclose(fp_Strain_S22);
        fclose(fp_Strain_S33);
        fclose(fp_Strain_S12);
        fclose(fp_Strain_S13);
        fclose(fp_Strain_S23);

        fclose(fp_Stress_S11);
        fclose(fp_Stress_S22);
        fclose(fp_Stress_S33);
        fclose(fp_Stress_S12);
        fclose(fp_Stress_S13);
        fclose(fp_Stress_S23);
    }
 /////////////////////////////////////////////////////
    free(filename_Strain_S11);
    free(filename_Strain_S22);
    free(filename_Strain_S33);
    free(filename_Strain_S12);
    free(filename_Strain_S13);
    free(filename_Strain_S23);

    free(filename_Stress_S11);
    free(filename_Stress_S22);
    free(filename_Stress_S33);
    free(filename_Stress_S12);
    free(filename_Stress_S13);
    free(filename_Stress_S23);

    free(Strain);
    free(Stress);

    return 1;
}

int Write_NMD_Magnetization(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                            double *MxV, double *MyV, double *MzV, Boundary_Shape *BS)
{
    int i, j, k, v, m, n, l, zc, Version;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV;

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

    sprintf(filename_MxV,"Version %d\\M3D\\MxV.dat",Version);
    sprintf(filename_MyV,"Version %d\\M3D\\MyV.dat",Version);
    sprintf(filename_MzV,"Version %d\\M3D\\MzV.dat",Version);
    fp_MxV = fopen(filename_MxV,"w+");
    fp_MyV = fopen(filename_MyV,"w+");
    fp_MzV = fopen(filename_MzV,"w+");

    for(v=0;v<BS->Virtual_Points;++v)
    {
        fprintf(fp_MxV,"%0.8f\t",MxV[v]);
        fprintf(fp_MyV,"%0.8f\t",MyV[v]);
        fprintf(fp_MzV,"%0.8f\t",MzV[v]);
    }
    fclose(fp_MxV);
    fclose(fp_MyV);
    fclose(fp_MzV);

    for(k=0;k<l;k++)
    {
        sprintf(filename_Mx,"Version %d\\M3D\\z=%d\\Mx.dat",Version,k-zc);
        sprintf(filename_My,"Version %d\\M3D\\z=%d\\My.dat",Version,k-zc);
        sprintf(filename_Mz,"Version %d\\M3D\\z=%d\\Mz.dat",Version,k-zc);
        fp_Mx = fopen(filename_Mx,"w+");
        fp_My = fopen(filename_My,"w+");
        fp_Mz = fopen(filename_Mz,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    fprintf(fp_Mx,"%0.8f\t",Mx[i][j][k]);
                    fprintf(fp_My,"%0.8f\t",My[i][j][k]);
                    fprintf(fp_Mz,"%0.8f\t",Mz[i][j][k]);
                }else
                {
                    fprintf(fp_Mx,"%0.8f",Mx[i][j][k]);
                    fprintf(fp_My,"%0.8f",My[i][j][k]);
                    fprintf(fp_Mz,"%0.8f",Mz[i][j][k]);
                }
            }
            if(i<m-1)
            {
                fprintf(fp_Mx,"\n");
                fprintf(fp_My,"\n");
                fprintf(fp_Mz,"\n");
            }
        }
        fclose(fp_Mx);
        fclose(fp_My);
        fclose(fp_Mz);
    }
//////////////////////////////////////////////////////////////////////////////////////
    for(k=0;k<l;k++)
    {
        sprintf(filename_Mx,"Version %d\\For Printing\\z=%d\\Mx.dat",Version,k-zc);
        sprintf(filename_My,"Version %d\\For Printing\\z=%d\\My.dat",Version,k-zc);
        sprintf(filename_Mz,"Version %d\\For Printing\\z=%d\\Mz.dat",Version,k-zc);
        fp_Mx = fopen(filename_Mx,"w+");
        fp_My = fopen(filename_My,"w+");
        fp_Mz = fopen(filename_Mz,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        fprintf(fp_Mx,"%0.8f\t",Mx[i][j][k]);
                        fprintf(fp_My,"%0.8f\t",My[i][j][k]);
                        fprintf(fp_Mz,"%0.8f\t",Mz[i][j][k]);
                    }else
                    {
                        fprintf(fp_Mx,"\t");
                        fprintf(fp_My,"\t");
                        fprintf(fp_Mz,"\t");
                    }
                }else
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        fprintf(fp_Mx,"%0.8f",Mx[i][j][k]);
                        fprintf(fp_My,"%0.8f",My[i][j][k]);
                        fprintf(fp_Mz,"%0.8f",Mz[i][j][k]);
                    }
                }
            }

            if(i<m-1)
            {
                fprintf(fp_Mx,"\n");
                fprintf(fp_My,"\n");
                fprintf(fp_Mz,"\n");
            }
        }
        fclose(fp_Mx);
        fclose(fp_My);
        fclose(fp_Mz);
    }
//////////////////////////////////////////////////////////////////////////////
    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    return 1;
}

int Write_NMD_Magnetization_SKN(Input_Parameter_NMD *IPN, double ***Mx, double ***My, double ***Mz,
                                double *MxV, double *MyV, double *MzV, Boundary_Shape *BS, char *dir_SKN)
{
    int i, j, k, v, m, n, l, zc, Version;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV, *String_Temp;

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

    sprintf(String_Temp,"Version %d\\M3D\\MxV.dat",Version);
    sprintf(filename_MxV, "%s\\%s", dir_SKN, String_Temp);

    sprintf(String_Temp,"Version %d\\M3D\\MyV.dat",Version);
    sprintf(filename_MyV, "%s\\%s", dir_SKN, String_Temp);

    sprintf(String_Temp,"Version %d\\M3D\\MzV.dat",Version);
    sprintf(filename_MzV, "%s\\%s", dir_SKN, String_Temp);

    fp_MxV = fopen(filename_MxV,"w+");
    fp_MyV = fopen(filename_MyV,"w+");
    fp_MzV = fopen(filename_MzV,"w+");

    for(v=0;v<BS->Virtual_Points;++v)
    {
        fprintf(fp_MxV,"%0.8f\t",MxV[v]);
        fprintf(fp_MyV,"%0.8f\t",MyV[v]);
        fprintf(fp_MzV, "%0.8f\t", MzV[v]);
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

        fp_Mx = fopen(filename_Mx,"w+");
        fp_My = fopen(filename_My,"w+");
        fp_Mz = fopen(filename_Mz,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    fprintf(fp_Mx,"%0.8f\t",Mx[i][j][k]);
                    fprintf(fp_My,"%0.8f\t",My[i][j][k]);
                    fprintf(fp_Mz,"%0.8f\t",Mz[i][j][k]);
                }else
                {
                    fprintf(fp_Mx,"%0.8f",Mx[i][j][k]);
                    fprintf(fp_My,"%0.8f",My[i][j][k]);
                    fprintf(fp_Mz,"%0.8f",Mz[i][j][k]);
                }
            }
            if(i<m-1)
            {
                fprintf(fp_Mx,"\n");
                fprintf(fp_My,"\n");
                fprintf(fp_Mz,"\n");
            }
        }
        fclose(fp_Mx);
        fclose(fp_My);
        fclose(fp_Mz);
    }
//////////////////////////////////////////////////////////////////////////////////////
    for(k=0;k<l;k++)
    {        
        sprintf(String_Temp,"Version %d\\For Printing\\z=%d\\Mx.dat",Version,k-zc);
        sprintf(filename_Mx, "%s\\%s", dir_SKN, String_Temp);

        sprintf(String_Temp,"Version %d\\For Printing\\z=%d\\My.dat",Version,k-zc);
        sprintf(filename_My, "%s\\%s", dir_SKN, String_Temp);

        sprintf(String_Temp,"Version %d\\For Printing\\z=%d\\Mz.dat",Version,k-zc);
        sprintf(filename_Mz, "%s\\%s", dir_SKN, String_Temp);

        fp_Mx = fopen(filename_Mx,"w+");
        fp_My = fopen(filename_My,"w+");
        fp_Mz = fopen(filename_Mz,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        fprintf(fp_Mx,"%0.8f\t",Mx[i][j][k]);
                        fprintf(fp_My,"%0.8f\t",My[i][j][k]);
                        fprintf(fp_Mz,"%0.8f\t",Mz[i][j][k]);
                    }else
                    {
                        fprintf(fp_Mx,"\t");
                        fprintf(fp_My,"\t");
                        fprintf(fp_Mz,"\t");
                    }
                }else
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        fprintf(fp_Mx,"%0.8f",Mx[i][j][k]);
                        fprintf(fp_My,"%0.8f",My[i][j][k]);
                        fprintf(fp_Mz,"%0.8f",Mz[i][j][k]);
                    }
                }
            }

            if(i<m-1)
            {
                fprintf(fp_Mx,"\n");
                fprintf(fp_My,"\n");
                fprintf(fp_Mz,"\n");
            }
        }
        fclose(fp_Mx);
        fclose(fp_My);
        fclose(fp_Mz);
    }
//////////////////////////////////////////////////////////////////////////////
    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    free(String_Temp);
    return 1;
}

int Write_NMD_Energy(Input_Parameter_NMD *IPN, double ***w)
{
    int i, j, k, m, n, l, zc, Version;
    FILE *fp;
    char *filename;

    filename  = Make1DString(500);
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

    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\M3D\\z=%d\\w.dat",Version,k-zc);
        fp = fopen(filename,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    fprintf(fp,"%0.8f\t",w[i][j][k]);
                }else
                {
                    fprintf(fp,"%0.8f",w[i][j][k]);
                }
            }
            fprintf(fp,"\n");
        }
        fclose(fp);
    }
    free(filename);
    return 1;
}

int Write_NMD_SKN_D(Input_Parameter_NMD *IPN, double ***SKN_D)
{
    int i, j, k, m, n, l, zc, Version;
    FILE *fp, *fp_SKN_Print;
    char *filename, *filename_SKN_Print;
    Boundary_Shape BS;
    BS.xn = IPN->xn;    BS.yn = IPN->yn;
    BS.zn = IPN->zn;

    filename  = Make1DString(500);
    filename_SKN_Print = Make1DString(500);
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

    BS.Boundary=Make3DArrayinteger(m,n,l);
    BS.LocalEnergy=Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,&BS);

    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\M3D\\z=%d\\SKN_D.dat",Version,k-zc);
        sprintf(filename_SKN_Print,"Version %d\\M3D\\z=%d\\SKN_Print.dat",Version,k-zc);
        fp = fopen(filename,"w+");
        fp_SKN_Print = fopen(filename_SKN_Print,"w+");

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if(j<n-1)
                {
                    fprintf(fp,"%0.8f\t",SKN_D[i][j][k]);
                }else
                {
                    fprintf(fp,"%0.8f",SKN_D[i][j][k]);
                }

                if(BS.Boundary[i][j][k]!=0)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_SKN_Print,"%0.8f\t",SKN_D[i][j][k]);
                    }else
                    {
                        fprintf(fp_SKN_Print,"%0.8f",SKN_D[i][j][k]);
                    }
                }else
                {
                    if(j<n-1)
                    {
                        fprintf(fp_SKN_Print,"\t");
                    }
                }
            }
            fprintf(fp,"\n");
            fprintf(fp_SKN_Print,"\n");
        }
        fclose(fp);
        fclose(fp_SKN_Print);
    }
    free(filename);
    free(filename_SKN_Print);
    free3DArrayinteger(BS.Boundary,m,n);
    free4DArrayinteger(BS.LocalEnergy,m,n,l);
    return 1;
}

int Write_NMD_Result(Input_Parameter_NMD *IPN)
{
    char *filename;
    FILE *fp;
    filename = Make1DString(500);
    sprintf(filename,"Version %d\\Result.dat",IPN->Version);
    fp=fopen(filename,"w+");
    fprintf(fp,"hex\t\t t\n");
    fprintf(fp,"%g\t\t %g\n", IPN->h, IPN->t);
    fprintf(fp,"Magnetic Field Mode=%d\n",IPN->MagneticFieldMode);
    if(IPN->MagneticFieldMode==0)
    {
        fprintf(fp,"Homogeneous steady magnetic field\n");
    }else if(IPN->MagneticFieldMode==1)
    {
        fprintf(fp,"b=%gsin(2pi*t/%d)\n", IPN->h, IPN->MagneticFieldPeriod);
    }

    fprintf(fp,"xn\t\t yn\t\t zn\n");
    fprintf(fp,"%d\t\t %d\t\t %d\n", IPN->xn, IPN->yn, IPN->zn);
    fprintf(fp,"Boundary Shape Type=%d\n",IPN->Boundary_Shape_Type);
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"W_Average=%0.8f\n",IPN->W);
    fprintf(fp,"\n");
    fprintf(fp,"SKN=%0.8f\t\t SKN Area=%0.8f\n",IPN->SKN,IPN->SKN_Area);
    fprintf(fp,"\n");
    fprintf(fp,"Time=%0.8f\n",IPN->Time);
    fprintf(fp,"\n");
    if(IPN->CalculateMode==2||IPN->CalculateMode==3)
    {
        fprintf(fp,"dM/dt=%0.8f\n",IPN->DerivativeM);
    }
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"Elastic Mode=%d\n",IPN->ElasticFieldMode);
    if(IPN->ElasticFieldMode==1)
    {
        fprintf(fp,"Uniaxial Strain\n");
    }else if(IPN->ElasticFieldMode==2)
    {
        fprintf(fp,"Uniaxial Stress\n");
    }else
    {
        fprintf(fp,"No Uniaxial Field\n");
    }


    fprintf(fp,"S11\t\t %0.8f\n",IPN->CST_Uniaxial[0]);
    fprintf(fp,"S22\t\t %0.8f\n",IPN->CST_Uniaxial[1]);
    fprintf(fp,"S33\t\t %0.8f\n",IPN->CST_Uniaxial[2]);
    fprintf(fp,"S12\t\t %0.8f\n",IPN->CST_Uniaxial[3]);
    fprintf(fp,"S13\t\t %0.8f\n",IPN->CST_Uniaxial[4]);
    fprintf(fp,"S23\t\t %0.8f\n",IPN->CST_Uniaxial[5]);
    fprintf(fp,"\n");
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"a=%g\t b=%g\t FC=%g\n",IPN->a,IPN->b,IPN->FC);
    fprintf(fp,"----------------------------------\n");
    fprintf(fp,"dx\t\t dy\t\t dz\n");
    fprintf(fp,"%0.8f\t %0.8f\t %0.8f\n", IPN->dxy[0], IPN->dxy[1], IPN->dxy[2]);
    fclose(fp);

    free(filename);
    return 1;
}

int WriteVectorM_3D_Monolayer(int m, int n, double ***Mx, double ***My, double ***Mz)
{
    int i, j;
    FILE *fp1, *fp2, *fp3;
    fp1=fopen("Mx.dat","w+");
    fp2=fopen("My.dat","w+");
    fp3=fopen("Mz.dat","w+");

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            fprintf(fp1,"%0.8f\t",Mx[i][j][0]);
            fprintf(fp2,"%0.8f\t",My[i][j][0]);
            fprintf(fp3,"%0.8f\t",Mz[i][j][0]);
        }
        fprintf(fp1,"\n");
        fprintf(fp2,"\n");
        fprintf(fp3,"\n");
    }
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);

    return 1;
}

int Simplify_Matrix_3D(Input_Parameter_NMD *IPN, Boundary_Shape *BS)
{
    FILE *fp_Mx,*fp_My, *fp_Mz, *fp_x, *fp_y, *fp_z;
    int k, i, j, m, n, r, l, zc;
    double ***Mx, ***My, ***Mz, *MxV, *MyV, *MzV, x, y, z;
    char *filename, *filename_Mx, *filename_My, *filename_Mz, *filename_x, *filename_y, *filename_z;

    filename     = Make1DString(500);
    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_x   = Make1DString(500);
    filename_y   = Make1DString(500);
    filename_z   = Make1DString(500);

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

    r=m/24;

    Mx=Make3DArray(m,n,l);
    My=Make3DArray(m,n,l);
    Mz=Make3DArray(m,n,l);
    MxV=Make1DArray(BS->Virtual_Points);
    MyV=Make1DArray(BS->Virtual_Points);
    MzV=Make1DArray(BS->Virtual_Points);

    InitialM_3D(IPN,Mx,My,Mz,MxV,MyV,MzV,BS);
    sprintf(filename,"SIM");
    mkdir(filename);
    sprintf(filename_Mx,"SIM\\Mx_Vector.dat");
    sprintf(filename_My,"SIM\\My_Vector.dat");
    sprintf(filename_Mz,"SIM\\Mz_Vector.dat");
    sprintf(filename_x,"SIM\\x.dat");
    sprintf(filename_y,"SIM\\y.dat");
    sprintf(filename_z,"SIM\\z.dat");

    fp_x=fopen(filename_x,"w+");
    fp_y=fopen(filename_y,"w+");
    fp_z=fopen(filename_z,"w+");

    fp_Mx=fopen(filename_Mx,"w+");
    fp_My=fopen(filename_My,"w+");
    fp_Mz=fopen(filename_Mz,"w+");

    for(k=0;k<l;k++)
    {
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if((i-IPN->xn-1)%r==0 && (j-IPN->yn-1)%r==0)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        x=(i-IPN->xn-1)*IPN->dxy[0];
                        y=(j-IPN->yn-1)*IPN->dxy[1];
                        z=(k-zc)*IPN->dxy[2];
                        fprintf(fp_x,"%lf\n",x);
                        fprintf(fp_y,"%lf\n",y);
                        fprintf(fp_z,"%lf\n",z);

                        fprintf(fp_Mx,"%lf\n",Mx[i][j][k]);
                        fprintf(fp_My,"%lf\n",My[i][j][k]);
                        fprintf(fp_Mz,"%lf\n",Mz[i][j][k]);
                    }
                }
            }
        }
    }


    fclose(fp_x);
    fclose(fp_y);
    fclose(fp_z);
    fclose(fp_Mx);
    fclose(fp_My);
    fclose(fp_Mz);

    free3DArray(Mx,m,n);
    free3DArray(My,m,n);
    free3DArray(Mz,m,n);
    free(MxV);
    free(MyV);
    free(MzV);
    free(filename);
    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_x);
    free(filename_y);
    free(filename_z);

    return 1;
}

int Write_Inplane_M(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ***Mx, double ***My)
{
    FILE *fp_M_Inplane;
    int k, i, j, m, n, l, zc, Inplane_Density, D;
    double x, y, Mx_Offset, My_Offset, x_max, x_min, y_max, y_min;
    double Mx0, My0, M_Module, M_Appropriate, Ratio, M_min;
    char *filename;

    Inplane_Density = IPN->Inplane_Density;
    filename     = Make1DString(500);

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

    if(Inplane_Density == 0)
    {
        Inplane_Density = 24;
    }
    D = m/Inplane_Density;

    if(D<=1)
    {
        D=1;
    }

    M_Appropriate = D*fabs(IPN->dxy[0])*0.8;
    M_min         = D*fabs(IPN->dxy[0])*0.02;

    x_max = fabs(IPN->dxy[0])*(m-1)/2;   x_min = -fabs(IPN->dxy[0])*(m-1)/2;
    y_max = fabs(IPN->dxy[1])*(n-1)/2;   y_min = -fabs(IPN->dxy[1])*(n-1)/2;
    for(k=0;k<l;k++)
    {
        sprintf(filename,"Version %d\\In-Plane Components\\z=%d\\M_InPlane.dat",IPN->Version,k-zc);
        fp_M_Inplane = fopen(filename,"w+");
        fprintf(fp_M_Inplane,"y\tx\tMy_Offset\tMx_Offset\tMy\tMx\n");
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if((i-IPN->xn-1)%D==0 && (j-IPN->yn-1)%D==0)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        x=(i-IPN->xn-1)*IPN->dxy[0];
                        y=(j-IPN->yn-1)*IPN->dxy[1];
                        Mx0 = Mx[i][j][k];
                        My0 = My[i][j][k];

                        M_Module = sqrt(Mx0*Mx0+My0*My0);
                        if(M_Module>M_Appropriate)
                        {
                            Ratio = M_Appropriate/M_Module;
                            Mx0 = Mx0*Ratio;
                            My0 = My0*Ratio;
                        }else if(M_Module<M_min)
                        {
                            Mx0 = 0;
                            My0 = 0;
                        }

                        Mx_Offset = Mx0+x;
                        My_Offset = My0+y;
                        if(Mx_Offset>=x_max || Mx_Offset<=x_min || My_Offset>=y_max || My_Offset<=y_min)
                        {
                            Mx_Offset = x;
                            My_Offset = y;

                        }

                        fprintf(fp_M_Inplane,"%0.2f\t%0.2f\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",y,x,My_Offset,Mx_Offset,My[i][j][k],Mx[i][j][k]);

                    }
                }
            }
        }
        fclose(fp_M_Inplane);
    }

    free(filename);

    return 1;
}

int Write_Inplane_M_SKN(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ***Mx, double ***My, char *dir_SKN)
{
    FILE *fp_M_Inplane;
    int k, i, j, m, n, l, zc, Inplane_Density, D;
    double x, y, Mx_Offset, My_Offset, x_max, x_min, y_max, y_min;
    double Mx0, My0, M_Module, M_Appropriate, Ratio, M_min;
    char *filename, *String_Temp;

    Inplane_Density = IPN->Inplane_Density;
    filename     = Make1DString(500);
    String_Temp  = Make1DString(500);

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

    if(Inplane_Density == 0)
    {
        Inplane_Density = 24;
    }
    D = m/Inplane_Density;

    if(D<=1)
    {
        D=1;
    }

    M_Appropriate = D*fabs(IPN->dxy[0])*0.8;
    M_min         = D*fabs(IPN->dxy[0])*0.02;

    x_max = fabs(IPN->dxy[0])*(m-1)/2;   x_min = -fabs(IPN->dxy[0])*(m-1)/2;
    y_max = fabs(IPN->dxy[1])*(n-1)/2;   y_min = -fabs(IPN->dxy[1])*(n-1)/2;
    for(k=0;k<l;k++)
    {
        sprintf(String_Temp,"Version %d\\In-Plane Components\\z=%d\\M_InPlane.dat",IPN->Version,k-zc);
        sprintf(filename, "%s\\%s", dir_SKN, String_Temp);
        fp_M_Inplane = fopen(filename,"w+");
        fprintf(fp_M_Inplane,"y\tx\tMy_Offset\tMx_Offset\tMy\tMx\n");
        for(i=0;i<m;i++)
        {
            for(j=0;j<n;j++)
            {
                if((i-IPN->xn-1)%D==0 && (j-IPN->yn-1)%D==0)
                {
                    if(BS->Boundary[i][j][k]!=0)
                    {
                        x=(i-IPN->xn-1)*IPN->dxy[0];
                        y=(j-IPN->yn-1)*IPN->dxy[1];
                        Mx0 = Mx[i][j][k];
                        My0 = My[i][j][k];

                        M_Module = sqrt(Mx0*Mx0+My0*My0);
                        if(M_Module>M_Appropriate)
                        {
                            Ratio = M_Appropriate/M_Module;
                            Mx0 = Mx0*Ratio;
                            My0 = My0*Ratio;
                        }else if(M_Module<M_min)
                        {
                            Mx0 = 0;
                            My0 = 0;
                        }

                        Mx_Offset = Mx0+x;
                        My_Offset = My0+y;
                        if(Mx_Offset>=x_max || Mx_Offset<=x_min || My_Offset>=y_max || My_Offset<=y_min)
                        {
                            Mx_Offset = x;
                            My_Offset = y;

                        }

                        fprintf(fp_M_Inplane,"%0.2f\t%0.2f\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",y,x,My_Offset,Mx_Offset,My[i][j][k],Mx[i][j][k]);

                    }
                }
            }
        }
        fclose(fp_M_Inplane);
    }

    free(filename);
    free(String_Temp);

    return 1;
}

int Write_Judge_Phase_Curve(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****CST)
{
    FILE *fp_Judge_Phase_Curve, *fp_Stress_Difference, *fp_Stress_Bulk, *fp_Stress_C, *fp_Curve_1D;
    int k, i, j, m, n, l, zc, r, **Curve_Array;
    double *Stress, *Strain, Stress_Bulk, Stress_Critical;
    double **Stress_Difference, Ratio;
    char *filename, *filename_Difference, *filename_Stress_Bulk, *filename_Stress_Critical;
    char *filename_Curve_1D;

    filename                 = Make1DString(500);
    filename_Difference      = Make1DString(500);
    filename_Stress_Bulk     = Make1DString(500);
    filename_Stress_Critical = Make1DString(500);
    filename_Curve_1D        = Make1DString(500);

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
    Stress_Critical = (0.9-IPN->t)/3/q_Stress;

    sprintf(filename_Stress_Critical,"Version %d\\For Printing\\Stress Critical.dat",IPN->Version);
    fp_Stress_C = fopen(filename_Stress_Critical,"w+");
    fprintf(fp_Stress_C,"%0.8f\n",Stress_Critical);
    fclose(fp_Stress_C);

    Curve_Array = Make2DArrayinteger(m,n);
    Stress_Difference = Make2DArray(m,n);
    Stress = Make1DArray(6);
    Strain = Make1DArray(6);
    for(k=0;k<l;++k)
    {
        sprintf(filename,"Version %d\\For Printing\\z=%d\\Judge Phase Curve.dat",IPN->Version,k-zc);
        sprintf(filename_Difference,"Version %d\\For Printing\\z=%d\\Stress Difference.dat",IPN->Version,k-zc);
        sprintf(filename_Stress_Bulk,"Version %d\\For Printing\\z=%d\\Bulk Stress.dat",IPN->Version,k-zc);
        sprintf(filename_Curve_1D,"Version %d\\For Printing\\z=%d\\Judge Phase Curve 1D.dat",IPN->Version,k-zc);
        fp_Judge_Phase_Curve = fopen(filename,"w+");
        fp_Stress_Difference = fopen(filename_Difference,"w+");
        fp_Stress_Bulk       = fopen(filename_Stress_Bulk,"w+");
        fp_Curve_1D          = fopen(filename_Curve_1D,"w+");

        fprintf(fp_Curve_1D,"j\ti\n");

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                }
                Stress_Bulk = (Stress[0]+Stress[1]+Stress[2])/3;
                if(BS->Boundary[i][j][k]!=0)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_Stress_Bulk,"%0.8f\t",Stress_Bulk);
                    }else
                    {
                        fprintf(fp_Stress_Bulk,"%0.8f",Stress_Bulk);
                    }
                }else
                {
                    if(j<n-1)
                    {
                        fprintf(fp_Stress_Bulk,"\t");
                    }
                }
                Ratio = fabs((Stress_Bulk - Stress_Critical)/Stress_Critical);
                Stress_Difference[i][j] = Stress_Bulk-Stress_Critical;

                if(Ratio<=0.01)
                {
                    Stress_Difference[i][j] = 0;
                }

            }
            fprintf(fp_Stress_Bulk,"\n");
        }

        for(i=0;i<m;i++)
        {
            for(j=0;j<n;++j)
            {
                Curve_Array[i][j]=0;
                if(Stress_Difference[i][j]>=0)
                {
                    if(i-1>=0)
                    {
                        if(Stress_Difference[i-1][j]<0)
                        {
                            Curve_Array[i][j]=1;
                        }
                    }

                    if(i+1<m)
                    {
                        if(Stress_Difference[i+1][j]<0)
                        {
                            Curve_Array[i][j]=1;
                        }
                    }

                    if(j-1>=0)
                    {
                        if(Stress_Difference[i][j-1]<0)
                        {
                            Curve_Array[i][j]=1;
                        }
                    }

                    if(j+1<n)
                    {
                        if(Stress_Difference[i][j+1]<0)
                        {
                            Curve_Array[i][j]=1;
                        }
                    }
                }
                if(j<n-1)
                {
                    if(Curve_Array[i][j]==1)
                    {
                        fprintf(fp_Judge_Phase_Curve,"%d\t",Curve_Array[i][j]);
                    }else
                    {
                        fprintf(fp_Judge_Phase_Curve,"\t");
                    }

                    fprintf(fp_Stress_Difference,"%0.8f\t",Stress_Difference[i][j]);
                }else
                {
                    if(Curve_Array[i][j]==1)
                    {
                        fprintf(fp_Judge_Phase_Curve,"%d",Curve_Array[i][j]);
                    }
                    fprintf(fp_Stress_Difference,"%0.8f",Stress_Difference[i][j]);
                }

                if(Curve_Array[i][j]==1)
                {
                    fprintf(fp_Curve_1D,"%d\t%d\n",j,i);
                }
            }
            fprintf(fp_Judge_Phase_Curve,"\n");
            fprintf(fp_Stress_Difference,"\n");
        }
        fclose(fp_Judge_Phase_Curve);
        fclose(fp_Stress_Difference);
        fclose(fp_Stress_Bulk);
        fclose(fp_Curve_1D);
    }

    free(filename);
    free(filename_Curve_1D);
    free(filename_Stress_Critical);
    free(filename_Difference);
    free(filename_Stress_Bulk);
    free(Stress);
    free(Strain);
    free2DArray(Stress_Difference,m);
    free2DArrayinteger(Curve_Array,m);
    return 1;
}

int Read_NMD_SKN_D(Input_Parameter_NMD *IPN, double **SKN_D)
{
    int i, j, m, n, Version;
    FILE *fp_SKN_D;
    char *filename_SKN_D;
    enum Initial_Mode_NMD{Read_M_3D, Read_M_2D, Initial_M_FE, Initial_M_CC, Read_M_Process};
    Initial_Mode_NMD IMN;

    filename_SKN_D  = Make1DString(500);
    Version = IPN->Version;

    m = 2 * IPN->xn + 3;    n = 2 * IPN->yn + 3;

    IMN = (Initial_Mode_NMD) IPN->InitialMode;
    switch(IMN)
    {
        case Read_M_3D :
        {
            sprintf(filename_SKN_D,"Version %d\\M3D\\z=%d\\SKN_D.dat",Version,0);
            fp_SKN_D = fopen(filename_SKN_D,"r+");

            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    fscanf(fp_SKN_D,"%lf",&SKN_D[i][j]);

                }
            }
            fclose(fp_SKN_D);
            break;
        }

        default :
            break;
    }
    free(filename_SKN_D);

    return 1;
}

int SKN_Boundary_Found(Input_Parameter_NMD *IPN)
{
    int i, j, k, m, n, l, Version;
    double SKN_B, SKN_B_U=0, SKN_B_L=0, SKN_Positif=0, SKN_Area=0, SKN_Target;
    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l  = 2*IPN->zn+3;
    }

    double **SKN_D;
    SKN_D = Make2DArray(m,n);
    Read_NMD_SKN_D(IPN,SKN_D);
    SKN_B = 0;
    //SKN_Target = 0.999;
    SKN_Target = 1;
    printf("Enter the SKN Target\n");
    scanf("%lf",&SKN_Target);
    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            if(SKN_D[i][j]>SKN_B)
            {
                SKN_Positif = SKN_Positif+SKN_D[i][j];
                //printf("%0.8f\t %0.8f\n",SKN_D[i][j],SKN_Positif);
            }
        }
    }
    printf("%0.8f\n",SKN_Positif);

    if(SKN_Positif>SKN_Target)
    {
        //SKN_B = (SKN_Positif-1)/m/n;
        SKN_B = 0.1;
        //printf("%g\n",SKN_B);
        do{
            SKN_Positif = 0;
            //printf("%g\t %g\n",SKN_B_U,SKN_B_L);
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    if(SKN_D[i][j]>SKN_B)
                    {
                        SKN_Positif = SKN_Positif+SKN_D[i][j];
                    }
                }
            }
            //printf("%0.8f\n",SKN_Positif);
            if(SKN_Positif<SKN_Target)
            {
                SKN_B_U = SKN_B;
            }else
            {
                SKN_B_L = SKN_B;
            }
            SKN_B = SKN_B_L*0.618+SKN_B_U*0.382;
        }while(fabs(SKN_Positif-SKN_Target)>0.0001);
    }

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            if(SKN_D[i][j]>SKN_B)
            {
                SKN_Area = SKN_Area+IPN->dxy[0]*IPN->dxy[1];
            }
        }
    }

    printf("%0.8f\t %0.8f\t %0.8f\n",SKN_Positif,SKN_Area,SKN_B);
    return 1;
}

int Stress_Analysis(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, zc, r, Version, D;
    double *Strain, *Stress, *StressAnalysis;
    Strain = Make1DArray(6);
    Stress = Make1DArray(6);
    StressAnalysis = Make1DArray(9);
    char *filename_MaxPS, *filename_MaxPS_x, *filename_MaxPS_y,
         *filename_MinPS, *filename_MinPS_x, *filename_MinPS_y,
         *filename_BulkStress;
    FILE *fp_MaxPS, *fp_MaxPS_x, *fp_MaxPS_y, *fp_MinPS, *fp_MinPS_x, *fp_MinPS_y, *fp_BulkStress;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;
    if (m>=n)
        D=n/20;
    else
        D=m/20;
    //D=min(m,n)/10;

    filename_MaxPS      = Make1DString(500);
    filename_MaxPS_x    = Make1DString(500);
    filename_MaxPS_y    = Make1DString(500);
    filename_MinPS      = Make1DString(500);
    filename_MinPS_x    = Make1DString(500);
    filename_MinPS_y    = Make1DString(500);
    filename_BulkStress = Make1DString(500);

    Version = IPN->Version;
    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }

    for(k=0;k<l;++k)
    {
        sprintf(filename_MaxPS,     "Version %d\\For Printing\\z=%d\\MaxPrincipalStress.dat",  Version,k-zc);
        sprintf(filename_MaxPS_x,   "Version %d\\Elastic Field\\Stress\\z=%d\\MaxPrincipalStress_x.dat",Version,k-zc);
        sprintf(filename_MaxPS_y,   "Version %d\\Elastic Field\\Stress\\z=%d\\MaxPrincipalStress_y.dat",Version,k-zc);
        sprintf(filename_MinPS,     "Version %d\\For Printing\\z=%d\\MinPrincipalStress.dat",  Version,k-zc);
        sprintf(filename_MinPS_x,   "Version %d\\Elastic Field\\Stress\\z=%d\\MinPrincipalStress_x.dat",Version,k-zc);
        sprintf(filename_MinPS_y,   "Version %d\\Elastic Field\\Stress\\z=%d\\MinPrincipalStress_y.dat",Version,k-zc);
        sprintf(filename_BulkStress,"Version %d\\Elastic Field\\Stress\\z=%d\\BulkStress.dat",          Version,k-zc);

        fp_MaxPS      = fopen(filename_MaxPS,"w+");
        fp_MaxPS_x    = fopen(filename_MaxPS_x,"w+");
        fp_MaxPS_y    = fopen(filename_MaxPS_y,"w+");
        fp_MinPS      = fopen(filename_MinPS,"w+");
        fp_MinPS_x    = fopen(filename_MinPS_x,"w+");
        fp_MinPS_y    = fopen(filename_MinPS_y,"w+");
        fp_BulkStress = fopen(filename_BulkStress,"w+");

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                }
                StressAnalysis[0]=(Stress[0]+Stress[1])/2+sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//max principal stress
                StressAnalysis[1]=(Stress[0]+Stress[1])/2-sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//min principal stress
                StressAnalysis[2]=(StressAnalysis[0]+StressAnalysis[1])/3;//bulk stress
                StressAnalysis[3]=atan2(Stress[3],(Stress[1]-StressAnalysis[0]));
                StressAnalysis[4]=atan2(Stress[3],(Stress[1]-StressAnalysis[1]));

                StressAnalysis[5]=StressAnalysis[0]*cos(StressAnalysis[3]);//max principal stress along x
                StressAnalysis[6]=StressAnalysis[0]*sin(StressAnalysis[3]);//max principal stress along y
                StressAnalysis[7]=StressAnalysis[1]*cos(StressAnalysis[4]);//min principal stress along x
                StressAnalysis[8]=StressAnalysis[1]*sin(StressAnalysis[4]);//min principal stress along y

                if(BS->Boundary[i][j][k]!=0)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_MaxPS,     "%0.8f\t",StressAnalysis[0]);
                        fprintf(fp_MaxPS_x,   "%0.8f\t",StressAnalysis[5]);
                        fprintf(fp_MaxPS_y,   "%0.8f\t",StressAnalysis[6]);
                        fprintf(fp_MinPS,     "%0.8f\t",StressAnalysis[1]);
                        fprintf(fp_MinPS_x,   "%0.8f\t",StressAnalysis[7]);
                        fprintf(fp_MinPS_y,   "%0.8f\t",StressAnalysis[8]);
                        fprintf(fp_BulkStress,"%0.8f\t",StressAnalysis[2]);
                    }else
                    {
                        fprintf(fp_MaxPS,     "%0.8f",StressAnalysis[0]);
                        fprintf(fp_MaxPS_x,   "%0.8f",StressAnalysis[5]);
                        fprintf(fp_MaxPS_y,   "%0.8f",StressAnalysis[6]);
                        fprintf(fp_MinPS,     "%0.8f",StressAnalysis[1]);
                        fprintf(fp_MinPS_x,   "%0.8f",StressAnalysis[7]);
                        fprintf(fp_MinPS_y,   "%0.8f",StressAnalysis[8]);
                        fprintf(fp_BulkStress,"%0.8f",StressAnalysis[2]);
                    }
                }else
                {
                    if(j<n-1)
                    {
                        fprintf(fp_MaxPS,     "\t");
                        fprintf(fp_MaxPS_x,   "\t");
                        fprintf(fp_MaxPS_y,   "\t");
                        fprintf(fp_MinPS,     "\t");
                        fprintf(fp_MinPS_x,   "\t");
                        fprintf(fp_MinPS_y,   "\t");
                        fprintf(fp_BulkStress,"\t");
                    }
                }
/*
                if(i%D==0 && j%D==0)
                {
                    if(i-D<0||i+D>m||j-D<0||j+D>m)
                    {
                        fprintf(fp_MaxPS,"%d\t%d\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",i,j,0.,0.,0.,0.);
                        fprintf(fp_MinPS,"%d\t%d\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",i,j,0.,0.,0.,0.);
                    }
                    else
                    {
                        fprintf(fp_MaxPS,"%d\t%d\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",i,j,StressAnalysis[3],StressAnalysis[0],StressAnalysis[3]+PI,StressAnalysis[0]);
                        fprintf(fp_MinPS,"%d\t%d\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",i,j,StressAnalysis[4],StressAnalysis[1],StressAnalysis[4]+PI,StressAnalysis[1]);
                    }
                }*/
            }
            fprintf(fp_MaxPS,     "\n");
            fprintf(fp_MaxPS_x,   "\n");
            fprintf(fp_MaxPS_y,   "\n");
            fprintf(fp_MinPS,     "\n");
            fprintf(fp_MinPS_x,   "\n");
            fprintf(fp_MinPS_y,   "\n");
            fprintf(fp_BulkStress,"\n");
        }
        fclose(fp_MaxPS);
        fclose(fp_MaxPS_x);
        fclose(fp_MaxPS_y);
        fclose(fp_MinPS);
        fclose(fp_MinPS_x);
        fclose(fp_MinPS_y);
        fclose(fp_BulkStress);
    }
    free(filename_MaxPS);
    free(filename_MaxPS_x);
    free(filename_MaxPS_y);
    free(filename_MinPS);
    free(filename_MinPS_x);
    free(filename_MinPS_y);
    free(filename_BulkStress);

    free(Strain);
    free(Stress);
    free(StressAnalysis);

    return 1;
}

int Write_Principal_Stress(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, zc, r, Version, ibegin, jbegin;
    double *Strain, *Stress, *StressAnalysis, MaxPositive, MaxNegative, MinPositive,
           MinNegative, ScaleMaxPositive, ScaleMaxNegative, ScaleMinPositive, ScaleMinNegative, MaxEndX, MaxEndY, MinEndX, MinEndY;
    MaxPositive=0;  MaxNegative=0; MinPositive=0; MinNegative=0;
    Strain = Make1DArray(6);
    Stress = Make1DArray(6);
    StressAnalysis = Make1DArray(9);

    double S_Appropriate, S_Min, S_Module, S_i, S_j, Ratio;
    int Inplane_Density, D;
    char *filename_MaxPS, *filename_MinPS, *filename_BulkStress;
    FILE *fp_MaxPS, *fp_MinPS, *fp_BulkStress;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    Inplane_Density = 12;
    if (m>=n)
    {
        D=n/Inplane_Density;
    }else
    {
        D=m/Inplane_Density;
    }

    ibegin=((m-1)/2)%D;
    jbegin=((n-1)/2)%D;

    S_Appropriate = D*0.4;
    S_Min         = D*0.02;

    filename_MaxPS          = Make1DString(500);
    filename_MinPS          = Make1DString(500);
    filename_BulkStress     = Make1DString(500);

    Version = IPN->Version;
    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }
    for(k=0;k<l;++k)
    {
        sprintf(filename_MaxPS, "Version %d\\For Printing\\z=%d\\MaxPrincipalStress.dat",  Version,k-zc);
        sprintf(filename_MinPS,  "Version %d\\For Printing\\z=%d\\MinPrincipalStress.dat",  Version,k-zc);
        sprintf(filename_BulkStress,"Version %d\\Elastic Field\\Stress\\z=%d\\BulkStress.dat",  Version,k-zc);

        fp_MaxPS         = fopen(filename_MaxPS,"w+");
        fp_MinPS         = fopen(filename_MinPS,"w+");
        fp_BulkStress = fopen(filename_BulkStress,"w+");
        for(i=ibegin;i<m;i+=D)
        {
            for(j=jbegin;j<n;j+=D)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                }
                StressAnalysis[0]=(Stress[0]+Stress[1])/2+sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//max principal stress
                if (StressAnalysis[0]>=0)
                {
                    if (StressAnalysis[0]>=MaxPositive)
                    {
                        MaxPositive=fabs(StressAnalysis[0]);
                    }
                }
                else
                {
                    if (fabs(StressAnalysis[0])>=MaxNegative)
                    {
                        MaxNegative=fabs(StressAnalysis[0]);
                    }
                }
                StressAnalysis[1]=(Stress[0]+Stress[1])/2-sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//min principal stress
                if (StressAnalysis[1]>=0)
                {
                    if (StressAnalysis[1]>=MinPositive)
                    {
                        MinPositive=fabs(StressAnalysis[1]);
                    }
                }
                else
                {
                    if (fabs(StressAnalysis[1])>=MinNegative)
                    {
                        MinNegative=fabs(StressAnalysis[1]);
                    }
                }
            }
        }

        ScaleMaxPositive=D/1.5/MaxPositive;
        ScaleMaxNegative=D/1.5/MaxNegative;
        ScaleMinPositive=D/1.5/MinPositive;
        ScaleMinNegative=D/1.5/MinNegative;

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                    StrainToStress(Stress,Strain);
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                }
                StressAnalysis[0]=(Stress[0]+Stress[1])/2+sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//max principal stress
                StressAnalysis[1]=(Stress[0]+Stress[1])/2-sqrt(((Stress[0]-Stress[1])/2)*((Stress[0]-Stress[1])/2)+Stress[3]*Stress[3]);//min principal stress
                StressAnalysis[2]=(StressAnalysis[0]+StressAnalysis[1])/3;//bulk stress
                StressAnalysis[3]=atan2(Stress[3],(Stress[1]-StressAnalysis[0]));
                StressAnalysis[4]=atan2(Stress[3],(Stress[1]-StressAnalysis[1]));
                if (StressAnalysis[0]>=0)
                {
                    StressAnalysis[5]=ScaleMaxPositive*(StressAnalysis[0]*cos(StressAnalysis[3]));//max principal stress along x
                    StressAnalysis[6]=ScaleMaxPositive*(StressAnalysis[0]*sin(StressAnalysis[3]));//max principal stress along y
                }
                else
                {
                    StressAnalysis[5]=ScaleMaxNegative*(StressAnalysis[0]*cos(StressAnalysis[3]));//max principal stress along x
                    StressAnalysis[6]=ScaleMaxNegative*(StressAnalysis[0]*sin(StressAnalysis[3]));//max principal stress along y
                }

                if (StressAnalysis[1]>=0)
                {
                    StressAnalysis[7]=ScaleMinPositive*(StressAnalysis[1]*cos(StressAnalysis[4]));//min principal stress along x
                    StressAnalysis[8]=ScaleMinPositive*(StressAnalysis[1]*sin(StressAnalysis[4]));//min principal stress along y
                }
                else
                {
                    StressAnalysis[7]=ScaleMinNegative*(StressAnalysis[1]*cos(StressAnalysis[4]));//min principal stress along x
                    StressAnalysis[8]=ScaleMinNegative*(StressAnalysis[1]*sin(StressAnalysis[4]));//min principal stress along y
                }

                if(BS->Boundary[i][j][k]!=0)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_BulkStress,"%0.8f\t",StressAnalysis[2]);
                    }
                    else
                    {
                        fprintf(fp_BulkStress,"%0.8f",StressAnalysis[2]);
                    }
                }else
                {
                    if(j<n-1)
                    {
                        fprintf(fp_BulkStress,"\t");
                    }
                }


                if(i%D==ibegin && j%D==jbegin)
                {
                    S_Module = sqrt(StressAnalysis[5]*StressAnalysis[5]+StressAnalysis[6]*StressAnalysis[6]);
                    S_i   = StressAnalysis[5];
                    S_j   = StressAnalysis[6];
                    if(S_Module>S_Appropriate)
                    {
                        Ratio = S_Appropriate/S_Module;
                        S_i   = StressAnalysis[5]*Ratio;
                        S_j   = StressAnalysis[6]*Ratio;
                    }else if(S_Module<S_Min)
                    {
                        S_i = 0;
                        S_j = 0;
                    }

                    MaxEndX = i+S_i;  MaxEndY = j+S_j;
                    MinEndX = i-S_i;  MinEndY = j-S_j;
                    if(MaxEndX<0||MaxEndX>m-1||MinEndX<0||MinEndX>m-1||MaxEndY<0||MaxEndY>n-1||MinEndY<0||MinEndY>n-1)
                    {
                        MaxEndX=i;
                        MinEndX=i;
                        MaxEndY=j;
                        MinEndY=j;
                    }

                    if(BS->Boundary[i][j][k]!=0)
                    {
                        if (StressAnalysis[0]>=0)
                        {
                            fprintf(fp_MaxPS,"%d\t%d\t%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\n",
                                    j, i, MaxEndY, MaxEndX, j, i, MinEndY, MinEndX);
                        }
                        else
                        {
                            fprintf(fp_MaxPS,"%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\t%d\t%d\n",
                                    MaxEndY, MaxEndX, j, i, MinEndY, MinEndX, j, i);
                        }
                    }

                    S_Module = sqrt(StressAnalysis[7]*StressAnalysis[7]+StressAnalysis[8]*StressAnalysis[8]);
                    S_i   = StressAnalysis[7];
                    S_j   = StressAnalysis[8];
                    if(S_Module>S_Appropriate)
                    {
                        Ratio = S_Appropriate/S_Module;
                        S_i   = StressAnalysis[7]*Ratio;
                        S_j   = StressAnalysis[8]*Ratio;
                    }else if(S_Module<S_Min)
                    {
                        S_i = 0;
                        S_j = 0;
                    }

                    MaxEndX = i+S_i;   MaxEndY = j+S_j;
                    MinEndX = i-S_i;   MinEndY = j-S_j;
                    if(MaxEndX<0||MaxEndX>m-1||MinEndX<0||MinEndX>m-1||MaxEndY<0||MaxEndY>n-1||MinEndY<0||MinEndY>n-1)
                    {
                        MaxEndX=i;
                        MinEndX=i;
                        MaxEndY=j;
                        MinEndY=j;
                    }

                    if(BS->Boundary[i][j][k]!=0)
                    {
                        if (StressAnalysis[1]>=0)
                        {
                            fprintf(fp_MinPS,"%d\t%d\t%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\n",
                                    j, i, MaxEndY, MaxEndX, j, i, MinEndY, MinEndX);
                        }
                        else
                        {
                            fprintf(fp_MinPS,"%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\t%d\t%d\n",
                                    MaxEndY, MaxEndX, j, i, MinEndY, MinEndX, j, i);
                        }
                    }
                }
            }
            fprintf(fp_BulkStress,"\n");
        }
        fclose(fp_MaxPS);
        fclose(fp_MinPS);
        fclose(fp_BulkStress);
    }
    free(filename_MaxPS);
    free(filename_MinPS);
    free(filename_BulkStress);

    free(Strain);
    free(Stress);
    free(StressAnalysis);

    return 1;
}

int Write_Principal_Strain(Input_Parameter_NMD *IPN, double ****CST, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, zc, r, Version, ibegin, jbegin;
    double *Strain, *Stress, *StressAnalysis, MaxPositive, MaxNegative, MinPositive,
           MinNegative, ScaleMaxPositive, ScaleMaxNegative, ScaleMinPositive, ScaleMinNegative, MaxEndX, MaxEndY, MinEndX, MinEndY;
    MaxPositive=0;  MaxNegative=0; MinPositive=0; MinNegative=0;
    Strain = Make1DArray(6);
    Stress = Make1DArray(6);
    StressAnalysis = Make1DArray(9);

    double S_Appropriate, S_Min, S_Module, S_i, S_j, Ratio;
    int Inplane_Density, D;
    char *filename_MaxPS, *filename_MinPS, *filename_BulkStress;
    FILE *fp_MaxPS, *fp_MinPS, *fp_BulkStress;

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    Inplane_Density = 12;
    if (m>=n)
    {
        D=n/Inplane_Density;
    }else
    {
        D=m/Inplane_Density;
    }

    ibegin=((m-1)/2)%D;
    jbegin=((n-1)/2)%D;

    S_Appropriate = D*0.4;
    S_Min         = D*0.02;

    filename_MaxPS          = Make1DString(500);
    filename_MinPS          = Make1DString(500);
    filename_BulkStress     = Make1DString(500);

    Version = IPN->Version;
    if(IPN->zn==0)
    {
        l=1;
        zc=0;
    }else
    {
        l  = 2*IPN->zn+3;
        zc = IPN->zn-1;
    }
    for(k=0;k<l;++k)
    {
        sprintf(filename_MaxPS, "Version %d\\For Printing\\z=%d\\MaxPrincipalStrain.dat",  Version,k-zc);
        sprintf(filename_MinPS,  "Version %d\\For Printing\\z=%d\\MinPrincipalStrain.dat",  Version,k-zc);
        sprintf(filename_BulkStress,"Version %d\\Elastic Field\\Stress\\z=%d\\BulkStrain.dat",  Version,k-zc);

        fp_MaxPS         = fopen(filename_MaxPS,"w+");
        fp_MinPS         = fopen(filename_MinPS,"w+");
        fp_BulkStress = fopen(filename_BulkStress,"w+");
        for(i=ibegin;i<m;i+=D)
        {
            for(j=jbegin;j<n;j+=D)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                    StressToStrain(Stress,Strain);
                }
                StressAnalysis[0]=(Strain[0]+Strain[1])/2+sqrt(((Strain[0]-Strain[1])/2)*((Strain[0]-Strain[1])/2)+Strain[3]*Strain[3]);//max principal strain
                if (StressAnalysis[0]>=0)
                {
                    if (StressAnalysis[0]>=MaxPositive && BS->Boundary[i][j][k]!=0)
                    {
                        MaxPositive=fabs(StressAnalysis[0]);
                    }
                }
                else
                {
                    if (fabs(StressAnalysis[0])>=MaxNegative && BS->Boundary[i][j][k]!=0)
                    {
                        MaxNegative=fabs(StressAnalysis[0]);
                    }
                }
                StressAnalysis[1]=(Strain[0]+Strain[1])/2-sqrt(((Strain[0]-Strain[1])/2)*((Strain[0]-Strain[1])/2)+Strain[3]*Strain[3]);//min principal strain
                if (StressAnalysis[1]>=0)
                {
                    if (StressAnalysis[1]>=MinPositive && BS->Boundary[i][j][k]!=0)
                    {
                        MinPositive=fabs(StressAnalysis[1]);
                    }
                }
                else
                {
                    if (fabs(StressAnalysis[1])>=MinNegative && BS->Boundary[i][j][k]!=0)
                    {
                        MinNegative=fabs(StressAnalysis[1]);
                    }
                }
            }
        }

        ScaleMaxPositive=D/1.5/MaxPositive;
        ScaleMaxNegative=D/1.5/MaxNegative;
        ScaleMinPositive=D/1.5/MinPositive;
        ScaleMinNegative=D/1.5/MinNegative;

        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(CST[6][i][j][k]==0)//Stain Type
                {
                    for(r=0;r<6;++r)
                    {
                        Strain[r]=CST[r][i][j][k];
                    }
                }else//Stress Type
                {
                    for(r=0;r<6;++r)
                    {
                        Stress[r]=CST[r][i][j][k];
                    }
                    StressToStrain(Stress,Strain);
                }
                StressAnalysis[0]=(Strain[0]+Strain[1])/2+sqrt(((Strain[0]-Strain[1])/2)*((Strain[0]-Strain[1])/2)+Strain[3]*Strain[3]);//max principal stress
                StressAnalysis[1]=(Strain[0]+Strain[1])/2-sqrt(((Strain[0]-Strain[1])/2)*((Strain[0]-Strain[1])/2)+Strain[3]*Strain[3]);//min principal stress
                StressAnalysis[2]=(StressAnalysis[0]+StressAnalysis[1])/3;//bulk stress
                StressAnalysis[3]=atan2(Strain[3],(Strain[1]-StressAnalysis[0]));
                StressAnalysis[4]=atan2(Strain[3],(Strain[1]-StressAnalysis[1]));
                if (StressAnalysis[0]>=0)
                {
                    StressAnalysis[5]=ScaleMaxPositive*(StressAnalysis[0]*cos(StressAnalysis[3]));//max principal stress along x
                    StressAnalysis[6]=ScaleMaxPositive*(StressAnalysis[0]*sin(StressAnalysis[3]));//max principal stress along y
                }
                else
                {
                    StressAnalysis[5]=ScaleMaxNegative*(StressAnalysis[0]*cos(StressAnalysis[3]));//max principal stress along x
                    StressAnalysis[6]=ScaleMaxNegative*(StressAnalysis[0]*sin(StressAnalysis[3]));//max principal stress along y
                }

                if (StressAnalysis[1]>=0)
                {
                    StressAnalysis[7]=ScaleMinPositive*(StressAnalysis[1]*cos(StressAnalysis[4]));//min principal stress along x
                    StressAnalysis[8]=ScaleMinPositive*(StressAnalysis[1]*sin(StressAnalysis[4]));//min principal stress along y
                }
                else
                {
                    StressAnalysis[7]=ScaleMinNegative*(StressAnalysis[1]*cos(StressAnalysis[4]));//min principal stress along x
                    StressAnalysis[8]=ScaleMinNegative*(StressAnalysis[1]*sin(StressAnalysis[4]));//min principal stress along y
                }

                if(BS->Boundary[i][j][k]!=0)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_BulkStress,"%0.8f\t",StressAnalysis[2]);
                    }
                    else
                    {
                        fprintf(fp_BulkStress,"%0.8f",StressAnalysis[2]);
                    }
                }else
                {
                    if(j<n-1)
                    {
                        fprintf(fp_BulkStress,"\t");
                    }
                }


                if(i%D==ibegin && j%D==jbegin)
                {
                    S_Module = sqrt(StressAnalysis[5]*StressAnalysis[5]+StressAnalysis[6]*StressAnalysis[6]);
                    S_i   = StressAnalysis[5];
                    S_j   = StressAnalysis[6];
                    if(S_Module>S_Appropriate)
                    {
                        Ratio = S_Appropriate/S_Module;
                        S_i   = StressAnalysis[5]*Ratio;
                        S_j   = StressAnalysis[6]*Ratio;
                    }else if(S_Module<S_Min)
                    {
                        S_i = 0;
                        S_j = 0;
                    }

                    MaxEndX = i+S_i;  MaxEndY = j+S_j;
                    MinEndX = i-S_i;  MinEndY = j-S_j;
                    if(MaxEndX<0||MaxEndX>m-1||MinEndX<0||MinEndX>m-1||MaxEndY<0||MaxEndY>n-1||MinEndY<0||MinEndY>n-1)
                    {
                        MaxEndX=i;
                        MinEndX=i;
                        MaxEndY=j;
                        MinEndY=j;
                    }

                    if(BS->Boundary[i][j][k]!=0)
                    {
                        if (StressAnalysis[0]>=0)
                        {
                            fprintf(fp_MaxPS,"%d\t%d\t%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\n",
                                    j, i, MaxEndY, MaxEndX, j, i, MinEndY, MinEndX);
                        }
                        else
                        {
                            fprintf(fp_MaxPS,"%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\t%d\t%d\n",
                                    MaxEndY, MaxEndX, j, i, MinEndY, MinEndX, j, i);
                        }
                    }

                    S_Module = sqrt(StressAnalysis[7]*StressAnalysis[7]+StressAnalysis[8]*StressAnalysis[8]);
                    S_i   = StressAnalysis[7];
                    S_j   = StressAnalysis[8];
                    if(S_Module>S_Appropriate)
                    {
                        Ratio = S_Appropriate/S_Module;
                        S_i   = StressAnalysis[7]*Ratio;
                        S_j   = StressAnalysis[8]*Ratio;
                    }else if(S_Module<S_Min)
                    {
                        S_i = 0;
                        S_j = 0;
                    }

                    MaxEndX = i+S_i;   MaxEndY = j+S_j;
                    MinEndX = i-S_i;   MinEndY = j-S_j;
                    if(MaxEndX<0||MaxEndX>m-1||MinEndX<0||MinEndX>m-1||MaxEndY<0||MaxEndY>n-1||MinEndY<0||MinEndY>n-1)
                    {
                        MaxEndX=i;
                        MinEndX=i;
                        MaxEndY=j;
                        MinEndY=j;
                    }

                    if(BS->Boundary[i][j][k]!=0)
                    {
                        if (StressAnalysis[1]>=0)
                        {
                            fprintf(fp_MinPS,"%d\t%d\t%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\n",
                                    j, i, MaxEndY, MaxEndX, j, i, MinEndY, MinEndX);
                        }
                        else
                        {
                            fprintf(fp_MinPS,"%0.8f\t%0.8f\t%d\t%d\t%0.8f\t%0.8f\t%d\t%d\n",
                                    MaxEndY, MaxEndX, j, i, MinEndY, MinEndX, j, i);
                        }
                    }
                }
            }
            fprintf(fp_BulkStress,"\n");
        }
        fclose(fp_MaxPS);
        fclose(fp_MinPS);
        fclose(fp_BulkStress);
    }
    free(filename_MaxPS);
    free(filename_MinPS);
    free(filename_BulkStress);

    free(Strain);
    free(Stress);
    free(StressAnalysis);

    return 1;
}

double Determinant_4D(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4)
{
    double Det;

    Det=(x1*(x4*(y2-y3)*(y1-y4)-x3*(y1-y3)*(y2-y4)+x2*(y1-y2)*(y3-y4))+
         x2*x3*(y2-y3)*(y1-y4)-x2*x4*(y1-y3)*(y2-y4)+x3*x4*(y1-y2)*(y3-y4));

    return Det;
}

int Calculate_Interpolation(double **Interpolate_Coordinate, double *Interpolate_Values)
{
    double x1, y1, x2, y2, x3, y3, x4, y4, Det;
    //double x0, y0;
    //x0 = Interpolate_Coordinate[0][0];  y0 = Interpolate_Coordinate[0][1];
    x1 = Interpolate_Coordinate[1][0];  y1 = Interpolate_Coordinate[1][1];
    x2 = Interpolate_Coordinate[2][0];  y2 = Interpolate_Coordinate[2][1];
    x3 = Interpolate_Coordinate[3][0];  y3 = Interpolate_Coordinate[3][1];
    x4 = Interpolate_Coordinate[4][0];  y4 = Interpolate_Coordinate[4][1];

    Det = Determinant_4D(x1,y1,x2,y2,x3,y3,x4,y4);



    if(Det==0)
    {
        Interpolate_Values[0] = -1000;
        Interpolate_Values[1] = -1000;
        Interpolate_Values[2] = -1000;
        Interpolate_Values[3] = -1000;
        printf("DET=0\n");
    }else
    {
        Interpolate_Values[0] = (x2*x4*y3*(y2-y4)+x2*x3*y4*(y3-y2)+x3*x4*y2*(y4-y3))/
                                (x1*(x4*(y2-y3)*(y1-y4)-x3*(y1-y3)*(y2-y4)+x2*(y1-y2)*(y3-y4))+
                                 x2*x3*(y2-y3)*(y1-y4)-x2*x4*(y1-y3)*(y2-y4)+x3*x4*(y1-y2)*(y3-y4));

        Interpolate_Values[1] = (x1*x4*y3*(y1-y4)+x1*x3*y4*(y3-y1)+x3*x4*y1*(y4-y3))/
                                (x2*(x4*(y1-y3)*(y2-y4)-x3*(y2-y3)*(y1-y4)+x1*(y2-y1)*(y3-y4))+
                                 x1*x3*(y1-y3)*(y2-y4)-x1*x4*(y2-y3)*(y1-y4)+x3*x4*(y2-y1)*(y3-y4));

        Interpolate_Values[2] = (x2*x4*y1*(y2-y4)+x2*x1*y4*(y1-y2)+x1*x4*y2*(y4-y1))/
                                (x3*(x4*(y2-y1)*(y3-y4)-x1*(y3-y1)*(y2-y4)+x2*(y3-y2)*(y1-y4))+
                                 x2*x1*(y2-y1)*(y3-y4)-x2*x4*(y3-y1)*(y2-y4)+x1*x4*(y3-y2)*(y1-y4));

        Interpolate_Values[3] = (x2*x1*y3*(y2-y1)+x2*x3*y1*(y3-y2)+x3*x1*y2*(y1-y3))/
                                (x4*(x1*(y2-y3)*(y4-y1)-x3*(y4-y3)*(y2-y1)+x2*(y4-y2)*(y3-y1))+
                                 x2*x3*(y2-y3)*(y4-y1)-x2*x1*(y4-y3)*(y2-y1)+x3*x1*(y4-y2)*(y3-y1));
    }

    return 1;
}

int Interpolation_NMD_2D(int N_Interpolate_Points, double **Interpolate_Points_Parameter, double *Interpolate_Coefficient)
{
    enum Interpolate_Point_Case{IP_1=1, IP_2, IP_3, IP_4};
    Interpolate_Point_Case IPC;
    double x1, y1, x2, y2, x3, y3, x4, y4, f4, xt, yt, c1, c2, c3, c4;

    IPC = (Interpolate_Point_Case)N_Interpolate_Points;
    switch (IPC)
    {
        /////////////////////////// 1 Point Case ///////////////////////////////////////////
        case IP_1:
        {
            c1 = 1.;    c2 = 0.;    c3 = 0.;    c4 = 0.;
            break;
        }
        /////////////////////////// 2 Points Case ///////////////////////////////////////////
        case IP_2:
        {
            xt = Interpolate_Points_Parameter[4][0];         yt = Interpolate_Points_Parameter[4][1];
            x1 = Interpolate_Points_Parameter[0][0] - xt;    y1 = Interpolate_Points_Parameter[0][1] - yt;
            x2 = Interpolate_Points_Parameter[1][0] - xt;    y2 = Interpolate_Points_Parameter[1][1] - yt;
            ////////////// Line Parameter: Ax+By+C=0 ////////////////
            double Line_A, Line_B, Line_C, Distance_1, Distance_2, Distance_TargetToLine, Distance_1_1D, Distance_2_1D;
            Line_A = y1 - y2;
            Line_B = x2 - x1;
            Line_C = x1 * y2 - x2 * y1;
            //Distance_TargetToLine = fabs(Line_A * xt + Line_B * yt + Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
            //Distance_1 = sqrt(pow((x1 - xt), 2) + pow((y1 - yt), 2));
            //Distance_2 = sqrt(pow((x2 - xt), 2) + pow((y2 - yt), 2));
            Distance_TargetToLine = fabs(Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
            Distance_1 = sqrt(pow((x1), 2) + pow((y1), 2));
            Distance_2 = sqrt(pow((x2), 2) + pow((y2), 2));
            Distance_1_1D = sqrt(pow(Distance_1, 2) - pow(Distance_TargetToLine, 2));
            Distance_2_1D = sqrt(pow(Distance_2, 2) - pow(Distance_TargetToLine, 2));

            c1 = Distance_2_1D / (Distance_1_1D + Distance_2_1D);
            c2 = Distance_1_1D / (Distance_1_1D + Distance_2_1D);
            c3 = 0.;
            c4 = 0.;
            break;
        }
        /////////////////////////// 3 Points Case ///////////////////////////////////////////
        case IP_3:
        {
            xt = Interpolate_Points_Parameter[4][0];         yt = Interpolate_Points_Parameter[4][1];
            x1 = Interpolate_Points_Parameter[0][0] - xt;    y1 = Interpolate_Points_Parameter[0][1] - yt;
            x2 = Interpolate_Points_Parameter[1][0] - xt;    y2 = Interpolate_Points_Parameter[1][1] - yt;
            x3 = Interpolate_Points_Parameter[2][0] - xt;    y3 = Interpolate_Points_Parameter[2][1] - yt;
            ////////////// Line Parameter: Ax+By+C=0 ////////////////
            ////////////// Line 12: Cross 1 and 2 Point /////////////
            double Line_A_12, Line_B_12, Line_C_12;
            Line_A_12 = (y1 - y2) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            Line_B_12 = (x2 - x1) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            Line_C_12 = (x1 * y2 - x2 * y1) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            ////////////// Line 13: Cross 1 and 3 Point /////////////
            double Line_A_13, Line_B_13, Line_C_13;
            Line_A_13 = (y1 - y3) / (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1));
            Line_B_13 = (x3 - x1) / (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1));
            Line_C_13 = (x1 * y3 - x3 * y1) / (x1 * (y3 - y2) + x2 * (y1 - y3) + x3 * (y2 - y1));
            ////////////// Line 23: Cross 2 and 3 Point /////////////
            double Line_A_23, Line_B_23, Line_C_23;
            Line_A_23 = (y2 - y3) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            Line_B_23 = (x3 - x2) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            Line_C_23 = (x2 * y3 - x3 * y2) / (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
            /////////////////////////////////// Justification for Triangle Inside //////////////////////////////////////
            //c1 = Line_A_23 * xt + Line_B_23 * yt + Line_C_23;
            //c2 = Line_A_13 * xt + Line_B_13 * yt + Line_C_13;
            //c3 = Line_A_12 * xt + Line_B_12 * yt + Line_C_12;
            //c1 = Line_A_23 * xt + Line_B_23 * yt + Line_C_23;
            c1 = Line_C_23;
            c2 = Line_C_13;
            c3 = Line_C_12;
            c4 = 0.;

            if (c1 < 0 || c2 < 0 || c3 < 0)
            {
                ////////////// Line Parameter: Ax+By+C=0 ////////////////
                double Line_A, Line_B, Line_C, Distance_1, Distance_2, Distance_TargetToLine, Distance_1_1D, Distance_2_1D;
                if (fabs(c1) <= fabs(c2) && fabs(c1) <= fabs(c3))
                {
                    ///////////////////// First Point=Point 2, Second Point=Point 3 /////////////////////////
                    Line_A = Line_A_23;
                    Line_B = Line_B_23;
                    Line_C = Line_C_23;
                    //Distance_TargetToLine = fabs(Line_A * xt + Line_B * yt + Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    //Distance_1 = sqrt(pow((x2 - xt), 2) + pow((y2 - yt), 2));
                    //Distance_2 = sqrt(pow((x3 - xt), 2) + pow((y3 - yt), 2));
                    Distance_TargetToLine = fabs(Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    Distance_1 = sqrt(pow((x2), 2) + pow((y2), 2));
                    Distance_2 = sqrt(pow((x3), 2) + pow((y3), 2));
                    Distance_1_1D = sqrt(pow(Distance_1, 2) - pow(Distance_TargetToLine, 2));
                    Distance_2_1D = sqrt(pow(Distance_2, 2) - pow(Distance_TargetToLine, 2));

                    c1 = 0.;
                    c2 = Distance_2_1D / (Distance_1_1D + Distance_2_1D);
                    c3 = Distance_1_1D / (Distance_1_1D + Distance_2_1D);
                    c4 = 0.;
                }else if(fabs(c2) < fabs(c1) && fabs(c2) <= fabs(c3))
                {
                    ///////////////////// First Point=Point 1, Second Point=Point 3 /////////////////////////
                    Line_A = Line_A_13;
                    Line_B = Line_B_13;
                    Line_C = Line_C_13;
                    //Distance_TargetToLine = fabs(Line_A * xt + Line_B * yt + Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    //Distance_1 = sqrt(pow((x1 - xt), 2) + pow((y1 - yt), 2));
                    //Distance_2 = sqrt(pow((x3 - xt), 2) + pow((y3 - yt), 2));
                    Distance_TargetToLine = fabs(Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    Distance_1 = sqrt(pow((x1), 2) + pow((y1), 2));
                    Distance_2 = sqrt(pow((x3), 2) + pow((y3), 2));
                    Distance_1_1D = sqrt(pow(Distance_1, 2) - pow(Distance_TargetToLine, 2));
                    Distance_2_1D = sqrt(pow(Distance_2, 2) - pow(Distance_TargetToLine, 2));

                    c1 = Distance_2_1D / (Distance_1_1D + Distance_2_1D);
                    c2 = 0.;
                    c3 = Distance_1_1D / (Distance_1_1D + Distance_2_1D);
                    c4 = 0.;
                }else if(fabs(c3) < fabs(c1) && fabs(c3) < fabs(c2))
                {
                    ///////////////////// First Point=Point 1, Second Point=Point 2 /////////////////////////
                    Line_A = Line_A_12;
                    Line_B = Line_B_12;
                    Line_C = Line_C_12;
                    //Distance_TargetToLine = fabs(Line_A * xt + Line_B * yt + Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    //Distance_1 = sqrt(pow((x1 - xt), 2) + pow((y1 - yt), 2));
                    //Distance_2 = sqrt(pow((x2 - xt), 2) + pow((y2 - yt), 2));
                    Distance_TargetToLine = fabs(Line_C) / sqrt(Line_A * Line_A + Line_B * Line_B);
                    Distance_1 = sqrt(pow((x1), 2) + pow((y1), 2));
                    Distance_2 = sqrt(pow((x2), 2) + pow((y2), 2));
                    Distance_1_1D = sqrt(pow(Distance_1, 2) - pow(Distance_TargetToLine, 2));
                    Distance_2_1D = sqrt(pow(Distance_2, 2) - pow(Distance_TargetToLine, 2));

                    c1 = Distance_2_1D / (Distance_1_1D + Distance_2_1D);
                    c2 = Distance_1_1D / (Distance_1_1D + Distance_2_1D);
                    c3 = 0.;
                    c4 = 0.;
                }
            }
            break;
        }
        /////////////////////////// 4 Points Case ///////////////////////////////////////////
        case IP_4:
        {
            xt = Interpolate_Points_Parameter[4][0];         yt = Interpolate_Points_Parameter[4][1];
            x1 = Interpolate_Points_Parameter[0][0] - xt;    y1 = Interpolate_Points_Parameter[0][1] - yt;
            x2 = Interpolate_Points_Parameter[1][0] - xt;    y2 = Interpolate_Points_Parameter[1][1] - yt;
            x3 = Interpolate_Points_Parameter[2][0] - xt;    y3 = Interpolate_Points_Parameter[2][1] - yt;
            x4 = Interpolate_Points_Parameter[3][0] - xt;    y4 = Interpolate_Points_Parameter[3][1] - yt;
            ////////////// Curve Parameter: Ax+By+Cxy+D=0 ////////////////
            c1 = (x2 * x4 * y3 * (y2 - y4) + x2 * x3 * y4 * (y3 - y2) + x3 * x4 * y2 * (y4 - y3)) /
                 (x1 * (x4 * (y2 - y3) * (y1 - y4) - x3 * (y1 - y3) * (y2 - y4) + x2 * (y1 - y2) * (y3 - y4)) +
                  x2 * x3 * (y2 - y3) * (y1 - y4) - x2 * x4 * (y1 - y3) * (y2 - y4) + x3 * x4 * (y1 - y2) * (y3 - y4));

            c2 = (x1 * x4 * y3 * (y1 - y4) + x1 * x3 * y4 * (y3 - y1) + x3 * x4 * y1 * (y4 - y3)) /
                 (x2 * (x4 * (y1 - y3) * (y2 - y4) - x3 * (y2 - y3) * (y1 - y4) + x1 * (y2 - y1) * (y3 - y4)) +
                  x1 * x3 * (y1 - y3) * (y2 - y4) - x1 * x4 * (y2 - y3) * (y1 - y4) + x3 * x4 * (y2 - y1) * (y3 - y4));

            c3 = (x2 * x4 * y1 * (y2 - y4) + x2 * x1 * y4 * (y1 - y2) + x1 * x4 * y2 * (y4 - y1)) /
                 (x3 * (x4 * (y2 - y1) * (y3 - y4) - x1 * (y3 - y1) * (y2 - y4) + x2 * (y3 - y2) * (y1 - y4)) +
                  x2 * x1 * (y2 - y1) * (y3 - y4) - x2 * x4 * (y3 - y1) * (y2 - y4) + x1 * x4 * (y3 - y2) * (y1 - y4));

            c4 = (x2 * x1 * y3 * (y2 - y1) + x2 * x3 * y1 * (y3 - y2) + x3 * x1 * y2 * (y1 - y3)) /
                 (x4 * (x1 * (y2 - y3) * (y4 - y1) - x3 * (y4 - y3) * (y2 - y1) + x2 * (y4 - y2) * (y3 - y1)) +
                  x2 * x3 * (y2 - y3) * (y4 - y1) - x2 * x1 * (y4 - y3) * (y2 - y1) + x3 * x1 * (y4 - y2) * (y3 - y1));

/*             printf("%0.2f\t%0.2f\n", x1, y1);
            printf("%0.2f\t%0.2f\n", x2, y2);
            printf("%0.2f\t%0.2f\n", x3, y3);
            printf("%0.2f\t%0.2f\n", x4, y4);
            printf("%0.2f\t%0.2f\t%0.2f\t%0.2f\n", c1, c2, c3, c4); */

            break;
        }
    
        default:
        {
            c1 = 1.;    c2 = 0.;    c3 = 0.;    c4 = 0.;
            break;
        }
    }

    Interpolate_Coefficient[0] = c1;
    Interpolate_Coefficient[1] = c2;
    Interpolate_Coefficient[2] = c3;
    Interpolate_Coefficient[3] = c4;

    return 1;
}


int FE_Neighbour_Points_Find(double **Grid, double **Coordinate, int ***Node_4Points, int NodeMax, int i, int j)
{
    int k, Sort_Numbers, *Node_Resort_Array, p, Index_R_Min, N, Count;
    int N1, N2, N3, N4, N_Permutation=1, **Extract_4Points, End_Cycles;
    double r_min, **Interpolate_Coordinate, *Interpolate_Values;

    r_min = Grid[0][0];
    Node_Resort_Array = Make1DArrayinteger(NodeMax);
    Interpolate_Coordinate = Make2DArray(5,2);
    Interpolate_Values     = Make1DArray(4);
    Sort_Numbers = 7;

    for(p=0;p<Sort_Numbers;++p)
    {
        N_Permutation = (Sort_Numbers-p)*N_Permutation;
    }

    N_Permutation = N_Permutation/4/3/2;

    for(p=4;p<Sort_Numbers;++p)
    {
        N_Permutation = N_Permutation/(Sort_Numbers-p);
    }


    Extract_4Points = Make2DArrayinteger(N_Permutation,4);

/////////////////////////////////////////////////////////////////
    for(p=0;p<NodeMax;++p)
    {
        Node_Resort_Array[p]=p;
    }

    for(k=0;k<Sort_Numbers;++k)
    {
        Count = 0;
        for(p=k;p<NodeMax;++p)
        {
            //printf("1\n");
            N = Node_Resort_Array[p];
            if(Coordinate[N][3]==1.)
            {
                Count = Count+1;
            }

            if(Count ==1)
            {
                r_min = Grid[N][0];
            }

            if(Grid[N][0]<=r_min && Coordinate[N][3]==1)
            {
                r_min = Grid[N][0];
                Index_R_Min = p;
                //printf("%0.6f\n",r_min);
            }

        }
        N = Node_Resort_Array[Index_R_Min];
        Node_Resort_Array[Index_R_Min] = Node_Resort_Array[k];
        Node_Resort_Array[k]           = N;
        //printf("%0.6f\n",r_min);
    }
    for(k=0;k<Sort_Numbers;++k)
    {
        //printf("%0.6f\t",Grid[Node_Resort_Array[k]][0]);
    }
    //printf("\n");

////////////////////////////////////////////////////////
    k = 0;
    for(N1=0;N1<=Sort_Numbers-4;++N1)
    {
        for(N2=N1+1;N2<=Sort_Numbers-3;++N2)
        {
            for(N3=N2+1;N3<=Sort_Numbers-2;++N3)
            {
                for(N4=N3+1;N4<=Sort_Numbers-1;++N4)
                {
                    Extract_4Points[k][0] = N1;
                    Extract_4Points[k][1] = N2;
                    Extract_4Points[k][2] = N3;
                    Extract_4Points[k][3] = N4;
                    k = k+1;
                }
            }
        }
    }



    for(k=0;k<N_Permutation;++k)
    {
        N1 = Node_Resort_Array[Extract_4Points[k][0]];
        N2 = Node_Resort_Array[Extract_4Points[k][1]];
        N3 = Node_Resort_Array[Extract_4Points[k][2]];
        N4 = Node_Resort_Array[Extract_4Points[k][3]];

        Interpolate_Coordinate[1][0]= Grid[N1][1];    Interpolate_Coordinate[1][1]= Grid[N1][2];
        Interpolate_Coordinate[2][0]= Grid[N2][1];    Interpolate_Coordinate[2][1]= Grid[N2][2];
        Interpolate_Coordinate[3][0]= Grid[N3][1];    Interpolate_Coordinate[3][1]= Grid[N3][2];
        Interpolate_Coordinate[4][0]= Grid[N4][1];    Interpolate_Coordinate[4][1]= Grid[N4][2];

        Calculate_Interpolation(Interpolate_Coordinate,Interpolate_Values);
        End_Cycles = 1;
        for(p=0;p<4;++p)
        {
            if(Interpolate_Values[p]<0 || Interpolate_Values[p]>1)
            {
                End_Cycles = 0;
            }
        }

        if(End_Cycles)
        {
            break;
        }

    }


//////////////////////////////////////////////////////////////
    Node_4Points[i][j][4] = 0;
    Node_4Points[i][j][0] = N1; Node_4Points[i][j][1] = N2;
    Node_4Points[i][j][2] = N3; Node_4Points[i][j][3] = N4;

    if(k>30)
    {
        printf("k=%d\n",k);
        Node_4Points[i][j][4] = 1;
        N1 = Node_Resort_Array[0];  N2 = Node_Resort_Array[1];
        N3 = Node_Resort_Array[2];  N4 = Node_Resort_Array[3];

        Node_4Points[i][j][0] = N1; Node_4Points[i][j][1] = N2;
        Node_4Points[i][j][2] = N3; Node_4Points[i][j][3] = N4;
    }

    free(Node_Resort_Array);
    free2DArray(Interpolate_Coordinate,5);
    free(Interpolate_Values);
    free2DArrayinteger(Extract_4Points,N_Permutation);

    return 1;

}

int ReadAnsys(Input_Parameter_NMD *IPN, double ****CST)
{
    int i, j, ss, Count=0, Count_Min, NodeMax=0, Node, m, n, l;
    int xc, yc, k, ***Node_4Points, N1, N2, N3, N4;
    double **Stress, x, y, z, thx, thy, thz, **Coordinate, **Grid;
    double x1, x2, x3, x4, y1, y2, y3, y4;
    double r, CST1, CST2, CST3, CST4;
    FILE *fp, *fp_Coordinate;
    char s[1000], s_SXZ[1000], s_5Stars[1000], s_MINIMUN[1000], s_THZX[1000], s_NODE[1000];

    m=2*IPN->xn+3;   n=2*IPN->yn+3;

    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l  = 2*IPN->zn+3;
    }
    xc=(m-1)/2; yc=(n-1)/2;

    sprintf(s_THZX,"THZX");
    sprintf(s_NODE,"NODE");
    fp_Coordinate=fopen("NLIST.lis","r+");

    while(!feof(fp_Coordinate))
    {
        ss=fscanf(fp_Coordinate,"%s",s);

        if(strcmp(s,s_THZX)==0)
        {
            //fscanf(fp,)
            while(!feof(fp_Coordinate))
            {
                ss=fscanf(fp_Coordinate,"%s",s);
                //printf("%s\n",s);
                if(strcmp(s,s_NODE)==0 || ss==EOF)
                {
                    break;
                }else
                {
                    Node=atoi(s);
                    fscanf(fp_Coordinate,"%lf%lf%lf%lf%lf%lf",&x,&y,&z,&thx,&thy,&thz);
                    if(Node>=NodeMax)
                    {
                        NodeMax=Node;
                    }

                }
            }

        }
        //printf("%s\n",s);
        if(ss==EOF)
        {
            break;
        }
    }
    fclose(fp_Coordinate);
    printf("%d\n",NodeMax);
///////////////////////////////////////////////////////////////
    Coordinate=Make2DArray(NodeMax,4);
    for(k=0;k<NodeMax;++k)
    {
        Coordinate[k][3]=0;
    }

    fp_Coordinate=fopen("NLIST.lis","r+");
    while(!feof(fp_Coordinate))
    {
        ss=fscanf(fp_Coordinate,"%s",s);

        if(strcmp(s,s_THZX)==0)
        {
            //fscanf(fp,)
            while(!feof(fp_Coordinate))
            {
                ss=fscanf(fp_Coordinate,"%s",s);
                //printf("%s\n",s);
                if(strcmp(s,s_NODE)==0 || ss==EOF)
                {
                    break;
                }else
                {
                    Node=atoi(s);
                    fscanf(fp_Coordinate,"%lf%lf%lf%lf%lf%lf",&x,&y,&z,&thx,&thy,&thz);
                    Coordinate[Node-1][0]=x;
                    Coordinate[Node-1][1]=y;
                    Coordinate[Node-1][2]=z;
                    Coordinate[Node-1][3]=0.;
                    //printf("%d\t%0.2f\t%0.2f\t%0.2f\n",Node-1,x,y,z);
               }
            }

        }
        //printf("%s\n",s);
        if(ss==EOF)
        {
            break;
        }
    }
    fclose(fp_Coordinate);

////////////////////////////////////////////////////

    Stress=Make2DArray(NodeMax,6);
    fp=fopen("PRNSOL.lis","r+");
    sprintf(s_SXZ,"SXZ");
    while(!feof(fp))
    {
        ss=fscanf(fp,"%s",s);

        if(strcmp(s,s_SXZ)==0)
        {
            sprintf(s_5Stars,"*****");
            sprintf(s_MINIMUN,"MINIMUM");
            while(!feof(fp))
            {
                fscanf(fp,"%s",s);
                //printf("%s\n",s);
                if(strcmp(s,s_5Stars)==0 || strcmp(s,s_MINIMUN)==0)
                {
                    break;
                }else
                {
                    Node=atoi(s);
                    Coordinate[Node-1][3]=1.;
                    fscanf(fp,"%lf%lf%lf%lf%lf%lf",&Stress[Node-1][0],&Stress[Node-1][1],&Stress[Node-1][2],
                                                   &Stress[Node-1][3],&Stress[Node-1][5],&Stress[Node-1][4]);
                    //printf("%0.2f\t%0.2f\t%0.2f\n",Stress[Node-1][0],Stress[Node-1][1],Stress[Node-1][3]);

                }
            }

        }
        if(ss==EOF)
        {
            break;
        }
    }
    fclose(fp);

////////////////////////////////////////////////////////////////
    double Grid_x_max, Grid_x_min, Grid_y_max, Grid_y_min;

    Grid_x_max=-fabs(Coordinate[0][0]*10);    Grid_x_min=fabs(Coordinate[0][0]*10);
    Grid_y_max=-fabs(Coordinate[0][1]*10);    Grid_y_min=fabs(Coordinate[0][1]*10);

    for(k=0;k<NodeMax;k++)
    {
        if(Coordinate[k][3]!=0.)
        {
            if(Coordinate[k][0]>=Grid_x_max)
            {
                Grid_x_max=Coordinate[k][0];
            }

            if(Coordinate[k][0]<=Grid_x_min)
            {
                Grid_x_min=Coordinate[k][0];
            }

            if(Coordinate[k][1]>=Grid_y_max)
            {
                Grid_y_max=Coordinate[k][1];
            }

            if(Coordinate[k][1]<=Grid_y_min)
            {
                Grid_y_min=Coordinate[k][1];
            }
        }
        //printf("%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",k,Grid_x_max,Grid_x_min,Grid_y_max,Grid_y_min);
    }

    printf("%0.2f\t%0.2f\t%0.2f\t%0.2f\n",Grid_x_max,Grid_x_min,Grid_y_max,Grid_y_min);

    double x_interval, y_interval;
    x_interval=(Grid_x_max-Grid_x_min)/m;
    y_interval=(Grid_y_max-Grid_y_min)/n;
///////////////////////////////////////////////////
    //Theta_Grid  =Make3DArray(m,n,NodeMax);
    Node_4Points=Make3DArrayinteger(m,n,5);
    Grid=Make2DArray(NodeMax,3);

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            x =(i-xc)*x_interval;
            y =(j-yc)*y_interval;
            for(k=0;k<NodeMax;k++)
            {
                if(Coordinate[k][3]==1.)
                {
                    r=sqrt((x-Coordinate[k][0])*(x-Coordinate[k][0])
                          +(y-Coordinate[k][1])*(y-Coordinate[k][1]));
                    Grid[k][0]=r;
                    Grid[k][1]=Coordinate[k][0]-x;
                    Grid[k][2]=Coordinate[k][1]-y;

                    //Theta_Grid[i][j][k]=atan2(y-Coordinate[k][1],x-Coordinate[k][0]);
                }

            }

            FE_Neighbour_Points_Find(Grid,Coordinate,Node_4Points,NodeMax,i,j);


            //N1 = Node_4Points[i][j][0]; N2 = Node_4Points[i][j][1];
            //N3 = Node_4Points[i][j][2]; N4 = Node_4Points[i][j][3];
            //printf("%d\t%d\t%d\t%d\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\n",Node_4Points[i][j][0],
                 //  Node_4Points[i][j][1],Node_4Points[i][j][2],Node_4Points[i][j][3],Coordinate[N1][0],
                   // Coordinate[N1][1],Coordinate[N2][0],Coordinate[N2][1],
                   // Coordinate[N3][0],Coordinate[N3][1],Coordinate[N4][0],Coordinate[N4][1],x,y);

            for(k=0;k<6;k++)
            {
                N1=Node_4Points[i][j][0];
                N2=Node_4Points[i][j][1];
                N3=Node_4Points[i][j][2];
                N4=Node_4Points[i][j][3];

                CST1=Stress[N1][k];     x1=Grid[N1][1];     y1=Grid[N1][2];
                CST2=Stress[N2][k];     x2=Grid[N2][1];     y2=Grid[N2][2];
                CST3=Stress[N3][k];     x3=Grid[N3][1];     y3=Grid[N3][2];
                CST4=Stress[N4][k];     x4=Grid[N4][1];     y4=Grid[N4][2];



                if(fabs(x1)<=(x_interval)/4 && fabs(y1)<=(y_interval)/4)
                {
                    CST[k][i][j][0]=CST1;
                }else if(fabs(x2)<=(x_interval)/4 && fabs(y2)<=(y_interval)/4)
                {
                    CST[k][i][j][0]=CST2;
                }else if(fabs(x3)<=(x_interval)/4 && fabs(y3)<=(y_interval)/4)
                {
                    CST[k][i][j][0]=CST3;
                }else if(fabs(x4)<=(x_interval)/4 && fabs(y4)<=(y_interval)/4)
                {
                    CST[k][i][j][0]=CST4;
                }else if(Node_4Points[i][j][4]==1)
                {
                    CST[k][i][j][0]=CST1;
                }else
                {
                    CST[k][i][j][0]=CST1*(x2*x4*y3*(y2-y4)+x2*x3*y4*(y3-y2)+x3*x4*y2*(y4-y3))/
                                         (x1*(x4*(y2-y3)*(y1-y4)-x3*(y1-y3)*(y2-y4)+x2*(y1-y2)*(y3-y4))+
                                          x2*x3*(y2-y3)*(y1-y4)-x2*x4*(y1-y3)*(y2-y4)+x3*x4*(y1-y2)*(y3-y4))+
                                    CST2*(x1*x4*y3*(y1-y4)+x1*x3*y4*(y3-y1)+x3*x4*y1*(y4-y3))/
                                         (x2*(x4*(y1-y3)*(y2-y4)-x3*(y2-y3)*(y1-y4)+x1*(y2-y1)*(y3-y4))+
                                          x1*x3*(y1-y3)*(y2-y4)-x1*x4*(y2-y3)*(y1-y4)+x3*x4*(y2-y1)*(y3-y4))+
                                    CST3*(x2*x4*y1*(y2-y4)+x2*x1*y4*(y1-y2)+x1*x4*y2*(y4-y1))/
                                         (x3*(x4*(y2-y1)*(y3-y4)-x1*(y3-y1)*(y2-y4)+x2*(y3-y2)*(y1-y4))+
                                          x2*x1*(y2-y1)*(y3-y4)-x2*x4*(y3-y1)*(y2-y4)+x1*x4*(y3-y2)*(y1-y4))+
                                    CST4*(x2*x1*y3*(y2-y1)+x2*x3*y1*(y3-y2)+x3*x1*y2*(y1-y3))/
                                         (x4*(x1*(y2-y3)*(y4-y1)-x3*(y4-y3)*(y2-y1)+x2*(y4-y2)*(y3-y1))+
                                          x2*x3*(y2-y3)*(y4-y1)-x2*x1*(y4-y3)*(y2-y1)+x3*x1*(y4-y2)*(y3-y1));
                }


                //printf("%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",x1,y1,x2,y2,x3,y3,x4,y4);
                if(CST[0][i][j][0]<0)
                {
                    printf("%d\t%d\t%d\t%d\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\n",N1,N2,N3,N4,CST1,CST2,CST3,CST4,CST[k][i][j][0]);
                }

                CST[k][i][j][0]=CST[k][i][j][0]*1.e6;

                //printf("%0.8f\n",CST[i][j][k]);
            }
        }
        printf("%d\n",i);
    }

    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            CST[6][i][j][0]=1.;
        }
    }

    FILE *fp_S11, *fp_S22, *fp_S33, *fp_S12, *fp_S13, *fp_S23;
    char *filename, *filename_S11, *filename_S22, *filename_S33, *filename_S12, *filename_S13, *filename_S23;
    filename_S11 = Make1DString(500); filename_S22 = Make1DString(500);
    filename_S33 = Make1DString(500); filename_S12 = Make1DString(500);
    filename_S13 = Make1DString(500); filename_S23 = Make1DString(500);
    filename     = Make1DString(500);

    sprintf(filename,"Elastic Field of ANSYS");
    mkdir(filename);
    sprintf(filename_S11,"Elastic Field of ANSYS\\S11.dat");
    sprintf(filename_S22,"Elastic Field of ANSYS\\S22.dat");
    sprintf(filename_S33,"Elastic Field of ANSYS\\S33.dat");
    sprintf(filename_S12,"Elastic Field of ANSYS\\S12.dat");
    sprintf(filename_S13,"Elastic Field of ANSYS\\S13.dat");
    sprintf(filename_S23,"Elastic Field of ANSYS\\S23.dat");
    fp_S11=fopen(filename_S11,"w+"); fp_S22=fopen(filename_S22,"w+");
    fp_S33=fopen(filename_S33,"w+"); fp_S12=fopen(filename_S12,"w+");
    fp_S13=fopen(filename_S13,"w+"); fp_S23=fopen(filename_S23,"w+");
    for(i=0;i<m;i++)
    {
        for(j=0;j<n;j++)
        {
            if(j<n-1)
            {
                fprintf(fp_S11,"%0.8f\t",CST[0][i][j][0]);
                fprintf(fp_S22,"%0.8f\t",CST[1][i][j][0]);
                fprintf(fp_S33,"%0.8f\t",CST[2][i][j][0]);
                fprintf(fp_S12,"%0.8f\t",CST[3][i][j][0]);
                fprintf(fp_S13,"%0.8f\t",CST[4][i][j][0]);
                fprintf(fp_S23,"%0.8f\t",CST[5][i][j][0]);
            }else
            {
                fprintf(fp_S11,"%0.8f",CST[0][i][j][0]);
                fprintf(fp_S22,"%0.8f",CST[1][i][j][0]);
                fprintf(fp_S33,"%0.8f",CST[2][i][j][0]);
                fprintf(fp_S12,"%0.8f",CST[3][i][j][0]);
                fprintf(fp_S13,"%0.8f",CST[4][i][j][0]);
                fprintf(fp_S23,"%0.8f",CST[5][i][j][0]);
            }

        }
        fprintf(fp_S11,"\n");
        fprintf(fp_S22,"\n");
        fprintf(fp_S33,"\n");
        fprintf(fp_S12,"\n");
        fprintf(fp_S13,"\n");
        fprintf(fp_S23,"\n");
    }
    fclose(fp_S11);
    fclose(fp_S22);
    fclose(fp_S33);
    fclose(fp_S12);
    fclose(fp_S13);
    fclose(fp_S23);

    free(filename);
    free(filename_S11);
    free(filename_S22);
    free(filename_S33);
    free(filename_S12);
    free(filename_S13);
    free(filename_S23);
    free2DArray(Stress,NodeMax);
    free2DArray(Coordinate,NodeMax);
    free2DArray(Grid,NodeMax);
    free3DArrayinteger(Node_4Points,m,n);



    return 1;

}

int Write_Image(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****Mx_Image, 
                double ****My_Image, double ****Mz_Image, double **MxV_Image, double **MyV_Image, double **MzV_Image)
{
    int i, j, k, v, m, n, l, zc, Version, N_Image, i_Image, Inplane_Density, D;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV, *fp_M_Inplane;
    double x, y, Mx_Offset, My_Offset, x_max, x_min, y_max, y_min;
    double Mx0, My0, M_Module, M_Appropriate, Ratio, M_min, Scale;
    char *filename_Dir,*filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV, *filename_M_Inplane;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);
    filename_Dir = Make1DString(500);
    filename_M_Inplane  = Make1DString(500);
    Version = IPN->Version;
    N_Image = IPN->N_Image;

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

    Inplane_Density = IPN->Inplane_Density;
    if(Inplane_Density == 0)
    {
        Inplane_Density = 24;
    }
    D = m/Inplane_Density;

    if(D<=1)
    {
        D=1;
    }

    Scale = 1.;
    M_Appropriate = D*fabs(IPN->dxy[0])*0.8;
    M_min         = D*fabs(IPN->dxy[0])*0.02;

    x_max = fabs(IPN->dxy[0])*(m-1)/2;   x_min = -fabs(IPN->dxy[0])*(m-1)/2;
    y_max = fabs(IPN->dxy[1])*(n-1)/2;   y_min = -fabs(IPN->dxy[1])*(n-1)/2;
    
    sprintf(filename_Dir,"Images");
    mkdir(filename_Dir);
    sprintf(filename_Dir,"Images\\Data");
    mkdir(filename_Dir);
    sprintf(filename_Dir,"Images\\For Printing");
    mkdir(filename_Dir);
    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        sprintf(filename_Dir,"Images\\Data\\Image %d",i_Image+1);
        mkdir(filename_Dir);
        for (k = 0; k < l; ++k)
        {
            sprintf(filename_Dir, "Images\\Data\\Image %d\\z=%d", i_Image + 1, k - zc);
            mkdir(filename_Dir);
        }

        sprintf(filename_Dir, "Images\\For Printing\\Image %d", i_Image + 1);
        mkdir(filename_Dir);
        for (k = 0; k < l; ++k)
        {
            sprintf(filename_Dir, "Images\\For Printing\\Image %d\\z=%d", i_Image + 1, k - zc);
            mkdir(filename_Dir);
        }
    }

    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        sprintf(filename_MxV, "Images\\Data\\Image %d\\MxV.dat", i_Image+1);
        sprintf(filename_MyV, "Images\\Data\\Image %d\\MyV.dat", i_Image+1);
        sprintf(filename_MzV, "Images\\Data\\Image %d\\MzV.dat", i_Image+1);
        fp_MxV = fopen(filename_MxV,"w+");
        fp_MyV = fopen(filename_MyV,"w+");
        fp_MzV = fopen(filename_MzV,"w+");

        for(v=0;v<BS->Virtual_Points;++v)
        {
            fprintf(fp_MxV,"%0.8f\t",MxV_Image[i_Image][v]);
            fprintf(fp_MyV,"%0.8f\t",MyV_Image[i_Image][v]);
            fprintf(fp_MzV,"%0.8f\t",MzV_Image[i_Image][v]);
        }
        fclose(fp_MxV);
        fclose(fp_MyV);
        fclose(fp_MzV);
    }


    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        for(k=0;k<l;k++)
        {
            sprintf(filename_Mx,"Images\\Data\\Image %d\\z=%d\\Mx.dat",i_Image+1,k-zc);
            sprintf(filename_My,"Images\\Data\\Image %d\\z=%d\\My.dat",i_Image+1,k-zc);
            sprintf(filename_Mz,"Images\\Data\\Image %d\\z=%d\\Mz.dat",i_Image+1,k-zc);
            
            fp_Mx = fopen(filename_Mx,"w+");
            fp_My = fopen(filename_My,"w+");
            fp_Mz = fopen(filename_Mz,"w+");

            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    if(j<n-1)
                    {
                        fprintf(fp_Mx,"%0.8f\t",Mx_Image[i_Image][i][j][k]);
                        fprintf(fp_My,"%0.8f\t",My_Image[i_Image][i][j][k]);
                        fprintf(fp_Mz,"%0.8f\t",Mz_Image[i_Image][i][j][k]);
                    }else
                    {
                        fprintf(fp_Mx,"%0.8f",Mx_Image[i_Image][i][j][k]);
                        fprintf(fp_My,"%0.8f",My_Image[i_Image][i][j][k]);
                        fprintf(fp_Mz,"%0.8f",Mz_Image[i_Image][i][j][k]);
                    }
                }
                if(i<m-1)
                {
                    fprintf(fp_Mx,"\n");
                    fprintf(fp_My,"\n");
                    fprintf(fp_Mz,"\n");
                }
            }
            fclose(fp_Mx);
            fclose(fp_My);
            fclose(fp_Mz);
        }
    }
//////////////////////////////////////////////////////////////////////////////////////
    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        for(k=0;k<l;k++)
        {
            sprintf(filename_Mx,"Images\\For Printing\\Image %d\\z=%d\\Mx.dat",i_Image+1,k-zc);
            sprintf(filename_My,"Images\\For Printing\\Image %d\\z=%d\\My.dat",i_Image+1,k-zc);
            sprintf(filename_Mz,"Images\\For Printing\\Image %d\\z=%d\\Mz.dat",i_Image+1,k-zc);
            sprintf(filename_M_Inplane,"Images\\For Printing\\Image %d\\z=%d\\M_InPlane.dat",i_Image+1,k-zc);
            
            fp_M_Inplane = fopen(filename_M_Inplane,"w+");
            fp_Mx = fopen(filename_Mx,"w+");
            fp_My = fopen(filename_My,"w+");
            fp_Mz = fopen(filename_Mz,"w+");

            fprintf(fp_M_Inplane,"y\tx\tMy_Offset\tMx_Offset\tMy\tMx\n");
            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    if(j<n-1)
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Mx,"%0.8f\t",Mx_Image[i_Image][i][j][k]);
                            fprintf(fp_My,"%0.8f\t",My_Image[i_Image][i][j][k]);
                            fprintf(fp_Mz,"%0.8f\t",Mz_Image[i_Image][i][j][k]);
                        }else
                        {
                            fprintf(fp_Mx,"\t");
                            fprintf(fp_My,"\t");
                            fprintf(fp_Mz,"\t");
                        }
                    }else
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            fprintf(fp_Mx,"%0.8f",Mx_Image[i_Image][i][j][k]);
                            fprintf(fp_My,"%0.8f",My_Image[i_Image][i][j][k]);
                            fprintf(fp_Mz,"%0.8f",Mz_Image[i_Image][i][j][k]);
                        }
                    }

                    if((i-IPN->xn-1)%D==0 && (j-IPN->yn-1)%D==0)
                    {
                        if(BS->Boundary[i][j][k]!=0)
                        {
                            x=(i-IPN->xn-1)*IPN->dxy[0];
                            y=(j-IPN->yn-1)*IPN->dxy[1];
                            Mx0 = Mx_Image[i_Image][i][j][k] * Scale;
                            My0 = My_Image[i_Image][i][j][k] * Scale;

                            M_Module = sqrt(Mx0*Mx0+My0*My0);
                            if(M_Module>M_Appropriate)
                            {
                                Ratio = M_Appropriate/M_Module;
                                Mx0 = Mx0*Ratio;
                                My0 = My0*Ratio;
                            }else if(M_Module<M_min)
                            {
                                Mx0 = 0;
                                My0 = 0;
                            }

                            Mx_Offset = Mx0+x;
                            My_Offset = My0+y;
                            if(Mx_Offset>=x_max || Mx_Offset<=x_min || My_Offset>=y_max || My_Offset<=y_min)
                            {
                                Mx_Offset = x;
                                My_Offset = y;
                            }
                        fprintf(fp_M_Inplane,"%0.2f\t%0.2f\t%0.8f\t%0.8f\t%0.8f\t%0.8f\n",y,x,My_Offset,Mx_Offset,My_Image[i_Image][i][j][k],Mx_Image[i_Image][i][j][k]);
                        }
                    }
                }

                if(i<m-1)
                {
                    fprintf(fp_Mx,"\n");
                    fprintf(fp_My,"\n");
                    fprintf(fp_Mz,"\n");
                }
            }
            fclose(fp_M_Inplane);
            fclose(fp_Mx);
            fclose(fp_My);
            fclose(fp_Mz);
        }
    }
//////////////////////////////////////////////////////////////////////////////
    free(filename_M_Inplane);
    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    free(filename_Dir);
    return 1;
}

int Read_Image(Input_Parameter_NMD *IPN, Boundary_Shape *BS, double ****Mx_Image, 
                double ****My_Image, double ****Mz_Image, double **MxV_Image, double **MyV_Image, double **MzV_Image)
{
    int i, j, k, v, m, n, l, zc, N_Image, i_Image;
    FILE *fp_Mx, *fp_My, *fp_Mz, *fp_MxV, *fp_MyV, *fp_MzV;
    char *filename_Mx, *filename_My, *filename_Mz, *filename_MxV, *filename_MyV, *filename_MzV;

    filename_Mx  = Make1DString(500);
    filename_My  = Make1DString(500);
    filename_Mz  = Make1DString(500);
    filename_MxV = Make1DString(500);
    filename_MyV = Make1DString(500);
    filename_MzV = Make1DString(500);
    N_Image = IPN->N_Image;

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

    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        sprintf(filename_MxV, "Images\\Data\\Image %d\\MxV.dat", i_Image+1);
        sprintf(filename_MyV, "Images\\Data\\Image %d\\MyV.dat", i_Image+1);
        sprintf(filename_MzV, "Images\\Data\\Image %d\\MzV.dat", i_Image+1);
        fp_MxV = fopen(filename_MxV,"r+");
        fp_MyV = fopen(filename_MyV,"r+");
        fp_MzV = fopen(filename_MzV,"r+");

        for(v=0;v<BS->Virtual_Points;++v)
        {
            fscanf(fp_MxV,"%lf",&MxV_Image[i_Image][v]);
            fscanf(fp_MyV,"%lf",&MyV_Image[i_Image][v]);
            fscanf(fp_MzV,"%lf",&MzV_Image[i_Image][v]);
        }
        fclose(fp_MxV);
        fclose(fp_MyV);
        fclose(fp_MzV);
    }

    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        for(k=0;k<l;k++)
        {
            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    Mx_Image[i_Image][i][j][k] = 0;
                    My_Image[i_Image][i][j][k] = 0;
                    Mz_Image[i_Image][i][j][k] = 0;
                }
            }
        }
    }

    for (i_Image = 0; i_Image < N_Image; ++i_Image)
    {
        for(k=0;k<l;k++)
        {
            sprintf(filename_Mx,"Images\\Data\\Image %d\\z=%d\\Mx.dat",i_Image+1,k-zc);
            sprintf(filename_My,"Images\\Data\\Image %d\\z=%d\\My.dat",i_Image+1,k-zc);
            sprintf(filename_Mz,"Images\\Data\\Image %d\\z=%d\\Mz.dat",i_Image+1,k-zc);
            fp_Mx = fopen(filename_Mx,"r+");
            fp_My = fopen(filename_My,"r+");
            fp_Mz = fopen(filename_Mz,"r+");

            for(i=0;i<m;i++)
            {
                for(j=0;j<n;j++)
                {
                    fscanf(fp_Mx,"%lf",&Mx_Image[i_Image][i][j][k]);
                    fscanf(fp_My,"%lf",&My_Image[i_Image][i][j][k]);
                    fscanf(fp_Mz,"%lf",&Mz_Image[i_Image][i][j][k]);
                }
            }
            fclose(fp_Mx);
            fclose(fp_My);
            fclose(fp_Mz);
        }
    }

    free(filename_Mx);
    free(filename_My);
    free(filename_Mz);
    free(filename_MxV);
    free(filename_MyV);
    free(filename_MzV);
    return 1;
}

int Write_Image_Result(double *Energy_Image, double *SKN_Image, int N_Image)
{
    int i_Image;
    double Energy_Difference_Max=0, Energy_Image_Initial;
    Energy_Image_Initial = Energy_Image[0];
    for (i_Image = 0; i_Image < N_Image + 2; ++i_Image)
    {
        if(fabs(Energy_Image[i_Image]-Energy_Image_Initial)>=Energy_Difference_Max)
        {
            Energy_Difference_Max = Energy_Image[i_Image] - Energy_Image_Initial;
        }
    }

    FILE *fp;
    char *filename;

    filename  = Make1DString(500);
    //sprintf(filename, "Images_Result.dat");
    sprintf(filename, "Images_Result_%d.dat", i_Angle_Rotation);
    fp = fopen(filename,"w+");

    fprintf(fp,"Energy_Difference_Max= %0.8f\n", Energy_Difference_Max);
    fprintf(fp,"---------------------------------------------\n");
    fprintf(fp,"Image\t Energy\t\t SKN\n");
    for (i_Image = 0; i_Image < N_Image + 2; ++i_Image)
    {
        fprintf(fp, "%d\t %0.8f\t %0.8f\n", i_Image, Energy_Image[i_Image], SKN_Image[i_Image]);
    }

    fclose(fp);
    free(filename);
    return 1;
}