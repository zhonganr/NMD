#include <GUI.h>
#include <ui_gui.h>
#include <Variable.h>
#include <QDebug>
#include <windows.h>
#include <Modify_MC.h>
#include <Modify_NL.h>
#include <LLG.h>
#include <GNEB.h>
#include <InOut.h>
#include <Initial_State.h>
#include <Elastic_Field.h>
#include <BoundaryShape.h>
#include <Array.h>
#include <pthread.h>

extern Input_Parameter_NMD *IPN;
extern double ***Mx_Global, ***My_Global, ***Mz_Global, *MxV_Global, *MyV_Global, *MzV_Global;
extern double ***Mx_OPGL, ***My_OPGL, ***Mz_OPGL, *MxV_OPGL, *MyV_OPGL, *MzV_OPGL;
GUI_Parameter GPM;
Pthread_State_Run_Parameter PSRP;
extern int N_Energy_Terms;

void GUI::SetPage_CalculateMode(int Index_Page_CalculateMode)
{
    int n_Count = ui->stackedWidget_Calculate_Mode->count();
    if (Index_Page_CalculateMode >= n_Count)
    {
        Index_Page_CalculateMode = 0;
    }
    ui->stackedWidget_Calculate_Mode->setCurrentIndex(Index_Page_CalculateMode);
}

void GUI::SetPage_ElasticField(int Index_Page_ElasticField)
{
    int n_Count = ui->stackedWidget_ElasticField->count();
    if (Index_Page_ElasticField >= n_Count)
    {
        Index_Page_ElasticField = 0;
    }
    ui->stackedWidget_ElasticField->setCurrentIndex(Index_Page_ElasticField);
}

void GUI::GUI_Initial(Input_Parameter_NMD *IPN)
{
    int n_Count, Index, i;
    Energy_Terms ET;
    ui->TextOut_Xn->setText(QString::number(IPN->xn, 10));
    ui->TextOut_Yn->setText(QString::number(IPN->yn, 10));
    ui->TextOut_Zn->setText(QString::number(IPN->zn, 10));
    ui->TextOut_dx->setText(QString::number(IPN->dxy[0],'g',6));
    ui->TextOut_dy->setText(QString::number(IPN->dxy[1],'g',6));
    ui->TextOut_dz->setText(QString::number(IPN->dxy[2],'g',6));

    Index = IPN->Boundary_Shape_Type;
    n_Count = ui->comboBox_Shape->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_Shape->setCurrentIndex(Index);

    Index = IPN->Boundary_Type_x;
    n_Count = ui->comboBox_ShapeType_x->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_ShapeType_x->setCurrentIndex(Index);

    Index = IPN->Boundary_Type_y;
    n_Count = ui->comboBox_ShapeType_y->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_ShapeType_y->setCurrentIndex(Index);

    Index = IPN->Boundary_Type_z;
    n_Count = ui->comboBox_ShapeType_z->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_ShapeType_z->setCurrentIndex(Index);

    Index = IPN->CalculateMode;
    n_Count = ui->comboBox_Calculate_Mode->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_Calculate_Mode->setCurrentIndex(Index);
    SetPage_CalculateMode(Index);
    ui->TextOut_MPx_SOT->setText(QString::number(IPN->MpX,'g',6));
    ui->TextOut_MPy_SOT->setText(QString::number(IPN->MpY,'g',6));
    ui->TextOut_MPz_SOT->setText(QString::number(IPN->MpZ,'g',6));
    ui->TextOut_CurrentDensity_SOT->setText(QString::number(IPN->CurrentDensity,'g',6));
    ui->TextOut_dT_SOT->setText(QString::number(IPN->dT,'g',6));

    ui->TextOut_MPx_STT->setText(QString::number(IPN->MpX,'g',6));
    ui->TextOut_MPy_STT->setText(QString::number(IPN->MpY,'g',6));
    ui->TextOut_MPz_STT->setText(QString::number(IPN->MpZ,'g',6));
    ui->TextOut_CurrentDensity_STT->setText(QString::number(IPN->CurrentDensity,'g',6));
    ui->TextOut_dT_STT->setText(QString::number(IPN->dT,'g',6));

    Index = IPN->Continuity_Modify_Mode;
    n_Count = ui->comboBox_Continue_Mode->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_Continue_Mode->setCurrentIndex(Index);
    ui->TextOut_Continue_Mode_Factor->setText(QString::number(IPN->Continuity_Modify_Coefficient,'g',6));

    Index = IPN->InitialMode + 1;
    n_Count = ui->comboBox_InitialMode->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_InitialMode->setCurrentIndex(Index);
    
    Index = IPN->ElasticFieldMode;
    n_Count = ui->comboBox_ElasticField->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_ElasticField->setCurrentIndex(Index);
    SetPage_ElasticField(Index);
    ui->TextOut_Strain_e11->setText(QString::number(IPN->CST_Uniaxial[0],'g',6));
    ui->TextOut_Strain_e22->setText(QString::number(IPN->CST_Uniaxial[1],'g',6));
    ui->TextOut_Strain_e33->setText(QString::number(IPN->CST_Uniaxial[2],'g',6));
    ui->TextOut_Strain_e12->setText(QString::number(IPN->CST_Uniaxial[3],'g',6));
    ui->TextOut_Strain_e13->setText(QString::number(IPN->CST_Uniaxial[4],'g',6));
    ui->TextOut_Strain_e23->setText(QString::number(IPN->CST_Uniaxial[5],'g',6));

    ui->TextOut_Stress_e11->setText(QString::number(IPN->CST_Uniaxial[1],'g',6));
    ui->TextOut_Stress_e22->setText(QString::number(IPN->CST_Uniaxial[2],'g',6));
    ui->TextOut_Stress_e33->setText(QString::number(IPN->CST_Uniaxial[3],'g',6));
    ui->TextOut_Stress_e12->setText(QString::number(IPN->CST_Uniaxial[4],'g',6));
    ui->TextOut_Stress_e13->setText(QString::number(IPN->CST_Uniaxial[5],'g',6));
    ui->TextOut_Stress_e23->setText(QString::number(IPN->CST_Uniaxial[6],'g',6));

    ui->TextOut_Edge_Dislocation_b->setText(QString::number(IPN->b,'g',2));
    ui->TextOut_Screw_Dislocation_b->setText(QString::number(IPN->b,'g',2));
    ui->TextOut_Crack_a->setText(QString::number(IPN->a,'g',2));
    ui->TextOut_Concentrated_Force_FC->setText(QString::number(IPN->FC,'g',2));
    ui->TextOut_Beam_Uniform_Force_FU->setText(QString::number(IPN->FU_Beam,'g',2));
    ui->TextOut_Beam_Concentrated_Force_FC->setText(QString::number(IPN->FC_Beam,'g',2));
    ui->TextOut_Beam_Bending_M->setText(QString::number(IPN->B_M,'g',2));
    ui->TextOut_Annulus_Torsion_T->setText(QString::number(IPN->B_M,'g',2));
    ui->TextOut_Ellipse_Torsion_T->setText(QString::number(IPN->B_M,'g',2));
    ui->TextOut_Read_Ansys_SF->setText(QString::number(IPN->ElasticFieldScaleFactor,'g',2));

    ui->TextOut_NCycle->setText(QString::number(IPN->NCycles,10));
    ui->TextOut_Version->setText(QString::number(IPN->Version,10));
    if (IPN->VersionMode == VM_Plus)
    {
        ui->radioButton_VersionPlus->setChecked(true);
    }else if (IPN->VersionMode == VM_Minus)
    {
        ui->radioButton_VersionMinus->setChecked(true);
    }
    ui->TextOut_ThreadNumber->setText(QString::number(IPN->Threads_Number,10));


    ui->TextOut_H_Amplitude->setText(QString::number(IPN->h,'g',3));
    Index = IPN->MagneticFieldMode;
    n_Count = ui->comboBox_MagneticField->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_MagneticField->setCurrentIndex(Index);
    ui->TextOut_H_Period->setText(QString::number(IPN->MagneticFieldPeriod, 10));

    ui->TextOut_T_Amplitude->setText(QString::number(IPN->t,'g',3));
    Index = IPN->TemperatureFieldMode;
    n_Count = ui->comboBox_TemperatureField->count();
    if (Index >= n_Count)
    {
        Index = 0;
    }
    ui->comboBox_TemperatureField->setCurrentIndex(Index);
    ui->TextOut_T_Period->setText(QString::number(IPN->MagneticFieldPeriod, 10));
    ui->TextOut_N_Image->setText(QString::number(IPN->N_Image, 10));

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch (ET)
        {
            case Exchange:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_Ex->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_Ex->setCheckState(Qt::Unchecked);
                }
                break;
            }

            case DM_Bloch:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_DM_Bloch->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_DM_Bloch->setCheckState(Qt::Unchecked);
                }
                break;
            }
            
            case DM_Neel:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_DM_Neel->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_DM_Neel->setCheckState(Qt::Unchecked);
                }
                break;
            }

            case Zeeman:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_ZM->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_ZM->setCheckState(Qt::Unchecked);
                }
                break;
            }

            case Landau:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_Landau->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_Landau->setCheckState(Qt::Unchecked);
                }
                break;
            }

            case Magnetoelastic:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_ME->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_ME->setCheckState(Qt::Unchecked);
                }
                break;
            }

            case Axial_Anisotropy:
            {
                if (IPN->Energy_Coefficient[i] == 1)
                {
                    ui->checkBox_AA->setCheckState(Qt::Checked);
                }else if (IPN->Energy_Coefficient[i] == 0)
                {
                    ui->checkBox_AA->setCheckState(Qt::Unchecked);
                }
                break;
            }

            default:
            break;
        }
    }
    
    if (IPN->OPGL_Mode == OPGL_Mode_On)
    {
        ui->checkBox_OPGL->setCheckState(Qt::Checked);
    }else if (IPN->OPGL_Mode == OPGL_Mode_Off)
    {
        ui->checkBox_OPGL->setCheckState(Qt::Unchecked);
    }
}

void GUI::Pthread_State_Run_Part_Prepare()
{
    int i, N_Thread, N_Thread_Modify, State_Checkbox;
    QString str;
    //N_Thread = Run thread + Synchronize thread + OPGL thread
    N_Thread = IPN->Threads_Number + 2;
    Energy_Terms ET;

    str = ui->TextOut_Xn->toPlainText();
    IPN->xn = str.toInt();
    str = ui->TextOut_Yn->toPlainText();
    IPN->yn = str.toInt();
    str = ui->TextOut_Zn->toPlainText();
    IPN->zn = str.toInt();

    str = ui->TextOut_dx->toPlainText();
    IPN->dxy[0] = str.toDouble();
    str = ui->TextOut_dy->toPlainText();
    IPN->dxy[1] = str.toDouble();
    str = ui->TextOut_dz->toPlainText();
    IPN->dxy[2] = str.toDouble();

    IPN->Boundary_Shape_Type = ui->comboBox_Shape->currentIndex();
    IPN->Boundary_Type_x = ui->comboBox_ShapeType_x->currentIndex();
    IPN->Boundary_Type_y = ui->comboBox_ShapeType_y->currentIndex();
    IPN->Boundary_Type_z = ui->comboBox_ShapeType_z->currentIndex();
    /////////////////////////////////////////////////////////////////////////////
    str = ui->TextOut_Version->toPlainText();
    IPN->Version = str.toInt();
    PSRP.Version = IPN->Version;
    str = ui->TextOut_ThreadNumber->toPlainText();
    IPN->Threads_Number = str.toInt();
    N_Thread_Modify = IPN->Threads_Number + 2;

    //////////////////////////  Signal Parameter /////////////////////////////////////
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_cond_destroy(&GPM.Pthread_GUI_Cond[i]);
        pthread_mutex_destroy(&GPM.Pthread_GUI_Mutex[i]);
    }

    free(GPM.Pthread_GUI_Cond);
    free(GPM.Pthread_GUI_Mutex);
    GPM.Pthread_GUI_Cond  = (pthread_cond_t *)malloc(sizeof(pthread_cond_t) * N_Thread_Modify);
    GPM.Pthread_GUI_Mutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t) * N_Thread_Modify);

    for (i = 0; i < N_Thread_Modify; ++i)
    {
        pthread_cond_init(&GPM.Pthread_GUI_Cond[i], NULL);
        pthread_mutex_init(&GPM.Pthread_GUI_Mutex[i], 0);
    }

    free2DArrayinteger(GPM.Thread_State, N_Thread);
    GPM.Thread_State = Make2DArrayinteger(N_Thread_Modify, 2);
    /////////////////////////////////////////////////////////////////////////////////////////////
    str = ui->TextOut_NCycle->toPlainText();
    IPN->NCycles = str.toInt();

    IPN->InitialMode = ui->comboBox_InitialMode->currentIndex() - 1;

    IPN->CalculateMode = ui->comboBox_Calculate_Mode->currentIndex();
    Calculate_Mode_GUI CMG;
    CMG = (Calculate_Mode_GUI)IPN->CalculateMode;
    switch (CMG)
    {
        case CMG_LLG_SOT:
        {
            str = ui->TextOut_MPx_SOT->toPlainText();
            IPN->MpX = str.toDouble();
            str = ui->TextOut_MPy_SOT->toPlainText();
            IPN->MpY = str.toDouble();
            str = ui->TextOut_MPz_SOT->toPlainText();
            IPN->MpZ = str.toDouble();
            str = ui->TextOut_CurrentDensity_SOT->toPlainText();
            IPN->CurrentDensity = str.toDouble();
            str = ui->TextOut_dT_SOT->toPlainText();
            IPN->dT = str.toDouble();
            break;
        }

        case CMG_LLG_STT:
        {
            str = ui->TextOut_MPx_STT->toPlainText();
            IPN->MpX = str.toDouble();
            str = ui->TextOut_MPy_STT->toPlainText();
            IPN->MpY = str.toDouble();
            str = ui->TextOut_MPz_STT->toPlainText();
            IPN->MpZ = str.toDouble();
            str = ui->TextOut_CurrentDensity_STT->toPlainText();
            IPN->CurrentDensity = str.toDouble();
            str = ui->TextOut_dT_STT->toPlainText();
            IPN->dT = str.toDouble();
            break;
        }

        default:
        break;
    }

    for (i = 0; i < N_Energy_Terms; ++i)
    {
        ET = (Energy_Terms)i;
        switch (ET)
        {
            case Exchange:
            {
                State_Checkbox = ui->checkBox_Ex->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            case DM_Bloch:
            {
                State_Checkbox = ui->checkBox_DM_Bloch->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }
            
            case DM_Neel:
            {
                State_Checkbox = ui->checkBox_DM_Neel->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            case Zeeman:
            {
                State_Checkbox = ui->checkBox_ZM->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            case Landau:
            {
                State_Checkbox = ui->checkBox_Landau->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            case Magnetoelastic:
            {
                State_Checkbox = ui->checkBox_ME->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            case Axial_Anisotropy:
            {
                State_Checkbox = ui->checkBox_AA->checkState();
                if (State_Checkbox == Qt::Checked)
                {
                    IPN->Energy_Coefficient[i] = 1;
                }else if (State_Checkbox == Qt::Unchecked)
                {
                    IPN->Energy_Coefficient[i] = 0;
                }
                break;
            }

            default:
            break;
        }
    }
    
    State_Checkbox = ui->checkBox_OPGL->checkState();
    if (State_Checkbox == Qt::Checked)
    {
        IPN->OPGL_Mode = OPGL_Mode_On;
    }else if (State_Checkbox == Qt::Unchecked)
    {
        IPN->OPGL_Mode = OPGL_Mode_Off;
    }
    
    IPN->Continuity_Modify_Mode = ui->comboBox_Continue_Mode->currentIndex();
    str = ui->TextOut_Continue_Mode_Factor->toPlainText();
    IPN->Continuity_Modify_Coefficient = str.toDouble();

    IPN->MagneticFieldMode = ui->comboBox_MagneticField->currentIndex();
    str = ui->TextOut_H_Amplitude->toPlainText();
    IPN->h = str.toDouble();
    str = ui->TextOut_H_Period->toPlainText();
    IPN->MagneticFieldPeriod = str.toInt();

    IPN->TemperatureFieldMode = ui->comboBox_TemperatureField->currentIndex();
    str = ui->TextOut_T_Amplitude->toPlainText();
    IPN->t = str.toDouble();
    str = ui->TextOut_T_Period->toPlainText();
    IPN->TemperatureFieldPeriod = str.toInt();

    IPN->ElasticFieldMode = ui->comboBox_ElasticField->currentIndex();
    Elastic_Field_Mode EFM;
    EFM = (Elastic_Field_Mode)IPN->ElasticFieldMode;
    switch (EFM)
    {   
        case EFM_UniaxialStrain:
        {
            str = ui->TextOut_Strain_e11->toPlainText();
            IPN->CST_Uniaxial[0] = str.toDouble();
            str = ui->TextOut_Strain_e22->toPlainText();
            IPN->CST_Uniaxial[1] = str.toDouble();
            str = ui->TextOut_Strain_e33->toPlainText();
            IPN->CST_Uniaxial[2] = str.toDouble();
            str = ui->TextOut_Strain_e12->toPlainText();
            IPN->CST_Uniaxial[3] = str.toDouble();
            str = ui->TextOut_Strain_e13->toPlainText();
            IPN->CST_Uniaxial[4] = str.toDouble();
            str = ui->TextOut_Strain_e23->toPlainText();
            IPN->CST_Uniaxial[5] = str.toDouble();
            break;
        }

        case EFM_UniaxialStress:
        {
            str = ui->TextOut_Stress_e11->toPlainText();
            IPN->CST_Uniaxial[0] = str.toDouble();
            str = ui->TextOut_Stress_e22->toPlainText();
            IPN->CST_Uniaxial[1] = str.toDouble();
            str = ui->TextOut_Stress_e33->toPlainText();
            IPN->CST_Uniaxial[2] = str.toDouble();
            str = ui->TextOut_Stress_e12->toPlainText();
            IPN->CST_Uniaxial[3] = str.toDouble();
            str = ui->TextOut_Stress_e13->toPlainText();
            IPN->CST_Uniaxial[4] = str.toDouble();
            str = ui->TextOut_Stress_e23->toPlainText();
            IPN->CST_Uniaxial[5] = str.toDouble();
            break;
        }

        case EFM_EdgeDislocation:
        {
            str = ui->TextOut_Edge_Dislocation_b->toPlainText();
            IPN->b = str.toDouble();
            break;
        }

        case EFM_ScrewDislocation:
        {
            str = ui->TextOut_Screw_Dislocation_b->toPlainText();
            IPN->b = str.toDouble();
            break;
        }

        case EFM_CrackI:
        {
            str = ui->TextOut_Crack_a->toPlainText();
            IPN->a = str.toDouble();
            break;
        }

        case EFM_ConcentratedForce:
        {
            str = ui->TextOut_Concentrated_Force_FC->toPlainText();
            IPN->FC = str.toDouble();
            break;
        }

        case EFM_BeamUniformForce:
        {
            str = ui->TextOut_Beam_Uniform_Force_FU->toPlainText();
            IPN->FU_Beam = str.toDouble();
            break;
        }

        case EFM_BeamConcentratedForce:
        {
            str = ui->TextOut_Beam_Concentrated_Force_FC->toPlainText();
            IPN->FC_Beam = str.toDouble();
            break;
        }

        case EFM_BeamBending:
        {
            str = ui->TextOut_Beam_Bending_M->toPlainText();
            IPN->B_M = str.toDouble();
            break;
        }

        case EFM_ReadAnsys:
        {
            str = ui->TextOut_Read_Ansys_SF->toPlainText();
            IPN->ElasticFieldScaleFactor = str.toDouble();
            break;
        }

        case EFM_HollowTorsion:
        {
            str = ui->TextOut_Annulus_Torsion_T->toPlainText();
            IPN->B_M = str.toDouble();
            break;
        }

        case EFM_EllipseTorsion:
        {
            str = ui->TextOut_Ellipse_Torsion_T->toPlainText();
            IPN->B_M = str.toDouble();
            break;
        }

        default :
            break;
    }

    str = ui->TextOut_N_Image->toPlainText();
    IPN->N_Image = str.toInt();
///////////////////////////////////////////////////////////////////////////////////////////////////////
    Boundary_Shape *BS;
    Calculate_Points_Parameters *CPP;
    BS = &PSRP.BS;
    CPP = &PSRP.CPP;
    int m, n, l;
    //GUI Parameter
    for (i = 0; i < N_Thread_Modify; ++i)
    {
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[i]);
        GPM.Thread_State[i][0] = Signal_GUI_Run;
        GPM.Thread_State[i][1] = Signal_GUI_Continue;
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[i]);
    }

    m=2*IPN->xn+3;   n=2*IPN->yn+3;
    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l=2*IPN->zn+3;
    }
    Mx_Global=Make3DArray(m,n,l);
	My_Global=Make3DArray(m,n,l);
	Mz_Global=Make3DArray(m,n,l);

    Mx_OPGL=Make3DArray(m,n,l);
	My_OPGL=Make3DArray(m,n,l);
	Mz_OPGL=Make3DArray(m,n,l);

    BS->xn=IPN->xn;  BS->yn=IPN->yn;
    BS->zn=IPN->zn;
    BS->Boundary=Make3DArrayinteger(m,n,l);
    BS->LocalEnergy=Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    printf("%d\n",BS->Calculate_Points);
    
    CPP->Position_3DTo1D=Make3DArrayinteger(m,n,l);
    CPP->Local_Points_1D=Make3DArrayinteger(BS->Calculate_Points,7,7);
    CPP->Virtual_Position=Make1DArrayinteger(BS->Virtual_Points);
    CPP->Local_Points_1D_LocalEnergy=Make4DArrayinteger(BS->Calculate_Points,7,7,7);
    CalculPointsParameters_Initial(CPP,BS);

    MxV_Global=Make1DArray(BS->Virtual_Points);
    MyV_Global=Make1DArray(BS->Virtual_Points);
    MzV_Global=Make1DArray(BS->Virtual_Points);

    MxV_OPGL=Make1DArray(BS->Virtual_Points);
    MyV_OPGL=Make1DArray(BS->Virtual_Points);
    MzV_OPGL=Make1DArray(BS->Virtual_Points);
    for (i = 0; i < BS->Virtual_Points; ++i)
    {
        MxV_Global[i]=0.;
        MyV_Global[i]=0.;
        MzV_Global[i]=1;
    }
      
    InitialM_3D(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, BS);
    Trim_Boundary_M3D_BoundaryToVirtual(BS,CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);

    
    if (ui->radioButton_VersionPlus->isChecked() == true)
    {
        IPN->VersionMode = VM_Plus;
    }else if (ui->radioButton_VersionMinus->isChecked() == true)
    {
        IPN->VersionMode = VM_Minus;
    }

    switch (IPN->VersionMode)
    {
        case VM_Plus:
        {
            IPN->Version = IPN->Version + 1;
            break;
        }

        case VM_Minus:
        {
            IPN->Version = IPN->Version - 1;
            break;
        }
        
    
        default:
        break;
    }
    ui->TextOut_Version->setText(QString::number(IPN->Version,10));
    Make_NMD_Result_Dir(IPN);
    
}

void GUI::Pthread_State_Run_Part_Run()
{
    Calculate_Mode_GUI CMG;
    CMG = (Calculate_Mode_GUI)IPN->CalculateMode;
    switch (CMG)
    {
        case CMG_MC:
        {
            DiscreatPoint_Run_MC(IPN);
            break;
        }

        case CMG_NL:
        {
            DiscreatPoint_Run_NL(IPN);
            break;
        }

        case CMG_LLG_STT:
        {
            DiscreatPoint_Run_LLG(IPN);
            break;
        }

        case CMG_LLG_SOT:
        {
            DiscreatPoint_Run_LLG(IPN);
            break;
        }

        case CMG_GNEB:
        {
            DiscreatPoint_Run_GNEB(IPN);
            break;
        }
        
    
        default:
        break;
    }
    pthread_t tid_Part_Free;
    int err = 0;
    err = pthread_create(&tid_Part_Free, NULL, Pthread_State_Run_Part_WriteResult_Transfer, NULL);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
}

void Pthread_State_Run_Part_Run_Strain_Thread(Pthread_State_Run_Strain_Parameter *PSRSP)
{
    int State_Pthread_Strain;
    printf("%d\n", PSRSP->Index_Pthread);
    if (PSRSP->Index_Pthread != 0)
    {
        do
        {
            Sleep(100);
            pthread_mutex_lock(&PSRP.State_Pthread_Strain_Mutex);
            State_Pthread_Strain = PSRP.State_Pthread_Strain[PSRSP->Index_Pthread - 1];
            pthread_mutex_unlock(&PSRP.State_Pthread_Strain_Mutex);
        } while (State_Pthread_Strain == SPS_Prepare);

        pthread_join(PSRP.tid_Part_Free[PSRSP->Index_Pthread - 1], NULL);
    }

    InitialM_3D_SKN(IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global, &PSRP.BS, PSRSP->dir_SKN);
    Trim_Boundary_M3D_BoundaryToVirtual(&PSRP.BS,&PSRP.CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);
    
    Calculate_Mode_GUI CMG;
    CMG = (Calculate_Mode_GUI)IPN->CalculateMode;
    switch (CMG)
    {
        case CMG_MC:
        {
            DiscreatPoint_Run_MC(IPN);
            break;
        }

        case CMG_NL:
        {
            DiscreatPoint_Run_NL(IPN);
            break;
        }

        case CMG_LLG_STT:
        {
            DiscreatPoint_Run_LLG(IPN);
            break;
        }

        case CMG_LLG_SOT:
        {
            DiscreatPoint_Run_LLG(IPN);
            break;
        }

        case CMG_GNEB:
        {
            DiscreatPoint_Run_GNEB(IPN);
            break;
        }

        default:
            break;
    }

    int err = 0;
    err = pthread_create(&PSRP.tid_Part_Free[PSRSP->Index_Pthread], NULL, Pthread_State_Run_Part_WriteResult_Strain_Transfer, PSRSP);
    if(err != 0)
    {
        printf("create thread error\n");
        Sleep(1000);
    }
}

void *Pthread_State_Run_Part_Run_Strain_Thread_Transfer(void *arg)
{
    Pthread_State_Run_Strain_Parameter *PSRSP;
    PSRSP = (Pthread_State_Run_Strain_Parameter *)arg;
    Pthread_State_Run_Part_Run_Strain_Thread(PSRSP);
    return NULL;
}

void GUI::Pthread_State_Run_Part_Run_Strain()
{
    FILE *fp_Note;
    char *file_SKN;
    file_SKN = Make1DString(500);
    int Version, i, j, err = 0, N;

    Version = PSRP.Version;

    N = 11;
    pthread_t *tid_Part_Run, *tid_Part_Free;
    tid_Part_Run  = (pthread_t *)malloc(sizeof(pthread_t) * N);
    tid_Part_Free = (pthread_t *)malloc(sizeof(pthread_t) * N);
    PSRP.tid_Part_Run  = tid_Part_Run;
    PSRP.tid_Part_Free = tid_Part_Free;
    PSRP.file_SKN     = file_SKN;

    Pthread_State_Run_Strain_Parameter *PSRSP;
    PSRSP = (Pthread_State_Run_Strain_Parameter *)malloc(sizeof(Pthread_State_Run_Strain_Parameter) * N);
    PSRP.State_Pthread_Strain = Make1DArrayinteger(N);
    pthread_mutex_init(&PSRP.State_Pthread_Strain_Mutex, 0);
    for (j = 1; j < 2; ++j)
    {
        //h = 2.7 - j * 0.05;
        //sprintf(file_SKN, "h=%0.2f\\Note.dat", h);
        sprintf(file_SKN, "Note.dat");
        fp_Note = fopen(file_SKN, "w+");
        PSRP.fp_Note = fp_Note;
        fprintf(fp_Note,"Number\t Energy\t SKN\n");
        for (i = 1; i < N + 1; ++i)
        {
            PSRSP[i - 1].dir_SKN = Make1DString(500);
            //IPN->h = h;
            if(i==0)
            {
                //sprintf(dir_SKN, "h=%0.2f\\FE", h);
            }else
            {
                //sprintf(dir_SKN, "h=%0.2f\\%d", h, i);
                sprintf(PSRSP[i - 1].dir_SKN, "%d", i);
            }
            PSRSP[i - 1].Index_Pthread = i - 1;
            PSRSP[i - 1].N_Pthread = N;
            PSRP.State_Pthread_Strain[i - 1] = SPS_Prepare;

            IPN->Version = Version;
            switch (IPN->VersionMode)
            {
                case VM_Plus:
                {
                    IPN->Version = IPN->Version + 1;
                    break;
                }

                case VM_Minus:
                {
                    IPN->Version = IPN->Version - 1;
                    break;
                }
    
                default:
                break;
            }

            Make_NMD_Result_Dir_SKN(IPN,PSRSP[i - 1].dir_SKN);

            err = pthread_create(&tid_Part_Run[i - 1], NULL, Pthread_State_Run_Part_Run_Strain_Thread_Transfer, &PSRSP[i - 1]);
            if(err != 0)
            {
                printf("create thread error\n");
                Sleep(1000);
            }

            Sleep(100);
        }    
    }
    WriteInput_NMD(IPN);

}

void Pthread_State_Run_Part_WriteResult_Strain(Pthread_State_Run_Strain_Parameter *PSRSP)
{
    printf("WriteResult Thread Run %d\n", PSRSP->Index_Pthread);

    pthread_mutex_lock(&PSRP.State_Pthread_Strain_Mutex);
    PSRP.State_Pthread_Strain[PSRSP->Index_Pthread] = SPS_Run;
    pthread_mutex_unlock(&PSRP.State_Pthread_Strain_Mutex);

    int State_GUI_Run;
    
    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    State_GUI_Run = GPM.State_GUI_Run;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    if (State_GUI_Run != TSG_Finish)
    {
        pthread_mutex_lock(&GPM.Run_GUI_Mutex);
        pthread_cond_wait(&GPM.Run_GUI_Cond, &GPM.Run_GUI_Mutex);
        pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    }

    printf("%d\t %0.8f\t %0.8f\n",PSRSP->Index_Pthread,IPN->W,IPN->SKN);
    fprintf(PSRP.fp_Note,"%d\t %0.8f\t %0.8f\n",PSRSP->Index_Pthread,IPN->W,IPN->SKN);
    Sleep(100);
    Trim_Boundary_M3D_VirtualToBoundary(&PSRP.BS,&PSRP.CPP,IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global);
    Write_NMD_Magnetization_SKN(IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global,&PSRP.BS,PSRSP->dir_SKN);
    Write_Inplane_M_SKN(IPN, &PSRP.BS, Mx_Global, My_Global, PSRSP->dir_SKN);
    free(PSRSP->dir_SKN);

    if (PSRSP->Index_Pthread == PSRSP->N_Pthread - 1)
    {
        int m, n, l;
        m=2*IPN->xn+3;   n=2*IPN->yn+3;
        if(IPN->zn==0)
        {
            l=1;
        }else
        {
            l=2*IPN->zn+3;
        }

        WriteInput_NMD(IPN);
    
    
        free3DArray(Mx_Global,m,n);
        free3DArray(My_Global,m,n);
        free3DArray(Mz_Global,m,n);
        free3DArray(Mx_OPGL, m, n);
        free3DArray(My_OPGL, m, n);
        free3DArray(Mz_OPGL, m, n);
        free(MxV_Global);
        free(MyV_Global);
        free(MzV_Global);
        free(MxV_OPGL);
        free(MyV_OPGL);
        free(MzV_OPGL);
        free3DArrayinteger(PSRP.BS.Boundary,m,n);
        free4DArrayinteger(PSRP.BS.LocalEnergy,m,n,l);
        free(PSRP.CPP.Virtual_Position);
        free4DArrayinteger(PSRP.CPP.Local_Points_1D_LocalEnergy,PSRP.BS.Calculate_Points,7,7);
        free3DArrayinteger(PSRP.CPP.Local_Points_1D,PSRP.BS.Calculate_Points,7);
        free3DArrayinteger(PSRP.CPP.Position_3DTo1D,m,n);
        free(PSRP.file_SKN);
        free(PSRP.tid_Part_Run);
        free(PSRP.tid_Part_Free);
        free(PSRP.State_Pthread_Strain);
        fclose(PSRP.fp_Note);
    }

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Prepare;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);

    pthread_mutex_lock(&PSRP.State_Pthread_Strain_Mutex);
    PSRP.State_Pthread_Strain[PSRSP->Index_Pthread] = SPS_Finish;
    pthread_mutex_unlock(&PSRP.State_Pthread_Strain_Mutex);
    printf("Write result thread %d Done\n", PSRSP->Index_Pthread);
}

void *Pthread_State_Run_Part_WriteResult_Strain_Transfer(void *arg)
{
    Pthread_State_Run_Strain_Parameter *PSRSP;
    PSRSP = (Pthread_State_Run_Strain_Parameter *)arg;
    Pthread_State_Run_Part_WriteResult_Strain(PSRSP);
    return NULL;
}

void Pthread_State_Run_Part_WriteResult()
{
    int State_GUI_Run;
    
    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    State_GUI_Run = GPM.State_GUI_Run;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    if (State_GUI_Run != TSG_Finish)
    {
        pthread_mutex_lock(&GPM.Run_GUI_Mutex);
        pthread_cond_wait(&GPM.Run_GUI_Cond, &GPM.Run_GUI_Mutex);
        pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    }

    Trim_Boundary_M3D_VirtualToBoundary(&PSRP.BS, &PSRP.CPP, IPN, Mx_Global, My_Global, Mz_Global, MxV_Global, MyV_Global, MzV_Global);
    Write_NMD_Magnetization(IPN,Mx_Global,My_Global,Mz_Global,MxV_Global,MyV_Global,MzV_Global,&PSRP.BS);
    double ****CST;
    int m, n, l;
    m=2*IPN->xn+3;   n=2*IPN->yn+3;
    if(IPN->zn==0)
    {
        l=1;
    }else
    {
        l=2*IPN->zn+3;
    }
    CST=Make4DArray(7,m,n,l);
    Elastic_Field_Initial(CST,IPN);
    
    //ReadAnsys(IPN,CST);

    Write_NMD_ElastiField(IPN,CST,&PSRP.BS);
    Write_Principal_Stress(IPN,CST,&PSRP.BS);
    Write_Principal_Strain(IPN,CST,&PSRP.BS);


    Write_NMD_Result(IPN);
    Write_Inplane_M(IPN,&PSRP.BS,Mx_Global,My_Global);
    Write_Judge_Phase_Curve(IPN,&PSRP.BS,CST);
    free4DArray(CST,7,m,n);
    WriteInput_NMD(IPN);
    
    
    free3DArray(Mx_Global,m,n);
    free3DArray(My_Global,m,n);
    free3DArray(Mz_Global,m,n);
    free3DArray(Mx_OPGL, m, n);
    free3DArray(My_OPGL, m, n);
    free3DArray(Mz_OPGL, m, n);
    free(MxV_Global);
    free(MyV_Global);
    free(MzV_Global);
    free(MxV_OPGL);
    free(MyV_OPGL);
    free(MzV_OPGL);
    free3DArrayinteger(PSRP.BS.Boundary,m,n);
    free4DArrayinteger(PSRP.BS.LocalEnergy,m,n,l);
    free(PSRP.CPP.Virtual_Position);
    free4DArrayinteger(PSRP.CPP.Local_Points_1D_LocalEnergy,PSRP.BS.Calculate_Points,7,7);
    free3DArrayinteger(PSRP.CPP.Local_Points_1D,PSRP.BS.Calculate_Points,7);
    free3DArrayinteger(PSRP.CPP.Position_3DTo1D,m,n);

    pthread_mutex_lock(&GPM.Run_GUI_Mutex);
    GPM.State_GUI_Run = TSG_Prepare;
    pthread_mutex_unlock(&GPM.Run_GUI_Mutex);
    printf("Done\n");
}

void *Pthread_State_Run_Part_WriteResult_Transfer(void *arg)
{
    Pthread_State_Run_Part_WriteResult();
    return NULL;
}

void GUI::Pthread_State_Run()
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
            Pthread_State_Run_Part_Prepare();
            Pthread_State_Run_Part_Run();
            //Pthread_State_Run_Part_Run_Strain();
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
}

void GUI::Pthread_State_Stop()
{
    int i, N_Thread;
    //N_Thread = Run thread + Synchronize thread + OPGL thread
    N_Thread = IPN->Threads_Number + 2;
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[i]);
        GPM.Thread_State[i][0] = Signal_GUI_Stop;
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[i]);
    }
}

void GUI::Pthread_State_Continue()
{
    int i, N_Thread, Signal_ContinueAndBreak;
    //N_Thread = Run thread + Synchronize thread + OPGL thread
    N_Thread = IPN->Threads_Number + 2;
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[i]);
        Signal_ContinueAndBreak = GPM.Thread_State[i][1];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[i]);
        if(Signal_ContinueAndBreak != Signal_GUI_Break)
        {
            return;
        }
    }

    for (i = 0; i < N_Thread; ++i)
    {
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[i]);
        GPM.Thread_State[i][1] = Signal_GUI_Continue;
        pthread_cond_signal(&GPM.Pthread_GUI_Cond[i]);
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[i]);
    }
}

void GUI::Pthread_State_Break()
{
    int i, N_Thread;
    //N_Thread = Run thread + Synchronize thread + OPGL thread
    N_Thread = IPN->Threads_Number + 2;
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[i]);
        GPM.Thread_State[i][1] = Signal_GUI_Break;
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[i]);
    }
}


GUI::GUI(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::GUI)
{
    ui->setupUi(this);
    GPM.State_GUI_Run = TSG_Prepare;
    int i, N_Thread;
    //N_Thread = Run thread + Synchronize thread
    N_Thread = IPN->Threads_Number + 2;
    //////////////////////////  Signal Parameter /////////////////////////////////////
    GPM.Pthread_GUI_Cond  = (pthread_cond_t *)malloc(sizeof(pthread_cond_t) * N_Thread);
    GPM.Pthread_GUI_Mutex = (pthread_mutex_t *)malloc(sizeof(pthread_mutex_t) * N_Thread);
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_cond_init(&GPM.Pthread_GUI_Cond[i], NULL);
        pthread_mutex_init(&GPM.Pthread_GUI_Mutex[i], 0);
    }
    pthread_mutex_init(&GPM.Run_GUI_Mutex, 0);
    pthread_cond_init(&GPM.Run_GUI_Cond, 0);

    GPM.Thread_State = Make2DArrayinteger(N_Thread, 2);
    Sleep(100);
    GUI_Initial(IPN);
    connect(ui->comboBox_Calculate_Mode,SIGNAL(currentIndexChanged(int)),this,SLOT(SetPage_CalculateMode(int)));
    connect(ui->comboBox_ElasticField,SIGNAL(currentIndexChanged(int)),this,SLOT(SetPage_ElasticField(int)));
    connect(ui->pushButton_Run, SIGNAL(clicked()), this, SLOT(Pthread_State_Run()));
    connect(ui->pushButton_Stop, SIGNAL(clicked()), this, SLOT(Pthread_State_Stop()));
    connect(ui->pushButton_Continue, SIGNAL(clicked()), this, SLOT(Pthread_State_Continue()));
    connect(ui->pushButton_Break, SIGNAL(clicked()), this, SLOT(Pthread_State_Break()));
}

GUI::~GUI()
{
    delete ui;
    int i, N_Thread;
    N_Thread = IPN->Threads_Number + 2;
    //////////////////////////  Signal Parameter /////////////////////////////////////
    for (i = 0; i < N_Thread; ++i)
    {
        pthread_cond_destroy(&GPM.Pthread_GUI_Cond[i]);
        pthread_mutex_destroy(&GPM.Pthread_GUI_Mutex[i]);
    }
    pthread_mutex_destroy(&GPM.Run_GUI_Mutex);
    pthread_cond_destroy(&GPM.Run_GUI_Cond);
    free(GPM.Pthread_GUI_Cond);
    free(GPM.Pthread_GUI_Mutex);
}

