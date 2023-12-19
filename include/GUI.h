#ifndef GUI_H
#define GUI_H

#include <QMainWindow>
#include <Variable.h>

QT_BEGIN_NAMESPACE
namespace Ui { class GUI; }
QT_END_NAMESPACE

class GUI : public QMainWindow
{
    Q_OBJECT

public:
    GUI(QWidget *parent = nullptr);
    ~GUI();

private slots:
    void SetPage_CalculateMode(int Index_Page_CalculateMode);
    void SetPage_ElasticField(int Index_Page_ElasticField);
    void Pthread_State_Run();
    void Pthread_State_Run_Part_Prepare();
    void Pthread_State_Run_Part_Run();
    void Pthread_State_Run_Part_Run_Strain();
    void Pthread_State_Stop();
    void Pthread_State_Continue();
    void Pthread_State_Break();

private:
    Ui::GUI *ui;
    void GUI_Initial(Input_Parameter_NMD *IPN);
};

void Pthread_State_Run_Part_Run_Strain_Thread(Pthread_State_Run_Strain_Parameter *PSRSP);
void *Pthread_State_Run_Part_Run_Strain_Thread_Transfer(void *arg);
void Pthread_State_Run_Part_WriteResult();
void *Pthread_State_Run_Part_WriteResult_Transfer(void *arg);
void Pthread_State_Run_Part_WriteResult_Strain(Pthread_State_Run_Strain_Parameter *PSRSP);
void *Pthread_State_Run_Part_WriteResult_Strain_Transfer(void *arg);
#endif // GUI_H
