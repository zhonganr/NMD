#ifndef OPGL_H_INCLUDED
#define OPGL_H_INCLUDED

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <Variable.h>
void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);
int Pastel_RGB(int **RGB_Matrix, int Level);
int Background_Matrix_Assignment(double ***Mz, int **Background_Matrix, int m, int n);
int Scale_Matrix_Assignment(double ***Mz, int **Scale_Matrix, int m, int n, int Level);
int Update_Background_Rectangle_Vertex(double ****Background_Matrix_Vertex_Coordinate, float *Background_Rectangle_Vertex, int i, int j);
int Update_Scale_Rectangle_Vertex(double ****Scale_Matrix_Vertex_Coordinate, float *Scale_Rectangle_Vertex, int i, int j);
int Background_Matrix_Vertex_Coordinate_Assignment(double ****Background_Matrix_Vertex_Coordinate, int m, int n, int Grid_m_Max, int Grid_m_Min, int Grid_n_Max, int Grid_n_Min);
int Scale_Matrix_Vertex_Coordinate_Assignment(double ****Scale_Matrix_Vertex_Coordinate, int Level);
int OPGL_InPlaneComponent_Parameter(int m, int n, int Inplane_Density, int *OI_Parameter);
int OPGL_InPlaneComponent_Coordiante(double ***Mx, double ***My, double ***Mz, int D,
									int *OI_Parameter, int Inplane_Density, double ****Arrow_Head_Coordinate, double ****Arrow_Body_Coordinate,
									int Grid_m_Max, int Grid_m_Min, int Grid_n_Max, int Grid_n_Min);
int Update_Arrow_Vertex(double ****Arrow_Head_Coordinate, double ****Arrow_Body_Coordinate, float *Arrow_Head_Vertex, float *Arrow_Body_Vertex, int i, int j);
void SetGridSize(int m, int n, double *GridSize);
int OPGL_NMD(OPGL_Parameters *OPM);
int OPGL_GNEB(Input_Parameter_NMD *IPN);
int Get_Max_Min(double ***Mz, float *MzMaxMin, int m, int n, int Scale);

#endif // OPGL_H_INCLUDED