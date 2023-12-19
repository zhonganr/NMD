#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <Array.h>
#include <malloc.h>
#include <math.h>
#include <iostream>
#include <map>
#include <string>
#include <stdio.h>
#include <pthread.h>
#include <OPGL.h>
#include <BoundaryShape.h>
#include <Variable.h>
#include <Energy.h>
#include <windows.h>
#include <ft2build.h>
#include FT_FREETYPE_H 
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <learnopengl/shader.h>

#define PI      3.14159265359

pthread_mutex_t OPGL_Mutex;
double Grid_m_Max = 0.5, Grid_m_Min = -0.9, Grid_n_Max = 0.7, Grid_n_Min = -0.7;//0.112
double Scale_m_Max = 0.7, Scale_m_Min = 0.6, Scale_n_Max = 0.7, Scale_n_Min = -0.7;
OPGL_Parameters_GUI OPGL_GUI;
extern int N_Energy_Terms;
extern GUI_Parameter GPM;

void SetGridSize(int m, int n, double *GridSize)//Grid_m_Max, int Grid_m_Min, int Grid_n_Max, int Grid_n_Min
{

	//int MaxSize;
	if (m >= n)
	{
		GridSize[0] = 0.5;
		GridSize[1] = -0.9;
		GridSize[2] = 1.4 * n / m / 2;
		GridSize[3] = -1.4 * n / m / 2;
	}
	else
	{
		GridSize[0] = 0.7;
		GridSize[1] = -0.7;
		GridSize[2] = -0.2+1.4 * m / n / 2;
		GridSize[3] = -0.2-1.4 * m / n / 2;
	}
}

// double Mz_Max=0., Mz_Min=0.;
// Mz_Max = 0.;
// Mz_Min = 0.;
 
void RenderText(Shader &shader, std::string text, float x, float y, float scale, glm::vec3 color, unsigned int VAOId_Text, unsigned int VBOId_Text); 
// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 1000;


/// Holds all state information relevant to a character as loaded using FreeType

struct Character {
    unsigned int TextureID; // ID handle of the glyph texture
    glm::ivec2   Size;      // Size of glyph
    glm::ivec2   Bearing;   // Offset from baseline to left/top of glyph
    unsigned int Advance;   // Horizontal offset to advance to next glyph
};

std::map<GLchar, Character> Characters;

const char *vertexShaderSource = "#version 330 core\n"
    "layout (location = 0) in vec3 aPos;\n"
    "void main()\n"
    "{\n"
    "   gl_Position = vec4(aPos.x, aPos.y, aPos.z, 1.0);\n"
    "}\0";
const char *fragmentShader1Source = "#version 330 core\n"
    "out vec4 FragColor;\n"
    "void main()\n"
    "{\n"
    "   FragColor = vec4(1.0f, 0.5f, 0.2f, 1.0f);\n"
    "}\n\0";
const char *fragmentShaderSource2 = "#version 330 core\n"
    "out vec4 FragColor;\n"
    "void main()\n"
    "{\n"
    "   FragColor = vec4(1.0f, 1.0f, 0.0f, 1.0f);\n"
    "}\n\0";

const char *fragmentShaderSource = "#version 330 core\n"
 	 "out vec4 FragColor;\n"
 	 "uniform vec4 ourColor;\n"
 	 "void main()\n"
 	 "{\n"
 	 "   FragColor = ourColor;\n"
 	 "}\n\0";



// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------
void processInput(GLFWwindow *window)
{
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);
}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and 
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}


int Background_Matrix_Vertex_Coordinate_Assignment(double ****Background_Matrix_Vertex_Coordinate, int m, int n, double Grid_m_Max, double Grid_m_Min, double Grid_n_Max, double Grid_n_Min)
{
	int i, j;
	double Grid_Interval_m, Grid_Interval_n;
	Grid_Interval_m = (Grid_m_Max - Grid_m_Min) / m;	Grid_Interval_n = (Grid_n_Max - Grid_n_Min) / n;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			//bottom left
			Background_Matrix_Vertex_Coordinate[i][j][0][0] = Grid_m_Min + i * Grid_Interval_m;
			Background_Matrix_Vertex_Coordinate[i][j][0][1] = Grid_n_Min + j * Grid_Interval_n;
			Background_Matrix_Vertex_Coordinate[i][j][0][2] = 0.;

			//bottom right
			Background_Matrix_Vertex_Coordinate[i][j][1][0] = Grid_m_Min + (i + 1) * Grid_Interval_m;
			Background_Matrix_Vertex_Coordinate[i][j][1][1] = Grid_n_Min + j * Grid_Interval_n;
			Background_Matrix_Vertex_Coordinate[i][j][1][2] = 0.;

			//top right
			Background_Matrix_Vertex_Coordinate[i][j][2][0] = Grid_m_Min + (i + 1) * Grid_Interval_m;
			Background_Matrix_Vertex_Coordinate[i][j][2][1] = Grid_n_Min + (j + 1) * Grid_Interval_n;
			Background_Matrix_Vertex_Coordinate[i][j][2][2] = 0.;

			//top left
			Background_Matrix_Vertex_Coordinate[i][j][3][0] = Grid_m_Min + i * Grid_Interval_m;
			Background_Matrix_Vertex_Coordinate[i][j][3][1] = Grid_n_Min + (j + 1) * Grid_Interval_n;
			Background_Matrix_Vertex_Coordinate[i][j][3][2] = 0.;
		}
	}

	return 1;
}

int Scale_Matrix_Vertex_Coordinate_Assignment(double ****Scale_Matrix_Vertex_Coordinate, int Level)
{
	int i, j;
	double Scale_Interval_m, Scale_Interval_n;
	Scale_Interval_m = (Scale_m_Max - Scale_m_Min);	Scale_Interval_n = (Scale_n_Max - Scale_n_Min) / Level;
	for (i = 0; i < 1; ++i)
	{
		for (j = 0; j < Level; ++j)
		{
			//bottom left
            Scale_Matrix_Vertex_Coordinate[i][j][0][0] = Scale_m_Min + i * Scale_Interval_m;
			Scale_Matrix_Vertex_Coordinate[i][j][0][1] = Scale_n_Min + j * Scale_Interval_n;
			Scale_Matrix_Vertex_Coordinate[i][j][0][2] = 0.;

			//bottom right
			Scale_Matrix_Vertex_Coordinate[i][j][1][0] = Scale_m_Min + (i + 1) * Scale_Interval_m;
			Scale_Matrix_Vertex_Coordinate[i][j][1][1] = Scale_n_Min + j * Scale_Interval_n;
			Scale_Matrix_Vertex_Coordinate[i][j][1][2] = 0.;

			//top right
			Scale_Matrix_Vertex_Coordinate[i][j][2][0] = Scale_m_Min + (i + 1) * Scale_Interval_m;
			Scale_Matrix_Vertex_Coordinate[i][j][2][1] = Scale_n_Min + (j + 1) * Scale_Interval_n;
			Scale_Matrix_Vertex_Coordinate[i][j][2][2] = 0.;

			//top left
			Scale_Matrix_Vertex_Coordinate[i][j][3][0] = Scale_m_Min + i * Scale_Interval_m;
			Scale_Matrix_Vertex_Coordinate[i][j][3][1] = Scale_n_Min + (j + 1) * Scale_Interval_n;
			Scale_Matrix_Vertex_Coordinate[i][j][3][2] = 0.;
		}
	}

	return 1;
}

int Update_Background_Rectangle_Vertex(double ****Background_Matrix_Vertex_Coordinate, float *Background_Rectangle_Vertex, int i, int j)
{
	int k, r;
	for (k = 0; k < 4; ++k)
	{
		for (r = 0; r < 3; ++r)
		{
			Background_Rectangle_Vertex[3*k+r] = Background_Matrix_Vertex_Coordinate[i][j][k][r];
		}	
	}

	return 1;
}

int Update_Scale_Rectangle_Vertex(double ****Scale_Matrix_Vertex_Coordinate, float *Scale_Rectangle_Vertex, int i, int j)
{
	int k, r;
	for (k = 0; k < 4; ++k)
	{
		for (r = 0; r < 3; ++r)
		{
			Scale_Rectangle_Vertex[3*k+r] = Scale_Matrix_Vertex_Coordinate[i][j][k][r];
		}	
	}

	return 1;
}

int Background_Matrix_Assignment(double ***Mz, int **Background_Matrix, int m, int n)
{
	int i, j, Level, Level_Actual;
	double Max, Min, Level_Interval;
	Max = Mz[1][1][0];	Min = Mz[1][1][0];
	Level = 45;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (Mz[i][j][0] >= Max)
			{
				Max = Mz[i][j][0];
			}

			if (Mz[i][j][0] <= Min)
			{
				Min = Mz[i][j][0];
			}
		}
	}

	Level_Interval = (Max - Min) / Level;

	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			Level_Actual = round((Mz[i][j][0] - Min) / Level_Interval);
			if(Level_Actual>=Level)
			{
				Level_Actual = Level-1;
			}
			Background_Matrix[i][j] = Level_Actual;
		}
	}

	return 1;
}

int Scale_Matrix_Assignment(double ***Mz, int **Scale_Matrix, int m, int n, int Level)
{
	int i, j, Level_Actual;
	double Max, Min;
	//double Level_Interval;
	Max = Mz[1][1][0];	Min = Mz[1][1][0];
	//Level = 5;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (Mz[i][j][0] >= Max)
			{
				Max = Mz[i][j][0];
			}

			if (Mz[i][j][0] <= Min)
			{
				Min = Mz[i][j][0];
			}
		}
	}

	//Level_Interval = (Max - Min) / Level;

	for (i = 0; i < 1; ++i)
	{
		for (j = 0; j < Level; ++j)
		{
			Level_Actual = j;
			Scale_Matrix[i][j] = Level_Actual;
		}
	}

	return 1;
}

int Get_Max_Min(double ***Mz, float *MzMaxMin, int m, int n, int Scale)
{
	int i, j;
	float Interval;
	MzMaxMin[0] = Mz[1][1][0];	MzMaxMin[Scale] = Mz[1][1][0];
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if (Mz[i][j][0] >= MzMaxMin[0])
			{
				MzMaxMin[0] = Mz[i][j][0];
			}

			if (Mz[i][j][0] <= MzMaxMin[Scale])
			{
				MzMaxMin[Scale] = Mz[i][j][0];
			}
		}
	}
	Interval = (MzMaxMin[0] - MzMaxMin[Scale]) / Scale;
	for (i = 0; i < Scale+1; ++i)
	{
		MzMaxMin[i] = MzMaxMin[0] - i * Interval;
	}
	return 1;
}

int Pastel_RGB(int **RGB_Matrix_Get, int Level)
{
	int i,j;
	int **RGB_Matrix;
	RGB_Matrix = Make2DArrayinteger(45, 3);
	RGB_Matrix[0][0]  = 203;	RGB_Matrix[0][1]  = 189;	RGB_Matrix[0][2]  = 223;
	RGB_Matrix[1][0]  = 207;	RGB_Matrix[1][1]  = 201;	RGB_Matrix[1][2]  = 227;
	RGB_Matrix[2][0]  = 210;	RGB_Matrix[2][1]  = 203;	RGB_Matrix[2][2]  = 228;
	RGB_Matrix[3][0]  = 209;	RGB_Matrix[3][1]  = 199;	RGB_Matrix[3][2]  = 226;
	RGB_Matrix[4][0]  = 218;	RGB_Matrix[4][1]  = 203;	RGB_Matrix[4][2]  = 227;
	RGB_Matrix[5][0]  = 226;	RGB_Matrix[5][1]  = 208;	RGB_Matrix[5][2]  = 229;
	RGB_Matrix[6][0]  = 230;	RGB_Matrix[6][1]  = 213;	RGB_Matrix[6][2]  = 232;
	RGB_Matrix[7][0]  = 237;	RGB_Matrix[7][1]  = 213;	RGB_Matrix[7][2]  = 231;
	RGB_Matrix[8][0]  = 237;	RGB_Matrix[8][1]  = 207;	RGB_Matrix[8][2]  = 221;
	RGB_Matrix[9][0]  = 226;	RGB_Matrix[9][1]  = 187;	RGB_Matrix[9][2]  = 196;
	RGB_Matrix[10][0] = 235;	RGB_Matrix[10][1] = 175;	RGB_Matrix[10][2] = 181;
	RGB_Matrix[11][0] = 243;	RGB_Matrix[11][1] = 171;	RGB_Matrix[11][2] = 175;
	RGB_Matrix[12][0] = 248;	RGB_Matrix[12][1] = 176;	RGB_Matrix[12][2] = 179;
	RGB_Matrix[13][0] = 248;	RGB_Matrix[13][1] = 172;	RGB_Matrix[13][2] = 177;
	RGB_Matrix[14][0] = 246;	RGB_Matrix[14][1] = 174;	RGB_Matrix[14][2] = 179;
	RGB_Matrix[15][0] = 238;	RGB_Matrix[15][1] = 189;	RGB_Matrix[15][2] = 191;
	RGB_Matrix[16][0] = 244;	RGB_Matrix[16][1] = 187;	RGB_Matrix[16][2] = 172;
	RGB_Matrix[17][0] = 249;	RGB_Matrix[17][1] = 195;	RGB_Matrix[17][2] = 167;
	RGB_Matrix[18][0] = 249;	RGB_Matrix[18][1] = 218;	RGB_Matrix[18][2] = 192;
	RGB_Matrix[19][0] = 251;	RGB_Matrix[19][1] = 233;	RGB_Matrix[19][2] = 197;
	RGB_Matrix[20][0] = 252;	RGB_Matrix[20][1] = 241;	RGB_Matrix[20][2] = 187;
	RGB_Matrix[21][0] = 250;	RGB_Matrix[21][1] = 243;	RGB_Matrix[21][2] = 154;
	RGB_Matrix[22][0] = 253;	RGB_Matrix[22][1] = 248;	RGB_Matrix[22][2] = 171;
	RGB_Matrix[23][0] = 252;	RGB_Matrix[23][1] = 249;	RGB_Matrix[23][2] = 190;
	RGB_Matrix[24][0] = 236;	RGB_Matrix[24][1] = 240;	RGB_Matrix[24][2] = 189;
	RGB_Matrix[25][0] = 212;	RGB_Matrix[25][1] = 229;	RGB_Matrix[25][2] = 154;
	RGB_Matrix[26][0] = 190;	RGB_Matrix[26][1] = 220;	RGB_Matrix[26][2] = 130;
	RGB_Matrix[27][0] = 180;	RGB_Matrix[27][1] = 220;	RGB_Matrix[27][2] = 175;
	RGB_Matrix[28][0] = 181;	RGB_Matrix[28][1] = 221;	RGB_Matrix[28][2] = 192;
	RGB_Matrix[29][0] = 187;	RGB_Matrix[29][1] = 222;	RGB_Matrix[29][2] = 192;
	RGB_Matrix[30][0] = 185;	RGB_Matrix[30][1] = 225;	RGB_Matrix[30][2] = 210;
	RGB_Matrix[31][0] = 179;	RGB_Matrix[31][1] = 221;	RGB_Matrix[31][2] = 214;
	RGB_Matrix[32][0] = 173;	RGB_Matrix[32][1] = 215;	RGB_Matrix[32][2] = 211;
	RGB_Matrix[33][0] = 160;	RGB_Matrix[33][1] = 215;	RGB_Matrix[33][2] = 215;
	RGB_Matrix[34][0] = 151;	RGB_Matrix[34][1] = 215;	RGB_Matrix[34][2] = 222;
	RGB_Matrix[35][0] = 147;	RGB_Matrix[35][1] = 217;	RGB_Matrix[35][2] = 231;
	RGB_Matrix[36][0] = 156;	RGB_Matrix[36][1] = 234;	RGB_Matrix[36][2] = 238;
	RGB_Matrix[37][0] = 159;	RGB_Matrix[37][1] = 233;	RGB_Matrix[37][2] = 244;
	RGB_Matrix[38][0] = 160;	RGB_Matrix[38][1] = 220;	RGB_Matrix[38][2] = 249;
	RGB_Matrix[39][0] = 186;	RGB_Matrix[39][1] = 214;	RGB_Matrix[39][2] = 243;
	RGB_Matrix[40][0] = 194;	RGB_Matrix[40][1] = 210;	RGB_Matrix[40][2] = 238;
	RGB_Matrix[41][0] = 184;	RGB_Matrix[41][1] = 206;	RGB_Matrix[41][2] = 233;
	RGB_Matrix[42][0] = 141;	RGB_Matrix[42][1] = 207;	RGB_Matrix[42][2] = 238;
	RGB_Matrix[43][0] = 142;	RGB_Matrix[43][1] = 209;	RGB_Matrix[43][2] = 238;
	RGB_Matrix[44][0] = 201;	RGB_Matrix[44][1] = 214;	RGB_Matrix[44][2] = 231;
	if (Level==45)
	{
		for (i = 0; i < 45; ++i)
		{
			for (j = 0; j < 3; j++)
			{
				RGB_Matrix_Get[i][j] = RGB_Matrix[i][j];
			}
		}
	}
	else
	{
		for (i = 0; i < Level; ++i)
		{
			j = round(44 / Level * (i + 0.5));
			RGB_Matrix_Get[i][0]=RGB_Matrix[j][0];
			RGB_Matrix_Get[i][1]=RGB_Matrix[j][1];
			RGB_Matrix_Get[i][2]=RGB_Matrix[j][2];
		}
	}
	
	return 1;
}

int OPGL_InPlaneComponent_Parameter(int m, int n, int Inplane_Density, int *OI_Parameter)
{
	int i, j, D, xc, yc, N_i, N_j, N_Total;

	xc = (m - 1) / 2;
	yc = (n - 1) / 2;


	if(Inplane_Density == 0)
    {
        Inplane_Density = 24;
    }
    D = m/Inplane_Density;

    if(D<=1)
    {
        D=1;
    }

	N_Total = 0;
	N_i = 0;
	N_j = 0;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			if ((i - xc) % D == 0 && (j - yc) % D == 0)
			{
				N_Total = N_Total + 1;
			}
			
		}
	}

	for (i = 0; i < m; ++i)
	{
		if ((i - xc) % D == 0)
		{
			N_i = N_i + 1;
		}
	}

	for (j = 0; j < n; ++j)
	{
		if ((j - yc) % D == 0)
		{
			N_j = N_j + 1;
		}
	}

	OI_Parameter[0] = N_Total;
	OI_Parameter[1] = N_i;
	OI_Parameter[2] = N_j;
	return D;
}

int OPGL_InPlaneComponent_Coordiante(double ***Mx, double ***My, double ***Mz, int D,
									int *OI_Parameter, int Inplane_Density, double ****Arrow_Head_Coordinate, double ****Arrow_Body_Coordinate,
									double Grid_m_Max, double Grid_m_Min, double Grid_n_Max, double Grid_n_Min)
{
	int i, j, i_G, j_G;
	double Mx0, My0, M_Module, Phi, M_Appropriate, M_min, Ratio;

	double Grid_Interval_m, Grid_Interval_n;
	double Arrow_Body_Width_Factor, Arrow_Head_Angle, Arrow_Head_length_Factor, Arrow_Head_Bottom;

	Arrow_Body_Width_Factor = 0.2;
	Arrow_Head_Angle = PI / 6;
	Arrow_Head_length_Factor = 2;
	
	Grid_Interval_m = (Grid_m_Max - Grid_m_Min) / OI_Parameter[1];	
	Grid_Interval_n = (Grid_n_Max - Grid_n_Min) / OI_Parameter[2];

	M_Appropriate = (Grid_m_Max - Grid_m_Min)*0.8/Inplane_Density;
    M_min         = (Grid_m_Max - Grid_m_Min)*0.02/Inplane_Density;

	for (i = 0; i < OI_Parameter[1]; ++i)
	{
		for (j = 0; j < OI_Parameter[2]; ++j)
		{
			i_G = D * i;
			j_G = D * j;
			Mx0 = Mx[i_G][j_G][0];
			My0 = My[i_G][j_G][0];
			M_Module = sqrt(Mx0 * Mx0 + My0 * My0);
			Phi = atan2(My0, Mx0);
			if(M_Module>M_Appropriate)
            {
                Ratio = M_Appropriate/M_Module;
                Mx0 = Mx0*Ratio;
                My0 = My0*Ratio;
				M_Module = M_Appropriate;
			}else if(M_Module<M_min)
            {
                Mx0 = 0;
                My0 = 0;
				M_Module = 0;
			}
			//Arrow Body
			//bottom left
			Arrow_Body_Coordinate[i][j][0][0] = Grid_m_Min + i * Grid_Interval_m + cos(PI / 2 + Phi) * 0.01 * Arrow_Body_Width_Factor;
			Arrow_Body_Coordinate[i][j][0][1] = Grid_n_Min + j * Grid_Interval_n + sin(PI / 2 + Phi) * 0.01 * Arrow_Body_Width_Factor;
			Arrow_Body_Coordinate[i][j][0][2] = 0.;

			//bottom right
			Arrow_Body_Coordinate[i][j][1][0] = Grid_m_Min + i * Grid_Interval_m + cos(Phi - PI / 2) * 0.01 * Arrow_Body_Width_Factor;
			Arrow_Body_Coordinate[i][j][1][1] = Grid_n_Min + j * Grid_Interval_n + sin(Phi - PI / 2) * 0.01 * Arrow_Body_Width_Factor;
			Arrow_Body_Coordinate[i][j][1][2] = 0.;

			//top right
			if(M_Module==0)
			{
				Arrow_Body_Coordinate[i][j][2][0] = Grid_m_Min + i * Grid_Interval_m;
				Arrow_Body_Coordinate[i][j][2][1] = Grid_n_Min + j * Grid_Interval_n;
			}else
			{
				Arrow_Body_Coordinate[i][j][2][0] = Grid_m_Min + i * Grid_Interval_m + cos(Phi - PI / 2) * 0.01 * Arrow_Body_Width_Factor + (M_Module - 0.01 * Arrow_Head_length_Factor) * cos(Phi);
				Arrow_Body_Coordinate[i][j][2][1] = Grid_n_Min + j * Grid_Interval_n + sin(Phi - PI / 2) * 0.01 * Arrow_Body_Width_Factor + (M_Module - 0.01 * Arrow_Head_length_Factor) * sin(Phi);
			}
			Arrow_Body_Coordinate[i][j][2][2] = 0.;

			//top left
			if (M_Module == 0)
			{
				Arrow_Body_Coordinate[i][j][3][0] = Grid_m_Min + i * Grid_Interval_m;
				Arrow_Body_Coordinate[i][j][3][1] = Grid_n_Min + j * Grid_Interval_n;
			}else
			{
				Arrow_Body_Coordinate[i][j][3][0] = Grid_m_Min + i * Grid_Interval_m + cos(PI / 2 + Phi) * 0.01 * Arrow_Body_Width_Factor + (M_Module - 0.01 * Arrow_Head_length_Factor) * cos(Phi);
				Arrow_Body_Coordinate[i][j][3][1] = Grid_n_Min + j * Grid_Interval_n + sin(PI / 2 + Phi) * 0.01 * Arrow_Body_Width_Factor + (M_Module - 0.01 * Arrow_Head_length_Factor) * sin(Phi);
			}	
			Arrow_Body_Coordinate[i][j][3][2] = 0.;

			//Arrow Head
			//Top
			if(M_Module==0)
			{
				Arrow_Head_Coordinate[i][j][0][0] = 0.;
				Arrow_Head_Coordinate[i][j][0][1] = 0.;
				Arrow_Head_Coordinate[i][j][0][2] = 0.;
			}else
			{
				Arrow_Head_Coordinate[i][j][0][0] = Grid_m_Min + i * Grid_Interval_m + M_Module * cos(Phi);
				Arrow_Head_Coordinate[i][j][0][1] = Grid_n_Min + j * Grid_Interval_n + M_Module * sin(Phi);
				Arrow_Head_Coordinate[i][j][0][2] = 0.;
			}
			
			Arrow_Head_Bottom = 0.01 * Arrow_Head_length_Factor * sin(Arrow_Head_Angle/2);
			//Bottom left
			if(M_Module==0)
			{
				Arrow_Head_Coordinate[i][j][1][0] = 0.;
				Arrow_Head_Coordinate[i][j][1][1] = 0.;
				Arrow_Head_Coordinate[i][j][1][2] = 0.;
			}else
			{
				Arrow_Head_Coordinate[i][j][1][0] = Grid_m_Min + i * Grid_Interval_m + (M_Module - 0.01 * Arrow_Head_length_Factor) * cos(Phi) + Arrow_Head_Bottom * cos(PI / 2 + Phi);
				Arrow_Head_Coordinate[i][j][1][1] = Grid_n_Min + j * Grid_Interval_n + (M_Module - 0.01 * Arrow_Head_length_Factor) * sin(Phi) + Arrow_Head_Bottom * sin(PI / 2 + Phi);
				Arrow_Head_Coordinate[i][j][1][2] = 0.;
			}

			//Bottom Right
			if(M_Module==0)
			{
				Arrow_Head_Coordinate[i][j][2][0] = 0.;
				Arrow_Head_Coordinate[i][j][2][1] = 0.;
				Arrow_Head_Coordinate[i][j][2][2] = 0.;
			}else
			{
				Arrow_Head_Coordinate[i][j][2][0] = Grid_m_Min + i * Grid_Interval_m + (M_Module - 0.01 * Arrow_Head_length_Factor) * cos(Phi) + Arrow_Head_Bottom * cos(-PI / 2 + Phi);
				Arrow_Head_Coordinate[i][j][2][1] = Grid_n_Min + j * Grid_Interval_n + (M_Module - 0.01 * Arrow_Head_length_Factor) * sin(Phi) + Arrow_Head_Bottom * sin(-PI / 2 + Phi);
				Arrow_Head_Coordinate[i][j][2][2] = 0.;
			}
		}
			
	}

	return 1;
}

int Update_Arrow_Vertex(double ****Arrow_Head_Coordinate, double ****Arrow_Body_Coordinate, float *Arrow_Head_Vertex, float *Arrow_Body_Vertex, int i, int j)
{
	int k, r;
	for (k = 0; k < 3; ++k)
	{
		for (r = 0; r < 3; ++r)
		{
			Arrow_Head_Vertex[3*k+r] = Arrow_Head_Coordinate[i][j][k][r];
		}	
	}

	for (k = 0; k < 4; ++k)
	{
		for (r = 0; r < 3; ++r)
		{
			Arrow_Body_Vertex[3*k+r] = Arrow_Body_Coordinate[i][j][k][r];
		}	
	}

	return 1;
}

int OPGL_NMD(OPGL_Parameters *OPM)
{

	int m, n, l, Inplane_Density, i, j, k, N_Total = 0, Level, N_Total_Arrow, D, Scale;
	double ***Mx, ***My, ***Mz;
	double Energy, SKN;
	int Thread_Index_OPGL, Signal_RunAndStop;
	Thread_Index_OPGL = OPM->IPN->Threads_Number + 1;

	m = 2 * OPM->IPN->xn + 3;
	n = 2 * OPM->IPN->yn + 3;
	if(OPM->IPN->zn==0)
    {
        l=1;
    }else
    {
        l = 2*OPM->IPN->zn+3;
    }
	Inplane_Density = OPM->IPN->Inplane_Density;
	Mx = OPM->Mx_OPGL;
	My = OPM->My_OPGL;
	Mz = OPM->Mz_OPGL;

	///////////////////////// Boundary Parameters Initial /////////////////////////
    Boundary_Shape BS0, *BS;
    Calculate_Points_Parameters CPP0, *CPP;

    BS=&BS0;
    BS->xn=OPM->IPN->xn;  BS->yn=OPM->IPN->yn;
    BS->zn=OPM->IPN->zn;
    BS->Boundary    = Make3DArrayinteger(m,n,l);
    BS->LocalEnergy = Make4DArrayinteger(m,n,l,8);
    Boundary_Initial_Type(OPM->IPN,BS);
    LocalEnergy_Initial(BS);
    CalculPoints_Numbers(BS);
    //printf("%d\n",Calculate_Points);

    CPP=&CPP0;
    CPP->Position_3DTo1D  = Make3DArrayinteger(m,n,l);
    CPP->Local_Points_1D  = Make3DArrayinteger(BS->Calculate_Points,7,7);
    CPP->Virtual_Position = Make1DArrayinteger(BS->Virtual_Points);
    CPP->Local_Points_1D_LocalEnergy       = Make4DArrayinteger(BS->Calculate_Points,7,7,7);
    CPP->Neighbour_Points_3Dto1D_Array     = Make2DArrayinteger(14,2);
    CPP->Sort_Calculate_Points             = Make1DArrayinteger(BS->Calculate_Points);
    CalculPointsParameters_Initial(CPP,BS);

	free3DArrayinteger(BS->Boundary,m,n);
    free4DArrayinteger(BS->LocalEnergy,m,n,l);
    free(CPP->Virtual_Position);
    free4DArrayinteger(CPP->Local_Points_1D_LocalEnergy,BS->Calculate_Points,7,7);
    free2DArrayinteger(CPP->Neighbour_Points_3Dto1D_Array,14);

	Energy_Variables_Distribution EVD;
    EVD.CST = OPM->CST;    EVD.dxy = OPM->IPN->dxy;
    EVD.h   = OPM->IPN->h; EVD.t   = OPM->IPN->t;
    EVD.m   = m;      EVD.n   = n;
    EVD.l   = l;
    EVD.Mx    = Mx;
    EVD.My    = My;
    EVD.Mz    = Mz;
    EVD.MxV   = OPM->MxV_OPGL;
    EVD.MyV   = OPM->MyV_OPGL;
    EVD.MzV   = OPM->MzV_OPGL;
    EVD.w     = Make3DArray(m,n,l);
    EVD.SKN_D = Make3DArray(m,n,l);
    EVD.Local_Points_1D=CPP->Local_Points_1D;
    EVD.Position_3DTo1D=CPP->Position_3DTo1D;
	EVD.Energy_Coefficient = Make1DArray(N_Energy_Terms);

	for (i = 0; i < N_Energy_Terms; ++i)
	{
		EVD.Energy_Coefficient[i] = OPM->IPN->Energy_Coefficient[i];
	}
	/////////////////////////////////////////////////////////////////////////////////////////
	double Grid_m_Max, Grid_m_Min, Grid_n_Max, Grid_n_Min, *GridSize;
	GridSize = Make1DArray(4);
	SetGridSize(m, n, GridSize);
	Grid_m_Max = GridSize[0];
	Grid_m_Min = GridSize[1];
	Grid_n_Max = GridSize[2];
	Grid_n_Min = GridSize[3];
	//printf("%lf %lf %lf %lf\n", GridSize[0], GridSize[1], GridSize[2], GridSize[3]);

	Level = 45;
	Scale = 5;
	float *MzMaxMin;
    MzMaxMin = Make1DArray_float(Scale+1);
	for (i = 0; i <= Scale;++i)
	{
		MzMaxMin[i] = 0.;
	}
	float *ScalePositionMaxMin;
    ScalePositionMaxMin = Make1DArray_float(Scale+1);
	for (i = 0; i <= Scale;++i)
	{
		ScalePositionMaxMin[i] = 860.-i*(860.-120.)/Scale;
	}

	float RGB_Red, RGB_Green, RGB_Blue;
	std::string s1, s2, s3, s4, s5, s6, s7, s8;

	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			N_Total = N_Total + 1;
		}
	}

	// glfw: initialize and configure
	// ------------------------------
	glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    // glfw window creation
    // --------------------
    GLFWwindow* window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

    // glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

////////////////////////////////////////////////////////////////////////////////////
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

//////////////////////////////////////Freetype//////////////////////////////////////
	Shader shader("text.vs", "text.fs");
    glm::mat4 projection = glm::ortho(0.0f, static_cast<float>(SCR_WIDTH), 0.0f, static_cast<float>(SCR_HEIGHT));
    shader.use();
    // FreeType
    FT_Library ft;
    // All functions return a value different than 0 whenever an error occurred
    if (FT_Init_FreeType(&ft))
    {
        std::cout << "ERROR::FREETYPE: Could not init FreeType Library" << std::endl;
        return -1;
    }
	FT_Face face;

    if (FT_New_Face(ft, "fonts/arial.ttf", 0, &face)) {
        std::cout << "ERROR::FREETYPE: Failed to load font" << std::endl;
        return -1;
    }
    else {
        // set size to load glyphs as
        FT_Set_Pixel_Sizes(face, 0, 48);
        // disable byte-alignment restriction
        glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
        // load first 128 characters of ASCII set
        for (unsigned char c = 0; c < 128; c++)
        {
            // Load character glyph 
            if (FT_Load_Char(face, c, FT_LOAD_RENDER))
            {
                std::cout << "ERROR::FREETYTPE: Failed to load Glyph" << std::endl;
                continue;
            }
            // generate texture
            unsigned int texture;
            glGenTextures(1, &texture);
            glBindTexture(GL_TEXTURE_2D, texture);
            glTexImage2D(
                GL_TEXTURE_2D,
                0,
                GL_RED,
                face->glyph->bitmap.width,
                face->glyph->bitmap.rows,
                0,
                GL_RED,
                GL_UNSIGNED_BYTE,
                face->glyph->bitmap.buffer
            );
            // set texture options
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
            // now store character for later use
            Character character = {
                texture,
                glm::ivec2(face->glyph->bitmap.width, face->glyph->bitmap.rows),
                glm::ivec2(face->glyph->bitmap_left, face->glyph->bitmap_top),
                face->glyph->advance.x
            };
            Characters.insert(std::pair<char, Character>(c, character));
        }
        glBindTexture(GL_TEXTURE_2D, 0);
    }
	    // destroy FreeType once we're finished
    FT_Done_Face(face);
    FT_Done_FreeType(ft);
////////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////////
    // build and compile our shader program
    // ------------------------------------
    // vertex shader
    int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    glCompileShader(vertexShader);
    // check for shader compile errors
    int success;
    char infoLog[512];
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // fragment shader
    int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    glCompileShader(fragmentShader);
    // check for shader compile errors
    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    }
    // link shaders
    int shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glLinkProgram(shaderProgram);
    // check for linking errors
    glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    if (!success) {
        glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
        std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    }
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);


    unsigned int indices[] = {  // note that we start from 0!
        0, 1, 2,  // first Triangle
        0, 2, 3   // second Triangle
    };

	unsigned int indices_Arrow_Head[] = {  // note that we start from 0!
        0, 1, 2  // Triangle
    };

    float *Background_Rectangle_Vertex;
    Background_Rectangle_Vertex = Make1DArray_float(12);
	int **Background_Matrix;
	Background_Matrix = Make2DArrayinteger(m, n);
	double ****Background_Matrix_Vertex_Coordinate;
	Background_Matrix_Vertex_Coordinate = Make4DArray(m, n, 4, 3);
	Background_Matrix_Vertex_Coordinate_Assignment(Background_Matrix_Vertex_Coordinate, m, n, Grid_m_Max, Grid_m_Min, Grid_n_Max, Grid_n_Min);

	float *Scale_Rectangle_Vertex;
    Scale_Rectangle_Vertex = Make1DArray_float(12);
	int **Scale_Matrix;
	Scale_Matrix = Make2DArrayinteger(1, Scale);
	double ****Scale_Matrix_Vertex_Coordinate;
	Scale_Matrix_Vertex_Coordinate = Make4DArray(1, Scale, 4, 3);
	Scale_Matrix_Vertex_Coordinate_Assignment(Scale_Matrix_Vertex_Coordinate, Scale);


	int **RGB_Matrix;
	RGB_Matrix = Make2DArrayinteger(Level, 3);
	Pastel_RGB(RGB_Matrix, Level);

	int **RGB_Matrix_Scale;
	RGB_Matrix_Scale = Make2DArrayinteger(Scale, 3);
	Pastel_RGB(RGB_Matrix_Scale, Scale);

	int *OI_Parameter;
	OI_Parameter = Make1DArrayinteger(3);
	D = OPGL_InPlaneComponent_Parameter(m, n, Inplane_Density, OI_Parameter);
	//printf("%d\t%d\t%d\t%d\n", D, OI_Parameter[0], OI_Parameter[1], OI_Parameter[2]);

	double ****Arrow_Head_Coordinate, ****Arrow_Body_Coordinate;
	Arrow_Head_Coordinate = Make4DArray(OI_Parameter[1], OI_Parameter[2], 3, 3);
	Arrow_Body_Coordinate = Make4DArray(OI_Parameter[1], OI_Parameter[2], 4, 3);
	OPGL_InPlaneComponent_Coordiante(Mx, My, Mz, D, OI_Parameter, Inplane_Density, Arrow_Head_Coordinate, Arrow_Body_Coordinate, Grid_m_Max, Grid_m_Min, Grid_n_Max, Grid_n_Min);

	float *Arrow_Head_Vertex, *Arrow_Body_Vertex;
	Arrow_Head_Vertex = Make1DArray_float(9);
	Arrow_Body_Vertex = Make1DArray_float(12);


	unsigned int *VBOId_Background, *VAOId_Background, *EBOId_Background, *VBOId_Arrow_Body, *VAOId_Arrow_Body,*EBOId_Arrow_Body;
	unsigned int *VBOId_Arrow_Head, *VAOId_Arrow_Head,*EBOId_Arrow_Head, *VBOId_Scale, *VAOId_Scale, *EBOId_Scale;
	unsigned int VAOId_Text, VBOId_Text;
	VBOId_Background = (unsigned int *)malloc(sizeof(unsigned int) * N_Total);
	VAOId_Background = (unsigned int *)malloc(sizeof(unsigned int) * N_Total);
	EBOId_Background = (unsigned int *)malloc(sizeof(unsigned int) * N_Total);

	VBOId_Arrow_Body = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);
	VAOId_Arrow_Body = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);
	EBOId_Arrow_Body = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);

	VBOId_Arrow_Head = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);
	VAOId_Arrow_Head = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);
	EBOId_Arrow_Head = (unsigned int *)malloc(sizeof(unsigned int) * OI_Parameter[0]);

	VBOId_Scale = (unsigned int *)malloc(sizeof(unsigned int) * Scale);
	VAOId_Scale = (unsigned int *)malloc(sizeof(unsigned int) * Scale);
	EBOId_Scale = (unsigned int *)malloc(sizeof(unsigned int) * Scale);

	glGenVertexArrays(N_Total, VAOId_Background); // we can also generate multiple VAOs or buffers at the same time
	glGenBuffers(N_Total, VBOId_Background);
	glGenBuffers(N_Total, EBOId_Background);

	k = 0;
	for (i = 0; i < m; ++i)
	{
		for (j = 0; j < n; ++j)
		{
			Update_Background_Rectangle_Vertex(Background_Matrix_Vertex_Coordinate, Background_Rectangle_Vertex, i, j);
			glBindVertexArray(VAOId_Background[k]);
    		glBindBuffer(GL_ARRAY_BUFFER, VBOId_Background[k]);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Background[k]);
			glBufferData(GL_ARRAY_BUFFER, 12*sizeof(float), NULL, GL_STREAM_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, 12*sizeof(float), Background_Rectangle_Vertex);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STREAM_DRAW);
   			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    		glEnableVertexAttribArray(0);
    
    		glBindBuffer(GL_ARRAY_BUFFER, 0);

			k = k + 1;
		}
	}

//////////////////////////////////////////ScaleText/////////////////////////////////////////////////////////////
    glGenVertexArrays(1, &VAOId_Text);
    glGenBuffers(1, &VBOId_Text);
    glBindVertexArray(VAOId_Text);
    glBindBuffer(GL_ARRAY_BUFFER, VBOId_Text);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 6 * 4, NULL, GL_DYNAMIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, 4 * sizeof(float), 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
    glUniformMatrix4fv(glGetUniformLocation(shader.ID, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

//////////////////////////////////////////Scale/////////////////////////////////////////////////////////////

    glGenVertexArrays(Scale, VAOId_Scale);
	glGenBuffers(Scale, VBOId_Scale);
	glGenBuffers(Scale, EBOId_Scale);

    k = 0;
	for (i = 0; i < 1; ++i)
	{
		for (j = 0; j < Scale; ++j)
		{
			Update_Scale_Rectangle_Vertex(Scale_Matrix_Vertex_Coordinate, Scale_Rectangle_Vertex, i, j);
			glBindVertexArray(VAOId_Scale[k]);
    		glBindBuffer(GL_ARRAY_BUFFER, VBOId_Scale[k]);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Scale[k]);
			glBufferData(GL_ARRAY_BUFFER, 12*sizeof(float), NULL, GL_STREAM_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, 12*sizeof(float), Scale_Rectangle_Vertex);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STREAM_DRAW);
   			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    		glEnableVertexAttribArray(0);
    		glBindBuffer(GL_ARRAY_BUFFER, 0);

			k = k + 1;
		}
	}
//////////////////////////////////////////Arrow/////////////////////////////////////////////////////////////

	glGenVertexArrays(OI_Parameter[0], VAOId_Arrow_Head); // we can also generate multiple VAOs or buffers at the same time
	glGenBuffers(OI_Parameter[0], VBOId_Arrow_Head);
	glGenBuffers(OI_Parameter[0], EBOId_Arrow_Head);

	glGenVertexArrays(OI_Parameter[0], VAOId_Arrow_Body); // we can also generate multiple VAOs or buffers at the same time
	glGenBuffers(OI_Parameter[0], VBOId_Arrow_Body);
	glGenBuffers(OI_Parameter[0], EBOId_Arrow_Body);


	k = 0;
	for (i = 0; i < OI_Parameter[1]; ++i)
	{
		for (j = 0; j < OI_Parameter[2]; ++j)
		{
			Update_Arrow_Vertex(Arrow_Head_Coordinate, Arrow_Body_Coordinate, Arrow_Head_Vertex, Arrow_Body_Vertex, i, j);
			//Arrow_Head
			glBindVertexArray(VAOId_Arrow_Head[k]);
    		glBindBuffer(GL_ARRAY_BUFFER, VBOId_Arrow_Head[k]);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Arrow_Head[k]);
			glBufferData(GL_ARRAY_BUFFER, 9*sizeof(float), NULL, GL_STREAM_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, 9*sizeof(float), Arrow_Head_Vertex);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices_Arrow_Head), indices_Arrow_Head, GL_STREAM_DRAW);
   			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    		glEnableVertexAttribArray(0);
    		glBindBuffer(GL_ARRAY_BUFFER, 0);

			//Arrow_Body
			glBindVertexArray(VAOId_Arrow_Body[k]);
    		glBindBuffer(GL_ARRAY_BUFFER, VBOId_Arrow_Body[k]);
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Arrow_Body[k]);
			glBufferData(GL_ARRAY_BUFFER, 12*sizeof(float), NULL, GL_STREAM_DRAW);
			glBufferSubData(GL_ARRAY_BUFFER, 0, 12*sizeof(float), Arrow_Body_Vertex);
			glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STREAM_DRAW);
   			glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    		glEnableVertexAttribArray(0);
    		glBindBuffer(GL_ARRAY_BUFFER, 0);

			k = k + 1;
		}
	}

    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window);

		pthread_mutex_lock(&OPGL_Mutex);
        Background_Matrix_Assignment(Mz, Background_Matrix, m, n);
		pthread_mutex_unlock(&OPGL_Mutex);
        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);


		k = 0;
		for (i = 0; i < m; ++i)
		{
			for (j = 0; j < n; ++j)
			{

				RGB_Red   = RGB_Matrix[Background_Matrix[i][j]][0] / 255.f;
				RGB_Green = RGB_Matrix[Background_Matrix[i][j]][1] / 255.f;
				RGB_Blue  = RGB_Matrix[Background_Matrix[i][j]][2] / 255.f;

				glUseProgram(shaderProgram);
				int vertexColorLocation = glGetUniformLocation(shaderProgram, "ourColor");
				glBindVertexArray(VAOId_Background[k]);
       			glUniform4f(vertexColorLocation, RGB_Red, RGB_Green, RGB_Blue, 1.0f);
       			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
				k = k + 1;

			}
		}

        pthread_mutex_lock(&OPGL_Mutex);
        Scale_Matrix_Assignment(Mz, Scale_Matrix, m, n, Scale);
		pthread_mutex_unlock(&OPGL_Mutex);

		k = 0;
		for (i = 0; i < 1; ++i)
		{
			for (j = 0; j < Scale; ++j)
			{

				RGB_Red   = RGB_Matrix_Scale[Scale_Matrix[i][j]][0] / 255.f;
				RGB_Green = RGB_Matrix_Scale[Scale_Matrix[i][j]][1] / 255.f;
				RGB_Blue  = RGB_Matrix_Scale[Scale_Matrix[i][j]][2] / 255.f;

				glUseProgram(shaderProgram);
				int vertexColorLocation = glGetUniformLocation(shaderProgram, "ourColor");
				glBindVertexArray(VAOId_Scale[k]);
       			glUniform4f(vertexColorLocation, RGB_Red, RGB_Green, RGB_Blue, 1.0f);
       			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
				k = k + 1;

			}
		}

		pthread_mutex_lock(&OPGL_Mutex);
		OPGL_InPlaneComponent_Coordiante(Mx, My, Mz, D, OI_Parameter, Inplane_Density, Arrow_Head_Coordinate, Arrow_Body_Coordinate, Grid_m_Max, Grid_m_Min, Grid_n_Max, Grid_n_Min);
		pthread_mutex_unlock(&OPGL_Mutex);
		k = 0;
		for (i = 0; i < OI_Parameter[1]; ++i)
		{
			for (j = 0; j < OI_Parameter[2]; ++j)
			{

				RGB_Red   = 0.f;
				RGB_Green = 0.f;
				RGB_Blue  = 0.f;

				Update_Arrow_Vertex(Arrow_Head_Coordinate, Arrow_Body_Coordinate, Arrow_Head_Vertex, Arrow_Body_Vertex, i, j);
				//Arrow_Head
				glBindVertexArray(VAOId_Arrow_Head[k]);
    			glBindBuffer(GL_ARRAY_BUFFER, VBOId_Arrow_Head[k]);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Arrow_Head[k]);
				glBufferData(GL_ARRAY_BUFFER, 9*sizeof(float), NULL, GL_STREAM_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, 9*sizeof(float), Arrow_Head_Vertex);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices_Arrow_Head), indices_Arrow_Head, GL_STREAM_DRAW);
   				glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    			glEnableVertexAttribArray(0);
    			glBindBuffer(GL_ARRAY_BUFFER, 0);

				//Arrow_Body
				glBindVertexArray(VAOId_Arrow_Body[k]);
    			glBindBuffer(GL_ARRAY_BUFFER, VBOId_Arrow_Body[k]);
				glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOId_Arrow_Body[k]);
				glBufferData(GL_ARRAY_BUFFER, 12*sizeof(float), NULL, GL_STREAM_DRAW);
				glBufferSubData(GL_ARRAY_BUFFER, 0, 12*sizeof(float), Arrow_Body_Vertex);
				glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STREAM_DRAW);
   				glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0); // because the vertex data is tightly packed we can also specify 0 as the vertex attribute's stride to let OpenGL figure it out
    			glEnableVertexAttribArray(0);
    			glBindBuffer(GL_ARRAY_BUFFER, 0);

				glUseProgram(shaderProgram);
				int vertexColorLocation = glGetUniformLocation(shaderProgram, "ourColor");
				glBindVertexArray(VAOId_Arrow_Head[k]);
       			glUniform4f(vertexColorLocation, RGB_Red, RGB_Green, RGB_Blue, 1.0f);
       			glDrawElements(GL_TRIANGLES, 3, GL_UNSIGNED_INT, 0);

				glBindVertexArray(VAOId_Arrow_Body[k]);
       			glUniform4f(vertexColorLocation, RGB_Red, RGB_Green, RGB_Blue, 1.0f);
       			glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
				k = k + 1;

			}
		}

		pthread_mutex_lock(&OPGL_Mutex);
		//Energy_Distribution(&EVD);
		//Energy = EVD.W_Average;
		//SKN = EVD.SKN;
		Energy = OPGL_GUI.Energy;
		SKN    = OPGL_GUI.SKN;
		Get_Max_Min(Mz, MzMaxMin, m, n, Scale);
		pthread_mutex_unlock(&OPGL_Mutex);

		s1 = std::to_string(MzMaxMin[0]);
		s2 = std::to_string(MzMaxMin[1]);
		s3 = std::to_string(MzMaxMin[2]);
		s4 = std::to_string(MzMaxMin[3]);
		s5 = std::to_string(MzMaxMin[4]);
		s6 = std::to_string(MzMaxMin[5]);
		s7 = "Energy: " + std::to_string(Energy);
		s8 = "SKN: " + std::to_string(SKN);

		//RenderText(shader, "Mz", 1000.0f, 1000.0f, 1.0f, glm::vec3(0.5, 0.8f, 0.2f));
        RenderText(shader, "Mz", 825.0f, 890.0f, 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s1, 870.0f, ScalePositionMaxMin[0], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s2, 870.0f, ScalePositionMaxMin[1], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s3, 870.0f, ScalePositionMaxMin[2], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s4, 870.0f, ScalePositionMaxMin[3], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s5, 870.0f, ScalePositionMaxMin[4], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s6, 870.0f, ScalePositionMaxMin[5], 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s7, 100.0f, 890.0f, 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
		RenderText(shader, s8, 500.0f, 890.0f, 0.5f, glm::vec3(0.3, 0.7f, 0.9f), VAOId_Text, VBOId_Text);
        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
    	glfwSwapBuffers(window);
        glfwPollEvents();

		/////////////////////////// Suspend thread by GUI /////////////////////////////
		/*
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Thread_Index_OPGL]);
        while (GPM.Thread_State[Thread_Index_OPGL][1] == Signal_GUI_Break)//Signal_ContinueAndBreak
        {
            pthread_cond_wait(&GPM.Pthread_GUI_Cond[Thread_Index_OPGL], &GPM.Pthread_GUI_Mutex[Thread_Index_OPGL]);
        }
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Thread_Index_OPGL]);
		*/

		/////////////////////////// Stop thread by GUI /////////////////////////////
        pthread_mutex_lock(&GPM.Pthread_GUI_Mutex[Thread_Index_OPGL]);
        Signal_RunAndStop = GPM.Thread_State[Thread_Index_OPGL][0];
        pthread_mutex_unlock(&GPM.Pthread_GUI_Mutex[Thread_Index_OPGL]);
        if (Signal_RunAndStop == Signal_GUI_Stop)
        {
            break;
        }

		/////////////////////////// Stop thread by Synchronize thread /////////////////////////////
		pthread_mutex_lock(&OPGL_Mutex);
		if (OPGL_GUI.State_Syn == StateSyn_Stop)
		{
            break;
        }
		pthread_mutex_unlock(&OPGL_Mutex);

		//Sleep(20);
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------	
	glDeleteVertexArrays(N_Total, VAOId_Background);
    glDeleteBuffers(N_Total, VBOId_Background);
	glDeleteBuffers(N_Total, EBOId_Background);

    glDeleteVertexArrays(Scale, VAOId_Scale);
    glDeleteBuffers(Scale, VBOId_Scale);
	glDeleteBuffers(Scale, EBOId_Scale);

	glDeleteVertexArrays(OI_Parameter[0], VAOId_Arrow_Body);
    glDeleteBuffers(OI_Parameter[0], VBOId_Arrow_Body);
	glDeleteBuffers(OI_Parameter[0], EBOId_Arrow_Body);

	glDeleteVertexArrays(OI_Parameter[0], VAOId_Arrow_Head);
    glDeleteBuffers(OI_Parameter[0], VBOId_Arrow_Head);
	glDeleteBuffers(OI_Parameter[0], EBOId_Arrow_Head);
    glDeleteProgram(shaderProgram);

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

	free(VBOId_Background);
	free(VAOId_Background);
	free(EBOId_Background);
    free(VBOId_Scale);
	free(VAOId_Scale);
	free(EBOId_Scale);
	free(VBOId_Arrow_Head);
	free(VAOId_Arrow_Head);
	free(EBOId_Arrow_Head);
	free(VBOId_Arrow_Body);
	free(VAOId_Arrow_Body);
	free(EBOId_Arrow_Body);
	free(Background_Rectangle_Vertex);
    free(Scale_Rectangle_Vertex);
	free(Arrow_Head_Vertex);
	free(Arrow_Body_Vertex);
	free(GridSize);
	free2DArrayinteger(RGB_Matrix,Level);
    free4DArray(Background_Matrix_Vertex_Coordinate,m,n,4);
    free4DArray(Scale_Matrix_Vertex_Coordinate,1,Scale,4);
	free4DArray(Arrow_Head_Coordinate, OI_Parameter[1], OI_Parameter[2], 3);
	free4DArray(Arrow_Body_Coordinate, OI_Parameter[1], OI_Parameter[2], 4);
	free(MzMaxMin);
	free(ScalePositionMaxMin);
	free3DArray(EVD.w,m,n);
    free3DArray(EVD.SKN_D,m,n);
	free(EVD.Energy_Coefficient);
	free3DArrayinteger(CPP->Local_Points_1D,BS->Calculate_Points,7);
    free3DArrayinteger(CPP->Position_3DTo1D,m,n);
	return 1;
}

void RenderText(Shader &shader, std::string text, float x, float y, float scale, glm::vec3 color, unsigned int VAOId_Text, unsigned int VBOId_Text)
{
    // activate corresponding render state	
    shader.use();
    glUniform3f(glGetUniformLocation(shader.ID, "textColor"), color.x, color.y, color.z);
    glActiveTexture(GL_TEXTURE0);
    glBindVertexArray(VAOId_Text);
    // iterate through all characters
    std::string::const_iterator c;
    for (c = text.begin(); c != text.end(); c++) 
    {
        Character ch = Characters[*c];
        float xpos = x + ch.Bearing.x * scale;
        float ypos = y - (ch.Size.y - ch.Bearing.y) * scale;
        float w = ch.Size.x * scale;
        float h = ch.Size.y * scale;
        // update VBO for each character
        float vertices[6][4] = {
            { xpos,     ypos + h,   0.0f, 0.0f },            
            { xpos,     ypos,       0.0f, 1.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos,     ypos + h,   0.0f, 0.0f },
            { xpos + w, ypos,       1.0f, 1.0f },
            { xpos + w, ypos + h,   1.0f, 0.0f }           
        };
        // render glyph texture over quad
        glBindTexture(GL_TEXTURE_2D, ch.TextureID);
        // update content of VBO memory
        glBindBuffer(GL_ARRAY_BUFFER, VBOId_Text);
        glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(vertices), vertices); // be sure to use glBufferSubData and not glBufferData
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        // render quad
        glDrawArrays(GL_TRIANGLES, 0, 6);
        // now advance cursors for next glyph (note that advance is number of 1/64 pixels)
        x += (ch.Advance >> 6) * scale; // bitshift by 6 to get value in pixels (2^6 = 64 (divide amount of 1/64th pixels by 64 to get amount of pixels))
    }
    glBindVertexArray(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// ////////////////////////////////GUI/////////////////////////////////////
// void TW_CALL CB_Set_Run( const void *value, void *clientData )
// {
// 	Play = *( int *)value;
//     if (Play!=0){
//         pthread_mutex_lock(&culc_mutex);
//             ENGINE_MUTEX=DO_IT;
//             SleepTime=100;
//         pthread_mutex_unlock(&culc_mutex);
//     }else{
//         pthread_mutex_lock(&culc_mutex);
//             ENGINE_MUTEX=WAIT;
//             SleepTime=3000;
//         pthread_mutex_unlock(&culc_mutex);  
//     }
// }

// void TW_CALL CB_Get_Run(void *value, void *clientData)
// {
//     *(float *)value = Play;
// }