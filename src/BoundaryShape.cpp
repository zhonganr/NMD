#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Array.h>
#include <BoundaryShape.h>
#include <Variable.h>

#define pi   3.14159265358

using namespace std;
//This function is to initialize the Circle/Cylinder Boundary parameters.
//rn is the radius of Circle/Cylinder.
int Boundary_Ellipse(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    double D, an, bn;
    //D: Distance between origin(xn, yn, zn) and the selected point (i, j, k).
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    an=yn;
    bn=xn;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    int xc=(m-1)/2, yc=(n-1)/2;//Circle Center Radius.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                D=(i-xc)*(i-xc)/bn/bn+(j-yc)*(j-yc)/an/an;
                if(D<=1.)
                {
                    BS->Boundary[i][j][k]=1;//1 is inner side.
                    if(l==1)//2 Dimensions Case
                    {
                        BS->Boundary[i][j][0]=1;//1 is inner side.
                    }else   //3 Dimensions Case, the upper and the lower are both virtual points.
                    {
                        BS->Boundary[i][j][0]  =0;//0 represents exterior side.
                        BS->Boundary[i][j][l-1]=0;//0 represents exterior side.
                    }
                }else
                {
                    BS->Boundary[i][j][k]=0;//0 represents exterior side.
                }
            }
        }
    }


    for(k=0;k<l;++k)
    {
        for(i=1;i<m-1;++i)
        {
            for(j=1;j<n-1;++j)
            {
                if(BS->Boundary[i][j][k]==1)
                {
                    if(BS->Boundary[i+1][j][k]*BS->Boundary[i-1][j][k]*BS->Boundary[i][j-1][k]*BS->Boundary[i][j+1][k]==0)
                    {
                        BS->Boundary[i][j][k]=2;//2 is on the boundary.
                    }

                    if(l>1)//3 Dimensions Case
                    {
                        BS->Boundary[i][j][1]  =2;//2 is on the boundary.
                        BS->Boundary[i][j][l-2]=2;//2 is on the boundary.
                    }
                }
            }
        }
    }

    return 1;

}

int Boundary_Circle(int rn, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, rn_Min;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    double D;
    //D: Distance between origin(xn, yn, zn) and the selected point (i, j, k).
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    //Obtain the minimum radius.
    if(xn>yn)
    {
        rn_Min=yn;
    }else
    {
        rn_Min=xn;
    }

    if(rn<rn_Min)
    {
        rn_Min=rn;
    }

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    int xc=(m-1)/2, yc=(n-1)/2;//Circle Center Radius.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                D=pow((i-xc)*(i-xc)+(j-yc)*(j-yc),0.5);
                if(D<=rn_Min)
                {
                    BS->Boundary[i][j][k]=1;//1 is inner side.
                    if(l==1)//2 Dimensions Case
                    {
                        BS->Boundary[i][j][0]=1;//1 is inner side.
                    }else   //3 Dimensions Case, the upper and the lower are both virtual points.
                    {
                        BS->Boundary[i][j][0]  =0;//0 represents exterior side.
                        BS->Boundary[i][j][l-1]=0;//0 represents exterior side.
                    }
                }else
                {
                    BS->Boundary[i][j][k]=0;//0 represents exterior side.
                }
            }
        }
    }


    for(k=0;k<l;++k)
    {
        for(i=1;i<m-1;++i)
        {
            for(j=1;j<n-1;++j)
            {
                if(BS->Boundary[i][j][k]==1)
                {
                    if(BS->Boundary[i+1][j][k]*BS->Boundary[i-1][j][k]*BS->Boundary[i][j-1][k]*BS->Boundary[i][j+1][k]==0)
                    {
                        BS->Boundary[i][j][k]=2;//2 is on the boundary.
                    }

                    if(l>1)//3 Dimensions Case
                    {
                        BS->Boundary[i][j][1]  =2;//2 is on the boundary.
                        BS->Boundary[i][j][l-2]=2;//2 is on the boundary.
                    }
                }
            }
        }
    }


    return 1;

}

int Boundary_Circle_Hole(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, rn_Min, rn_Max;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    double D;
    //D: Distance between origin(xn, yn, zn) and the selected point (i, j, k).
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    //Obtain the minimum radius.
    rn_Max=xn;
    rn_Min=xn*0.5;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    int xc=(m-1)/2, yc=(n-1)/2;//Circle Center Radius.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                D=pow((i-xc)*(i-xc)+(j-yc)*(j-yc),0.5);
                if(D<=rn_Max&&D>=rn_Min)
                {
                    BS->Boundary[i][j][k]=1;//1 is inner side.
                    if(l==1)//2 Dimensions Case
                    {
                        BS->Boundary[i][j][0]=1;//1 is inner side.
                    }else   //3 Dimensions Case, the upper and the lower are both virtual points.
                    {
                        BS->Boundary[i][j][0]  =0;//0 represents exterior side.
                        BS->Boundary[i][j][l-1]=0;//0 represents exterior side.
                    }
                }else
                {
                    BS->Boundary[i][j][k]=0;//0 represents exterior side.
                }
            }
        }
    }


    for(k=0;k<l;++k)
    {
        for(i=1;i<m-1;++i)
        {
            for(j=1;j<n-1;++j)
            {
                if(BS->Boundary[i][j][k]==1)
                {
                    if(BS->Boundary[i+1][j][k]*BS->Boundary[i-1][j][k]*BS->Boundary[i][j-1][k]*BS->Boundary[i][j+1][k]==0)
                    {
                        BS->Boundary[i][j][k]=2;//2 is on the boundary.
                    }

                    if(l>1)//3 Dimensions Case
                    {
                        BS->Boundary[i][j][1]  =2;//2 is on the boundary.
                        BS->Boundary[i][j][l-2]=2;//2 is on the boundary.
                    }
                }
            }
        }
    }


    return 1;

}

int Regular_Polygon_Plane_Shape(int rn, int N_Polygon, double **Line_Polygon)
{
    int i, N1, N2;
    double **Vertex_Polar, **Vertex_Cartesian;

    Vertex_Polar     = Make2DArray(N_Polygon,3);
    Vertex_Cartesian = Make2DArray(N_Polygon,3);
    for(i=0;i<N_Polygon;++i)
    {
        Vertex_Polar[i][0] = rn;//Circumradius
        Vertex_Polar[i][1] = 0; //Theta
        Vertex_Polar[i][2] = -pi+2*pi/N_Polygon*i;//Phi

        Vertex_Cartesian[i][0] = Vertex_Polar[i][0]*cos(Vertex_Polar[i][1])*cos(Vertex_Polar[i][2]);//x
        Vertex_Cartesian[i][1] = Vertex_Polar[i][0]*cos(Vertex_Polar[i][1])*sin(Vertex_Polar[i][2]);//y
        Vertex_Cartesian[i][2] = Vertex_Polar[i][0]*sin(Vertex_Polar[i][1]);//z
    }
    ////////////// Line Parameter: Ax+By+C=0 ////////////////
    for(i=0;i<N_Polygon;++i)
    {
        N1 = i;
        N2 = i+1;
        if(N2==N_Polygon)
        {
            N2=0;
        }
        Line_Polygon[i][0] = Vertex_Cartesian[N1][1]-Vertex_Cartesian[N2][1];//A
        Line_Polygon[i][1] = Vertex_Cartesian[N1][0]-Vertex_Cartesian[N2][0];//B
        Line_Polygon[i][2] = Vertex_Cartesian[N1][0]*Vertex_Cartesian[N2][1]-Vertex_Cartesian[N2][0]*Vertex_Cartesian[N1][1];//C
    }

    free2DArray(Vertex_Polar,N_Polygon);
    free2DArray(Vertex_Cartesian,N_Polygon);
    return 1;
}

int Boundary_Rectangle(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }

    return 1;

}

int Boundary_Rectangle_20Grid_Crack(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc;
    int Grids=200;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }

    for(i=-Grids;i<=Grids;++i)
    {
        BS->Boundary[xc+Grids][yc][0]=0;
        BS->Boundary[xc][yc-1][0]=2;
        BS->Boundary[xc][yc+1][0]=2;
    }
    BS->Boundary[xc-Grids-1][yc][0]=2;
    BS->Boundary[xc+Grids+1][yc][0]=2;

    return 1;

}

int Boundary_Rectangle_Circle(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc, rn, r;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    rn = 100;
    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            r = sqrt((i-xc)*(i-xc)+(j-yc)*(j-yc));
            if(r<rn)
            {
                BS->Boundary[i][j][0]=0;
            }
        }
    }

    for(i=xc-rn-1;i<xc+rn+1;++i)
    {
        for(j=yc-rn-1;j<yc+rn+1;++j)
        {
            if(BS->Boundary[i][j][0]!=0)
            {
                if(BS->Boundary[i-1][j][0]*BS->Boundary[i+1][j][0]*BS->Boundary[i][j-1][0]*BS->Boundary[i][j+1][0]==0)
                {
                    BS->Boundary[i][j][0]=2;
                }
            }
        }
    }


    return 1;

}

int Boundary_Rectangle_Ellipse(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc;
    double r, r_a, r_b, an, bn;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    an = 100;   bn=50;
    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            r_a = (i-xc)*(i-xc)/(bn*bn);
            r_b = (j-yc)*(j-yc)/(an*an);
            r = r_a+r_b;
            //printf("%0.3f\n",r);
            if(r<1)
            {
                BS->Boundary[i][j][0]=0;
            }
        }
    }

    for(i=xc-bn-1;i<xc+bn+1;++i)
    {
        for(j=yc-an-1;j<yc+an+1;++j)
        {
            if(BS->Boundary[i][j][0]!=0)
            {
                if(BS->Boundary[i-1][j][0]*BS->Boundary[i+1][j][0]*BS->Boundary[i][j-1][0]*BS->Boundary[i][j+1][0]==0)
                {
                    BS->Boundary[i][j][0]=2;
                }
            }
        }
    }


    return 1;

}

int Boundary_Rectangle_4Circle_Holes(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc, rn, p;
    int xc1, yc1, xc2, yc2, xc3, yc3, xc4, yc4, **Center;
    double r;
    Center=Make2DArrayinteger(4,2);
    xc1 = 151; yc1 = 151;
    xc2 = 452; yc2 = 151;
    xc3 = 151; yc3 = 452;
    xc4 = 452; yc4 = 452;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    Center[0][0] = xc1; Center[0][1] = yc1;
    Center[1][0] = xc2; Center[1][1] = yc2;
    Center[2][0] = xc3; Center[2][1] = yc3;
    Center[3][0] = xc4; Center[3][1] = yc4;

    rn = 50;
    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }



    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(p=0;p<4;++p)
            {
                r = sqrt((i-Center[p][0])*(i-Center[p][0])+(j-Center[p][1])*(j-Center[p][1]));
                if(r<rn)
                {
                    BS->Boundary[i][j][0]=0;
                }
            }
        }
    }

    for(p=0;p<4;++p)
    {
        for(i=Center[p][0]-rn-1;i<Center[p][0]+rn+1;++i)
        {
            for(j=Center[p][1]-rn-1;j<Center[p][1]+rn+1;++j)
            {
                if(BS->Boundary[i][j][0]!=0)
                {
                    if(BS->Boundary[i-1][j][0]*BS->Boundary[i+1][j][0]*BS->Boundary[i][j-1][0]*BS->Boundary[i][j+1][0]==0)
                    {
                        BS->Boundary[i][j][0]=2;
                    }
                }
            }
        }
    }

    free2DArrayinteger(Center,4);
    return 1;

}

int Boundary_Rectangle_Rectangle_Hole(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, zc, Edge;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    Edge = 100;
    if(zn==0)//2 Dimensions Grid
    {
        l=1;
        zc=0;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
        zc=zn+1;
    }

    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }else
                {
                    if(abs(i-xc)<=xn && abs(j-yc)<=yn && abs(k-zc)<=zn)
                    {
                        if(abs(i-xc)<xn && abs(j-yc)<yn && abs(k-zc)<zn)
                        {
                            BS->Boundary[i][j][k]=1;//1 is inner side.
                        }else if(abs(i-xc)==xn || abs(j-yc)==yn || abs(k-zc)==zn)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }

                    }else
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }
                }
            }
        }
    }

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            if(abs(i-xc)<Edge && abs(j-yc)<Edge)
            {
                BS->Boundary[i][j][0]=0;
            }
        }
    }

    for(i=xc-Edge-1;i<xc+Edge+1;++i)
    {
        for(j=yc-Edge-1;j<yc+Edge+1;++j)
        {
            if(BS->Boundary[i][j][0]!=0)
            {
                if(BS->Boundary[i-1][j][0]*BS->Boundary[i+1][j][0]*BS->Boundary[i][j-1][0]*BS->Boundary[i][j+1][0]==0)
                {
                    BS->Boundary[i][j][0]=2;
                }
            }
        }
    }


    return 1;

}


int Boundary_Polygon(Boundary_Shape *BS, int N_Polygon)
{
    int i, j, k, m, n, l, xn, yn, zn, xc, yc, rn, p;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    double **Line_Polygon, Line_Zone;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    //Obtain the minimum radius.
    if(xn>yn)
    {
        rn=yn;
    }else
    {
        rn=xn;
    }

    Line_Polygon = Make2DArray(N_Polygon,3);
    Regular_Polygon_Plane_Shape(rn,N_Polygon,Line_Polygon);
    xc=(m-1)/2, yc=(n-1)/2;//Origin position.

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(l==1)//2 Dimensions Case
                {
                    BS->Boundary[i][j][k]=1;//1 is inner side.
                    for(p=0;p<N_Polygon;++p)
                    {
                        Line_Zone = Line_Polygon[p][2]*(Line_Polygon[p][0]*(i-xc)+Line_Polygon[p][1]*(j-yc)+Line_Polygon[p][2]);
                        if(Line_Zone<0)
                        {
                            BS->Boundary[i][j][k]=0;//0 represents exterior side.
                        }
                    }
                }else
                {
                    BS->Boundary[i][j][k]=1;//1 is inner side.
                    for(p=0;p<N_Polygon;++p)
                    {
                        Line_Zone = Line_Polygon[p][2]*(Line_Polygon[p][0]*(i-xc)+Line_Polygon[p][1]*(j-yc)+Line_Polygon[p][2]);
                        if(Line_Zone<0)
                        {
                            BS->Boundary[i][j][k]=0;//0 represents exterior side.
                        }
                    }

                    if(k==0 || k==l-1)
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.

                    }
                }
            }
        }
    }

    Shape_Boundary_Assignment(BS);
    free2DArray(Line_Polygon,N_Polygon);

    return 1;

}

int Shape_Boundary_Assignment(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, Boundary_Jude;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    if(l==1)//2 Dimensions Case
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    if(i==0 || j==0 || i== m-1 || j==n-1)
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }

                    if(BS->Boundary[i][j][k]!=0)//Inner side or on Boundary
                    {
                        Boundary_Jude = BS->Boundary[i-1][j][k]*BS->Boundary[i+1][j][k]*BS->Boundary[i][j-1][k]*BS->Boundary[i][j+1][k];
                        if(Boundary_Jude==0)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }
                }
            }
        }
    }else//3 Dimensions Case
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    if(i==0 || j==0 || i==m-1 || j==n-1 || k==0 || k==l-1)
                    {
                        BS->Boundary[i][j][k]=0;//0 represents exterior side.
                    }

                    if(BS->Boundary[i][j][k]!=0)//Inner side or on Boundary
                    {
                        Boundary_Jude = BS->Boundary[i-1][j][k]*BS->Boundary[i+1][j][k]*BS->Boundary[i][j-1][k]*
                                        BS->Boundary[i][j+1][k]*BS->Boundary[i][j][k-1]*BS->Boundary[i][j][k+1];
                        if(Boundary_Jude==0)
                        {
                            BS->Boundary[i][j][k]=2;//2 is on the boundary.
                        }
                    }
                }
            }
        }
    }

    return 1;
}

int LocalEnergy_Initial(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;
    int p;

    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    if(l==1)//2 Dimensions Case
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                for(p=0;p<7;++p)
                {
                    LPP=(Local_Point_Postion) p;
                    ///////////////////////////////////////////////////////
                    switch(LPP)
                    {
                        ///////////////////////////////////////////////////
                        case Origin :
                        {
                            if(BS->Boundary[i][j][0]==0)
                            {
                                BS->LocalEnergy[i][j][0][p]=0;//0 represents no energy on this position.

                            }else
                            {
                                BS->LocalEnergy[i][j][0][p]=1;//1 represents that there is energy on this position.
                            }
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case North :
                        {
                            if(i==0)//Case North Boundary
                            {
                                BS->LocalEnergy[i][j][0][p]=0;
                            }else
                            {
                                if(BS->Boundary[i-1][j][0]==0)
                                {
                                    BS->LocalEnergy[i][j][0][p]=0;
                                }else
                                {
                                    BS->LocalEnergy[i][j][0][p]=1;
                                }
                            }
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case South :
                        {
                            if(i==m-1)//Case South Boundary
                            {
                                BS->LocalEnergy[i][j][0][p]=0;
                            }else
                            {
                                if(BS->Boundary[i+1][j][0]==0)
                                {
                                    BS->LocalEnergy[i][j][0][p]=0;
                                }else
                                {
                                    BS->LocalEnergy[i][j][0][p]=1;
                                }
                            }
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case West :
                        {
                            if(j==0)//Case West Boundary
                            {
                                BS->LocalEnergy[i][j][0][p]=0;
                            }else
                            {
                                if(BS->Boundary[i][j-1][0]==0)
                                {
                                    BS->LocalEnergy[i][j][0][p]=0;
                                }else
                                {
                                    BS->LocalEnergy[i][j][0][p]=1;
                                }
                            }
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case East :
                        {
                            if(j==n-1)//Case West Boundary
                            {
                                BS->LocalEnergy[i][j][0][p]=0;
                            }else
                            {
                                if(BS->Boundary[i][j+1][0]==0)
                                {
                                    BS->LocalEnergy[i][j][0][p]=0;
                                }else
                                {
                                    BS->LocalEnergy[i][j][0][p]=1;
                                }
                            }
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case Lower :
                        {
                            BS->LocalEnergy[i][j][0][p]=1;
                            break;
                        }
                        ///////////////////////////////////////////////////
                        case Upper :
                        {
                            BS->LocalEnergy[i][j][0][p]=1;
                            break;
                        }
                        ///////////////////////////////////////////////////
                        default :
                            break;
                    }
                }
            }
        }

    }else  //3 Dimensions Grid
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    for(p=0;p<7;++p)
                    {
                        LPP=(Local_Point_Postion) p;
                        ///////////////////////////////////////////////////////
                        switch(LPP)
                        {
                            ///////////////////////////////////////////////////
                            case Origin :
                            {
                                if(BS->Boundary[i][j][k]==0)
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;//0 represents no energy on this position.

                                }else
                                {
                                    BS->LocalEnergy[i][j][k][p]=1;//1 represents that there is energy on this position.
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case North :
                            {
                                if(i==0)//Case North Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i-1][j][k]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case South :
                            {
                                if(i==m-1)//Case South Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i+1][j][k]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case West :
                            {
                                if(j==0)//Case West Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i][j-1][k]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case East :
                            {
                                if(j==n-1)//Case East Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i][j+1][k]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case Lower :
                            {
                                if(k==0)//Case Lower Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i][j][k-1]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            case Upper :
                            {
                                if(k==l-1)//Case Upper Boundary
                                {
                                    BS->LocalEnergy[i][j][k][p]=0;
                                }else
                                {
                                    if(BS->Boundary[i][j][k+1]==0)
                                    {
                                        BS->LocalEnergy[i][j][k][p]=0;
                                    }else
                                    {
                                        BS->LocalEnergy[i][j][k][p]=1;
                                    }
                                }
                                break;
                            }
                            ///////////////////////////////////////////////////
                            default:
                                break;
                        }
                    }
                }
            }
        }
    }

    return 1;
}

int CalculPoints_Numbers(Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, p;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    int Boundary_Points=0, Inner_Points=0, Virtual_Points=0, Calculate_Points;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP;
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;
    //2 is on the boundary.
    //1 is inner side.
    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }
////////////////////////  Obtain the total points calculated  /////////////////////
    if(l==1)//2 Dimensions Case
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
                if(BS->Boundary[i][j][0]==1)//1 is inner side.
                {
                    Inner_Points=Inner_Points+1;
                }

                if(BS->Boundary[i][j][0]==2)//2 is on the boundary.
                {
                    Boundary_Points=Boundary_Points+1;
                }

                if(BS->Boundary[i][j][0]==2)//2 is on the boundary.
                {
                    for(p=1;p<5;++p)
                    {
                        LPP=(Local_Point_Postion) p;
                        if(BS->LocalEnergy[i][j][0][p]==0)
                        {
                            switch(LPP)
                            {
                                case North:
                                {
                                    if(BS->Boundary_Type_x==0)
                                    {
                                        Virtual_Points=Virtual_Points+1;
                                    }
                                    break;
                                }

                                case South:
                                {
                                    if(BS->Boundary_Type_x==0)
                                    {
                                        Virtual_Points=Virtual_Points+1;
                                    }
                                    break;
                                }

                                case West:
                                {
                                    if(BS->Boundary_Type_y==0)
                                    {
                                        Virtual_Points=Virtual_Points+1;
                                    }
                                    break;
                                }

                                case East:
                                {
                                    if(BS->Boundary_Type_y==0)
                                    {
                                        Virtual_Points=Virtual_Points+1;
                                    }
                                    break;
                                }

                                default:
                                    break;
                            }
                        }
                    }
                }
            }
        }
    }else//3 Dimensions Case
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
                    if(BS->Boundary[i][j][k]==1)//1 is inner side.
                    {
                        Inner_Points=Inner_Points+1;
                    }

                    if(BS->Boundary[i][j][k]==2)//2 is on the boundary.
                    {
                        Boundary_Points=Boundary_Points+1;
                    }

                    if(BS->Boundary[i][j][k]==2)//2 is on the boundary.
                    {
                        for(p=1;p<7;++p)
                        {
                            LPP=(Local_Point_Postion) p;
                            if(BS->LocalEnergy[i][j][k][p]==0)
                            {
                                switch(LPP)
                                {
                                    case North:
                                    {
                                        if(BS->Boundary_Type_x==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    case South:
                                    {
                                        if(BS->Boundary_Type_x==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    case West:
                                    {
                                        if(BS->Boundary_Type_y==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    case East:
                                    {
                                        if(BS->Boundary_Type_y==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    case Lower:
                                    {
                                        if(BS->Boundary_Type_z==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    case Upper:
                                    {
                                        if(BS->Boundary_Type_z==0)
                                        {
                                            Virtual_Points=Virtual_Points+1;
                                        }
                                        break;
                                    }

                                    default:
                                        break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    Virtual_Points=Virtual_Points+1;
    Calculate_Points=Inner_Points+Boundary_Points+Virtual_Points;

    BS->Inner_Points     = Inner_Points;
    BS->Boundary_Points  = Boundary_Points;
    BS->Virtual_Points   = Virtual_Points;
    BS->Calculate_Points = Calculate_Points;

    //printf("%d\n",Calculate_Points);
    return 1;
}

int CalculPointsParameters_Initial(Calculate_Points_Parameters *CPP, Boundary_Shape *BS)
{
    int i, j, k, m, n, l, xn, yn, zn, p, q, N, N1, r, s, v, PointNumbers, VirtualPointsPositions;
    //x-direction: i, m, xn; y-direction: j, n, yn; z-direction: k, l, zn;
    int Calculate_Points;
    int ***Local_Points_1D, ***Boundary, ****LocalEnergy, ***Position_3DTo1D, *Virtual_Position;
    int ****Local_Points_1D_LocalEnergy;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP, LPP1;
    PointNumbers=1;//0 is the common one.
    VirtualPointsPositions=1;//0 is the common one.
    //////////////////Input Parameter////////////////////////
    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }

    Boundary         = BS->Boundary;
    LocalEnergy      = BS->LocalEnergy;
    Calculate_Points = BS->Calculate_Points;
    Position_3DTo1D  = CPP->Position_3DTo1D;
    Virtual_Position = CPP->Virtual_Position;
    Local_Points_1D  = CPP->Local_Points_1D;
    Local_Points_1D_LocalEnergy=CPP->Local_Points_1D_LocalEnergy;

    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                Position_3DTo1D[i][j][k]=0;
            }
        }
    }

    //////////////////Initial Part////////////////////////////
    if(l==1)//2 Dimensions Case
    {
        for(i=0;i<m;++i)
        {
            for(j=0;j<n;++j)
            {
///////////////////////////////////////////   Inner side   /////////////////////////////////////////////////////////////
                if(Boundary[i][j][0]==1)//1 is inner side.
                {
                    Position_3DTo1D[i][j][0]=PointNumbers;
                    for(p=0;p<7;++p)
                    {
                        LPP=(Local_Point_Postion) p;
                        switch(LPP)
                        {
                            case Origin:
                            {
                                Local_Points_1D[PointNumbers][p][0]=1;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case North:
                            {
                                Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][0];//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i-1;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case South:
                            {
                                Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][0];//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i+1;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case West:
                            {
                                Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][0];//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j-1;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case East:
                            {
                                Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][0];//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j+1;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case Lower:
                            {
                                Local_Points_1D[PointNumbers][p][0]=1;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case Upper:
                            {
                                Local_Points_1D[PointNumbers][p][0]=1;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            default:
                                break;
                        }
                    }
                    PointNumbers=PointNumbers+1;
                }
///////////////////////////////////////////   Boundary points   /////////////////////////////////////////////////////////////
                if(Boundary[i][j][0]==2)//2 is on the boundary.
                {
                    Position_3DTo1D[i][j][0]=PointNumbers;
                    for(p=0;p<7;++p)
                    {
                        LPP=(Local_Point_Postion) p;
                        switch(LPP)
                        {
                            case Origin:
                            {
                                Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case North:
                            {
                                if(LocalEnergy[i][j][0][p]==1)
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][0];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i-1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                }else
                                {
                                    if(BS->Boundary_Type_x==0)//Free Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][0];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=0;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=0;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                        VirtualPointsPositions=VirtualPointsPositions+1;
                                    }else//Periodic Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=m-1-i;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    }
                                }
                                break;
                            }

                            case South:
                            {
                                if(LocalEnergy[i][j][0][p]==1)
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][0];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i+1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                }else
                                {
                                    if(BS->Boundary_Type_x==0)//Free Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][0];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=0;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=0;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                        VirtualPointsPositions=VirtualPointsPositions+1;
                                    }else//Periodic Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=m-1-i;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    }
                                }
                                break;
                            }

                            case West:
                            {
                                if(LocalEnergy[i][j][0][p]==1)
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][0];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j-1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                }else
                                {
                                    if(BS->Boundary_Type_y==0)//Free Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][0];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                        VirtualPointsPositions=VirtualPointsPositions+1;
                                    }else//Periodic Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=n-1-j;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    }
                                }
                                break;
                            }

                            case East:
                            {
                                if(LocalEnergy[i][j][0][p]==1)
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][0];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j+1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                }else
                                {
                                    if(BS->Boundary_Type_y==0)//Free Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][0];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][0][p];//5 is the energy proportion.
                                        VirtualPointsPositions=VirtualPointsPositions+1;
                                    }else//Periodic Boundary
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=n-1-j;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    }
                                }
                                break;
                            }

                            case Lower:
                            {
                                Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            case Upper:
                            {
                                Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                break;
                            }

                            default:
                                break;
                        }
                    }
                    PointNumbers=PointNumbers+1;
                }
///////////////////////////////////////////   Virtual points   /////////////////////////////////////////////////////////////
                if(Boundary[i][j][0]==2)//2 is on the boundary.
                {
                    N=Position_3DTo1D[i][j][0];
                    for(p=0;p<7;++p)
                    {
                        LPP=(Local_Point_Postion) p;
                        switch(LPP)
                        {
                            case North:
                            {
                                if(LocalEnergy[i][j][0][p]==0 && BS->Boundary_Type_x==0)//Free Boundary
                                {
                                    for(q=0;q<7;++q)
                                    {
                                        LPP1=(Local_Point_Postion) q;
                                        switch(LPP1)
                                        {
                                            case Origin:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case South:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                break;
                                            }

                                            case Lower:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case Upper:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            default:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }
                                        }
                                    }
                                    Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                    PointNumbers=PointNumbers+1;
                                }
                                break;
                            }

                            case South:
                            {
                                if(LocalEnergy[i][j][0][p]==0 && BS->Boundary_Type_x==0)//Free Boundary
                                {
                                    for(q=0;q<7;++q)
                                    {
                                        LPP1=(Local_Point_Postion) q;
                                        switch(LPP1)
                                        {
                                            case Origin:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case North:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                break;
                                            }

                                            case Lower:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case Upper:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            default:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }
                                        }
                                    }
                                    Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                    PointNumbers=PointNumbers+1;
                                }
                                break;
                            }

                            case West:
                            {
                                if(LocalEnergy[i][j][0][p]==0 && BS->Boundary_Type_y==0)//Free Boundary
                                {
                                    for(q=0;q<7;++q)
                                    {
                                        LPP1=(Local_Point_Postion) q;
                                        switch(LPP1)
                                        {
                                            case Origin:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case East:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                break;
                                            }

                                            case Lower:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case Upper:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            default:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }
                                        }
                                    }
                                    Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                    PointNumbers=PointNumbers+1;
                                }
                                break;
                            }

                            case East:
                            {
                                if(LocalEnergy[i][j][0][p]==0 && BS->Boundary_Type_y==0)//Free Boundary
                                {
                                    for(q=0;q<7;++q)
                                    {
                                        LPP1=(Local_Point_Postion) q;
                                        switch(LPP1)
                                        {
                                            case Origin:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case West:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                break;
                                            }

                                            case Lower:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            case Upper:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }

                                            default:
                                            {
                                                Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                break;
                                            }
                                        }
                                    }
                                    Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                    PointNumbers=PointNumbers+1;
                                }
                                break;
                            }

                            default:
                                break;
                        }
                    }
                }
            }
        }
    }else///////////////////////////////////////////// 3D case //////////////////////////////////////////////////////////////////////////
    {
        for(k=0;k<l;++k)
        {
            for(i=0;i<m;++i)
            {
                for(j=0;j<n;++j)
                {
///////////////////////////////////////////   Inner side   /////////////////////////////////////////////////////////////
                    if(Boundary[i][j][k]==1)//1 is inner side.
                    {
                        Position_3DTo1D[i][j][k]=PointNumbers;
                        for(p=0;p<7;++p)
                        {
                            LPP=(Local_Point_Postion) p;
                            switch(LPP)
                            {
                                case Origin:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=1;//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case North:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][k];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i-1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case South:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][k];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i+1;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case West:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][k];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j-1;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case East:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][k];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j+1;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case Lower:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k-1];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k-1;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=0;//5 is the energy proportion.
                                    break;
                                }

                                case Upper:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k+1];//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k+1;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=0;//5 is the energy proportion.
                                    break;
                                }

                                default:
                                    break;
                            }
                        }
                        PointNumbers=PointNumbers+1;
                    }
///////////////////////////////////////////   Boundary points   /////////////////////////////////////////////////////////////
                    if(Boundary[i][j][k]==2)//2 is on the boundary.
                    {
                        Position_3DTo1D[i][j][k]=PointNumbers;
                        for(p=0;p<7;++p)
                        {
                            LPP=(Local_Point_Postion) p;
                            switch(LPP)
                            {
                                case Origin:
                                {
                                    Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                    Local_Points_1D[PointNumbers][p][1]=i;//1 is the x-coordinate.
                                    Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                    Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                    Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                    Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                    break;
                                }

                                case North:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][k];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i-1;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_x==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i-1][j][k];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else//Periodic Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=m-1-i;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                case South:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][k];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i+1;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_x==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i+1][j][k];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else//Periodic Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=m-1-i;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=j;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                case West:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][k];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j-1;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_y==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i][j-1][k];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else//Periodic Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=n-1-j;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                case East:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][k];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j+1;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_y==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i][j+1][k];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else//Periodic Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=n-1-j;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                case Lower:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k-1];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k-1;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_z==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k-1];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else//Periodic Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=j;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=l-1-k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                case Upper:
                                {
                                    if(LocalEnergy[i][j][k][p]==1)
                                    {
                                        Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k+1];//0 is the boundary type.
                                        Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                        Local_Points_1D[PointNumbers][p][2]=j;//1 is the x-coordinate.
                                        Local_Points_1D[PointNumbers][p][3]=k+1;//3 is the z-coordinate.
                                        Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                        Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                    }else
                                    {
                                        if(BS->Boundary_Type_z==0)//Free Boundary
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=Boundary[i][j][k+1];//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=0;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=0;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=0;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=VirtualPointsPositions;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=LocalEnergy[i][j][k][p];//5 is the energy proportion.
                                            VirtualPointsPositions=VirtualPointsPositions+1;
                                        }else
                                        {
                                            Local_Points_1D[PointNumbers][p][0]=2;//0 is the boundary type.
                                            Local_Points_1D[PointNumbers][p][1]=i;//2 is the y-coordinate.
                                            Local_Points_1D[PointNumbers][p][2]=j;//1 is the x-coordinate.
                                            Local_Points_1D[PointNumbers][p][3]=l-1-k;//3 is the z-coordinate.
                                            Local_Points_1D[PointNumbers][p][4]=0;//4 is the virtual point position.
                                            Local_Points_1D[PointNumbers][p][5]=1;//5 is the energy proportion.
                                        }
                                    }
                                    break;
                                }

                                default:
                                    break;
                            }
                        }
                        PointNumbers=PointNumbers+1;
                    }
///////////////////////////////////////////   Virtual points   /////////////////////////////////////////////////////////////
                    if(Boundary[i][j][k]==2)//2 is on the boundary.
                    {
                        N=Position_3DTo1D[i][j][k];
                        for(p=0;p<7;++p)
                        {
                            LPP=(Local_Point_Postion) p;
                            switch(LPP)
                            {
                                case North:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_x==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case South:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                case South:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_x==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case North:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                case West:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_y==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case East:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                case East:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_y==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case West:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                case Lower:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_z==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case Upper:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                case Upper:
                                {
                                    if(LocalEnergy[i][j][k][p]==0 && BS->Boundary_Type_z==0)//Free Boundary
                                    {
                                        for(q=0;q<7;++q)
                                        {
                                            LPP1=(Local_Point_Postion) q;
                                            switch(LPP1)
                                            {
                                                case Origin:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=Local_Points_1D[N][p][4];//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }

                                                case Lower:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=2;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=i;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=j;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=k;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=1;//5 is the energy proportion.
                                                    break;
                                                }

                                                default:
                                                {
                                                    Local_Points_1D[PointNumbers][q][0]=0;//0 is the boundary type.
                                                    Local_Points_1D[PointNumbers][q][1]=0;//1 is the x-coordinate.
                                                    Local_Points_1D[PointNumbers][q][2]=0;//2 is the y-coordinate.
                                                    Local_Points_1D[PointNumbers][q][3]=0;//3 is the z-coordinate.
                                                    Local_Points_1D[PointNumbers][q][4]=0;//4 is the virtual point position.
                                                    Local_Points_1D[PointNumbers][q][5]=0;//5 is the energy proportion.
                                                    break;
                                                }
                                            }
                                        }
                                        Virtual_Position[Local_Points_1D[N][p][4]]=PointNumbers;
                                        PointNumbers=PointNumbers+1;
                                    }
                                    break;
                                }

                                default:
                                    break;
                            }
                        }
                    }
                }
            }
        }
    }

    Virtual_Position[0]=0;
    for(p=0;p<7;++p)
    {
        for(s=0;s<7;++s)
        {
            Local_Points_1D[0][p][s]=0;
        }
    }

///////////////////////////////////////// Local Energy ///////////////////////////////////
    for(N=1;N<Calculate_Points;++N)
    {
        for(r=0;r<7;++r)
        {
            for(s=0;s<7;++s)
            {
                for(p=0;p<7;++p)
                {
                    if(Local_Points_1D[N][p][0]==0)
                    {
                        v = Local_Points_1D[N][p][4];
                        N1= Virtual_Position[v];
                    }else
                    {
                        i = Local_Points_1D[N][p][1];
                        j = Local_Points_1D[N][p][2];
                        k = Local_Points_1D[N][p][3];
                        N1= Position_3DTo1D[i][j][k];
                    }
                    Local_Points_1D[N][p][6]=N1;//Neighbor points in 1D.
                }
            }
        }
    }

    for(N=1;N<Calculate_Points;++N)
    {
        for(r=0;r<7;++r)
        {
            for(s=0;s<7;++s)
            {
                for(p=0;p<7;++p)
                {
                    N1=Local_Points_1D[N][p][6];
                    Local_Points_1D_LocalEnergy[N][p][r][s]=Local_Points_1D[N1][r][s];
                }
            }
        }
    }

    for(r=0;r<7;++r)
    {
        for(s=0;s<7;++s)
        {
            for(p=0;p<7;++p)
            {
                Local_Points_1D_LocalEnergy[0][p][r][s]=0;
            }
        }
    }

    //printf("%d\t%d\t%d\t%d\n",PointNumbers-1,Calculate_Points,VirtualPointsPositions,Virtual_Points);
/*
    for(N=1;N<Calculate_Points;++N)
    {
        printf("%d\n",N);
        for(p=0;p<7;++p)
        {
            for(r=0;r<7;++r)
            {
                printf("%d\t",Local_Points_1D[N][p][r]);
            }
            printf("\n");
        }
        printf("\n");
    }
*/
    return 1;
}

int Neighbour_Points_Local_Initial(Calculate_Points_Parameters *CPP)
{
    int ***Neighbour_Points_Local_Coordinate, p, q, i, j, k;
    Neighbour_Points_Local_Coordinate=CPP->Neighbour_Points_Local_Coordinate;

    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    Local_Point_Postion LPP, LPP1;

    i=2;    j=2;    k=2;

    for(p=0; p<7; ++p)
    {
        LPP=(Local_Point_Postion) p;
        switch(LPP)
        {
            case Origin:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i;//0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j;//1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k;//2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case North:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i-1;//0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j;//1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k;//2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case South:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i+1;//0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j;  //1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k;  //2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case West:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i;  //0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j-1;//1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k;  //2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case East:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i;  //0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j+1;//1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k;  //2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case Lower:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i;  //0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j;  //1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k-1;//2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            case Upper:
            {
                Neighbour_Points_Local_Coordinate[0][p][0]=i;  //0 is the x-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][1]=j;  //1 is the y-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][2]=k+1;//2 is the z-coordinate.
                Neighbour_Points_Local_Coordinate[0][p][3]=1;//Repetition Index, 1 is not repeated.
                break;
            }

            default:
                break;
        }
    }

    for(p=0; p<7; ++p)
    {
        LPP=(Local_Point_Postion) p;
        for(q=0; q<7; ++q)
        {
            LPP1=(Local_Point_Postion) q;
            switch(LPP)
            {
                case Origin:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0];//0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1];//1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2];//2 is the z-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                    break;
                }

                case North:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0]-1;//0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1];  //1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2];  //2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case North:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                case South:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0]+1;//0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1];  //1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2];  //2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case South:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                case West:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0];  //0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1]-1;//1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2];  //2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case West:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                case East:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0];  //0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1]+1;//1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2];  //2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case East:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                case Lower:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0];  //0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1];  //1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2]-1;//2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case Lower:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                case Upper:
                {
                    Neighbour_Points_Local_Coordinate[p][q][0]=Neighbour_Points_Local_Coordinate[0][q][0];  //0 is the x-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][1]=Neighbour_Points_Local_Coordinate[0][q][1];  //1 is the y-coordinate.
                    Neighbour_Points_Local_Coordinate[p][q][2]=Neighbour_Points_Local_Coordinate[0][q][2]+1;//2 is the z-coordinate.
                    switch(LPP1)
                    {
                        case Upper:
                        Neighbour_Points_Local_Coordinate[p][q][3]=1;//Repetition Index, 1 is not repeated.
                        break;

                        default:
                        Neighbour_Points_Local_Coordinate[p][q][3]=0;//Repetition Index, 1 is not repeated.
                        break;
                    }
                    break;
                }

                default:
                    break;
            }
        }
    }

    return 1;
}

int Neighbour_Points_3Dto1D(Calculate_Points_Parameters *CPP)
{
    int **Neighbour_Points_3Dto1D_Array, p, q;
    Neighbour_Points_3Dto1D_Array=CPP->Neighbour_Points_3Dto1D_Array;
    enum Local_Point_Postion{Origin, North, South, West, East, Lower, Upper};
    enum Local_Point_Postion_L{Null, O, N, S, W, E, L, U, NN, SS, WW, EE, LL, UU};
    //Origin(i,j,k), North(i-1,j,k), South(i+1,j,k), West(i,j-1,k), East(i,j+1,k), Lower(i,j,k-1), Upper(i,j,k+1).
    //NN(i-2,j,k),SS(i+2,j,k), WW(i,j-2,k), EE(i,j+2,k), LL(i,j,k-2), UU(i,j,k+2);
    Local_Point_Postion LPP;
    Local_Point_Postion_L LPPL;

    for(p=0;p<14;++p)
    {
        LPPL=(Local_Point_Postion_L) p;
        for(q=0;q<7;++q)
        {
            LPP=(Local_Point_Postion) q;
            switch(LPPL)
            {
///////////////////////////////////////////////////////////////////////////
                case O :
                {
                    switch(LPP)
                    {
                        case Origin:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case N :
                {
                    switch(LPP)
                    {
                        case North:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case S :
                {
                    switch(LPP)
                    {
                        case South:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case W :
                {
                    switch(LPP)
                    {
                        case West:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case E :
                {
                    switch(LPP)
                    {
                        case East:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case L :
                {
                    switch(LPP)
                    {
                        case Lower:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }

///////////////////////////////////////////////////////////////////////////
                case U :
                {
                    switch(LPP)
                    {
                        case Upper:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=0;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case NN :
                {
                    switch(LPP)
                    {
                        case North:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case SS :
                {
                    switch(LPP)
                    {
                        case South:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case WW :
                {
                    switch(LPP)
                    {
                        case West:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case EE :
                {
                    switch(LPP)
                    {
                        case East:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case LL :
                {
                    switch(LPP)
                    {
                        case Lower:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                case UU :
                {
                    switch(LPP)
                    {
                        case Upper:
                        {
                            Neighbour_Points_3Dto1D_Array[p][0]=q;
                            Neighbour_Points_3Dto1D_Array[p][1]=q;
                            break;
                        }

                        default:
                            break;
                    }
                    break;
                }
///////////////////////////////////////////////////////////////////////////
                default:
                {
                    Neighbour_Points_3Dto1D_Array[p][0]=0;
                    Neighbour_Points_3Dto1D_Array[p][1]=0;
                    break;
                }
            }
        }
    }

    for(p=0;p<14;++p)
    {
        //printf("%d\t%d\n",Neighbour_Points_3Dto1D_Array[p][0],Neighbour_Points_3Dto1D_Array[p][1]);
    }

    return 1;
}

int ReSort_Calcualte_Points(Calculate_Points_Parameters *CPP,Boundary_Shape *BS)
{
    int i=0, p;
    for(p=0;p<BS->Calculate_Points;++p)
    {
        if(CPP->Local_Points_1D[p][0][0]==0)
        {
            CPP->Sort_Calculate_Points[i]=p;
            i=i+1;
            //printf("%d\t%d\n",i,p);
        }
    }

    for(p=0;p<BS->Calculate_Points;++p)
    {
        if(CPP->Local_Points_1D[p][0][0]!=0)
        {
            CPP->Sort_Calculate_Points[i]=p;
            i=i+1;
            //printf("%d\t%d\n",i,p);
        }
    }
    return 1;
}

int ReSort_Calcualte_Points_ByThreads(Calculate_Points_Parameters *CPP,Boundary_Shape *BS, int Threads_Number)
{
    int i=0, p, q, CalculatePointsThreads_Lower, CalculatePointsThreads_Upper;
    for(p=0;p<Threads_Number;++p)
    {
        CalculatePointsThreads_Lower=floor(BS->Calculate_Points/Threads_Number)*p;
        if(p==Threads_Number-1)
        {
            CalculatePointsThreads_Upper = BS->Calculate_Points;
        }else
        {
            CalculatePointsThreads_Upper = floor(BS->Calculate_Points/Threads_Number)*(p+1);
        }

        for(q=CalculatePointsThreads_Lower;q<CalculatePointsThreads_Upper;++q)
        {
            if(CPP->Local_Points_1D[q][0][0]==0)
            {
                CPP->Sort_Calculate_Points[i]=q;
                i=i+1;
                //printf("%d\t%d\n",i,q);
            }
        }

        for(q=CalculatePointsThreads_Lower;q<CalculatePointsThreads_Upper;++q)
        {
            if(CPP->Local_Points_1D[q][0][0]!=0)
            {
                CPP->Sort_Calculate_Points[i]=q;
                i=i+1;
                //printf("%d\t%d\n",i,q);
            }
        }
    }

    return 1;
}

int ReSort_Calcualte_Points_ByOrders(Calculate_Points_Parameters *CPP,Boundary_Shape *BS)
{
    int p;
/*
    for(p=1;p<BS->Calculate_Points;++p)
    {
        CPP->Sort_Calculate_Points[p]=BS->Calculate_Points-p;
    }
    CPP->Sort_Calculate_Points[0]=0;
*/
    for(p=0;p<BS->Calculate_Points;++p)
    {
        CPP->Sort_Calculate_Points[p]=p;
    }
    return 1;
}

int ReSort_Calcualte_Points_BySymmetry_Xaxis(Calculate_Points_Parameters *CPP,Boundary_Shape *BS)
{
    int p, i, j, j_sym, k, v, xn, yn, zn, m, n, l, N, q, q_sym;
    int *Position;
    Position=Make1DArrayinteger(7);
    Position[0]=0;  Position[1]=1;
    Position[2]=2;  Position[3]=4;
    Position[4]=3;  Position[5]=5;
    Position[6]=6;

    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }
    p = 1;
    for(i=0;i<m;++i)
    {
        for(j=0;j<n;++j)
        {
            j_sym = n-j-1;
            for(k=0;k<l;++k)
            {
                if(BS->Boundary[i][j_sym][k]!=0)
                {
                    N = CPP->Position_3DTo1D[i][j_sym][k];
                    CPP->Sort_Calculate_Points[p] = N;
                    p = p+1;
                    if(CPP->Local_Points_1D[N][0][0]==2)
                    {
                        for(q=0;q<7;++q)
                        {
                            q_sym = Position[q];
                            if(CPP->Local_Points_1D[N][q_sym][5]==0)
                            {
                                v = CPP->Local_Points_1D[N][q_sym][4];
                                if(v!=0)
                                {
                                    CPP->Sort_Calculate_Points[p]=CPP->Local_Points_1D[N][q_sym][6];
                                    p = p+1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    CPP->Sort_Calculate_Points[0]=0;
    free(Position);
    //printf("%d\t%d\n",p-1,BS->Calculate_Points);
    return 1;
}

int ReSort_Calcualte_Points_BySymmetry_Yaxis(Calculate_Points_Parameters *CPP,Boundary_Shape *BS)
{
    int p, i, j, i_sym, k, v, xn, yn, zn, m, n, l, N, q, q_sym;
    int *Position;
    Position=Make1DArrayinteger(7);
    Position[0]=0;  Position[1]=2;
    Position[2]=1;  Position[3]=3;
    Position[4]=4;  Position[5]=5;
    Position[6]=6;

    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }
    p = 1;
    for(i=0;i<m;++i)
    {
        i_sym = m-i-1;
        for(j=0;j<n;++j)
        {
            for(k=0;k<l;++k)
            {
                if(BS->Boundary[i_sym][j][k]!=0)
                {
                    N = CPP->Position_3DTo1D[i_sym][j][k];
                    CPP->Sort_Calculate_Points[p] = N;
                    p = p+1;
                    if(CPP->Local_Points_1D[N][0][0]==2)
                    {
                        for(q=0;q<7;++q)
                        {
                            q_sym = Position[q];
                            if(CPP->Local_Points_1D[N][q_sym][5]==0)
                            {
                                v = CPP->Local_Points_1D[N][q_sym][4];
                                if(v!=0)
                                {
                                    CPP->Sort_Calculate_Points[p]=CPP->Local_Points_1D[N][q_sym][6];
                                    p = p+1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    CPP->Sort_Calculate_Points[0]=0;
    free(Position);
    //printf("%d\t%d\n",p-1,BS->Calculate_Points);
    return 1;
}

int ReSort_Calcualte_Points_ByRow(Calculate_Points_Parameters *CPP,Boundary_Shape *BS)
{
    int p, i, j, k, v, xn, yn, zn, m, n, l, N, q, q_row;
    int *Position;
    Position=Make1DArrayinteger(7);
    Position[0]=0;  Position[1]=3;
    Position[2]=4;  Position[3]=1;
    Position[4]=2;  Position[5]=5;
    Position[6]=6;

    xn=BS->xn;      yn=BS->yn;
    zn=BS->zn;      m =xn*2+3;
    n =yn*2+3;

    if(zn==0)//2 Dimensions Grid
    {
        l=1;
    }else    //3 Dimensions Grid
    {
        l=zn*2+3;
    }
    p = 1;

    for(j=0;j<n;++j)
    {
        for(i=0;i<m;++i)
        {
            for(k=0;k<l;++k)
            {
                if(BS->Boundary[i][j][k]!=0)
                {
                    N = CPP->Position_3DTo1D[i][j][k];
                    CPP->Sort_Calculate_Points[p] = N;
                    p = p+1;
                    if(CPP->Local_Points_1D[N][0][0]==2)
                    {
                        for(q=0;q<7;++q)
                        {
                            q_row = Position[q];
                            if(CPP->Local_Points_1D[N][q_row][5]==0)
                            {
                                v = CPP->Local_Points_1D[N][q_row][4];
                                if(v!=0)
                                {
                                    CPP->Sort_Calculate_Points[p]=CPP->Local_Points_1D[N][q_row][6];
                                    p = p+1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    CPP->Sort_Calculate_Points[0]=0;
    free(Position);
    //printf("%d\t%d\n",p-1,BS->Calculate_Points);
    return 1;
}

int Boundary_Initial_Type(Input_Parameter_NMD *IPN, Boundary_Shape *BS)
{
    int xn, yn, rn;

    xn=BS->xn;  yn=BS->yn;
    if(xn>yn)
    {
        rn=yn;
    }else
    {
        rn=xn;
    }
    BS->Boundary_Type_x = IPN->Boundary_Type_x;
    BS->Boundary_Type_y = IPN->Boundary_Type_y;
    BS->Boundary_Type_z = IPN->Boundary_Type_z;

    switch (IPN->Boundary_Shape_Type)
    {
        case BITN_Rectangle:
        {
            Boundary_Rectangle(BS);
            break;
        }

        case BITN_Circle:
        {
            Boundary_Circle(rn,BS);
            break;
        }

        case BITN_Circle_Hole:
        {
            Boundary_Circle_Hole(BS);
            break;
        }

        case BITN_Ellipse:
        {
            Boundary_Ellipse(BS);
            break;
        }

        case BITN_Triangle:
        {
            Boundary_Polygon(BS,3);
            break;
        }

        default:
        break;
    }
/*
    if(IPN->Boundary_Shape_Type==0)
    {
        Boundary_Rectangle(BS);
    }else if(IPN->Boundary_Shape_Type==1)
    {
        Boundary_Circle(rn,BS);
    }else if(IPN->Boundary_Shape_Type==-1)
    {
        Boundary_Circle_Hole(BS);
    }else if(IPN->Boundary_Shape_Type==-2)
    {
        Boundary_Ellipse(BS);
    }else if(IPN->Boundary_Shape_Type==2)
    {
        Boundary_Rectangle_20Grid_Crack(BS);
        //Boundary_Rectangle_Circle(BS);
        //Boundary_Rectangle_4Circle_Holes(BS);
        //Boundary_Rectangle_Rectangle_Hole(BS);
        //Boundary_Rectangle_Ellipse(BS);
    }else if(IPN->Boundary_Shape_Type==3)
    {
        Boundary_Polygon(BS,3);
    }else if(IPN->Boundary_Shape_Type==4)
    {
        Boundary_Polygon(BS,4);
    }else if(IPN->Boundary_Shape_Type==5)
    {
        Boundary_Polygon(BS,5);
    }else if(IPN->Boundary_Shape_Type==6)
    {
        Boundary_Polygon(BS,6);
    }else if(IPN->Boundary_Shape_Type==7)
    {
        Boundary_Polygon(BS,7);
    }else if(IPN->Boundary_Shape_Type==8)
    {
        Boundary_Polygon(BS,8);
    }else if(IPN->Boundary_Shape_Type==9)
    {
        Boundary_Polygon(BS,9);
    }else if(IPN->Boundary_Shape_Type==10)
    {
        Boundary_Polygon(BS,10);
    }else if(IPN->Boundary_Shape_Type==11)
    {
        Boundary_Polygon(BS,11);
    }else if(IPN->Boundary_Shape_Type==12)
    {
        Boundary_Polygon(BS,12);
    }else
    {
        Boundary_Polygon(BS,IPN->Boundary_Shape_Type);
    }
*/
    return 1;
}
