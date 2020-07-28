// cpp file


/***********************************************

      FileName: wavefront.cpp

        Author: stj
   Description: ---
 First  Create: 2019-07-10 15:33:56
 Last Modified: 2020-07-23 11:36:40

***********************************************/

#include"wavefront.h"
#include<iostream>
#include<cmath>
#include<cstdlib>
using namespace std;

void wavefront::getv(const Array2D<float> velocity)
{
    int a=floor(x);
    int b=floor(z);
   //cout<<a<<"  "<<b<<endl;
    float dx=(x-a)/k;
    float dz=(z-b)/h;
	if((x<2) && (z<2))
    {
        v=velocity[1][1];
    }
    else if((x>nx-3) && (z>nz-3))
    {
        v=velocity[nx-2][nz-2];
    }
    else if(x<2)
    {
        v=velocity[1][b];

    }
    else if(z<2)
    {
        
        v=velocity[a][1];
    }
    else if(x>nx-3)
    {
       // cout<<"yes"<<endl;
        v=velocity[nx-2][b];
    }
    else if(z>nz-3)
    {
        v=velocity[a][nz-2];
    }// border ,fuck!!!
    else
    {
    Array1D<float> X;
    Array1D<float> Z;
    X.Resize(4);
    Z.Resize(4);
    X.Zero();
    Z.Zero();
   
    X[0]= -1*dz*dz*dz+2*dz*dz-dz;
    X[1]=  3*dz*dz*dz-5*dz*dz+2;
    X[2]= -3*dz*dz*dz+4*dz*dz+dz;
    X[3]= dz*dz*dz-dz*dz;
   
    Z[0]= -dx*dx*dx+2*dx*dx-dx;
    Z[1]= 3*dx*dx*dx-5*dx*dx+2;
    Z[2]= -3*dx*dx*dx+4*dx*dx+dx;
    Z[3]= dx*dx*dx-dx*dx;

    Array2D<float> H;
    H.Resize(4,4);
    H.Zero();
    H[0][0]=velocity[a-1][b-1]*X[0]; 
    H[0][1]=velocity[a][b-1]*X[0];
    H[0][2]=velocity[a+1][b-1]*X[0];
    H[0][3]=velocity[a+2][b-1]*X[0];
    H[1][0]=velocity[a-1][b]*X[1];   H[1][1]=velocity[a][b]*X[1];   H[1][2]=velocity[a+1][b]*X[1];   H[1][3]=velocity[a+2][b]*X[1];
    H[2][0]=velocity[a-1][b+1]*X[2]; H[2][1]=velocity[a][b+1]*X[2]; H[2][2]=velocity[a+1][b+1]*X[2]; H[2][3]=velocity[a+2][b+1]*X[2];
    H[3][0]=velocity[a-1][b+2]*X[3]; H[3][1]=velocity[a][b+2]*X[3]; H[3][2]=velocity[a+1][b+2]*X[3]; H[3][3]=velocity[a+2][b+2]*X[3];


    v=   (H[0][0]+H[1][0]+H[2][0]+H[3][0])*Z[0]*0.25
        +(H[0][1]+H[1][1]+H[2][1]+H[3][1])*Z[1]*0.25
        +(H[0][2]+H[1][2]+H[2][2]+H[3][2])*Z[2]*0.25
        +(H[0][3]+H[1][3]+H[2][3]+H[3][3])*Z[3]*0.25;
//cout << (H[0][0] + H[1][0] + H[2][0] + H[3][0])*Z[0] * 0.25 << endl;
//	cout << (H[0][1] + H[1][1] + H[2][1] + H[3][1])*Z[1] * 0.25 << endl;
//	cout << (H[0][2] + H[1][2] + H[2][2] + H[3][2])*Z[2] * 0.25 << endl;
//	cout << (H[0][3] + H[1][3] + H[2][3] + H[3][3])*Z[3] * 0.25 << endl;
//	cout << v << endl;
    }
}



void wavefront::getvx(Array2D<float> velocityX)
{

    int a=floor(x);
    int b=floor(z);
    float dx=(x-a)/k;
    float dz=(z-b)/h;
    //cout<<a<<"  "<<b<<endl;
	/*if (a >= nx - 1)
	{
		a = a - 2;
	}
	if (b >= nz - 1)
	{
		b = b - 2;
	}
    if((x<2)||(z<2))
    {
        float t=0;
        float u=0;
        t=(a+1-x)*velocityX[a][b]+(x-a)*velocityX[a][b+1];
        u=(a+1-x)*velocityX[a+1][b]+(x-a)*velocityX[a+1][b+1];
        vx=(b+1-z)*t+(z-b)*u;
    

    }
    else if((x+1>nx-1)||(z+1>nz-1))
    {
        float t=0;
        float u=0;
        t=(a+1-x)*velocityX[a][b]+(x-a)*velocityX[a+1][1];
        u=(a+1-x)*velocityX[a][b+1]+(x-a)*velocityX[a+1][b+1];
        vx=(b+1-z)*t+(z-b)*u;
    
    }*/
    if((x<2) && (z<2))
    {
        vx=velocityX[1][1];
    }
    else if((x>nx-3) && (z>nz-3))
    {
        vx=velocityX[nx-2][nz-2];
    }
    else if(x<2)
    {
        vx=velocityX[1][b];

    }
    else if(z<2)
    {
        
        vx=velocityX[a][1];
    }
    else if(x>nx-3)
    {
        vx=velocityX[nx-2][b];
    }
    else if(z>nz-3)
    {
        vx=velocityX[a][nz-2];
    }// border ,fuck!!!
    else
    {
    Array1D<float> X;
    Array1D<float> Z;
    X.Resize(4);
    Z.Resize(4);
    X.Zero();
    Z.Zero();
    
    X[0]= -1*dz*dz*dz+2*dz*dz-dz;
    X[1]=  3*dz*dz*dz-5*dz*dz+2;
    X[2]= -3*dz*dz*dz+4*dz*dz+dz;
    X[3]= dz*dz*dz-dz*dz;
    
    Z[0]= -dx*dx*dx+2*dx*dx-dx;
    Z[1]= 3*dx*dx*dx-5*dx*dx+2;
    Z[2]= -3*dx*dx*dx+4*dx*dx+dx;
    Z[3]= dx*dx*dx-dx*dx;

    Array2D<float> H;
    H.Resize(4,4);
    H.Zero();
    H[0][0]=velocityX[a-1][b-1]*X[0]; H[0][1]=velocityX[a][b-1]*X[0]; H[0][2]=velocityX[a+1][b-1]*X[0]; H[0][3]=velocityX[a+2][b-1]*X[0];
    H[1][0]=velocityX[a-1][b]*X[1];   H[1][1]=velocityX[a][b]*X[1];   H[1][2]=velocityX[a+1][b]*X[1];   H[1][3]=velocityX[a+2][b]*X[1];
    H[2][0]=velocityX[a-1][b+1]*X[2]; H[2][1]=velocityX[a][b+1]*X[2]; H[2][2]=velocityX[a+1][b+1]*X[2]; H[2][3]=velocityX[a+2][b+1]*X[2];
    H[3][0]=velocityX[a-1][b+2]*X[3]; H[3][1]=velocityX[a][b+2]*X[3]; H[3][2]=velocityX[a+1][b+2]*X[3]; H[3][3]=velocityX[a+2][b+2]*X[3];


    vx=   (H[0][0]+H[1][0]+H[2][0]+H[3][0])*Z[0]*0.25
        +(H[0][1]+H[1][1]+H[2][1]+H[3][1])*Z[1]*0.25
        +(H[0][2]+H[1][2]+H[2][2]+H[3][2])*Z[2]*0.25
        +(H[0][3]+H[1][3]+H[2][3]+H[3][3])*Z[3]*0.25;
    }




}

void wavefront::getvz(Array2D<float> velocityZ)
{

    int a=floor(x);
    int b=floor(z);
    float dx=(x-a)/k;
    float dz=(z-b)/h;
	/*if (a >= nx - 1)
	{
		a = a - 2;
	}
	if (b >= nz - 1)
	{
		b = b - 2;
	}
    if((x<1)||(z<1))
    {
        float t=0;
        float u=0;
        t=(a+1-x)*velocityZ[a][b]+(x-a)*velocityZ[a][b+1];
        u=(a+1-x)*velocityZ[a+1][b]+(x-a)*velocityZ[a+1][b+1];
        vz=(b+1-z)*t+(z-b)*u;
    

    }
   else if((x+1>nx-1)||(z+1>nz-1)) 
    {
        float t=0;
        float u=0;
        t=(a+1-x)*velocityZ[a][b]+(x-a)*velocityZ[a+1][1];
        u=(a+1-x)*velocityZ[a][b+1]+(x-a)*velocityZ[a+1][b+1];
        vz=(b+1-z)*t+(z-b)*u;
    }*/

    if((x<2) && (z<2))
    {
        vz=velocityZ[1][1];
    }
    else if((x>nx-3) && (z>nz-3))
    {
        vz=velocityZ[nx-2][nz-2];
    }
    else if(x<2)
    {
        vz=velocityZ[1][b];

    }
    else if(z<2)
    {
        
        vz=velocityZ[a][1];
    }
    else if(x>nx-3)
    {
        vz=velocityZ[nx-2][b];
    }
    else if(z>nz-3)
    {
        vz=velocityZ[a][nz-2];
    }// border ,fuck!!!
    else
    {

    Array1D<float> X;
    Array1D<float> Z;
    X.Resize(4);
    Z.Resize(4);
    X.Zero();
    Z.Zero();
    
    X[0]= -1*dz*dz*dz+2*dz*dz-dz;
    X[1]=  3*dz*dz*dz-5*dz*dz+2;
    X[2]= -3*dz*dz*dz+4*dz*dz+dz;
    X[3]= dz*dz*dz-dz*dz;
    
    Z[0]= -dx*dx*dx+2*dx*dx-dx;
    Z[1]= 3*dx*dx*dx-5*dx*dx+2;
    Z[2]= -3*dx*dx*dx+4*dx*dx+dx;
    Z[3]= dx*dx*dx-dx*dx;

    Array2D<float> H;
    H.Resize(4,4);
    H.Zero();
    H[0][0]=velocityZ[a-1][b-1]*X[0]; H[0][1]=velocityZ[a][b-1]*X[0]; H[0][2]=velocityZ[a+1][b-1]*X[0]; H[0][3]=velocityZ[a+2][b-1]*X[0];
    H[1][0]=velocityZ[a-1][b]*X[1];   H[1][1]=velocityZ[a][b]*X[1];   H[1][2]=velocityZ[a+1][b]*X[1];   H[1][3]=velocityZ[a+2][b]*X[1];
    H[2][0]=velocityZ[a-1][b+1]*X[2]; H[2][1]=velocityZ[a][b+1]*X[2]; H[2][2]=velocityZ[a+1][b+1]*X[2]; H[2][3]=velocityZ[a+2][b+1]*X[2];
    H[3][0]=velocityZ[a-1][b+2]*X[3]; H[3][1]=velocityZ[a][b+2]*X[3]; H[3][2]=velocityZ[a+1][b+2]*X[3]; H[3][3]=velocityZ[a+2][b+2]*X[3];


    vz=   (H[0][0]+H[1][0]+H[2][0]+H[3][0])*Z[0]*0.25
        +(H[0][1]+H[1][1]+H[2][1]+H[3][1])*Z[1]*0.25
        +(H[0][2]+H[1][2]+H[2][2]+H[3][2])*Z[2]*0.25
        +(H[0][3]+H[1][3]+H[2][3]+H[3][3])*Z[3]*0.25;


    }



}

wavefront::~wavefront()
{

}


