// cpp file


/***********************************************

      FileName: Raytrace.cpp

        Author: stj
   Description: ---
 First  Create: 2020-01-02 15:03:55
 Last Modified: 2020-07-21 10:49:58

***********************************************/
#include<cstdio>
#include<iostream>
#include<string>
#include<fstream>
#include<iomanip>
#include<cmath>
#include"Array.hpp"
#include"wavefront.h"
#include"operator.h"
#define PI 3.1415926

void rule(wavefront now, wavefront &nex) 
{

    float v2 = now.v*now.v;
    nex.x = v2 * now.x_slow;
    //cout<<"rule: now"<<now.x_slow<<endl;
    
    nex.z = v2 * now.z_slow;
    nex.x_slow = -now.vx/now.v;
    //cout<<"rule: nex"<<nex.x_slow<<endl;
    //cout<<now.vx<<endl;
    //cout<<now.v<<endl;
    nex.z_slow = -now.vz/now.v;

}

void rk4(const wavefront now, wavefront &nex,float dx,float dz,float nx,float nz,float dt,Array2D<float> velocity,Array2D<float> velocityX, Array2D<float> velocityZ)
{
    wavefront tmp1(dx,dz,nx,nz);
    rule(now,tmp1);
    wavefront tmpp1(dx,dz,nx,nz);

    tmpp1.x=now.x+0.5*dt*tmp1.x;
    tmpp1.z=now.z+0.5*dt*tmp1.z;
    tmpp1.x_slow=now.x_slow+0.5*dt*tmp1.x_slow;
    tmpp1.z_slow=now.z_slow+0.5*dt*tmp1.z_slow;
	//cout <<"tmpp1:"<< tmpp1.z_slow << endl;
    tmpp1.getv(velocity);
    tmpp1.getvx(velocityX);
    tmpp1.getvz(velocityZ);
    wavefront tmp2(dx,dz,nx,nz);
    rule(tmpp1,tmp2);
    wavefront tmpp2(dx,dz,nx,nz);
    tmpp2.x=now.x+0.5*dt*tmp2.x;
    tmpp2.z=now.z+0.5*dt*tmp2.z;
    tmpp2.x_slow=now.x_slow+0.5*dt*tmp2.x_slow;
    tmpp2.z_slow=now.z_slow+0.5*dt*tmp2.z_slow;
	//cout << "tmpp2:" << tmpp2.z_slow << endl;
    tmpp2.getv(velocity);
    tmpp2.getvx(velocityX);
    tmpp2.getvz(velocityZ);

    wavefront tmp3(dx,dz,nx,nz);

    rule(tmpp2,tmp3);
    wavefront tmpp3(dx,dz,nx,nz);
    tmpp3.x=now.x+dt*tmp3.x;
    tmpp3.z=now.z+dt*tmp3.z;

    tmpp3.x_slow=now.x_slow+dt*tmp3.x_slow;
    tmpp3.z_slow=now.z_slow+dt*tmp3.z_slow;
    //cout << "tmpp3:" << tmpp3.z_slow << endl;
    //cout << "tmpp3:" << tmpp3.x<<"  "<<tmpp3.z << endl;
    tmpp3.getv(velocity);
    //cout<<"tmpp3.v="<<tmpp3.v<<endl;
    tmpp3.getvx(velocityX);
    tmpp3.getvz(velocityZ);
    wavefront tmp4(dx,dz,nx,nz);
    rule(tmpp3,tmp4);
    //cout<<"tmp4"<<tmp4.x_slow<<endl;
    nex.x=now.x+(dt/6.0)*(tmp1.x+2.0*tmp2.x+2.0*tmp3.x+tmp4.x);
    nex.z=now.z+(dt/6.0)*(tmp1.z+2.0*tmp2.z+2.0*tmp3.z+tmp4.z);

    nex.x_slow=now.x_slow+(dt/6.0)*(tmp1.x_slow+2.0*tmp2.x_slow+2.0*tmp3.x_slow+tmp4.x_slow);
    //cout<<"rk: x_slow"<<nex.x_slow<<endl;
    //cout<<tmp1.x_slow<<"  "<<tmp2.x_slow<<"  "<<tmp3.x_slow<<"  "<<"  "<<tmp4.x_slow<<endl;
    nex.z_slow=now.z_slow+(dt/6.0)*(tmp1.z_slow+2.0*tmp2.z_slow+2.0*tmp3.z_slow+tmp4.z_slow);
   
}

Array2D<float> DerivateX(Array2D<float> velocity,int nx,int nz,int dx)
{
    Array2D<float> speedX;
    speedX.Resize(nx,nz);
    speedX.Zero();
    
    for(int i=0;i<nz;i++)
    {
        speedX[0][i]=(velocity[1][i]-velocity[0][i])*(1/dx);
        speedX[nx-1][i]=(velocity[nx-1][i]-velocity[nx-2][i])*(1/dx);
    }
    
    for(int i=0;i<nz;i++)
        for(int j=1;j<nx-1;j++)
        {
            speedX[j][i]=(velocity[j+1][i]-velocity[j-1][i])*(0.5/dx);
        }
    return speedX;
    
/*

    if(x==0)
    {

        der_x=(velocity[x+1][z]-velocity[x][z])*(0.5/dx);
    }
    else if(x==nx-1)
        der_x=(velocity[x][z]-velocity[x-1][z])*(0.5/dx);
    else
    der_x=(velocity[x+1][z]-velocity[x-1][z])*(0.5/dx);

    return der_x;
    */
}

Array2D<float> DerivateZ(Array2D<float> velocity,int nx,int nz,int dz)
{   
    Array2D<float> speedZ;
    speedZ.Resize(nx,nz);
    speedZ.Zero();
    
    for(int i=0;i<nx;i++)
    {
        speedZ[i][0]=(velocity[i][1]-velocity[i][0])*(1/dz);
        speedZ[i][nz-1]=(velocity[i][nz-1]-velocity[i][nz-2])*(1/dz);
    }
    
    for(int i=0;i<nx;i++)
        for(int j=1;j<nz-1;j++)
        {
            speedZ[i][j]=(velocity[i][j+1]-velocity[i][j-1])*(0.5/dz);
        }

    return speedZ;
/*
    if(z==0)
        der_z=(velocity[x][z+1]-velocity[x][z])*(0.5/dz);
    else if(z==nz-1)
        der_z=(velocity[x][z]-velocity[x][z-1])*(0.5/dz);
    else
    der_z=(velocity[x][z+1]-velocity[x][z-1])*(0.5/dz);

    return der_z;
    */
}

void Raytrace(ray &path,int detector,wavefront now, wavefront nex,
        float dx,float dz,float nx,float nz,
        float dt,
        Array2D<float> velocity,Array2D<float> velocityX, Array2D<float> velocityZ)
{
    int a=floor(path.Xposition)-1;
    if(a<0)
    {
        a=0;
    }
	//cout << a << endl;
    float hole=(velocity[a][detector]*dt+0.02)/2.0;
    //cout<<hole<<endl;
    //float hole=1.0;
    now.x=path.Xposition;
    now.z=path.Zposition;
    now.getv(velocity);
    now.getvx(velocityX);
    now.getvz(velocityZ);

    now.x_slow=(1.0/now.v)*cos(path.angle*PI/180);
    now.z_slow=(1.0/now.v)*sin(path.angle*PI/180);
    double t = 0;
    for(int i=0;i<666666;i++)
    {
        //cout<<i<<endl;
		
        if(now.x<=2.4||now.x>=nx-2)
            break;
        if(now.z<detector||now.z-detector<hole)
        {
			//cout <<"time is"<< i*5 << endl;
			//cout << "x=" << now.x << " " << "z=" << now.z << endl;
			//float a = sqrt((now.x - 20) * (now.x - 20) + (now.z - 60) * (now.z - 60));
			//cout << "dis(20,60)=" << a*5/2 << endl;
			//cout << "t=" << t/dt << endl;
			
			path.gettime(t/dt);
            break;
        }
        rk4(now,nex,dx,dz,nx,nz,dt,velocity,velocityX,velocityZ);
		
		//cout << "v=" << now.v << "   " << "xslow=" << now.x_slow <<"   "<<"zslow=" << now.z_slow << endl;
		//cout << "t=" << t / dt / 10 << endl;
		//cout << nex.x-now.x << "  " << nex.z-now.z << endl;
		t = t + sqrt(dx * (nex.x - now.x) * dx * (nex.x - now.x) + dz * (nex.z - now.z) * dz * (nex.z - now.z)) / now.v;
        now.x=nex.x;
        now.z=nex.z;
        path.pushpath(now.x,now.z);
        now.x_slow=nex.x_slow;
        //cout<<"fuck"<<now.x_slow<<endl;
        now.z_slow=nex.z_slow;
        now.getv(velocity);
        now.getvx(velocityX);
        now.getvz(velocityZ);
	//cout << "x=" << now.x << " " << "z=" << now.z << " " << i << endl;
    //cout << "v=" << now.v << "   " << "xslow=" << now.x_slow <<"   "<<"zslow=" << now.z_slow << endl;
    //cout<<endl;
		/*
		if (i % 100 == 0)
		{
			cout << "x=" << now.x << " " << "z=" << now.z << " " << i << endl;
			cout << "v=" << now.v << "   " << "xslow=" << now.x_slow << "   " << "zslow=" << now.z_slow << endl;
			cout << "t=" << t / dt / 10 << endl;
		}
		*/
    }

}

