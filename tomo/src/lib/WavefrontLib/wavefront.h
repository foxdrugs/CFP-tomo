// h file


/***********************************************

      FileName: wavefront.h

        Author: stj
   Description: ---
 First  Create: 2019-07-08 15:24:16
 Last Modified: 2020-07-10 22:14:42

***********************************************/

#ifndef WAVEFRONT_H
#define WAVEFRONT_H

#include<iostream>
#include"Array.hpp"
#include <string>
#include <cstdio>
using namespace std;
class wavefront
{
public:
    wavefront(double dx,double dz,double Nx,double Nz)
    {
        k=dx;
        h=dz;
        nx=Nx;
        nz=Nz;
        
    }
    void getv(Array2D<float> velocity);
    void getvx(Array2D<float> velocityX);
    void getvz(Array2D<float> velocityZ);

    double k,h;//distance of the cube
    int nz,nx;//numbers of the cube
    double x,z,t;//time and position;
    double x_slow, z_slow;//slowness
    double v, vx,vz;//velocity and derivative
    double angle;
    ~wavefront();

};



#endif 

