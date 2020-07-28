// h file


/***********************************************

      FileName: Raytrace.h

        Author: stj
   Description: ---
 First  Create: 2020-01-02 14:45:39
 Last Modified: 2020-07-10 22:13:44

***********************************************/

#ifndef RAYTRACE_H 
#define RAYTRACE_H
#include"wavefront.h"
#include <Ray.hpp>
void rule(wavefront now, wavefront &nex);
void rk4(const wavefront now, wavefront &nex,
        float dx,float dz,float nx,float nz,
        float dt,
        Array2D<float> velocity,Array2D<float> velocityX, Array2D<float> velocityZ);
Array2D<float> DerivateX(Array2D<float> velocity,int nx,int nz,int dx);
Array2D<float> DerivateZ(Array2D<float> velocity,int nx,int nz,int dz);
void Raytrace(ray& path, int detector, wavefront now, wavefront nex,
	float dx, float dz, float nx, float nz,
	float dt,
	Array2D<float> velocity, Array2D<float> velocityX, Array2D<float> velocityZ);


#endif 


