// h file


/***********************************************

      FileName: Sirt.h

        Author: stj
   Description: ---
 First  Create: 2020-07-12 21:04:50
 Last Modified: 2020-07-16 15:45:02

***********************************************/



#ifndef SIRT_H 
#define SIRT_H 
#include "Array.hpp"
#include<vector>
#include<string>
#include"Point.hpp"
using std::vector;

Array2D<float> sirt(vector<Array2D<float>> velocity, int nx, int nz);
Array2D<float> constraint(Array2D<float> velocity,vector<float> l_up, vector<float> l_down,int nx, int nz);

vector<float> interpola(int nx, vector<point> layer);
Array2D<float> smooth(Array2D<float> velocity, int nx, int nz,int n);

#endif 
