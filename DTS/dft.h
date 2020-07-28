// h file


/***********************************************

      FileName: dft.h

        Author: stj
   Description: ---
 First  Create: 2020-01-02 14:45:39
 Last Modified: 2020-01-03 16:39:57

***********************************************/



#ifndef DFT_H 
#define DFT_H
#include "Array.hpp"
#include<vector>
#include<string>
#include"Ricker.h"
#include"cov.h"
using std::vector;
#define PI 3.14159265358
void read(vector<Array2D<float>> &record, vector<std::string> filename,int x,int z);
Array2D<float> focous_1(Array2D<float> record,Array2D<float> operators, int x,int z);
Array2D<float> integrate(int x,int z,int ixsrc,Array2D<float> CFP);
Array2D<float> mkcfp(vector<Array2D<float>> record,Array2D<float> operators,int x,int z, vector<int> ixsrc);
Array2D<float> mkdts(Array2D<float> CFP,Array2D<float> operators,int x,int z);
Array2D<float> idts(Array2D<float> CFP, Array2D<float> DTS, int x, int z);
Array2D<float> mkoperator(Array2D<float> deconvolution, int x, int z, int* ixsrc);
vector<float> insert_hermite(int* x, int* z, int nx);
void ricker_wavelet(int nt, float dt, float fpeak);
Array2D<float> mk_focusop(Array2D<float> operators, int x, int z);
Array2D<float> same_convolution(Array2D<float> operators, Array2D<float> record, int x, int z);
Array2D<float> same_dts(Array2D<float> operators, Array2D<float> record, int x, int z);
Array2D<float> full_convolution(Array2D<float> operators, Array2D<float> record, int x, int z);
Array2D<float> pad_zero(Array2D<float> demo, int x, int z);

#endif 
