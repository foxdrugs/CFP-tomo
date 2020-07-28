// cpp file


/***********************************************

      FileName: changeM.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-12 15:21:10
 Last Modified: 2020-07-12 15:31:04

***********************************************/
#include <iostream>
#include <fstream>
#include "Array.hpp"
#include <sstream>
#include<string>
using namespace std;
int main()
{
    int nx=348;
    int nz=648;
    Array2D<float> velocity;
    velocity.Resize(nx,nz);
    velocity.Zero();

    FILE* oper;
	oper = fopen("velocity.bin", "rb");
	fread(&velocity[0][0], sizeof(float), nx*nz, oper);
	fclose(oper);

    for(int i=0;i<nx;++i)
    {
        for(int j=0;j<nz;++j)
        {
            velocity[i][j]=velocity[100][10];
        }
    }
    FILE* fp;
	fp = fopen("ch_velocity.bin", "wb");
	fwrite(&(velocity[0][0]), sizeof(float), nx * nz, fp);
	fclose(fp);
}
