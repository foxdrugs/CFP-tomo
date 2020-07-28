// cpp file


/***********************************************

      FileName: smooth.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-12 15:53:51
 Last Modified: 2020-07-13 19:51:19

***********************************************/

#include<iostream>
#include"Array.hpp"
#include <cstdio>
Array2D<float> smooth(Array2D<float> velocity, int nx, int nz, int n, int border)
{
	//float w = 0.5;//光滑算子;
	Array2D<float> result;
	result.Resize(nx - 2 * border, nz - 2 * border);
	result.Zero();

	Array2D<float> tmp;
	tmp.Resize(nx, nz);
	tmp.Zero();
	tmp = velocity;

	for (int k = 0; k < n; k++)
	{
		std::cout << k << std::endl;
		for (int i = border-5; i < nx - border+5; i++)
		{
			for (int j = border-5; j < nz - border+5; ++j)
			{
				tmp[i][j] = (velocity[i][j - 1] + velocity[i][j] + velocity[i][j + 1]+ velocity[i][j + 2]+ velocity[i][j -2])/5.0;
				
			}
		}
		velocity = tmp;
        
        for (int i = border-5; i < nx - border+5; i++)
		{
			for (int j = border-5; j < nz - border+5; ++j)
			{
				tmp[i][j] = (velocity[i-1][j] + velocity[i][j] + velocity[i+1][j]+ velocity[i+2][j]+ velocity[i-2][j])/5.0;
				
			}
		}
        velocity=tmp;

	}

	for (int i = border; i < nx - border; i++)
	{
		for (int j = border; j < nz - border; ++j)
		{
			result[i - border][j - border] = velocity[i][j];
		}
	}
	return result;
}
int main()
{
	int nx=348;
	int nz=648;
	int border = 24;
	Array2D<float> velocity;
	velocity.Resize(nx, nz);
	velocity.Zero();
	FILE* s;
	s = fopen("velocity.bin", "rb");
	fread(&(velocity[0][0]), sizeof(float), nz * nx, s);
	fclose(s);

	Array2D<float> result;
	result.Resize(nx - 2*border, nz - 2*border);
	result.Zero();
	result = smooth(velocity, nx, nz, 50, border);

	FILE* c;
	c = fopen("smooth_velocity.bin", "wb");
	fwrite(&(result[0][0]), sizeof(float), (nx - 2 * border)* (nz - 2 * border), c);
	fclose(c);



}

