// cpp file


/***********************************************

      FileName: model.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-15 13:41:43
 Last Modified: 2020-07-26 16:30:49

***********************************************/

#include <iostream>
#include <fstream>
#include "Array.hpp"
using namespace std;
int main()
{
	int n1 = 300;
	int n2 = 600;
	int bclenth = 24;
	n2 = n2 + 2 * bclenth;
	n1 = n1 + 2 * bclenth;
	Array2D<float> velocity;
	velocity.Resize(n1, n2);
	velocity.Zero();
	cout << "n1=" << n1 << "  n2=" << n2 << endl;
	

	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2; j++)
		{
			
			velocity[i][j] = 1500;
		}
		
	}

	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < 240; j++)
		{

			velocity[i][j] = 1000;
		}
		cout << "i=" << i << "  j=" << 174 + 0.3 * i << endl;

	}


	for (int i = 0; i < n1; i++)
	{
		for (int j = 474; j < n2; j++)
		{

			velocity[i][j] = 2100;
		}
		

	}


	FILE* fp;
	fp = fopen("velocity.bin", "wb");
	fwrite(&(velocity[0][0]), sizeof(float), n1 * n2, fp);
	fclose(fp);


	return 0;
}




