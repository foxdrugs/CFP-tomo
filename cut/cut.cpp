// cpp file


/***********************************************

      FileName: cut.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-12 15:16:50
 Last Modified: 2020-07-12 15:31:08

***********************************************/


#include <iostream>
#include <fstream>
#include "Array.hpp"
#include <sstream>
#include<string>
using namespace std;
int main()
{
	int n1 = 300;
	int n2 = 2000;
	
	Array2D<float> offset;
	offset.Resize(n1, n2);
	offset.Zero();


	Array2D<float> direct;
	direct.Resize(n1, n2);
	direct.Zero();

	Array2D<float> result;
	result.Resize(n1, n2);
	result.Zero();

	int izsrc = 0;

	for (int ixsrc = 0; ixsrc < n1; ixsrc += 10)
	{
		result.Zero();
		direct.Zero();
		offset.Zero();

		ostringstream off;
		off<<"offset" << "(" << ixsrc << "," << izsrc << ")" << ".bin";
		string Filename = off.str();
		cout << Filename << endl;
		const char* filename = Filename.data();
		FILE* oper;
		oper = fopen(filename, "rb");
		fread(&offset[0][0], sizeof(float), n1*n2, oper);
		fclose(oper);

		ostringstream direc;
		direc << "direct" << "(" << ixsrc << "," << izsrc << ")" << ".bin";
		string dir_Filename = direc.str();
		cout << dir_Filename << endl;
		const char* dir_filename = dir_Filename.data();
		FILE* dir;
		dir = fopen(dir_filename, "rb");
		fread(&direct[0][0], sizeof(float), n1 * n2, dir);
		fclose(dir);

		for(int i=0;i<n1;i++)
			for (int j = 0; j < n2; ++j)
			{
				result[i][j] = offset[i][j] - direct[i][j];
			}

		ostringstream re;
		re << "off" << "(" << ixsrc << "," << izsrc << ")" << ".bin";
		string re_Filename = re.str();
		cout << re_Filename << endl;
		const char* re_filename = re_Filename.data();
		FILE* fp;
		fp = fopen(re_filename, "wb");
		fwrite(&(result[0][0]), sizeof(float), n1 * n2, fp);
		fclose(fp);

	}
	


	


	return 0;
}
