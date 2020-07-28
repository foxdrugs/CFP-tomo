// cpp file


/***********************************************

      FileName: main.cpp

        Author: stj
   Description: ---
 First  Create: 2019-12-17 20:38:04
 Last Modified: 2019-12-29 21:13:13

***********************************************/
#include "dft.h"
#include<iostream>
#include<string>
#include <sstream>
#include"cov.h"
using std::cout;
using std::endl;
using std::vector;
int main(int argc,char *argv[])
{
	int x = 600;
	int z = 1500;
	int zz = 2 * z - 1;
	vector<Array2D<float>>	offset;
	int izsrc = 0;
	
	vector<std::string> filename;
	vector<int> IXSRC;
	for (int ixsrc = 10; ixsrc < 501; ixsrc += 10)
	{
		std::ostringstream off;
		off << "off" << "(" << ixsrc << "," << izsrc << ")" << ".bin";
		std::string Filename = off.str();
		cout << Filename << endl;
		filename.push_back(Filename);

		Array2D<float> record;
		record.Resize(x, z);
		record.Zero();
		offset.push_back(record);

		IXSRC.push_back(ixsrc);
	}
	
	
	read(offset, filename, x, z);
	Array2D<float> operators;
	operators.Resize(x, z);
	operators.Zero();
	FILE* a;
	a = fopen("operator(306,145).bin", "rb");
	fread(&(operators[0][0]), sizeof(float), x * z, a);
	fclose(a);

	Array2D<float> CFP;
	CFP.Resize(x, zz);
	CFP.Zero();
	CFP = mkcfp(offset,mk_focusop(operators,x, z),x, z,IXSRC);


	FILE* cfp;
	cfp = fopen("cfp_mess.bin", "wb");//积分后的道集
	fwrite(&(CFP[0][0]), sizeof(float),x* zz, cfp);
	fclose(cfp);

	/*
	Array2D<float> DTS;
	DTS.Resize(x, 2*z-1);
	DTS.Zero();
	for (int i = 0; i < x; i++)
	{
		Array1D<float> CFP1D;
		CFP1D.Resize(2 * z - 1);
		CFP1D.Zero();
		for (int j = 0; j < 2 * z - 1; j++)
		{
			CFP1D[j] = CFP[i][2 * z - 1 - j];
		}
		convolve_cwp(z, 0, operators[i], 2 * z - 1,0,CFP1D,2*z-1,0,DTS[i]);
	}
	FILE* w;
	w = fopen("dts.bin", "wb");//dts面板
	fwrite(&(DTS[0][0]), sizeof(float), x * (2*z-1), w);
	fclose(w);
	*/
	
	Array2D<float> DTS;
	DTS.Resize(x, zz);
	DTS.Zero();
	DTS = mkdts(pad_zero(operators, x, z), CFP, x, zz);
	FILE* w;
	w = fopen("dts.bin", "wb");//dts面板
	fwrite(&(DTS[0][0]), sizeof(float), x * zz, w);
	fclose(w);
	


	/*Array2D<float> iop;
	iop.Resize(x, z);
	iop.Zero();
	iop = idts(CFP, DTS, x, z);


	Array2D<float> test_op;
	test_op.Resize(x, z);
	test_op.Zero();
	test_op = mkoperator(iop,x, z,ixsrc);

	

	FILE* aa;
	aa = fopen("test_op.bin", "wb");//插值后的结果
	fwrite(&(test_op[0][0]),sizeof(float),x* z,aa);
	fclose(aa);

	FILE* bb;
	bb = fopen("iop.bin", "wb");//反褶积后的结果
	fwrite(&(iop[0][0]), sizeof(float), x * z, bb);
	fclose(bb);*/
	return 0;
}






int __main()
{
	int x = 600;
	int z = 1500;
	Array2D<float> operators;
	operators.Resize(x, z);
	operators.Zero();
	FILE* a;
	a = fopen("operator(306,145).bin", "rb");
	fread(&(operators[0][0]), sizeof(float), x * z, a);
	fclose(a);
	
	Array2D<float> record;
	record.Resize(x, z);
	record.Zero();
	FILE* b;
	b = fopen("off(500,0).bin", "rb");
	fread(&(record[0][0]), sizeof(float), x * z, b);
	fclose(b);

	Array2D<float> CFP;
	CFP.Resize(x, z);
	CFP.Zero();

	/*Array2D<float> test;
	test.Resize(20, 10);
	test.Zero();
	for (int i = 0; i < 20; ++i)
		for (int j = 0; j < 10; ++j)
			test[i][j] = 1;*/

	CFP = same_convolution(mk_focusop(operators,x,z), record, x, z);
	

	FILE* aa;
	aa = fopen("focus_test.bin", "wb");//第一次聚焦的结果
	fwrite(&(CFP[0][0]), sizeof(float), x * z, aa);
	fclose(aa);
	return 0;


}



