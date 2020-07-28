// cpp file


/***********************************************

      FileName: dft.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-13 11:28:59
 Last Modified: 2020-07-13 11:31:23

***********************************************/

#include<stdio.h>
#include <iostream>
#include <cmath>
#include<complex>
#include<fftw3.h>
#include"dft.h"
#include<string>
#include<omp.h>
#include"cov.h"

using std::cout;
using std::endl;
using std::vector;

int to_int(float x)
{
	if (x - floor(x) < 0.5)
		return (int)x;
	else
		return (int)x + 1;

}


void read(vector<Array2D<float>> &record, vector<std::string> filename,int x,int z)
{
	int n1=record.size();
	int n2=filename.size();
	if(n1!=n2)
	{
		cout<<"两个数组都不求一样"<<endl;
		exit(0);
	}
	for(int i=0;i<n1;++i)
	{
		const char* a = filename[i].c_str();
		FILE *fp;
		fp=fopen(a,"rb");
		//cout << filename[i] << endl;
		
		if(fp==NULL)
		{
			cout<<"劳资打不开这个文件"<<endl;
			exit(0);
		}
		
		fread(&(record[i][0][0]),sizeof(float),x*z,fp);

		fclose(fp);

	}
	cout << "读取完毕" << endl;
}

Array2D<float> pad_zero(Array2D<float> demo, int x, int z)
{
	Array2D<float> result;
	result.Resize(x, 2 * z - 1);
	result.Zero();
	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < z; j++)
		{
			result[i][j] = demo[i][j];
		}
		for (int j = z; j < 2 * z - 1; j++)
		{
			result[i][j] = 0;
		}
	}
	return result;

}


Array2D<float> focous_1(Array2D<float> record,Array2D<float> operators, int x,int z)
{
	cout << "一次聚焦" << endl;
    Array1D<double> offset1D;
    Array1D<double> operator1D;
    Array1D<double> CFP1D;
    offset1D.Resize(z);
    operator1D.Resize(z);
    CFP1D.Resize(z);
    offset1D.Zero();
    operator1D.Zero();
    CFP1D.Zero();

    Array2D<float> CFP;
    CFP.Resize(x,z);
    CFP.Zero();

    fftw_complex *f_out,*g_out,*c_out;
    f_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);
    g_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);
    c_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);
	

	std::complex<double>* f_out_w = (std::complex<double> *) f_out;
	std::complex<double>* g_out_w = (std::complex<double>*) g_out;
	std::complex<double>* c_out_w = (std::complex<double>*) c_out;

    fftw_plan p1, p2, p3;
    p1 = fftw_plan_dft_r2c_1d(z, offset1D, f_out,FFTW_ESTIMATE);
    p2 = fftw_plan_dft_r2c_1d(z, operator1D, g_out,FFTW_ESTIMATE);
    p3 = fftw_plan_dft_c2r_1d(z, c_out, CFP1D,FFTW_ESTIMATE);//有待提高。一次plan
	//std::complex<double>* a = new std::complex<double>[z];

    for(int i=0;i<x;i++)
    {
		
        for(int j=0;j<z;j++)
        {
			offset1D[j]=record[i][j];
		
            operator1D[j]=operators[i][j];//聚焦算子
        }
        fftw_execute(p1);
        fftw_execute(p2);
        for(int j=0;j<z;j++)
        {
            if(j<z/2+1)
            {
                c_out_w[j]=f_out_w[j]*g_out_w[j];
            }
            else
            {
                c_out_w[j]=conj(f_out_w[z-i])*conj(g_out_w[z-i]);//fftw变成复数域的时候，输出为共轭对称
            }
        }
        fftw_execute(p3);
        for(int j=0;j<z;j++)
        {
            CFP[i][j]=CFP1D[j]/z;//fftw变成实数域时没有缩放
        }
    }

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_cleanup();
	
	cout << "聚焦完毕" << endl;

    return CFP;

}

Array2D<float> integrate(int x,int z,int ixsrc,Array2D<float> CFP)
{
	cout << "开始积分" << endl;
	Array2D<float> integ;
	integ.Resize(x,z);
	integ.Zero();
	for(int j=0;j<z;++j)
	{
		for(int i=0;i<x;++i)
		{
			integ[ixsrc][j]+=CFP[i][j];
		}
	}
	cout << "积分完毕" << endl;
	return integ;
}

Array2D<float> mkcfp(vector<Array2D<float>> record,Array2D<float> operators,int x,int z,vector<int> ixsrc )
{
	cout << "cfp" << endl;
	int n=record.size();
	int zz = 2 * z - 1;
	vector<Array2D<float>> sum(n);
	Array2D<float> CFP;
	//CFP.Resize(x,z);//same卷积所要的大小
	CFP.Resize(x,zz);//全卷积的大小
	CFP.Zero();

	//#pragma omp parallel for
	for(int i=0;i<n;i++)
	{
		
		cout << "第" << i << "次" << endl;
		//CFP=same_convolution(operators,record[i],x,z);//same卷积
		//CFP = full_convolution(operators, record[i], x, z);
		CFP = focous_1(pad_zero(record[i], x, z), pad_zero(operators, x, z),x,zz);
		
		sum[i]=integrate(x,zz,ixsrc[i],CFP);
		
		std::string file;
		std::string str = std::to_string(i);
		file = "first" + str + ".bin";
		const char* a = file.c_str();
		FILE* aa;
		aa = fopen(a, "wb");//第一次聚焦的结果
		fwrite(&(CFP[0][0]), sizeof(float), x * zz, aa);
		fclose(aa);
		
	}
	CFP.Zero();
	for(int i=0;i<n;++i)
	{
		for(int j=0;j<x;++j)
			for(int k=0;k<zz;++k)
				CFP[j][k] = sum[i][j][k] + CFP[j][k];
	}
	cout << "cfp over" << endl;
	return CFP;
}

Array2D<float> mkdts(Array2D<float> CFP,Array2D<float> operators,int x,int z)
{
	Array1D<double> DTS1D;
    Array1D<double> operator1D;
    Array1D<double> CFP1D;
    DTS1D.Resize(z);
    operator1D.Resize(z);
    CFP1D.Resize(z);
    DTS1D.Zero();
    operator1D.Zero();
    CFP1D.Zero();

    Array2D<float> DTS;
    DTS.Resize(x,z);
    DTS.Zero();

    fftw_complex *f_out,*g_out,*c_out;
    f_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);
    g_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);
    c_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * z);

	std::complex<double>* f_out_w = (std::complex<double>*) f_out;
	std::complex<double>* g_out_w = (std::complex<double>*) g_out;
	std::complex<double>* c_out_w = (std::complex<double>*) c_out;

    fftw_plan p1, p2, p3;
    p1 = fftw_plan_dft_r2c_1d(z, CFP1D, f_out,FFTW_ESTIMATE);
    p2 = fftw_plan_dft_r2c_1d(z, operator1D, g_out,FFTW_ESTIMATE);
    p3 = fftw_plan_dft_c2r_1d(z, c_out, DTS1D,FFTW_ESTIMATE);//有待提高。一次plan

    for(int i=0;i<x;i++)
    {
        for(int j=0;j<z;j++)
        {
            CFP1D[j]=CFP[i][j];
            operator1D[j]=operators[i][j];//逆时算子
        }
        fftw_execute(p1);
        fftw_execute(p2);
        for(int j=0;j<z;j++)
        {
			
			
            if(j<z/2+1)
            {
                c_out_w[j]=f_out_w[j]*conj(g_out_w[j]);
            }
            else
            {
                c_out_w[j]=conj(f_out_w[z-i])*g_out_w[z-i];//fftw变成复数域的时候，输出为共轭对称
            }
        }
        fftw_execute(p3);
        for(int j=0;j<z;j++)
        {
            DTS[i][j]=DTS1D[j]/z;//fftw变成实数域时没有缩放
        }
    }

    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_destroy_plan(p3);
    fftw_cleanup();

    return DTS;
}

Array2D<float> idts(Array2D<float> CFP, Array2D<float> DTS, int x, int z)
{
	cout << "开始反特么的卷积" << endl;
	Array1D<double> operators1D;
	Array1D<double> DTS1D;
	Array1D<double> CFP1D;
	operators1D.Resize(z);
	DTS1D.Resize(z);
	CFP1D.Resize(z);
	operators1D.Zero();
	DTS1D.Zero();
	CFP1D.Zero();

	Array2D<float> operators;
	operators.Resize(x, z);
	operators.Zero();

	fftw_complex* f_out, * g_out, * c_out;
	f_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * z);
	g_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * z);
	c_out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * z);

	std::complex<double>* f_out_w = (std::complex<double>*) f_out;
	std::complex<double>* g_out_w = (std::complex<double>*) g_out;
	std::complex<double>* c_out_w = (std::complex<double>*) c_out;

	fftw_plan p1, p2, p3;
	p1 = fftw_plan_dft_r2c_1d(z, CFP1D, f_out, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_r2c_1d(z, DTS1D, g_out, FFTW_ESTIMATE);
	p3 = fftw_plan_dft_c2r_1d(z, c_out, operators1D, FFTW_ESTIMATE);//有待提高。一次plan

	for (int i = 0; i < x; i++)
	{
		for (int j = 0; j < z; j++)
		{
			CFP1D[j] = CFP[i][j];
			DTS1D[j] = DTS[i][j];
		}
		fftw_execute(p1);
		fftw_execute(p2);
		for (int j = 0; j < z; j++)
		{
			if (j < z / 2 + 1)
			{
				c_out_w[j] = conj(g_out_w[j]/f_out_w[j]);
			}
			else
			{
				c_out_w[j] = g_out_w[j] / f_out_w[j];//fftw变成复数域的时候，输出为共轭对称
			}
		}
		fftw_execute(p3);
		for (int j = 0; j < z; j++)
		{
			operators[i][j] = operators1D[j] / z;//fftw变成实数域时没有缩放
		}
	}

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	fftw_cleanup();
	cout << "反卷积完成" << endl;
	return operators;
}

Array2D<float> mkoperator(Array2D<float> deconvolution, int x, int z, int* ixsrc)
{
	const int n = 3;
	cout << "n=" << n << endl;
	int maxz[n] = {0};
	for(int i=0;i<n;i++)
		for (int j = 0; j < z; j++)
		{
			
			if (deconvolution[ixsrc[i]][j] > maxz[i])
			{
				maxz[i] = j;
			}
		}

	for (int i = 0; i < n; i++)
		cout << "z=" << maxz[i] << endl;

	vector<float> offset(x);
	offset = insert_hermite(ixsrc, maxz, x);

	Array2D<float> test_op;
	test_op.Resize(x, z);
	test_op.Zero();
	for (int i = 0; i < x; i++)
	{
		test_op[i][to_int(offset[i])] = 1;

	}

	return test_op;
}


vector<float> insert_hermite(int* x,int* z, int nx)
{

	if ((sizeof(x)/sizeof(x[0]))!= sizeof(z) / sizeof(z[0]))
	{
		cout << "insert lost information" << endl;
	}
	int n = 3;

	
	vector<float> grad(n);//存导数，利用有限差分求取.

	vector<float> offset(nx);

	for (int i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			grad[i] = (z[i + 1] - z[i]) / (x[i + 1] - x[i]);
		}
		else if (i == n - 1)
		{
			grad[i] = (z[i] - z[i - 1]) / (x[i] - x[i - 1]);
		}
		else
		{
			grad[i] = (z[i + 1] - z[i - 1]) / (x[i + 1] - x[i - 1]);
		}
	}//一阶差分计算导数

	for (int i = 0; i < n - 1; ++i)
	{
		float x00 = x[i];
		float x11 = x[i + 1];
		float p;
		float a1, a0, b1, b0;
		float value;
		for (p = x[i]; p < x[i + 1]; p += 1.0)
		{
			float temp = (x11 - x00) * (x11 - x00);
			a0 = (x11 - 3 * x00 + 2 * p) * (x11 - p) * (x11 - p) / temp / (x11 - x00);
			a1 = (3 * x11 - x00 - 2 * p) * (p - x00) * (p - x00) / temp / (x11 - x00);
			b0 = (p - x00) * (p - x11) * (p - x11) / temp;
			b1 = (p - x00) * (p - x00) * (p - x11) / temp;
			value = a0 * z[i] + a1 * z[i + 1] + b0 * grad[i] + b1 * grad[i + 1];
			offset[p] = value;
		}
	}



	/*Array2D<float> a;
	a.Resize(nx, 2301);
	a.Zero();
	for (int i = 0; i < nx; i++)
	{
		//cout << offset[i] << endl;
		a[i][int(offset[i])] = 1;


	}
	FILE* F;
	F = fopen("hermite_test.bin", "wb");
	fwrite(&(a[0][0]), sizeof(float), nx * 2301, F);
	fclose(F);*/





	return offset;

}

Array2D<float> mk_focusop(Array2D<float> operators, int x, int z)
{
	Array2D<float> focusop;
	focusop.Resize(x, z);
	focusop.Zero();
	for (int i = 0; i < x; i++)
	{
		//cout << i << endl;
		for (int j = 0; j < z; j++)
		{
			focusop[i][j] = operators[i][z - j - 1];

			//if (i == 299)
				//cout << j << "  " << z - j - 1 << endl;


		}
	}
	FILE* cfp;
	cfp = fopen("focusop.bin", "wb");
	fwrite(&(focusop[0][0]), sizeof(float), x * z, cfp);
	fclose(cfp);
	return focusop;
}//得到聚焦算子





/////////////////////////////////

//   以下部分为卷积的尝试     //

///////////////////////////////
Array2D<float> same_convolution(Array2D<float> operators, Array2D<float> record, int x, int z)
{
	int n = 0;
	if (z % 2 == 0)
		n = (z + 2) / 2;
	else
		n = (z + 1) / 2;//定位卷积核

	Array2D<float> result;
	result.Resize(x, z);
	result.Zero();

	Array1D<float> operators1D;
	operators1D.Resize(z);
	Array1D<float> record1D;
	record1D.Resize(z);

	for (int i =0; i < x; ++i)
	{
		operators1D.Zero();
		record1D.Zero();
		for (int j = 0; j < z; ++j)
		{
			//cout << record1D[j] << endl;
			operators1D[j] = operators[i][j];
			record1D[j] = record[i][z-j-1];
			
			
		}//提取一维的向量进行卷积

		int tmp = n;
		
		for (int j = 0; j < n; ++j)
		{
			//cout <<"j="<< j << endl;
			for (int k=tmp-1;k<z;k++)
			{
				result[i][j] += record1D[k] * operators1D[k-tmp+1];
				//cout << k << ": "<<operators1D[k] << "*" << k - tmp +1<<" :"<< record1D[k - tmp + 1] << endl;
			}
			//cout << result[i][j] << endl;
			tmp--;
		}
		
		//cout << "hou/////////////////////////////////////////////////////////////////" << endl;
		//cout << tmp << endl;
		for (int j = n; j < z; ++j)
		{
			//cout << "j=" << j << endl;
			for (int k = 0; k < z-tmp-1; k++)
			{
				result[i][j] += record1D[k] * operators1D[tmp+k+1];
				//cout << k <<" :"<< operators1D[k] << "*" << tmp+k +1<<" :"<< record1D[tmp + k + 1]<< endl;
				//cout << result[i][j] << endl;
				
			}
			//cout << endl;
			//cout << result[i][j] << endl;
		
			
			tmp++;
		}

	}
	return result;


}


Array2D<float> same_dts(Array2D<float> operators, Array2D<float> record, int x, int z)
{
	int n = 0;
	if (z % 2 == 0)
		n = (z + 2) / 2;
	else
		n = (z + 1) / 2;//定位卷积核

	Array2D<float> result;
	result.Resize(x, z);
	result.Zero();

	Array1D<float> operators1D;
	operators1D.Resize(z);
	Array1D<float> record1D;
	record1D.Resize(z);

	for (int i = 0; i < x; ++i)
	{
		operators1D.Zero();
		record1D.Zero();
		for (int j = 0; j < z; ++j)
		{
			//cout << record1D[j] << endl;
			operators1D[j] = operators[i][j];
			record1D[j] = record[i][j ];


		}//提取一维的向量进行卷积

		int tmp = n;

		for (int j = 0; j < n; ++j)
		{
			//cout <<"j="<< j << endl;
			for (int k = tmp - 1; k < z; k++)
			{
				result[i][j] += record1D[k] * operators1D[k - tmp + 1];
				//cout << k << ": "<<operators1D[k] << "*" << k - tmp +1<<" :"<< record1D[k - tmp + 1] << endl;
			}
			//cout << result[i][j] << endl;
			tmp--;
		}

		//cout << "hou/////////////////////////////////////////////////////////////////" << endl;
		//cout << tmp << endl;
		for (int j = n; j < z; ++j)
		{
			//cout << "j=" << j << endl;
			for (int k = 0; k < z - tmp - 1; k++)
			{
				result[i][j] += record1D[k] * operators1D[tmp + k + 1];
				//cout << k <<" :"<< operators1D[k] << "*" << tmp+k +1<<" :"<< record1D[tmp + k + 1]<< endl;
				//cout << result[i][j] << endl;

			}
			//cout << endl;
			//cout << result[i][j] << endl;


			tmp++;
		}

	}
	return result;


}

Array2D<float> full_convolution(Array2D<float> operators, Array2D<float> record, int x, int z)
{
	Array2D<float> result;
	result.Resize(x, 2*z-1);
	result.Zero();
	for (int i = 0; i < x; ++i)
	{
		
		convolve_cwp(z,0, record[i], z, 0, operators[i],2*z-1,0,result[i]);

	}
	return result;
}

//g++ dft.cpp -o test -I/opt/OSF/FFTW/include/ -L/opt/OSF/FFTW/lib/ -lfftw3 -lm



