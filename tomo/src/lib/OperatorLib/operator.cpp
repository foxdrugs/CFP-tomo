// cpp file


/***********************************************

      FileName: operator.cpp

        Author: stj
   Description: ---
 First  Create: 2019-12-27 15:44:40
 Last Modified: 2020-07-26 15:32:05

***********************************************/
#include"operator.h"
#include<math.h>
#define PI 3.14159265358979
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

vector<float> operators::insert_hermite(vector<ray*> tmp_ray)
{
	if (n != (int)tmp_ray.size())
	{
		cout << "insert lost information" << endl;
	}

	vector<float> t(n);
	vector<float> x(n);
	vector<float> grad(n);//存导数，利用有限差分求取.

	for (int i = 0; i < n; ++i)
	{

		t[i] = tmp_ray[i]->time;
		//cout << "t=" << t[i] << endl;
		x[i] = tmp_ray[i]->getx();
		//cout << "x=" << x[i] << endl;

	}
	vector<float> offset(nx);

	for (int i = 0; i < n; ++i)
	{
		if (i == 0)
		{
			grad[i] = (t[i + 1] - t[i]) / (x[i + 1] - x[i]);
		}
		else if (i == n - 1)
		{
			grad[i] = (t[i] - t[i - 1]) / (x[i] - x[i - 1]);
		}
		else
		{
			grad[i] = (t[i + 1] - t[i - 1]) / (x[i + 1] - x[i - 1]);
		}
	}//一阶差分计算导数

	for (int i = 0; i < n-1; ++i)
	{
		float x00 = x[i];
		float x11 = x[i + 1];
		float p;
		float a1, a0, b1, b0;
		float value;
		for (p = x[i]; p < x[i + 1]; p += 1.0)
		{
			float temp = (x11 - x00) * (x11 - x00);
			a0 = (x11-3*x00+2*p)*(x11-p)*(x11-p)/temp/(x11-x00);
			a1=(3*x11-x00-2*p)*(p-x00)*(p-x00)/temp / (x11 - x00);
			b0 = (p - x00) * (p - x11) * (p - x11) / temp;
			b1 = (p - x00) * (p - x00) * (p - x11) / temp;
			value = a0 * t[i] + a1 * t[i + 1] + b0 * grad[i] + b1 * grad[i + 1];
			offset[p] = value;
		}
	}
	


	Array2D<float> a;
	a.Resize(nx, 2000);
	a.Zero();
	for (int i = 0; i < nx; i++)
	{
		a[i][int(offset[i])] = 1;


	}
	FILE* F;
	F = fopen("hermite_test.bin", "wb");
	fwrite(&(a[0][0]),sizeof(float),nx*2000,F);
	fclose(F);



	

	return offset;

}

operators::operators(float xposition, float zposition,int Nx,int Nz,float Dt)
{
    Xposition=xposition;
    Zposition=zposition;
	n=0;
	nx = Nx;
	nz = Nz;
	dt = Dt;
}

operators::operators(operators &p)
{
    Xposition=p.Xposition;
    Zposition=p.Zposition;
    ray_group=p.ray_group;
	n=p.n;
	nx = p.nx;
	nz = p.nz;
}

void operators::push_back(ray path)
{
    ray *tmp=new ray(path);
    ray_group.push_back(tmp);
	n++;

}

vector<ray*> operators::getray_group()
{
    return ray_group;
}

vector<float> operators::insert(vector<ray*> tmp_ray)
{
    if(n!= (int)tmp_ray.size())
	{
		cout<<"insert lost information"<<endl;
	}
	
    vector<float> t(n);
    vector<float> x(n);
	
    for(int i=0;i<n;++i)
    {
		
        t[i]= tmp_ray[i]->time;
		//cout << "t=" << t[i] << endl;
        x[i]= tmp_ray[i]->getx();
		//cout << "x=" << x[i] << endl;
		
    }

	vector<float> offset(nx);
	
	for(int i=1;i<n;i++)
	{
		float grad=(t[i]-t[i-1])/(x[i]-x[i-1]);
		
		int num = 1;
		for(int j=x[i-1];j<x[i];j++)
		{
			offset[j]=t[i-1]+num*grad;
			num++;
			
		}
		
		
	}


	/*Array2D<float> a;
	a.Resize(nx, 2301);
	a.Zero();
	for (int i = 0; i < nx; i++)
	{
		a[i][int(offset[i])] = 1;
		

	}
	FILE* F;
	F = fopen("test.bin", "wb");
	fwrite(&(a[0][0]),sizeof(float),nx*2301,F);
	fclose(F);*/
	return offset;
	
}



vector<float> operators::perturb(int left, int right, int up, int down,Array2D<float> velocity,Array2D<float> record,int dx,int dz)
{


    if(Xposition-left<0)
    {
        left=Xposition;
    }
    if(Xposition+right>nx)
    {
        right=nx-Xposition-1;
    }
    if(Zposition-up<0)
    {
        up=Zposition;
    }
    if(Zposition+down>nz)
    {
        down=nz-Zposition-1;
    }

	vector<float> offset(nx);
	float amplitude=-10000000;//叠加振幅
	
	vector<ray*> tmp_ray(n);//射线组的中间变量
	for (int i = 0; i < n; i++)
	{
		tmp_ray[i] = new ray(*ray_group[i]);
	}

	int a = 0, b = 0;//坐标的偏移量
	vector<float> dt_v(n);
	float dt_h;
	int flag_h=0;
	int flag_v=up;

    
	for(int j=0;j<up+down+1;++j)
	{
		
		//dt_v=0;
		for(int k=0;k<n;++k)
		{
			dt_v[k]=-(flag_v*dz*sin((360-tmp_ray[k]->angle)*PI/180))/velocity[(int)Xposition][(int)Zposition-1]/dt;
		}

		flag_h=left;
		
		for(int i=0;i<left+right+1;++i)
		{
			float tmp_amplitude = 0;//临时变量，每次循环都归零

			



			//cout << endl;
			//cout <<"////////////////////////////////////////////////////////////////////"<< endl;
			//cout << "x=" << Xposition - flag_h << " " << "z=" << Zposition - flag_v << endl;
			for(int k=0;k<n;++k)
			{
				dt_h=(float(flag_h)*float(dx)*cos((360.0-tmp_ray[k]->angle)*PI/180))/velocity[(int)Xposition][(int)Zposition-1]/dt;
				tmp_ray[k]->time=ray_group[k]->time+dt_v[k]+dt_h;
				//cout << "ray" << k << "=" << tmp_ray[k]->time<<"  x="<< tmp_ray[k]->getx() << " 原射线"<<ray_group[k]->time<<"  角度："<< 360-ray_group[k]->angle<<endl;
				//cout << "扰动前-扰动后=" << ray_group[k]->time - tmp_ray[k]->time << "  角度：" << 360 - ray_group[k]->angle << endl;
				//cout << endl;
				//cout << "水平扰动量=" << dt_h << "  角度：" << 360 - ray_group[k]->angle << endl;
				//cout << "竖直扰动量=" << dt_v[k] <<  endl;
				//cout << "cos=" << cos((360 - tmp_ray[k]->angle) * PI / 180) << endl;
				//cout << "总扰动量=" << dt_v[k] + dt_h << endl;
				
			}//扰动后求取时间
			
			offset= insert_hermite(tmp_ray);//插值函数得到炮记录
			

			for(int i=0;i<nx;++i)
			{
				tmp_amplitude += record[i][to_int(offset[i])];
				//cout << "i=" << i << "  t=" << to_int(offset[i]) << endl;
				
			}//叠加振幅

			
			
			if(tmp_amplitude>amplitude)
			{
				amplitude=tmp_amplitude;
				//ray_group=tmp_ray;
				a=flag_h;
				b=flag_v;
			}//比较振幅大小
			//cout<< "  amplitude:" << tmp_amplitude << endl;
			
			flag_h--;
		}
		flag_v--;
	}
	
	
	for (int i = 0; i < n; i++)
	{
		delete(tmp_ray[i]);

	}
	

	vector<float> mes(3);
	mes[0] = Xposition - a;
	mes[1] = Zposition - b;
	mes[2] = amplitude;//返回临时的坐标与振幅
	//cout << "finalx=" << mes[0]<< " " << "finalz=" << mes[1]<< "  :" << mes[2] << endl;
	
	return mes;
}

void operators::perturb_veclocity(float low, float high,float step ,int left, int right, int up, int down,Array2D<float> velocity,Array2D<float> record,int dx,int dz)
{
	
	float rate = 1;
	
	
	vector<float> tmp_mes(3);
	tmp_mes=perturb(left,right,up,down,velocity,record,dx,dz);

	
	vector<ray*> tmp_ray(n);//射线组的中间变量
	for (int i = 0; i < n; i++)
	{
		tmp_ray[i] = new ray(*ray_group[i]);
	}//扰动前的射线束拷贝

	
	vector<float> mes;
	

	Array2D<float> tmp_velocity;
	tmp_velocity.Resize(nx,nz);
	tmp_velocity.Zero();
	tmp_velocity=velocity;
	

	vector<float> offset(nx);
	for (float r = low; r <high; r += step)
	{
		for (int i = 0; i < n; i++)
		{
			ray_group[i]->time= tmp_ray[i]->time  ;
		}
		//cout << low << endl;
		//cout << high << endl;
		//cout << step << endl;
		//cout << r << endl;
		//rate = r;
		//cout << "velocity=" << r *velocity[20][20]<< endl;
		//cout << "r=" << r;
		
		for (int i = 0; i < nx; i++)
		{
			for (int j = 0; j < nz; j++)
			{
				tmp_velocity[i][j] = velocity[i][j] * r;
				
			}
		}
		
		for (int i = 0; i < n; i++)
		{
			ray_group[i]->time = (tmp_ray[i]->time) / r;
			//cout << "tmp: " << tmp_ray[i]->time << endl;
			//cout << "time: " << ray_group[i]->time << endl;
		}
		/*
		for (int i = 0; i < n; ++i)
		{
			cout << "x=" << ray_group[i]->getx() << " z=" << ray_group[i]->getz() << " t=" << ray_group[i]->time << endl;
		}*/

		mes=perturb(left,right,up,down,tmp_velocity,record,dx,dz);
		
		//cout << "xposition=" << mes[0] << "  zposition=" << mes[1] << "  " << mes[2] << endl;

		
		if (mes[2] > tmp_mes[2])
		{
			tmp_mes[0]= mes[0];
			 tmp_mes[1]= mes[1] ;
			tmp_mes[2] =mes[2];
			//cout << "xposition=" << tmp_mes[0] << "  zposition=" << tmp_mes[1] << "  " << tmp_mes[2] << endl;
			rate = r;
		}
		//cout << endl;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	Xposition=tmp_mes[0];
    Zposition=tmp_mes[1];
    amp=tmp_mes[2];
    velocity_rate=rate;
	cout << "xposition=" << tmp_mes[0] << "  zposition=" << tmp_mes[1] <<"  "<< tmp_mes[2]<<"  rate="<<rate<< endl;
    cout << "velocity=" << rate*velocity[20][20]<< endl;
	for (int i = 0; i < n; i++)
	{
		ray_group[i] = tmp_ray[i];
	}
	//cout << "over" << endl;
	
	for (int i = 0; i < n; i++)
	{
		delete(tmp_ray[i]);

	}

}
Array2D<float> operators::change_velocity(Array2D<float> velocity)
{
    /*
	vector<point> first_line;
	first_line=ray_group[0]->insert_line();
	int n_f=first_line.size();
	if(first_line[n_f-1].z>detector)
	{
		first_line.pop_back();
	}
	else if(first_line[n_f-1].z<detector)
	{
		point tmp(first_line[n_f-1].x,detector);
		first_line.push_back(tmp);
	}
	n_f=first_line.size();

	vector<point> last_line;
	last_line=ray_group[n-1]->insert_line();
	int n_l=last_line.size();
	if(last_line[n_l-1].z>detector)
	{
		last_line.pop_back();
	}
	else if(last_line[n_l-1].z<detector)
	{
		point tmp2(last_line[n_l-1].x,detector);
		last_line.push_back(tmp2);
	}
	n_l=last_line.size();//两射线在检波面处对齐
	if(n_l!=n_f)
	{
		cout<<"fuck tow line"<<endl;
	}//检查对齐结果

	 Array2D<float> ratio;
	 ratio.Resize(nx, nz);
	 ratio.Zero();
	for(int i=1;i<n_f;++i)
	{
		int a=first_line[i].x;
		int b=last_line[i].z;
		for(int j=a;j<b+1;++j)
		{
			velocity[j][i]=velocity[j][i]*velocity_rate;
			ratio[j][i]=1;
		}
	}
	vector<Array2D<float>> a(2);
	a[0]=velocity;
	a[1]=ratio;

	return a;*/

    vector<int> l_1(Zposition);
    vector<int> l_2(Zposition);
    float grad1=Xposition/Zposition;
    float grad2=(nx-Xposition)/Zposition;
    
    for(int i=0;i<Zposition;++i)
    {
        l_1[i]=to_int(i*grad1);
        l_2[i]=to_int(nx-i*grad2);
    }
    //cout<<"////////////////"<<endl;
    //cout<<Zposition<<endl;

    for(int i=0;i<Zposition;i++)
    {
        for(int j=l_1[i];j<l_2[i];j++)
        {
            velocity[j][i]=velocity[j][i]*velocity_rate;
        }
    }
    //cout<<"change is over"<<endl;
    return velocity;
    
}
/*
Array2D<float> operators::change_layervelocity(Array2D<float> velocity, vector<float> up)
{


    vector<int> l_1(Zposition);
    vector<int> l_2(Zposition);
    float grad1=Xposition/Zposition;
    float grad2=(nx-Xposition)/Zposition;
    
    for(int i=0;i<Zposition;++i)
    {
        l_1[i]=to_int(i*grad1);
        l_2[i]=to_int(nx-i*grad2);
    }
    //cout<<"////////////////"<<endl;
    //cout<<Zposition<<endl;

    for(int i=l_1[0];i<l_2[0];i++)
    {
        if(i<Xposition)
        {
            for(int j=up[i];j<l_1[i];j++)
            {
                velocity[i][j]=velocity[j][i]*velocity_rate;
            }
        }
        else if(i>Xposition)
        {
            for(int j=up[i];j<l_1[i];j++)
            {
                velocity[i][j]=velocity[j][i]*velocity_rate;
            }
        }
    }
    //cout<<"change is over"<<endl;
    return velocity;
    
}*/


void operators::sort_ray()
{
	if (n != (int)ray_group.size())
	{
		cout << "n is erro" << endl;
		n = ray_group.size();
	}
	for(int i=0;i< n-1;++i)
		for (int j = 0; j < n - i-1; ++j)
		{
			
			if (ray_group[j]->angle > ray_group[j + 1]->angle)
			{
				std::swap(ray_group[j], ray_group[j+1]);
			}
		}
}




