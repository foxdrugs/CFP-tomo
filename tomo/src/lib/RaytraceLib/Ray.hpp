// h file


/***********************************************

      FileName: Ray.h

        Author: stj
   Description: ---
 First  Create: 2019-12-27 13:53:40
 Last Modified: 2019-12-31 15:10:01

***********************************************/

#ifndef RAY_H 
#define RAY_H
#include<Point.hpp>
#include<iostream>
#include<vector>
#include<cmath>
using std::vector;
using std::cout;
using std::endl;
using std::swap;
class ray
{
    public:
        float Xposition, Zposition;
        point star;
        float time;
        float angle;
        vector<point> path;

        ray(float xposition, float zposition,float Angle)
        {
            Xposition=xposition;
            Zposition=zposition;
            star.voluation(xposition,zposition);
            path.resize(1);
            path[0]=star;
            angle=Angle;
            time=0;
        }

        ray(const ray &tmp)
        {
            Xposition=tmp.Xposition;
            Zposition=tmp.Zposition;
            star.voluation(tmp.Xposition,tmp.Zposition);
            time=tmp.time;
            angle=tmp.angle;
			path = tmp.path;
           // swap(path,tmp.path);
            
        }
		
        void print()
        {
            int n=path.size();
            for(int i=0;i<n;i++)
            {
                cout<<"x="<<path[i].x<<" "<<"z="<<path[i].z<<endl;
            }
        }

        void pushpath(float x,float z)
        {
            point tmp(x,z);
            path.push_back(tmp);

        }

        vector<point> getpath()
        {
            return path;
        }

        float gettime(float t)
        {
			time = t;
            return t;
        }

        float getx()
        {
            return path[(path.size()-1)].x;
        }
		float getz()
		{
			return path[(path.size()-1)].z;
		}
	    
        vector<point> insert_line()
		{
			int n=path.size();
			vector<point> line;
			point a(Xposition,Zposition);
			line.push_back(a);
			for(int i=0;i<n-1;i++)
			{
				float grad=(path[n+1].x-path[n].x)/(path[n+1].z-path[n].z);
				int num=floor(path[n+1].z-path[n].z);
				int fz=ceil(path[n].z);
				float dz=fz-path[n].z;
				float fx=Xposition+dz*grad;
				point first(fx,fz);
				line.push_back(first);

				for(int j=1;j<num;++j)
				{
					point c(fx+j*grad,fz+j);
					line.push_back(c);
				}
			}
			return line;
			
		}
		
};


#endif 



