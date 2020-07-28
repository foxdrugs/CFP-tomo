// h file


/***********************************************

      FileName: operator.h

        Author: stj
   Description: ---
 First  Create: 2019-12-27 15:44:35
 Last Modified: 2020-07-21 13:43:12

***********************************************/
#ifndef OPERATOR_H 
#define OPERATOR_H
#include<iostream>
#include<vector>
#include"Array.hpp"
#include"Ray.hpp"
#include<cmath>
class operators 
{
    public:
        int n;
		int nx;
		int nz;
		float dt;
		float velocity_rate;//速度扰动率
        float Xposition, Zposition;
        float amp;
        vector<ray*> ray_group;
        operators(float Xposition,float Zposition,int Nx, int Nz,float Dt);
        operators(operators &p);
        void push_back(ray tmp);
		void sort_ray();
		vector<float> perturb(int left, int right, int up, int down, Array2D<float> velocity, Array2D<float> record, int dx, int dz);
		void perturb_veclocity(float low, float high, float step, int left, int right, int up, int down, Array2D<float> velocity, Array2D<float> record, int dx, int dz);
        vector<ray*> getray_group();
		Array2D<float> change_velocity(Array2D<float> velocity);
		vector<float> insert(vector<ray*> tmp_ray);
		vector<float> insert_hermite(vector<ray*> tmp_ray);
        


};


#endif 
