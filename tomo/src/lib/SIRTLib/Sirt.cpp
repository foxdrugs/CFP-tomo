// cpp file


/***********************************************

      FileName: Sirt.cpp

        Author: stj
   Description: ---
 First  Create: 2020-07-13 10:31:13
 Last Modified: 2020-07-25 22:25:14

***********************************************/

#include"Sirt.h"
#include<vector>
using std::vector;

Array2D<float> sirt(vector<Array2D<float>> velocity, int nx, int nz)
{
    int n=velocity.size();
    Array2D<float> result;
    result.Resize(nx,nz);
    result.Zero();
    
    for(int i=0;i<n;++i)
    {
        for(int j=0;j<nx;++j)
        {
            for(int k=0;k<nz;++k)
            {
                result[j][k]+=velocity[i][j][k];
            }
        }
    }
    for(int j=0;j<nx;++j)
    {
        for(int k=0;k<nz;++k)
        {
            result[j][k]=result[j][k]/n;
        }
    }
    

    return result;
}

Array2D<float> constraint(Array2D<float> velocity,vector<float> l_up, vector<float> l_down,int nx,int nz)
{
    float sum=0;
    int num=0;

    Array2D<float> result;
    result.Resize(nx,nz);
    result.Zero();
    result=velocity;
    for(int i=0;i<nx;i++)
    {
        for(int j=l_up[i];j<l_down[i];++j)
        {
            sum+=velocity[i][j];
            num++;
        }
    }

    float v=sum/num;
    for(int i=0;i<nx;i++)
    {
        for(int j=l_up[i];j<l_down[i];++j)
        {
            result[i][j]=v;    
        }
    }

    
    return result;

}


vector<float> interpola(int nx, vector<point> layer)
{
    vector<float> result(nx);

    int n=layer.size();
    float grad[n-1];
    for(int i=1;i<n;++i)
    {
        if(layer[i].x==layer[i-1].x)
        {
            grad[i-1]=0;
        }
        else
        {
        grad[i-1]=(layer[i].z-layer[i-1].z)/(layer[i].x-layer[i-1].x);
        }
        int num=0;
        for(int j=layer[i-1].x;j<layer[i].x;++j)
        {
            result[j]=layer[i-1].z+num*grad[i-1];
            num++;
        }
    }
    
    float g1=0;
    float g2=0;
    if(n>6)
    {
        g1=(grad[0]+grad[1]+grad[3])/3;
        g2=(grad[n-2]+grad[n-3]+grad[n-4])/3;
    }
    else
    {
        g1=grad[0];
        g2=grad[n-2];
    }
    int num=0;
    for(int j=layer[0].x;j>0;j=j-1)
    {
        result[j]=layer[0].z-g1*num;
        num++;
    }
    
    for(int j=layer[n-1].x;j<nx-1;++j)
    {
        result[j]=layer[n-1].z+g2*(j-layer[n-1].x);
    }
    result[0]=result[1];
    result[nx-1]=result[nx-2];
    return result;

}

Array2D<float> smooth(Array2D<float> velocity, int nx, int nz,int n)
{
    Array2D<float> tmp;
    tmp.Resize(nx,nz);
    tmp.Zero();
    for(int k=0;k<n;k++)
    {
        for(int i=0;i<nx;++i)
        {
            for(int j=2;j<nz-2;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i][j-1]+velocity[i][j+1]+velocity[i][j-2]+velocity[i][j+2])/5;
            }
        }

        for(int i=0;i<nx;++i)
        {
            for(int j=0;j<2;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i][j+1]+velocity[i][j+2])/3;
            }
        }

        for(int i=0;i<nx;++i)
        {
            for(int j=nz-2;j<nz;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i][j-1]+velocity[i][j-2])/3;
            }
        }//border smooth  z demension

        velocity=tmp;

        for(int i=2;i<nx-2;++i)
        {
            for(int j=0;j<nz;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i-1][j]+velocity[i+1][j]+velocity[i-2][j]+velocity[i+2][j])/5;
            }
        }

        for(int i=0;i<2;++i)
        {
            for(int j=0;j<nz;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i+1][j]+velocity[i+2][j])/3;
            }
        }

        for(int i=nx-2;i<nx;++i)
        {
            for(int j=0;j<nz;++j)
            {
                tmp[i][j]=(velocity[i][j]+velocity[i-1][j]+velocity[i-2][j])/3;
            }
        }//border  x demension

        velocity=tmp;

    }
    return velocity;

}
