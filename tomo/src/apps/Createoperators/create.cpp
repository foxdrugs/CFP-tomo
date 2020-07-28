// cpp file


/***********************************************

      FileName: create.cpp

        Author: stj
   Description: ---
 First  Create: 2020-01-03 14:52:16
 Last Modified: 2020-07-26 16:32:41

***********************************************/
#include<iostream>
#include<omp.h>
#include "Raytrace.h"
#include"Array.hpp"
#include"wavefront.h"
#include"mpi.h"
#include"operator.h"
#include"Sirt.h"
#include <sstream>
#include<fstream>

using std::cout;
using std::endl;


void creat_op(int nx,int nz,int dx, int dz,float dt,Array2D<float> velocity,operators &test)
{	
	wavefront now(dx, dz, nx, nz);
	wavefront nex(dx, dz, nx, nz);
	Array2D<float> velocityX;
	velocityX.Resize(nx, nz);
	velocityX.Zero();
	Array2D<float> velocityZ;
	velocityZ.Resize(nx, nz);
	velocityZ.Zero();
	velocityX = DerivateX(velocity, nx, nz, dx);
	velocityZ = DerivateZ(velocity, nx, nz, dz);
    
	//#pragma omp parallel for
	for (int i = 10 + 180; i < 170 + 180; i=i+2)
	{

		ray path(test.Xposition,test.Zposition, i);
		//cout << "角度为" << 360 - i << "射线追踪" << endl;
		Raytrace(path, 0, now, nex, dx, dz, nx, nz, dt, velocity, velocityX, velocityZ);
		if ((path.time-0)>0.001)
		{
			test.push_back(path);
                //cout<<path.time<<endl;
		}

    }
	test.sort_ray();

}

vector<Array2D<float>> to_vec(float*a,int nx,int nz,int n)
{
    vector<Array2D<float>> result;
    int num=0;
    Array2D<float> tmp;
    for(int k=0;k<n;k++)
    {
        tmp.Resize(nx,nz);
        tmp.Zero();
        for(int x=0;x<nx;x++)
            for(int y=0;y<nz;y++)
            {
                tmp[x][y]=a[num];
                num++;
            }

        result.push_back(tmp);
    
    }
    return result;
}

vector<point> to_point(float a[][2], int n)
{
    vector<point> result;
    for(int i=0;i<n;i++)
    {
        point tmp(a[i][0],a[i][1]);
        result.push_back(tmp);
    }
    return result;
}

vector<point> read_node(const char* filename, int n)
{
    int x=0;
    int y=0;
    vector<point> result;
    ifstream infile;
    infile.open(filename,ios::in);
    if(!infile.is_open())
        cout<<"open file failure"<<endl;
    for(int i=0;i<n;i++)
    {
        infile>>x>>y;
        point tmp(x,y);
        result.push_back(tmp);
    }
    infile.close();
    return result;
}


Array2D<float> read_record(point p, int nx, int nt)
{
    Array2D<float> result;
    result.Resize(nx,nt);
    ostringstream name;
    name<<"/home/tshen/CFP-self/CFP-tomo/infile/"<<"operator"<<"("<<p.x<<","<<p.z<<")"<<".bin";
    string Filename=name.str();
    const char* filename=Filename.data();
    cout<<filename<<endl;
    FILE* f;
	f = fopen(filename, "rb");
	fread(&(result[0][0]), sizeof(float), nx * nt, f);
	fclose(f);
    return result;

}

int _main(int argc, char** argv)
{
    int rank,size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	float nx = 300, nz = 600, dx = 5, dz = 5;
	int nt = 2000;
    float dt = 0.001;
    int n=30;
    int n1=15;

	Array2D<float> velocity;
	velocity.Resize(nx, nz);
	velocity.Zero();
	Array2D<float> record;

	FILE* s;
	s = fopen("smooth_velocity.bin", "rb");
	fread(&(velocity[0][0]), sizeof(float), nz * nx, s);
	fclose(s);

    vector<point> focus=read_node("/home/tshen/CFP-self/CFP-tomo/infile/node.txt",n);
    record=read_record(focus[rank],nx,nt);
    float P[n][2];
    for(int i=0;i<n;++i)
    {
        P[i][0]=focus[i].x;
        P[i][1]=focus[i].z-5;
    }
    float sum_amp[n];
    float local_ampitude=0;
    float global_sumamp=0;
    cout<<"n="<<size<<endl;

    float *velocity_group= new float[int(n*nx*nz)]; 
    cout<<"get in loop"<<endl;
/////////////////////////////////////////////////////////////////////
for(int loop=0;loop<30;loop++)
{
	
    cout<<"x="<<P[rank][0]<<"  "<<"z="<<P[rank][1]<<"  rank="<<rank<<endl;

	operators test(P[rank][0], P[rank][1], int(nx), int(nz),dt);
    cout<<"Raytrace"<<endl;
	creat_op(nx,nz,dx,dz,dt,velocity,test);
    int left=0;
    int right=0;
    if(rank==0||rank==n1)
    {
        left=P[rank][0]-1;
        right=P[rank+1][0]-P[rank][0]-1;
    }
    else if(rank==(n1-1)||(rank==n-1))
    {
        left=P[rank][0]-P[rank-1][0]-1;
        right=nx-P[rank][0]-1;
    }
    else
    {
        left=P[rank][0]-P[rank-1][0]-1;
        right=P[rank+1][0]-P[rank][0]-1;
    }

    cout<<"pertub"<<"  rank="<<rank<<endl;   
	test.perturb_veclocity(0.5,2,0.2,left, right, 5, 5, velocity, record, dx, dz);

    Array2D<float> v_change;
    v_change.Resize(nx,nz);
    v_change.Zero();
    v_change=test.change_velocity(velocity);
    
    
    float local_p[2];
    local_p[0]=test.Xposition;
    local_p[1]=test.Zposition;
    local_ampitude=test.amp;

    float tran[(int)nx][(int)nz];
    for(int i=0;i<nx;i++)
        for(int j=0;j<nz;j++)
        {
            tran[i][j]=v_change[i][j];
        }
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Allgather((void*)&local_p,2,MPI_INT,(void *)&P,2,MPI_INT,MPI_COMM_WORLD);
    MPI_Allgather((void*)&local_ampitude,1,MPI_FLOAT,&sum_amp,1,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Allgather(tran,nx*nz,MPI_FLOAT,velocity_group,nx*nz,MPI_FLOAT,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    vector<Array2D<float>> out(n);
    out=to_vec(velocity_group,nx,nz,n);

    float p1[n1][2];
    float p2[n-n1][2];
    for(int i=0;i<n1;i++)
    {
        p1[i][0]=P[i][0];
        p1[i][1]=P[i][1];
        p2[i][0]=P[n1+i][0];
        p2[i][1]=P[n1+i][1];
    }

    vector<float> l_1(nx);
    vector<float> l_2(nx);
    vector<float> l_3(nx);
    vector<float> l_4(nx);
    l_2=interpola(nx,to_point(p1,n1));
    l_3=interpola(nx,to_point(p2,n1));
    if(rank==0)
    {
        for(int i=0;i<n1;i++)
        {
            cout<<p1[i][0]<<"   "<<p1[i][1]<<endl;
            cout<<p2[i][0]<<"   "<<p2[i][1]<<endl;
        }
    }

    for(int i=0;i<nx;i++)
        l_4[i]=nz-1;
    cout<<"sirt"<<endl;
    Array2D<float> result;
    result.Resize(nx,nz);
    result.Zero();
    Array2D<float> sirted;
    sirted.Resize(nx,nz);
    sirted.Zero();
    sirted=sirt(out,nx,nz);
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==0)
    {
        for(int i=0;i<nx;i++)
        {
        cout<<l_2[i]<<endl;
        cout<<l_3[i]<<endl;
        }
    }
    cout<<"constraint"<<endl;
    result=constraint(sirted,l_1,l_2,nx,nz);
    result=constraint(result,l_2,l_3,nx,nz);
    result=constraint(result,l_3,l_4,nx,nz);
    cout<<"constraint finish"<<endl;
    float local_sumamp=0;
    for(int i=0;i<n;++i)
    {
            local_sumamp+=sum_amp[i];
    }
    
    //if(local_sumamp>global_sumamp)
    //{
      //  global_sumamp=local_sumamp;
       if(rank==0)
        {
            ostringstream name;
            name<<"/home/tshen/CFP-self/CFP-tomo/outfile/"<<"modelmpi"<<loop<<"a"<<local_sumamp<<".bin";
            string Filename=name.str();
            const char* filename=Filename.data();
            cout<<filename<<endl;
            FILE* F;
	        F = fopen(filename, "wb");
	        fwrite(&(result[0][0]),sizeof(float),nx*nz,F);
	        fclose(F);
            cout<<"amp="<<local_sumamp<<endl;

        }
                
    //}

    velocity=smooth(result,nx,nz,50);
}
    delete[]velocity_group;
    MPI_Finalize();
    return 0;

}

int main(int argc, char** argv)
{
	int rank,size;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

	float nx = 300, nz = 600, dx = 5, dz = 5;
	int nt = 2000;
    float dt = 0.001;
    int n=30;
    int n1=15;

	Array2D<float> velocity;
	velocity.Resize(nx, nz);
	velocity.Zero();
	Array2D<float> record;

	FILE* s;
	s = fopen("smooth_velocity.bin", "rb");
	fread(&(velocity[0][0]), sizeof(float), nz * nx, s);
	fclose(s);

    vector<point> focus=read_node("/home/tshen/CFP-self/CFP-tomo/infile/node.txt",n);
    
    float P[n][2];
    float tmp_P[n][2];
    for(int i=0;i<n;++i)
    {
        P[i][0]=focus[i].x;
        P[i][1]=focus[i].z-5;
    }
    float sum_amp[n/2];
    float local_ampitude=0;
    float global_sumamp=0;
    cout<<"n="<<size<<endl;

    float *velocity_group= new float[int(n*nx*nz)/2]; 
    cout<<"get in loop"<<endl;
/////////////////////////////////////////////////////////////////////
	for(int loop=0;loop<20;loop++)
	{
		vector<float> l_1(nx);
		vector<float> l_2(nx);
        vector<float> l_0(nx);
        float local_sumamp=0;
        Array2D<float> result;
		result.Resize(nx,nz);
		result.Zero();
		for(int k=0;k<1;k++)
		{
            record=read_record(focus[rank+k*n1],nx,nt);
			cout<<"x="<<P[rank+k*n1][0]<<"  "<<"z="<<P[rank+k*n1][1]<<"  rank="<<rank<<endl;

			operators test(P[rank+k*n1][0], P[rank+k*n1][1], int(nx), int(nz),dt);
			cout<<"Raytrace"<<endl;
			creat_op(nx,nz,dx,dz,dt,velocity,test);
			int left=0;
			int right=0;
			if(rank==0)
			{
				left=P[rank][0];
				right=P[rank+1][0]-P[rank][0];
			}
			else if(rank==(n1-1))
			{
				left=P[rank][0]-P[rank-1][0];
				right=nx-P[rank][0];
			}
			else
			{
				left=P[rank][0]-P[rank-1][0];
				right=P[rank+1][0]-P[rank][0];
			}

			cout<<"pertub"<<"  rank="<<rank<<endl;   
			test.perturb_veclocity(0.5,2,0.1,left, right, 5, 5, velocity, record, dx, dz);

			Array2D<float> v_change;
			v_change.Resize(nx,nz);
			v_change.Zero();
			v_change=test.change_velocity(velocity);
    
    
			float local_p[2];
			local_p[0]=test.Xposition;
			local_p[1]=test.Zposition;
			local_ampitude=test.amp;

			float tran[(int)nx][(int)nz];
			for(int i=0;i<nx;i++)
				for(int j=0;j<nz;j++)
				{
					tran[i][j]=v_change[i][j];
				}
    
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Allgather((void*)&local_p,2,MPI_INT,(void *)&P,2,MPI_INT,MPI_COMM_WORLD);
			MPI_Allgather((void*)&local_ampitude,1,MPI_FLOAT,&sum_amp,1,MPI_FLOAT,MPI_COMM_WORLD);
			MPI_Allgather(tran,nx*nz,MPI_FLOAT,velocity_group,nx*nz,MPI_FLOAT,MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			vector<Array2D<float>> out(n/2);
			
			out=to_vec(velocity_group,nx,nz,n/2);

			float p1[n1][2];
			
			for(int i=0;i<n1;i++)
			{
				p1[i][0]=P[i][0];
				p1[i][1]=P[i][1];
                tmp_P[i+k*n1][0]=P[i][0];
                tmp_P[i+k*n1][1]=P[i][1];
				
			}

			
			l_2=interpola(nx,to_point(p1,n1));
			
			cout<<"sirt"<<endl;
			
			Array2D<float> sirted;
			sirted.Resize(nx,nz);
			sirted.Zero();
			sirted=sirt(out,nx,nz);
			MPI_Barrier(MPI_COMM_WORLD);
			

			cout<<"constraint"<<endl;
            
			result=constraint(sirted,l_0,l_1,nx,nz);
            result=constraint(result,l_1,l_2,nx,nz);
			
			for(int i=0;i<n/2;++i)
			{
					local_sumamp+=sum_amp[i];
			}
    
			velocity=smooth(result,nx,nz,50);
            l_0=l_1;
            l_1=l_2;
                
		}
        if(rank==0)
        {
            ostringstream name;
            name<<"/home/tshen/CFP-self/CFP-tomo/outfile/"<<"test"<<loop<<"one layer"<<local_sumamp<<".bin";
            string Filename=name.str();
            const char* filename=Filename.data();
            cout<<filename<<endl;
            FILE* F;
	        F = fopen(filename, "wb");
	        fwrite(&(result[0][0]),sizeof(float),nx*nz,F);
	        fclose(F);
            cout<<"amp="<<local_sumamp<<endl;

        }
        for(int i=0;i<n;i++)
        {
            P[i][0]=tmp_P[i][0];
            P[i][1]=tmp_P[i][1];
        }

		
	}
    return 0;
}




///////////////////////////////////////////////////down is un-mpi program  ///////////////////////////////////////////////////////////////////////////////////////////////
int __main(int argc, char** argv)
{
	float nx = 300, nz = 600, dx = 5, dz = 5;
	int nt = 2000;
	Array2D<float> velocity;
	velocity.Resize(nx, nz);
	velocity.Zero();
	vector<Array2D<float>> record;
    vector<Array2D<float>> record_l;

	vector<point> focus;
    //vector<point> focus_l;

	FILE* s;
	s = fopen("smooth_velocity.bin", "rb");
	fread(&(velocity[0][0]), sizeof(float), nz * nx, s);
	fclose(s);

    float dt = 0.001;

    int num =0;
    for(int k=16;k<300;k+=20)
    {    
        //num=8;
        ostringstream name;
        name<<"operator"<<"("<<k<<","<<162+num*6<<")"<<".bin";
        string Filename=name.str();
        const char* filename=Filename.data();
        cout<<filename<<endl;
        Array2D<float> w;
        w.Resize(nx, nt);
	    w.Zero();
	    FILE* f;
	    f = fopen(filename, "rb");
	    fread(&(w[0][0]), sizeof(float), nx * nt, f);
	    fclose(f);
        record.push_back(w);

        int z=150+num*5;
        point ff(k,z);
        focus.push_back(ff);
        num++;
    }
    //num=0;
/*
     for(int k=16;k<300;k+=20)
    { 

        ostringstream name_l;
        name_l<<"operator"<<"("<<k<<","<<450<<")"<<".bin";
        string Filename_l=name_l.str();
        const char* filename_l=Filename_l.data();
        cout<<filename_l<<endl;
        Array2D<float> w_l;
        w_l.Resize(nx, nt);
	    w_l.Zero();
	    FILE* fl;
	    fl = fopen(filename_l, "rb");
	    fread(&(w_l[0][0]), sizeof(float), nx * nt, fl);
	    fclose(fl);
        record.push_back(w_l);
        


        point f_l(k,450);
        focus_l.push_back(f_l);
    }
*/

    int n1=focus.size();
    //int n2=focus_l.size();
    /*
    for(int i=0;i<n2;i++)
    {
        focus.push_back(focus_l[i]);
    }*/
    int nn=n1;//+n2;


    if(nn!=(int)focus.size())
    {
        cout<<"erro ////////////////////////////////////////"<<endl;
    }

    float step=0.2;
//cout<<n1<<" "<<n2<<endl;
for(int loop=0;loop<1;loop++)
{
    

    float amptitud=0;
    vector<Array2D<float>> velocity_group;
    vector<point> layer;
    vector<point> layer_l;
    

    for(int k=0;k<nn;k++)
    {    
        cout<<"x_begin="<<focus[k].x<<"  z_begin="<<focus[k].z<<endl;

	    /*
	    float amp = 0;
	    for (int i = 0; i < nx; ++i)
	    {
		    float a = 0;
		    int z = 0;
		    a = record[i][0];
		    for (int j = 0; j < nt-1; ++j)
		    {
			
			    if (record[i][j + 1] > a)
			    {
				    a = record[i][j + 1];
				    z = j + 1;
			    }
		    }
		    cout << "x=" << i << "  t=" << z << "  a="<<a<<endl;
		    amp += a;
	    }
	    cout << "amp=" << amp;
	    */
	    wavefront now(dx, dz, nx, nz);
	    wavefront nex(dx, dz, nx, nz);
	    Array2D<float> velocityX;
	    velocityX.Resize(nx, nz);
	    velocityX.Zero();
	    Array2D<float> velocityZ;
	    velocityZ.Resize(nx, nz);
	    velocityZ.Zero();
	    velocityX = DerivateX(velocity, nx, nz, dx);
	    velocityZ = DerivateZ(velocity, nx, nz, dz);
	    operators test(focus[k].x, focus[k].z, int(nx), int(nz),dt);
        //operators test(35, 165, int(nx), int(nz),dt);
    
        cout << "star Raytrace" <<k<< endl;
	    #pragma omp parallel for
	    for (int i = 10 + 180; i < 170 + 180; i += 2)
	    {

		    ray path(focus[k].x,focus[k].z, i);
            //ray path(35,165, i);
		    //cout << "角度为" << 360 - i << "射线追踪" << endl;
		    Raytrace(path, 0, now, nex, dx, dz, nx, nz, dt, velocity, velocityX, velocityZ);
		    if ((path.time-0)>0.001)
		    {
			    test.push_back(path);
                //cout<<path.time<<endl;
		    }
	    }

	    cout << "star pertub" << endl;
	    test.sort_ray();
        int left=0;
        int right=0;
        if(k==0)
        {
            left=5;
            right=focus[k+1].x-focus[k].x;
        }
        else if(k==(nn-1))
        {
            left=focus[k].x-focus[k-1].x;
            right=5;
        }
        else
        {
            left=focus[k].x-focus[k-1].x;
            right=focus[k+1].x-focus[k].x;
        }

        
	    test.perturb_veclocity(0.5,2,step,left, right, 5, 10, velocity, record[k], dx, dz);
        
        cout<<"pertub finish"<<endl;
        velocity_group.push_back(test.change_velocity(velocity));
        
        if(k<n1)
        {
            point tmp(test.Xposition,test.Zposition);
            layer.push_back(tmp);
            focus[k]=tmp;
        }
        else
        {
            point tmp(test.Xposition,test.Zposition);
            layer_l.push_back(tmp);
            focus[k]=tmp;
        }
        cout<<"one point finish"<<endl;
        amptitud+=test.amp;
        

    }	
    cout<<"sum_ampitude="<<amptitud<<endl;

    Array2D<float> test_velocity;
    test_velocity.Resize(nx,nz);
    test_velocity.Zero();

    vector<float> l_1(nx);
    vector<float> l_2(nx);
    //vector<float> l_3(nx);
    vector<float> l_4(nx);

    l_2=interpola(nx,layer);
    //l_3=interpola(nx,layer_l);
    cout<<"star sirt"<<endl;
    test_velocity=constraint(sirt(velocity_group,nx,nz),l_1,l_2,nx,nz);
    //test_velocity=constraint(test_velocity,l_2,l_3,nx,nz);
    cout<<"sirt finish"<<endl;



    ostringstream name;
    name<<"model"<<loop<<".bin";
    string Filename=name.str();
    const char* filename=Filename.data();
    FILE* F;
	F = fopen(filename, "wb");
	fwrite(&(test_velocity[0][0]),sizeof(float),nx*nz,F);
	fclose(F);

    velocity=smooth(test_velocity,nx,nz,50);
    cout<<"loop ="<<loop<<endl;
    cout<<endl;
    //step=step-0.05;
    cout<<"step="<<step<<endl;
}

	return 0;
}

//////////////////////////debug test//////////////////////////////////////
int test1(int argc, char** argv)
{
	float nx = 300, nz = 600, dx = 5, dz = 5;
	//int nt = 2000;
	Array2D<float> velocity;
	velocity.Resize(nx, nz);
	velocity.Zero();

    FILE* s;
	s = fopen("model15.bin", "rb");
	fread(&(velocity[0][0]), sizeof(float), nz * nx, s);
	fclose(s);

    float dt = 0.001;

    wavefront now(dx, dz, nx, nz);
	wavefront nex(dx, dz, nx, nz);
	Array2D<float> velocityX;
	velocityX.Resize(nx, nz);
	velocityX.Zero();
	Array2D<float> velocityZ;
	velocityZ.Resize(nx, nz);
	velocityZ.Zero();
	velocityX = DerivateX(velocity, nx, nz, dx);
	velocityZ = DerivateZ(velocity, nx, nz, dz);
	operators test(0, 173, int(nx), int(nz),dt);
    //cout << "star Raytrace" <<k<< endl;
	    //#pragma omp parallel for
	    for (int i = 10 + 180; i < 170 + 180; i += 2)
	    {

		    ray path(0,173, i);
            //ray path(35,165, i);
		    //cout << "角度为" << 360 - i << "射线追踪" << endl;
		    Raytrace(path, 0, now, nex, dx, dz, nx, nz, dt, velocity, velocityX, velocityZ);
		    if ((path.time-0)>0.001)
		    {
			    test.push_back(path);
                //cout<<path.time<<endl;
		    }
	    }
        return 0;

}




