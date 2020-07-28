// cpp file


/***********************************************

      Filename: Propagator.cpp

        Author: ShenTianJing
   Description: ---
 First  Create: 2018-10-30 17:20:45
 Last Modified: 2020-07-17 21:32:04

***********************************************/


#include"Propagator.h"
#include<iostream>
#include"CommandLineParser.h"
#include<fstream>
#include<assert.h>
#include"FDCoeffs.h"
#include<vector>
#include<omp.h>
#include"CPMLBoundary.h"

using namespace std;

Propagator::Propagator(CommandLineParser &cmlParser, int Ixsrc,int Izsrc,const string Outfilename): m_cmlParser(cmlParser)
{
     npml=16;
    ixsrc=Ixsrc;
    izsrc=Izsrc;
    outfilename=Outfilename;
	a_ratio.Resize(npml);
    b_ratio.Resize(npml);
    a_ratio.Zero();
    b_ratio.Zero();
    
	InitializeParameters();
	CreatRickerClass();
    CreatDerivativeCoefficients();
    CreatCPMLratio();
//////////cout<<"construct correct"<<endl;
}



void Propagator::InitializeParameters()
{


   
	
    string inputfilename;
    m_cmlParser.getIntValue("waveletNt", m_waveletNt);
    //////////cout<<"Ricker_t="<<m_waveletNt<<endl;
    
    m_cmlParser.getStringValue("inputfile", inputfilename);
    //////////cout<<"Ricker_t="<<m_waveletNt<<endl;

	m_cmlParser.getIntValue("nStep", m_nStep);
    //////////cout<<"nStep="<<m_nStep<<endl;

    m_cmlParser.getIntValue("order", order);
    //////////cout<<"order="<<order<<endl;


	m_cmlParser.getFloatValue("dt", m_dt);
    //////////cout<<"dt="<<m_dt<<endl;

    m_cmlParser.getFloatValue("dx", m_dx);
    //////////cout<<"dx="<<m_dx<<endl;


    

	m_cmlParser.getIntValue("nx", m_nx);
    //////////cout<<"nx="<<m_nx<<endl;

	m_cmlParser.getIntValue("nz", m_nz);
    //////////cout<<"nz="<<m_nz<<endl;

	m_cmlParser.getIntValue("peakF", m_peakF);
    //////////cout<<"peakF="<<m_peakF<<endl;

   	m_speed.Resize(m_nx, m_nz);

    const char* filename=inputfilename.data();
    int a=0;
    FILE *FP;
    FP=fopen(filename,"rb");
    a=fread(&m_speed[0][0],sizeof(float),m_nx*m_nz,FP);
    fclose(FP);
    if(a==m_nx*m_nz)
    {
        cout<<"read correct"<<endl;
        //cout<<"speed size is "<<a<<endl;
    }







}


void Propagator::CreatRickerClass()
{

	m_ricker = new Ricker(m_waveletNt, m_dt, m_peakF);

}


void Propagator::CreatDerivativeCoefficients()
{

SecondDerivativeCoefficients(order,SecondDerivative);
//for(int i=0;i<SecondDerivative.size();i++)
   // cout<<SecondDerivative[i]<<endl;

FirstDerivativeCoefficients(order,FirstDerivative);
//for(int i=0;i<FirstDerivative.size();i++)
   //cout<<FirstDerivative[i]<<endl;


}


void Propagator::FireOneShot()
{   
    bclength=npml+order;
	cout<<"bclength="<<bclength<<endl;
    const char* filename=outfilename.data();

    Array2D<float> L_X_psi;
    Array2D<float> R_X_psi;  
    Array2D<float> L_X_eta;  
    Array2D<float> R_X_eta;

    Array2D<float> D_Z_psi;
    Array2D<float> UP_Z_psi;  
    Array2D<float> D_Z_eta;  
    Array2D<float> UP_Z_eta;

    L_X_psi.Resize(bclength,m_nz);  
    R_X_psi.Resize(bclength,m_nz);  
    L_X_eta.Resize(bclength,m_nz);  
    R_X_eta.Resize(bclength,m_nz);
    offset.Resize(m_nx-bclength-bclength,m_nStep);

    UP_Z_psi.Resize(m_nx,bclength);  
    D_Z_psi.Resize(m_nx,bclength);  
    UP_Z_eta.Resize(m_nx,bclength);  
    D_Z_eta.Resize(m_nx,bclength); 

    int detector=bclength;  
    //m_cmlParser.getIntValue("detector", detector);
    cout<<"detector="<<detector<<endl;
    offset.Zero();
    L_X_psi.Zero();
    R_X_eta.Zero();
    R_X_psi.Zero();
    L_X_eta.Zero();

    D_Z_psi.Zero();
    UP_Z_eta.Zero();
    D_Z_psi.Zero();
    UP_Z_eta.Zero();

  
    m_befor.Resize(m_nx,m_nz);
    m_now.Resize(m_nx, m_nz);
	m_nex.Resize(m_nx, m_nz);
    m_befor.Zero();
    m_now.Zero();
    m_nex.Zero();

	ixsrc += bclength;
	izsrc += bclength;
   cout<<"xsrc="<<ixsrc<<endl;
   cout<<"zrc="<<izsrc<<endl;

    cout<<filename<<endl;
    string file="offset.bin";
///////////////////////////////////////////////////////////////////////////
//change the model 
	/*
    if(outfilename!=file)
    {
        cout<<"change the velocity"<<endl;
        for(int i=0;i<m_nx;i++)
            for(int j=izsrc-1;j<m_nz;j++)
            {
                m_speed[i][j]=m_speed[ixsrc][izsrc-2];
            }
    }
    else
    {
        cout<<"no change"<<endl;
    }*/
//this part can change according to the actual situation
///////////////////////////////////////////////////////////////////////////

int secondderivative_numbers=SecondDerivative.size();
int number=secondderivative_numbers-1;

////////////cout<<"secondderivative_numbers="<<secondderivative_numbers<<endl;
int firstderivative_numbers=FirstDerivative.size();

////////////cout<<"firstderivative_numbers="<<firstderivative_numbers<<endl;

float dt2;
dt2=m_dt*m_dt;

float dx;
dx=1.0/(m_dx*m_dx) ;

int distance;
distance=m_nx-order-npml;
int aaa=m_nz-m_nx;
//////////cout<<"distance"<<distance<<endl;


//////////cout<<"bclength="<<bclength<<endl;
//////////cout<<"oeder="<<order<<endl;


for (int istep = 0; istep < m_nStep; istep++)
{
    
        
    if( istep< m_waveletNt)
    {
        m_now[ixsrc][izsrc] = m_ricker->RickerWavelet()[istep];
		//cout << m_ricker->RickerWavelet()[istep] << "  " << istep << endl;
    }



        #pragma omp parallel for collapse(2)  
        
        for (int i = number; i < m_nx - number; i++)                                                                   
		    for (int j = number; j < m_nz - number; j++)
		     {
			
                float P_X_2_derivative = SecondDerivative[0] * m_now[i][j];
                float P_Z_2_derivative = SecondDerivative[0] * m_now[i][j];
				
                
                for (int k = 1; k < secondderivative_numbers; k++)
			        {
			     	    P_X_2_derivative= P_X_2_derivative+SecondDerivative[k] * (m_now[i][j - k] + m_now[i][j + k]);
			            P_Z_2_derivative=P_Z_2_derivative +SecondDerivative[k] * (m_now[i - k][j] + m_now[i + k][j]);			     
                    }
                
				//  1:calculate second derivative of u int x direction; part1
				//  2:calculate second derivative of u int z direction; part2
                    m_nex[i][j] = 2 * m_now[i][j] - m_befor[i][j]+(P_X_2_derivative +P_Z_2_derivative )*m_speed[i][j]*m_speed[i][j] *dt2*dx;
                   
			    
            }   //calculate middle field

      #pragma omp parallel for collapse(2)  
     for (int i = number; i < bclength-number; i++)
         {                                                                   
		    for (int j = 0; j < m_nz; j++)
		     {
                float L_Boundary_X_1_derivative=0;
                float R_Boundary_X_1_derivative=0;

                for (int k = 1; k < firstderivative_numbers+1; k++)
			        {
                        
                        L_Boundary_X_1_derivative += FirstDerivative[k-1] * (m_now[i+k][j] - m_now[i-k][j]);
                        R_Boundary_X_1_derivative += FirstDerivative[k-1] * (m_now[i+distance+k][j] - m_now[i+distance-k][j]);
				     }    
                            L_X_psi[i][j]= L_X_psi[i][j]*a_ratio[i-number]+b_ratio[i-number]*L_Boundary_X_1_derivative/m_dx;  //calculate psi field 
                            R_X_psi[i][j]= R_X_psi[i][j]*a_ratio[npml-1-(i-number)]+b_ratio[npml-1-(i-number)]*R_Boundary_X_1_derivative/m_dx;     
            }
        }
#pragma omp parallel for collapse(2)
    for (int i = number; i < bclength-number; i++)
    {                                                                   
	    for (int j = 0; j < m_nz; j++)
		{
            float L_part1 = 0.f;
            float L_part2 = 0.f;
            float R_part1 = 0.f;
            float R_part2 = 0.f; 
         
            L_part2 = SecondDerivative[0] * m_now[i][j];//normsl first and second derivative
            R_part2 = SecondDerivative[0] * m_now[i+distance][j];

            for (int k = 1; k < secondderivative_numbers; k++)  
            {   
               L_part2 += SecondDerivative[k] * (m_now[i-k][j] + m_now[i+k][j]);
               R_part2 += SecondDerivative[k] * (m_now[i+distance-k][j] + m_now[i+distance+k][j]);
            }

            for (int k=1;k<firstderivative_numbers+1;k++)
            {    
                L_part1 += FirstDerivative[k-1]*(L_X_psi[i+k][j]-L_X_psi[i-k][j]);  //ratio first derivative  
                R_part1 += FirstDerivative[k-1]*(R_X_psi[i+k][j]-R_X_psi[i-k][j]);
            }

            L_X_eta[i][j] = a_ratio[i-number]*L_X_eta[i][j] + b_ratio[i-number]*(L_part1/m_dx + L_part2*dx);//calculate eta field
            R_X_eta[i][j] = a_ratio[npml-1-(i-number)]*R_X_eta[i][j] + b_ratio[npml-1-(i-number)]*(R_part1/m_dx + R_part2*dx);


            m_nex[i][j] += (L_part1/m_dx + L_X_eta[i][j])*dt2*m_speed[i][j]*m_speed[i][j];
            m_nex[i+distance][j] += (R_part1/m_dx + R_X_eta[i][j])*dt2*m_speed[i+distance][j]*m_speed[i+distance][j];

        }
    }                          //save x direction boundary  

 #pragma omp parallel for collapse(2)
for (int i = 0; i < m_nx; i++)
         {                                                                   
		    for (int  j = number; j < bclength-number; j++)
		     {         
                float UP_Boundary_Z_1_derivative=0;
                float D_Boundary_Z_1_derivative=0;
                
                for (int k = 1; k < firstderivative_numbers+1; k++)
			        {
                        
                        UP_Boundary_Z_1_derivative += FirstDerivative[k-1] * (m_now[i][j+k] - m_now[i][j-k]);     
                        D_Boundary_Z_1_derivative += FirstDerivative[k-1] * (m_now[i][j+distance+aaa+k] - m_now[i][j+distance+aaa-k]);
                    }    
                
                            UP_Z_psi[i][j]= UP_Z_psi[i][j]*a_ratio[j-number]+b_ratio[j-number]*UP_Boundary_Z_1_derivative/m_dx;  //calculate psi field 
                            D_Z_psi[i][j]= D_Z_psi[i][j]*a_ratio[npml-1-(j-number)]+b_ratio[npml-1-(j-number)]*D_Boundary_Z_1_derivative/m_dx;
            }
        }
#pragma omp parallel for collapse(2)
    for (int i = 0; i < m_nx; i++)   
    {                                                                   
	    for (int j=number; j < bclength-number; j++)
		{
            float UP_part1 = 0.f;
            float UP_part2 = 0.f;
            float D_part1 = 0.f;
            float D_part2 = 0.f; 

            UP_part2 = SecondDerivative[0] * m_now[i][j];//normsl first and second derivative
            D_part2 = SecondDerivative[0] * m_now[i][j+distance+aaa];

            for (int k = 1; k < secondderivative_numbers; k++)  
            {   
               UP_part2 += SecondDerivative[k] * (m_now[i][j-k] + m_now[i][j+k]);
               D_part2 += SecondDerivative[k] * (m_now[i][j+distance+aaa-k] + m_now[i][j+distance+aaa+k]);
            }

            for (int k=1;k<firstderivative_numbers+1;k++)
            {    
                UP_part1 += FirstDerivative[k-1]*(UP_Z_psi[i][j+k]-UP_Z_psi[i][j-k]);  //ratio first derivative  
                D_part1 += FirstDerivative[k-1]*(D_Z_psi[i][j+k]-D_Z_psi[i][j-k]);
            }

            UP_Z_eta[i][j] = a_ratio[j-number]*UP_Z_eta[i][j] + b_ratio[j-number]*(UP_part1/m_dx + UP_part2*dx);//calculate eta field
            D_Z_eta[i][j] = a_ratio[npml-1-(j-number)]*D_Z_eta[i][j] + b_ratio[npml-1-(j-number)]*(D_part1/m_dx + D_part2*dx);
            
            m_nex[i][j] += (UP_part1/m_dx + UP_Z_eta[i][j])*dt2*m_speed[i][j]*m_speed[i][j];
            m_nex[i][j+distance+aaa] += (D_part1/m_dx + D_Z_eta[i][j])*dt2*m_speed[i][j+distance+aaa]*m_speed[i][j+distance+aaa];
        }
    }                          //save z direction boundary 
	

m_befor=m_now;
m_now=m_nex;                  //exchange the time
for(int i=bclength;i<m_nx-bclength;i++)
    offset[i-bclength][istep]=m_now[i][detector];
               
}

//cout<<"bclength="<<bclength<<endl;
//cout<<"m_nz="<<m_nz<<endl;
                FILE *fp;
                fp=fopen("Propagator.bin","wb");
                fwrite(&(m_now[0][0]),sizeof(float),m_nx*m_nz, fp);
                fclose(fp);
cout<<"n1(Propagator.bin)="<<m_nz<<endl;
cout<<"n2(Propagator.bin)="<<m_nx<<endl;
                FILE *Fp;
                Fp=fopen(filename,"wb");
                fwrite(&(offset[0][0]),sizeof(float),(m_nx-bclength-bclength)*m_nStep, fp);
                fclose(Fp);
cout<<"n1("<<filename<<")="<<m_nStep<<endl;
cout<<"n2("<<filename<<")="<<m_nx-bclength-bclength<<endl;

               

}




void Propagator::CreatCPMLratio()
{   

    
    CPMLBoundary ratio(npml, 5, m_dt, 25, 3000, 1e-5);
    for (int i=0;i < npml;i++)
    {
        a_ratio[i]=ratio.A()[i];
        b_ratio[i]=ratio.B()[i];
////////////cout<<"a="<<a_ratio[i]<<endl;
////////////cout<<"b="<<b_ratio[i]<<endl;

    }




}
