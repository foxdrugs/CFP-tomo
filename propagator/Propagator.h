// h file


/***********************************************

      Filename: Propagator.h

        Author: ShenTianJing
   Description: ---
 First  Create: 2018-10-27 14:51:07
 Last Modified: 2019-11-13 16:53:04

***********************************************/



#ifndef PROPAGATOR_H
#define  PROPAGATOR_H

#include"Array.hpp"
#include"CommandLineParser.h"
#include"Ricker.h"
#include<vector>
using namespace std;

class Propagator
{
public:
	Propagator(CommandLineParser &cmlParser,int Ixsrc,int Izsrc,const string Outfilename);
	~Propagator() {};

public:
	void InitializeParameters();
	void CreatRickerClass();
	void FireOneShot();
    void CreatDerivativeCoefficients();
    void CreatCPMLratio();
    int getnx(){return m_nx;}
    int getnz(){return m_nz;}
    Array2D<float> getoffset(){return offset;}
    int get_nstep(){return m_nStep;}
    int get_n2(){return m_nx-bclength-bclength;}
private:
	int m_nx;
	int m_nz;
	float m_dx;
	int m_dz;
    int order;
    int npml;
	Array2D<float> m_now;
    Array2D<float> offset;
	Array2D<float> m_nex;

	Array2D<float> m_befor;
    int bclength;
	Array2D<float> m_speed;
    Array1D<float> a_ratio;
    Array1D<float> b_ratio;

	int m_waveletNt;
	float m_dt;
    vector<float> SecondDerivative;
    vector<float> FirstDerivative;
	int m_peakF;
	float m_modelingTime;
	int m_nStep;
    CommandLineParser m_cmlParser;

    Ricker *m_ricker;
    int ixsrc;
    int izsrc;
    string outfilename;

    
};

#endif 



