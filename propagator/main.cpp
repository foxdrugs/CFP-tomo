// cpp file


/***********************************************

      Filename: main.cpp

        Author: ShenTianJing
   Description: ---
 First  Create: 2018-10-28 19:03:15
 Last Modified: 2020-07-20 21:45:42

***********************************************/

#include"Array.hpp"
#include"Propagator.h"
#include"Ricker.h"
#include"iostream"
#include"CommandLineParser.h"
//#include "mpi.h"
#include <ctime>
#include <sstream>
#include<omp.h>
using namespace std;

/*
int _main(int argc, char **argv)
{
   
    int rank,size;
    cout<<"get in the MPI"<<endl;
	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);

    CommandLineParser m_cmlParser(argc,argv);
    float ixsrc=0;
    float izsrc=0;
   
    //m_cmlParser.getFloatValue("ixsrc", ixsrc);
    //m_cmlParser.getFloatValue("izsrc", izsrc);
    
   	
    ixsrc=ixsrc+rank*10;
    ostringstream name;
    name<<"offset"<<" "<<ixsrc<<","<<izsrc<<" "<<".bin";
    string Filename=name.str();
    Propagator a(m_cmlParser,ixsrc,izsrc,Filename);
    a.FireOneShot();
cout<<rank<<endl;

	MPI_Finalize();
    


    return 0;

}*/


int main(int argc, char **argv)
{
    int ixsrc=0;
    int izsrc=0;

    CommandLineParser m_cmlParser(argc,argv);
    m_cmlParser.getIntValue("ixsrc", ixsrc);
    m_cmlParser.getIntValue("izsrc", izsrc);
    ostringstream name;
    name<<"operator"<<"("<<ixsrc<<","<<izsrc<<")"<<".bin";
    string Filename=name.str();
    cout<<"Filename="<<Filename<<endl;
    Propagator a(m_cmlParser,ixsrc,izsrc,Filename);
    a.FireOneShot();

//#pragma omp parallel for
/*    for(ixsrc=26;ixsrc<27;ixsrc+=10)
    {
        ostringstream name;
        name<<"operator"<<"("<<ixsrc<<","<<izsrc<<")"<<".bin";
        string Filename=name.str();
        Propagator a(m_cmlParser,ixsrc,izsrc,Filename);
        a.FireOneShot();
    }*/
    return 0;
}

