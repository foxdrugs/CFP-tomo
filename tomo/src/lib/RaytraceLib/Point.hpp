// h file


/***********************************************

      FileName: Point.h

        Author: stj
   Description: ---
 First  Create: 2019-12-29 16:25:31
 Last Modified: 2019-12-29 21:13:10

***********************************************/

#ifndef POINT_H 
#define POINT_H
#include<iostream>
using std::cout;
using std::endl;
class point
{
    public:
        float x,z;
    point(float a,float b)
    {
        x=a;
        z=b;
    }
	 point(const point &tmp)
	{
		x = tmp.x;
		z = tmp.z;
	}
    point()
    {
        x=0;
        z=0;
    }
    point(float a)
    {
        x=a;
        z=a;
    }
    void voluation(float a, float b)
    {
        x=a;
        z=b;
    }
    void display() const
    {
        cout<<"x="<<x<<" "<<"z="<<z<<endl;
    }
};



#endif 
