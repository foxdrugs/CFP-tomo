// cpp file


/***********************************************

      Filename: Ricker.cpp

        Author: ShenTianJing
   Description: ---
 First  Create: 2018-10-27 14:54:51
 Last Modified: 2018-10-28 19:02:54

***********************************************/


#include"Ricker.h"
#include<cstdio>
#include<cmath>

#define PI 3.1415926

Ricker::Ricker(int Nt, float Dt, float Fpeak)
{
	N_set(Nt);
	D_set(Dt);
	F_set(Fpeak);

	ricker_wavelet(Nt, Dt, Fpeak);


}


void Ricker::ricker_wavelet(int nt, float dt, float fpeak)
{
	int it;
	float   t1, t0;

	t0 = 1.0 / fpeak;
	m_wavelet.Resize(nt);

	for (it = 0; it < nt; it++) {
		t1 = it * dt;
		m_wavelet[it] = exp(-PI * PI*fpeak*fpeak*(t1 - t0)*(t1 - t0))*(1.0 - 2.*PI*PI*fpeak*fpeak*(t1 - t0)*(t1 - t0));
	}
}



