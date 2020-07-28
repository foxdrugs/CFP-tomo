// h file


/***********************************************

      Filename: Ricker.h

        Author: ShenTianJing
   Description: ---
 First  Create: 2018-10-27 14:54:01
 Last Modified: 2018-10-28 19:04:36

***********************************************/
#ifndef RICKER_H
#define RICKER_H
#include"Array.hpp"

class Ricker {
public:
	void N_set(int Nt) { nt = Nt; }
	void D_set(int Dt) { dt = Dt; }
	void F_set(int Fpeak) { fpeak = Fpeak; }

	Array1D<float>& RickerWavelet() { return m_wavelet; }


public:
	Ricker(int Nt, float Dt, float Fpeak);
	~Ricker() { ; }

public:
	void ricker_wavelet(int nt, float dt, float fpeak);

private:
	int nt;
	float dt;
	float fpeak;
	Array1D<float> m_wavelet;

};
#endif






