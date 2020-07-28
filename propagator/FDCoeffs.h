#ifndef FINITE_DIFFERENCE_COEFFICIENTS_H
#define FINITE_DIFFERENCE_COEFFICIENTS_H
// Yezheng (Jim) HU
// History:
// Date: May 7 2015 first version
// July 20 2015 Add interface for coeffs to be pointer and changed return type to bool
//
#include<vector>

// Routines below only handles center point finite difference

// Routine: FirstdDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side exclude the central point
// coeffs[0]: weight for x0+h
// coeffs[1]: weight for point x0+2h
// ....
bool FirstDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs);
bool FirstDerivativeCoefficients(const int & npoint, float *&coeffs);

// Routine: SecondDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side include the central point
// coeffs[0]: weight for x0
// coeffs[1]: weight for point x0+h
// ....
bool SecondDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs);
bool SecondDerivativeCoefficients(const int & npoint, float *&coeffs);

// Routine: ThirddDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side exclude the central point
// coeffs[0]: weight for x0+h
// coeffs[1]: weight for point x0+2h
// ....
bool ThirdDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs);
bool ThirdDerivativeCoefficients(const int & npoint, float *&coeffs);

// Routine: FourthDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side include the central point
// coeffs[0]: weight for x0
// coeffs[1]: weight for point x0+h
// ....
bool FourthDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs);
bool FourthDerivativeCoefficients(const int & npoint, float *&coeffs);

void CenteredFirstDerivativeCoefficients(const int &spaceOrder, const float &delta, std::vector<float> &coeffs);
void CenteredSecondDerivativeCoefficients(const int &spaceOrder, const float &delta, std::vector<float> &coeffs);

#endif
