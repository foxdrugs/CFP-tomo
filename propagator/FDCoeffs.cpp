// Yezheng (Jim) HU
// History:
// Date: May 7 2015 first version
// July 20 2015 Add interface for coeffs to be pointer and changed return type to bool

#include "FDCoeffs.h"
#include<iostream>
#include<vector>
#include <eigen3/Eigen/Dense>
using namespace Eigen;

// Routines below only handles center point finite difference

// Routine: FirstdDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side exclude the central point
// coeffs[0]: weight for x0+h
// coeffs[1]: weight for point x0+2h
// ....
void ComputeFirstDerivativeCoefficients(float *coeffs, const int &size)
{
    MatrixXd A(size, size);
    A.setZero();

    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            A(row, col) = pow(col+1.0, 2.0 * row + 1.0);
        }
    }
    VectorXd b(size);
    b.setZero();
    b(0) = 0.5;
    
    VectorXd x = A.fullPivLu().solve(b);
    for(int iv=0; iv<size; ++iv) coeffs[iv] = x(iv);
}

bool FirstDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs)
{
    int nv = npoint / 2;
    if (nv < 1) {
       std::cerr << "Too few points !" << std::endl;
       return false;  
    }
    coeffs.resize(nv);
    ComputeFirstDerivativeCoefficients(&coeffs[0], nv);
    return true;
}

bool FirstDerivativeCoefficients(const int & npoint, float *&coeffs)
{
   int nv = npoint / 2;
    if (nv < 1) {
       std::cerr << "Too few points !" << std::endl;
       return false;
    }
    coeffs = (float *) malloc(nv * sizeof(float));
    ComputeFirstDerivativeCoefficients(coeffs, nv);
    return true;
}

// Routine: SecondDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side include the central point
// coeffs[0]: weight for x0
// coeffs[1]: weight for point x0+h
// ....
void ComputeSecondDerivativeCoefficients(float *coeffs, const int &size)
{
    MatrixXd A(size, size);
    A.setZero();
    
    float fact = 1.f;
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (row == 0) {
               A(row,col) = (col > 0) ? 2.0 : 1.0;
            } else {
               A(row,col) = pow(1.0 * col, 2.0 * row) / fact;
            }
        }
        fact *= (2.0 * row+1.0) * (2.0 * row + 2.0);
    }
    VectorXd b(size);
    b.setZero();
    b(1) = 0.5;
    VectorXd x = A.fullPivLu().solve(b);
    for(int iv=0; iv<size; ++iv) coeffs[iv] = x(iv);
}

bool SecondDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs)
{
    int nv = npoint / 2 + 1;
    if (nv < 1) {
       std::cerr << "Too few points !" << std::endl;
       return false;  
    }
    coeffs.resize(nv);
    ComputeSecondDerivativeCoefficients(&coeffs[0], nv);
    return true;  
}

bool SecondDerivativeCoefficients(const int & npoint, float *&coeffs)
{
    int nv = npoint / 2 + 1;
    if (nv < 1) {
       std::cerr << "Too few points !" << std::endl;
       return false;  
    }
    coeffs = (float *) malloc(nv * sizeof(float));
    ComputeSecondDerivativeCoefficients(&coeffs[0], nv);
    return true;  
}

// Routine: ThirddDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side exclude the central point
// coeffs[0]: weight for x0+h
// coeffs[1]: weight for point x0+2h
// ....
void ComputeThirdDerivativeCoefficients(float *coeffs, const int & size)
{
    MatrixXd A(size, size);
    A.setZero();

    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            A(row, col) = pow(col+1.0, 2.0 * row + 1.0);
        }
    }
    VectorXd b(size);
    b.setZero();
    b(1) = 3.0;
    
    VectorXd x = A.fullPivLu().solve(b);
    for(int iv=0; iv<size; ++iv) coeffs[iv] = x(iv);
}

bool ThirdDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs)
{
    int nv = npoint / 2;
    if (nv < 2) {
       std::cerr << "Too few points !" << std::endl;
       return false;  
    }
    coeffs.resize(nv);

    ComputeThirdDerivativeCoefficients(&coeffs[0], nv);
    return true;  
}

bool ThirdDerivativeCoefficients(const int & npoint, float *&coeffs)
{
    int nv = npoint / 2;
    if (nv < 2) {
       std::cerr << "Too few points !" << std::endl;
       return false;  
    }
    coeffs = (float *) malloc(nv * sizeof(float));
    ComputeThirdDerivativeCoefficients(&coeffs[0], nv);
    return true;  
}

// Routine: FourthDerivativeCoefficients:
// Input:
// npoint:  no of points used
// output:
// coeffs: coefficients for one side include the central point
// coeffs[0]: weight for x0
// coeffs[1]: weight for point x0+h
// ....
void ComputeFourthDerivativeCoefficients(float *coeffs, const int &size)
{
    MatrixXd A(size, size);
    A.setZero();
    
    float fact = 1.f;
    for (int row = 0; row < size; row++) {
        for (int col = 0; col < size; col++) {
            if (row == 0) {
               A(row,col) = (col > 0) ? 2.0 : 1.0;
            } else {
               A(row,col) = pow(1.f * col, 2.0 * row) / fact;
            }
        }
        fact *= (2.0 * row+1.f) * (2.0 * row + 2.0);
    }
    VectorXd b(size);
    b.setZero();
    b(2) = 0.5;
    VectorXd x = A.fullPivLu().solve(b);
    for(int iv=0; iv<size; ++iv) coeffs[iv] = x(iv);
}

bool FourthDerivativeCoefficients(const int & npoint, std::vector<float> &coeffs)
{
    int nv = npoint / 2 + 1;
    if (nv < 2) {
       std::cerr << "Too few points !" << std::endl;
       return false;
    }
    coeffs.resize(nv);
    ComputeFourthDerivativeCoefficients(&coeffs[0], nv);
    return true;
}

bool FourthDerivativeCoefficients(const int & npoint, float *&coeffs)
{
    int nv = npoint / 2 + 1;
    if (nv < 2) {
       std::cerr << "Too few points !" << std::endl;
       return false;
    }
    coeffs = (float *) malloc(nv * sizeof(float));
    ComputeFourthDerivativeCoefficients(&coeffs[0], nv);
    return true;
}

/*
 * Following two functions compute center-pointed first and second derivative coefficients
 * include interval.The code is adapted from Changhua Zhang's code
 */
void CenteredFirstDerivativeCoefficients(const int &spaceOrder, const float &delta, std::vector<float> &coeffs)
{
    int radius = spaceOrder / 2;
    coeffs.resize(radius);
    coeffs[0] = 0.f;
    for (int m=1; m<=radius; m++) {
         coeffs[m-1] = 0.5f / (m * delta);
         if (m % 2 == 0) coeffs[m-1] = -coeffs[m-1]; 
         for (int n=1; n<=radius; n++) {
             if (n != m) {
	        float temp=fabsf(1.0f/(1.0f-float(m*m)/float(n*n)));
	        coeffs[m-1] *= temp;
             }
         }
    }
}

void CenteredSecondDerivativeCoefficients(const int &spaceOrder, const float &delta, std::vector<float> &coeffs)
{
    int radius = spaceOrder / 2;
    coeffs.resize(radius + 1);
    for (int m=1; m<=radius; m++) {
        float offset = m * delta;
        float offset2 = offset * offset;
        coeffs[m] = 1.0f / offset2;
        if (m%2==0) coeffs[m] = -coeffs[m]; 
        for (int n=1; n<=radius; n++) {
            if (n != m) {
               float temp=fabs(1.0f/(1.0f-float(m*m)/float(n*n)));
	       coeffs[m]*=temp;
            }
        }
        coeffs[0] -= 2.0f * coeffs[m];
    }
}

#define TEST_FD 0
#if TEST_FD
int main(int argc, char **argv)
{
    std::vector<float> coeff;

    if (argc != 3) {
        std::cout << "Please provide derivative level and no of poinst" << std::endl;
    }
   
    int level = atoi(argv[1]);
    int np    = atoi(argv[2]);
   std::cout << "level = "<<level <<std::endl;
   std::cout << "npoint = "<<np <<std::endl;
    if ( level == 1) {
      FirstDerivativeCoefficients(np, coeff);
    } else if (level == 2) {
      SecondDerivativeCoefficients(np, coeff);
    } else if (level == 3) {
      ThirdDerivativeCoefficients(np, coeff);
    } else if (level == 4) {
      FourthDerivativeCoefficients(np, coeff);
    } else {
        std::cout << "The "<<level<<"-th derivative not implemented!" << std::endl;
    }
    std::cout << "The weight coefficients for "<<level<<"-th derivative: " << std::endl;
    for (int i=0; i<coeff.size(); ++i) {
        std::cout << coeff[i] << std::endl;
    } 
    return 0;
}
#endif
