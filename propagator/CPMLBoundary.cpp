#include<cmath>

#include "CPMLBoundary.h"

// Yezheng (Jim) Hu
// May 12 2015
//

/* Papers used to implement this class
 * (1) Collino, F. and Tsogka, C., 2001, Application of the perfectly matched absorbing layer model to the linear
 *     elastodynamic problem in anisotropic heterogeneous media: Geophysics, Vol 66 No. 1 P.294-307
 * (2) Komatitsch, D. and Martin, R., 2007, An unsplit convolutional perfectly matched layer improved at grazing
 *     incidence for the seismic wave equation: Geophysics, Vol 72 No. 5 P.SM155-SM167
 * (3) Yang, P., 2014, A numerical tour of wave propagation
 */

CPMLBoundary::CPMLBoundary(void) : AbsorbingBoundary()
{
}

/* 
 * Input:
 * nlayer: No of boundary layers
 * dh:     of space interval
 * dt:     time step
 * pf:     principal frequency
 * vmean:  Average velocity
 * damping factor: Damping factor
 */
CPMLBoundary::CPMLBoundary(const int nlayer, const float dh, const float dt, const float maxFreq, const float vmean, const float dampingRate)
{
    //const float kmax = 1.0;
    setNumberOfLayers(nlayer);
    m_A.Resize(nlayer);
    m_B.Resize(nlayer);

    // See Page (9( of Ref (3) for meaning of R
    float R = dampingRate;
    if (R > 1E-4) R = 1E-4;
    if (R < 1E-6) R = 1E-6;
    float delta = nlayer * dh;
    // See Eq (21) of ref (1)
    float d0 = log(1.0 / R) * (3.0 * vmean / (delta + delta) );
    float alpha_max_pml = M_PI * maxFreq;
    for (int i=0; i<nlayer; ++i) {
        float frac = 1. - (i+0.)/ nlayer;
        float frac2 = frac * frac;
        float d = d0 * frac2;
        float alpha = alpha_max_pml * (1. - frac) + 0.1 * alpha_max_pml;
        //float alpha = alpha_max_pml * (1. - frac);
        //float kx = 1.0 + (kmax - 1.0) * frac2;
        // See eq. (25) of ref. (2), Coefficients A and B are swaped, one is chosen for k
        m_A[i] = exp(-(d + alpha) * dt);
        m_B[i] = d / (d + alpha) * (m_A[i]-1);
    }
}
