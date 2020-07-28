#ifndef CPML_BOUNDARY_H
#define CPML_BOUNDARY_H

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


#include "Array.hpp"
#include "AbsorbingBoundary.h"

class CPMLBoundary : public AbsorbingBoundary 
{
   public:
      CPMLBoundary(void);
      CPMLBoundary(const int nlayer, const float dh, const float dt, const float pf, const float vmean, const float dampingRate);
      const Array1D<float> & A(void) const {return m_A; }
      const Array1D<float> & B(void) const {return m_B; }
   private:
      // A, B parameters are attenuation parameters defined in Eq (25) of ref (2)
      Array1D<float> m_A, m_B;
};

#endif
