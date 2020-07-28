#ifndef GENERAL_ABSORBING_BOUNDARY_H
#define GENERAL_ABSORBING_BOUNDARY_H

// Yezheng (Jim) Hu
// May 11 2015
// TODO
// This class is empty for now and will be implenetd as I have some clear ideas
class AbsorbingBoundary
{
    public:
      AbsorbingBoundary(void): m_nlayer(0) {}
      virtual ~AbsorbingBoundary(void) {;}
      int NumberOfLayers(void) const {return m_nlayer; }

      void setNumberOfLayers(const int nlayer) {m_nlayer = nlayer;}
    private:
      int m_nlayer;
};

#endif
