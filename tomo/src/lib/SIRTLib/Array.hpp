#ifndef MULTIDIMENSIONAL_ARRAY_H
#define MULTIDIMENSIONAL_ARRAY_H

#include <cstring>
#include <vector>

// Author: Yezheng (Jim) Hu
// Date:   11 May 2015
//
// Template class for array from one to five dimensions
// May 20 2015 Add overload = operator to assign the same value for all array elements

// ### 1D array class ###
template<class T>
class Array1D 
{
public:
    Array1D(void) {}
    Array1D(size_t n) {
        Resize(n);
    }
    Array1D(const Array1D<T> &orig) : m_pool(orig.m_pool) {}
    const Array1D<T> &operator= (T rhs) { m_pool.assign (m_pool.size(), rhs); return *this; }
    const Array1D<T> &operator=(const Array1D<T> &rhs);
    const size_t N1(void) const {
        return m_pool.size();
    }
    void Resize(const size_t n) {
      std::vector<T> x (n);
      m_pool.swap (x);
    }
    void Zero(void) {
        memset ((void *)&m_pool[0], 0, m_pool.size() * sizeof(m_pool[0]));
    }
    operator T*() {
        return &m_pool[0];
    }
    operator const T*() const {
        return &m_pool[0];
    }
private:
    std::vector<T> m_pool;
};

template <class T>
const Array1D<T> &Array1D<T>::operator= (const Array1D<T> &rhs)
{
    if (this != &rhs) {
        m_pool = rhs.m_pool;
    }
    return (*this);
}

// ### 2D array class ###
template <class T>
class Array2D
{
public:
    Array2D(void);
    Array2D(const size_t n1, const size_t n2);
    Array2D (const Array2D<T> & orig);
    const Array2D<T> &operator= (T rhs) { m_pool.assign (m_pool.size(), rhs); return *this; }
    const Array2D<T> &operator=(const Array2D<T> &rhs);
    void Resize(const size_t n1, const size_t n2);
    
    void Zero(void) {
        memset ((void *)&m_pool[0], 0, m_pool.size() * sizeof(m_pool[0]));
    }
    const size_t N1(void) const {
        return m_n1;
    }
    const size_t N2(void) const {
        return m_n2;
    }
    operator T**() {
        return &m_row[0];
    }
    operator const T* const*() const {
        return &m_row[0];
    }
private:
    void MapIndex(void);
    void Copy(const Array2D<T> &rhs);
private:
    size_t m_n1;
    size_t m_n2;

    std::vector<T> m_pool;
    std::vector<T *> m_row;
};

template <class T>
Array2D<T>::Array2D(void)
{
    m_n1 =  0;
    m_n2 = 0;
}

template <class T>
Array2D<T>::Array2D (const Array2D<T> &orig)
        : m_n1 (orig.m_n1),
        m_n2 (orig.m_n2),
        m_pool (orig.m_pool)
{
    MapIndex();
}

template <class T>
Array2D<T>::Array2D (const size_t n1, const size_t n2)
{
    Resize (n1, n2);
}

template <class T>
void Array2D<T>::Resize (const size_t n1, const size_t n2)
{
    m_n1 = n1;
    m_n2 = n2;
    
    std::vector<T> x(n1 * n2);
    m_pool.swap (x);
    MapIndex();
}

template <class T>
void Array2D<T>::MapIndex (void)
{
  std::vector<T *> x (m_n1);
  m_row.swap(x);
  for (size_t i = 0; i < m_n1; ++i) {
    m_row[i] = &m_pool[i * m_n2];
  }
}

template<class T>
const Array2D<T> &Array2D<T>::operator=(const Array2D<T> &rhs)
{
  size_t n1 = rhs.N1();
  size_t n2 = rhs.N2();
  Resize (n1, n2);
  Copy(rhs);
  return *this;
}

template<class T>
void Array2D<T>::Copy(const Array2D<T> &rhs)
{
  Array2D<T> &out = *this;
  for (size_t i = 0; i < m_n1; i++) {
    for (size_t j = 0; j < m_n2; j++) {
      T val = rhs[i][j];
      out[i][j] = val;
    }
  }
}

// ### 3D array class ###
template <class T>
class Array3D
{
public:
    Array3D(void);
    Array3D(const size_t n1, const size_t n2, const size_t n3);
    Array3D (const Array3D<T> & orig);
    const Array3D<T> &operator= (T rhs) { m_pool.assign (m_pool.size(), rhs); return *this; }
    const Array3D<T> &operator=(const Array3D<T> &rhs);
    void Resize(const size_t n1, const size_t n2, const size_t n3);
    void Zero(void) {
        memset ((void *)&m_pool[0], 0, m_pool.size() * sizeof(m_pool[0]));
    }
    const size_t N1(void) const {
        return m_n1;
    }
    const size_t N2(void) const {
        return m_n2;
    }
    const size_t N3(void) const {
        return m_n3;
    }
    operator T***() {
        return &m_dim1P[0];
    }
    operator const T* const* const*() const {
        return &m_dim1P[0];
    }
  
private:
    void MapIndex(void);
    void Copy(const Array3D<T> &rhs);
private:
    size_t m_n1;
    size_t m_n2;
    size_t m_n3;

    std::vector<T> m_pool;
    std::vector<T **> m_dim1P;
    std::vector<T *> m_dim2P;
};

template <class T>
Array3D<T>::Array3D(void)
{
    m_n1 = 0;
    m_n2 = 0;
    m_n3 = 0;
    m_pool.resize(0);
    m_dim1P.resize(0);
    m_dim2P.resize(0);
}

template <class T>
Array3D<T>::Array3D (const Array3D<T> &orig)
        : m_n1 (orig.m_n1),
        m_n2 (orig.m_n2),
        m_n3 (orig.m_n3),
        m_pool (orig.m_pool)
{
    MapIndex();
}

template <class T>
Array3D<T>::Array3D (const size_t n1, const size_t n2, const size_t n3)
{
    Resize (n1, n2, n3);
}

template<class T>
const Array3D<T> &Array3D<T>::operator=(const Array3D<T> &rhs)
{
  size_t n1 = rhs.N1();
  size_t n2 = rhs.N2();
  size_t n3 = rhs.N3();
  Resize (n1, n2, n3);
  Copy(rhs);
  return *this;
}

template<class T>
void Array3D<T>::Copy(const Array3D<T> &rhs)
{
  Array3D<T> &out = *this;
  for (size_t i = 0; i < m_n1; i++) {
    for (size_t j = 0; j < m_n2; j++) {
      for (size_t k = 0; k < m_n3; k++) {
        T val = rhs[i][j][k];
        out[i][j][k] = val;
      }
    }
  }
}

template <class T>
void Array3D<T>::Resize (const size_t n1, const size_t n2, const size_t n3)
{
    m_n1 = n1;
    m_n2 = n2;
    m_n3 = n3;
    std::vector<T> x (n1 * n2 * n3);
    m_pool.swap(x);

    MapIndex();
}

template <class T>
void Array3D<T>::MapIndex (void)
{
  std::vector<T **> x1 (m_n1);
  std::vector<T *>  x2 (m_n1 * m_n2);
  m_dim1P.swap(x1);
  m_dim2P.swap(x2);

  for (size_t i = 0; i < m_n1; ++i) {
    m_dim1P[i] = &m_dim2P[m_n2 * i];
    for (size_t j = 0; j < m_n2; ++j) {
      m_dim1P[i][j] = &m_pool[(j + m_n2 * i) * m_n3];
    }
  }
}

// ### 4D array class ###
template <class T>
class Array4D
{
public:
    Array4D(void);
    Array4D(const size_t n1, const size_t n2, const size_t n3, const size_t n4);
    Array4D (const Array4D<T> & orig);
    const Array4D<T> &operator= (T rhs) { m_pool.assign (m_pool.size(), rhs); return *this; }
    const Array4D<T> &operator=(const Array4D<T> &rhs);
    void Resize(const size_t n1, const size_t n2, const size_t n3, const size_t n4);
    void Zero(void) {
        memset ((void *)&m_pool[0], 0, m_pool.size() * sizeof(m_pool[0]));
    }
    const size_t N1(void) const {
        return m_n1;
    }
    const size_t N2(void) const {
        return m_n2;
    }
    const size_t N3(void) const {
        return m_n3;
    }
    const size_t N4(void) const {
        return m_n4;
    }
    operator T****() {
        return &m_dim1P[0];
    }
    operator const T* const* const* const*() const {
        return &m_dim1P[0];
    }
private:
    void MapIndex(void);
private:
    size_t m_n1;
    size_t m_n2;
    size_t m_n3;
    size_t m_n4;

    std::vector<T> m_pool;
    std::vector<T ***> m_dim1P;
    std::vector<T **> m_dim2P;
    std::vector<T *> m_dim3P;
};

template <class T>
Array4D<T>::Array4D(void)
{
    m_n1 = 0;
    m_n2 = 0;
    m_n3 = 0;
    m_n4 = 0;
}

template <class T>
Array4D<T>::Array4D (const Array4D<T> &orig)
        : m_n1 (orig.m_n1),
        m_n2 (orig.m_n2),
        m_n3 (orig.m_n3),
        m_n4 (orig.m_n4),
        m_pool (orig.m_pool)
{
    MapIndex();
}

template <class T>
Array4D<T>::Array4D (const size_t n1, const size_t n2, const size_t n3, const size_t n4)
{
    Resize (n1, n2, n3, n4);
}

template <class T>
void Array4D<T>::Resize (const size_t n1, const size_t n2, const size_t n3, const size_t n4)
{
    m_n1 = n1;
    m_n2 = n2;
    m_n3 = n3;
    m_n4 = n4;
    std::vector<T> x (n1 * n2 * n3 * n4);
    m_pool.swap (x);

    MapIndex();
}

template <class T>
void Array4D<T>::MapIndex(void)
{
  std::vector<T ***> x1 (m_n1);
  std::vector<T **> x2 (m_n1 * m_n2);
  std::vector<T *> x3 (m_n1 * m_n2 * m_n3);
  
  m_dim1P.swap(x1);
  m_dim2P.swap(x2);
  m_dim3P.swap(x3);
  for (size_t i = 0; i < m_n1; ++i) {
    m_dim1P[i] = &m_dim2P[m_n2 * i];
    for (size_t j = 0; j < m_n2; ++j) {
      m_dim1P[i][j] = &m_dim3P[(j + m_n2 * i) * m_n3];
      for (size_t k=0; k<m_n3; ++k) {
        m_dim1P[i][j][k] = &m_pool[((i* m_n2 + j) * m_n3 + k) * m_n4];
      }
    }
  }
}

// ### 5D array class ###
template <class T>
class Array5D
{
public:
    Array5D(void);
    Array5D(const size_t n1, const size_t n2, const size_t n3, const size_t n4, const size_t n5);
    Array5D (const Array5D<T> & orig);
    const Array5D<T> &operator= (T rhs) { m_pool.assign (m_pool.size(), rhs); return *this; }
    const Array5D<T> &operator=(const Array5D<T> &rhs);
    void Resize(const size_t n1, const size_t n2, const size_t n3, const size_t n4, const size_t n5);
    void Zero(void) {
        memset ((void *)&m_pool[0], 0, m_pool.size() * sizeof(m_pool[0]));
    }
    const size_t N1(void) const {
        return m_n1;
    }
    const size_t N2(void) const {
        return m_n2;
    }
    const size_t N3(void) const {
        return m_n3;
    }
    const size_t N4(void) const {
        return m_n4;
    }
    const size_t N5(void) const {
        return m_n5;
    }
    operator T*****() {
        return &m_dim1P[0];
    }
    operator const T* const* const* const* const*() const {
        return &m_dim1P[0];
    }
private:
    void MapIndex(void);
private:
    size_t m_n1;
    size_t m_n2;
    size_t m_n3;
    size_t m_n4;
    size_t m_n5;

    std::vector<T> m_pool;
    std::vector<T ****> m_dim1P;
    std::vector<T ***> m_dim2P;
    std::vector<T **> m_dim3P;
    std::vector<T *> m_dim4P;
};

template <class T>
Array5D<T>::Array5D(void)
{
    m_n1 = 0;
    m_n2 = 0;
    m_n3 = 0;
    m_n4 = 0;
    m_n5 = 0;
}

template <class T>
Array5D<T>::Array5D (const Array5D<T> &orig)
        : m_n1 (orig.m_n1),
        m_n2 (orig.m_n2),
        m_n3 (orig.m_n3),
        m_n4 (orig.m_n4),
        m_n5 (orig.m_n5),
        m_pool (orig.m_pool)
{
    MapIndex();
}

template <class T>
Array5D<T>::Array5D (const size_t n1, const size_t n2, const size_t n3, const size_t n4, const size_t n5)
{
    Resize (n1, n2, n3, n4, n5);
}

template <class T>
void Array5D<T>::Resize (const size_t n1, const size_t n2, const size_t n3, const size_t n4, const size_t n5)
{
    m_n1 = n1;
    m_n2 = n2;
    m_n3 = n3;
    m_n4 = n4;
    m_n5 = n5;
    std::vector<T> x(m_n1 * m_n2 * m_n3 * m_n4 * m_n5);
    m_pool.swap (x);

    MapIndex();
}

template <class T>
void Array5D<T>::MapIndex(void)
{
  std::vector<T ****> x1 (m_n1);
  std::vector<T ***> x2 (m_n1 * m_n2);
  std::vector<T **> x3 (m_n1 * m_n2 * m_n3);
  std::vector<T *> x4 (m_n1 * m_n2 * m_n3 * m_n4);
  
  m_dim1P.swap(x1);
  m_dim2P.swap(x2);
  m_dim3P.swap(x3);
  m_dim4P.swap(x4);
  for (size_t i = 0; i < m_n1; ++i) {
    m_dim1P[i] = &m_dim2P[m_n2 * i];
    for (size_t j = 0; j < m_n2; ++j) {
      m_dim1P[i][j] = &m_dim3P[(j + m_n2 * i) * m_n3];
      for (size_t k=0; k<m_n3; ++k) {
        m_dim1P[i][j][k] = &m_dim4P[(k + m_n3 * (j + m_n2 * i)) * m_n4];
        for (size_t l=0; l<m_n4; ++l) {
          m_dim1P[i][j][k][l] = &m_pool[(l + m_n4 * (k + m_n3 * (j + m_n2 * i))) * m_n5];
        }
      }
    }
  }
}

template<class T>
T **allocData (int nr, int ns) {
  T ** pp = new T*[nr];
  for (int i = 0; i < nr; i++) {
    pp[i] = new T[ns];
  }
  return pp;
}

template<class T>
void freeData (T **data, int nr) {
  T ** pp = data;
  for (int i = 0; i < nr; i++) {
    delete[] pp[i];
  }
  delete[] pp;
}

#endif
