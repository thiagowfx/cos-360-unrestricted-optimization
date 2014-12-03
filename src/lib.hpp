#ifndef _LIB_H_
#define _LIB_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

template<class T = double>
class Matrix {
  public:
    Matrix();
    Matrix(unsigned rows, unsigned cols, T value = 0.0);
    Matrix(const Matrix&);
    Matrix(const vector<T>&);
    Matrix(const vector<vector<T> >&);
    unsigned getCols() const;
    unsigned getRows() const;
    T get(unsigned i) const;
    T get(unsigned i, unsigned j) const;
    void set(unsigned i, unsigned j, T value);
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const; 
    Matrix transpose() const;
    bool isVector() const;
    T det2() const;
    T mod() const;
    T x1() const;
    T x2() const;
    std::string toString() const;
    unsigned length() const;
    template<class Y> friend Matrix<Y> operator*(double, const Matrix<Y>&);
  private:
    /** Dimensions of the Matrix.
     * m = number of lines
     * n = number of columns
     */
    unsigned m, n;

    /// Internal representation of the matrix.
    vector< vector<T> > v;
};

template<class T>
Matrix<T> eye(unsigned n) {
  Matrix<T> w(n,n,0.0);
  for (unsigned i = 1; i <= n; ++i)
    w.set(i,i,1.0);
  return w;
}

template<class T>
Matrix<T>::Matrix() :
  m(0), 
  n(0) {
  v.clear();
}

template<class T>
Matrix<T>::Matrix(unsigned rows, unsigned cols, T value) :
  m(rows),
  n(cols) {
  v = vector< vector<T> >(rows, vector<T>(cols, value));
}

template<class T>
Matrix<T>::Matrix(const Matrix& o) :
  m(o.m),
  n(o.n),
  v(o.v) {
}

template<class T>
Matrix<T>::Matrix(const vector<T>& w) :
  m(w.size()),
  n(1) {
    for (unsigned i = 0; i < w.size(); ++i)
      v.push_back(vector<T>(1, w[i]));
}

template<class T>
T Matrix<T>::det2() const {
  if (m != 2 && n != 2)
    throw std::invalid_argument("Can't apply det2 to a non 2x2 matrix");
  return get(1,1) * get(2,2) - get(1,2) * get(2,1); 
}

template<class T>
Matrix<T>::Matrix(const vector<vector<T> >& w) :
  m(w.size()),
  n(w[0].size()),
  v(w) {
}

template<class T>
unsigned Matrix<T>::getCols() const {
  return n;
}

template<class T>
unsigned Matrix<T>::getRows() const {
  return m;
}

template<class T>
T Matrix<T>::get(unsigned i) const {
  return v.at((i - 1) % m).at((i-1) / m);
}

template<class T>
T Matrix<T>::get(unsigned i, unsigned j) const {
  return v.at(i-1).at(j-1);
}

template<class T>
void Matrix<T>::set(unsigned i, unsigned j, T value) {
  v.at(i-1).at(j-1) = value;
}
    
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix& o) const {
  Matrix<T> a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, get(i,j) + o.get(i,j));
  return a;
}

template<class T>
Matrix<T> Matrix<T>::operator-(const Matrix& o) const {
  Matrix<T> a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, get(i,j) - o.get(i,j));
  return a;
}

template<class T>
Matrix<T> Matrix<T>::operator*(double s) const {
  Matrix<T> a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, s * get(i,j));
  return a;
}

template<class T>
Matrix<T> Matrix<T>::operator/(double s) const {
  Matrix<T> a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, s / get(i,j));
  return a;
}

template<class T>
Matrix<T> Matrix<T>::transpose() const {
  Matrix w(getCols(), getRows());
  for (unsigned i = 1; i <= getRows(); ++i)
    for (unsigned j = 1; j <= getCols(); ++j)
      w.set(j,i,get(i,j));
  return w;
}

template<class T>
bool Matrix<T>::isVector() const {
  return m == 1 || n == 1;
}

template<class T>
T Matrix<T>::mod() const {
  T sum = 0.0;
  for (unsigned i = 1; i <= length(); ++i)
    sum += get(i) * get(i);
  return sqrt(sum);
}

template<class T>
T Matrix<T>::x1() const {
  if (m == 2 && n == 1)
    return get(1,1);  
  else
    throw std::invalid_argument("Not a (2,1) column vector");
}

template<class T>
T Matrix<T>::x2() const {
  if (m == 2 && n == 1)
    return get(2,1);  
  else
    throw std::invalid_argument("Not a (2,1) column vector");
}

template<class T>
unsigned Matrix<T>::length() const {
  return m * n;
}

template<class T>
Matrix<T> Matrix<T>::operator*(const Matrix& o) const {
  Matrix<T> w(getRows(), o.getCols());
  if (getCols() != o.getRows())
    throw std::invalid_argument("Invalid matrix multiplication");
  for (unsigned i = 1; i <= getRows(); ++i)
    for (unsigned j = 1; j <= o.getCols(); ++j) {
      T sum = 0.0;
      for (unsigned k = 1; k <= getCols(); ++k) {
        sum += get(i,k) * o.get(k,j);
      }
      w.set(i,j,sum);
    }
  return w;
}

template<class Y>
Matrix<Y> operator*(double s, const Matrix<Y>& o) {
  return o * s;
}

template<class T>
std::string Matrix<T>::toString() const {
  std::cout << "INFO: Matrix debug" << std::endl;
  std::cout << "\t" << "#rows=" << m << ", #cols=" << "n = " << n << std::endl;
  for (unsigned i = 1; i <= m; ++i) {
    std::cout << "\t";
    for (unsigned j = 1; j <= n; ++j)
      std::cout << get(i,j) << " ";
    std::cout << std::endl;
  }
}

template<class T>
double armijo(double s, double beta, double sigma, const Matrix<T>& xk, double (*f)(Matrix<T>), Matrix<T> (*grad)(Matrix<T>), const Matrix<T>& dk) {
  std::cout << "INFO: Armijo run" << std::endl;
  unsigned iter = 0;
  while (((*f)(xk + s * pow(beta, iter) * dk) - (*f)(xk)) < sigma * s * pow(beta, iter) * ((*grad)(xk) * dk.transpose()).get(1,1))
    ++iter;
  double param = s * pow(beta, iter);
  std::cout << "\t#iter=" << iter+1 << ", parameter=" << param << std::endl;
  return param;
}

#endif
