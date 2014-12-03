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
    Matrix(const vector<vector<T> >& v);
    unsigned getCols() const;
    unsigned getRows() const;
    T get(unsigned i, unsigned j) const;
    void set(unsigned i, unsigned j, T value);
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const; 
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
Matrix<T>::Matrix(const vector<vector<T> >& v) :
  m(v.size()),
  n(v[0].size()),
  v(v) {
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
T Matrix<T>::get(unsigned i, unsigned j) const {
  return v.at(i-1).at(j-1);
}

template<class T>
void Matrix<T>::set(unsigned i, unsigned j, T value) {
  v.at(i-1).at(j-1) = value;
}
    
template<class T>
Matrix<T> Matrix<T>::operator+(const Matrix& o) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, a.get(i,j) + get(i,j));
  return a;
}

#endif
