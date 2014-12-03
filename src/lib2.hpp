#ifndef _LIB_H_
#define _LIB_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace std;

class Matrix {
  public:
  Matrix(const vector<double>&);
  Matrix(const vector< vector<double> >&);
  Matrix operator+(const Matrix&) const;
  Matrix operator-(const Matrix&) const;
  Matrix operator*(const Matrix&) const;
  Matrix operator*(const double) const;
  Matrix operator/(const double) const;
  friend Matrix operator*(const double, const Matrix&);
  Matrix eye_inplace(const int);
  double x() const;
  double x1() const;
  double x2() const;
  double mod() const;
  double det2() const;
  void transpose_inplace();
  Matrix transpose() const;
  void debug() const;
  vector< vector< double > > v;
  unsigned m;
  unsigned n;
};
Matrix eye(const int);
enum {PURE, HARMONIC, ARMIJO};
enum {POSTO1, POSTO2};
enum {GRADIENT, NEWTON, QUASINEWTON};
Matrix method_grad(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, const int search = PURE, const double epsilon = 1e-4);
Matrix method_newton(double (*f)(Matrix), Matrix (*grad)(Matrix), Matrix (*inv_hess)(Matrix), const Matrix x0, const int search = PURE, const double epsilon = 1e-4);
Matrix method_quasinewton(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, Matrix B0, const int search = PURE, const int posto = POSTO1, const double epsilon = 1e-4);
inline double armijo(double s, double beta, double sigma, const Matrix& xk, double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix& dk);

#endif
