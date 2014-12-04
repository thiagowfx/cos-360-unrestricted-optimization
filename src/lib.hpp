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
    /// Construct a empty (0x0) matrix.
    Matrix();

    /// Construct a matrix object with the specified dimensions (rows x cols), initialized to value.
    Matrix(unsigned rows, unsigned cols, double value = 0.0);

    /// Copy a existing matrix.
    Matrix(const Matrix&);

    /// Construct a column vector from a vector.
    Matrix(const vector<double>&);

    /// Construct a matrix from a vector of vectors.
    Matrix(const vector<vector<double> >&);

    /// Return the number of columns of this matrix.
    unsigned getCols() const;

    /// Return the number of rows of this matrix.
    unsigned getRows() const;

    /// Get the value of the ith element. Row-major.
    double get(unsigned i) const;

    /// Get the value of the Aij element.
    double get(unsigned i, unsigned j) const;

    /// Set Aij to value.
    void set(unsigned i, unsigned j, double value);

    // Usual matrix operations.
    Matrix operator+(const Matrix&) const;
    Matrix operator-(const Matrix&) const;
    Matrix operator*(const Matrix&) const;
    Matrix operator*(double) const;
    Matrix operator/(double) const; 
    friend Matrix operator*(double, const Matrix&);

    /// Return the transpose of this matrix.
    Matrix transpose() const;

    /// Return true if this is a vector (i.e., a row vector or a column vector).
    bool isVector() const;

    /// IF this is a 2x2 matrix, return its determinant.
    double det2() const;

    /// If this is a vector, return its module.
    double mod() const;

    /// If this is a 2x1 vector, return its first element.
    double x1() const;

    /// If this is a 2x1 vector, return its second element.
    double x2() const;

    /// Return a string representation of this matrix.
    std::string toString() const;

    /// Return the number of elements of this matrix.
    unsigned length() const;
  private:
    /** Dimensions of the Matrix.
     * m = number of lines
     * n = number of columns
     */
    unsigned m, n;

    /// Internal representation of the matrix.
    vector< vector<double> > v;
};

Matrix eye(unsigned n) {
  Matrix w(n,n,0.0);
  for (unsigned i = 1; i <= n; ++i)
    w.set(i,i,1.0);
  return w;
}

Matrix::Matrix() :
  m(0), 
  n(0) {
  v.clear();
}

Matrix::Matrix(unsigned rows, unsigned cols, double value) :
  m(rows),
  n(cols) {
  v = vector< vector<double> >(rows, vector<double>(cols, value));
}

Matrix::Matrix(const Matrix& o) :
  m(o.m),
  n(o.n),
  v(o.v) {
}

Matrix::Matrix(const vector<double>& w) :
  m(w.size()),
  n(1) {
    for (unsigned i = 0; i < w.size(); ++i)
      v.push_back(vector<double>(1, w[i]));
}

double Matrix::det2() const {
  if (m != 2 && n != 2)
    throw std::invalid_argument("Can't apply det2 to a non 2x2 matrix");
  return get(1,1) * get(2,2) - get(1,2) * get(2,1); 
}

Matrix::Matrix(const vector<vector<double> >& w) :
  m(w.size()),
  n(w[0].size()),
  v(w) {
}

unsigned Matrix::getCols() const {
  return n;
}

unsigned Matrix::getRows() const {
  return m;
}

double Matrix::get(unsigned i) const {
  return v.at((i - 1) % m).at((i-1) / m);
}

double Matrix::get(unsigned i, unsigned j) const {
  return v.at(i-1).at(j-1);
}

void Matrix::set(unsigned i, unsigned j, double value) {
  v.at(i-1).at(j-1) = value;
}
    
Matrix Matrix::operator+(const Matrix& o) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, get(i,j) + o.get(i,j));
  return a;
}

Matrix Matrix::operator-(const Matrix& o) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, get(i,j) - o.get(i,j));
  return a;
}

Matrix Matrix::operator*(double s) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, s * get(i,j));
  return a;
}

Matrix Matrix::operator/(double s) const {
  Matrix a(m, n);
  for (unsigned i = 1; i <= m; ++i)
    for (unsigned j = 1; j <= n; ++j)
      a.set(i, j, s / get(i,j));
  return a;
}

Matrix Matrix::transpose() const {
  Matrix w(getCols(), getRows());
  for (unsigned i = 1; i <= getRows(); ++i)
    for (unsigned j = 1; j <= getCols(); ++j)
      w.set(j,i,get(i,j));
  return w;
}

bool Matrix::isVector() const {
  return m == 1 || n == 1;
}

double Matrix::mod() const {
  double sum = 0.0;
  for (unsigned i = 1; i <= length(); ++i)
    sum += get(i) * get(i);
  return sqrt(sum);
}

double Matrix::x1() const {
  if (m == 2 && n == 1)
    return get(1,1);  
  else
    throw std::invalid_argument("Not a (2,1) column vector");
}

double Matrix::x2() const {
  if (m == 2 && n == 1)
    return get(2,1);  
  else
    throw std::invalid_argument("Not a (2,1) column vector");
}

unsigned Matrix::length() const {
  return m * n;
}

Matrix Matrix::operator*(const Matrix& o) const {
  Matrix w(getRows(), o.getCols());
  if (getCols() != o.getRows())
    throw std::invalid_argument("Invalid matrix multiplication");
  for (unsigned i = 1; i <= getRows(); ++i)
    for (unsigned j = 1; j <= o.getCols(); ++j) {
      double sum = 0.0;
      for (unsigned k = 1; k <= getCols(); ++k) {
        sum += get(i,k) * o.get(k,j);
      }
      w.set(i,j,sum);
    }
  return w;
}

Matrix operator*(double s, const Matrix& o) {
  return o * s;
}

std::string Matrix::toString() const {
  std::cout << "INFO: Matrix debug" << std::endl;
  std::cout << "\t" << "#rows=" << m << ", #cols=" << "n = " << n << std::endl;
  for (unsigned i = 1; i <= m; ++i) {
    std::cout << "\t";
    for (unsigned j = 1; j <= n; ++j)
      std::cout << get(i,j) << " ";
    std::cout << std::endl;
  }
}

double armijo(double s, double beta, double sigma, const Matrix& xk, double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix& dk) {
  std::cout << "INFO: Armijo run" << std::endl;
  unsigned iter = 0;
  while (((*f)(xk + s * pow(beta, iter) * dk) - (*f)(xk)) < sigma * s * pow(beta, iter) * (((*grad)(xk)).transpose() * dk).get(1,1))
    ++iter;
  double param = s * pow(beta, iter);
  std::cout << "\t#iter=" << iter+1 << ", parameter=" << param << std::endl;
  return param;
}

#endif
