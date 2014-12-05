#ifndef _LIB_H_
#define _LIB_H_

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <functional>
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

    /// Get the value of the ith element. Column-major.
    double get(unsigned i) const;

    /// Get the value of the Aij element.
    double get(unsigned i, unsigned j) const;

    /// Set the value of the ith element. Column-major.
    void set(unsigned i, double value);

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

    /// Transpose alias.
    Matrix t() const;

    /// Return true if this is a vector (i.e., a row vector or a column vector).
    bool isVector() const;

    /// IF this is a 2x2 matrix, return its determinant.
    double det2() const;

    /// If this is a vector, return its module.
    double mod() const;

    /// If this is a 1x1 vector, return its only element.
    double x() const;

    /// If this is a 2x1 vector, return its first element.
    double x1() const;

    /// If this is a 2x1 vector, return its second element.
    double x2() const;

    /// Return a string representation of this matrix.
    void debug() const;

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

void Matrix::set(unsigned i, double value) {
  v.at((i - 1) % m).at((i-1) / m) = value;
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

Matrix Matrix::t() const {
  return this->transpose();
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

double Matrix::x() const {
  if (m == 1 && n == 1)
    return get(1,1);
  else
    throw std::invalid_argument("Not a 1x1 Matrix");
}

double Matrix::x1() const {
  if (m == 2 && n == 1)
    return get(1,1);  
  else
    throw std::invalid_argument("Not a 2x1 column vector");
}

double Matrix::x2() const {
  if (m == 2 && n == 1)
    return get(2,1);  
  else
    throw std::invalid_argument("Not a 2x1 column vector");
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

void Matrix::debug() const {
  std::cout << "INFO: Matrix debug" << std::endl;
  std::cout << "\t" << "#rows=" << m << ", #cols=" <<  n << std::endl;
  for (unsigned i = 1; i <= m; ++i) {
    std::cout << "\t";
    for (unsigned j = 1; j <= n; ++j)
      std::cout << get(i,j) << " ";
    std::cout << std::endl;
  }
}

/**
 * Classe para contar o tempo de um método. 
 * Como usar: 
 *    Timer timer;
 *    <rotina demorada>;
 *    timer.elapsed();
 * Retorna o tempo em segundos.
 */
class Timer {
  public:
    Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }
    double elapsed() {
      clock_gettime(CLOCK_REALTIME, &end_);
      return end_.tv_sec - beg_.tv_sec +
        (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
    }
    void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }
  private:
    timespec beg_, end_;
};

/// Retorna um inteiro aleatório entre -limit e +limit (não-incluso).
int rand_int(unsigned limit) { 
  int signal = ((rand() % 2) == 0) ? (-1) : (+1);
  return signal * (rand() % limit); 
};

/// Retorna um número aleatório entre -limit e +limit.
int rand_double(unsigned limit) {
  int signal = ((rand() % 2) == 0) ? (-1) : (+1);
  double decimals = (rand() % 1000) / 1000.0;
  return signal * ((rand() % limit) + decimals); 
};

/// fa
double fa(Matrix x) {
  return pow(x.x1(), 2) + pow(exp(x.x1()) - x.x2(), 2.0);
}

/// gradiente de fa
Matrix gradfa(Matrix x) {
  Matrix w(2,1);
  w.set(1, (2 * x.x1()) + 2 * ( exp(x.x1()) - x.x2() ) * exp(x.x1()));
  w.set(2,                2 * ( exp(x.x1()) - x.x2() ) * (-1));
  return w;
}

/// fb
double fb(Matrix x) {
  return sqrt(fa(x));
}

/// gradiente de fb
Matrix gradfb(Matrix x) {
  Matrix w(2,1);
  w.set(1, (1.0/(2 * fb(x))) * gradfa(x).x1());
  w.set(2, (1.0/(2 * fb(x))) * gradfa(x).x2());
  return w;
}

/// fc
double fc(Matrix x) {
  return log(1.0 + fa(x));
}

/// gradiente de fc
Matrix gradfc(Matrix x) {
  Matrix w(2,1);
  w.set(1, (1.0/fc(x)) * gradfa(x).x1());
  w.set(2, (1.0/fc(x)) * gradfa(x).x2());
  return w;
}

/**
 * d do subproblema.
 * x: a variável
 * xkk: o ponto anterior
 */
double d(Matrix x, Matrix xkk) {
  return
    pow(x.x1() - xkk.x1(), 2.0) +
    pow((x.x2() - xkk.x2()) - (exp(x.x1()) - exp(xkk.x1())), 2.0);
}

/// gradiente de d
Matrix gradd(Matrix x, Matrix xkk) {
  Matrix w(2,1);

  w.set(1,
      2 * (x.x1() - xkk.x1()) +
      2 * (-exp(x.x1())) * ((x.x2() - xkk.x2()) - (exp(x.x1()) - exp(xkk.x1())))
      );

  w.set(2,
      2 * ((x.x2() - xkk.x2()) - (exp(x.x1()) - exp(xkk.x1())))
      );

  return w;
}

/**
 * g do subproblema
 * g = f + (lambdak / 2) * d
 */
double g(
    std::function<double(Matrix)> f,
    double lambdak,
    Matrix x,
    Matrix xkk
    )
{
  return f(x) + (((lambdak/2.0) * d(x,xkk)));
}

/// gradiente de g
Matrix gradg(
    std::function<Matrix(Matrix)> gradf,
    double lambdak,
    Matrix x,
    Matrix xkk
    )
{
  return gradf(x) + ((lambdak/2.0) * gradd(x,xkk));
}

/**
 * Regra de Armijo.
 * Encontrar um t = sb^m tal que
 * f(x + sb^m d) - f(x) <= o sb^m gradf(x)' d ==>
 * f(x + td) - f(x) <= o t gradf(x)' d
 *
 * Constraints:
 *    m >= 0
 *    0 < o < 1 (sigma)
 *    0 < b < 1 (beta)
 *    0 << t < 1 (bs^m)
 */
double armijo_call(
    double s,
    double beta,              // 0 < b < 1
    double sigma,             // 0 < o < 1
    std::function<double(Matrix)> f,
    std::function<Matrix(Matrix)> gradf,
    const Matrix& x,
    const Matrix& d
    ) {
  std::cout << "\t" << "INFO: armijo_call run" << std::endl;

  // Skipping right through the test means 1 iteration.
  unsigned iter = 0;

  while (true) {
    if ( (f(x) - f(x + s * pow(beta, iter) * d)) >= -sigma * s * pow(beta, iter) * ((gradf(x)).t() * d).x() )
      break;
    ++iter;
  }

  double t = s * pow(beta, iter);
  std::cout << "\t\t" << "#iter=" << iter+1 << ", t=" << t << std::endl;
  return t;
}

Matrix gradient_method(
    std::function<double(Matrix)> f,
    std::function<Matrix(Matrix)> gradf,
    Matrix x0,
    double epsilon
    ) {
  std::cout << "INFO: gradient_method run" << std::endl;
  std::cout << "initial point: " << "(" << x0.x1() << ", " << x0.x2() << ")" << std::endl;

  Timer timer;
  Matrix dk;                  // descida (o gradiente)
  Matrix xk = x0;             // x atual
  unsigned iter = 0;          // Iteração
  unsigned n_call_armijo = 0; // Número de chamadas de Armijo.

  while(true) {
    // Critério de parada.
    if ((gradf(xk)).mod() < epsilon) 
      break;
    ++iter;

    dk = (-1) * gradf(xk);

    double ak = armijo_call(0.8, 0.8, 0.8, f, gradf, xk, dk);
    ++n_call_armijo;

    // Atualização do xk.
    xk = xk + ak * dk;

    std::cout << "\tINFO: gradient_method iter " << iter << std::endl;
    std::cout << "\t\t" << "dk: " << "(" << dk.x1() << ", " << dk.x2() << ")" << std::endl;
    std::cout << "\t\t" << "xk: " << "(" << xk.x1() << ", " << xk.x2() << ")" << std::endl;
    std::cout << "\t\t" << "f(xk): " << f(xk) << std::endl;
  }

  std::cout << "\t" << "elapsed time:" << timer.elapsed() << "s" << std::endl;
  std::cout << "\t" << "initial point: " << "(" << x0.x1() << ", " << x0.x2() << ")" << std::endl;
  std::cout << "\t" << "epsilon: " << epsilon << std::endl;
  std::cout << "\t" << "n_iterations: " << iter + 1 << std::endl;
  std::cout << "\t" << "n_call_armijo: " << n_call_armijo << std::endl;
  std::cout << "\t" << "optimal point: " << "(" << xk.x1() << ", " << xk.x2() << ")" << std::endl;
  std::cout << "\t" << "optimal value: " << f(xk) << std::endl;

  return xk;
}

Matrix newton_method(
    std::function<double(Matrix)> f,
    std::function<Matrix(Matrix)> gradf,
    std::function<Matrix(Matrix)> invhessf,
    Matrix x0,
    double epsilon
    ) {
  std::cout << "INFO: newton_method run" << std::endl;
  std::cout << "initial point: " << "(" << x0.x1() << ", " << x0.x2() << ")" << std::endl;

  Timer timer;
  Matrix dk;
  Matrix xk = x0;
  unsigned iter = 0;
  unsigned n_call_armijo = 0;

  while(true) {
    // Critério de parada.
    if ((gradf(xk)).mod() < epsilon)
      break;
    ++iter;

    dk = (-1) * invhessf(xk) * gradf(xk);

    double ak = armijo_call(0.8, 0.8, 0.8, f, gradf, xk, dk);
    ++n_call_armijo;

    // Atualização do xk.
    xk = xk + ak * dk;

    std::cout << "\tINFO: newton_method iter " << iter << std::endl;
    std::cout << "\t\t" << "dk: " << "(" << dk.x1() << ", " << dk.x2() << ")" << std::endl;
    std::cout << "\t\t" << "xk: " << "(" << xk.x1() << ", " << xk.x2() << ")" << std::endl;
    std::cout << "\t\t" << "f(xk): " << f(xk) << std::endl;
  }

  std::cout << "\t" << "elapsed time:" << timer.elapsed() << "s" << std::endl;
  std::cout << "\t" << "initial point: " << "(" << x0.x1() << ", " << x0.x2() << ")" << std::endl;
  std::cout << "\t" << "epsilon: " << epsilon << std::endl;
  std::cout << "\t" << "n_iterations: " << iter + 1 << std::endl;
  std::cout << "\t" << "n_call_armijo: " << n_call_armijo << std::endl;
  std::cout << "\t" << "optimal point: " << "(" << xk.x1() << ", " << xk.x2() << ")" << std::endl;
  std::cout << "\t" << "optimal value: " << f(xk) << std::endl;

  return xk;
}

enum {POSTO1, POSTO2};

Matrix quasinewton_method(
    double (*f)(Matrix),
    Matrix (*gradf)(Matrix),
    Matrix x0,
    Matrix B0,
    int posto,
    double epsilon
    ) {
  // TODO.
  return Matrix();
}

enum {FA, FB, FC};
enum {GRADIENT, NEWTON, QUASINEWTON};

Matrix solve_it(
    int function,
    Matrix x0sub,
    int limitx0,
    double epsilonSub,
    double epsilonMeth,
    int method
    ) 
{
  std::cout << "INFO: solve_it run" << std::endl;
  std::cout << "initial point: " << "(" << x0sub.x1() << ", " << x0sub.x2() << ")" << std::endl;

  Timer timer;
  Matrix xk = x0sub;
  Matrix xnext;
  unsigned iter = 0;

  while(true) {
    ++iter;

    // Inicialização do problema de otimização
    Matrix x0(2,1);
    x0.set(1, rand_double(limitx0));
    x0.set(2, rand_double(limitx0));

    // Lambdak
    double lambdak = 1.0/iter;

    // Resolver um problema de otimização.
    if (method == GRADIENT) {
      switch(function) {
        case FA:
          xnext = gradient_method(
              [lambdak,xk](Matrix x) -> double { return g(fa, lambdak, x, xk); },
              [lambdak,xk](Matrix x) -> Matrix { return gradg(gradfa, lambdak, x, xk); },
              x0,
              epsilonMeth
              );
          // xnext = gradient_method( 
          //<g(fa, lambdak, _1, xk)>,
          //<gradg(gradfa, lambdak, _1, xk)>,
          //x0,
          //epsilonMeth);
          break;
        case FB:
          xnext = gradient_method(
              [lambdak,xk](Matrix x) -> double { return g(fb, lambdak, x, xk); },
              [lambdak,xk](Matrix x) -> Matrix { return gradg(gradfb, lambdak, x, xk); },
              x0,
              epsilonMeth
              );
          break;
        case FC:
          xnext = gradient_method(
              [lambdak,xk](Matrix x) -> double { return g(fc, lambdak, x, xk); },
              [lambdak,xk](Matrix x) -> Matrix { return gradg(gradfc, lambdak, x, xk); },
              x0,
              epsilonMeth
              );
          break;
      }
    }

    // Critério de parada 1.
    if ((xnext - xk).mod() < epsilonSub) {
      std::cout << "\t" << "Finished solve_it. Reason: xnext near to xk" << std::endl;
      std::cout << "\t\t" << "xnext (optimal point):" << "(" << xnext.x1() << ", " << xnext.x2() << ")" << std::endl;
      std::cout << "\t\t" << "xk (previous one): " << "(" << xk.x1() << ", " << xk.x2() << ")" << std::endl;
      break;
    }

    // Critério de parada 2.
    bool terminate = false;
    switch(function) {
      case FA:
        if (((*gradfa)(xnext)).mod() < epsilonSub)
          terminate = true;
        break;
      case FB:
        if (((*gradfb)(xnext)).mod() < epsilonSub)
          terminate = true;
        break;
      case FC:
        if (((*gradfc)(xnext)).mod() < epsilonSub)
          terminate = true;
        break;
    }
    if (terminate) {
      std::cout << "\t" << "Finished solve_it. Reason: f(xnext) near to zero" << std::endl;
      std::cout << "\t\t" << "xnext (optimal point):" << "(" << xnext.x1() << ", " << xnext.x2() << ")" << std::endl;
      break;
    }

    xk = xnext;
  }

  std::cout << "\t" << "elapsed time:" << timer.elapsed() << "s" << std::endl;
  std::cout << "\t" << "initial point: " << "(" << x0sub.x1() << ", " << x0sub.x2() << ")" << std::endl;
  std::cout << "\t" << "epsilon: " << epsilonSub << std::endl;
  std::cout << "\t" << "n_iterations: " << iter << std::endl;
  std::cout << "\t" << "optimal point: " << "(" << xnext.x1() << ", " << xnext.x2() << ")" << std::endl;

  switch(function) {
    case FA:
      std::cout << "\t" << "optimal value: " << (*fa)(xnext) << std::endl;
      break;
    case FB:
      std::cout << "\t" << "optimal value: " << (*fb)(xnext) << std::endl;
      break;
    case FC:
      std::cout << "\t" << "optimal value: " << (*fc)(xnext) << std::endl;
      break;
  }
  return xnext;
}

#endif
