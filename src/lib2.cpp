#include "lib.hpp"

Matrix::Matrix() :
  m(0),n(0) {
  v.clear();
}

Matrix::Matrix(const Matrix& o) :
  m(o.m),n(o.n) {
  v = vector< vector<double> >(m, vector<double>(n, 0.0));
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < n; ++j) 
      v[i][j] = o.v[i][j];
}

Matrix::Matrix(unsigned m, unsigned n, double value) :
  m(m), n(n) {
  v = vector< vector<double> >(m, vector<double>(n, value));
}

Matrix::Matrix(const vector<double>& o) :
  m(1), n(o.size()) {
  v = vector < vector<double > >(1, vector<double>(n, 0.0));
  for (unsigned i = 0; i < o.size(); ++i)
    v[0][i] = o[i];
}

Matrix::Matrix(const vector< vector<double> >& o) :
  m(o.size()), n(o[0].size()) {
  v = vector< vector<double> >(m, vector<double>(n, 0.0));
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < n; ++j) 
      v[i][j] = o[i][j];
}

Matrix Matrix::operator+(const Matrix& o) const {
  if (!(m == o.m && n == o.n))
    throw std::invalid_argument("matrices should have the same dimensions for this operation");

  Matrix w(m, n);
  for (unsigned i = 0; i < m; ++i)
    for(unsigned j = 0; j < n; ++j)
      w.v[i][j] = v[i][j] + o.v[i][j];
  return w;
}

Matrix Matrix::operator-(const Matrix& o) const {
  if (!(m == o.m && n == o.n))
    throw std::invalid_argument("matrices should have the same dimensions for this operation");
  Matrix w(m, n);
  for (unsigned i = 0; i < m; ++i)
    for(unsigned j = 0; j < n; ++j)
      w.v[i][j] = v[i][j] - o.v[i][j];
  return w;
}

Matrix Matrix::operator*(const Matrix& o) const {
  if (! (n == o.m))
    throw std::invalid_argument("invalid matrix multiplication");
  Matrix w(m, o.n);
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < o.n; ++j) {
      double tmp = 0.0;
      for (unsigned k = 0; k < n; ++k)
        tmp += v[i][k] * o.v[k][j];
      w.v[i][j] = tmp;
    }
  return w;
}

Matrix Matrix::operator*(const double a) const {
  Matrix w(*this);
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0 ; j < n; ++j)
      w.v[i][j] *= a;
  return w;
}

Matrix Matrix::operator/(const double a) const {
  Matrix w(*this);
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < n; ++j)
      w.v[i][j] /= a;
  return w;
}

Matrix Matrix::eye_inplace(const int m) {
  this->m = this->n = m;
  v = vector< vector<double> >(m, vector<double>(m, 0.0));
  for (unsigned i = 0; i < m; ++i)
    v[i][i] = 1.0;
  return v;
}

void Matrix::debug() const {
  cout << "INFO: Matrix debug" << endl;
  cout << "\t" << "m = " << m << " x " << "n = " << n << endl;
  for (unsigned i = 0; i < m; ++i) {
    cout << "\t";
    for (unsigned j = 0; j < n; ++j)
      cout << v[i][j] << " ";
    cout << endl;
  }
}

double Matrix::x() const {
  if ( !(m == 1 && n == 1) )
    throw std::logic_error("x() must be applied only on a 1x1 vector");
  return v[0][0];
}
double Matrix::x1() const {
  if ( !(m == 1 && n == 2) )
    throw std::logic_error("x1() must be applied only on a 1x2 vector");
  return v[0][0];
}
double Matrix::x2() const {
  if ( !(m == 1 && n == 2) )
    throw std::logic_error("x2() must be applied only on a 1x2 vector");
  return v[0][1];
}

double Matrix::mod() const {
  if ( !(m == 1 || n == 1)) 
    throw std::logic_error("mod operation must be applied only on a vector");
  double ret = 0.0;
  if (m == 1)
    for (unsigned i = 0; i < n; ++i)
      ret += v[0][i] * v[0][i];
  else
    for (unsigned i = 0; i < m; ++i)
      ret += v[i][0] * v[i][0];
  return sqrt(ret);
}

double Matrix::det2() const {
  if ( !(m == 2 || n == 2)) 
    throw std::logic_error("det2() must be applied on a 2x2 matrix");
  return v[0][0] * v[1][1] - v[0][1] * v[1][0];  
}

void Matrix::transpose_inplace() {
  Matrix o(*this);
  v = vector< vector<double> >(n, vector<double>(m, 0.0));
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < n; ++j)
      v[j][i] = o.v[i][j];
  swap(m, n);
}

Matrix Matrix::transpose() const {
  Matrix w(n, m);
  for (unsigned i = 0; i < m; ++i)
    for (unsigned j = 0; j < n; ++j)
      w.v[j][i] = v[i][j];
  return w;
}

Matrix eye(const int m) {
  Matrix w(m, m, 0.0);
  for (unsigned i = 0; i < m; ++i)
    w.v[i][i] = 1.0;
  return w;
}

Matrix operator*(const double a, const Matrix& m) {
  Matrix w(m * a);
  return w;
}

Matrix method_grad(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, const int search, const double epsilon) {
  Matrix xk = x0, xprev = x0;
  Matrix dk;
  unsigned k = 0, call_armijo = 0;
  while( ((*grad)(xprev)).mod() >= epsilon ) {
    cout << "INFO: grad iter" << endl;
    printf("\txprev: (%.6lf, %.6lf)\n", xprev.x1(), xprev.x2());
    dk = (*grad)(xprev);
    printf("\tdk: (%.6lf, %.6lf)\n", dk.x1(), dk.x2());
    if (search == PURE)
      xk = xprev - dk;
    else if (search == HARMONIC) 
      xk= xprev - (1.0/(k+1)) * dk;
    else {
      ++call_armijo;
      double lambdak = armijo(0.01, 0.8, 0.8, xprev, f, grad, dk);
      xk = xprev - lambdak * dk;
    }
    ++k;
    xprev = xk;
  }
  Matrix opt = xk;
  cout << "INFO: gradient method" << endl;
  printf("\tx0: (%.6lf, %.6lf)\n", x0.v[0][0], x0.v[0][1]); 
  cout << "\titer: " << k << endl;
  cout << "\tcall_armijo: " << call_armijo << endl;
  printf("\topt_point: (%.6lf, %.6lf)\n", opt.v[0][0], opt.v[0][1]);
  printf("\topt_value: %.6lf\n", (*f)(opt));
  return opt;
}

Matrix method_newton(double (*f)(Matrix), Matrix (*grad)(Matrix), Matrix (*inv_hess)(Matrix), const Matrix x0, const int search, const double epsilon) {
  Matrix xk = x0, xprev = x0;
  Matrix dk;
  unsigned k = 0, call_armijo = 0;
  while( ((*grad)(xprev)).mod() >= epsilon ) {
    dk = ((*inv_hess)(xprev) * (*grad)(xprev).transpose()).transpose();
    printf("\tdk: (%.6lf, %.6lf)\n", dk.v[0][0], dk.v[0][1]);
    ++call_armijo;
    double lambdak = armijo(0.8, 0.8, 0.8, xprev, f, grad, dk);
    xk = xprev - lambdak * dk;
    ++k;
    xprev = xk;
  }
  Matrix opt = xk;
  cout << "INFO: Newton method" << endl;
  printf("\tx0: (%.6lf, %.6lf)\n", x0.v[0][0], x0.v[0][1]); 
  cout << "\titer: " << k << endl;
  cout << "\tcall_armijo: " << call_armijo << endl;
  printf("\topt_point: (%.6lf, %.6lf)\n", opt.v[0][0], opt.v[0][1]);
  printf("\topt_value: %.6lf\n", (*f)(opt));
  return opt;
}

Matrix method_quasinewton(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, Matrix B0, const int search, const int posto, const double epsilon) {
  Matrix xk = x0, xprev = x0;
  Matrix bk = B0, bprev = B0;
  Matrix dk, sk, yk;
  unsigned k = 0, call_armijo = 0;
  while(((*grad)(xprev)).mod() >= epsilon) {
    cout << "INFO: quasi newton iter" << endl;
    printf("\txprev: (%.6lf, %.6lf)\n", xprev.x1(), xprev.x2());
    dk = (-1.0) * (bprev * ((*grad)(xprev)).transpose()).transpose();
    printf("\tdk: (%.6lf, %.6lf)\n", dk.x1(), dk.x2());
    if (search == PURE)
      xk = xprev + dk;
    else if (search == HARMONIC)
      xk = xprev + (1/(k+1)) * dk;
    else {
      ++call_armijo;
      double lambdak = armijo(0.8, 0.8, 0.8, xprev, f, grad, dk);
      xk = xprev + lambdak * dk;
    }
    sk = xk - xprev;
    yk = (*grad)(xk) - (*grad)(xprev);
    if (posto == POSTO1) {
        double alguma1_den = (((yk.transpose() - (bprev * sk.transpose())).transpose()) * sk.transpose()).x();
        Matrix alguma1 = (yk.transpose() - bprev * sk.transpose()) * ((yk.transpose() - bprev * sk.transpose()).transpose());
        bk = bprev + (alguma1 / alguma1_den);
    }
    else { // POSTO2
        Matrix alguma2_1 = (-1) * ((bprev * sk.transpose()) * ((bprev * sk.transpose())).transpose()) / ((sk * (bprev * sk.transpose())).x());
        Matrix alguma2_2 = (yk.transpose() * yk) / ((yk * sk.transpose()).x());
        Matrix alguma2 = alguma2_1 + alguma2_2;
        bk = bprev + alguma2;
    }
    ++k;
    xprev = xk;
    bprev = bk;
  }
  Matrix opt = xk;
  cout << "INFO: Quasi Newton method" << endl;
  printf("\tx0: (%.6lf, %.6lf)\n", x0.v[0][0], x0.v[0][1]); 
  cout << "\titer: " << k << endl;
  cout << "\tcall_armijo: " << call_armijo << endl;
  printf("\topt_point: (%.6lf, %.6lf)\n", opt.v[0][0], opt.v[0][1]);
  printf("\topt_value: %.6lf\n", (*f)(opt));
  return opt;
}

double armijo(double s, double beta, double sigma, const Matrix& xk, double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix& dk) {
  unsigned m = 0;
  while ( ((*f)(xk + s * pow(beta, m) * dk) - (*f)(xk)) < sigma * s * pow(beta, m) * ((*grad)(xk) * dk.transpose()).v[0][0])
    ++m;
  cout << "INFO: armijo" << endl;
  cout << "\tIter: m = " << m << endl;
  cout << "\tLambdak = " << s * pow(beta, m) << endl;
  return s * pow(beta, m); /* lambda_k */
}


