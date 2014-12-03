#ifndef _LIB_H_
#define _LIB_H_

enum {PURE, HARMONIC, ARMIJO};
enum {POSTO1, POSTO2};
enum {GRADIENT, NEWTON, QUASINEWTON};
Matrix method_grad(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, const int search = PURE, const double epsilon = 1e-4);
Matrix method_newton(double (*f)(Matrix), Matrix (*grad)(Matrix), Matrix (*inv_hess)(Matrix), const Matrix x0, const int search = PURE, const double epsilon = 1e-4);
Matrix method_quasinewton(double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix x0, Matrix B0, const int search = PURE, const int posto = POSTO1, const double epsilon = 1e-4);
inline double armijo(double s, double beta, double sigma, const Matrix& xk, double (*f)(Matrix), Matrix (*grad)(Matrix), const Matrix& dk);

#endif
