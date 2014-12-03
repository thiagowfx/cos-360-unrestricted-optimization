#include "gtest/gtest.h"
#include "lib.hpp"
using namespace std;

TEST(MatrixTest, EmptyConstructor) {
  Matrix<double> m;
  EXPECT_EQ(m.getRows(), 0);
  EXPECT_EQ(m.getCols(), 0);
}

TEST(MatrixTest, StandardConstructor) {
  Matrix<double> m(2, 1);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
}

TEST(MatrixTest, StandardConstructor2) {
  Matrix<double> m(2, 1, 10.0);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 10.0);
  EXPECT_EQ(m.get(2,1), 10.0);
}

//TEST(MatrixTest, modTest) {
  //vector< double > v = {3, 4};
  //Matrix m(v);
  //EXPECT_DOUBLE_EQ(m.mod(), 5.0);
  //m.transpose_inplace();
  //EXPECT_DOUBLE_EQ(m.mod(), 5.0);

  //Matrix n(3, 3, 0.0);
  //EXPECT_ANY_THROW(n.mod());

  //Matrix z(3, 1, 0.0);
  //EXPECT_DOUBLE_EQ(z.mod(), 0.0);
//}

//TEST(MatrixTest, xx1x2Test) {
  //Matrix m(vector<double>{1.0, 2.0});
  //EXPECT_DOUBLE_EQ(m.x1(), 1.0);
  //EXPECT_DOUBLE_EQ(m.x2(), 2.0);

  //m.transpose_inplace();
  //EXPECT_ANY_THROW(m.x1());
  //EXPECT_ANY_THROW(m.x2());

  //Matrix n1(3, 1);
  //EXPECT_ANY_THROW(n1.x1());
  //EXPECT_ANY_THROW(n1.x2());

  //Matrix n2(1, 3);
  //EXPECT_ANY_THROW(n2.x1());
  //EXPECT_ANY_THROW(n2.x2());

  //Matrix nn(3, 3);
  //EXPECT_ANY_THROW(nn.x1());
  //EXPECT_ANY_THROW(nn.x2());

  //Matrix ss(vector<double>{1.0});
  //EXPECT_DOUBLE_EQ(ss.x(), 1.0);
  //EXPECT_ANY_THROW(nn.x());
//}

//TEST(MatrixTest, scalarMultiplication) {
  //Matrix m(vector<double>{2.0, 1.0});
  //Matrix n(2 * m);
  //EXPECT_DOUBLE_EQ(n.v[0][0], 4.0);
  //EXPECT_DOUBLE_EQ(n.v[0][1], 2.0);

  //Matrix nn(m * 2);
  //EXPECT_DOUBLE_EQ(nn.v[0][0], 4.0);
  //EXPECT_DOUBLE_EQ(nn.v[0][1], 2.0);
//}

//TEST(MatrixTest, scalarDivision) {
  //Matrix m(vector<double>{2.0, 1.0});
  //Matrix n(m / 2.0);
  //EXPECT_DOUBLE_EQ(n.v[0][0], 1.0);
  //EXPECT_DOUBLE_EQ(n.v[0][1], 0.5);
//}

//TEST(MatrixTest, transposeTest) {
  //vector< vector<double > > v;
  //v = vector< vector<double> >(2, vector<double>(3, 0.0));
  //v[0][0] = 1.0;
  //v[0][1] = 2.0;
  //v[0][2] = 3.0;
  //v[1][0] = 4.0;
  //v[1][1] = 5.0;
  //v[1][2] = 6.0;

  //Matrix m(v);
  //EXPECT_EQ(m.m, 2);
  //EXPECT_EQ(m.n, 3);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);
  //EXPECT_DOUBLE_EQ(m.v[0][1], 2.0);
  //EXPECT_DOUBLE_EQ(m.v[0][2], 3.0);
  //EXPECT_DOUBLE_EQ(m.v[1][0], 4.0);
  //EXPECT_DOUBLE_EQ(m.v[1][1], 5.0);
  //EXPECT_DOUBLE_EQ(m.v[1][2], 6.0);

  //m.transpose_inplace();
  //EXPECT_EQ(m.m, 3);
  //EXPECT_EQ(m.n, 2);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);
  //EXPECT_DOUBLE_EQ(m.v[0][1], 4.0);
  //EXPECT_DOUBLE_EQ(m.v[1][0], 2.0);
  //EXPECT_DOUBLE_EQ(m.v[1][1], 5.0);
  //EXPECT_DOUBLE_EQ(m.v[2][0], 3.0);
  //EXPECT_DOUBLE_EQ(m.v[2][1], 6.0);

  //Matrix n(vector<double>{3.0, 4.0});
  //EXPECT_EQ(n.m, 1);
  //EXPECT_EQ(n.n, 2);
  //EXPECT_DOUBLE_EQ(n.v[0][0], 3.0);
  //EXPECT_DOUBLE_EQ(n.v[0][1], 4.0);

  //Matrix zz(n.transpose());
  //EXPECT_EQ(zz.m, 2);
  //EXPECT_EQ(zz.n, 1);
  //EXPECT_DOUBLE_EQ(zz.v[0][0], 3.0);
  //EXPECT_DOUBLE_EQ(zz.v[1][0], 4.0);

  //n.transpose_inplace();
  //EXPECT_EQ(n.m, 2);
  //EXPECT_EQ(n.n, 1);
  //EXPECT_DOUBLE_EQ(n.v[0][0], 3.0);
  //EXPECT_DOUBLE_EQ(n.v[1][0], 4.0);
//}

//TEST(MatrixTest, eyeTest) {
  //Matrix m;

  //m.eye_inplace(1);
  //EXPECT_EQ(m.m, 1);
  //EXPECT_EQ(m.n, 1);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);

  //m.eye_inplace(2);
  //EXPECT_EQ(m.m, 2);
  //EXPECT_EQ(m.n, 2);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);
  //EXPECT_DOUBLE_EQ(m.v[0][1], 0.0);
  //EXPECT_DOUBLE_EQ(m.v[1][0], 0.0);
  //EXPECT_DOUBLE_EQ(m.v[1][1], 1.0);

  //m.eye_inplace(1);
  //EXPECT_EQ(m.m, 1);
  //EXPECT_EQ(m.n, 1);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);
//}

//TEST(eyeTest, eyeTest) {
  //Matrix m(eye(2));
  //EXPECT_EQ(m.m, 2);
  //EXPECT_EQ(m.n, 2);
  //EXPECT_DOUBLE_EQ(m.v[0][0], 1.0);
  //EXPECT_DOUBLE_EQ(m.v[0][1], 0.0);
  //EXPECT_DOUBLE_EQ(m.v[1][0], 0.0);
  //EXPECT_DOUBLE_EQ(m.v[1][1], 1.0);
  
  //Matrix n(eye(1));
  //EXPECT_EQ(n.m, 1);
  //EXPECT_EQ(n.n, 1);
  //EXPECT_DOUBLE_EQ(n.v[0][0], 1.0);
//}

////TEST(fTest, valueTest) {
  ////vector<double> v = {0.0, 0.0};
  ////EXPECT_DOUBLE_EQ(f(Matrix(v)), 0.0);
  ////v = {1.0, 0.0};
  ////EXPECT_DOUBLE_EQ(f(Matrix(v)), 0.0);
  ////v = {1.0, 1.0};
  ////EXPECT_DOUBLE_EQ(f(Matrix(v)), -1.0);
  ////v = {2.0, 1.0};
  ////EXPECT_DOUBLE_EQ(f(Matrix(v)), -sqrt(2.0));
////}

////TEST(fTest, exceptionTest) {
  ////vector<double> v = {1.0, -1.0};
  ////EXPECT_ANY_THROW(f(Matrix(v)));
////}

//TEST(MatrixTest, operatorPlusMinusTest) {
  //Matrix m(vector<double>{1.0, 2.0});
  //Matrix n(vector<double>{1.0, -0.5});
  //Matrix ans(m + n);
  //EXPECT_EQ(ans.m, 1);
  //EXPECT_EQ(ans.n, 2);
  //EXPECT_DOUBLE_EQ(ans.v[0][0], 2.0);
  //EXPECT_DOUBLE_EQ(ans.v[0][1], 1.5);
  //Matrix ans2(m - n);
  //EXPECT_EQ(ans2.m, 1);
  //EXPECT_EQ(ans2.n, 2);
  //EXPECT_DOUBLE_EQ(ans2.v[0][0], 0.0);
  //EXPECT_DOUBLE_EQ(ans2.v[0][1], 2.5);
//}

//TEST(MatrixTest, operatorTimesTest) {
  //Matrix m(vector<double>{1.0, 2.0});
  //Matrix n(2, 1);
  //n.v[0][0] = 3.0;
  //n.v[1][0] = 4.0;
  //Matrix ans(m * n);
  //EXPECT_EQ(ans.m, 1);
  //EXPECT_EQ(ans.n, 1);
  //EXPECT_DOUBLE_EQ(ans.v[0][0], 11);
  //Matrix ans2(n * m);
  //EXPECT_EQ(ans2.m, 2);
  //EXPECT_EQ(ans2.n, 2);
  //EXPECT_DOUBLE_EQ(ans2.v[0][0], 3);
  //EXPECT_DOUBLE_EQ(ans2.v[0][1], 6);
  //EXPECT_DOUBLE_EQ(ans2.v[1][0], 4);
  //EXPECT_DOUBLE_EQ(ans2.v[1][1], 8);
//}

//TEST(MatrixTest, det2Test) {
  //Matrix m(2, 2);
  //m.v[0][0] = 1.0;
  //m.v[1][1] = 3.0;
  //m.v[1][0] = 2.0;
  //m.v[0][1] = 6.0;
  //EXPECT_DOUBLE_EQ(m.det2(), -9);
//}
