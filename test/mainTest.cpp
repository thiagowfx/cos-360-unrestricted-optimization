#include "gtest/gtest.h"
#include "lib.hpp"
using namespace std;

TEST(MatrixTest, EmptyConstructor) {
  Matrix m;
  EXPECT_EQ(m.getRows(), 0);
  EXPECT_EQ(m.getCols(), 0);
}

TEST(MatrixTest, StandardConstructor) {
  Matrix m(2, 1);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, StandardConstructor2) {
  Matrix m(2, 1, 10.0);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 10.0);
  EXPECT_EQ(m.get(2,1), 10.0);
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, AliasConstructor) {
  Matrix n(2, 1, 2.0);
  Matrix m(n);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 2.0);
  EXPECT_EQ(m.get(2,1), 2.0);
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, VectorConstructor) {
  vector<double> v;
  v.push_back(1.0);
  v.push_back(2.0);
  v.push_back(3.0);
  Matrix m(v);
  EXPECT_EQ(m.getRows(), 3);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_DOUBLE_EQ(m.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(m.get(2,1), 2.0);
  EXPECT_DOUBLE_EQ(m.get(3,1), 3.0);
}

TEST(MatrixTest, VectorVectorConstructor) {
  vector<vector<double> > v;
  v.push_back(vector<double>(1, 2.0));
  v.push_back(vector<double>(1, 2.0));
  Matrix m(v);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 2.0);
  EXPECT_EQ(m.get(2,1), 2.0);
}

TEST(MatrixTest, get) {
  Matrix m(2,1);
  EXPECT_THROW(m.get(1,2), std::exception);
  EXPECT_THROW(m.get(2,2), std::exception);
}

TEST(MatrixTest, get2) {
  Matrix m(1,2);
  EXPECT_THROW(m.get(2,1), std::exception);
  EXPECT_THROW(m.get(2,2), std::exception);
}

TEST(MatrixTest, get3) {
  Matrix m(2,2);
  m.set(1,1,1.0);
  m.set(2,1,2.0);
  m.set(1,2,3.0);
  m.set(2,2,4.0);
  EXPECT_DOUBLE_EQ(m.get(1), 1.0);
  EXPECT_DOUBLE_EQ(m.get(2), 2.0);
  EXPECT_DOUBLE_EQ(m.get(3), 3.0);
  EXPECT_DOUBLE_EQ(m.get(4), 4.0);
}

TEST(MatrixTest, set) {
  Matrix m(2,1);
  m.set(1,1,1.0);
  m.set(2,1,2.0);
  EXPECT_EQ(m.get(1,1), 1.0);
  EXPECT_EQ(m.get(2,1), 2.0);
}

TEST(MatrixTest, operatorPlusMinus) {
  Matrix m(2, 1, 2.0);
  Matrix n(2, 1, 3.0);
  Matrix s = m + n;
  Matrix d = m - n;
  EXPECT_DOUBLE_EQ(s.get(1,1), 5.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 5.0);
  EXPECT_DOUBLE_EQ(d.get(1,1), -1.0);
  EXPECT_DOUBLE_EQ(d.get(2,1), -1.0);
}

TEST(MatrixTest, operatorMultScalar) {
  Matrix m(2,1, 2.0);
  Matrix s = m * 2;
  EXPECT_DOUBLE_EQ(s.get(1,1), 4.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 4.0);
}

TEST(MatrixTest, operatorDivScalar) {
  Matrix m(2,1, 2.0);
  Matrix s = m / 2;
  EXPECT_DOUBLE_EQ(s.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 1.0);
}

TEST(MatrixTest, isVector) {
  Matrix a(1,1);
  Matrix b(2,1);
  Matrix c(1,2);
  Matrix d(2,2);
  EXPECT_TRUE(a.isVector());
  EXPECT_TRUE(b.isVector());
  EXPECT_TRUE(c.isVector());
  EXPECT_FALSE(d.isVector());
}

TEST(MatrixTest, x) {
  Matrix m(1,1,4.0);
  EXPECT_DOUBLE_EQ(m.x(), 4.0);

  Matrix n(2,1);
  Matrix o(1,2);
  Matrix p(2,2);
  EXPECT_THROW(n.x(), std::exception);
  EXPECT_THROW(o.x(), std::exception);
  EXPECT_THROW(p.x(), std::exception);
}

TEST(MatrixTest, x1x2) {
  Matrix m(2, 1, 3.0);
  EXPECT_DOUBLE_EQ(m.x1(), 3.0);
  EXPECT_DOUBLE_EQ(m.x2(), 3.0);

  Matrix  n(1, 2);
  EXPECT_THROW(n.x1(), std::exception);
  EXPECT_THROW(n.x2(), std::exception);

  Matrix  o(1, 1);
  EXPECT_THROW(o.x1(), std::exception);
  EXPECT_THROW(o.x2(), std::exception);

  Matrix p(2, 2);
  EXPECT_THROW(p.x1(), std::exception);
  EXPECT_THROW(p.x2(), std::exception);
}

TEST(MatrixTest, length) {
  Matrix a(1,1);
  Matrix b(2,1);
  Matrix c(1,2);
  Matrix d(2,2);
  EXPECT_EQ(a.length(), 1);
  EXPECT_EQ(b.length(), 2);
  EXPECT_EQ(c.length(), 2);
  EXPECT_EQ(d.length(), 4);
}

TEST(MatrixTest, mod) {
  Matrix m(2,1);
  m.set(1, 1, 3.0);
  m.set(2, 1, 4.0);
  EXPECT_DOUBLE_EQ(m.mod(), 5.0);
}

TEST(MatrixTest, operatorMultMatrix) {
  Matrix a(2,2,2.0);
  a.set(1,1,1.0);
  Matrix b(2,1,3.0);
  b.set(1,1,4.0);
  Matrix c = a * b;
  EXPECT_EQ(c.getRows(), 2);
  EXPECT_EQ(c.getCols(), 1);
  EXPECT_DOUBLE_EQ(c.get(1,1), 10.0);
  EXPECT_DOUBLE_EQ(c.get(2,1), 14.0);
}

TEST(MatrixTest, transpose) {
  Matrix m(2,1);
  m.set(1,1,1.0);
  m.set(2,1,2.0);
  Matrix t = m.transpose();
  EXPECT_EQ(t.getRows(), 1);
  EXPECT_EQ(t.getCols(), 2);
  EXPECT_DOUBLE_EQ(t.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(t.get(1,2), 2.0);
  EXPECT_THROW(t.get(2,1), std::exception);
  EXPECT_THROW(t.get(2,2), std::exception);
}

TEST(MatrixTest, t) {
  Matrix m(2,1);
  m.set(1,1,1.0);
  m.set(2,1,2.0);
  Matrix t = m.t();
  EXPECT_EQ(t.getRows(), 1);
  EXPECT_EQ(t.getCols(), 2);
  EXPECT_DOUBLE_EQ(t.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(t.get(1,2), 2.0);
  EXPECT_THROW(t.get(2,1), std::exception);
  EXPECT_THROW(t.get(2,2), std::exception);
}

TEST(MatrixTest, det2) {
  Matrix m(2,2,2.0);
  m.set(2,1,3.0);
  m.set(1,2,3.0);
  EXPECT_DOUBLE_EQ(m.det2(), -5);
}

TEST(MatrixTest, operatorMultScalar2) {
  Matrix m(2,1, 2.0);
  Matrix s = 2 * m;
  EXPECT_DOUBLE_EQ(s.get(1,1), 4.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 4.0);
}

TEST(eyeTest, eyeTest) {
  Matrix m = eye(2);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 2);
  EXPECT_DOUBLE_EQ(m.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(m.get(2,2), 1.0);
  EXPECT_DOUBLE_EQ(m.get(1,2), 0.0);
  EXPECT_DOUBLE_EQ(m.get(2,1), 0.0);
}

