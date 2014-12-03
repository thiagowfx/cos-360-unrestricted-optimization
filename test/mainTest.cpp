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
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, StandardConstructor2) {
  Matrix<double> m(2, 1, 10.0);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 10.0);
  EXPECT_EQ(m.get(2,1), 10.0);
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, AliasConstructor) {
  Matrix<double> n(2, 1, 2.0);
  Matrix<double> m(n);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 2.0);
  EXPECT_EQ(m.get(2,1), 2.0);
  EXPECT_THROW(m.get(1,2), std::exception);
}

TEST(MatrixTest, VectorVectorConstructor) {
  vector<vector<double> > v;
  v.push_back(vector<double>(1, 2.0));
  v.push_back(vector<double>(1, 2.0));
  Matrix<double> m(v);
  EXPECT_EQ(m.getRows(), 2);
  EXPECT_EQ(m.getCols(), 1);
  EXPECT_EQ(m.get(1,1), 2.0);
  EXPECT_EQ(m.get(2,1), 2.0);
}

TEST(MatrixTest, get) {
  Matrix<double> m(2,1);
  EXPECT_THROW(m.get(1,2), std::exception);
  EXPECT_THROW(m.get(2,2), std::exception);
}

TEST(MatrixTest, get2) {
  Matrix<double> m(1,2);
  EXPECT_THROW(m.get(2,1), std::exception);
  EXPECT_THROW(m.get(2,2), std::exception);
}

TEST(MatrixTest, set) {
  Matrix<double> m(2,1);
  m.set(1,1,1.0);
  m.set(2,1,2.0);
  EXPECT_EQ(m.get(1,1), 1.0);
  EXPECT_EQ(m.get(2,1), 2.0);
}

TEST(MatrixTest, operatorPlusMinus) {
  Matrix<double> m(2, 1, 2.0);
  Matrix<double> n(2, 1, 3.0);
  Matrix<double> s = m + n;
  Matrix<double> d = m - n;
  EXPECT_DOUBLE_EQ(s.get(1,1), 5.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 5.0);
  EXPECT_DOUBLE_EQ(d.get(1,1), -1.0);
  EXPECT_DOUBLE_EQ(d.get(2,1), -1.0);
}

TEST(MatrixTest, operatorMultScalar) {
  Matrix<double> m(2,1, 2.0);
  Matrix<double> s = m * 2;
  EXPECT_DOUBLE_EQ(s.get(1,1), 4.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 4.0);
}

TEST(MatrixTest, operatorDivScalar) {
  Matrix<double> m(2,1, 2.0);
  Matrix<double> s = m / 2;
  EXPECT_DOUBLE_EQ(s.get(1,1), 1.0);
  EXPECT_DOUBLE_EQ(s.get(2,1), 1.0);
}

TEST(MatrixTest, isVector) {
  Matrix<double> a(1,1);
  Matrix<double> b(2,1);
  Matrix<double> c(1,2);
  Matrix<double> d(2,2);
  EXPECT_TRUE(a.isVector());
  EXPECT_TRUE(b.isVector());
  EXPECT_TRUE(c.isVector());
  EXPECT_FALSE(d.isVector());
}

TEST(MatrixTest, x1x2) {
  Matrix<double> m(2, 1, 3.0);
  EXPECT_DOUBLE_EQ(m.x1(), 3.0);
  EXPECT_DOUBLE_EQ(m.x2(), 3.0);

  Matrix <double> n(1, 2);
  EXPECT_THROW(n.x1(), std::exception);
  EXPECT_THROW(n.x2(), std::exception);

  Matrix <double> o(1, 1);
  EXPECT_THROW(o.x1(), std::exception);
  EXPECT_THROW(o.x2(), std::exception);

  Matrix<double> p(2, 2);
  EXPECT_THROW(p.x1(), std::exception);
  EXPECT_THROW(p.x2(), std::exception);
}

//TEST(MatrixTest, operatorMultScalar2) {
  //Matrix<double> m(2,1, 2.0);
  //Matrix<double> s = 2 * m;
  //EXPECT_DOUBLE_EQ(s.get(1,1), 4.0);
  //EXPECT_DOUBLE_EQ(s.get(2,1), 4.0);
//}

