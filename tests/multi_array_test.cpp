
#include <gtest/gtest.h>

#include "NDArray.hpp"

#include <cmath>
#include <complex>

template<typename T>
using matrix = lsms::NDArray<T, 2>;



TEST(MultiArrayTest, Constructor) {

  lsms::NDArray<double, 1> c1(2);
  lsms::NDArray<double, 2> c2(2, 4);
  lsms::NDArray<double, 3> c3(2, 4, 6);

}

TEST(MultiArrayTest, CopyConstructor) {

  lsms::NDArray<double, 2> c2(3, 4);
  c2 = 3.0;

  auto c1(c2);
  EXPECT_EQ(c1.size(), c2.size());

  for (auto i = 0; i < 3; i++) {
    for (auto j = 0; j < 4; j++) {
      EXPECT_NEAR(c1(i, j), c2(i, j), 1e-12);
    }
  }

}


TEST(MultiArrayTest, CopyAssigment) {


  lsms::NDArray<double, 2> c1(2, 5);
  lsms::NDArray<double, 2> c2(2, 3);

  c1 = 1.0;
  c2 = 3.0;

  // Copy-Assigment constructor
  EXPECT_NE(c1.size(), c2.size());
  c1 = c2;
  EXPECT_EQ(c1.size(), c2.size());

  for (auto i = 0; i < 2; i++) {
    for (auto j = 0; j < 3; j++) {
      EXPECT_NEAR(c1(i, j), c2(i, j), 1e-12);
    }
  }

}

TEST(MultiArrayTest, MoveAssigment) {


  lsms::NDArray<double, 2> c1(2, 5);
  lsms::NDArray<double, 2> c2(2, 3);
  lsms::NDArray<double, 2> c3(2, 3);

  c1 = 1.0;
  c2 = 3.0;
  c3 = c2;

  EXPECT_NE(c1.size(), c2.size());
  c1 = std::move(c2);
  EXPECT_EQ(c1.size(), c2.size());

  for (auto i = 0; i < 2; i++) {
    for (auto j = 0; j < 3; j++) {
      EXPECT_NEAR(c3(i, j), c1(i, j), 1e-12);
    }
  }

}

TEST(MultiArrayTest, MoveConstructor) {

  lsms::NDArray<double, 2> c1(2, 5);
  lsms::NDArray<double, 2> c2(2, 3);

  c1 = 1.0;
  c2 = 3.0;

  c1 = c2;
  EXPECT_EQ(c1.size(), c1.size());

  auto c3(std::move(c2));
  EXPECT_EQ(c3.size(), c1.size());

  for (auto i = 0; i < 2; i++) {
    for (auto j = 0; j < 3; j++) {
      EXPECT_NEAR(c3(i, j), c1(i, j), 1e-12);
    }
  }

}


TEST(MultiArrayTest, RandomAccessAndOrder) {

  constexpr auto cols = 4;
  constexpr auto rows = 3;

  lsms::NDArray<double, 2> c(cols, rows);

  constexpr auto ref_size = cols * rows;

  EXPECT_EQ(ref_size, c.size());

  for (auto col = 0; col < cols; col++) {
    for (auto row = 0; row < rows; row++) {
      c(col, row) = col * 10 + row;
    }
  }

  EXPECT_EQ(c(1, 0), 10);
  EXPECT_EQ(c(1, 0), c[1]);

  EXPECT_EQ(c(2, 1), 21);
  EXPECT_EQ(c(2, 1), c[6]);

  lsms::NDArray<double, 3> cc;
  auto shape = cc.shape();
  EXPECT_EQ(shape[0], 0);
  EXPECT_EQ(shape[1], 0);
  EXPECT_EQ(shape[2], 0);

}

TEST(MultiArrayTest, ListInitializer1) {

  lsms::NDArray<double, 1> c {1, 3, 4, 5, 7.7, 8.7};

  EXPECT_NEAR(c(0), 1, 1e-12);
  EXPECT_NEAR(c(1), 3, 1e-12);
  EXPECT_NEAR(c(5), 8.7, 1e-12);

  EXPECT_EQ(c.size(), 6);

}

TEST(MultiArrayTest, ListInitializer2) {

  lsms::NDArray<double, 2> c2 { {1, 2}, {3, 4}};

  EXPECT_NEAR(c2(0, 0), 1, 1e-12);
  EXPECT_NEAR(c2(0, 1), 2, 1e-12);
  EXPECT_NEAR(c2(1, 0), 3, 1e-12);
  EXPECT_NEAR(c2(1, 1), 4, 1e-12);

  EXPECT_EQ(c2.size(), 4);

}
