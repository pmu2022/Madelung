
#include <gtest/gtest.h>

#include "NDArray.hpp"

#include <cmath>
#include <complex>


TEST(MultiArrayTest, TestAccess) {


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


  lsms::NDArray<double, 2> ccc(c);
  lsms::NDArray<double, 2> cccc = c;
  lsms::NDArray<double, 2> cm(std::move(c));


}

