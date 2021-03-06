//
// Created by F.Moitzi on 13.01.2022.
//

#ifndef SRC_MADELUNG_COMMON_HPP
#define SRC_MADELUNG_COMMON_HPP

#include "Array3d.hpp"
#include "Matrix.hpp"

#include "NDArray.hpp"

constexpr auto Y0inv = 2.0 * 1.772453850905514;

namespace lsms {

template <typename T>
using matrix = NDArray<T,2>;

template <typename T>
using array3d = NDArray<T, 3>;

}  // namespace lsms

#endif  // SRC_MADELUNG_COMMON_HPP
