// Copyright (C) 2020 Anders Logg, Anton J Olsson
// Licensed under the MIT License

#include <cassert>
#include <cstddef>

#ifndef DTCC_UTILS_H
#define DTCC_UTILS_H

namespace DTCC_BUILDER
{

/// Index utility used by native raster sampling.
class Utils
{
public:
  /// Crop integer x to interval [0, n - k). Requires care
  /// due to involving both signed and unsigned integers.
  ///
  /// @param x Signed integer
  /// @param n Unsigned integer (length of array)
  /// @param k Unsigned integer (maring at end of array)
  /// @return Unsigned integer within specified range
  static size_t crop(long int x, size_t n, size_t k = 0)
  {
    assert(n > 0);
    assert(k < n);
    return x < 0 ? 0 : (x + k >= n ? n - 1 - k : x);
  }


};
} // namespace DTCC_BUILDER

#endif
