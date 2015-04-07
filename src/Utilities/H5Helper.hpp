/// \file H5Helper.hpp
/// \author jlippuner
/// \since Feb 17, 2015
///
/// \brief
///
///

#ifndef SKYNET_UTILITIES_H5HELPER_HPP_
#define SKYNET_UTILITIES_H5HELPER_HPP_

#include <H5Cpp.h>

template<typename T>
H5::DataType GetH5DataType() {
  if (std::is_same<T, double>::value)
    return H5::DataType(H5::PredType::NATIVE_DOUBLE);
  else if (std::is_same<T, int>::value)
    return H5::DataType(H5::PredType::NATIVE_INT);
  else
    throw std::invalid_argument("Unknown type");
}

#endif // SKYNET_UTILITIES_H5HELPER_HPP_
