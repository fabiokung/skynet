/// \file BuildInfo.hpp
/// \author jlippuner
/// \since Aug 26, 2014
///
/// \brief
///
///

#ifndef SKYNET_BUILDINFO_HPP_
#define SKYNET_BUILDINFO_HPP_

#include <string>

namespace BuildInfo {

extern const std::string GitHash;

extern const bool GitClean;

extern const std::string BuildDateTime;

extern const std::string MatrixSolver;

} // namespace BuildInfo

extern const std::string SkyNetRoot;

#endif // SKYNET_BUILDINFO_HPP_
