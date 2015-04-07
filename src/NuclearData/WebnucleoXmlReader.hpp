/// \file WebnucleoXmlReader.hpp
/// \author jlippuner
/// \since Jul 3, 2014
///
/// \brief
///
///

#ifndef SKYNET_NUCLEARDATA_WEBNUCLEOXMLREADER_HPP_
#define SKYNET_NUCLEARDATA_WEBNUCLEOXMLREADER_HPP_

#include <string>
#include <unordered_map>

#include "tinyxml2/tinyxml2.h"

#include "NuclearData/Nuclide.hpp"

class WebnucleoXmlReader {
public:
  static std::unordered_map<std::string, Nuclide> GetAllNuclides(
      const std::string& pathToWebnucleoXmlFile);

private:
  // don't allow this class to be constructed
  WebnucleoXmlReader();

  static Nuclide ReadNuclide(const tinyxml2::XMLNode& node);
};

#endif // SKYNET_NUCLEARDATA_WEBNUCLEOXMLREADER_HPP_
