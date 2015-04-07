/// \file WebnucleoXmlReader.cpp
/// \author jlippuner
/// \since Jul 3, 2014
///
/// \brief
///
///

#include <NuclearData/WebnucleoXmlReader.hpp>
#include <Utilities/File.hpp>

#include <array>


std::unordered_map<std::string, Nuclide> WebnucleoXmlReader::GetAllNuclides(
    const std::string& pathToWebnucleoXmlFile) {
  if (!File::Exists(pathToWebnucleoXmlFile)) {
    throw std::invalid_argument("The webnucleo XML file '"
        + pathToWebnucleoXmlFile + "' does not exist.");
  }

  tinyxml2::XMLDocument doc;
  doc.LoadFile(pathToWebnucleoXmlFile.c_str());

  std::unordered_map<std::string, Nuclide> nuclides;

  auto root = doc.RootElement();
  for (auto node = root->FirstChildElement("nuclide"); node != nullptr;
      node = node->NextSiblingElement("nuclide")) {
    auto nuclide = ReadNuclide(*node);
    auto pair = nuclides.insert({nuclide.Name(), nuclide});
    if (!pair.second)
      throw std::invalid_argument("Duplicate nuclide with name "
          + nuclide.Name());
  }

  return nuclides;
}

Nuclide WebnucleoXmlReader::ReadNuclide(const tinyxml2::XMLNode& node) {
  // read Z
  auto pZ = node.FirstChildElement("z");
  if (!pZ)
    throw std::invalid_argument("No z element");
  int z = std::stoi(pZ->GetText());

  // read A
  auto pA = node.FirstChildElement("a");
  if (!pA)
    throw std::invalid_argument("No a element");
  int a = std::stoi(pA->GetText());

  // read source
//  auto pSource = node.FirstChildElement("source");
//  if (!pSource)
//    throw std::invalid_argument("No source element");
//  std::string source = pSource->GetText();

  // read mass excess
  auto pMassExcess = node.FirstChildElement("mass_excess");
  if (!pMassExcess)
    throw std::invalid_argument("No mass_excess element");
  double massExcess = std::stod(pMassExcess->GetText());

  // read spin
  auto pSpin = node.FirstChildElement("spin");
  if (!pSpin)
    throw std::invalid_argument("No spin element");
  double spin = std::stod(pSpin->GetText());

  // read partition function
  unsigned int counter = 0;
  std::array<double, Nuclide::NumPartitionFunctionEntries> partfLog10;

  auto pPartfTable = node.FirstChildElement("partf_table");
  if (!pPartfTable)
    throw std::invalid_argument("No partf_table element");
  // loop over points
  for (auto point = pPartfTable->FirstChildElement("point"); point != nullptr;
      point = point->NextSiblingElement("point")) {
    auto pT9 = point->FirstChildElement("t9");
    if (!pT9)
      throw std::invalid_argument("No t9 element");
    double t9 = std::stod(pT9->GetText());

    auto pLog10 = point->FirstChildElement("log10_partf");
    if (!pLog10)
      throw std::invalid_argument("No log10_partf element");
    double log10 = std::stod(pLog10->GetText());

    // check T9 is as expected
    if (t9 != Nuclide::PartitionFunctionT9[counter])
      throw std::invalid_argument("Unexpected value for T9 in partition "
          "function");

    if (counter >= Nuclide::NumPartitionFunctionEntries)
      throw std::invalid_argument("Too many partition function entries");

    partfLog10[counter++] = log10;
  }

  if (counter != Nuclide::NumPartitionFunctionEntries) {
    throw std::invalid_argument("Wrong number of partition function entries");
  }

//  std::cout << "Z = " << z << ", A = " << a
//      //<< ", source = '" << source
//      << "', mass excess = " << massExcess
//      //<< ", spin = " << spin
//      << std::endl;
//
//  for (unsigned int i = 0; i < Nuclide::PartitionFunctionT9.size(); ++i) {
//    std::cout << Nuclide::PartitionFunctionT9[i] << " "
//        << partfLog10[i] << std::endl;
//  }

  return Nuclide(z, a, massExcess, spin, partfLog10, "");
}

