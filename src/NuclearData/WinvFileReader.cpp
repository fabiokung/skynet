/// \file WinvFileReader.cpp
/// \author jlippuner
/// \since Jul 28, 2014
///
/// \brief
///
///

#include <Utilities/File.hpp>
#include "NuclearData/WinvFileReader.hpp"

#include <fstream>

#include "Utilities/FloatingPointComparison.hpp"

std::unordered_map<std::string, Nuclide> WinvFileReader::GetAllNuclides(
      const std::string& pathToWinvFile) {
  if (!File::Exists(pathToWinvFile)) {
    throw std::invalid_argument("The Winv file '"
        + pathToWinvFile + "' does not exist.");
  }

  std::ifstream ifs(pathToWinvFile.c_str(), std::ifstream::in);
  std::string line;

  // the first line may or may not contain the total number of nuclides, so just
  // ignore it
  std::getline(ifs, line);

  // The next line contains the temperatures of the partition function tables.
  // There are 24 entries and each entry is 3 characters long, with no
  // separating character in between. Each entry is converted to an integer and
  // then multiplied by a multiplier, which starts at 0.01. The entries are
  // strictly monotonically increasing and if an entry would be smaller than the
  // previous one, the multiplier is multiplied by 10.
  std::getline(ifs, line);

  double multiplier = 0.01;
  double prevT9 = 0.0;

  for (int i = 0; i < 24; ++i) {
    std::string entry = line.substr(i * 3, 3);
    double T9 = (double)std::stoi(entry) * multiplier;

    if (T9 <= prevT9) {
      T9 *= 10.0;
      multiplier *= 10.0;
    }

    if (T9 <= prevT9)
      throw std::invalid_argument("T9 is not monotonically increasing in winv "
          "file");

    prevT9 = T9;

    // check that this is the expected value for T9
    if (fcmp(T9, Nuclide::PartitionFunctionT9[i]) != 0)
      throw std::invalid_argument("Unexpected T9 value " + std::to_string(T9)
          + " instead of " + std::to_string(Nuclide::PartitionFunctionT9[i])
          + " in winv file");
  }

  // now read the names
  int numNuclides = 0;
  std::getline(ifs, line);

  while (line.size() == 5) {
    ++numNuclides;
    std::getline(ifs, line);
  }

  //printf("reading %i nuclides...\n", numNuclides);
  std::unordered_map<std::string, Nuclide> nuclidesMap;

  for (int i = 0; i < numNuclides; ++i) {
    // read first line
    char name[5];
    double A;
    int Z;
    int N;
    double spin;
    double massExcess;
    if (sscanf(line.c_str(), "%5s%lf%i%i%lf%lf", name, &A, &Z, &N, &spin,
        &massExcess) != 6)
      throw std::invalid_argument("Could not parse first line of nuclide "
          "record: " + line);

    if ((double)(Z + N) != A)
      throw std::invalid_argument("Z + N != A");

    // the next 3 lines contain the partition functions (plain values, not log)
    std::array<double, 24> partitionFunctionLog10;

    for (int j = 0; j < 3; ++j) {
      std::getline(ifs, line);
      double * p = partitionFunctionLog10.data() + j * 8;

      if (sscanf(line.c_str(), "%lf%lf%lf%lf%lf%lf%lf%lf", p + 0, p + 1, p + 2,
          p + 3, p + 4, p + 5, p + 6, p + 7) != 8)
        throw std::invalid_argument("Expected 8 partition function values");
    }

    // take log
    for (unsigned int j = 0; j < partitionFunctionLog10.size(); ++j)
      partitionFunctionLog10[j] = log10(partitionFunctionLog10[j]);

    // create nuclide
    Nuclide nuclide(Z, Z + N, massExcess, spin, partitionFunctionLog10, "");
    nuclidesMap.insert({nuclide.Name(), nuclide});

    // read header line for next nuclide
    std::getline(ifs, line);
  }

  ifs.close();

  return nuclidesMap;
}

