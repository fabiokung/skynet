#include <cstdio>
#include <fstream>
#include <regex>
#include <cmath>

#include "Network/NetworkOutput.hpp"

int main(int argc, char** argv) {
  if (argc != 3) {
    printf("Usage: calc_error <h5 file> <NSE file>\n");
    return 1;
  }

  // read NSE
  std::ifstream ifs(argv[2]);
  std::vector<double> nse;

  int Z, N;
  double Y;

  while (true) {
    ifs >> Z;
    ifs >> N;
    ifs >> Y;

    if (ifs.eof())
      break;

    nse.push_back(Y);
  }

  ifs.close();

  // read SkyNet results
  auto out = NetworkOutput::ReadFromFile(std::string(argv[1]));
  auto& ts = out.Times();
  auto& Ys = out.YVsTime();

  if (nse.size() != Ys[0].size()) {
    printf("Mismatch in number of abundances\n");
    printf("nse: %lu, h5: %lu\n", nse.size(), Ys[0].size());
    return 1;
  }

  printf("# %16s  %18s\n", "Time [s]", "Abundance error");
  for (size_t i = 0; i < ts.size(); ++i) {
    double error = 0.0;
    for (size_t j = 0; j < Ys[i].size(); ++j) {
      error = std::max(error, fabs(nse[j] - Ys[i][j]));
    }

    printf("%18.10E  %18.10E\n", ts[i], error);
  }

  return 0;
}
