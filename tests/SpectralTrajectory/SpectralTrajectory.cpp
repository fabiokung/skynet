/// \file SpectralTrajectory.cpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "SpectralTrajectory.hpp"

namespace {

bool failed = false;

void Check(const bool ok, const std::string& what) {
  if (!ok) {
    std::cerr << "FAILED: " << what << std::endl;
    failed = true;
  }
}

double FermiDirac(const double e, const double T, const double eta) {
  return 1.0 / (exp(e / T - eta) + 1.0);
}

// Parse a whole file the way ReadTrajectory does: header lines through
// ParseSpectrumHeaderLine, data rows through ReadSpectrumRow.
SpectrumHeader ParseHeaderBlock(const std::string& text) {
  std::istringstream in(text);
  std::string line;
  SpectrumHeader header;

  while (std::getline(in, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) continue;
    if (line[first] == '#')
      ParseSpectrumHeaderLine(line, &header);
  }

  return header;
}

void TestHeaderParsing() {
  auto header = ParseHeaderBlock(
      "# SPECTRUM_GRID: nbins=3 width=0.5\n"
      "# SPECTRUM_SPECIES: NuE AntiNuE\n"
      "# SPECTRUM_KIND: number_flux\n");

  Check(header.Present, "grid key marks the header present");
  Check(header.NBins == 3 && header.Width == 0.5, "nbins and width are read");
  // nbins intervals -> nbins+1 edges at 0, width, ..., nbins*width
  Check(header.EnergiesMeV.size() == 4, "grid expands to nbins+1 edges");
  Check(header.EnergiesMeV[2] == 1.0, "grid edges are k*width");
  Check(header.EnergiesMeV.back() == 1.5, "grid ends at nbins*width");
  Check(header.Kind == "number_flux", "kind is read");

  // key order doesn't matter
  auto swapped = ParseHeaderBlock("# SPECTRUM_GRID: width=0.25 nbins=400\n");
  Check(swapped.NBins == 400 && swapped.Width == 0.25,
      "nbins and width parse in either order");
  Check(swapped.EnergiesMeV.back() == 100.0, "400 x 0.25 reaches 100 MeV");

  // a file without the grid key is a legacy 18-column trajectory
  auto legacy = ParseHeaderBlock("# just a comment\n# another\n");
  Check(!legacy.Present, "a file with no grid key is not a spectral file");
  Check(legacy.Kind == "occupation", "kind defaults to occupation");
}

void TestHeaderRejection() {
  auto rejects = [] (const std::string& line, const std::string& what) {
    bool threw = false;
    try {
      SpectrumHeader header;
      ParseSpectrumHeaderLine(line, &header);
    } catch (const std::runtime_error&) {
      threw = true;
    }
    Check(threw, what);
  };

  rejects("# SPECTRUM_KIND: bogus", "an unknown SPECTRUM_KIND is rejected");
  rejects("# SPECTRUM_KIND: occupation number_flux",
      "a SPECTRUM_KIND with two values is rejected");
  rejects("# SPECTRUM_SPECIES: NuE NuMu",
      "a non electron-flavor SPECTRUM_SPECIES is rejected");
  rejects("# SPECTRUM_SPECIES: NuE", "a single-species header is rejected");
  rejects("# SPECTRUM_GRID: nbins=400 0.25",
      "a SPECTRUM_GRID token without key=value is rejected");
  rejects("# SPECTRUM_GRID: nbins=400 step=0.25",
      "an unknown SPECTRUM_GRID key is rejected");
  rejects("# SPECTRUM_GRID: width=0.25", "a grid with no nbins is rejected");
  rejects("# SPECTRUM_GRID: nbins=400", "a grid with no width is rejected");
  rejects("# SPECTRUM_GRID: nbins=0 width=0.25",
      "a grid with non-positive nbins is rejected");
}

void TestRowReading() {
  SpectrumHeader header;
  header.Present = true;
  header.EnergiesMeV = { 0.0, 0.5, 1.0 };

  // a full row: 3 nu_e columns then 3 nu_e-bar columns
  std::istringstream full("1.0 2.0 3.0 4.0 5.0 6.0");
  auto row = ReadSpectrumRow(full, header);
  Check(row.size() == 2 && row[0].size() == 3,
      "a spectrum row is [species][energy]");
  Check(row[0][2] == 3.0 && row[1][0] == 4.0,
      "nu_e columns come before nu_e-bar columns");

  // a short row has fewer than 2*nE columns
  bool threw = false;
  try {
    std::istringstream shortRow("1.0 2.0 3.0 4.0 5.0");
    ReadSpectrumRow(shortRow, header);
  } catch (const std::runtime_error&) {
    threw = true;
  }
  Check(threw, "a row with too few spectrum columns is rejected");
}

// A number_flux column holds n(E) = E^2 f(E); the driver has to recover f(E).
// This round-trip is what lets the occupation and number_flux trajectories give
// the same yields.
void TestNumberFluxRoundTrip() {
  const double T = 2.66, eta = 2.1;
  std::vector<double> energies;
  for (int i = 0; i <= 200; ++i)
    energies.push_back(i * 0.5);

  std::vector<std::vector<std::vector<double>>> spectra(1,
      std::vector<std::vector<double>>(2, std::vector<double>(energies.size())));
  for (std::size_t e = 0; e < energies.size(); ++e) {
    const double f = FermiDirac(energies[e], T, eta);
    spectra[0][0][e] = f * energies[e] * energies[e];
    spectra[0][1][e] = f * energies[e] * energies[e];
  }

  const std::size_t dropped = ConvertNumberFluxToOccupation(&energies, &spectra);
  Check(dropped == 1, "the E = 0 bin is dropped");
  Check(energies.front() == 0.5, "the recovered grid starts at the first "
      "non-zero energy");

  double maxRel = 0.0;
  for (std::size_t e = 0; e < energies.size(); ++e) {
    const double exact = FermiDirac(energies[e], T, eta);
    maxRel = std::max(maxRel, std::fabs(spectra[0][0][e] / exact - 1.0));
  }
  Check(maxRel < 1.0e-12, "number_flux inverts back to the occupation number");

  // a grid that does not start at zero drops nothing
  std::vector<double> nonzero = { 0.5, 1.0, 1.5 };
  std::vector<std::vector<std::vector<double>>> tiny(1,
      std::vector<std::vector<double>>(2, std::vector<double>{ 1.0, 2.0, 3.0 }));
  Check(ConvertNumberFluxToOccupation(&nonzero, &tiny) == 0,
      "a grid without an E = 0 bin drops nothing");
}

} // namespace

int main(int, char**) {
  TestHeaderParsing();
  TestHeaderRejection();
  TestRowReading();
  TestNumberFluxRoundTrip();

  return failed ? EXIT_FAILURE : EXIT_SUCCESS;
}
