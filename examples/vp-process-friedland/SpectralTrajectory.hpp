/// \file SpectralTrajectory.hpp
/// \author Fabio Kung
/// \since Jul 2026
///
/// \brief
///
///

#ifndef VP_PROCESS_SPECTRALTRAJECTORY_HPP_
#define VP_PROCESS_SPECTRALTRAJECTORY_HPP_

#include <istream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Optional tabulated neutrino spectra appended to the 18-column trajectory,
// one column per energy bin (nu_e first, then nu_e-bar), announced by header
// keys:
//
//   # SPECTRUM_ENERGY_GRID_MEV: 0 0.25 ... 100
//   # SPECTRUM_SPECIES: NuE AntiNuE
//   # SPECTRUM_KIND: occupation | number_flux
//
// Only the shape matters; the scale is set by the L(t) the driver builds, the
// same one the (T, eta) path uses.
struct SpectrumHeader {
  std::vector<double> EnergiesMeV;
  std::string Kind = "occupation";
  bool Present = false;
};

// Reads "# KEY: v1 v2 ..." out of a comment line
inline bool ParseHeaderKey(const std::string& line, const std::string& key,
    std::vector<std::string> * const pValues) {
  const auto at = line.find(key + ":");
  if (at == std::string::npos) return false;

  std::istringstream iss(line.substr(at + key.size() + 1));
  std::string token;
  pValues->clear();
  while (iss >> token) pValues->push_back(token);

  return true;
}

// Applies one '#' comment line to the header, throwing on a malformed value.
inline void ParseSpectrumHeaderLine(const std::string& line,
    SpectrumHeader * const pHeader) {
  std::vector<std::string> tokens;

  if (ParseHeaderKey(line, "SPECTRUM_ENERGY_GRID_MEV", &tokens)) {
    pHeader->Present = true;
    pHeader->EnergiesMeV.clear();
    for (const auto& t : tokens)
      pHeader->EnergiesMeV.push_back(std::stod(t));
  } else if (ParseHeaderKey(line, "SPECTRUM_KIND", &tokens)) {
    if (tokens.size() != 1)
      throw std::runtime_error("SPECTRUM_KIND takes a single value");
    pHeader->Kind = tokens[0];
    if (pHeader->Kind != "occupation" && pHeader->Kind != "number_flux")
      throw std::runtime_error("SPECTRUM_KIND must be occupation or "
          "number_flux, got " + pHeader->Kind);
  } else if (ParseHeaderKey(line, "SPECTRUM_SPECIES", &tokens)) {
    if (tokens != std::vector<std::string>{"NuE", "AntiNuE"})
      throw std::runtime_error("SPECTRUM_SPECIES must be \"NuE AntiNuE\"; "
          "this driver only evolves electron-flavor neutrinos");
  }
}

// Reads the 2*nE spectrum columns of one data row from a stream positioned
// after the 18 trajectory columns; returns [species][energy], throws short.
inline std::vector<std::vector<double>> ReadSpectrumRow(std::istream& iss,
    const SpectrumHeader& header) {
  const std::size_t nE = header.EnergiesMeV.size();
  std::vector<std::vector<double>> row(2, std::vector<double>(nE));

  for (std::size_t s = 0; s < 2; ++s) {
    for (std::size_t e = 0; e < nE; ++e) {
      if (!(iss >> row[s][e]))
        throw std::runtime_error("Trajectory row has fewer spectrum columns "
            "than the " + std::to_string(2 * nE) + " the energy grid calls "
            "for");
    }
  }

  return row;
}

// number_flux n(E) = E^2 f(E): drop the E=0 bin (n(0)=0 carries no shape and
// can't be inverted) and divide out E^2 in place to recover the occupation
// number. Returns the number of leading bins dropped.
inline std::size_t ConvertNumberFluxToOccupation(
    std::vector<double> * const pEnergies,
    std::vector<std::vector<std::vector<double>>> * const pSpectra) {
  const std::size_t first = pEnergies->front() == 0.0 ? 1 : 0;

  pEnergies->erase(pEnergies->begin(), pEnergies->begin() + first);

  for (auto& row : *pSpectra) {
    for (auto& perSpecies : row) {
      perSpecies.erase(perSpecies.begin(), perSpecies.begin() + first);
      for (std::size_t e = 0; e < perSpecies.size(); ++e)
        perSpecies[e] /= (*pEnergies)[e] * (*pEnergies)[e];
    }
  }

  return first;
}

#endif // VP_PROCESS_SPECTRALTRAJECTORY_HPP_
