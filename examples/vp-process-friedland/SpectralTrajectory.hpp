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
//   # SPECTRUM_GRID: nbins=400 width=0.25
//   # SPECTRUM_SPECIES: NuE AntiNuE
//   # SPECTRUM_KIND: occupation | number_flux
//
// The grid is uniform: nbins intervals of `width` MeV starting at 0, sampled at
// the nbins+1 edges 0, width, ..., nbins*width (one column per edge, shared by
// both species). Only the shape matters; the scale is set by the L(t) the
// driver builds, the same one the (T, eta) path uses.
//
// A grid replaces the analytic pinched Fermi-Dirac shape outright, so whatever
// (T, eta) the driver would otherwise use goes unread. The grid is also
// frame-fixed -- nothing downstream can rescale it per radius -- so a GR
// blueshift has to be baked into the file, spectra and L alike, rather than
// applied on the way in.
struct SpectrumHeader {
  std::vector<double> EnergiesMeV;
  int NBins = 0;
  double Width = 0.0;
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

  if (ParseHeaderKey(line, "SPECTRUM_GRID", &tokens)) {
    int nbins = 0;
    double width = 0.0;
    for (const auto& t : tokens) {
      const auto eq = t.find('=');
      if (eq == std::string::npos)
        throw std::runtime_error("SPECTRUM_GRID tokens are key=value, got " + t);
      const std::string key = t.substr(0, eq), val = t.substr(eq + 1);
      if (key == "nbins") nbins = std::stoi(val);
      else if (key == "width") width = std::stod(val);
      else throw std::runtime_error("SPECTRUM_GRID: unknown key " + key);
    }
    if (nbins <= 0)
      throw std::runtime_error("SPECTRUM_GRID needs a positive nbins");
    if (!(width > 0.0))
      throw std::runtime_error("SPECTRUM_GRID needs a positive width");

    pHeader->Present = true;
    pHeader->NBins = nbins;
    pHeader->Width = width;
    pHeader->EnergiesMeV.clear();
    for (int k = 0; k <= nbins; ++k)
      pHeader->EnergiesMeV.push_back(k * width);
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
