/// \file NeutrinoHistoryBlackBody.hpp
/// \author jlippuner
/// \since Dec 22, 2015
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYBLACKBODY_HPP_
#define SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYBLACKBODY_HPP_

#include <cassert>

#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class NeutrinoHistoryBlackBody : public NeutrinoHistory {
public:
  static NeutrinoHistoryBlackBody CreateConstant(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& luminosities,
      const bool interpLogSpace) {
    return DoCreateConstant(times, radii, species, T9s, etas, luminosities,
        interpLogSpace);
  }

  static NeutrinoHistoryBlackBody CreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, radii, species, T9s, etas, luminosities,
        interpLogSpace);
  }

  // SWIG automatically converts a std::vector of NeutrinoSpecies to a
  // std::vector of ints
  static NeutrinoHistoryBlackBody CreateConstantSWIG(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<int>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& luminosities,
      const bool interpLogSpace) {
    return DoCreateConstant(times, radii,
        ConvertIntVecToNeutrinoSpeciesVec(species), T9s, etas, luminosities,
        interpLogSpace);
  }

  static NeutrinoHistoryBlackBody CreateTimeDependentSWIG(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<int>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, radii,
        ConvertIntVecToNeutrinoSpeciesVec(species), T9s, etas, luminosities,
        interpLogSpace);
  }

  std::shared_ptr<NeutrinoDistribution> operator()(const double time) const;

  std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const {
    return
        std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>(
            new NeutrinoHistoryBlackBody(*this));
  }

  std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const {
    return std::shared_ptr<NeutrinoHistory>(new NeutrinoHistoryBlackBody(*this));
  }

  void PrintInfo(NetworkOutput * const pOutput) const;

private:
  template<typename T>
  NeutrinoHistoryBlackBody(const bool isConstant, const std::vector<T>& species,
      const GeneralPiecewiseLinearFunction<double>& radiusVsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& T9VsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& etaVsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& lVsTime) :
  NeutrinoHistory(isConstant, species),
  mRadiusVsTime(radiusVsTime),
  mT9VsTime(T9VsTime),
  mEtaVsTime(etaVsTime),
  mLVsTime(lVsTime) {}

  template<typename T>
  static NeutrinoHistoryBlackBody DoCreateConstant(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<T>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& luminosities,
      const bool interpLogSpace);

  template<typename T>
  static NeutrinoHistoryBlackBody DoCreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<T>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace);

  double ConvertLuminosityAndRadiusToNormalization(const double luminosity,
      const double radius);

  GeneralPiecewiseLinearFunction<double> mRadiusVsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mT9VsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mEtaVsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mLVsTime;
};

template<typename T>
NeutrinoHistoryBlackBody NeutrinoHistoryBlackBody::DoCreateConstant(
    const std::vector<double>& times,
    const std::vector<double>& radii,
    const std::vector<T>& species,
    const std::vector<double>& T9s,
    const std::vector<double>& etas,
    const std::vector<double>& luminosities,
    const bool interpLogSpace) {
  assert(times.size() == radii.size());
  assert(species.size() == T9s.size());
  assert(species.size() == etas.size());
  assert(species.size() == luminosities.size());

  GeneralPiecewiseLinearFunction<double> radiusVsTime(times, radii,
      interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> T9VsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(T9s.data(), T9s.size())), interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> etaVsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(etas.data(), etas.size())), interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> lVsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(luminosities.data(), luminosities.size())),
          interpLogSpace);

  return NeutrinoHistoryBlackBody(true, species, radiusVsTime, T9VsTime,
      etaVsTime, lVsTime);
}

template<typename T>
NeutrinoHistoryBlackBody NeutrinoHistoryBlackBody::DoCreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<double>& radii,
      const std::vector<T>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& luminosities,
      const bool interpLogSpace) {
  assert(times.size() == radii.size());
  assert(times.size() == T9s.size());
  assert(times.size() == etas.size());
  assert(times.size() == luminosities.size());

  for (unsigned int i = 0; i < times.size(); ++i) {
    assert(species.size() == T9s[i].size());
    assert(species.size() == etas[i].size());
    assert(species.size() == luminosities[i].size());
  }

  GeneralPiecewiseLinearFunction<double> radiusVsTime(times, radii,
      interpLogSpace);

  std::vector<std::valarray<double>> T9VsTime(times.size(),
      std::valarray<double>(species.size()));
  std::vector<std::valarray<double>> etaVsTime(times.size(),
      std::valarray<double>(species.size()));
  std::vector<std::valarray<double>> lVsTime(times.size(),
      std::valarray<double>(species.size()));

  for (unsigned int i = 0; i < times.size(); ++i) {
    T9VsTime[i] = std::valarray<double>(T9s[i].data(), species.size());
    etaVsTime[i] = std::valarray<double>(etas[i].data(), species.size());
    lVsTime[i] = std::valarray<double>(luminosities[i].data(), species.size());
  }

  return NeutrinoHistoryBlackBody(false, species, radiusVsTime,
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, T9VsTime,
          interpLogSpace),
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, etaVsTime,
          interpLogSpace),
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, lVsTime,
          interpLogSpace));
}

#endif // SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYBLACKBODY_HPP_
