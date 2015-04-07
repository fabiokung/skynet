/// \file NeutrinoHistoryFermiDirac.hpp
/// \author jlippuner
/// \since Feb 4, 2016
///
/// \brief
///
///

#ifndef SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYFERMIDIRAC_HPP_
#define SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYFERMIDIRAC_HPP_

#include <cassert>

#include "EquationsOfState/NeutrinoHistory.hpp"
#include "EquationsOfState/NeutrinoDistributionFermiDirac.hpp"

#include "Utilities/Interpolators/PiecewiseLinearFunction.hpp"

class NeutrinoHistoryFermiDirac : public NeutrinoHistory {
public:
  static NeutrinoHistoryFermiDirac CreateConstant(
      const std::vector<double>& times,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& normalizations,
      const bool interpLogSpace) {
    return DoCreateConstant(times, species, T9s, etas, normalizations,
        interpLogSpace);
  }

  static NeutrinoHistoryFermiDirac CreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<NeutrinoSpecies>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& normalizations,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, species, T9s, etas, normalizations,
        interpLogSpace);
  }

  static NeutrinoHistoryFermiDirac CreateConstantSWIG(
      const std::vector<double>& times,
      const std::vector<int>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& normalizations,
      const bool interpLogSpace) {
    return DoCreateConstant(times, species, T9s, etas, normalizations,
        interpLogSpace);
  }

  static NeutrinoHistoryFermiDirac CreateTimeDependentSWIG(
      const std::vector<double>& times,
      const std::vector<int>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& normalizations,
      const bool interpLogSpace) {
    return DoCreateTimeDependent(times, species, T9s, etas, normalizations,
        interpLogSpace);
  }

  std::shared_ptr<NeutrinoDistribution> operator()(const double time) const {
    return NeutrinoDistributionFermiDirac::Create(mSpecies, mT9VsTime(time),
        mEtaVsTime(time), mNormVsTime(time));
  }

  std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const {
    return
        std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>(
            new NeutrinoHistoryFermiDirac(*this));
  }

  std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const {
    return std::shared_ptr<NeutrinoHistory>(
        new NeutrinoHistoryFermiDirac(*this));
  }

  void PrintInfo(NetworkOutput * const pOutput) const;

private:
  template<typename T>
  NeutrinoHistoryFermiDirac(const bool isConstant,
      const std::vector<T>& species,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& T9VsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& etaVsTime,
      const GeneralPiecewiseLinearFunction<std::valarray<double>>& normVsTime) :
  NeutrinoHistory(isConstant, species),
  mT9VsTime(T9VsTime),
  mEtaVsTime(etaVsTime),
  mNormVsTime(normVsTime) {}

  template<typename T>
  static NeutrinoHistoryFermiDirac DoCreateConstant(
      const std::vector<double>& times,
      const std::vector<T>& species,
      const std::vector<double>& T9s,
      const std::vector<double>& etas,
      const std::vector<double>& normalizations,
      const bool interpLogSpace);

  template<typename T>
  static NeutrinoHistoryFermiDirac DoCreateTimeDependent(
      const std::vector<double>& times,
      const std::vector<T>& species,
      const std::vector<std::vector<double>>& T9s,
      const std::vector<std::vector<double>>& etas,
      const std::vector<std::vector<double>>& normalizations,
      const bool interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> mT9VsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mEtaVsTime;
  GeneralPiecewiseLinearFunction<std::valarray<double>> mNormVsTime;
};

template<typename T>
NeutrinoHistoryFermiDirac NeutrinoHistoryFermiDirac::DoCreateConstant(
    const std::vector<double>& times,
    const std::vector<T>& species,
    const std::vector<double>& T9s,
    const std::vector<double>& etas,
    const std::vector<double>& normalizations,
    const bool interpLogSpace) {
  assert(species.size() == T9s.size());
  assert(species.size() == etas.size());
  assert(species.size() == normalizations.size());

  GeneralPiecewiseLinearFunction<std::valarray<double>> T9VsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(T9s.data(), T9s.size())), interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> etaVsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(etas.data(), etas.size())), interpLogSpace);

  GeneralPiecewiseLinearFunction<std::valarray<double>> normVsTime(times,
      std::vector<std::valarray<double>>(times.size(),
          std::valarray<double>(normalizations.data(), normalizations.size())),
      interpLogSpace);

  return NeutrinoHistoryFermiDirac(true, species, T9VsTime, etaVsTime,
      normVsTime);
}

template<typename T>
NeutrinoHistoryFermiDirac NeutrinoHistoryFermiDirac::DoCreateTimeDependent(
    const std::vector<double>& times,
    const std::vector<T>& species,
    const std::vector<std::vector<double>>& T9s,
    const std::vector<std::vector<double>>& etas,
    const std::vector<std::vector<double>>& normalizations,
    const bool interpLogSpace) {
  assert(times.size() == T9s.size());
  assert(times.size() == etas.size());
  assert(times.size() == normalizations.size());

  for (unsigned int i = 0; i < times.size(); ++i) {
    assert(species.size() == T9s[i].size());
    assert(species.size() == etas[i].size());
    assert(species.size() == normalizations[i].size());
  }

  std::vector<std::valarray<double>> T9VsTime(times.size(),
      std::valarray<double>(species.size()));
  std::vector<std::valarray<double>> etaVsTime(times.size(),
      std::valarray<double>(species.size()));
  std::vector<std::valarray<double>> normVsTime(times.size(),
      std::valarray<double>(species.size()));

  for (unsigned int i = 0; i < times.size(); ++i) {
    T9VsTime[i] = std::valarray<double>(T9s[i].data(), species.size());
    etaVsTime[i] = std::valarray<double>(etas[i].data(), species.size());
    normVsTime[i] = std::valarray<double>(normalizations[i].data(),
        species.size());
  }

  return NeutrinoHistoryFermiDirac(false, species,
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, T9VsTime,
          interpLogSpace),
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, etaVsTime,
          interpLogSpace),
      GeneralPiecewiseLinearFunction<std::valarray<double>>(times, normVsTime,
          interpLogSpace));
}

#endif // SRC_EQUATIONSOFSTATE_NEUTRINOHISTORYFERMIDIRAC_HPP_
