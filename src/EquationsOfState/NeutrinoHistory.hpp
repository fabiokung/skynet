/// \file NeutrinoHistory.hpp
/// \author lroberts
/// \since September 10, 2015
///
/// \brief
///
///

#ifndef SKYNET_EQUATIONSOFSTATE_NEUTRINOHISTORY_HPP_
#define SKYNET_EQUATIONSOFSTATE_NEUTRINOHISTORY_HPP_

#include <memory>

#include "EquationsOfState/NeutrinoDistribution.hpp"
#include "Network/NetworkOutput.hpp"
#include "Utilities/FunctionVsTime.hpp"

class NeutrinoHistory:
    public FunctionVsTime<std::shared_ptr<NeutrinoDistribution>> {
public:
  virtual std::shared_ptr<NeutrinoDistribution> operator()(
      const double time) const =0;

  virtual std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const =0;

  virtual std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const =0;

  virtual void PrintInfo(NetworkOutput * const pOutput) const =0;

  static std::vector<NeutrinoSpecies> ConvertIntVecToNeutrinoSpeciesVec(
      const std::vector<int>& vec) {
    std::vector<NeutrinoSpecies> res(vec.size());

    for (unsigned int i = 0; i < res.size(); ++i)
      res[i] = static_cast<NeutrinoSpecies>(vec[i]);

    return res;
  }

protected:
  NeutrinoHistory(const bool isConst,
      const std::vector<NeutrinoSpecies>& species) :
      mIsConst(isConst),
      mSpecies(species.data(), species.size()) {}

  NeutrinoHistory(const bool isConst, const std::vector<int>& species) :
    NeutrinoHistory(isConst, ConvertIntVecToNeutrinoSpeciesVec(species)) {}

  bool mIsConst;
  std::valarray<NeutrinoSpecies> mSpecies;
};

class DummyNeutrinoHistory: public NeutrinoHistory {
public:
  DummyNeutrinoHistory() :
    NeutrinoHistory(true, std::vector<NeutrinoSpecies>(0)) {}

  std::shared_ptr<NeutrinoDistribution> operator()(
      const double /*time*/) const {
    return std::shared_ptr<NeutrinoDistribution>(
        new DummyNeutrinoDistribution());
  }

  std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>
  MakeUniquePtr() const {
    return
        std::unique_ptr<FunctionVsTime<std::shared_ptr<NeutrinoDistribution>>>(
            new DummyNeutrinoHistory(*this));
  }

  std::shared_ptr<NeutrinoHistory> MakeSharedPtr() const {
    return std::shared_ptr<NeutrinoHistory>(new DummyNeutrinoHistory(*this));
  }

  void PrintInfo(NetworkOutput * const pOutput) const {
    pOutput->Log("# Neutrino History: Dummy (i.e. no neutrinos)\n");
  }
};


#endif // SKYNET_EQUATIONSOFSTATE_NEUTRINOHISTORY_HPP_
