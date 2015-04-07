
/// \file FFNEntry.hpp
/// \author lroberts
/// \since May 13, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_NEUTRINOENTRY_HPP_
#define SKYNET_REACTIONS_NEUTRINOENTRY_HPP_

#include <string>
#include <vector>

#include "Reactions/Reaction.hpp"

class NeutrinoEntry {
 public:

  NeutrinoEntry(const std::vector<int> ParentsAs,
      const std::vector<int> ParentZs, const std::vector<int> DaughterAs,
      const std::vector<int> DaughterZs, const double Q,
      const double matrixElement, const double Wm=0.0, const bool isNue=true);

  Reaction GetReaction(const NuclideLibrary& nuclib) const;
  double GetQ() const { return mQ;}
  double GetMatrixElement() const { return mMatrixElement;}
  double GetWm() const { return mWm;}
  double IsNueReaction() const { return mIsNue; }

 private:
  std::vector<std::string> mParentNames, mDaughterNames;
  double mQ;
  double mMatrixElement;
  double mWm;
  bool mIsNue;
};

#endif // SKYNET_REACTIONS_NEUTRINOENTRY_HPP_
