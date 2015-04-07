/// \file ReactionData.hpp
/// \author jlippuner
/// \since Jan 14, 2015
///
/// \brief
///
///

#ifndef SKYNET_REACTIONS_REACTIONDATA_HPP_
#define SKYNET_REACTIONS_REACTIONDATA_HPP_

#include <stdexcept>
#include <vector>

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

// so that we can loop over a list of ReactionData classes without needing to
// know their types
class GeneralReactionData {
public:
  virtual ~GeneralReactionData() {}

  virtual void RemoveData(const std::vector<bool>& removeFlags) =0;

  virtual void SetActive(const std::vector<bool>& activeFlags) =0;
};

template<typename T>
class ReactionData : public GeneralReactionData {
public:
  ReactionData() :
      mAllData(),
      mActiveData(),
      mActiveFlags() {}

  explicit ReactionData(const std::vector<T>& data) :
      mAllData(data),
      mActiveData(data),
      mActiveFlags(data.size(), true) {}

  std::size_t size() const {
    return mActiveData.size();
  }

  inline T operator[](const std::size_t idx) const {
    return mActiveData[idx];
  }

  ReactionData& operator=(const std::vector<T>& data) {
    mAllData = data;
    mActiveData = data;
    mActiveFlags = std::vector<bool>(mActiveData.size(), true);

    return *this;
  }

  const std::vector<bool>& ActiveFlags() const {
    return mActiveFlags;
  }

  void RemoveData(const std::vector<bool>& removeFlags) {
    if (removeFlags.size() != mAllData.size())
      throw std::invalid_argument("removeFlag must have same length as "
          "mAllData");

    std::vector<T> newAllData;
    std::vector<bool> newActiveFlags;

    for (unsigned int i = 0; i < mAllData.size(); i++) {
      if (!removeFlags[i]) {
        newAllData.push_back(mAllData[i]);
        newActiveFlags.push_back(mActiveFlags[i]);
      }
    }

    // set the new data
    mAllData = newAllData;
    SetActive(newActiveFlags);
  }

  void SetActive(const std::vector<bool>& activeFlags) {
    if (activeFlags.size() != mAllData.size())
      throw std::invalid_argument("activeFlags must have same length as "
          "mAllData");

    mActiveFlags = activeFlags;
    mActiveData = std::vector<T>();

    for (unsigned int i = 0; i < activeFlags.size(); ++i) {
      if (activeFlags[i])
        mActiveData.push_back(mAllData[i]);
    }
  }

  const std::vector<T>& AllData() const {
    return mAllData;
  }

  const std::vector<T>& ActiveData() const {
    return mActiveData;
  }

private:
  friend class boost::serialization::access;
  template<class AR>
  void serialize(AR& ar, const unsigned int /*version*/) {
    ar & mAllData;
    ar & mActiveData;
    ar & mActiveFlags;
  }

  std::vector<T> mAllData; // data for all reactions
  std::vector<T> mActiveData; // data for only active reactions
  std::vector<bool> mActiveFlags; // which reactions are active
};

#endif // SKYNET_REACTIONS_REACTIONDATA_HPP_
