/// \file NetworkOutputDataset.hpp
/// \author jlippuner
/// \since Aug 17, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_NETWORKOUTPUTDATASET_HPP_
#define SKYNET_NETWORK_NETWORKOUTPUTDATASET_HPP_

#include <stdexcept>
#include <string>
#include <vector>

#include <H5Cpp.h>

template<typename T>
class NetworkOutputDataset {
public:
  NetworkOutputDataset();

  void InitForWriting(H5::H5File * const pH5File, const hsize_t ChunkSize1D,
      const hsize_t ChunkSize2D, const bool reduceCache);

  void WriteToFile();

  void ReadFromFile(const H5::H5File& H5File);

  void Close();

  void SetName(const std::string& name) {
    mName = name;
  }

  void SetDescription(const std::string& description) {
    mDescription = description;
  }

  void SetUnits(const std::string& units) {
    mUnits = units;
  }

  void SetDim2Size(const std::size_t dim2Size);

  const std::vector<T>& Data() const {
    return mData;
  }

  std::vector<T>& Data() {
    return mData;
  }

  const std::string& Name() const {
    return mName;
  }

  bool Is2D() const {
    return (mRank == 2);
  }

  bool ReduceCache() const {
    return mReduceCache;
  }

private:
  void AddStringAttribute(const std::string& name, const std::string& value);

  std::string GetStringAttribute(const std::string& name);

  const H5::DataType mDataType;
  const std::size_t mRank;
  std::vector<hsize_t> mDim;

  bool mReduceCache;
  std::size_t mNumRemovedCacheEntries;

  std::string mName;
  std::string mDescription;
  std::string mUnits;

  H5::DataSet mH5Dataset;
  std::vector<T> mData;
};

#endif // SKYNET_NETWORK_NETWORKOUTPUTDATASET_HPP_
