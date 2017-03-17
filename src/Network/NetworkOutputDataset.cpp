/// \file NetworkOutputDataset.cpp
/// \author jlippuner
/// \since Aug 17, 2014
///
/// \brief
///
///

#include "Network/NetworkOutputDataset.hpp"

#include <cstring>

#include "Network/NetworkOutput.hpp"
#include "Utilities/H5Helper.hpp"

namespace { // unnamed so that these functions can only be used in this file

template<typename T>
struct GetDatasetType {
  typedef T type;
};

template<typename T>
struct GetDatasetType<std::vector<T>> {
  typedef T type;
};

template<typename T>
  std::size_t GetDatasetRank(const NetworkOutputDataset<T> * const) {
  return 1;
}

template<typename T>
  std::size_t GetDatasetRank(
      const NetworkOutputDataset<std::vector<T>> * const) {
  return 2;
}


template<typename T>
const std::vector<T>& GetDatasetInnermostDataVector(
    const NetworkOutputDataset<T> * const pDataset,
    const unsigned int /*idx*/) {
  return pDataset->Data();
}

template<typename T>
const std::vector<T>& GetDatasetInnermostDataVector(
    const NetworkOutputDataset<std::vector<T>> * const pDataset,
    const unsigned int idx) {
  return pDataset->Data()[idx];
}

template<typename T>
std::vector<T>& GetDatasetInnermostDataVector(
    NetworkOutputDataset<T> * const pDataset, const unsigned int) {
  return pDataset->Data();
}

template<typename T>
std::vector<T>& GetDatasetInnermostDataVector(
    NetworkOutputDataset<std::vector<T>> * const pDataset, const unsigned int idx) {
  return pDataset->Data()[idx];
}

} // namespace [unnamed]

template<typename T>
NetworkOutputDataset<T>::NetworkOutputDataset() :
    mDataType(GetH5DataType<typename GetDatasetType<T>::type>()),
    mRank(GetDatasetRank(this)),
    mDim(mRank),
    mReduceCache(false),
    mNumRemovedCacheEntries(0) {}

template<typename T>
void NetworkOutputDataset<T>::InitForWriting(H5::H5File * const pH5File,
    const hsize_t ChunkSize1D, const hsize_t ChunkSize2D,
    const bool reduceCache) {
  mDim[0] = 0;
  mReduceCache = reduceCache;

  hsize_t maxDim[mRank];
  maxDim[0] = H5S_UNLIMITED;
  if (Is2D()) maxDim[1] = mDim[1];

  H5::DataSpace dataSpace(mRank, mDim.data(), maxDim);

  H5::DSetCreatPropList props;
  props.setDeflate(5);

  hsize_t chunkDim[mRank];
  if (Is2D()) {
    chunkDim[0] = ChunkSize2D;
    chunkDim[1] = mDim[1];
  } else {
    chunkDim[0] = ChunkSize1D;
  }
  props.setChunk(mRank, chunkDim);

  mH5Dataset = pH5File->createDataSet(mName.c_str(), mDataType, dataSpace, props);

  AddStringAttribute("Description", mDescription);
  AddStringAttribute("Units", mUnits);
}

template<typename T>
void NetworkOutputDataset<T>::WriteToFile() {
  std::size_t totalNumRecs = mNumRemovedCacheEntries + mData.size();
  if (totalNumRecs < mDim[0])
    throw std::runtime_error("Total number of records is smaller than "
        "what has already been written to file");

  if (totalNumRecs == mDim[0]) {
    // there is no new data to write
    return;
  }

  std::size_t numEntriesToWrite = totalNumRecs - mDim[0];

  // dimension of the space by which we are extending
  hsize_t extendDim[mRank];
  extendDim[0] = numEntriesToWrite;
  if (Is2D()) extendDim[1] = mDim[1];

  // make new size
  hsize_t newSize[mRank];
  newSize[0] = mDim[0] + extendDim[0];
  if (Is2D()) newSize[1] = mDim[1];

  mH5Dataset.extend(newSize);

  // select hyperslap to write new data to
  hsize_t offset[mRank];
  offset[0] = mDim[0];
  if (Is2D()) offset[1] = 0;
  H5::DataSpace fileSpace = mH5Dataset.getSpace();
  fileSpace.selectHyperslab(H5S_SELECT_SET, extendDim, offset);

  H5::DataSpace memSpace(mRank, extendDim);

  if (Is2D()) {
    // we need to arrange the data contiguously
    // we must use new here, otherwise the memory will be allocated on the stack,
    // which is rather small and would lead to segfaults for large datasets
    std::size_t typeSize = sizeof(typename GetDatasetType<T>::type);
    char * const rawData = new char[extendDim[0] * extendDim[1] * typeSize];

    for (unsigned int i = 0; i < extendDim[0]; ++i) {
      auto dataPtr = GetDatasetInnermostDataVector(this,
          i + mDim[0] - mNumRemovedCacheEntries).data();
      memcpy(rawData + (i * extendDim[1] * typeSize), dataPtr,
          extendDim[1] * typeSize);
    }

    mH5Dataset.write(rawData, mDataType, memSpace, fileSpace);
    delete[] rawData;
  } else {
    mH5Dataset.write(mData.data() + (mDim[0] - mNumRemovedCacheEntries),
        mDataType, memSpace, fileSpace);
  }

  mDim[0] += numEntriesToWrite;
  mH5Dataset.flush(H5F_SCOPE_GLOBAL);

  if (mReduceCache && (mData.size() > NetworkOutput::CacheSize)) {
    std::size_t numToRemove = mData.size() - NetworkOutput::CacheSize;
    mData.erase(mData.begin(), mData.begin() + numToRemove);
    mNumRemovedCacheEntries += numToRemove;
  }
}

template<typename T>
void NetworkOutputDataset<T>::ReadFromFile(const H5::H5File& H5File) {
  mH5Dataset = H5File.openDataSet(mName.c_str());
  H5::DataSpace dataSpace = mH5Dataset.getSpace();

  if (mRank != (std::size_t)dataSpace.getSimpleExtentDims(mDim.data()))
    throw std::runtime_error("Unexpected number of dimensions in data set "
        + mName);

  if (!(mH5Dataset.getDataType() == mDataType))
    throw std::runtime_error("Wrong data type in data set " + mName);

  std::string desc = GetStringAttribute("Description");
  std::string units = GetStringAttribute("Units");

  if (desc != mDescription)
    throw std::runtime_error("Wrong description in data set " + mName);

  if (units != mUnits)
    throw std::runtime_error("Wrong units in data set " + mName);

  // set size of data vector
  mData.resize(mDim[0]);

  if (Is2D()) {
    // we read one row at a time, because each row goes into a separate
    // std::vector
    hsize_t readDim[mRank];
    readDim[0] = 1;
    readDim[1] = mDim[1];
    H5::DataSpace memSpace(mRank, readDim);

    hsize_t offset[mRank];
    offset[1] = 0;

    for (unsigned int i = 0; i < mDim[0]; ++i) {
      auto& dataVector = GetDatasetInnermostDataVector(this, i);
      dataVector.resize(mDim[1]);

      // select the row in the file
      offset[0] = i;
      H5::DataSpace fileSpace = dataSpace;
      fileSpace.selectHyperslab(H5S_SELECT_SET, readDim, offset);

      mH5Dataset.read(dataVector.data(), mDataType, memSpace, fileSpace);
    }
  } else {
    mH5Dataset.read(mData.data(), mDataType);
  }
}

template<typename T>
void NetworkOutputDataset<T>::Close() {
  mH5Dataset.close();
}

template<typename T>
void NetworkOutputDataset<T>::SetDim2Size(const std::size_t dim2Size) {
  if (!Is2D())
    throw std::invalid_argument("Tried to set dim 2 size of the 1-dim "
        "NetworkOutputDataset " + mName);
  else
    mDim[1] = dim2Size;
}

template<typename T>
void NetworkOutputDataset<T>::AddStringAttribute(const std::string& name,
    const std::string& value) {
  hsize_t dim[1] = { 1 };
  H5::DataSpace attrDataspace(1, dim);

  H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);

  H5::Attribute attr = mH5Dataset.createAttribute(name.c_str(), strType, attrDataspace);

  const char * valueCharPtr = value.c_str();
  attr.write(strType, &valueCharPtr);

  attr.close();
}

template<typename T>
std::string NetworkOutputDataset<T>::GetStringAttribute(
    const std::string& name) {
  H5::Attribute attr = mH5Dataset.openAttribute(name.c_str());

  hsize_t dim[1];
  H5::DataSpace attrSpace = attr.getSpace();
  if (attrSpace.getSimpleExtentNdims() != 1)
    throw std::runtime_error("String attribute space does not have rank 1");

  attrSpace.getSimpleExtentDims(dim);
  if (dim[0] != 1)
    throw std::runtime_error("String attribute space does not have size 1");

  H5::DataType dataType = attr.getDataType();
  H5::StrType strType(H5::PredType::C_S1, H5T_VARIABLE);
  if (!(dataType == strType))
    throw std::runtime_error("String attribute does not have variable string "
        "data type");

  // the space is allocated by HDF5, but we need to free it
  char *charPtr;
  attr.read(dataType, &charPtr);

  attr.close();

  std::string value(charPtr);
  free(charPtr);

  return value;
}

// explicit template instantiation
template class NetworkOutputDataset<double>;
template class NetworkOutputDataset<int>;
template class NetworkOutputDataset<std::vector<double>>;
