/// \file NetworkOutput.hpp
/// \author jlippuner
/// \since Aug 1, 2014
///
/// \brief
///
///

#ifndef SKYNET_NETWORK_NETWORKOUTPUT_HPP_
#define SKYNET_NETWORK_NETWORKOUTPUT_HPP_

#include <chrono>
#include <cstdio>
#include <ctime>
#include <utility>
#include <string>
#include <vector>

#include <H5Cpp.h>

#include "EquationsOfState/EOS.hpp"
#include "Network/NetworkOutputDataset.hpp"
#include "NuclearData/NuclideLibrary.hpp"
#include "Utilities/FlushStdout.hpp"

class NetworkOptions;

class NetworkOutput {
public:
  static NetworkOutput CreateNew(const std::string& filePrefix,
      const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput);

  static NetworkOutput CreateNew(const std::string& filePrefix,
      const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput,
      const std::vector<int>& massFractionColumnNuclideIds);

  static NetworkOutput ReadFromFile(const std::string& h5FilePath);

  static void MakeDatFile(const std::string& h5FilePath);

  static std::string GetDateTimeStr();

  NetworkOutput(const NetworkOutput& other);

  void Close();

  ~NetworkOutput();

  template<typename... Args>
  void Log(const char * format, Args... args);

  template<typename... Args>
  void LogWithTimestamp(const char * format, Args... args);

  void WriteColumnDescriptionToLog();

  void AddEntry(const double time,
      const std::pair<double, std::string> dtAndLimitReason, const int numFailed,
      const int numNR, const ThermodynamicState thermoState,
      const double heatingRate, const std::vector<double>& Y,
      const double totalMass);

  const std::vector<double>& DensityVsTime() const {
    return mDensityVsTime.Data();
  }

  const std::vector<double>& DtVsTime() const {
    return mDtVsTime.Data();
  }

  const std::vector<double>& EntropyVsTime() const {
    return mEntropyVsTime.Data();
  }

  const std::vector<double>& HeatingRateVsTime() const {
    return mHeatingRateVsTime.Data();
  }

  const std::vector<double>& TemperatureVsTime() const {
    return mTemperatureVsTime.Data();
  }

  const std::vector<double>& Times() const {
    return mTimes.Data();
  }

  const std::vector<double>& YeVsTime() const {
    return mYeVsTime.Data();
  }

  const std::vector<double>& EtaEVsTime() const {
    return mEtaEVsTime.Data();
  }

  const std::vector<std::vector<double>>& YVsTime() const {
    return mYVsTime.Data();
  }

  const std::vector<std::vector<double>>& YVsAVsTime() const {
    return mYVsAVsTime;
  }

  const std::vector<int>& A() const {
    return mA.Data();
  }

  const std::vector<int>& Z() const {
    return mZ.Data();
  }

  const std::vector<double>& BindingEnergy() const {
    return mBindingEnergy.Data();
  }

  const std::vector<double>& FinalY() const {
    return mYVsTime.Data()[mYVsTime.Data().size() - 1];
  }

  const std::vector<double>& FinalYVsA() const {
    return mYVsAVsTime[mYVsAVsTime.size() - 1];
  }

  bool IsReadOnly() const {
    return mIsReadOnly;
  }

  std::size_t NumEntries() const {
    return mNumEntries;
  }

  std::size_t NumNuclides() const {
    return mpNuclib->NumNuclides();
  }

  const NuclideLibrary& GetNuclideLibrary() const {
    return *mpNuclib;
  }

  // All the data will be written to disk when there are BufferSize unwritten
  // entries
  static const std::size_t BufferSize = 8;

  // For NetworkOutputDatasets that have the ReduceCache option enabled,
  // only this many records will be kept in memory
  static const std::size_t CacheSize = 32;

  static_assert(CacheSize >= BufferSize, "CacheSize must be at least as "
      "big as BufferSize");

private:
  NetworkOutput(const std::string& filePrefix,
      const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput,
      const std::vector<int>& massFractionColumnNuclideIds);

  NetworkOutput(const std::string& filePrefix);

  void OpenLogFile(const bool append);

  void SetMetadata();

  template<typename T>
  void InitH5Dataset(NetworkOutputDataset<T> * const pDataset,
      const bool reduceCache);

  void WriteToDisk();

  std::vector<double> MakeYVsA(const std::vector<double>& Y) const;

  std::vector<double> MakeYVsZ(const std::vector<double>& Y) const;

  // 1D data will be written in chunks of size ChunkSize1D
  static const hsize_t ChunkSize1D = 256;

  // 2D data will be written in chunks of size ChunkSize2D x dataDim[1]
  static const hsize_t ChunkSize2D = 8;

  // the column header mColumnHeader is repeated every
  // ColumnHeaderPrintFrequency entries
  static const int ColumnHeaderPrintFrequency = 50;

  // A NetworkOutput read from an H5 file is read-only, an empty newly created
  // H5 file is not read-only.
  bool mIsReadOnly;

  // the HDF5 output file will be mFilePrevix + ".h5" and the log file will be
  // mFilePrefix + ".log"
  std::string mFilePrefix;

  // HDF5 file for binary output
  H5::H5File mFile;

  // log file for text output
  FILE * mLogFile;

  // true if output to stdout is disabled
  bool mDisableStdoutOutput;

  // names and ids of the nuclides whose mass fractions should be written to the
  // log file
  std::vector<std::pair<std::string, int>> mMassFractionColumnNamesAndIds;

  // column header that is repeated every ColumnHeaderPrintFrequency entries
  std::string mColumnHeader;

  // copy of the nuclide library
  std::unique_ptr<NuclideLibrary> mpNuclib;

  // maximum value of A
  int mMaxA;
  int mMaxZ;

  // current total number of entries
  std::size_t mNumEntries;

  // current total number of entries that have been written to the H5 file
  std::size_t mNumEntriesWrittenToFile;

  NetworkOutputDataset<double> mTimes; // in s
  NetworkOutputDataset<double> mDtVsTime; // in s
  NetworkOutputDataset<double> mTemperatureVsTime; // in GK
  NetworkOutputDataset<double> mDensityVsTime; // in g / cm^3
  NetworkOutputDataset<double> mEntropyVsTime; // in k_b / baryon
  NetworkOutputDataset<double> mHeatingRateVsTime; // in erg / s / g
  NetworkOutputDataset<double> mYeVsTime; // dimensionless
  NetworkOutputDataset<double> mEtaEVsTime; // dimensionless
  NetworkOutputDataset<std::vector<double>> mYVsTime; // dimensionless

  NetworkOutputDataset<int> mA; // nuclide mass number
  NetworkOutputDataset<int> mZ; // nuclide atomic number
  NetworkOutputDataset<double> mBindingEnergy; // nuclide binding energy in MeV

  // this data is provided for convenience and not stored in the H5 file
  std::vector<std::vector<double>> mYVsAVsTime;
};

// these template member functions must be in the header file

template<typename... Args>
void NetworkOutput::Log(const char * format, Args... args) {
  if (mIsReadOnly)
    throw std::runtime_error("Tried to write to a read-only NetworkOutput.");

  fprintf(mLogFile, format, args...);
  fflush(mLogFile);

  if (!mDisableStdoutOutput) {
    printf(format, args...);
    // flush stdout, this seems to be necessary for big networks
    FlushStdout::Flush();
  }
}

template<typename... Args>
void NetworkOutput::LogWithTimestamp(const char * format, Args... args) {
  std::string formatStr = "# [" + GetDateTimeStr() + "]\n"
      + std::string(format);
  Log(formatStr.c_str(), args...);
}

#endif // SKYNET_NETWORK_NETWORKOUTPUT_HPP_
