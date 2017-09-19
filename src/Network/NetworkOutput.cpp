/// \file NetworkOutput.cpp
/// \author jlippuner
/// \since Aug 1, 2014
///
/// \brief
///
///

#include "Network/NetworkOutput.hpp"

#include <algorithm>

#include "BuildInfo.hpp"
#include "Utilities/File.hpp"

NetworkOutput NetworkOutput::CreateNew(const std::string& filePrefix,
    const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput) {
  std::vector<int> massFractionColumnNuclideIds;
  for (int i = 0; i < std::min(nuclideLibrary.NumNuclides(), 2); ++i)
    massFractionColumnNuclideIds.push_back(i);

  return NetworkOutput(filePrefix, nuclideLibrary, disableStdoutOutput,
      massFractionColumnNuclideIds);
}

NetworkOutput NetworkOutput::CreateNew(const std::string& filePrefix,
    const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput,
    const std::vector<int>& massFractionColumnNuclideIds) {
  return NetworkOutput(filePrefix, nuclideLibrary, disableStdoutOutput,
      massFractionColumnNuclideIds);
}

NetworkOutput NetworkOutput::ReadFromFile(const std::string& h5FilePath) {
  return NetworkOutput(h5FilePath);
}

std::string NetworkOutput::GetDateTimeStr() {
  auto now = std::chrono::system_clock::now();
  auto nowTime_t = std::chrono::system_clock::to_time_t(now);
  auto local = localtime(&nowTime_t);

  int hour = local->tm_hour;
  std::string AMPM = (hour < 12 ? "AM" : "PM");
  hour = hour % 12;
  if (hour == 0)
    hour = 12;

  char buf[16];
  sprintf(buf, "%04i-%02i-%02i", local->tm_year + 1900, local->tm_mon + 1,
      local->tm_mday);
  std::string dateStr(buf);

  sprintf(buf, "%02i:%02i:%02i %s", hour, local->tm_min, local->tm_sec,
      AMPM.c_str());
  std::string timeStr(buf);

  return dateStr + " @ " + timeStr;
}

void NetworkOutput::MakeDatFile(const std::string& h5FilePath) {
  std::string outFile = h5FilePath + ".dat";
  FILE * f = fopen(outFile.c_str(), "w");

  if (f == nullptr) {
    printf("ERROR: could not open '%s' for writing\n", outFile.c_str());
    return;
  }

  int nIdx = -1;
  int pIdx = -1;
  int aIdx = -1;

  NetworkOutput out(h5FilePath);

  for (unsigned int i = 0; i < out.A().size(); ++i) {
    if (out.A()[i] == 1) {
      if (out.Z()[i] == 0) {
        nIdx = i;
        continue;
      } else if (out.Z()[i] == 1) {
        pIdx = i;
        continue;
      }
    } else if ((out.A()[i] == 4) && (out.Z()[i] == 2)) {
      aIdx = i;
      continue;
    }

    if ((nIdx >= 0) && (pIdx >= 0) && (aIdx >= 0))
      break;
  }

  if (nIdx == -1) {
    printf("ERROR: could not find neutron in '%s'\n", h5FilePath.c_str());
    return;
  }
  if (pIdx == -1) {
    printf("ERROR: could not find proton in '%s'\n", h5FilePath.c_str());
    return;
  }
  if (aIdx == -1) {
    printf("ERROR: could not find alpha in '%s'\n", h5FilePath.c_str());
    return;
  }

  std::vector<double> Yn(out.Times().size(), 0.0);
  std::vector<double> Yp(out.Times().size(), 0.0);
  std::vector<double> Ya(out.Times().size(), 0.0);
  std::vector<double> X2to11(out.Times().size(), 0.0);
  std::vector<double> Y2to11(out.Times().size(), 0.0);
  std::vector<double> X12to119(out.Times().size(), 0.0);
  std::vector<double> Y12to119(out.Times().size(), 0.0);
  std::vector<double> X120to249(out.Times().size(), 0.0);
  std::vector<double> Y120to249(out.Times().size(), 0.0);
  std::vector<double> X250andUp(out.Times().size(), 0.0);
  std::vector<double> Y250andUp(out.Times().size(), 0.0);
  std::vector<double> Xlan(out.Times().size(), 0.0);
  std::vector<double> Ylan(out.Times().size(), 0.0);
  std::vector<double> Xact(out.Times().size(), 0.0);
  std::vector<double> Yact(out.Times().size(), 0.0);
  std::vector<double> Abar(out.Times().size(), 0.0);
  std::vector<double> Abar12andUp(out.Times().size(), 0.0);

  for (unsigned int i = 0; i < out.Times().size(); ++i) {
    Yn[i] = out.YVsTime()[i][nIdx];
    Yp[i] = out.YVsTime()[i][pIdx];
    Ya[i] = out.YVsTime()[i][aIdx];

    double sumY = 0.0;
    double sumX = 0.0;

    for (unsigned int j = 0; j < out.A().size(); ++j) {
      int a = out.A()[j];
      int z = out.Z()[j];
      double y = out.YVsTime()[i][j];
      double x = y * (double)a;

      if ((a >= 2) && (a <= 11)) {
        X2to11[i] += x;
        Y2to11[i] += y;
      } else if ((a >= 12) && (a <= 119)) {
        X12to119[i] += x;
        Y12to119[i] += y;
      } else if ((a >= 120) && (a <= 249)) {
        X120to249[i] += x;
        Y120to249[i] += y;
      } else if (a >= 250) {
        X250andUp[i] += x;
        Y250andUp[i] += y;
      }

      if ((z >= 58) && (z <= 71)) {
        Xlan[i] += x;
        Ylan[i] += y;
      } else if ((z >= 90) && (z <= 103)) {
        Xact[i] += x;
        Yact[i] += y;
      }

      sumX += x;
      sumY += y;
    }

    Abar[i] = sumX / sumY;
    Abar12andUp[i] = (X12to119[i] + X120to249[i] + X250andUp[i]) /
        (Y12to119[i] + Y120to249[i] + Y250andUp[i]);
  }

  fprintf(f, "# [ 1] = time [s]\n");
  fprintf(f, "# [ 2] = dt [s]\n");
  fprintf(f, "# [ 3] = temperature [GK]\n");
  fprintf(f, "# [ 4] = density [g / cm^3]\n");
  fprintf(f, "# [ 5] = entropy [kb / baryon]\n");
  fprintf(f, "# [ 6] = specific heating rate [erg / s / g]\n");
  fprintf(f, "# [ 7] = Ye []\n");
  fprintf(f, "# [ 8] = neutron abundance []\n");
  fprintf(f, "# [ 9] = proton abundance []\n");
  fprintf(f, "# [10] = alpha abundance []\n");
  fprintf(f, "# [11] = Average A []\n");
  fprintf(f, "# [12] = Average A for A >= 12 []\n");
  fprintf(f, "# [13] = total mass fraction of 2 <= A <= 11 []\n");
  fprintf(f, "# [14] = total mass fraction of 12 <= A <= 119 []\n");
  fprintf(f, "# [15] = total mass fraction of 120 <= A <= 249 []\n");
  fprintf(f, "# [16] = total mass fraction of 250 <= A []\n");
  fprintf(f, "# [17] = total abundance of 2 <= A <= 11 []\n");
  fprintf(f, "# [18] = total abundance of 12 <= A <= 119 []\n");
  fprintf(f, "# [19] = total abundance of 120 <= A <= 249 []\n");
  fprintf(f, "# [20] = total abundance of 250 <= A []\n");
  fprintf(f, "# [21] = total mass fraction of Lanthanides (58 <= Z <= 71) []\n");
  fprintf(f, "# [22] = total mass fraction of Actinides (90 <= Z <= 103) []\n");
  fprintf(f, "# [23] = total abundance of Lanthanides (58 <= Z <= 71) []\n");
  fprintf(f, "# [24] = total abundance of Actinides (90 <= Z <= 103) []\n");

  for (unsigned int i = 0; i < out.Times().size(); ++i) {
    fprintf(f, "%25.15e", out.Times()[i]);
    fprintf(f, " %25.15e", out.DtVsTime()[i]);
    fprintf(f, " %25.15e", out.TemperatureVsTime()[i]);
    fprintf(f, " %25.15e", out.DensityVsTime()[i]);
    fprintf(f, " %25.15e", out.EntropyVsTime()[i]);
    fprintf(f, " %25.15e", out.HeatingRateVsTime()[i]);
    fprintf(f, " %25.15e", out.YeVsTime()[i]);
    fprintf(f, " %25.15e", Yn[i]);
    fprintf(f, " %25.15e", Yp[i]);
    fprintf(f, " %25.15e", Ya[i]);
    fprintf(f, " %25.15e", Abar[i]);
    fprintf(f, " %25.15e", Abar12andUp[i]);
    fprintf(f, " %25.15e", X2to11[i]);
    fprintf(f, " %25.15e", X12to119[i]);
    fprintf(f, " %25.15e", X120to249[i]);
    fprintf(f, " %25.15e", X250andUp[i]);
    fprintf(f, " %25.15e", Y2to11[i]);
    fprintf(f, " %25.15e", Y12to119[i]);
    fprintf(f, " %25.15e", Y120to249[i]);
    fprintf(f, " %25.15e", Y250andUp[i]);
    fprintf(f, " %25.15e", Xlan[i]);
    fprintf(f, " %25.15e", Xact[i]);
    fprintf(f, " %25.15e", Ylan[i]);
    fprintf(f, " %25.15e", Yact[i]);
    fprintf(f, "\n");
  }

  fprintf(f, "\n\n\n");
  fprintf(f, "# [1] = mass number A\n");
  fprintf(f, "# [2] = final abundance (of A)\n");
  fprintf(f, "# [3] = initial abundance (of A)\n");
  fprintf(f, "# [4] = charge number Z\n");
  fprintf(f, "# [5] = final abundance (of Z)\n");
  fprintf(f, "# [6] = initial abundance (of Z)\n");

  auto finYAs = out.MakeYVsA(out.FinalY());
  auto finYZs = out.MakeYVsZ(out.FinalY());
  auto initYAs = out.MakeYVsA(out.YVsTime()[0]);
  auto initYZs = out.MakeYVsZ(out.YVsTime()[0]);

  auto num = std::max(finYAs.size(), finYZs.size());

  for (unsigned int i = 0; i < num; ++i) {
    unsigned int a, z;
    double finYA, finYZ, initYA, initYZ;

    if (i < finYAs.size()) {
      a = i;
      finYA = finYAs[a];
      initYA = initYAs[a];
    } else {
      a = 0;
      finYA = 0.0;
      initYA = 0.0;
    }

    if (i < finYZs.size()) {
      z = i;
      finYZ = finYZs[z];
      initYZ = initYZs[z];
    } else {
      z = 0;
      finYZ = 0.0;
      initYZ = 0.0;
    }

    fprintf(f, "%6u  %30.20E  %30.20E  %6u  %30.20E  %30.20E\n", a, finYA,
        initYA, z, finYZ, initYZ);
  }

  fclose(f);
}

NetworkOutput::NetworkOutput(const NetworkOutput& other):
  mIsReadOnly(other.mIsReadOnly),
  mFilePrefix(other.mFilePrefix),
  mFile(other.mFile),
  mLogFile(nullptr),
  mDisableStdoutOutput(other.mDisableStdoutOutput),
  mMassFractionColumnNamesAndIds(other.mMassFractionColumnNamesAndIds),
  mColumnHeader(other.mColumnHeader),
  mpNuclib(new NuclideLibrary(*other.mpNuclib)),
  mMaxA(other.mMaxA),
  mMaxZ(other.mMaxZ),
  mNumEntries(other.mNumEntries),
  mNumEntriesWrittenToFile(other.mNumEntriesWrittenToFile),
  mTimes(other.mTimes),
  mDtVsTime(other.mDtVsTime),
  mTemperatureVsTime(other.mTemperatureVsTime),
  mDensityVsTime(other.mDensityVsTime),
  mEntropyVsTime(other.mEntropyVsTime),
  mHeatingRateVsTime(other.mHeatingRateVsTime),
  mYeVsTime(other.mYeVsTime),
  mEtaEVsTime(other.mEtaEVsTime),
  mYVsTime(other.mYVsTime),
  mA(other.mA),
  mZ(other.mZ),
  mBindingEnergy(other.mBindingEnergy),
  mYVsAVsTime(other.mYVsAVsTime) {
  if (!mIsReadOnly)
    OpenLogFile(true);
}

void NetworkOutput::Close() {
  if (!mIsReadOnly) {
    // this is a file that we are writing to disk

    // write remaining entries to disk
    WriteToDisk();

    // close the datasets
    mTimes.Close();
    mDtVsTime.Close();
    mTemperatureVsTime.Close();
    mDensityVsTime.Close();
    mEntropyVsTime.Close();
    mHeatingRateVsTime.Close();
    mYeVsTime.Close();
    mEtaEVsTime.Close();
    mYVsTime.Close();

    // close the files
    mFile.close();
    fclose(mLogFile);

    mIsReadOnly = true;
  }

  // if this output is read-only, it was read from disk and we don't need to do
  // anything here
}

NetworkOutput::~NetworkOutput() {
  Close();
}

void NetworkOutput::WriteColumnDescriptionToLog() {
  Log("# Column descriptions:\n");
  Log("#\n");
  Log("# [ 1] = step number\n");
  Log("# [ 2] = time [s]\n");
  Log("# [ 3] = dt [s]\n");
  Log("# [ 4] = nuclide that limits dt (or other dt limit reason)\n");
  Log("# [ 5] = number of failed time step attempts\n");
  Log("# [ 6] = number of Newton-Raphson iterations in successful time step\n");
  Log("# [ 7] = Temperature [GK]\n");
  Log("# [ 8] = Density [g cm^{-3}]\n");
  Log("# [ 9] = Ye\n");
  Log("# [10] = deviation from total mass = 1\n");
  for (unsigned int i = 0; i < mMassFractionColumnNamesAndIds.size(); ++i) {
    Log("# [%i] = mass fraction of %s\n", i + 11,
        mMassFractionColumnNamesAndIds[i].first.c_str());
  }
  Log("#\n");
}

void NetworkOutput::AddEntry(const double time,
    const std::pair<double, std::string> dtAndLimitReason, const int numFailed,
    const int numNR, const ThermodynamicState thermoState,
    const double heatingRate, const std::vector<double>& Y,
    const double totalMass) {
  if (mIsReadOnly)
    throw std::runtime_error("Tried to add an entry to a read-only "
        "NetworkOutput.");

  if (Y.size() != (std::size_t)mpNuclib->NumNuclides())
    throw std::invalid_argument("Tried to add an entry with a wrong number of "
        "abundances.");

  mTimes.Data().push_back(time);
  mDtVsTime.Data().push_back(dtAndLimitReason.first);
  mTemperatureVsTime.Data().push_back(thermoState.T9());
  mDensityVsTime.Data().push_back(thermoState.Rho());
  mEntropyVsTime.Data().push_back(thermoState.S());
  mHeatingRateVsTime.Data().push_back(heatingRate);
  mYeVsTime.Data().push_back(thermoState.Ye());
  mEtaEVsTime.Data().push_back(thermoState.EtaElectron());
  mYVsTime.Data().push_back(Y);
  mYVsAVsTime.push_back(MakeYVsA(Y));

  ++mNumEntries;

  if (mNumEntries >= mNumEntriesWrittenToFile + BufferSize)
    WriteToDisk();

  if (mNumEntries == 1)
    Log("%s\n", mColumnHeader.c_str());

  if (mNumEntries % ColumnHeaderPrintFrequency == 0)
    LogWithTimestamp("%s\n", mColumnHeader.c_str());

  // write info to stdout and log file
  Log("% 5i% 12.4e% 12.4e%6s% 4i% 4i% 12.4e% 12.4e% 12.4e% 12.4e", mNumEntries,
      time, dtAndLimitReason.first, dtAndLimitReason.second.c_str(), numFailed,
      numNR, thermoState.T9(), thermoState.Rho(), thermoState.Ye(),
      totalMass - 1.0);
  for (unsigned int i = 0; i < mMassFractionColumnNamesAndIds.size(); ++i) {
    int id = mMassFractionColumnNamesAndIds[i].second;
    double massFraction = Y[id] * (double)mA.Data()[id];
    Log("% 12.4e", massFraction);
  }
  Log("\n");
}

NetworkOutput::NetworkOutput(const std::string& filePrefix,
    const NuclideLibrary& nuclideLibrary, const bool disableStdoutOutput,
    const std::vector<int>& massFractionColumnNuclideIds):
      mIsReadOnly(false),
      mFilePrefix(filePrefix),
      mFile(filePrefix + ".h5", H5F_ACC_TRUNC),
      mLogFile(nullptr),
      mDisableStdoutOutput(disableStdoutOutput),
      mMassFractionColumnNamesAndIds(0),
      mpNuclib(new NuclideLibrary(nuclideLibrary)),
      mNumEntries(0),
      mNumEntriesWrittenToFile(0) {
  for (int id : massFractionColumnNuclideIds) {
    mMassFractionColumnNamesAndIds.push_back({nuclideLibrary.Names()[id], id});
  }

  SetMetadata();

  OpenLogFile(false);

  // initialize H5 data sets
  InitH5Dataset(&mTimes, false);
  InitH5Dataset(&mDtVsTime, false);
  InitH5Dataset(&mTemperatureVsTime, false);
  InitH5Dataset(&mDensityVsTime, false);
  InitH5Dataset(&mEntropyVsTime, false);
  InitH5Dataset(&mHeatingRateVsTime, false);
  InitH5Dataset(&mYeVsTime, false);
  InitH5Dataset(&mEtaEVsTime, false);
  InitH5Dataset(&mYVsTime, true);
  InitH5Dataset(&mA, false);
  InitH5Dataset(&mZ, false);
  InitH5Dataset(&mBindingEnergy, false);

  // write nuclear data to file, since that doesn't change
  if ((std::size_t)mpNuclib->NumNuclides() != nuclideLibrary.As().size())
    throw std::runtime_error("Wrong number of A entries");
  if ((std::size_t)mpNuclib->NumNuclides() != nuclideLibrary.Zs().size())
    throw std::runtime_error("Wrong number of Z entries");
  if ((std::size_t)mpNuclib->NumNuclides()
      != nuclideLibrary.BindingEnergiesInMeV().size())
    throw std::runtime_error("Wrong number of binding energy entries");

  mA.Data() = nuclideLibrary.As();
  mA.WriteToFile();
  mA.Close();

  mZ.Data() = nuclideLibrary.Zs();
  mZ.WriteToFile();
  mZ.Close();

  mMaxA = *std::max_element(mA.Data().begin(), mA.Data().end());
  mMaxZ = *std::max_element(mZ.Data().begin(), mZ.Data().end());

  mBindingEnergy.Data() = nuclideLibrary.BindingEnergiesInMeV();
  mBindingEnergy.WriteToFile();
  mBindingEnergy.Close();

  char buffer[256];
  sprintf(buffer, "%5s%12s%12s%6s%4s%4s%12s%12s%12s%12s", "# num", "time [s]",
      "dt [s]", "limit", "#fl", "#NR", "temp [GK]", "rho [cgs]", "Ye",
      "mass - 1");
  for (unsigned int i = 0; i < mMassFractionColumnNamesAndIds.size(); ++i) {
    std::string str = "X(" + mMassFractionColumnNamesAndIds[i].first + ")";
    sprintf(buffer, "%s%12s", buffer, str.c_str());
  }
  mColumnHeader = std::string(buffer);

  // write log header
  Log("#\n");
  Log("#                                      ##\n");
  Log("#                                    ######\n");
  Log("#                                  ###    ###\n");
  Log("#                                ###   ##   ###\n");
  Log("#                              ###    ####    ###\n");
  Log("#                            ####      ##      ####\n");
  Log("#                          ###  ###          ###  ###\n");
  Log("#                        ###      ###      ###      ###\n");
  Log("#                      ###          ###  ###          ###\n");
  Log("#                    ###              ####              ###\n");
  Log("#                  ###___  _           ##   _  _       _  ###\n");
  Log("#                ### / __|| |__ _  _   ##  | \\| | ___ | |_  ###\n");
  Log("#              ###   \\__ \\| / /| || |  ##  | .` |/ -_)|  _|   ###\n");
  Log("#            ###     |___/|_\\_\\ \\_, |  ##  |_|\\_|\\___| \\__|     ###\n");
  Log("#          ###                  |__/   ##                         ###\n");
  Log("#        ##############################################################\n");
  Log("#\n");
  Log("#\n");
  Log("# SkyNet Nuclear Reaction Network\n");
  Log("# Copyright (c) 2014 - 2017, California Institute of Technology\n");
  Log("# Written by Jonas Lippuner and Luke Roberts\n");
  Log("# All rights reserved\n");
  Log("#\n");
  Log("# SkyNet is free software available at https://bitbucket.org/jlippuner/skynet\n");
  Log("# See the LICENSE file in the source code for conditions. There is NO warranty;\n");
  Log("# not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.\n");
  Log("#\n");
  if (BuildInfo::GitClean) {
    Log("# Git revision: %s\n", BuildInfo::GitHash.c_str());
  } else {
    Log("# Git revision: WORKING DIRECTORY DIRTY\n");
    Log("# There are unstaged changes, uncommitted changes, or untracked files,\n");
    Log("# thus the current Git revision (%s)\n", BuildInfo::GitHash.c_str());
    Log("# is not very meaningful.\n");
    Log("#\n");
  }
  Log("# Built on %s with %s\n", BuildInfo::BuildDateTime.c_str(),
      BuildInfo::MatrixSolver.c_str());
  Log("# SkyNetRoot: %s\n", SkyNetRoot.c_str());
  Log("#\n");
}

NetworkOutput::NetworkOutput(const std::string& h5FilePath) :
    mIsReadOnly(true),
    mFilePrefix(""),
    mFile(h5FilePath, H5F_ACC_RDONLY),
    mLogFile(nullptr),
    mDisableStdoutOutput(true),
    mMassFractionColumnNamesAndIds(0),
    mpNuclib(nullptr),
    mNumEntries(0),
    mNumEntriesWrittenToFile(0) {
  SetMetadata();

  mTimes.ReadFromFile(mFile);
  mDtVsTime.ReadFromFile(mFile);
  mTemperatureVsTime.ReadFromFile(mFile);
  mDensityVsTime.ReadFromFile(mFile);
  mEntropyVsTime.ReadFromFile(mFile);
  mHeatingRateVsTime.ReadFromFile(mFile);
  mYeVsTime.ReadFromFile(mFile);

  // EtaE was added later, so earlier output files do not contain it
  bool hasEtaE = true;
  try {
      mFile.openDataSet(mEtaEVsTime.Name());
  } catch(H5::FileIException& err) {
      hasEtaE = false;
  }

  if (hasEtaE)
    mEtaEVsTime.ReadFromFile(mFile);
  else
    mEtaEVsTime.Data() = std::vector<double>(mTimes.Data().size(), 0.0);

  mYVsTime.ReadFromFile(mFile);
  mA.ReadFromFile(mFile);
  mZ.ReadFromFile(mFile);
  mBindingEnergy.ReadFromFile(mFile);

  // reconstruct nuclide library
  std::vector<std::string> names(mA.Data().size());
  for (unsigned int i = 0; i < names.size(); ++i)
    names[i] = Nuclide::GetName(mZ.Data()[i], mA.Data()[i]);

  mpNuclib = std::unique_ptr<NuclideLibrary>(new NuclideLibrary(
      NuclideLibrary::CreateFromWebnucleoXML(SkyNetRoot +
          "/data/webnucleo_nuc_v2.0.xml", names)));

  mMaxA = *std::max_element(mA.Data().begin(), mA.Data().end());
  mMaxZ = *std::max_element(mZ.Data().begin(), mZ.Data().end());
  mNumEntries = mTimes.Data().size();

  if (mDtVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("dt has wrong number of entries");
  if (mTemperatureVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("temperature has wrong number of entries");
  if (mDensityVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("density has wrong number of entries");
  if (mEntropyVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("entropy has wrong number of entries");
  if (mHeatingRateVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("heating rate has wrong number of entries");
  if (mYeVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("Ye has wrong number of entries");
  if (mEtaEVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("EtaE has wrong number of entries");
  if (mYVsTime.Data().size() != mNumEntries)
    throw std::runtime_error("Y has wrong number of entries");

  for (unsigned int i = 0; i < mNumEntries; ++i) {
    if (mYVsTime.Data()[i].size() != (std::size_t)mpNuclib->NumNuclides())
      throw std::runtime_error("Y[" + std::to_string(i) + "] has wrong number "
          "of nuclides");
  }

  if (mZ.Data().size() != (std::size_t)mpNuclib->NumNuclides())
    throw std::runtime_error("Z has wrong number of nuclides");
  if (mBindingEnergy.Data().size() != (std::size_t)mpNuclib->NumNuclides())
    throw std::runtime_error("binding energy has wrong number of nuclices");

  // make YVsAVsTime
  mYVsAVsTime.resize(mNumEntries);
  for (unsigned int i = 0; i < mNumEntries; ++i) {
    mYVsAVsTime[i] = MakeYVsA(mYVsTime.Data()[i]);
  }

  mFile.close();
}

void NetworkOutput::OpenLogFile(const bool append) {
  if (mLogFile == nullptr) {
    std::string logFilePath = mFilePrefix + ".log";
    if (!File::Writable(logFilePath)) {
      throw std::invalid_argument("The output log file '"
          + logFilePath + "' is not writable.");
    }

    mLogFile = fopen(logFilePath.c_str(), append ? "a" : "w");
  }
}

void NetworkOutput::SetMetadata() {
  mTimes.SetName("Time");
  mTimes.SetDescription("final time of the time step");
  mTimes.SetUnits("second");

  mDtVsTime.SetName("Dt");
  mDtVsTime.SetDescription("time step size");
  mDtVsTime.SetUnits("second");

  mTemperatureVsTime.SetName("Temperature");
  mTemperatureVsTime.SetDescription(
      "temperature at the end of the time step (if the network is self-heating, "
      "the temperature is calculated from the entropy at the beginning of the "
      "time step and the density at the end of the time step)");
  mTemperatureVsTime.SetUnits("Giga Kelvin");

  mDensityVsTime.SetName("Density");
  mDensityVsTime.SetDescription("density at the end of the time step");
  mDensityVsTime.SetUnits("g / cm^3");

  mEntropyVsTime.SetName("Entropy");
  mEntropyVsTime.SetDescription("if the network is self-heating, this is the "
      "entropy at the end of the time step, otherwise it is the change in "
      "entropy during the time step");
  mEntropyVsTime.SetUnits("kb (Boltzmann constant) / baryon");

  mHeatingRateVsTime.SetName("HeatingRate");
  mHeatingRateVsTime.SetDescription(
      "amount of specific heating (due to nuclear reactions and neutrino "
      "heating/cooling) per unit time that occurred in the time step");
  mHeatingRateVsTime.SetUnits("erg / s / g");

  mYeVsTime.SetName("Ye");
  mYeVsTime.SetDescription("electron fraction at the end of the time step");
  mYeVsTime.SetUnits("dimensionless");

  mEtaEVsTime.SetName("EtaE");
  mEtaEVsTime.SetDescription("electron degeneracy factor at the end of the time step");
  mEtaEVsTime.SetUnits("dimensionless");

  mYVsTime.SetName("Y");
  mYVsTime.SetDescription(
      "abundances (number density of nuclear species / baryon number density) "
      "of all nuclides at the end of the time step, the first index is the "
      "time step and the second index is the nuclide index");
  mYVsTime.SetUnits("dimensionless");
  if (mpNuclib != nullptr)
    mYVsTime.SetDim2Size(mpNuclib->NumNuclides());
  else
    mYVsTime.SetDim2Size(0);

  mA.SetName("A");
  mA.SetDescription("mass number A");
  mA.SetUnits("dimensionless");

  mZ.SetName("Z");
  mZ.SetDescription("atomic number Z");
  mZ.SetUnits("dimensionless");

  mBindingEnergy.SetName("BindingEnergy");
  mBindingEnergy.SetDescription("nuclear binding energy");
  mBindingEnergy.SetUnits("MeV");
}

template<typename T>
void NetworkOutput::InitH5Dataset(NetworkOutputDataset<T> * const pDataset,
    const bool reduceCache) {
  pDataset->InitForWriting(&mFile, ChunkSize1D, ChunkSize2D, reduceCache);
}

void NetworkOutput::WriteToDisk() {
  mTimes.WriteToFile();
  mDtVsTime.WriteToFile();
  mTemperatureVsTime.WriteToFile();
  mDensityVsTime.WriteToFile();
  mEntropyVsTime.WriteToFile();
  mHeatingRateVsTime.WriteToFile();
  mYeVsTime.WriteToFile();
  mEtaEVsTime.WriteToFile();
  mYVsTime.WriteToFile();

  mNumEntriesWrittenToFile = mNumEntries;

  mFile.flush(H5F_SCOPE_GLOBAL);
}

std::vector<double> NetworkOutput::MakeYVsA(
    const std::vector<double>& Y) const {
  std::vector<double> yVsA(mMaxA + 1, 0.0);
  for (unsigned int i = 0; i < Y.size(); ++i)
    yVsA[mA.Data()[i]] += Y[i];

  return yVsA;
}

std::vector<double> NetworkOutput::MakeYVsZ(
    const std::vector<double>& Y) const {
  std::vector<double> yVsZ(mMaxZ + 1, 0.0);
  for (unsigned int i = 0; i < Y.size(); ++i)
    yVsZ[mZ.Data()[i]] += Y[i];

  return yVsZ;
}

