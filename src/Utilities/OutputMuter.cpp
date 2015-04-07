/// \file OutputMuter.cpp
/// \author jlippuner
/// \since Jul 11, 2014
///
/// \brief
///
///

#include "Utilities/OutputMuter.hpp"

#include <unistd.h>

void OutputMuter::Mute() {
  if (!mIsCurrentlyMuted) {
    // save position of current stdout
    fgetpos(stdout, &mOldStdoutPos);
    // duplicate the file descriptor of current stdout
    mOldStdoutFileDescriptor = dup(fileno(stdout));
    // redirect stdout to /dev/null
    freopen("/dev/null", "w", stdout);
    mIsCurrentlyMuted = true;
  }
}

void OutputMuter::Unmute() {
  if (mIsCurrentlyMuted) {
    // flush stdout to /dev/null
    fflush(stdout);
    // close output to /dev/null and restore original stdout
    dup2(mOldStdoutFileDescriptor, fileno(stdout));
    // deallocate the file descriptor that stored the original stdout
    close(mOldStdoutFileDescriptor);
    // clear errors from stdout
    clearerr(stdout);
    // restore original position in stdout
    fsetpos(stdout, &mOldStdoutPos);

    mIsCurrentlyMuted = false;
  }
}

bool OutputMuter::mIsCurrentlyMuted;
fpos_t OutputMuter::mOldStdoutPos;
int OutputMuter::mOldStdoutFileDescriptor;
