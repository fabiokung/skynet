%{
#include "BuildInfo.hpp"
%}

%rename(BuildInfo_GitHash) BuildInfo::GitHash;
%rename(BuildInfo_GitClean) BuildInfo::GitClean;
%rename(BuildInfo_BuildDateTime) BuildInfo::BuildDateTime;
%rename(BuildInfo_MatrixSolver) BuildInfo::MatrixSolver;

%include "BuildInfo.hpp"