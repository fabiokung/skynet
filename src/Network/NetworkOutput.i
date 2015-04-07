%{
#include "Network/NetworkOutput.hpp"
%}

%ignore Log;
%ignore LogWithTimestamp;
%ignore WriteColumnDescriptionToLog;
%ignore AddEntry;
%ignore BufferSize;

%include "Network/NetworkOutput.hpp"
