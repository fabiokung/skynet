%{
#include "Reactions/ReactionData.hpp"
%}

%ignore ReactionData::operator[];
%ignore ReactionData::operator=;

%include "Reactions/ReactionData.hpp"
