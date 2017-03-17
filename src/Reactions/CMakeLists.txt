set(Reactions_SOURCES
  ConstReactionLibrary.cpp
  Reaction.cpp
  ReactionLibraryBase.cpp
  REACLIB.cpp
  REACLIBEntry.cpp
  REACLIBReactionLibrary.cpp
  FFN.cpp
  FFNEntry.cpp
  FFNReactionLibrary.cpp
  Neutrino.cpp 
  NeutrinoEntry.cpp
  NeutrinoReactionLibrary.cpp
  SpecialReactionLibrary.cpp
  NeutrinoExternalRatesReactionLibrary.cpp
  ReactionPostProcess.cpp
)

add_SkyNet_library(Reactions "${Reactions_SOURCES}")

set(Reactions_SWIG_DEPS
  ConstReactionLibrary.i
  ConstReactionLibrary.hpp
  REACLIB.i
  REACLIB.hpp
  REACLIBEntry.i
  REACLIBEntry.hpp
  REACLIBReactionLibrary.i
  REACLIBReactionLibrary.hpp
  Reaction.i
  Reaction.hpp
  ReactionLibraryBase.hpp
  FFN.hpp
  FFNReactionLibrary.i
  FFNReactionLibrary.hpp
  Neutrino.i
  Neutrino.hpp
  NeutrinoExternalRatesReactionLibrary.i
  NeutrinoExternalRatesReactionLibrary.hpp
  NeutrinoReactionLibrary.i
  NeutrinoReactionLibrary.hpp
  ReactionPostProcess.i
  ReactionPostProcess.hpp
)

add_SWIG_dependencies("${Reactions_SWIG_DEPS}")
