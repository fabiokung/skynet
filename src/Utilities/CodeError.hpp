#ifndef SKYNET_UTILITIES_CODEERROR_HPP_
#define SKYNET_UTILITIES_CODEERROR_HPP_

#include <sstream>
#include <string>

// Function to handle ASSERT failure. Should never return.
void AssertFail(const char* expression, const char* file, const int line,
    std::string msg) __attribute__((noreturn));

// ASSERT signals that an expression should be true. If the preprocessor macro
// DEBUG is defined and the expression is false, an error is printed and the
// program exits.
// It contains a void cast to avoid an unused variable compiler warning when a
// passed parameter is only used in an ASSERT. This is inside an if(false)
// statement so that the code never actually gets executed and only placates
// the compiler.
#ifdef DEBUG
#define ASSERT(a,m) \
do { if (!(a)) { \
  std::ostringstream amsg; amsg << m; \
  AssertFail(#a, __FILE__, __LINE__, amsg.str()); \
} } while (0)
#else
#define ASSERT(a,m) do { if (false) (void)(a); } while (0)
#endif

/// Marker-class for code errors.  Currently does nothing.
class CodeError {
public:
  CodeError() { }
};

#endif // SKYNET_UTILITIES_CODEERROR_HPP_
