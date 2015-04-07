/// \file NucleiChartOptions.hpp
/// \author jonas
/// \since Jun 12, 2014
///
/// \brief
///
///

#ifndef SKYNET_MOVIEMAKER_NUCLEICHARTOPTIONS_HPP_
#define SKYNET_MOVIEMAKER_NUCLEICHARTOPTIONS_HPP_

#include <array>

namespace Opt {

static constexpr double fontSize = 14.0;
static constexpr const char* fontFamily = "sans-serif";

static constexpr double cellSize = 5.0;
static constexpr double lineWidth = 2.0;
static constexpr double lowerLeftBorder = 3.5; // in number of cells
static constexpr double upperRightBorder = 1.5; // in number of cells

static constexpr double colorScaleWidth = 20.0;
static constexpr double colorScaleLeftBorder = 20.0;
static constexpr double colorScaleTopBorder = 80.0;

static constexpr std::array<double, 3> backgroundRGB {{ 1.0, 1.0, 1.0 }};
static constexpr std::array<double, 3> borderRGB {{ 0.8, 0.8, 0.8 }};
static constexpr std::array<double, 3> stableBorderRGB {{ 0.0, 0.0, 0.0 }};
static constexpr std::array<double, 3> peakFillRGB = borderRGB;
static constexpr std::array<double, 4> missingFillRGBA {{ 0.3, 0.3, 0.3, 1.0 }};

// see ../stable_nuclides_notes on how these were found
static constexpr std::array<int, 289> stableN = {{ 0, 1, 1, 2, 3, 4, 5, 5, 6,
    6, 7, 7, 8, 8, 9, 10, 10, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16, 16,
    16, 17, 18, 20, 18, 20, 18, 20, 22, 20, 21, 22, 20, 22, 23, 24, 26, 28, 24,
    24, 25, 26, 27, 28, 27, 28, 26, 28, 29, 30, 30, 28, 30, 31, 32, 32, 30, 32,
    33, 34, 36, 34, 36, 34, 36, 37, 38, 40, 38, 40, 38, 40, 41, 42, 44, 42, 40,
    42, 43, 44, 46, 48, 44, 46, 42, 44, 46, 47, 48, 50, 48, 50, 46, 48, 49, 50,
    50, 50, 51, 52, 54, 56, 51, 52, 50, 52, 53, 54, 55, 56, 58, 52, 54, 55, 56,
    57, 58, 60, 58, 56, 58, 59, 60, 62, 64, 60, 62, 58, 60, 62, 63, 64, 65, 66,
    68, 64, 66, 62, 64, 65, 66, 67, 68, 69, 70, 72, 74, 70, 72, 68, 70, 71, 72,
    73, 74, 76, 78, 74, 70, 72, 74, 75, 76, 77, 78, 80, 82, 78, 74, 76, 78, 79,
    80, 81, 82, 81, 82, 78, 80, 82, 84, 82, 82, 83, 84, 85, 86, 88, 90, 82, 84,
    85, 86, 87, 88, 90, 92, 88, 90, 88, 90, 91, 92, 93, 94, 96, 94, 90, 92, 94,
    95, 96, 97, 98, 98, 94, 96, 98, 99, 100, 102, 100, 98, 100, 101, 102, 103,
    104, 106, 104, 105, 102, 104, 105, 106, 107, 108, 107, 108, 106, 108, 109,
    110, 112, 110, 112, 108, 110, 111, 112, 113, 114, 116, 114, 116, 112, 114,
    116, 117, 118, 120, 118, 116, 118, 119, 120, 121, 122, 124, 122, 124, 122,
    124, 125, 126, 126, 142, 143, 146, 150 }};

static constexpr std::array<int, 289> stableZ = {{ 1, 1, 2, 2, 3, 3, 4, 5, 5,
    6, 6, 7, 7, 8, 8, 8, 9, 10, 10, 10, 11, 12, 12, 12, 13, 14, 14, 14, 15, 16,
    16, 16, 16, 17, 17, 18, 18, 18, 19, 19, 19, 20, 20, 20, 20, 20, 20, 21, 22,
    22, 22, 22, 22, 23, 23, 24, 24, 24, 24, 25, 26, 26, 26, 26, 27, 28, 28, 28,
    28, 28, 29, 29, 30, 30, 30, 30, 30, 31, 31, 32, 32, 32, 32, 32, 33, 34, 34,
    34, 34, 34, 34, 35, 35, 36, 36, 36, 36, 36, 36, 37, 37, 38, 38, 38, 38, 39,
    40, 40, 40, 40, 40, 41, 41, 42, 42, 42, 42, 42, 42, 42, 44, 44, 44, 44, 44,
    44, 44, 45, 46, 46, 46, 46, 46, 46, 47, 47, 48, 48, 48, 48, 48, 48, 48, 48,
    49, 49, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 51, 51, 52, 52, 52, 52, 52,
    52, 52, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 56, 56, 56, 56, 56,
    56, 56, 57, 57, 58, 58, 58, 58, 59, 60, 60, 60, 60, 60, 60, 60, 62, 62, 62,
    62, 62, 62, 62, 62, 63, 63, 64, 64, 64, 64, 64, 64, 64, 65, 66, 66, 66, 66,
    66, 66, 66, 67, 68, 68, 68, 68, 68, 68, 69, 70, 70, 70, 70, 70, 70, 70, 71,
    71, 72, 72, 72, 72, 72, 72, 73, 73, 74, 74, 74, 74, 74, 75, 75, 76, 76, 76,
    76, 76, 76, 76, 77, 77, 78, 78, 78, 78, 78, 78, 79, 80, 80, 80, 80, 80, 80,
    80, 81, 81, 82, 82, 82, 82, 83, 90, 92, 92, 94 }};

// closed shell neutron numbers: 2, 8, 20, 28, 50, 82, 126
static constexpr std::array<int, 7> magicN {{ 2, 8, 20, 28, 50, 82, 126 }};
static constexpr int magicNLineOvershoot = 0;
static constexpr std::array<double, 3> magicNumberLineRGB {{ 0.0, 0.0, 0.0 }};

static constexpr std::array<int, 3> peakA {{ 80, 130, 194 }};
static constexpr int peakWidth = 4;
static constexpr int peakLineOvershoot = 0;
static constexpr std::array<double, 3> peakLineRGB {{ 0.0, 0.0, 0.0 }};

static constexpr std::array<double, 3> pathLineRGB {{ 0.0, 0.0, 0.0 }};

class str_const {  // constexpr string
private:
  const char* const p_;
  const std::size_t sz_;
public:
  template<std::size_t N>
  constexpr str_const(const char (&a)[N]) :  // ctor
      p_(a),
      sz_(N - 1) { }

  constexpr char operator[](std::size_t n) { // []
    return n < sz_ ? p_[n] : throw std::out_of_range(
        std::to_string(n) + " >= " + std::to_string(sz_));
  }

  constexpr std::size_t size() { return sz_; } // size()
};


// histogram
static constexpr double binWdith = 2.0;
static constexpr std::array<int, 10> leftBinEdges {{ 1, 2, 72, 94, 111, 145,
    170, 188, 206, 253 }};
static constexpr std::array<const char*, 10> binNames {{ "neutrons", "",
  "1^{st} peak", "", "2^{nd} peak  ", "  rare-earth", "", "3^{rd} peak", "",
  "fission material" }};
static constexpr std::array<std::array<double, 4>, 10> binColors {{
  {{ 1.0, 1.0, 1.0, 1.0 }}, // neutrons
  {{ 1.0, 1.0, 1.0, 1.0 }}, //
  {{ 0.8, 1.0, 0.8, 1.0 }}, // 1st peak
  {{ 1.0, 1.0, 1.0, 1.0 }}, //
  {{ 1.0, 1.0, 0.8, 1.0 }}, // 2nd peak
  {{ 1.0, 0.9, 0.8, 1.0 }}, // rare-earth
  {{ 1.0, 1.0, 1.0, 1.0 }}, //
  {{ 1.0, 0.85, 0.85, 1.0 }}, // 3rd peak
  {{ 1.0, 1.0, 1.0, 1.0 }}, //
  {{ 0.9, 0.9, 0.9, 1.0 }}  // fission material
}};

} // namespace Opt

#endif // SKYNET_MOVIEMAKER_NUCLEICHARTOPTIONS_HPP_
