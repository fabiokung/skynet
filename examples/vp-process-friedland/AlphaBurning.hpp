/// \file r-process.cpp
/// \author jlippuner
/// \since Sep 3, 2014
///
/// \brief
///
///
#ifndef ALPHA_BURNING_HPP_
#define ALPHA_BURNING_HPP_

#include <iostream> 
#include <math.h>
#include <vector>
#include "Reactions/Reaction.hpp" 
#include "Reactions/SpecialReactionLibrary.hpp" 

namespace AlphaBurning {  
  
  auto REACLIBFunc = [] (const std::vector<double>& aa, const double& T9) { 
    double T9onethird = pow(T9, 1.0/3.0); 
    double arg = 
       aa[0] 
     + aa[1] / T9
     + aa[2] / T9onethird
     + aa[3] * T9onethird
     + aa[4] * T9 
     + aa[5] * T9 * T9onethird * T9onethird
     + aa[6] * log(T9); 
    return exp(std::min(arg,1.e10)); 
  };

  double ReaclibRate(double /*time*/, ThermodynamicState state, 
      int nparticleEntrance, double doubleCount, 
      std::vector<std::vector<double>> params) { 
    double rate = 0.0;
    double T9 = std::min(std::max(state.T9(), 0.01), 10.0);  
    for (const auto &aa : params) 
      rate += REACLIBFunc(aa, T9);
    return rate * pow(state.Rho(), nparticleEntrance - 1)/doubleCount;
  }
   
  const std::vector<std::vector<double>> TripleAlpha = {{ 
    -2.435050e+01,-4.126560e+00,-1.349000e+01, 2.142590e+01,
    -1.347690e+00, 8.798160e-02,-1.316530e+01},{ 
    -1.178840e+01,-1.024460e+00,-2.357000e+01, 2.048860e+01,                      
    -1.298820e+01,-2.000000e+01,-2.166670e+00},{
    -9.710520e-01, 0.000000e+00,-3.706000e+01, 2.934930e+01,
    -1.155070e+02,-1.000000e+01,-1.333330e+00}};
  
  const std::vector<std::vector<double>> TripleAlphaInverse = {{ 
     4.577340e+01,-8.442270e+01,-3.706000e+01, 2.934930e+01,
    -1.155070e+02,-1.000000e+01, 1.666670e+00},{                                   
     2.239400e+01,-8.854930e+01,-1.349000e+01, 2.142590e+01,                     
    -1.347690e+00, 8.798160e-02,-1.016530e+01},{                                   
     3.495610e+01,-8.544720e+01,-2.357000e+01, 2.048860e+01,                     
    -1.298820e+01,-2.000000e+01, 8.333300e-01}};                                   

  //const std::vector<std::vector<double>> AlphaAlphaN = {{ 
  //  -8.004650e+00,-1.019530e+00,-1.476090e+00,-5.481100e+00,
  //   9.317100e-01,-7.803630e-02,-1.500000e+00},{
  //  -7.195380e+00, 0.000000e+00,-1.323870e+01, 7.156830e+00,                      
  //  -1.732780e+00,-3.000000e+00,-6.666670e-01}}; 
  //                                    
  //const std::vector<std::vector<double>> AlphaAlphaNInverse = {{ 
  //   3.749730e+01,-1.927920e+01,-1.476090e+00,-5.481100e+00,                      
  //   9.317100e-01,-7.803630e-02, 1.500000e+00},{                                   
  //   3.830660e+01,-1.825970e+01,-1.323870e+01, 7.156830e+00,                      
  //  -1.732780e+00,-3.000000e+00, 2.333330e+00}};                  
  

  // AC12 alphaalphaN 
  const std::vector<std::vector<double>> AlphaAlphaN = {{ 
    -8.228980e+00, 0.000000e+00,-1.333170e+01, 1.322370e+01,
    -9.063390e+00, 0.000000e+00,-6.666670e-01},{
    -6.811780e+00,-1.019530e+00,-1.566730e+00,-5.434970e+00, 
     6.738070e-01,-4.101400e-02,-1.500000e+00}};
  
  // AC12 be9 -> alpha alpha n
  const std::vector<std::vector<double>> AlphaAlphaNInverse = {{ 
     3.727300e+01,-1.825970e+01,-1.333170e+01, 1.322370e+01,
    -9.063390e+00, 0.000000e+00, 2.333330e+00},{
     3.869020e+01,-1.927920e+01,-1.566730e+00,-5.434970e+00, 
     6.738070e-01,-4.101400e-02, 1.500000e+00}};
  
  const std::vector<std::vector<double>> Be9Alpha = {{ 
     1.174400e+01,-4.179000e+00, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00,-1.500000e+00},{                                   
    -1.482810e+00,-1.834000e+00, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00,-1.500000e+00},{                                   
    -9.519590e+00,-1.184000e+00, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00,-1.500000e+00},{                                   
     3.146400e+01, 0.000000e+00,-2.387000e+01, 5.666980e-01,                      
     4.409570e+01,-3.142320e+02,-6.666670e-01,},{                                   
     1.929620e+01,-1.273200e+01, 0.000000e+00, 0.000000e+00,                   
     0.000000e+00, 0.000000e+00, 0.000000e+00}};
  
  const std::vector<std::vector<double>> Be9AlphaInverse = {{ 
     8.582560e-01,-6.799130e+01, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00,-1.500000e+00},{                                   
    -7.178520e+00,-6.734130e+01, 0.000000e+00, 0.000000e+00,                     
     0.000000e+00, 0.000000e+00,-1.500000e+00},{                                   
     3.380510e+01,-6.615730e+01,-2.387000e+01, 5.666980e-01,                      
     4.409570e+01,-3.142320e+02,-6.666670e-01},{                                   
     2.163730e+01,-7.888930e+01, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00, 0.000000e+00},{                                   
     1.408510e+01,-7.033630e+01, 0.000000e+00, 0.000000e+00,                      
     0.000000e+00, 0.000000e+00,-1.500000e+00}};                                   


  std::vector<std::pair<Reaction, RateFunc>> GetRates(const NuclideLibrary& nuclib, 
      bool enhanced, double Be9Fac,  
      double enhancePEnhance, double enhanceNEnhance) {
    std::vector<std::pair<Reaction, RateFunc>> rates; 
  
    auto rateTripleAlpha = std::pair<Reaction, RateFunc>(
        Reaction({"he4"},{"c12"},{3},{1}, false, false, false, 
        "Triple alpha", nuclib),
        [](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 3, 6.0, TripleAlpha);
        });
    rates.push_back(rateTripleAlpha); 
    
    auto rateTripleAlphaInverse = std::pair<Reaction, RateFunc>(
        Reaction({"c12"},{"he4"},{1},{3}, false, true, false, 
        "Triple alpha inverse", nuclib),
        [](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 1, 1.0, TripleAlphaInverse); 
        });
    rates.push_back(rateTripleAlphaInverse); 
  
    auto rateTripleAlphaNeutron = std::pair<Reaction, RateFunc>(
        Reaction({"n","he4"},{"n","c12"},{1,3},{1,1}, false, false, false, 
        "Neutron Catalyzed Triple alpha", nuclib),
        [enhanceNEnhance](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 3, 6.0, TripleAlpha)
	        * enhanceNEnhance
	        * state.Rho() * 1.e-6 * (75.1*exp(-state.T9()) + 88.7);
        });
    if (enhanced) rates.push_back(rateTripleAlphaNeutron); 
    
    auto rateTripleAlphaNeutronInverse = std::pair<Reaction, RateFunc>(
        Reaction({"n","c12"},{"n","he4"},{1,1},{1,3}, false, true, false, 
        "Triple alpha inverse", nuclib),
        [enhanceNEnhance](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 1, 1.0, TripleAlphaInverse)
	        * enhanceNEnhance
	        * state.Rho() * 1.e-6 * (75.1*exp(-state.T9()) + 88.7);
        });
    if (enhanced) rates.push_back(rateTripleAlphaNeutronInverse); 

    auto rateTripleAlphaProton = std::pair<Reaction, RateFunc>(
        Reaction({"p","he4"},{"p","c12"},{1,3},{1,1}, false, false, false, 
        "Neutron Catalyzed Triple alpha", nuclib),
        [enhancePEnhance](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 3, 6.0, TripleAlpha)
	        * enhancePEnhance
	        * state.Rho() * 1.e-6 
                *(0.03680 - 1.667*state.T9()+ 2.350*pow(state.T9(),2) 
                - 0.2911* pow(state.T9(),3)+ 0.01160*pow(state.T9(),4));
        });
    if (enhanced) rates.push_back(rateTripleAlphaProton); 
    
    auto rateTripleAlphaProtonInverse = std::pair<Reaction, RateFunc>(
        Reaction({"p","c12"},{"p","he4"},{1,1},{1,3}, false, true, false, 
        "Triple alpha inverse", nuclib),
        [enhancePEnhance](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 1, 1.0, TripleAlphaInverse)
	        * enhancePEnhance
                * state.Rho() * 1.e-6 
	        * (0.03680 - 1.667*state.T9()+ 2.350*pow(state.T9(),2) 
                - 0.2911* pow(state.T9(),3)+ 0.01160*pow(state.T9(),4));
        });
    if (enhanced) rates.push_back(rateTripleAlphaProtonInverse);

    auto rateAlphaAlphaN = std::pair<Reaction, RateFunc>(
        Reaction({"n","he4"},{"be9"},{1,2},{1}, false, false, false, 
        "2alpha + n", nuclib),
        [](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 3, 2.0, AlphaAlphaN); 
        });
    rates.push_back(rateAlphaAlphaN); 
  
    auto rateAlphaAlphaNInverse = std::pair<Reaction, RateFunc>(
        Reaction({"be9"}, {"n","he4"}, {1}, {1,2}, false, true, false, 
        "2alpha + n inverse", nuclib),
        [](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 1, 1.0, AlphaAlphaNInverse); 
        });
    rates.push_back(rateAlphaAlphaNInverse); 
    
    auto rateBe9Alpha = std::pair<Reaction, RateFunc>(
        Reaction({"he4","be9"},{"n","c12"},{1,1},{1,1}, false, false, false, 
        "Be9 + alpha", nuclib),
        [Be9Fac](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 2, 1.0, Be9Alpha)*Be9Fac; 
         });
    rates.push_back(rateBe9Alpha); 
  
    auto rateBe9AlphaInverse = std::pair<Reaction, RateFunc>(
        Reaction({"n","c12"},{"he4","be9"},{1,1},{1,1}, false, true, false, 
        "Be9 + alpha Inverse", nuclib),
        [Be9Fac](double time, ThermodynamicState state){
          return ReaclibRate(time, state, 2, 1.0, Be9AlphaInverse)*Be9Fac; 
        });
    rates.push_back(rateBe9AlphaInverse); 
    
    return rates;  
  }
  

}

#endif // ALPHA_BURNING_HPP_

  
  
  
  
  
  
  

  
  
  
  
  
  
  
  
  
  
  
