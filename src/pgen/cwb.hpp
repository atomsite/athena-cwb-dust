// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

//! \class Star
//  \brief Class used to describe a star in a CWB problem

class Star {
  public:
    Real x1,x2,x3;
    Real v1,v2,v3;
    Real mdot;
    Real vinf;
    Real twnd;
    Real X,Y,Z;

    int Init(ParameterInput *pin , std::string id);
};

int Star::Init(ParameterInput *pin, std::string id){
  // Get star position from model file
  x1 = pin->GetReal("problem","x1_"+id);
  x2 = pin->GetReal("problem","x2_"+id);
  x3 = pin->GetReal("problem","x3_"+id);
  // Initialise wind velocity
  v1 = 0.0;
  v2 = 0.0;
  v3 = 0.0;
  // Get wind properties from model file
  mdot = pin->GetReal("problem","mdot_"+id);
  vinf = pin->GetReal("problem","vinf_"+id);
  twnd = pin->GetReal("problem","twnd_"+id);
  // Get metallicities from model file
  X = pin->GetReal("problem","X_"+id);
  Y = pin->GetReal("problem","Y_"+id);
  Z = pin->GetReal("problem","Z_"+id);
}