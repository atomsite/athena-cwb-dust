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
    // Imported variables
    Real x1,x2,x3;   // Position (cm)
    Real v1,v2,v3;   // Velocity (cm s^-1)
    Real mdot;       // Mass loss rate (msol yr^-1)
    Real vinf;       // Terminal velocity (cm s^-1)
    Real twnd;       // Wind temperature (K)
    Real X,Y,Z;      // Metallicities
    Real remap;      // Remap radius (cells)
    // Calculated variables
    Real mdot_cgs;   // Mass loss rate (g s^-1)
    Real remap_cgs;  // Remap pradius (cm)
    Real remap_vol;  // Remap volume (cm^3)
    Real mu;         // Mean molecular mass
    Real avgm;       // Mean molecular mass (g)
    // Functions
    int Init(ParameterInput *pin , std::string id);
};

int Star::Init(ParameterInput *pin, std::string id){
  // Read from model file:
  // Get star position
  x1 = pin->GetReal("problem","x1_"+id);
  x2 = pin->GetReal("problem","x2_"+id);
  x3 = pin->GetReal("problem","x3_"+id);
  // Get wind properties
  mdot = pin->GetReal("problem","mdot_"+id);
  vinf = pin->GetReal("problem","vinf_"+id);
  twnd = pin->GetReal("problem","twnd_"+id);
  // Get metallicities
  X = pin->GetReal("problem","X_"+id);
  Y = pin->GetReal("problem","Y_"+id);
  Z = pin->GetReal("problem","Z_"+id);
  // Get remap radius
  remap = pin->GetReal("problem","remap_"+id);

  // Initialise star velocity
  v1 = 0.0;
  v2 = 0.0;
  v3 = 0.0;

  // Calculate and/or convert additional variables

  // Convert variables from input units to CGS
  // Calculate mass loss rate in CGS
  Real yr_to_sec = 3.154e+7;
  Real msol_to_g = 1.989e+33;
  mdot_cgs       = mdot * msol_to_g / yr_to_sec;
  // Calculate remap radius in CGS
  int  nx1      = pin->GetInteger("mesh","nx1");
  int  nlevs    = pin->GetInteger("mesh","numlevel");
  Real x1l      = pin->GetReal("mesh","x1min");
  Real x1r      = pin->GetReal("mesh","x1max");
  Real x1coarse = (x1r-x1l) / nx1;
  Real x1fine   = x1coarse / std::pow(2.0,nlevs-1);
  remap_cgs     = x1coarse * remap;
  remap_vol     = (4.0 / 3.0) * PI * std::pow(remap_cgs,3.0);
  // Calculate mean molecular mass
  Real massh = 1.67355e-24;
  mu         = 1.0 / ((2 * X) + (0.75 * Y) + (0.5 * Z));
  avgm       = mu*massh;
}