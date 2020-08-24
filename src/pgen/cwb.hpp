// C headers

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <sstream>
#include <stdexcept>
#include <string>
#include <fstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

#define CART 0
#define CYL  1

//! \class Star
//  \brief Class used to describe a star in a CWB problem

class Star {
  public:
    // Imported variables
    Real x[3];   // Position (cm)
    Real v[3];   // Velocity (cm s^-1)
    Real mass;
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
    // Scalars
    int  scal;       // Identifying scalar, 1 for WR wind, 0 for OB
    // Functions
    int Init(ParameterInput *pin , std::string id, int scal_id);
};

int Star::Init(ParameterInput *pin, std::string id, int scal_id){
  // Read from model file:
  // Get star properties
  x[0] = pin->GetReal("problem","x1_"+id);
  x[1] = pin->GetReal("problem","x2_"+id);
  x[2] = pin->GetReal("problem","x3_"+id);

  mass = pin->GetReal("problem","mass_"+id);
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
  v[0] = 0.0;
  v[1] = 0.0;
  v[2] = 0.0;

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
  // Set scalar value, which is used to identify wind
  scal = scal_id;
}

//! \class CoolCurve
//  \brief 

class CoolCurve {
  public:
    // Variables
    std::vector<Real> temp;   // Temperature (K)
    std::vector<Real> lambda;  // Energy loss constant
    Real tmin, tmax;  // Values of first and last index in temp
    Real gmin, gmax;  // Values of first and last index in lambda
    // Functions
    int  Init(std::string coolCurveFileName);
    Real FindLambda(Real t);
};

//! \fn Init
//  \brief Initialise cooling curve by reading in a space separated file
//  each line containing a pair of values, temperature and lambda.
//  Each line is then read into the appropriate vector.
//  Function checks for validity of file, as well as simple check to determine
//  if vectors are the same length.

int CoolCurve::Init(std::string coolCurveFileName) {
  std::ifstream file(coolCurveFileName);
  if (!file) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.hpp CoolCurve::Init" <<
           std::endl <<
           "File " <<
           coolCurveFileName <<
           " Cannot be read!" <<
           std::endl;
    ATHENA_ERROR(msg);
    return 1;
  }
  else {
    Real tbuf,lbuf;
    std::string line;
    file.seekg(0);
    while (!file.eof()) {
      file >> tbuf;
      file >> lbuf;
      if (lbuf == tbuf) {break;}
      temp.push_back(tbuf);
      lambda.push_back(lbuf);
    }
  }

  // Convert log(T) to T
  for (int n = 0; n < temp.size(); n++) {
    temp[n] = pow(10.0,temp[n]);
  }

  // Compare lengths of vectors, if they are not the same, import has failed
  if (temp.size() != lambda.size()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.hpp CoolCurve::Init" <<
           std::endl <<
           "Import has failed, mismatch in temperature and lambda vector size" <<
           std::endl;
    ATHENA_ERROR(msg);
  }

  tmin = temp.front();
  tmax = temp.back();
  gmin = lambda.front();
  gmax = lambda.back();

  return 0;
}

//! \fn FindLambda
// Perform a binary search to find first temperature in cooling curve > input temperature, tu
// Function then performs a linear interpolation between tu and the previous temperature in index, tl, to find the interpolated value for lambda
// Fast, but not particularly accurate at lower temperatures
// Uses C++ algorithm library
Real CoolCurve::FindLambda(Real t) {
  // Quickly return values that would cause search to fail
  if      (t <= tmin) {return 0.0;}
  else if (t >= tmax) {return gmax;}
  // Perform binary search
  auto upper = std::upper_bound(temp.begin(), temp.end(), t);
  // Assign indexes based on search
  int  iu    = std::distance(temp.begin(), upper);
  int  il    = iu - 1;
  // Perform linear interpolation
  // Assign values from each array from indexes found by binary search
  Real tl = temp[il];
  Real tu = temp[iu];
  Real gl = lambda[il];
  Real gu = lambda[iu];
  // Calculate interpolated value of lambda
  Real g  = gl + ((t - tl) * ((gu - gl) / (tu - tl)));
  return g;
}