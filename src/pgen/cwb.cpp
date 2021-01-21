//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cwb.cpp
//  \brief Problem generator for CWB problem
//
// The CWB problem consists of two supersonic winds that collide together
// Input parameters are:
//    - problem/mdot1  = mass-loss rate of star1 (Msol/yr)
//    - problem/mdot2  = mass-loss rate of star2 (Msol/yr)
//    - problem/vinf1  = terminal wind speed of star1 (cm/s)
//    - problem/vinf2  = terminal wind speed of star2 (cm/s)
//    - problem/xpos1  = coordinate-1 position of star1 (cm)
//    - problem/ypos1  = coordinate-2 position of star1 (cm)
//    - problem/zpos1  = coordinate-3 position of star1 (cm)
//    - problem/xpos2  = coordinate-1 position of star2 (cm)
//    - problem/ypos2  = coordinate-2 position of star2 (cm)
//    - problem/zpos2  = coordinate-3 position of star2 (cm)
//
// Orbital motion is only included for cartesian calculations (not cylindrically symmetric).
// Wind colour uses the first advected scalar (index 0).
// Dust uses the second and third advected scalars (indices 1 and 2).
//
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // upper_bound() used for fast binary search
#include <cmath>      // sqrt()
#include <fstream>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"   // JMP: needed for scalars

// Cube and power-4 macros, significantly faster than using pow
#define CUBE(x) ( (x)*(x)*(x) )
#define POW4(x) ( (x)*(x)*(x)*(x) )
// Scalar memory locations for dust grain parameters
#define CLOC 0  // Wind "colour"
#define ZLOC 1  // Dust-to-gas mass ratio
#define ALOC 2  // Grain radius
// Logical operators, enabled with pgen file
bool cool, dust;
bool dust_cool;
// Dust properties
Real initialDustToGasMassRatio;
Real initialGrainRadiusMicrons;
// Stellar properties
Real mass1, mass2;  // Mass, msol
Real mdot1, mdot2;  // Mass loss rates, msol yr^-1
Real vinf1, vinf2;  // Terminal veloicites, cm s^-1

// Array for storing mass fractions, first index is for star, 0 for WR 1 for OB
// Second index is for specific element in the form:
// 0: Hydrogen
// 1: Helium
// 2: Carbon
// 3: Nitrogen
// 4: Oxygen
// This assumes that 
Real massFrac[2][5];     // Mass fractions

Real remapRadius1, remapRadius2;  // Remap radii (cm)
// Orbital positions
Real xpos1, xpos2, ypos1, ypos2, zpos1, zpos2;  // Star positions from barycenter
Real xvel1, xvel2, yvel1, yvel2, zvel1, zvel2;  // Star velocities, cm s^-1
// Other orbital properties
Real orbitPhase;  // Orbital phase (t/P)
Real dsep;        // Star separation distance (cm)
Real rWR, rOB;    // Stagnation point distance (cm)
Real xiWR, xiOB;  // Cooling parameters (dimensionless)
Real period;      // orbit period (s)  
Real phaseoff;    // phase offset of orbit (from periastron) 
Real ecc;         // orbit eccentricity
// Thermodynamic properties
Real tmin,tmax;   // min/max temperature allowed on grid
// Refinement properties
int G0_level = 0;
// Physical and mathematical constants
const Real pi       = 2.0*asin(1.0);  // Take a wild guess
const Real Msol     = 1.9891e33;      // Solar mass to grams conversion
const Real yr       = 3.15569e7;      // Seconds to year converson
const Real boltzman = 1.380658e-16;   // Boltzman constant (CGS)
const Real massh    = 1.6726219e-24;  // Proton/hydrogen mass, g
const Real masse    = 9.1093837e-28;  // Electron mass, g
// User defined constants
const Real minimumDustToGasMassRatio = 1.0e-7;
const Real minimumGrainRadiusMicrons = 0.0001;
const Real grainBulkDensity = 3.0;                       // (g/cm^3)
const Real Twind = 1.0e4;                                // K
const Real avgmass = 1.0e-24;     // g
// Structure to store ga
// Structure to store gas cooling curve data
struct coolingCurve{
  int ntmax;  // Number of bins in cooling curve 
  std::string coolCurveFile;  // Filename
  std::vector<double> logt,lambdac,te,loglambda;
  double t_min,t_max,logtmin,logtmax,dlogt;
};

//! \class FreeElectronCurve

class FreeElectronCurve {
  public:
    // Variables
    std::vector<Real> T;
    std::vector<Real> e;
    Real TMin,TMax;
    Real eMin,eMax;
    int Init(std::string freeElectronFilename);
};

int FreeElectronCurve::Init(std::string freeElectronFilename) {
  std::ifstream file(freeElectronFilename);
  if (!file) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.cpp freeElectronCurve::Init" <<
           std::endl <<
           "File " <<
           freeElectronFilename <<
           " Cannot be read!" <<
           std::endl;
    ATHENA_ERROR(msg);
    return 1;
  }
  else {
    Real tbuf,ebuf;
    std::string line;
    file.seekg(0);
    while (true) {
      file >> tbuf >> ebuf;
      if(file.eof()) break;
      T.push_back(tbuf);
      e.push_back(ebuf);
    }
  }
  // Convert log(T) to T
  for (int n = 0; n < T.size(); n++) {
    T[n] = pow(10.0,T[n]);
  }
  // Compare lengths of vectors, if they are not the same, import has failed
  if (T.size() != e.size()) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.cpp FreeElectronCurve::Init" <<
           std::endl <<
           "Import has failed, mismatch in temperature and free elctron vector size" <<
           std::endl;
    ATHENA_ERROR(msg);
  }
  TMin = T.front();
  TMax = T.back();
  eMin = e.front();
  eMax = e.back();
  if (TMax < TMin) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.cpp FreeElectronCurve::Init" <<
           std::endl <<
           "Import has failed, Data does not appear to be sorted!" <<
           std::endl;
    ATHENA_ERROR(msg);
  }
  return 0;
}

FreeElectronCurve eCurve[2];

//! \class CoolCurve
//  \brief A self-initialising class containing cooling curves used in CWB problem
//         with cooling enabled.
//  - Cooling is stored as a b
//  - Binary search and linear interpolation is used to derive cooling parameter Lambda
//    which is then used to calculate energy loss using the formulae
//      dE/dt = nG^2 * Lambda(T)
//  - Where nG is the gas number density.

class CoolCurve {
  public:
    // Variables
    std::vector<Real> temp;    // Temperature (K)
    std::vector<Real> lambda;  // Energy loss constant
    Real tmin, tmax;  // Values of first and last index in temperature
    Real lmin, lmax;  // Values of first and last index in lambda
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
    msg << "### FATAL ERROR in cwb.cpp CoolCurve::Init" <<
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
    while (true) {
      file >> tbuf >> lbuf;
      if(file.eof()) break;
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
    msg << "### FATAL ERROR in cwb.cpp CoolCurve::Init" <<
           std::endl <<
           "Import has failed, mismatch in temperature and lambda vector size" <<
           std::endl;
    ATHENA_ERROR(msg);
  }
  tmin = temp.front();
  tmax = temp.back();
  lmin = lambda.front();
  lmax = lambda.back();
  if (tmax < tmin) {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.cpp CoolCurve::Init" <<
           std::endl <<
           "Import has failed, Data does not appear to be sorted!" <<
           std::endl;
    ATHENA_ERROR(msg);
  }
  return 0;
}

//! \fn FindLambda
// Perform a binary search to find first temperature in cooling curve > input temperature, tu
// Function then performs a linear interpolation between tu and the previous temperature in index, tl, to find the interpolated value for lambda
// Fast, but not particularly accurate at lower temperatures
// Uses C++ algorithm library
Real CoolCurve::FindLambda(Real t) {
  // Quickly return values that would cause search to fail
  if      (t <= 1.1*tmin) {return 0.0;}
  else if (t >= tmax) {return lmax;}
  // Perform binary search
  auto upper = std::upper_bound(temp.begin(),temp.end(),t);
  // Assign indexes based on search
  int iu = std::distance(temp.begin(), upper);
  int il = iu - 1;
  // Perform linear interpolation
  // Assign values from each array from indexes found by binary search
  Real tl = temp[il];
  Real tu = temp[iu];
  Real ll = lambda[il];
  Real lu = lambda[iu];
  // Calculate interpolated value of lambda
  Real l = ll + ((t - tl) * ((lu - ll) / (tu - tl)));
  return l;
}

// Define two cooling curves, one for each star
CoolCurve cc[2];

// JWE prototypes
Real Calc_h_e(Real x_e);  // Calculate electron heating efficiency
Real SearchAndInterpolate(Real x, std::vector<Real> xarr, std::vector<Real> yarr);
// JMP prototypes
void AdjustPressureDueToCooling(int is,int ie,int js,int je,int ks,int ke,Real gmma1,AthenaArray<Real> &dei,AthenaArray<Real> &cons);
// History user functions
Real UserHistoryFunction(MeshBlock *pmb, int iout);
Real DustCreationRateInWCR(MeshBlock *pmb, int iout);
Real ReturnOrbitalProperties(MeshBlock *pmb, int iout);
Real ReturnStag(MeshBlock *pmb, int iout);
Real ReturnXi(MeshBlock *pmb, int iout);
void EvolveDust(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
//void FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void OrbitCalc(Real t);
void PhysicalSources(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
void RadiateHeatCoolOld(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
int RefinementCondition(MeshBlock *pmb);
void RestrictCool(int is,int ie,int js,int je,int ks,int ke,int nd,Real gmma1,AthenaArray<Real> &dei,const AthenaArray<Real> &cons);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor. Also called when restarting.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // First, read in star properties
  // Read in mass loss rates and terminal velocities
  mdot1 = pin->GetReal("problem","mdot1");
  mdot2 = pin->GetReal("problem","mdot2");
  vinf1 = pin->GetReal("problem","vinf1");
  vinf2 = pin->GetReal("problem","vinf2");
  // Read in WR wind mass fractions
  massFrac[0][0] = pin->GetReal("problem","xH1");
  massFrac[0][1] = pin->GetReal("problem","xHe1");
  massFrac[0][2] = pin->GetReal("problem","xC1");
  massFrac[0][3] = pin->GetReal("problem","xN1");
  massFrac[0][4] = pin->GetReal("problem","xO1");
  // Read in OB wind mass fractions
  massFrac[1][0] = pin->GetReal("problem","xH2");
  massFrac[1][1] = pin->GetReal("problem","xHe2");
  massFrac[1][2] = pin->GetReal("problem","xC2");
  massFrac[1][3] = pin->GetReal("problem","xN2");
  massFrac[1][4] = pin->GetReal("problem","xO2");
  // Remap radius, typically set to 5-10x finest cell radius
  remapRadius1 = pin->GetReal("problem","remapRadius1");
  remapRadius2 = pin->GetReal("problem","remapRadius2");
  // Read in star positions
  // TODO: Get rid of these, as they aren't really needed!
  xpos1 = pin->GetReal("problem","xpos1");
  ypos1 = pin->GetReal("problem","ypos1");
  zpos1 = pin->GetReal("problem","zpos1");
  xpos2 = pin->GetReal("problem","xpos2");
  ypos2 = pin->GetReal("problem","ypos2");
  zpos2 = pin->GetReal("problem","zpos2");
  // Read in other orbital deteails
  mass1 = pin->GetReal("problem","mass1");
  mass2 = pin->GetReal("problem","mass2");
  ecc   = pin->GetReal("problem","ecc");
  period = pin->GetReal("problem","period");
  phaseoff = pin->GetReal("problem","phaseoff");
  // Determine if dust is being modelled
  std::string dusty = pin->GetString("problem","dust");
  if      (dusty == "on")  dust = true;
  else if (dusty == "off") dust = false;
  else{
    std::cout << "dust value not recognized: " << dust << "; Aborting!\n";
    exit(EXIT_SUCCESS);
  }
  // Determine if cooling is present
  std::string cooling = pin->GetString("problem","cooling");
  if      (cooling == "on")  cool = true;
  else if (cooling == "off") cool = false;
  else{
    std::cout << "cooling value not recognized: " << cooling << "; Aborting!\n";
    exit(EXIT_SUCCESS);
  }
  // Determine if dust cooling is present
  if (dust == true) {
    std::string dust_cooling = pin->GetString("problem","dust_cooling");
    if      (dust_cooling == "on")  dust_cool = true;
    else if (dust_cooling == "off") dust_cool = false;
    else{
      std::cout << "dust_cooling value not recognized: " << dust_cooling << "; Aborting!\n";
      exit(EXIT_SUCCESS);
    }
  }

  if (dust && NSCALARS < 3){
    // Scalars are: 0 = wind colour
    //              1 = dust to gas mass ratio
    //              2 = dust grain radius (microns)
    // Can also be called with:
    //           CLOC = Wind colour
    //           ZLOC = dust to gas mass ratio
    //           ALOC = dust grain radius (microns)
    std::cout << "Not enough scalars for dust modelling. NSCALARS = " << NSCALARS << ". Aborting!\n";
    exit(EXIT_SUCCESS);
  }
  if (dust){
    initialDustToGasMassRatio = pin->GetReal("problem","initialDustToGasMassRatio");
    initialGrainRadiusMicrons = pin->GetReal("problem","initialGrainRadiusMicrons");
  }

  // Convert values into CGS
  mdot1 *= Msol/yr; 
  mdot2 *= Msol/yr;

  // Read in cooling curves
  if (cool) {
    std::string WRCoolCurveFileName = pin->GetString("problem","WRCoolCurve");
    std::string OBCoolCurveFileName = pin->GetString("problem","OBCoolCurve");
    cc[0].Init(WRCoolCurveFileName);
    cc[1].Init(OBCoolCurveFileName);
  }

  if (dust_cool) {
    eCurve[0].Init("electron_WC");
    eCurve[1].Init("electron_solar");
  }

  // Note: it is only possible to have one source functions enrolled by the user.
  EnrollUserExplicitSourceFunction(PhysicalSources);

  if (adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);  

  // Add a user-defined global output (see https://github.com/PrincetonUniversity/athena-public-version/wiki/Outputs)
  if (dust){
    // Because of the way user history outputs work, this function needs
    // to be run 5 times across entire numerical grid, which isn't
    // particularly efficient, but still not too bad compared to entire
    // hydro grid
    AllocateUserHistoryOutput(15);
    // Orbital properties
    EnrollUserHistoryOutput(0, UserHistoryFunction, "phase",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(1, UserHistoryFunction, "x1",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(2, UserHistoryFunction, "y1",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(3, UserHistoryFunction, "x2",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(4, UserHistoryFunction, "y2",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(5, UserHistoryFunction, "dsep",
                            UserHistoryOperation::maxpm);
    // Stagnation points
    EnrollUserHistoryOutput(6, UserHistoryFunction, "rWR",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(7, UserHistoryFunction, "rOB",
                            UserHistoryOperation::maxpm);
    // Cooling parameters
    EnrollUserHistoryOutput(8, UserHistoryFunction, "xiWR",
                            UserHistoryOperation::maxpm);
    EnrollUserHistoryOutput(9, ReturnXi, "xiOB",
                            UserHistoryOperation::maxpm);
    // Basic dust parameters
    EnrollUserHistoryOutput(10, UserHistoryFunction, "dmdustdt_WCR");
    EnrollUserHistoryOutput(11, UserHistoryFunction, "dmdust_WCR_dt_created");
    EnrollUserHistoryOutput(12, UserHistoryFunction, "dmdust_WCR_dt_lost");
    EnrollUserHistoryOutput(13, UserHistoryFunction, "dust_WCR");
    EnrollUserHistoryOutput(14, UserHistoryFunction, "dust_TOTAL");
  }

  // No mesh blocks exist at this point...
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the CWB test. 
//         This is not called during a restart, so do not put any important stuff in here.
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  //std::cout << "[MeshBlock::ProblemGenerator]\n";

  Real gmma  = peos->GetGamma();
  Real gmma1 = gmma - 1.0;
  
  Real dsep = std::sqrt(std::pow(xpos1 - xpos2,2) + std::pow(ypos1 - ypos2,2) + std::pow(zpos1 - zpos2,2));
  Real eta = mdot2*vinf2/(mdot1*vinf1);                   // wind mtm ratio
  Real rob = (std::sqrt(eta)/(1.0 + std::sqrt(eta)))*dsep;//distance of stagnation point from star 1 (distance from star 0 is rwr)

  OrbitCalc(pmy_mesh->time);
  
  // Map on wind 1
  Real xmaxst1 = xpos1 + dsep - rob; // maximum x-extent of wind1
  for (int k=ks; k<=ke; k++) {
    Real zc = pcoord->x3v(k) - zpos1;
    Real zc2 = zc*zc;
    for (int j=js; j<=je; j++) {
      Real yc = pcoord->x2v(j) - ypos1;
      Real yc2 = yc*yc;
      for (int i=is; i<=ie; i++) {
        Real xc = pcoord->x1v(i) - xpos1;
	Real xc2 = xc*xc;
	Real r2 = xc2 + yc2 + zc2;
	Real r = std::sqrt(r2);
	Real xy = std::sqrt(xc2 + yc2);
	Real sinphi = xy/r;
	Real cosphi = zc/r;
        Real costhta = xc/xy;
	Real sinthta = yc/xy;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          if (pcoord->x1v(i) < xmaxst1) {
	    Real rho = mdot1/(4.0*pi*r2*vinf1);
	    Real u1 = vinf1*sinphi*costhta + xvel1;
	    Real u2 = vinf1*sinphi*sinthta + yvel1;
	    Real u3 = vinf1*cosphi;
	    Real pre = (rho/avgmass)*boltzman*Twind;
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	    if (NSCALARS > 0) {
              pscalars->s(0,k,j,i) = rho;
            }
	    if (dust){
              pscalars->s(1,k,j,i) = initialDustToGasMassRatio*rho;
              pscalars->s(2,k,j,i) = initialGrainRadiusMicrons*rho;
	    }
	  }
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) { // RPZ
          if (pcoord->x3v(k) < xmaxst1) { // x = R; y = P; z = Z
	    Real rho = mdot1/(4.0*pi*r2*vinf1);
	    Real u1 = vinf1*sinphi;
	    Real u2 = 0.0;
	    Real u3 = vinf1*cosphi;
	    Real pre = (rho/avgmass)*boltzman*Twind;
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*vinf1*vinf1;
	    if (NSCALARS > 0) {
              pscalars->s(0,k,j,i) = rho;
            }
	    if (dust){
              pscalars->s(1,k,j,i) = initialDustToGasMassRatio*rho;
              pscalars->s(2,k,j,i) = initialGrainRadiusMicrons*rho;
	    }
	  }	  
        }
      }
    }
  }

  // Map on wind 2
  for (int k=ks; k<=ke; k++) {
    Real zc = pcoord->x3v(k) - zpos2;
    Real zc2 = zc*zc;
    for (int j=js; j<=je; j++) {
      Real yc = pcoord->x2v(j) - ypos2;
      Real yc2 = yc*yc;
      for (int i=is; i<=ie; i++) {
        Real xc = pcoord->x1v(i) - xpos2;
	Real xc2 = xc*xc;
	Real r2 = xc2 + yc2 + zc2;
	Real r = std::sqrt(r2);
	Real xy = std::sqrt(xc2 + yc2);
	Real sinphi = xy/r;
	Real cosphi = zc/r;
        Real costhta = xc/xy;
	Real sinthta = yc/xy;
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          if (pcoord->x1v(i) >= xmaxst1) {
	    Real rho = mdot2/(4.0*pi*r2*vinf2);
	    Real u1 = vinf2*sinphi*costhta + xvel2;
	    Real u2 = vinf2*sinphi*sinthta + yvel2;
	    Real u3 = vinf2*cosphi;
	    Real pre = (rho/avgmass)*boltzman*Twind;
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	    if (NSCALARS > 0) {
              pscalars->s(0,k,j,i) = 0.0;
            }
	    if (dust){
              pscalars->s(1,k,j,i) = initialDustToGasMassRatio*rho;
              pscalars->s(2,k,j,i) = initialGrainRadiusMicrons*rho;
	    }
	  }
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) { // RPZ
          if (pcoord->x3v(k) >= xmaxst1) { // x = R; y = P; z = Z
	    Real rho = mdot2/(4.0*pi*r2*vinf2);
	    Real u1 = vinf2*sinphi;
	    Real u2 = 0.0;
	    Real u3 = vinf2*cosphi;
	    Real pre = (rho/avgmass)*boltzman*Twind;
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*vinf2*vinf2;	    
	    if (NSCALARS > 0) {
              pscalars->s(0,k,j,i) = 0.0;
            }
	    if (dust){
              pscalars->s(1,k,j,i) = initialDustToGasMassRatio*rho;
              pscalars->s(2,k,j,i) = initialGrainRadiusMicrons*rho;
	    }
	  }	  
        }
	

      }
    }
  }
  

  return;
}


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
// JMP: This is called after the end of all of the hydro steps (i.e. just before the simulation exits)
//========================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // Cycle over all MeshBlocks
  //MeshBlock *pmb = pblock;
  //while (pmb != nullptr) {
  //  pmb = pmb->next;
  //}
  return;
}


//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================
void MeshBlock::UserWorkInLoop() {

  Real gmma  = peos->GetGamma();
  Real gmma1 = gmma - 1.0;

  // Calculate the stellar positions
  OrbitCalc(pmy_mesh->time); 

  Real xpos[2]={xpos1,xpos2};
  Real ypos[2]={ypos1,ypos2};
  Real zpos[2]={zpos1,zpos2};
  Real xvel[2]={xvel1,xvel2};
  Real yvel[2]={yvel1,yvel2};
  Real mdot[2]={mdot1,mdot2};
  Real vinf[2]={vinf1,vinf2};
  Real scalar[2]={1.0,0.0};
  Real remapRadius[2] = {remapRadius1,remapRadius2};

  Real dx = pcoord->dx1v(0);
  int remapi[2] = {int(remapRadius1/dx),int(remapRadius2/dx)};
    
  // Remap winds
  for (int nw = 0; nw < 2; ++nw){ // Loop over each wind
    int istar = int((xpos[nw] - pcoord->x1f(0))/pcoord->dx1f(0));
    int jstar = int((ypos[nw] - pcoord->x2f(0))/pcoord->dx2f(0));
    int kstar = int((zpos[nw] - pcoord->x3f(0))/pcoord->dx3f(0));
    int istl = std::max(is,istar-remapi[nw]-2);
    int jstl = std::max(js,jstar-remapi[nw]-2);
    int kstl = std::max(ks,kstar-remapi[nw]-2);
    int istu = std::min(ie,istar+remapi[nw]+2);
    int jstu = std::min(je,jstar+remapi[nw]+2);
    int kstu = std::min(ke,kstar+remapi[nw]+2);
    if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      jstl = js; jstu = je;
    }
    for (int k=kstl; k<=kstu; k++) {
      Real zc = pcoord->x3v(k) - zpos[nw];
      Real zc2 = zc*zc;
      for (int j=jstl; j<=jstu; j++) {
        Real yc = pcoord->x2v(j) - ypos[nw];
        Real yc2 = yc*yc;
        for (int i=istl; i<=istu; i++) {
          Real xc = pcoord->x1v(i) - xpos[nw];
          Real xc2 = xc*xc;
          Real r2 = xc2 + yc2 + zc2;
          Real r = std::sqrt(r2);
          if (r < remapRadius[nw]) {
            Real xy = xy = std::sqrt(xc2 + yc2);
            Real sinphi = xy/r;
            Real cosphi = zc/r;
            Real costhta = xc/xy;
            Real sinthta = yc/xy;
            Real rho = mdot[nw]/(4.0*pi*r2*vinf[nw]);
            Real u1, u2, u3;
            Real pre = (rho/avgmass)*boltzman*Twind;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
              u1 = vinf[nw]*sinphi*costhta + xvel[nw];
              u2 = vinf[nw]*sinphi*sinthta + yvel[nw];
              u3 = vinf[nw]*cosphi;
            }
            else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) { // RPZ
              u1 = vinf1*sinphi;
              u2 = 0.0;
              u3 = vinf1*cosphi;
            }
            Real v2 = SQR(u1) + SQR(u2) + SQR(u3);
            phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*v2;
            // Set passive scalars
            if (NSCALARS > 0) {
              // wind "colour"
              pscalars->s(0,k,j,i) = scalar[nw]*rho;
              if (dust){
                pscalars->s(1,k,j,i) = initialDustToGasMassRatio*rho;
                pscalars->s(2,k,j,i) = initialGrainRadiusMicrons*rho;
              }
            }
          }
	      }
      }
    }
  }


  // Limit wind colour
  if (NSCALARS > 0) {  
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
	  // wind "colour"
          pscalars->s(0,k,j,i) = std::min(std::max(pscalars->s(0,k,j,i),0.0),1.0);
        }
      }
    }
  }


  
  return;
}



//========================================================================================
//! \fn void PhysicalSources(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Physical source terms
//========================================================================================
void PhysicalSources(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){

  if (cool) RadiateHeatCoolOld(pmb,dt,cons);
  if (dust) EvolveDust(pmb,dt,cons);
  return;
}

//! \fn Real SearchAndInterpolate(Real xarray[], Real yarray[], Real x)
// Function then performs a linear interpolation between two values in a lookup table
// Uses C++ algorithm library
Real SearchAndInterpolate(Real x, std::vector<Real> xarr, std::vector<Real> yarr) {
  // Perform binary search using <algorithm> routine
  auto upper = std::upper_bound(xarr.begin(),xarr.end(),x);
  // Assign indexes to upper bound, found by <algorithm> and lower bound
  int iu = std::distance(xarr.begin(),upper);
  int il = iu - 1;
  // Perform linear interpolation
  // Find upper and lower values for each axis
  Real xl = xarr[il];
  Real xu = xarr[iu];
  Real yl = yarr[il];
  Real yu = yarr[iu];
  // Interpolate
  Real y = yl + ((x - xl) * ((yu - yl) / (xu - xl)));
  // Finish up and return
  return y;
}

//! \fn Real Calc_h_e(Real x_e)
//  \brief Approximation of electron heating efficiency
//         - Electron heating is based on on temperature and grain radius
//         - This function uses a power law series to simplify a complex 
//           integral in order to calculate lambda due to dust on a
//           cell-by-cell basis
//         - Previously the integral form was used, but this massively reduced
//           performance, even with a simple trapezoidal rule integration
//         Derived from:
//         Dwek, E., & Werner, M. W. (1981).
//         The Infrared Emission From Supernova Condensates.
//         The Astrophysical Journal, 248, 138.
//         https://doi.org/10.1086/159138

Real Calc_h_e(Real x_e) {
  Real h_e;
       if (x_e > 4.5) h_e = 1;
  else if (x_e > 1.5) h_e = 0.37 * pow(x_e,0.62);
  else                h_e = 0.27 * pow(x_e,1.50);
  return h_e;
}

//! \fn Real CalcLambdaDust(Real nH, Real a, Real T)
//  \brief Calculate energy loss per dust grain, multiply by nD to calculate
//         cell cooling rate
//         - Energy lost from the gas flow due to dust is mainly due to
//           collisional heating of the dust particles from atoms and
//           electrons
//         - Efficiency losses can occur at high temperatures as particles
//           are so energetic they pass through one another
//         - This function approximates this effect
//         - Resultant value for single grain, to find energy loss in erg/s/cm^3
//           value must be multiplied by nD
//         Derived from:
//         Dwek, E., & Werner, M. W. (1981).
//         The Infrared Emission From Supernova Condensates.
//         The Astrophysical Journal, 248, 138.
//         https://doi.org/10.1086/159138

Real CalcGrainCoolRate(Real rhoG, Real a, Real T, int nc) {
  // === Setup ===
  // Shorten temp * boltzmann constant to kBT, as this is used a lot
  Real kBT = boltzman * T;
  // - Allocate a counter for the total amount of energy from each
  //   element collision
  // - Precalculated arrays for critical energy constant and atomic mass
  //   in CGS units for each type of element are declared:
  //   - For H atoms:   Ec = 133keV, m = 1.00797 AMU
  //   - For He atoms:  Ec = 222keV, m = 4.00260 AMU
  //   - For CNO atoms: Ec = 665keV, m = 12.011,14.0067,15.9994 AMU
  // - This assumes that carbon is dominant in CNO bin
  Real E_E[5] = {2.1308949e-07,3.5568321e-07,1.0654475e-06,1.0654475e-06,1.0654475e-06};
  Real m_E[5] = {1.6737736e-24,6.6464737e-24,1.9944735e-23,2.3258673e-23,2.6567629e-23};
  Real n_E[5] = {0.0,0.0,0.0,0.0,0.0};  // Elemental number density, cm^-3
  Real H_E[5] = {0.0,0.0,0.0,0.0,0.0};  // Heating for each element, erg s^-1
  // === Processing ===
  // Determine heating rate from individual grain due to colliding atoms
  Real H_coll = 0.0;  // Heating rate, erg s^-1
  Real nT     = 0.0;  // Total number density, cm^-3
  for (int n = 0 ; n < 5 ; n++) {
    // Calculate number density for element
    n_E[n] = rhoG * massFrac[nc][n] / m_E[n];
    // Calculate total number density
    nT += n_E[n];
    // Calculate the critical energy of incident hydrogen atoms
    Real EC = E_E[n] * a;  
    // Calculate h_n (Grain heating efficiency due to atoms)
    Real h_n = 1.0 - (1.0 + EC / (2.0 * kBT)) * exp(-EC/kBT);
    // Calculate heating rate of element, Eq 2 of DW81
    H_E[n]  = 1.26e-19 * SQR(a) * pow(T,1.5) * n_E[n] * h_n;
    H_E[n] /= sqrt(m_E[n]/massh);
    H_coll += H_E[n];  // Add to counter
  }
  // Determine heating rate due to colliding electrons
  // Approximation from Eq A12 of DW81 used for calculating electron-grain
  // "transparency", this can be off by up to 1%, but is nearly 2,000x
  // faster than a trap rule integration with 400 components and
  // uses almost no memory
  // Estimate the total number of free electrons per ion
  Real nFreeElectrons = 0.0;
       if (T < eCurve[nc].TMin) nFreeElectrons = std::min(1.0,eCurve[nc].eMin);
  else if (T > eCurve[nc].TMax) nFreeElectrons = eCurve[nc].eMax;
  else nFreeElectrons = SearchAndInterpolate(T,eCurve[nc].T,eCurve[nc].e);
  // Using this estimated value, calculate the electron number density
  Real ne = nT * nFreeElectrons;
  // Calculate the critical energy for electron to penetrate grain
  // This makes the assumption of an uncharged dust grain, Ee = Ec
  Real Ee  = 3.6850063e-08 * pow(a,2.0/3.0);  // DW81 Eq A6
  Real x_e = Ee/kBT;                          // DW81 Eq A11
  // Approximate electron-grain "transparency"
  Real h_e = Calc_h_e(x_e);                   // DW81 Eq A12
  // Calculate heating rate of electrons, using DW81 Eq 2
  Real H_el  = 1.26e-19 * SQR(a) * pow(T,1.5) * ne * h_e;
       H_el /= sqrt(masse/massh);
  // Summate heating rates and normalise by nd*np
  Real edotGrain = (H_coll + H_el);
  // Finish up and return!
  return edotGrain;
}

//! \fn void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons)
//  \brief Reduce energy in simulation to simulate cooling through gas and dust emission
//         - Gas cooling is handled through a cooling curve
//           - Cooling curve stored as lookup table, binary search and interpolation
//             are performed in order to improve accuracy
//         - Dust cooling is estimated through Dwek Werner prescription, see function
//           CalcGrainCoolRate()
//         - Simple checks are performed to ensure that results are sensible
//         - Adaptive cooling look is used to determine the change in total temperature
//         - This is then used to calculate average energy loss over hydro timestep
//         - This is then integrated and used to recalculate conserved arrays, modifying
//           the total cell energy

void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons) {
  const int ncool = 2;  // Number of cooling curves
  // Calculate gamma, shorten gamma-1 to a constant
  const Real gamma = pmb->peos->GetGamma();
  const Real g1    = gamma - 1.0;
  const Real tmin  = Twind;
  const Real tmax  = 3.0e9;

  // Build an array to store 
  AthenaArray<Real> dEintArray(pmb->ke+1,pmb->je+1,pmb->ie+1);
  // Iterate through meshblocks array, calculating cooling rates
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // Import conserved variables
        Real rho  = cons(IDN,k,j,i);        // Gas density (g cm^-3)
        Real u1   = cons(IM1,k,j,i) / rho;  // x velocity (cm s^-1)
        Real u2   = cons(IM2,k,j,i) / rho;  // y velocity (cm s^-1)
        Real u3   = cons(IM3,k,j,i) / rho;  // z velocity (cm s^-1)
        Real Et   = cons(IEN,k,j,i);        // Total energy (erg)
        // Import wind "colour"
        Real col  = pmb->pscalars->s(0,k,j,i) / rho;
        // Calculate gas parameters from conserved variables
        Real v    = std::sqrt(SQR(u1)+SQR(u2)+SQR(u3));  // Gas velocity (cm s^-1)
        Real kE   = 0.5 * rho * SQR(v);                  // Kinetic energy (erg)
        Real ie   = Et - kE;                             // Internal energy (erg)
        Real P    = g1 * ie;                             // Gas pressure (Ba)
        Real Ti   = P * avgmass / (rho*boltzman);      // Initial temperature (K)
        // Convert wind "colour" to an array form
        Real windCol[2] = {0.0,0.0};  // Initialise array
        windCol[0] = col;             // WR wind mass fraction
        windCol[1] = 1.0 - col;       // OB wind mass fraction
        // Set up cooling loop variables
        Real T      = Ti;             // Current temperature (K)
        Real rhomH  = rho / massh;
        Real rhomH2 = rhomH * rhomH;  // Plasma cooling curve constant (cm^-6)
        // Set up dust cooling variables
        Real a;     // Grain radius (micron)
        Real aCGS;  // Grain radius (cm)
        Real z;     // Dust-to-gas mass ratio
        Real nT;    // Total number density (cm^-3)
        Real nD;    // Dust number density (cm^-3)
        // Calculate dust cooling variables
        if (dust_cool) {
          z = pmb->pscalars->s(ZLOC,k,j,i) / rho;
          a = pmb->pscalars->s(ZLOC,k,j,i) / rho;
          // Convert a in microns to CGS units, cm
          aCGS = a * 1e-4;
          // Calculate dust grain number density 
          Real rhoD  = rho * z;
          Real volD  = (4.0/3.0) * PI * CUBE(aCGS);
          Real massD = grainBulkDensity * volD;
          // Finish nD calculation
          nD = rhoD / massD;
        }
        // If the temperature is close to the wind temperature
        // cool to the wind temperature to save time by bypassing cooling loop
        if (T < 1.5 * Twind) {
          T = Twind;
        }
        // Otherwise, begin the cooling loop
        else {
          // Initialise cooling variables
          Real dtInt = 0.0;  // Amount of timestep integrated (s)
          // Perform cooling step in loop
          while (dtInt < dt) {            
            // Initialise loop cooling arrays
            Real lambdaWind[2] = {0.0,0.0};  // Gas cooling parameter
            Real edotGas[2]    = {0.0,0.0};  // Gas cooling energy loss rate (erg/s)
            Real edotDust[2]   = {0.0,0.0};  // Dust cooling energy loss rate (erg/s)
            Real edotWind[2]   = {0.0,0.0};  // Total wind energy loss rate (erg/s)
            Real totalCool     = 0.0;        // Total cooling rate (erg/s)
            for (int nc = 0; nc < 2; ++nc) {
              // Calculate cooling contribution to plasma
              lambdaWind[nc] = cc[nc].FindLambda(T);
              edotGas[nc]    = rhomH2 * lambdaWind[nc];
              if (dust_cool) {
                // Calculate cooling contribution to 
                Real edotGrain = CalcGrainCoolRate(rho,a,T,nc);
                edotDust[nc]   = nD * edotGrain;
              }
              // Add individual contributions to total array
              edotWind[nc] += edotGas[nc];
              edotWind[nc] += edotDust[nc];
              // Multiply by wind colour to find wind contribution in cell
              edotWind[nc] *= windCol[nc];
              // Add to total energy lost in loop iteration
              totalCool += edotWind[nc];
            }

            if (totalCool > 0) std::cout << totalCool << " " << rhomH << " " << T <<  " " << lambdaWind[0] << "\n";

            // Calculate timestep interval
            Real Eint   = P / g1;            // Current internal energy
            Real tCool  = Eint / totalCool;  // Cooling time (s)
            Real dtCool = 0.1 * abs(tCool);  // Timestep, based on cooling time
            if ((dtInt + dtCool) > dt) {
              dtCool = dt - dtInt;  // 
            }
            dtInt += dtCool;  // Update total integrated time
            // Calculate new temperature and update pressure
            Real dEint = -totalCool * dtCool;
            Real TNew  = T * (Eint + dEint) / Eint;
            // Check to see if new temperature has reached lower or upper bound
            TNew = std::max(tmin,TNew);
            TNew = std::min(TNew,tmax);
            // Update pressure and temperature
            P *= TNew / T;
            T  = TNew;
          }
        }
        // Determine if cooling is sensible
        // If a solution cannot be found, abandon cooling!
        // Rejects NaN temperatures, infinite temperatures, or negative temperatures
        if (std::isnan(T) || std::isinf(T) || T < 0.0) {
          T = Ti;
        }
        // Calculate average rate of change in energy
        // by calculating 
        Real dEintAvg = (Ti - T) * rho;
        dEintArray(k,j,i) = dEintAvg;
      }
    }
  }
  // Restrict cooling at unresolve interfaces
  RestrictCool(pmb->is,pmb->ie,
               pmb->js,pmb->je,
               pmb->ks,pmb->ke,
               pmb->pmy_mesh->ndim,
               g1,
               dEintArray,
               cons);
  // Adjust pressure, performing the actual cooling
  AdjustPressureDueToCooling(pmb->is,pmb->ie,
                             pmb->js,pmb->je,
                             pmb->ks,pmb->ke,
                             g1,
                             dEintArray,
                             cons);
  // Perform garbage collection and finish
  dEintArray.DeleteAthenaArray();

  std::cout << "Completed another meshblock!\n";

  return;
}

//========================================================================================
//! \fn void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
//  \brief Calculating radiative heating and cooling
//========================================================================================
void RadiateHeatCoolOld(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons){

  const int ncool = 2; // number of cooling curves

  static coolingCurve cc[ncool];
  static bool firstHeatCool = true;
  bool restrictUnresolvedCooling = true;

  Real gmma  = pmb->peos->GetGamma();
  Real gmma1 = gmma - 1.0;

  tmin = Twind;
  tmax = 3.0e9;

  //const double a = grainRadius;                      // grain radius (cm)
  //const double dens_g = grainBulkDensity;            // (g/cm^3)
  //const double massD = (4.0/3.0)*pi*pow(a,3)*dens_g; // grain mass (g)

  double totEdotGas = 0.0;
  double totEdotDust = 0.0;

  if (firstHeatCool){
    // Specify the cooling curves and the number of temperature bins.
    // The first cooling curve should extend to the same or lower temperature
    // than the second.
    // The cooling_KI02_4.0_CLOUDY_7.6_MEKAL.txt gas cooling rate is lambda*(dens/mH)^2.
    // The dust_lambdaD_solar_a0.1.txt dust cooling rate is lambda_D*np*ne.
    cc[0].coolCurveFile = "cooling_curve_WC";
    cc[1].coolCurveFile = "cooling_curve_solar";
    cc[0].ntmax = 101;
    cc[1].ntmax = 101;
    //#ifdef DUST
    //    cc[1].coolCurveFile = "dust_lambdaD_solar_a0.1.txt";
    //    cc[1].ntmax = 101;
    //#endif

    // Read in each cooling curve. It is assumed that the temperature binning is
    // uniform in log space i.e. dlog(T) is a constant
    for (int nc = 0; nc < ncool; ++nc)
    {
      std::ifstream file(cc[nc].coolCurveFile.c_str());
      if (!file)
      {
        std::cerr << "radiateHeatCool: failed to open file " << cc[nc].coolCurveFile << "\n";
        exit(EXIT_FAILURE);
      }
      for (int n = 0; n < cc[nc].ntmax; ++n)
      {
        double logt, lambda;
        file >> logt >> lambda;
        cc[nc].logt.push_back(logt);
        cc[nc].lambdac.push_back(lambda);
        cc[nc].te.push_back(pow(10, logt));
        cc[nc].loglambda.push_back(log10(lambda));
        //if (Globals::my_rank == 0) std::cout << "nc = " << nc << "; n = " << n << "\t" << cc[nc].te[n] << "\t" << cc[nc].lambdac[n] << "\n";
      }
      file.close();

      cc[nc].t_min = cc[nc].te[0];                // K
      cc[nc].t_max = cc[nc].te[cc[nc].ntmax - 1]; // K
      cc[nc].logtmax = cc[nc].logt[cc[nc].ntmax - 1];
      cc[nc].logtmin = cc[nc].logt[0];
      cc[nc].dlogt = (cc[nc].logtmax - cc[nc].logtmin) / float(cc[nc].ntmax - 1);
    }

    // Check that tmin is not below the minimum temperature in the first
    // cooling curve
    //if (tmin < cc[0].te[0]){
    //  std::cout << "Minimum temperature is below that in cooling curve\n";
    //  std::cout << " tmin = " << tmin << " (K); te[0] = " << cc[0].te[0] << " (K). Aborting!\n";
    //  exit(EXIT_FAILURE);
    //}

    firstHeatCool = false;
    if (Globals::my_rank == 0){
      std::cout << "Finished firstHeatCool!\n";
    }
  }

  AthenaArray<Real> dei(pmb->ke+1,pmb->je+1,pmb->ie+1);
  
  // Now loop over cells, calculating the heating/cooling rate
  for (int k = pmb->ks; k <= pmb->ke; ++k){
    for (int j = pmb->js; j <= pmb->je; ++j){
      for (int i = pmb->is; i <= pmb->ie; ++i){
        Real rho  = cons(IDN, k, j, i);
        Real u1   = cons(IM1, k, j, i) / rho;
        Real u2   = cons(IM2, k, j, i) / rho;
        Real u3   = cons(IM3, k, j, i) / rho;
        Real v    = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
        Real ke   = 0.5 * rho * v * v;
        Real ie   = cons(IEN, k, j, i) - ke;
        Real pre  = gmma1 * ie;
        Real temp = pre * avgmass / (rho * boltzman);

        Real tempold = temp;
        Real logtemp = std::log10(temp);
        Real rhomh   = rho / massh;
        
        // Import wind "colour"
        Real col  = pmb->pscalars->s(0,k,j,i) / rho;
        // Convert wind "colour" to an array form
        Real windCol[2] = {0.0,0.0};  // Initialise array
        windCol[0] = col;             // WR wind mass fraction
        windCol[1] = 1.0 - col;       // OB wind mass fraction

        // Dust parameters, left empty unless dust cooling enabled
        Real a;      // Grain radius (micron) 
        Real z;      // Dust to gas mass ratio
        Real nH,nD;  // Hydrogen and dust number densities (cm^-3)
        if (dust_cool) {
          z = pmb->pscalars->s(ZLOC,k,j,i)/rho;
          a = pmb->pscalars->s(ALOC,k,j,i)/rho;
          // Calculate number density of gas
          nH = rho * (10.0/14.0)/massh;
          // Calculate dust grain number density
          Real rhoD  = rho * z;
          Real volD  = 4.0*PI*CUBE(a*1.0e-4);
          Real massD = grainBulkDensity * volD;
          // Finish nD calculation
          nD = rhoD / massD;
        }

        if (std::isnan(tempold) || std::isinf(tempold)){
          std::cout << "tempold is a NaN or an Inf. Aborting!\n";
          exit(EXIT_SUCCESS);
        }

        if (temp < 1.5 * Twind){
          // this should prevent wasting time in the unshocked gas
          // set temp to Twind
          temp = Twind;
        }
        else{       
          double dtint = 0.0;
          double lcool = 0.0;
          while (dtint < dt)
          { // && temp > 0.5*tempold){ // "resolve" cooling by restricting it to 10%

            //GammaHeat = 2.0e-26; // erg/s/cm^3

            // Loop over all cooling curves (e.g. gas plus dust)
            double lambda_cool_nc[ncool] = {0.0};
            double totalGasCool = 0.0;
            for (int nc = 0; nc < ncool; ++nc){
            
              double lambda_cool = 0.0;
              // If temp <= 1.1*te[0] there is zero cooling. We use te[0] instead of tmin,
              // because if tmin < te[0] it won't trigger if t > tmin but < te[0].
              // Cooling is only calculated if the temperature is
              // within the temperature range of the cooling curve.
              if (temp > 1.1 * cc[nc].te[0] && temp < cc[nc].te[cc[nc].ntmax - 1])
              {
                int kk = int((logtemp - cc[nc].logtmin) / cc[nc].dlogt);
                if (kk == cc[nc].ntmax - 1)
                  lambda_cool = cc[nc].lambdac[cc[nc].ntmax - 1];
                else
                {
                  // temperature lies between logt[kk] and logt[kk+1]
                  double grad = (cc[nc].loglambda[kk + 1] - cc[nc].loglambda[kk]) / cc[nc].dlogt;
                  lambda_cool = std::pow(10, cc[nc].loglambda[kk] + grad * (logtemp - cc[nc].logt[kk]));
                }
                //if ((temp > 1.0e5)&&(nc == 1)){
                //  cout << "temp = " << temp << "; kk = " << kk << "; lambda_cool = " << lambda_cool << "\n";
                //  cout << "te[0] = " << cc[nc].te[0] << "\n";
                //  cout << "te[ntmax-1] = " << cc[nc].te[cc[nc].ntmax-1] << "\n";
                //  cout << "logtmin = " << cc[nc].logtmin << "\n";
                //  cout << "dlogt = " << cc[nc].dlogt << "\n";
                //  quit();
                //}

                //if (Globals::my_rank == 0 && tempold > 1.0e7){
                //  std::cout << "kk = " << kk << "; ntmax = " << cc[nc].ntmax << "; lambda_cool = " << lambda_cool << "\n";
                //}
              }
              lambda_cool_nc[nc]  = lambda_cool;
              lambda_cool_nc[nc] *= windCol[nc];
              totalGasCool += lambda_cool_nc[nc];
            }

            // Calculate cooling rate due to gas
            double edotGas = rhomh * rhomh * totalGasCool; // erg/cm^-3/s
            double total_cool = edotGas;
            //if (Globals::my_rank == 0 && tempold > 1.0e7){
            //  std::cout << "edotGas = " << edotGas << "; lambda_cool_nc[0] = " << lambda_cool_nc[0] << "\n";
            //}

            // Calculate cooling rate due to dust
            // Uses a separate function as it is quite involved
            if (dust_cool) {
              for (int nc = 0; nc < 2; ++nc) {
                Real edotGrain = CalcGrainCoolRate(rho,a,temp,nc);
                Real edotDust  = nD * edotGrain;
                     edotDust *= windCol[nc];
                // Add result to cooling total
                total_cool += edotDust;
              }
            }

            double Eint = pre / gmma1;
            double t_cool = Eint / total_cool;
            double dtcool = 0.1 * std::abs(t_cool); // maximum timestep to sample the cooling curve
            if ((dtint + dtcool) > dt)
              dtcool = dt - dtint;
            dtint += dtcool;
            lcool += v*dtcool;

            // Calculate new temperature, and update pressure
            double dEint = -total_cool * dtcool;
            double tempnew = temp * (Eint + dEint) / Eint;
            tempnew = std::max(tmin, tempnew);
            tempnew = std::min(tempnew, tmax);
            pre *= tempnew / temp;
            temp = tempnew;
            /*
            if (temp <= 0.97 * tempold)
            { // maximum cooling allowed
              temp = 0.97 * tempold;
              break;
            }*/

            logtemp = std::log10(temp);

            //if (Globals::my_rank == 0 && tempold > 1.0e7){
            //  std::cout << "tempold = " << tempold << "; temp = " << temp << "\n";
            //  std::cout << "dt = " << dt << "; dtcool = " << dtcool << "; total_cool = " << total_cool << "; Eint = " << Eint << "; dEint = " << dEint << "\n";
            //}

            // Monitor total energy loss rate
            //totEdotGas += edotGas*vol;               // erg/s
#ifdef DUST
            //totEdotDust += edotDust*vol;
#endif
          } // end of cooling sub-cycling

	  /*
          // Check that cooling length is resolved and make sure that it is
          if (v != 0.0 && lcool < 10.0*dx){
            double lcool_new = 10.0*dx;
            double tcool = lcool_new/std::abs(v);
            double dtemp = dt*tempold/tcool;
            temp = tempold - dtemp;
            //std::cout << "v = " << v << "; dx = " << dx << "; lcool = " << lcool << "; dtemp = " << dtemp << "; tempold = " << tempold << "; temp = " << temp << "\n";
            //exit(EXIT_SUCCESS);
          } 
	  */


        }   // cool?

	// Calculate change in energy, dei = delta(temp) * density.
	// dei is positive if the gas is cooling.
        dei(k,j,i) = (tempold - temp)*rho;

	/*
        // Update conserved values
        temp = std::max(temp, tmin);
        temp = std::min(temp, tmax);
        if (std::isnan(temp) || std::isinf(temp))
        {
          std::cout << "temp is a NaN or an Inf. Aborting!\n";
          exit(EXIT_SUCCESS);
        }
        Real pnew = rho * temp * boltzman / avgmass;
        cons(IEN, k, j, i) = ke + pnew / gmma1;
	*/
	
        //if (Globals::my_rank == 0 && tempold > 1.0e7){
        //  std::cout << "pnew = " << pnew << "\n";
        //  exit(EXIT_SUCCESS);
        //}

        // Calculate change in energy, dei = delta(temp) * density
        //lg.dei[k][j][i] = (tempold - tempnew)*rho;

        //heatCoolRate = rhomh*(GammaHeat - rhomh*lambda_cool);

        //if (ncycle == 499 && i == 20 && procRank == 0){
        //std::cout << "Heatcool: rho = " << rho << "; pre = " << pre << "; tempold = " << tempold << "; temp = " << temp << "\n";
        //exit(EXIT_SUCCESS);
        //}

        //lg.zedot[k][j][i] += heatCoolRate;
        //if (i == 0 && y == 1.0 && procRank == 0){
        //  cout << "radiateHeatCool: "<<lambda_coolheatCoolRate<<" "<<lg.zedot[k][j][i] << "\n";
        //}
      }
    }
  }

  //if (ncycle%1000 == 0){
  //  cout << "totEdotGas = " << totEdotGas << "; totEdotDust = " << totEdotDust << " (erg/s)\n";
  //}

  // Restrict cooling at unresolved interfaces
  if (restrictUnresolvedCooling) RestrictCool(pmb->is,pmb->ie,pmb->js,pmb->je,pmb->ks,pmb->ke,pmb->pmy_mesh->ndim,gmma1,dei,cons);

  // Adjust pressure due to cooling
  AdjustPressureDueToCooling(pmb->is,pmb->ie,pmb->js,pmb->je,pmb->ks,pmb->ke,gmma1,dei,cons);

  dei.DeleteAthenaArray();
  return;
}

//========================================================================================
//! \fn void MeshBlock::FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Fix the wind temperature
//========================================================================================
void FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){

  Real g = pmb->peos->GetGamma();
  Real gmma1 = g - 1.0;

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
	
        Real temp = prim(IEN,k,j,i)*avgmass/(prim(IDN,k,j,i)*boltzman);
	if (temp < Twind){
	  Real rho = prim(IDN,k,j,i);
	  Real u1 = prim(IM1,k,j,i);
	  Real u2 = prim(IM2,k,j,i);
	  Real u3 = prim(IM3,k,j,i);
	  Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	  Real pre = (rho/avgmass)*boltzman*Twind;
          cons(IEN,k,j,i) = ke + pre/gmma1;
	}

	/*
        Real rho = cons(IDN,k,j,i);
	Real u1  = cons(IM1,k,j,i)/rho; 
	Real u2  = cons(IM2,k,j,i)/rho; 
	Real u3  = cons(IM3,k,j,i)/rho;
	Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	Real ie = cons(IEN,k,j,i) - ke;
	Real pre = gmma1*ie;
	Real temp = pre*avgmass/(rho*boltzman);
	if (temp < Twind){
	  pre = (rho/avgmass)*boltzman*Twind;
	  cons(IEN,k,j,i) = ke + pre/gmma1;
	}
	*/
      }
    }
  }
  return;
}


// Refinement condition. We look for the presence of a shock using the criteria from VH-1.
// Note that for a CWB with equal speed winds, there is no density or pressure
// jump at the stagnation point, so this will not be picked up as shock. This can lead to
// a lack or refinement at the stagnation point, despite refinement occuring off-axis.
// However, there is a strong velocity covergence at the stagnation point so this can
// also be looked for.
// In 2D RPZ geometry, i(1) is R, j(2) is P and k(3) is Z, and the stars are located along Z (k).
// Therefore there should be a velocity convergence in Z at the stagnation point.
int RefinementCondition(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> &r = pmb->pscalars->r;

  Real dx = pmb->pcoord->dx1f(0);
  Real dy = pmb->pcoord->dx2f(0);
  Real dz = pmb->pcoord->dx3f(0);
  //std::cout << "dx = " << dx << "; dy = " << dy << "; dz = " << dz << "\n";
  //exit(EXIT_SUCCESS);
  
  Real epsilon = 0.33;

  // Mesh contains: root_level, max_level, current_level;
  static bool first = true;
  if (first){
    G0_level = pmb->loc.level;
    first = false;
  }
  
  int levelAboveRoot = pmb->loc.level - G0_level;
  //std::cout << "currentLevel = " << currentLevel << "\n";
  //std::cout << "root_level = " << pmb->pmy_mesh->root_level << "\n";
  //exit(EXIT_SUCCESS);

  Real stagx, stagy, stagz;
  Real eta = mdot2*vinf2/(mdot1*vinf1);                   // wind mtm ratio
  Real dsep = std::sqrt(std::pow(xpos1 - xpos2,2) + std::pow(ypos1 - ypos2,2) + std::pow(zpos1 - zpos2,2));
  Real rob = (std::sqrt(eta)/(1.0 + std::sqrt(eta)))*dsep;
  Real fracrob = 1.0 - 1.0/(1.0 + std::sqrt(eta));
  
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
    stagx = xpos2 - fracrob*(xpos2 - xpos1);
    stagy = ypos2 - fracrob*(ypos2 - ypos1);
    stagz = 0.0;
  }
  else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    stagx = 0.0;
    stagy = 0.0;
    stagz = zpos2 - rob;
  }

  for(int k=pmb->ks-1; k<=pmb->ke+1; k++) {
    Real dz = (pmb->pcoord->x3v(k+1) - pmb->pcoord->x3v(k-1));
    for(int j=pmb->js-1; j<=pmb->je+1; j++) {
      Real dy = (pmb->pcoord->x2v(j+1) - pmb->pcoord->x2v(j-1));
      for(int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real dx = (pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i-1));
        //std::cout << "dx = " << dx << "; dy = " << dy << "; dz = " << dz << "\n";
        //exit(EXIT_SUCCESS);

	// Do not refine beyond G0 unless T > 1.0e7 K
        //Real temp = w(IPR,k,j,i)*avgmass/(w(IDN,k,j,i)*boltzman);
        //if (temp > 5.0e7 && levelAboveRoot < 3)      return 1;
	//else if (temp > 1.0e7 && levelAboveRoot < 2) return 1; // refine

	/*
	// Refine on overdensity
	Real rho = w(IDN,k,j,i);
        if (NSCALARS > 0){
	  if (r(0,k,j,i) > 0.5){ // In wind 1
            Real zc = pmb->pcoord->x3v(k) - zpos1;
            Real yc = pmb->pcoord->x2v(j) - ypos1;
            Real xc = pmb->pcoord->x1v(i) - xpos1;
    	    Real r2 = xc*xc + yc*yc + zc*zc;
	    Real rhoWind = mdot1/(4.0*pi*r2*vinf1);
	    if (rhoWind > 2.0*rho) return 1; // refine
	  }
	  else{
            Real zc = pmb->pcoord->x3v(k) - zpos2;
            Real yc = pmb->pcoord->x2v(j) - ypos2;
            Real xc = pmb->pcoord->x1v(i) - xpos2;
    	    Real r2 = xc*xc + yc*yc + zc*zc;
	    Real rhoWind = mdot2/(4.0*pi*r2*vinf2);
	    if (rhoWind > 2.0*rho) return 1; // refine
	  }
	}
	*/

	// Refine on divergence condition
	Real divV = 0.0;
        //if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
	  Real dudx = (w(IVX,k,j,i+1)-w(IVX,k,j,i-1))/dx;
	  Real dvdy = (w(IVY,k,j+1,i)-w(IVY,k,j-1,i))/dy;
	  Real dwdz = (w(IVZ,k+1,j,i)-w(IVZ,k-1,j,i))/dz;
	  divV = dudx + dvdy + dwdz;
	  /*
	}
	else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
	  Real r = pmb->pcoord->x1v(i);
	  Real rm = pmb->pcoord->x1v(i-1);
	  Real rp = pmb->pcoord->x1v(i+1);
	  Real divr = (rp*w(IVX,k,j,i+1) - rm*w(IVX,k,j,i-1))/r;
	  Real divp = 0.0;
	  Real divz = (w(IVZ,k+1,j,i)-w(IVZ,k-1,j,i))/dz;
	  divV = divr + divp + divz;
	}
	  */
	  
        //std::cout << "dudx = " << dudx << "; dvdy = " << dvdy << "; dwdz = " << dwdz << "; divV = " << divV << "\n";
        //exit(EXIT_SUCCESS);
	//if (divV < 0.0) return 1; // refine
	if (divV < 0.0){ // potentially refine
	  // Check to see if not too far away from the stagnation point.
	  Real x,y,z;
          if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
	    x = pmb->pcoord->x1v(i) - stagx;
	    y = pmb->pcoord->x2v(j) - stagy;
	    z = pmb->pcoord->x3v(k) - stagz;
	  }
	  else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
	    x = pmb->pcoord->x1v(i);
	    y = 0.0;
	    z = pmb->pcoord->x3v(k) - stagz;
	  }
	  Real r2 = x*x + y*y + z*z;
	  Real r = std::sqrt(r2);
	  int ri = int(r/dx);
	  //return 1;              // always refine
	  if (ri < 10) return 1; // refine (dx is twice as big as previously)
	  else if (ri < 30);     // do nothing
	  else return -1;	 // derefine 
	}

	  
	/*
        // Look for presence of a shock
	// Direction 1
	Real delp1 = w(IPR,k,j,i+1) - w(IPR,k,j,i-1);
	Real delp2 = w(IPR,k,j,i+2) - w(IPR,k,j,i-2);
        Real jump = std::abs(delp1)/std::min(w(IPR,k,j,i+1),w(IPR,k,j,i-1));
        Real shock = jump - epsilon;                       // Eq 46 in Fryxell etal
        shock = std::max(0.0,shock);
        if (shock > 0.0) shock = 1.0;
        if (w(IVX,k,j,i-1) < w(IVX,k,j,i+1)) shock = 0.0;  // Eq 47 in Fryxell etal
	if (shock == 1.0) return 1; // refine
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
	  // Also look for a reverse in the x-velocity (assumes that the stars are
	  // initially located along the x-xaxis)
	  if ((w(IVX,k,j,i-1) > 0.0) && (w(IVX,k,j,i+1) < 0.0)) return 1; // refine
	}
	
        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
	  // Direction 2
	  delp1 = w(IPR,k,j+1,i) - w(IPR,k,j-1);
	  delp2 = w(IPR,k,j+2,i) - w(IPR,k,j-2);
          jump = std::abs(delp1)/std::min(w(IPR,k,j+1,i),w(IPR,k,j-1,i));
          shock = jump - epsilon;                       // Eq 46 in Fryxell etal
          shock = std::max(0.0,shock);
          if (shock > 0.0) shock = 1.0;
          if (w(IVY,k,j-1,i) < w(IVY,k,j+1,i)) shock = 0.0;  // Eq 47 in Fryxell etal
	  if (shock == 1.0) return 1; // refine
	}

	// Direction 3
  	delp1 = w(IPR,k+1,j,i) - w(IPR,k-1,j,i);
        delp2 = w(IPR,k+2,j,i) - w(IPR,k-2,j,i);
        jump = std::abs(delp1)/std::min(w(IPR,k+1,j,i),w(IPR,k-1,j,i));
        shock = jump - epsilon;                       // Eq 46 in Fryxell etal
        shock = std::max(0.0,shock);
        if (shock > 0.0) shock = 1.0;
        if (w(IVZ,k-1,j,i) < w(IVZ,k+1,j,i)) shock = 0.0;  // Eq 47 in Fryxell etal
	if (shock == 1.0) return 1; // refine       
	if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
	  // Also look for a reverse in the z-velocity
	  //if ((w(IVZ,k-1,j,i) > 0.0) && (w(IVZ,k+1,j,i) < 0.0)) return 1; // refine
	}
	*/
      dontLookForShock: ;
      }
    }
  }

  // Set refinement based on distance to each star
  Real xpos[2]={xpos1,xpos2};
  Real ypos[2]={ypos1,ypos2};
  Real zpos[2]={zpos1,zpos2};

  bool derefineWind = true;
  
  for (int nw = 0; nw < 2; ++nw){ // Loop over each wind
    for (int k=pmb->ks; k<=pmb->ke; k++) {
      Real zc = pmb->pcoord->x3v(k) - zpos[nw];
      Real zc2 = zc*zc;
      for (int j=pmb->js; j<=pmb->je; j++) {
        Real yc = pmb->pcoord->x2v(j) - ypos[nw];
        Real yc2 = yc*yc;
        for (int i=pmb->is; i<=pmb->ie; i++) {
          Real xc = pmb->pcoord->x1v(i) - xpos[nw];
	  Real xc2 = xc*xc;
	  Real r2 = xc2 + yc2 + zc2;
	  Real r = std::sqrt(r2);
	  int ri = int(r/dx);
	  if (ri < 10) return 1; // refine
	  else if (ri < 20) derefineWind = false;
	}
      }
    }
  }
  
  if (derefineWind) return -1; // derefine
  return 0; // keep as is 
}


// Calculate the position and velocities of the stars based on the model time.
void OrbitCalc(Real t){
  float xdist,ydist,zdist;
  double time_offset,torbit,phase,phi,E,dE,cosE,sinE;
  double sii,coi,theta,ang;
  double rrel;      // radius vector
  double gamma_ang; // angle between velocity and radius vectors
  double a1,a2;     // semi-major axis of barycentric orbits
  double M;         // effective orbit mass
  double m1,m2,v1,v2;

  //std::cout << "dsep = " << dsep << "; t = " << t << "; phaseoff = " << phaseoff << "\n";

 
  time_offset = phaseoff*period;
  torbit = t + time_offset;  // time in seconds
  phase = torbit/period;

  //      write(*,*) T,torbit,phase
  //      quit()

  phi = 2.0*pi*phase;
  E = phi;           // first guess at eccentric anomaly

  // Now use Newton-Raphson solver to calculate E
  // (dsin, dcos, and dabs are double precision versions)

  cosE = std::cos(E);
  sinE = std::sin(E);
  
  dE = (phi - E + ecc*sinE)/(1.0 - ecc*cosE);
 iterate:
  if (std::abs(dE) > 1.0e-10){
    E = E + dE;
    cosE = std::cos(E);
    sinE = std::sin(E);
    dE = (phi - E + ecc*sinE)/(1.0 - ecc*cosE);
    goto iterate;
  }

  sii = (std::sqrt(1.0 - ecc*ecc))*sinE/(1.0 - ecc*cosE);
  coi = (cosE - ecc)/(1.0 - ecc*cosE);
  theta = std::atan2(sii,coi);
  if (theta < 0.0) theta= 2.0*pi + theta;   // 0 < theta < 2pi
         
  rrel = 1.0 - ecc*cosE;   // radius vector

  // Compute barycentric orbital positions and velocities of masses m1 and m2

  m1 = mass1;
  m2 = mass2;

  Real G = 6.67259e-8;

  M = (std::pow(m2,3)/std::pow((m1+m2),2)) * Msol;
  a1 = std::pow((G*M*period*period/(4.0*pi*pi)),(1.0/3.0));
  v1 = std::sqrt(G*M*(2.0/rrel/a1 - 1.0/a1));
  M = (std::pow(m1,3)/std::pow((m1+m2),2)) * Msol;
  a2 = std::pow((G*M*period*period/(4.0*pi*pi)),(1.0/3.0));
  v2 = std::sqrt(G*M*(2.0/rrel/a2 - 1.0/a2));

  // Compute angle between velocity and radius vectors

  gamma_ang = pi/2.0 + std::acos(std::sqrt((1.0 - ecc*ecc)/(rrel*(2.0 - rrel))));

  if (theta <= pi) ang = pi-gamma_ang+theta;
  else             ang = theta+gamma_ang;

  double sintheta = std::sin(theta);
  double costheta = std::cos(theta);
  double sinang = std::sin(ang);
  double cosang = std::cos(ang);

  // Orbit is clockwise in the xy plane, with the semi-major axis along the x-axis
  
  zpos1 = 0.0;
  zpos2 = 0.0;
  zvel1 = 0.0;
  zvel2 = 0.0;
  ypos1 =  a1*rrel*sintheta;
  xpos1 = -a1*rrel*costheta;
  ypos2 = -a2*rrel*sintheta;
  xpos2 =  a2*rrel*costheta;
  yvel1 =  v1*sinang;
  xvel1 = -v1*cosang;
  yvel2 = -v2*sinang;
  xvel2 =  v2*cosang;

  Real eta = (mdot2 * vinf2) / (mdot1 * vinf1);
  Real sqrteta = std::sqrt(eta);
  // Update global orbit phase, for use in history file
  orbitPhase = phase;
  // Update dsep and stagnation positions
  dsep = std::sqrt(SQR(xpos1 - xpos2) + (SQR(ypos1 - ypos2)));
  rWR  = (1 / (1 + sqrteta)) * dsep;
  rOB  = (sqrteta / (1 + sqrteta)) * dsep;
  // Calculate cooling parameter, xi, for each star
  xiWR = (POW4(vinf1 / 1e8) * (rWR / 1e12)) / ((mdot1 / (Msol/yr)) / 1e-7);
  xiOB = (POW4(vinf2 / 1e8) * (rOB / 1e12)) / ((mdot2 / (Msol/yr)) / 1e-7);
  
  return;
}

/*!  \brief Restrict the cooling rate at unresolved interfaces between hot 
 *          diffuse gas and cold dense gas.
 *
 *   Replace deltaE with minimum of neighboring deltaEs at the interface.
 *   Updates dei, which is positive if the gas is cooling.
 *
 *   \author Julian Pittard (Original version 13.09.11)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 13.09.11 (JMP)
 */
void RestrictCool(int is,int ie,int js,int je,int ks,int ke,int nd,Real gmma1,AthenaArray<Real> &dei,const AthenaArray<Real> &cons){

  AthenaArray<Real> pre(ke+1,je+1,ie+1);
  AthenaArray<Real> scrch(ie+1), dis(ie+1), drhox(ie+1), drhoy(ie+1), drhoz(ie+1);
  
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        Real rho = cons(IDN,k,j,i);
        Real u1  = cons(IM1,k,j,i)/rho; 
        Real u2  = cons(IM2,k,j,i)/rho; 
        Real u3  = cons(IM3,k,j,i)/rho;
        Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
        Real ie = cons(IEN,k,j,i) - ke;
        pre(k,j,i) = gmma1*ie;
      }
    }
  }
  
  for (int k = ks; k <= ke; k++) {
    int ktp = std::min(k+1,ke);
    int kbt = std::max(k-1,ks);
    for (int j = js; j <= je; j++) {
      // Locate cloud interfaces => dis = 1, otherwise, dis = -1
      int jtp = std::min(j+1,je);
      int jbt = std::max(j-1,js);
      for (int i = is+1; i < ie; i++) {
        scrch(i) = std::min(cons(IDN,k,j,i+1),cons(IDN,k,j,i-1));
        drhox(i) = cons(IDN,k,j,i+1) - cons(IDN,k,j,i-1);
        drhox(i) = std::copysign(drhox(i),(pre(k,j,i-1)/cons(IDN,k,j,i-1)
                   - pre(k,j,i+1)/cons(IDN,k,j,i+1)))/scrch(i) - 2;
        if (nd > 1) {
          scrch(i) = std::min(cons(IDN,k,jtp,i),cons(IDN,k,jbt,i));
          drhoy(i) = cons(IDN,k,jtp,i) - cons(IDN,k,jbt,i);
          drhoy(i) = std::copysign(drhoy(i),(pre(k,jbt,i)/cons(IDN,k,jbt,i)
                     - pre(k,jtp,i)/cons(IDN,k,jtp,i)))/scrch(i) - 2;
        }
        if (nd == 3) {
          scrch(i) = std::min(cons(IDN,ktp,j,i),cons(IDN,kbt,j,i));
          drhoz(i) = cons(IDN,ktp,j,i) - cons(IDN,kbt,j,i);
          drhoz(i) = std::copysign(drhoz(i),(pre(kbt,j,i)/cons(IDN,kbt,j,i)
                     - pre(ktp,j,i)/cons(IDN,ktp,j,i)))/scrch(i) - 2;
        }
        if      (nd == 1) dis(i) = drhox(i);
        else if (nd == 2) dis(i) = std::max(drhox(i),drhoy(i));
        else              dis(i) = std::max(drhox(i),std::max(drhoy(i),drhoz(i)));
      }
      dis(is) = -1.0;
      dis(ie) = -1.0;

      for (int i = is; i <= ie; i++) {
        int itp = std::min(i+1,ie);
        int ibt = std::max(i-1,is);
        if (dis(i) > 0.0) {
          if      (nd == 1) scrch(i) = std::min(dei(k,j,ibt),dei(k,j,itp));
          else if (nd == 2) scrch(i) = std::min(dei(k,j,ibt),std::min(dei(k,j,itp),
                                       std::min(dei(k,jtp,i),dei(k,jbt,i))));
          else              scrch(i) = std::min(dei(k,j,ibt),std::min(dei(k,j,itp),
                                       std::min(dei(k,jtp,i),std::min(dei(k,jbt,i),
                                       std::min(dei(ktp,j,i),dei(kbt,j,i))))));
          dei(k,j,i) = scrch(i);
        }
      }
    } // j loop
  } // k loop
  return;
}


/*!  \brief Adjust the pressure due to cooling.
 *
 *   Uses dei().
 *
 *   \author Julian Pittard (Original version 13.09.11)
 *   \version 1.0-stable (Evenstar):
 *   \date Last modified: 13.09.11 (JMP)
 */
void AdjustPressureDueToCooling(int is,int ie,int js,int je,int ks,int ke,Real gmma1,AthenaArray<Real> &dei,AthenaArray<Real> &cons) {
  for (int k = ks; k <= ke; k++){
    for (int j = js; j <= je; j++){
      for (int i = is; i <= ie; i++){
        Real rho = cons(IDN,k,j,i);
	      Real u1  = cons(IM1,k,j,i)/rho; 
	      Real u2  = cons(IM2,k,j,i)/rho; 
	      Real u3  = cons(IM3,k,j,i)/rho;
	      Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	      Real pre = (cons(IEN, k, j, i) - ke)*gmma1;
        //col  = prim[iqal0][k][j][i];
        // (=1.0 for solar, 0.0 for WC)
      	//if (col > 0.5) avgm = stavgm[0][0];
        //else           avgm = stavgm[1][0]; 
        Real avgm = avgmass;
        //if (col <= 0.5) dei[k][j][i] = 0.0;
        Real const_1 = avgm/boltzman;
        Real tmpold = const_1 * pre / rho;
        Real dtemp = dei(k,j,i)/rho;
        Real tmpnew = std::max((tmpold-dtemp),tmin);
        tmpnew = std::min(tmpnew,tmax);
        Real pnew = tmpnew * rho / const_1;
        // Update conserved values
        cons(IEN, k, j, i) = ke + pnew / gmma1;
      } // i loop
    }   // j loop
  }     // k loop
  return;
}


// Evolve the dust (e.g. due to sputtering and grain growth).
// Sputter the dust using the Draine & Salpeter (1979) prescription.
// For T > 1e6 K, the dust lifetime tau_d = 1e6 (a/n) yr (a in microns).
// where a is the grain radius and n the gas nucleon number density.
// See p156-157 of "WBB" book for further details. Sputtered atoms
// contribute to the gas density (and if are co-moving with the gas
// the gas maintains its bulk velocity and pressure).
// Grain growth occurs if T < 1.5e4 K (this gas is assumed to be cool enough
// to form dust). See p20 of CWB Book 2020.
// JMP 22/11/17 - Correctly working.
//
// *** NOTE: In this and DustCreationRateInWCR() I need to replace grainRadius
//           with an advected scalar value...  ****
//
void EvolveDust(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons){

  const Real dens_g = grainBulkDensity;               // (g/cm^3)
  const Real A = 12.0;    // Carbon dust grains
  const Real eps_a = 0.1; // probability of sticking

  Real gmma  = pmb->peos->GetGamma();
  Real gmma1 = gmma - 1.0;
  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real zc = pmb->pcoord->x3v(k) - zpos1;
    Real zc2 = zc*zc;
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real yc = pmb->pcoord->x2v(j) - ypos1;
      Real yc2 = yc*yc;
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        //cout << "i = " << i << "; j = " << j << "; k = " << k << "\n";
        Real xc = pmb->pcoord->x1v(i) - xpos1;
	      Real xc2 = xc*xc;
	      Real r2 = xc2 + yc2 + zc2;

	      Real rho = cons(IDN,k,j,i);
	      Real nH = rho*(10.0/14.0)/massh;  // solar with 10 H per 1 He, avg nucleon mass mu_nu = 14/11. 
        Real ntot = 1.1*nH;                        // total nucleon number density
        Real z = pmb->pscalars->s(1,k,j,i)/rho;    // dust mass fraction
        Real rhod = rho*z;                         // dust mass density (g/cm^3)

        Real a = (pmb->pscalars->s(2,k,j,i)/rho)*1.0e-4; // grain radius (cm)
        Real massD = (4.0/3.0)*pi*std::pow(a,3)*dens_g; // grain mass (g)

	      Real nD = rhod/massD;                      // grain number density (cm^-3)

	      Real u1  = cons(IM1,k,j,i)/rho; 
	      Real u2  = cons(IM2,k,j,i)/rho; 
	      Real u3  = cons(IM3,k,j,i)/rho;
	      Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	      Real pre = (cons(IEN, k, j, i) - ke)*gmma1;

        Real col = pmb->pscalars->s(0,k,j,i)/rho;  // wind colour (1.0 for primary wind, 0.0 for secondary wind)
        Real cols[2];
        cols[0] = col;
        cols[1] = 1.0 - col;
	
        //col  = prim[iqal0][k][j][i];   // (=1.0 for solar, 0.0 for WC)
      	//if (col > 0.5) avgm = stavgm[0][0];
        //else           avgm = stavgm[1][0]; 
        Real avgm = avgmass;     
	      Real temp = avgm * pre / (rho*boltzman);

      	// Determine overdensity from smooth WR wind
	      Real rho_smooth =  mdot1/(4.0*pi*r2*vinf1);
	
        Real rhod_dot = 0.0;
        Real dadt = 0.0;
	
        if (temp > 1.0e6 && z > 0.0){
          // Dust thermal sputtering
          Real tau_d = 3.156e17*a/ntot;            // grain destruction time (s)
          dadt = -a/tau_d;
          rhod_dot = -1.33e-17*dens_g*a*a*ntot*nD; // dust destruction rate (g cm^-3 s^-1)
        }
        else if (temp < 1.4e4){
          for (int nw = 0; nw < 2; nw++) {
            Real windFrac = cols[nw];
            Real carbonMassFraction = massFrac[nw][2];
            // Dust growth. Requires some grains to exist otherwise rhod_dot = 0.0	  
            Real wa  = std::sqrt(3.0*boltzman*temp/(A*massh));
            dadt     = 0.25*eps_a*rho*wa/dens_g;
            dadt    *= windFrac;
            dadt    *= carbonMassFraction;
            rhod_dot = 4.0*pi*a*a*dens_g*nD*dadt;    // dust growth rate (g cm^-3 s^-1)
          }
        }
        if (rhod_dot != 0.0){
          Real drhod = rhod_dot*dt;
          Real rhodnew = std::max(minimumDustToGasMassRatio*rho,rhod + drhod);    // new dust density
          Real rhonew = rho + (rhod - rhodnew);                                   // new gas  density (preserving total mass)
          cons(IDN,k,j,i) = rhonew;
          // Update the conserved wind colour (as the gas density may have changed)
          pmb->pscalars->s(0,k,j,i) = col*rhonew;
          // Update the dust to gas mass ratio
          pmb->pscalars->r(1,k,j,i) = rhodnew/rhonew;
          pmb->pscalars->s(1,k,j,i) = pmb->pscalars->r(1,k,j,i)*rhonew;
          // Update the grain radius
          Real anew = std::max(minimumGrainRadiusMicrons,(a + dadt*dt)/1.0e-4); // (microns)
          pmb->pscalars->r(2,k,j,i) = anew;
          pmb->pscalars->s(2,k,j,i) = anew*rhonew;
        }
      }
    }
  }
  
  return;
}

// \brief This function is used to link all the User History file functions
// While not particularly efficient, it's still better than nothing!

Real UserHistoryFunction(MeshBlock *pmb, int iout) {
  // indexes 0 through 7, orbital parameters
  // These have already been calculated every time calcOrbit() is called
  // Hence nothing has to be done, so running these first saves time
  if      (iout == 0)  {return orbitPhase;}
  else if (iout == 1)  {return xpos1;}
  else if (iout == 2)  {return ypos1;}
  else if (iout == 3)  {return xpos2;}
  else if (iout == 4)  {return ypos2;}
  else if (iout == 5)  {return dsep;}
  else if (iout == 6)  {return rWR;}
  else if (iout == 7)  {return rOB;}
  else if (iout == 8)  {return xiWR;}
  else if (iout == 9)  {return xiOB;}
  // Dust creation rate calculator, 
  else if (iout == 10) {return DustCreationRateInWCR(pmb,0);}
  else if (iout == 11) {return DustCreationRateInWCR(pmb,1);}
  else if (iout == 12) {return DustCreationRateInWCR(pmb,2);}
  else if (iout == 13) {return DustCreationRateInWCR(pmb,3);}
  else if (iout == 14) {return DustCreationRateInWCR(pmb,4);}
  // Somethings up, so return
  return 0.0;
}



Real DustCreationRateInWCR(MeshBlock *pmb, int iout){
  const Real dens_g = grainBulkDensity;               // (g/cm^3)
  const Real A = 12.0;    // Carbon dust grains
  const Real eps_a = 0.1; // probability of sticking
  
  Real gmma  = pmb->peos->GetGamma();
  Real gmma1 = gmma - 1.0;

  AthenaArray<Real>& cons = pmb->phydro->u;
  
  Real dmdustdt_WCR = 0.0;
  Real dmdust_WCR_dt_created   = 0.0;
  Real dmdust_WCR_dt_destroyed = 0.0;
  Real dust_WCR   = 0.0;
  Real dust_TOTAL = 0.0;

  Real a_TOTAL = 0.0;
  Real z_TOTAL = 0.0;
  Real vol_TOTAL = 0.0;
  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real zc = pmb->pcoord->x3v(k) - zpos1;
    Real zc2 = zc*zc;
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real yc = pmb->pcoord->x2v(j) - ypos1;
      Real yc2 = yc*yc;
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real xc = pmb->pcoord->x1v(i) - xpos1;
        Real xc2 = xc*xc;
        Real vol = pmb->pcoord->GetCellVolume(k, j, i); // Cell volume (cm)
        Real r2 = xc2 + yc2 + zc2;
        Real rho = cons(IDN,k,j,i);
        Real u1  = cons(IM1,k,j,i)/rho;
        Real u2  = cons(IM2,k,j,i)/rho;
        Real u3  = cons(IM3,k,j,i)/rho;
        Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
        Real pre = (cons(IEN, k, j, i) - ke)*gmma1;
        Real col = pmb->pscalars->s(0,k,j,i)/rho;  // wind colour (1.0 for primary wind, 0.0 for secondary wind)
        Real cols[2];
             cols[0] = col;
             cols[1] = 1.0 - col;
        //if (col > 0.5) avgm = stavgm[0][0];
        //else           avgm = stavgm[1][0]; 
        Real z   = pmb->pscalars->s(1,k,j,i)/rho;  // dust mass fraction
        Real a   = pmb->pscalars->s(2,k,j,i)/rho;  // Grain radius (micron)
             a  *= 1.0e-4;                         // Convert to cm
        Real nH = rho*(10.0/14.0)/massh;           // solar with 10 H per 1 He, avg nucleon mass mu_nu = 14/11. 
        Real ntot = 1.1*nH;                        // total nucleon number density
        Real rhod = rho*z;                         // dust mass density (g/cm^3)
        Real massD = (4.0/3.0)*pi*std::pow(a,3)*dens_g; // grain mass (g)
        Real nD = rhod/massD;                      // grain number density (cm^-3)
        Real avgm = avgmass;     
        Real temp = avgm * pre / (rho*boltzman);
        // Determine overdensity from smooth WR wind
        Real rho_smooth =  mdot1/(4.0*pi*r2*vinf1);
        Real rhod_dot = 0.0;

        Real r = std::sqrt(r2);

        // Dust only reasonably grows 
        if (r > remapRadius1) {
          if (temp < 1.4e4){
            for (int nw = 0; nw < 2; nw++) {
              // For each wind, determine the carbon content and wind mass fraction
              Real windFraction       = cols[nw];
              Real carbonMassFraction = massFrac[nw][2];
              // Dust growth. Requires some grains to exist otherwise rhod_dot = 0.0	  
              Real wa    = std::sqrt(3.0*boltzman*temp/(A*massh));
              Real dadt  = 0.25*eps_a*rho*wa/dens_g;
                   dadt *= windFraction;
                   dadt *= carbonMassFraction;
              rhod_dot   = 4.0*pi*a*a*dens_g*nD*dadt;    // dust growth rate (g cm^-3 s^-1)
            }
          }
          if (temp > 1.0e6 && z > 0.0) {
            Real tauD = 3.156e17 * a / ntot;
            Real dadt = -a / tauD;
            rhod_dot  = -1.33e-17 * dens_g * SQR(a) * ntot * nD;
          }
          if (rho > 2.0*rho_smooth && col > 0.5) {
            if (rhod_dot > 0.0) {
              dmdust_WCR_dt_created += rhod_dot * vol;
            }
            if (rhod_dot < 0.0) {
              dmdust_WCR_dt_destroyed += rhod_dot * vol;
            }
            dust_WCR += rhod*vol;
            dmdustdt_WCR += rhod_dot * vol;
          }
          dust_TOTAL += rhod*vol;
        }
        a_TOTAL   += a * vol;
        z_TOTAL   += z * vol;
        vol_TOTAL += vol;
      }
    }
  }

  Real a_avg = a_TOTAL / vol_TOTAL;
  Real z_avg = z_TOTAL / vol_TOTAL;
  
  // Quick check for 
  if (iout == 0) {
    if (std::isnan(dmdustdt_WCR) || std::isinf(dmdustdt_WCR)) {
      return 0;
    }
    else return dmdustdt_WCR;
  }
  else if (iout == 1) {
    if (std::isnan(dmdust_WCR_dt_created) || std::isinf(dmdust_WCR_dt_created)) {
      return 0;
    }
    else return dmdust_WCR_dt_created;
  }
  else if (iout == 2) {
    if (std::isnan(dmdust_WCR_dt_destroyed) || std::isinf(dmdust_WCR_dt_destroyed)) {
      return 0;
    }
    else return dmdust_WCR_dt_destroyed;
  }
  else if (iout == 3) {
    if (std::isnan(dust_WCR) || std::isinf(dust_WCR)) {
      return 0;
    }
    else return dust_WCR;
  }
  else if (iout == 4) {
    if (std::isnan(dust_TOTAL) || std::isinf(dust_TOTAL)) {
      return 0;
    }
    else return dust_TOTAL;
  }
  return 0;

}

Real ReturnOrbitalProperties(MeshBlock *pmb, int iout) {
       if (iout == 0) {return orbitPhase;}
  else if (iout == 1) {return xpos1;}
  else if (iout == 2) {return ypos1;}
  else if (iout == 3) {return xpos2;}
  else if (iout == 4) {return ypos2;}
  else if (iout == 5) {return dsep;}
  return 0;
}

Real ReturnStag(MeshBlock *pmb, int iout) {
       if (iout == 0) {return rWR;}
  else if (iout == 1) {return rOB;}
  return 0;
}

Real ReturnXi(MeshBlock *pmb, int iout) {
       if (iout == 0) {return xiWR;}
  else if (iout == 1) {return xiOB;}
  return 0;
}