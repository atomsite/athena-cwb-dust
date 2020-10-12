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

bool cool, dust;

Real initialDustToGasMassRatio;
Real initialGrainRadiusMicrons;

Real mass1, mass2;
Real mdot1, mdot2, vinf1, vinf2;
Real remapRadius1, remapRadius2;
Real xpos1, xpos2, ypos1, ypos2, zpos1, zpos2;
Real xvel1, xvel2, yvel1, yvel2, zvel1, zvel2;

Real period;  // orbit period (s)  
Real phaseoff;// phase offset of orbit (from periastron) 
Real ecc;     // orbit eccentricity

Real tmin,tmax; // min/max temperature allowed on grid

int G0_level = 0;

// Physical and mathematical constants
const Real pi = 2.0*asin(1.0);
const Real Msol = 1.9891e33;
const Real yr = 3.15569e7;
const Real boltzman = 1.380658e-16;
const Real massh = 1.67e-24;

// User defined constants
const Real minimumDustToGasMassRatio = 1.0e-6;
const Real minimumGrainRadiusMicrons = 0.01;
const Real grainBulkDensity = 3.0;                       // (g/cm^3)
const Real Twind = 1.0e4;                                // K
const Real avgmass = 1.0e-24;     // g


// Structure to store cooling curve data
struct coolingCurve{
  int ntmax;
  std::string coolCurveFile;
  std::vector<double> logt,lambdac,te,loglambda;
  double t_min,t_max,logtmin,logtmax,dlogt;
};

// JMP prototypes
void AdjustPressureDueToCooling(int is,int ie,int js,int je,int ks,int ke,Real gmma1,AthenaArray<Real> &dei,AthenaArray<Real> &cons);
Real DustCreationRateInWCR(MeshBlock *pmb, int iout);
void EvolveDust(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
//void FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void OrbitCalc(Real t);
void PhysicalSources(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
int RefinementCondition(MeshBlock *pmb);
void RestrictCool(int is,int ie,int js,int je,int ks,int ke,int nd,Real gmma1,AthenaArray<Real> &dei,const AthenaArray<Real> &cons);


//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor. Also called when restarting.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  //std::cout << "[Mesh::InitUserMeshData]\n";
  
  mdot1 = pin->GetReal("problem","mdot1");
  mdot2 = pin->GetReal("problem","mdot2");
  vinf1 = pin->GetReal("problem","vinf1");
  vinf2 = pin->GetReal("problem","vinf2");

  remapRadius1 = pin->GetReal("problem","remapRadius1");
  remapRadius2 = pin->GetReal("problem","remapRadius2");
  
  xpos1 = pin->GetReal("problem","xpos1");
  ypos1 = pin->GetReal("problem","ypos1");
  zpos1 = pin->GetReal("problem","zpos1");
  xpos2 = pin->GetReal("problem","xpos2");
  ypos2 = pin->GetReal("problem","ypos2");
  zpos2 = pin->GetReal("problem","zpos2");

  mass1 = pin->GetReal("problem","mass1");
  mass2 = pin->GetReal("problem","mass2");
  ecc   = pin->GetReal("problem","ecc");
  period = pin->GetReal("problem","period");
  phaseoff = pin->GetReal("problem","phaseoff");

  std::string cooling = pin->GetString("problem","cooling");
  if      (cooling == "on")  cool = true;
  else if (cooling == "off") cool = false;
  else{
    std::cout << "cooling value not recognized: " << cooling << "; Aborting!\n";
    exit(EXIT_SUCCESS);
  }

  std::string dusty = pin->GetString("problem","dust");
  if      (dusty == "on")  dust = true;
  else if (dusty == "off") dust = false;
  else{
    std::cout << "dust value not recognized: " << dust << "; Aborting!\n";
    exit(EXIT_SUCCESS);
  }
  if (dust && NSCALARS < 3){
    // Scalars are: 0 = wind colour
    //              1 = dust to gas mass ratio
    //              2 = dust grain radius (microns)
    std::cout << "Not enough scalars for dust modelling. NSCALARS = " << NSCALARS << ". Aborting!\n";
    exit(EXIT_SUCCESS);
  }
  if (dust){
    initialDustToGasMassRatio = pin->GetReal("problem","initialDustToGasMassRatio");
    initialGrainRadiusMicrons = pin->GetReal("problem","initialGrainRadiusMicrons");
  }

  mdot1 *= Msol/yr;
  mdot2 *= Msol/yr;

  // Note: it is only possible to have one source functions enrolled by the user.
  EnrollUserExplicitSourceFunction(PhysicalSources);

  if (adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);  

  // Add a user-defined global output (see https://github.com/PrincetonUniversity/athena-public-version/wiki/Outputs)
  if (dust){
    AllocateUserHistoryOutput(1);
    EnrollUserHistoryOutput(0, DustCreationRateInWCR, "dmdustdt_WCR");
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
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) { // RPZ
	      u1 = vinf1*sinphi;
	      u2 = 0.0;
	      u3 = vinf1*cosphi;
	    }
	    //std::cout << "nw = " << nw << "; r = " << r << "; rho = " << rho << "\n";
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	    
	    // Set passive scalars
	    if (NSCALARS > 0) {
	      // wind "colour"
              pscalars->s(0,k,j,i) = scalar[nw]*rho;
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

  if (cool) RadiateHeatCool(pmb,dt,cons);
  if (dust) EvolveDust(pmb,dt,cons);
  return;
}

//========================================================================================
//! \fn void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons);
//  \brief Calculating radiative heating and cooling
//========================================================================================
void RadiateHeatCool(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons){

  const int ncool = 1; // number of cooling curves

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
    cc[0].coolCurveFile = "cooling_KI02_4.0_CLOUDY_7.6_MEKAL.txt";
    cc[0].ntmax = 161;
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

        Real rho = cons(IDN, k, j, i);
        Real u1 = cons(IM1, k, j, i) / rho;
        Real u2 = cons(IM2, k, j, i) / rho;
        Real u3 = cons(IM3, k, j, i) / rho;
        Real v = std::sqrt(u1 * u1 + u2 * u2 + u3 * u3);
        Real ke = 0.5 * rho * v * v;
        Real ie = cons(IEN, k, j, i) - ke;
        Real pre = gmma1 * ie;
        Real temp = pre * avgmass / (rho * boltzman);

        Real tempold = temp;
        Real logtemp = std::log10(temp);
        Real rhomh = rho / massh;

#ifdef DUST
        // In this simplest implementation the dust moves with the gas and its
        // mass fraction is given by an advected scalar
        //	z = lg.P0[iqal0][k][j][i];            // dust mass fraction
        //rhod = rho*z;                        // dust mass density (g/cm^3)
        //nD = rhod/massD;                      // grain number density (cm^-3)
        //nH = rho*(10.0/14.0)/massh;          // solar with 10 H per 1 He, avg nucleon mass mu_nu = 14/11.
        ////ntot = 1.1*nH;                        // total nucleon number density
        ////ne = 1.2*nH;                          // electron number density
#endif

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
              lambda_cool_nc[nc] = lambda_cool;
            }

            // Calculate cooling rate due to gas
            double edotGas = rhomh * rhomh * lambda_cool_nc[0]; // erg/cm^-3/s
            double total_cool = edotGas;
            //if (Globals::my_rank == 0 && tempold > 1.0e7){
            //  std::cout << "edotGas = " << edotGas << "; lambda_cool_nc[0] = " << lambda_cool_nc[0] << "\n";
            //}

#ifdef DUST
            // Calculate cooling rate due to dust
            double edotDust = nH * nD * lambda_cool_nc[1]; // due to dust
            //if (temp > 1.0e5){
            //   cout << "nH = " << nH << "; nD = " << nD << "; lambda_D = " << lambda_cool_nc[1] << "\n";
            //   quit();
            //}
            total_cool += edotDust;
#endif

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
  
  // Update dsep and rob
  //xdist = xpos1 - xpos2;
  //ydist = ypos1 - ypos2;
  //zdist = zpos1 - zpos2;
  //dsep = std::sqrt(xdist*xdist + ydist*ydist + zdist*zdist);
  //eta = mdot2*vinf2/(mdot1*vinf1);                   // wind mtm ratio
  //rob = (1.0 - 1.0/(1.0 + std::sqrt(eta)))*dsep;     //distance of stagnation point from star 1 (distance from star 0 is rwr)
;
  
  //std::cout << "xpos1 = " << xpos1 << "; xpos2 = " << xpos2 << "; ypos1 = " << ypos1 << "; ypos2 = " << ypos2 << "; zpos1 = " << zpos1 << "; zpos2 = " << zpos2 << "\n";
  //std::cout << "dsep = " << dsep << "\n";
  //exit(EXIT_SUCCESS);
  
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
  
  for (int k = ks; k <= ke; k++){
    for (int j = js; j <= je; j++){
      for (int i = is; i <= ie; i++){
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
  
  for (int k = ks; k <= ke; k++){
    int ktp = std::min(k+1,ke);
    int kbt = std::max(k-1,ks);
    for (int j = js; j <= je; j++){

      // Locate cloud interfaces => dis = 1, otherwise, dis = -1

      int jtp = std::min(j+1,je);
      int jbt = std::max(j-1,js);
      for (int i = is+1; i < ie; i++){
	scrch(i) = std::min(cons(IDN,k,j,i+1),cons(IDN,k,j,i-1));
	drhox(i) = cons(IDN,k,j,i+1) - cons(IDN,k,j,i-1);
        drhox(i) = std::copysign(drhox(i),(pre(k,j,i-1)/cons(IDN,k,j,i-1)
			     - pre(k,j,i+1)/cons(IDN,k,j,i+1)))/scrch(i) - 2;
	if (nd > 1){
  	  scrch(i) = std::min(cons(IDN,k,jtp,i),cons(IDN,k,jbt,i));
	  drhoy(i) = cons(IDN,k,jtp,i) - cons(IDN,k,jbt,i);
          drhoy(i) = std::copysign(drhoy(i),(pre(k,jbt,i)/cons(IDN,k,jbt,i)
			     - pre(k,jtp,i)/cons(IDN,k,jtp,i)))/scrch(i) - 2;
        }
	if (nd == 3){
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

      for (int i = is; i <= ie; i++){
	int itp = std::min(i+1,ie);
	int ibt = std::max(i-1,is);
	if (dis(i) > 0.0){
  	  if      (nd == 1) scrch(i) = std::min(dei(k,j,ibt),dei(k,j,itp));
	  else if (nd == 2) scrch(i) = std::min(dei(k,j,ibt),std::min(dei(k,j,itp),
					  std::min(dei(k,jtp,i),dei(k,jbt,i))));
          else              scrch(i) = std::min(dei(k,j,ibt),std::min(dei(k,j,itp),
				     std::min(dei(k,jtp,i),std::min(dei(k,jbt,i),
				     std::min(dei(ktp,j,i),dei(kbt,j,i))))));
	  //std::cout << "k = " << k << "; j = " << j << "; i = " << i << "; scrch = " << scrch(i) << "; dei = " << dei(k,j,i) << "\n";
	  //exit(EXIT_SUCCESS);
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
void AdjustPressureDueToCooling(int is,int ie,int js,int je,int ks,int ke,Real gmma1,AthenaArray<Real> &dei,AthenaArray<Real> &cons){
  
  for (int k = ks; k <= ke; k++){
    for (int j = js; j <= je; j++){
      for (int i = is; i <= ie; i++){
        
        Real rho = cons(IDN,k,j,i);
	Real u1  = cons(IM1,k,j,i)/rho; 
	Real u2  = cons(IM2,k,j,i)/rho; 
	Real u3  = cons(IM3,k,j,i)/rho;
	Real ke = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
	Real pre = (cons(IEN, k, j, i) - ke)*gmma1;

        //col  = prim[iqal0][k][j][i];   // (=1.0 for solar, 0.0 for WC)

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
	Real nH = rho*(10.0/14.0)/massh;           // solar with 10 H per 1 He, avg nucleon mass mu_nu = 14/11. 
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
	  // Dust growth. Requires some grains to exist otherwise rhod_dot = 0.0	  
	  Real wa = std::sqrt(3.0*boltzman*temp/(A*massh));
	  dadt = 0.25*eps_a*rho*wa/dens_g;
	  rhod_dot = 4.0*pi*a*a*dens_g*nD*dadt;    // dust growth rate (g cm^-3 s^-1)
	  //if (rho > 2.0*rho_smooth && col > 0.5){
	    // In cool WCR
	    //std::cout << "Growing grains inside WCR...\n";
	    //std::cout << "dadt = " << dadt << "; nD = " << nD << " (cm^-3); rhod_dot = " << rhod_dot << "\n";
	    //std::cout << "Aborting!\n";
	    //exit(EXIT_SUCCESS);
	  //}
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



Real DustCreationRateInWCR(MeshBlock *pmb, int iout){
  const Real dens_g = grainBulkDensity;               // (g/cm^3)
  const Real A = 12.0;    // Carbon dust grains
  const Real eps_a = 0.1; // probability of sticking
  
  Real gmma  = pmb->peos->GetGamma();
  Real gmma1 = gmma - 1.0;

  AthenaArray<Real>& cons = pmb->phydro->u;
  
  Real dmdustdt_WCR = 0.0;
  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    Real zc = pmb->pcoord->x3v(k) - zpos1;
    Real zc2 = zc*zc;
    for (int j=pmb->js; j<=pmb->je; ++j) {
      Real yc = pmb->pcoord->x2v(j) - ypos1;
      Real yc2 = yc*yc;
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real xc = pmb->pcoord->x1v(i) - xpos1;
	Real xc2 = xc*xc;
	Real r2 = xc2 + yc2 + zc2;

	Real rho = cons(IDN,k,j,i);
	Real nH = rho*(10.0/14.0)/massh;           // solar with 10 H per 1 He, avg nucleon mass mu_nu = 14/11. 
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
	
	//if (col > 0.5) avgm = stavgm[0][0];
        //else           avgm = stavgm[1][0]; 
        Real avgm = avgmass;     
	Real temp = avgm * pre / (rho*boltzman);

	// Determine overdensity from smooth WR wind
	Real rho_smooth =  mdot1/(4.0*pi*r2*vinf1);


	
	Real rhod_dot = 0.0;
	
	
        if (temp < 1.4e4){
	  // Dust growth. Requires some grains to exist otherwise rhod_dot = 0.0	  
	  Real wa = std::sqrt(3.0*boltzman*temp/(A*massh));
	  Real dadt = 0.25*eps_a*rho*wa/dens_g;
	  rhod_dot = 4.0*pi*a*a*dens_g*nD*dadt;    // dust growth rate (g cm^-3 s^-1)
	  Real r = std::sqrt(r2);
	  if (rho > 2.0*rho_smooth && col > 0.5 && r > remapRadius1){
	    // In cool WCR
            Real vol = pmb->pcoord->GetCellVolume(k, j, i);
	    dmdustdt_WCR += rhod_dot*vol;
	  }
	}
      }
    }
  }
  return dmdustdt_WCR;
}
