//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file default_pgen.cpp
//  \brief Problem generator for Colliding Wind Binary problem.
//  Works in catersian and cylindrical coordinates.

//  CWB PROBLEM:

//  DUST PHYSICS:
//  Dust physics are calculated on passive scalars, dust is stored in the form:
//  a: grain radius in cm
//  z: dust density, g/cm^3
//  Dust is non-interactive with gas and assumed to be co-moving

//  COOLING:
//  Plasma cooling curve is used for stellar wind, while dust is calculated
//  based on the grain radius
//  System can read in a Plasma cooling curve in the form:
//  log10(T) | Lambda

//  ORBITS:
//  Orbits are calculated "on-rails" by changing the remap zones per timestep

// C headers
#include "stdio.h"

// C++ headers
#include <algorithm>
#include <cmath>
#include <cstdio>     // fopen(), fprintf(), freopen()
#include <cstring>    // strcmp()
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <chrono>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"
#include "../scalars/scalars.hpp"

// CWB problem header
#include "cwb.hpp"

// Function declaration


//  CWB global namespace
//  - Globally accessible variables that are needed for the programme to work
//    without passing a ridiculous number of variables through each function
//  - These are typically constants, aside from certain properties within
//    the WR and 
//  - A copy is kept per CPU, all updating calculations are performed across
//    all threads, this prevents race conditions for arising at a neglible cost
//    to performance, especially when compared to the hydro loop as a whole

namespace {
  // Define objects
  // Create star objects, see cwb.hpp
  Star WR;  // Object describing properties of Wolf-Rayet star
  Star OB;  // Object describing properties of OB star
  // Vector containing pointers to both stars
  std::vector<Star*> Stars;

  // Create cooling curve class
  CoolCurve CCurve;  // Object containing the cooling curve
  // Star orbital parameters
  Real dsep;          // Star orbital separation (cm)
  Real eta;           // Wind momentum ratio, OB numerators
  Real rob;           // OB distance from stagnation point (cm)
  Real rwr;           // WR distnace from stagnation point (cm)
  Real phase_offset;  // Orbital phase offset
  Real period;        // Orbital period (s)
  Real ecc;           // Orbital eccentricity

  // Dust variables
  Real a_min;      // Minimum grain size (cm)
  Real z_min;      // Minimum dust-to-gas mass ratio
  Real bulk_dens;  // Grain bulk density (g cm^-3)
  Real stick_eff;  // Grain sticking efficiency (fraction)
  Real nuc_temp;   // Grain nucleation temperature (K)

  // Refinement parameters
  int G0_level = 0;  // Coarse cell leve, leave at zero for most things

  // Simulation parameters
  int  remap     = 3;      // Remap radius in cell widths
  bool orbit     = false;  // Turn on orbital dynamics
  bool cool      = false;  // Turn on cooling
  bool dust      = false;  // Turn on dust evolution
  bool dust_cool = false;  // Turn on dust cooling

  // Simulation variables
  int coord; // Coordinate system, 0 for cartesian, 1 for cylindrical

  // Constants, all values are in CGS
  Real massh  = 1.6735575e-24;  // Hydrogen atom mass (g)
  Real kboltz = 1.380649e-16;   // Boltzmann constant (erg K^-1)
  Real G      = 6.67259e-8;     // Gravitational constant (dyn cm^2 g^-2)
  Real Msol   = 1.9891e33;      // 1 Solar mass (g)
  Real yr     = 3.15569e7;      // Seconds in a year (s)
  // Thermodynamic properties
  Real g  = 5.0 / 3.0;  // Ratio of specific heats, gamma
  Real g1 = g - 1;     // Gamma-1
}

void CalcOrbit(Real t) {
  Real time_offset;  // Time offset (s)
  Real t_orbit;      // Orbital time, inclusive of offset (s)
  Real phase;        // Orbital phase, fraction

  // Keplerian orbital parameters
  Real phi;          // True anomaly (rad)
  Real E;            // Eccentric anomaly
  Real sinE, cosE;   // Sine and cosine components of E
  Real dE;           // Change in eccentric anomaly

  Real sii,coi;
  Real a1,a2;
  Real v1,v2;
  Real theta,ang;
  Real rrel;
  Real M,gamma_ang;
  Real dx,dy; 

  time_offset = phase_offset * period;
  t_orbit     = t + time_offset;
  phase       = t_orbit / period;

  phi  = 2.0 * PI * phase; // Calculate true anomaly
  E    = phi;              // First guess at eccentric anomaly
  sinE = std::sin(E);      // Calculate sine component
  cosE = std::cos(E);      // Calculate cosine component
  dE   = (phi - E + ecc * sinE)/(1.0 - ecc*cosE);

  while (std::abs(dE) > 1.0e-10) {
    E = E + dE;
    sinE = std::sin(E);
    cosE = std::cos(E);
    dE = (phi - E + ecc * sinE) / (1.0 - ecc * cosE);
  }

  sii = (std::sqrt(1.0 - SQR(ecc))) * sinE / (1.0 - ecc*cosE);
  coi = (cosE - ecc) / (1.0 - ecc*cosE);
  theta = std::atan2(sii,coi);
  // Ensure that theta is between 0 and 2*PI
  if (theta < 0.0) {theta = 2.0 * PI + theta;}
  // Calc radius vector
  rrel = 1.0 - ecc*cosE;

  // Compute barycentric masses
  M  = (std::pow(OB.mass,3)/std::pow((WR.mass+OB.mass),2)) * Msol;
  a1 = std::pow((G*M*period*period/(4.0*PI*PI)),(1.0/3.0));
  v1 = std::sqrt(G*M*(2.0/rrel/a1 - 1.0/a1));
  M  = (std::pow(WR.mass,3)/std::pow((WR.mass+OB.mass),2)) * Msol;
  a2 = std::pow((G*M*period*period/(4.0*PI*PI)),(1.0/3.0));
  v2 = std::sqrt(G*M*(2.0/rrel/a2 - 1.0/a2));

  // Compute angle between velocity and radius vectors

  gamma_ang = PI/2.0 + std::acos(std::sqrt((1.0 - ecc*ecc)/(rrel*(2.0 - rrel))));

  if (theta <= PI) {ang = PI - gamma_ang + theta;}
  else             {ang = theta + gamma_ang;}

  // Orbit is in the xy plane, with the semi-major axis along the y-axis
  // Hence, z is not updated, and is zeroed on init
  WR.x[0] =  a1 * rrel * std::sin(theta);
  WR.x[1] = -a1 * rrel * std::cos(theta);
  OB.x[0] = -a2 * rrel * std::sin(theta);
  OB.x[1] =  a2 * rrel * std::cos(theta);
  WR.v[0] =  v1 * std::sin(ang);
  WR.v[1] = -v1 * std::cos(ang);
  OB.v[0] = -v2 * std::sin(ang);
  OB.v[1] =  v2 * std::cos(ang);
  // Update dsep and rob
  dx   = WR.x[0] - OB.x[0];
  dy   = WR.x[1] - OB.x[1];
  dsep = std::sqrt(SQR(dx) + SQR(dy));
  rob  = (std::sqrt(eta)/(1.0 + std::sqrt(eta))) * dsep;
  rwr  = dsep - rob;
  return;
}

//! \fn RestrictCool()
//  \brief Restrict the cooling rate at unresolved interfaces between hot 
//         diffuse gas and cold dense gas.
// 
//   Replace deltaE with minimum of neighboring deltaEs at the interface.
//   Updates dei, which is positive if the gas is cooling.
// 
//   \author Julian Pittard (Original version 13.09.11)
//   \version 1.0-stable (Evenstar):
//   \date Last modified: 13.09.11 (JMP)

void RestrictCool(int is,int ie,int js,int je,int ks,int ke,int nd,AthenaArray<Real> &dei,const AthenaArray<Real> &cons){

  AthenaArray<Real> pre(ke+1,je+1,ie+1);
  AthenaArray<Real> scrch(ie+1), dis(ie+1), drhox(ie+1), drhoy(ie+1), drhoz(ie+1);
  
  for (int k = ks; k <= ke; k++){
    for (int j = js; j <= je; j++){
      for (int i = is; i <= ie; i++){
        Real rho    = cons(IDN,k,j,i);
        Real u1     = cons(IM1,k,j,i)/rho;
        Real u2     = cons(IM2,k,j,i)/rho;
        Real u3     = cons(IM3,k,j,i)/rho;
        Real ke     = 0.5*rho*(u1*u1 + u2*u2 + u3*u3);
        Real ie     = cons(IEN,k,j,i) - ke;
        pre (k,j,i) = g1*ie;
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
          if      (nd == 1) scrch(i) = std::min(dei(k,j,ibt),
                                      dei(k,j,itp));
          else if (nd == 2) scrch(i) = std::min(dei(k,j,ibt),
                                      std::min(dei(k,j,itp),
                                      std::min(dei(k,jtp,i),
                                      dei(k,jbt,i))));
          else              scrch(i) = std::min(dei(k,j,ibt),
                                      std::min(dei(k,j,itp),
                                      std::min(dei(k,jtp,i),
                                      std::min(dei(k,jbt,i),
                                      std::min(dei(ktp,j,i),
                                      dei(kbt,j,i))))));
          dei(k,j,i) = scrch(i);
        }
      }
    } // j loop
  } // k loop
  
  return;
}

// Finish up cooling step by updating conserved value of internal energy
// This is accomplished by calculat
// This separate step (and the additional memory overhead of dei array) 
// is needed in order to restrict cooling at unresolved interfaces
// using the restrictcool() function

void AdjustPressureDueToCooling(int is, int ie,
                                int js, int je,
                                int ks, int ke,
                                AthenaArray<Real> &dei,
                                AthenaArray<Real> &cons) {
  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        Real rho = cons(IDN,k,j,i);
        Real u1  = cons(IM1,k,j,i) / rho;
        Real u2  = cons(IM2,k,j,i) / rho;
        Real u3  = cons(IM3,k,j,i) / rho;
        
        Real v2  = SQR(u1) + SQR(u2) + SQR(u3);
        Real ke  = 0.5 * rho * v2;
        Real pre = (cons(IEN,k,j,i) - ke) * g1;
        // Calculate change in temperature
        Real const_1 = WR.avgm / kboltz;
        Real TOld    = const_1 * pre / rho;
        Real dT      = dei(k,j,i) / rho;
        // Calculate new temperature
        Real TNew;
        TNew = std::max((TOld - dT),CCurve.tmin);
        TNew = std::min(TNew,CCurve.tmax);
        // Calculate new pressure
        Real preNew = TNew * rho / const_1;

        // Update conserved values
        cons(IEN,k,j,i) = ke + preNew / g1;
      }
    }
  }
  return;
}

//!  \fn CalcHe(Real x_e)
//  \brief Perform an integration to find h_e, the effective grain heating
//  efficiency due to electrons.
//  - Efficiency is dependent on grain charge, but in this instance it
//    is assumed that all grains are neutral.
//  Inputs:
//  - x_e: E*/kBT, where E* is the critical energy required for an
//    electron to penetrate a dust grain, D&W (1983), eq. A6
//  Returns:
//  - He: effective grain heating factor due to electons

Real CalcHe(Real x_e){
  Real intf,zpxe,phi,Ixstar,hat;
  Real f1;
  const Real x_ep1p5 = pow(x_e,1.5);
  int nzmax = 400;
  Real logzmin = -2.0;
  Real logzmax = 2.0;
  Real dlogz = (logzmax - logzmin)/(Real)(nzmax);
  Real f[nzmax+1];
  bool first = true;
  Real z[nzmax+1];
  Real expmz[nzmax+1];
  Real dz[nzmax];
  // Create z bins
  if (first){
    for (int n = 0; n <= nzmax; ++n){
      z[n] = pow(10,logzmin+(Real)(n)*dlogz);
      expmz[n] = exp(-z[n]);
    }
    for (int n = 0; n < nzmax; ++n){
      dz[n] = z[n+1] - z[n];
    }
    first = false;
  }
  // Evaluate the function at each z position, and integrate using the
  // trapezium rule
  for (int n = 0; n <= nzmax; ++n){
    zpxe = z[n] + x_e;
    f1 = pow(zpxe,1.5) - x_ep1p5;
    f[n] = zpxe*pow(f1,2.0/3.0)*expmz[n]; //integral in Dwek and Werner Eq A11 (note typo in TT2013)
  }
  intf = 0.0;
  for (int n = 0; n < nzmax; ++n){
    intf += 0.5*(f[n+1] + f[n])*dz[n];
  }
  Ixstar = 0.5*exp(-x_e)*intf;  // Eq A11 in Dwek & Werner
  return 1.0 - Ixstar;          // Eq A10 in DW81
  phi = 0.15;
  hat = exp(-phi)*(1.0 - Ixstar); // Eq A10 in DW1981
  // CURRENT STATUS - neither of these return the "correct" result - the resulting curve is different to the approximation of Eq A13 in DW1981
  return 0.0;
}

//! \fn void CalcDustLambda(rho,T,aCM,z)
//  \brief Function to calculate the cooling function, Lambda
//  for dust in a cell, given in the units (erg cm^-3 s^-1)
//  - Function is normalised using the electron and proton
//    number densities, is based on proton and electron collisions
//    with dust grains.
//  - Collisional energy is then emitted in the form
//    of black-body emission. This is not caclulated, and is assumed that
//    the time scale of emission is << timescale of simulation step
//  - Code is based off of a cooling curve function provided by Julian,
//    itself based on Dwek & Werner prescription:
//      Dwek, E., & Werner, M. W. (1981).
//      The Infrared Emission From Supernova Condensates.
//      The Astrophysical Journal, 248, 138.
//      https://doi.org/10.1086/159138
//  - A cooling curve is not used for calculating dust emission as
//    dust parameters z and a can very on a per-cell basis, unlike
//    plasma cooling, where abundances are assumed to be fairly constant
//  Input parameters:
//  - rho, density in g cm^-3
//  - T, temperature in Kelvin
//  - aCM, mean grain radius in cm, a is reserved for value in microns
//  - z, dust-to-gas mass ratio in cell, mD/mG
//  Returns:
//  - LambdaNorm, Normalised cooling function due to collisional heating
//    (erg cm^-3 s^-1) 

Real CalcDustLambda(Real rho, Real T, Real aCM, Real z) {
  // Constants from model
  const Real grainBulkDensity = bulk_dens;
  // Fundamental constants
  const Real electronMass = 9.1093837e-28;
  const Real hydrogenMass = 1.6735575e-24;  // Hydrogen atom mass (g)
  // Unit conversions
  const Real KeVtoErg   = 1.6021766e-09;
  const Real cmtoMicron = 1e4;
  // Most calculations done in microns
  Real a   = aCM * cmtoMicron; // Convert a from cm to microns
  Real kBT = kboltz * T;       // Shorthand to reduce fp ops
  // Calculate hydrogen number density
  Real nH = (10.0 / 14.0) * rho / massh; // Hydrogen number density
  // Calculate number densities dependent on nH
  Real nE   = 1.2 * nH;  // Electron number density
  Real nP   = nE;        // Proton number density
  // Caclulate density, grain mass and number density of dust
  Real rhoD = rho * z;  // Dust density (g cm^-3)
  Real volD = (4.0/3.0) * PI * CUBE(aCM);  // Dust grain volume (cm^3)
  Real massD = grainBulkDensity * volD;    // Dust grain mass (g)
  Real nD = rhoD / massD;
  // Calculate heating rate due to incident protons
  // First, calculate critical energy for penetration (KeV)
  Real eH  = 133.0 * a;
       eH *= KeVtoErg;   // Convert to ergs 
  Real hN = 1.0 - (1.0 + eH/(2.0 * kBT)) *
            std::exp(-eH / kBT);
  Real Hcoll = 1.26e-19 * SQR(a) * pow(T,1.5) * nH * hN;
  // Calculate heating rate due to incident electrons
  Real Ee  = 23.0 * pow(a, 2.0/3.0);
       Ee *= KeVtoErg;   // Calculate threshold energy
  Real xe = Ee/kBT;
  Real he = CalcHe(xe);
  Real Hel = 1.26e-19 * SQR(a) * pow(T,1.5) * nE * he /
             std::sqrt(electronMass / hydrogenMass);
  // Caclulate lambda
  Real lambda = (Hcoll + Hel) / nH;
  // Wrap up and return!
  return lambda;
}

//! \fn CoolingFunction()
//  \brief A function used to simulate the energy lost to the system through
//         radiation
//  - Cooling is simulated by the loss of cell energy by modification of the
//    cons(IEN) variable
//  - As cooling timescale may be smaller than the timestep, dt, and modifying
//    this timestep such it is < tau_cool would result in significantly more
//    computational load from hydro and other functions in the code, substepping
//    is used
//    - Energy loss integration is based on cooling time, which is updated per
//      substep, until dt has been reached
//  - Energy loss is then stored in an additional array
//    - This is done as additional checks are required to restrict cooling along
//      unresolved interfaces (see function RestrictCool) through interpolation
//    - Overall more accurate solution and prevents values from going wonky
//  - Additional checks are made in order to prevent failures in the Riemann
//    solver, as cooling can render this unstable.
//    - Main check is for div/0 errors, this seems to be enough in most cases/
//    - If you are finding the program crashes due to this, decrease the courant
//      number of your problem
//
//  ++ TO-DO ++
//  - There are a number of optimisations to be made in the dust code
//    - nE, nP and nD are calculated at every step in loop, calculate
//      initially and instead forward to function as arguments.
//  - Function assumes that WR is absolutely dominant, as such properties
//    such as wind temperature and average mass.
//    - Using an average based on wind momentum ratio might be an idea,
//      but overall it's probably not worth pursuing.

void CoolingFunction(MeshBlock *pmb, const Real dt,
                     AthenaArray<Real> &cons)
{
  AthenaArray<Real> dei(pmb->ke+1,pmb->je+1,pmb->ie+1);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Read conserved variables
        Real rho = cons(IDN,k,j,i);        // Denisty (g cm^-3)
        Real te  = cons(IEN,k,j,i);        // Total energy (ergs)
        Real u1  = cons(IM1,k,j,i) / rho;  // U_1 (cm s^-1)
        Real u2  = cons(IM2,k,j,i) / rho;  // U_2 (cm s^-1)
        Real u3  = cons(IM3,k,j,i) / rho;  // U_2 (cm s^-1)
        // Read scalars
        Real a = 0.0;
        Real z = 0.0;
        if (dust_cool) {
          z = pmb->pscalars->r(ZLOC,k,j,i);  // Dust mass fraction
          a = pmb->pscalars->r(ALOC,k,j,i);  // Dust avg. radius (cm)
        }
        // Calculate fluid velocity
        Real v2  = SQR(u1) + SQR(u2) + SQR(u3);  // Square of scalar velocity
        Real v   = std::sqrt(v2);                // Scalar velocity (cm/s)
        // Calculate average mass of particle in cell
        Real avgmass = WR.avgm;
        // Calculate average unshocked wind temperature (K)
        Real TWind   = WR.twnd; 
        // Calculate internal energy of cell
        Real ke     = 0.5 * rho * v2;  // Bulk kinetic energy of fluid (erg)
        Real ie     = te - ke;         // Internal energy of cell (erg)
        Real pre    = g1 * ie;         // Running pressure (Ba)
        Real preOld = pre;             // Original pressure (Ba)
        // Calculate Initial temperature
        Real T = (pre * avgmass) / (rho * kboltz);  // Running temperature (K)
        Real TOld = T;  // Original temperature (K)

        // Calculate cooling constants for plasma and dust cooling
        // - These are parts of the equation that are independent of pressure
        //   and temperature, so can be calculated once and used in the
        //   cooling loop
        Real edotConstGas  = SQR(rho)/SQR(massh);
        Real edotConstDust = 0.0;
        if (dust_cool && a > 0.0 && z > 0.0) {
          Real rhoD  = rho * z;
          Real massD = (4.0/3.0) * PI * CUBE(a);
          Real nD    = rhoD / massD;
          edotConstDust = (rho / massh) * nD;
        }

        // Check for most egregious error, if temperature found to be invalid
        // - Div/0 and negative pressures result in invalid values that are not
        //   checked on by Athena++ 
        // - This prevents the simulation from running on indefinitely
        //   Even after it has utterly failed
        // - This may cause the simulation to fail, if that is the case, 
        //   typically lowering Courant number fixes this!
        //   - Time step will be too large, resulting in invalid values, which
        //     then propagate to next step
        //   - It's also better than coming back to your workstation to find a
        //     bunch of gibberish, trust me!

        if (std::isnan(TOld) || std::isinf(TOld)) {
          printf("CRASH:\n");
          printf("GRID_POS %05d %05d %05d\n",i,j,k);
          printf("RHO = %.3e\n",rho);
          printf("PRE = %.3e\n",pre);

          std::stringstream msg;
          msg <<
          "### FATAL ERROR in cwb.cpp CoolingFunction" <<
          std::endl <<
          "Temperature is invalid: " <<
          TOld <<
          " Either NaN or Inf! Aborting!" <<
          std::endl;
          ATHENA_ERROR(msg);
        }
        // This should prevent wasting time in the unshocked gas
        // Forces temp to Twind, Pseudo heating effect
        if (T < 1e3) {
          T = 1e3; 
        }
        else {
          Real dtInt = 0.0;  // Internal dt, for substepping
          while (dtInt < dt) {
            Real Eint    = pre/g1;
            Real lambda  = CCurve.FindLambda(T);

            Real edotGas  = edotConstGas * lambda;
            Real edotDust = 0.0;
            if (dust_cool && z > 0.0 && a > 0.0) {
              // First, calculate the normalised cooling function due to dust
              Real lambdaDust = CalcDustLambda(rho,T,a,z);
              // Then calculate the energy lost (erg s^-1)
              edotDust = edotConstDust * lambdaDust;
            }

            Real totalCoolRate  = 0.0;
                 totalCoolRate += edotGas;
                 totalCoolRate += edotDust;
            Real t_cool = Eint / totalCoolRate;
            Real dtCool = 0.1 * std::abs(t_cool);
            if ((dtInt + dtCool) > dt) {
              dtCool = dt - dtInt;
            }
            dtInt += dtCool;
            // Calculate new temperature, update pressure
            Real dEint = -totalCoolRate * dtCool;
            Real TNew  = T * (Eint + dEint) / Eint;
            TNew = std::max(CCurve.tmin, TNew);
            TNew = std::min(TNew, CCurve.tmax);
            // Recalculate pressure based on ratio, update running temp
            pre *= TNew / T;
            T    = TNew;
            if (TNew < 1.1 * CCurve.tmin) break;
            if (TNew > 0.9 * CCurve.tmax) break;
          }
        }
        dei(k,j,i) = (TOld - T) * rho;
      }
    }
  }

  // Restrict cooling along unresolved interfaces
  RestrictCool(pmb->is,pmb->ie,
               pmb->js,pmb->je,
               pmb->ks,pmb->ke,
               pmb->pmy_mesh->ndim,
               dei,cons);
  // Adjust pressure due to cooling
  // - This is where the dei() array is actually used to change cons(IEN)
  AdjustPressureDueToCooling(pmb->is,pmb->ie,
                             pmb->js,pmb->je,
                             pmb->ks,pmb->ke,
                             dei,cons);
  return;
}

//! \fn void EvolveDust
//  \brief Function to evolve dust due to sputtering and grain growth
//  Grain sputtering using Draine & Salpeter prescription (1979) prescription
//    Ref: Draine, B. T., & Salpeter, E. E. (1979).
//         On the physics of dust grains in hot gas.
//         The Astrophysical Journal, 231, 77–94.
//         https://doi.org/10.1086/157165
//  - The dust lifetime is of the order 1e6 * (a / n) years if T > 1e6K
//    where a is the grain radius in microns, and n is the gas nucleon no.
//    density (eq. 44, fig 7.)
//  - Grain growth occurs at T < 1.5e4K, the governing equations can be found
//    at page 208 of: 
//    Spitzer Jr., L. (2008).
//    Physical Processes in the Interstellar Medium.
//    - Grain growth occurs at a rate:
//      da/dt = (w_a rho_a eps_a)/(4 rho_s)
//      Where w_a is the grain RMS velocity, rho_a is the interstellar mass
//      density of atoms, rho_s is the density of the grains and eps is
//      the sticking probability of atoms onto grains

void EvolveDust(MeshBlock *pmb, const Real dt, AthenaArray<Real> &cons) {
  // Problem dependent dust variables
  const Real grainDensity = bulk_dens;  // Grain bulk density (g cm^-3)
  const Real minGrainRadius = a_min;
  const Real minDustToGasMassRatio = z_min;
  // Hard-coded dust variables
  const Real A            = 12.0;  // Molecular weight of carbon grains
  const Real eps_a        = 0.1;   // Sticking probability

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Import and calculate gas parameters of cell
        Real rho = cons(IDN,k,j,i);  // Gas density (g cm^-3)
        Real te  = cons(IEN,k,j,i);  // Gas total energy (ergs)
        // Converved momentum must be converted into velocities
        Real u1  = cons(IM1,k,j,i) / rho;  // Gas velocity x direction (cm s^-1)
        Real u2  = cons(IM2,k,j,i) / rho;  // Gas velocity y direction (cm s^-1)
        Real u3  = cons(IM3,k,j,i) / rho;  // Gas velocity z direction (cm s^-1)
        // Import scalars from pscalars
        Real col = pmb->pscalars->s(CLOC,k,j,i) / rho;  // Wind colour
        Real a   = pmb->pscalars->s(ALOC,k,j,i) / rho;  // Grain radius (cm)
        Real z   = pmb->pscalars->s(ZLOC,k,j,i) / rho;  // Dust mass fraction 
        // Calculate temperature
        Real v2   = SQR(u1) + SQR(u2) + SQR(u3);
        Real ke   = 0.5 * rho * v2;
        Real pre  = (te - ke)*g1;
        Real temp = WR.avgm * pre / (rho * kboltz);
        // Replace a with minimum nucleation size and z if needed
        a = std::max(minGrainRadius,a);
        z = std::max(minDustToGasMassRatio,z);
        // Estimate gas number density (Hydrogen and ions)
        Real nH = rho * (10.0/14.0) / massh;  // Hyrogen number density (cm^-1)
        Real nTot = 1.1 * nH;                 // Total number density (cm^-1)
        // Calculate additional dust properties
        Real a2    = SQR(a); 
        Real a3    = CUBE(a); 
        Real massD = (4.0 / 3.0) * PI * a3 * grainDensity;
        Real rhoD  = rho * z;  // Dust density (g cm^-3)
        Real nD    = rhoD / massD;  // Dust number density (cm^-3)
        // Check to see if dust should grow or shrink in cell
        Real rhoD_dot = 0.0;
        Real dadt     = 0.0;

        if (temp > 1.0e6 && z > 0.0) {
          // Dust destruction through thermal sputtering occurs
          Real tau_D    = 3.156e17 * a / nTot;
          Real dadt     = -a / tau_D;
          Real rhoD_dot = -1.33e-17 * grainDensity * a2 * nTot * nD;
        }
        else if (temp < 1.5e4) {
          // Dust growth occurs
          Real wa   = std::sqrt(3.0*kboltz*temp/(A*massh));
          Real dadt = 0.25 * eps_a * rho * wa / grainDensity;
          rhoD_dot  = 4.0 * PI * a2 * grainDensity * nD * dadt;
        }
        if (rhoD_dot != 0.0) {
          Real dRhoD   = rhoD_dot * dt;  // Integrate to find total change
          // Calculate new dust density
          Real minRhoD = minDustToGasMassRatio * rho;
          Real rhoDNew = std::max(minRhoD, rhoD + dRhoD);
          // Calculate new gas density
          Real rhoNew  = rho + (rhoD - rhoDNew);
          // Calculate new dust-to-gass mass ratio
          Real zNew = rhoDNew / rhoNew;
          // Calculate new grain density
          Real da   = dadt * dt;
          Real aNew = (a + da);
          // Rewrite conserved arrays
          cons(IDN,k,j,i) = rhoNew;
          // Rewrite scalars
          // Update the conserved wind colour, as gas density has changed
          pmb->pscalars->s(CLOC,k,j,i) = col * rhoNew;
          // Update dust scalars
          pmb->pscalars->r(ZLOC,k,j,i) = zNew;
          pmb->pscalars->r(ALOC,k,j,i) = aNew;
          pmb->pscalars->s(ZLOC,k,j,i) = zNew * rhoNew;
          pmb->pscalars->s(ALOC,k,j,i) = aNew * rhoNew;
        }
      }
    }
  }

  return;  
}


//! \fn void CWBExplicitSourceFunction
//  \brief Function containing smaller functions to be enrolled into meshblock class
//  As only one explicit source function can be enrolled into Athena, and
//  using a single long source function would be clunky, this function
//  calls other smaller functions that handle:
//  - Cooling due to plasma and dust emission
//  - Dust growth and desctruction
//  This does mean that these events do not happen simultaneously, but they are
//  not particularly dependent on one another.
//  Initially wind remap was performed here, 

void CWBExplicitSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                               const AthenaArray<Real> &prim,
                               const AthenaArray<Real> &bcc,
                               AthenaArray<Real> &cons) {
  // Run cooling function
  if (cool) {CoolingFunction(pmb,dt,cons);}
  // Run dust evolution function
  if (dust) {EvolveDust(pmb,dt,cons);}
  // Finish up
  return;
}

int CWBRefinementCondition(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  AthenaArray<Real> &r = pmb->pscalars->r;

  Real dx = pmb->pcoord->dx1f(0);
  Real dy = pmb->pcoord->dx2f(0);
  Real dz = pmb->pcoord->dx3f(0);

  // Refine around stagnation point
  
  Real epsilon = 0.33;
  Real avgmass = 1.0e-24;
  Real boltzman = 1.380658e-16;
  Real pi = 2.0*asin(1.0);

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

  Real dxx = WR.x[0] - OB.x[0];
  Real dyy = WR.x[1] - OB.x[1];
  Real theta = std::atan2(dyy,dxx);

  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
    stagx = OB.x[0] + (rob * std::cos(theta));
    stagy = OB.x[1] + (rob * std::sin(theta));
    stagz = 0.0;
  }

  for(int k=pmb->ks-1; k<=pmb->ke+1; k++) {
    Real dz = (pmb->pcoord->x3v(k+1) - pmb->pcoord->x3v(k-1));
    for(int j=pmb->js-1; j<=pmb->je+1; j++) {
      Real dy = (pmb->pcoord->x2v(j+1) - pmb->pcoord->x2v(j-1));
      for(int i=pmb->is-1; i<=pmb->ie+1; i++) {
        Real dx = (pmb->pcoord->x1v(i+1) - pmb->pcoord->x1v(i-1));
	      // Refine on divergence condition
	      Real divV = 0.0;
        //if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
        Real dudx = (w(IVX,k,j,i+1)-w(IVX,k,j,i-1))/dx;
        Real dvdy = (w(IVY,k,j+1,i)-w(IVY,k,j-1,i))/dy;
        Real dwdz = (w(IVZ,k+1,j,i)-w(IVZ,k-1,j,i))/dz;
        divV = dudx + dvdy + dwdz;
	 
        if (divV < 0.0){ // potentially refine
          // Check to see if not too far away from the stagnation point.
          //**** CHECK ALL OF THIS!!!!****
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
      }
    }
  }




  // Refine based on proximity to remap zone
  // Check to see if cell is within 10 cells of remap zone
  bool derefineWind = true;
  for (auto Star : Stars) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      Real zc = pmb->pcoord->x3v(k) - Star->x[2];
      for (int j=pmb->js; j<=pmb->je; ++j) {
        Real yc = pmb->pcoord->x2v(j) - Star->x[1];
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real xc  = pmb->pcoord->x1v(i) - Star->x[0];
          Real rad = std::sqrt(SQR(xc) + SQR(yc) + SQR(zc));
          int ri = int(rad / dx);
          if      (ri < 10) {return 1;}
          else if (ri < 20) {derefineWind = false;}
        }
      }
    }
  }
  if (derefineWind) {return -1;}
  // All checks passed, no change
  return 0;
}




// int RefineRemaps(MeshBlock *pmb) {
//   return ;
// } 
//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Get boolean values from problem file
  orbit = pin->GetBoolean("problem","orbit");
  cool  = pin->GetBoolean("problem","cool");
  dust  = pin->GetBoolean("problem","dust");
  if (dust) {dust_cool = pin->GetBoolean("problem","dust_cool");}

  // Initialise stars
  WR.Init(pin,"wr",1);
  OB.Init(pin,"ob",0);
  // Initialise star vector
  Stars.push_back(&WR);
  Stars.push_back(&OB);
  // Initialise cooling curve
  if (cool) {
    std::string coolcurvefilename = pin->GetString("problem","coolcurve");
    CCurve.Init(coolcurvefilename);
  }

  // Initialise dust
  a_min = pin->GetReal("problem","a_min");
  z_min = pin->GetReal("problem","z_min");

  // Get orbit parameters from problem file
  period       = pin->GetReal("problem","period");
  ecc          = pin->GetReal("problem","ecc");
  phase_offset = pin->GetReal("problem","phase_offset");

  // Get dust parameters
  bulk_dens = pin->GetReal("problem","bulk_dens");
  stick_eff = pin->GetReal("problem","stick_eff");
  nuc_temp  = pin->GetReal("problem","nuc_temp");

  // Enroll functions based on boolean values
  // Refinement functions
  if (adaptive) {
    EnrollUserRefinementCondition(CWBRefinementCondition);
  }
  // Source functions
  // EnrollUserExplicitSourceFunction(RemapWinds);
  EnrollUserExplicitSourceFunction(CWBExplicitSourceFunction);
  
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void __attribute__((weak)) Mesh::UserWorkInLoop() {
  // do nothing
  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief Function called after main loop is finished for user-defined work.
//========================================================================================

void __attribute__((weak)) Mesh::UserWorkAfterLoop(ParameterInput *pin) {
  // do nothing
  return;
}

// 4x members of MeshBlock class:

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//  used to initialize variables which are global to other functions in this file.
//  Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {

  // do nothing
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief CWB problem generator
//========================================================================================

// void MeshBlock::ProblemGenerator(ParameterInput *pin) {

// }
 

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real d_amb = pin->GetReal("problem","d_amb");
  Real p_amb = pin->GetReal("problem","p_amb");

  eta = (OB.mdot*OB.vinf)/(WR.mdot*WR.vinf);

  // If orbits are enabled, recalculate positions and parameters
  if (orbit) {
    for (auto Star : Stars) {Star->x[2] = 0.0;}
    CalcOrbit(0.0); // Find orbits at t = 0
  }


  


  // Check to see if in a valid co-ordinate system, otherwise, exit with errors
  // Coordinate described with integer value from here on out
  // This is done as strcmp() is quite costly, especially within loops
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {coord = CART;}
  else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {coord = CYL;}
  else {
    std::stringstream msg;
    msg << "### FATAL ERROR in cwb.cpp ProblemGenerator" << std::endl
        << "Unrecognized COORDINATE_SYSTEM=" << COORDINATE_SYSTEM << std::endl;
    ATHENA_ERROR(msg);
  }
  // This makes the assumption that WR is the primary star!!!
  Real dx = WR.x[0] - OB.x[0];
  Real dy = WR.x[1] - OB.x[1];
  Real theta = std::atan2(dy,dx);
  Real stagx = OB.x[0] + (rob * std::cos(theta));
  Real stagy = OB.x[1] + (rob * std::sin(theta));
  // Run through grid in meshblock, finding the side the cell is on,
  // then rewriting properties based on single wind outflow from star
  // the cell is associated with
  for (int k=ks; k<=ke; k++) {
    Real z  = pcoord->x3v(k);
    Real z2 = SQR(z);
    for (int j=js; j<=je; j++) {
      Real y  = pcoord->x2v(j);
      Real y2 = SQR(y); 
      for (int i=is; i<=ie; i++) {
        Real x  = pcoord->x1v(i);
        Real x2 = SQR(x);
        // Translate coordinates to centre around stagnation point
        Real xp = x - stagx;
        Real yp = y - stagy;
        // Rotate coordinate system such that stars run along 
        Real xpp = (xp * std::cos(theta)) + (yp * std::sin(theta));
        Real ypp = (yp * std::cos(theta)) - (xp * std::sin(theta));
        // In rotated coordinate system, determine which side of the
        // stagnation point the cell is on
        int sid;
        if      (xpp >= 0) {sid = 0;} // Wind is on WR side
        else if (xpp <= 0) {sid = 1;} // Wind is on OB side
        // Find the distance relative to the star in question
        Real xc  = x - Stars[sid]->x[0];
        Real yc  = y - Stars[sid]->x[1];
        Real zc  = z - Stars[sid]->x[2];
        Real xc2 = SQR(xc);
        Real yc2 = SQR(yc);
        Real zc2 = SQR(zc);
        Real r2  = xc2 + yc2 + zc2;       // Square of distance from star
        Real r   = std::sqrt(r2);         // Distance from star
        Real xy  = std::sqrt(xc2 + yc2);  // xy axis distance from star
        // Calculate angles, names should be self-explanitory
        Real sinphi = xy/r;
        Real cosphi = zc/r;
        Real costht = xc/xy;
        Real sintht = yc/xy;
        // Copy star properties from class
        Real mdot = Stars[sid]->mdot_cgs;  // Mass loss rate (g s^-1)
        Real vinf = Stars[sid]->vinf;      // Terminal velocity (cm s^-1)
        Real twnd = Stars[sid]->twnd;      // Wind temperature (K)
        // Calculate cell properties based on distance from star
        // Density, pressure, kinetic energy and velocities in all dimensions
        Real rho = Stars[sid]->mdot_cgs / (4.0 * PI * r2 * WR.vinf);
        Real pre = (rho/WR.avgm) * kboltz * WR.twnd;
        Real ke  = 0.5 * rho * SQR(vinf);
        Real u1  = vinf * sinphi * costht;
        Real u2  = vinf * sinphi * sintht;
        Real u3  = vinf * cosphi;
        // Modify properties of cell
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * u1;
        phydro->u(IM2,k,j,i) = rho * u2;
        phydro->u(IM3,k,j,i) = rho * u3;
        phydro->u(IEN,k,j,i) = pre/g1 + ke;

        if (NSCALARS > 0) {
          pscalars->s(0,k,j,i) = Stars[sid]->scal * rho;
          pscalars->r(0,k,j,i) = Stars[sid]->scal;
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.

//  Function handles:
//  - Updating "remap zone" positions based on Keplerian "track"
//  - Remapping a supersonic wind within a small region called a "remap zone"
//    this is done in order to emulate a star in that area
//  Most other functions are handled using a source function, however 
//  remapping is performed after each timestep in order to "force" scalars
//  to a specific value
//========================================================================================

void MeshBlock::UserWorkInLoop() {

  // Update star positions based on Keplerian orbit
  if (orbit) {CalcOrbit(pmy_mesh->time);}

  // Map wind onto stars, this requires a lot of mess to spin out
  // into a function, so I'll just include it here
  // First, calculate the remap radius in CGS units, rather than cell widths
  Real remap_radius = Real(remap) * pcoord->dx1v(0);

  // As a meshblock is refined down to the appropriate level it
  // needs to be in an adaptive mesh, this means that the mesh
  // itself is ordered, therefore rather than looping over
  // the entire mesh, only a cubic region around the remap zone
  // is looped through.
  // i/j/kstar: Closest index on axis to star position
  // i/j/kstl:  "Lower" array to access
  // i/j/kstu:  "Upper" array to access

  for (auto Star : Stars) {
    Star->ncellsremapped = 0;

    int istar = int((Star->x[0] - pcoord->x1f(0))/pcoord->dx1f(0));
    int jstar = int((Star->x[1] - pcoord->x2f(0))/pcoord->dx2f(0));
    int kstar = int((Star->x[2] - pcoord->x3f(0))/pcoord->dx3f(0));
    int istl  = std::max(is,istar-remap-2);
    int jstl  = std::max(js,jstar-remap-2);
    int kstl  = std::max(ks,kstar-remap-2);
    int istu  = std::min(ie,istar+remap+2);
    int jstu  = std::min(je,jstar+remap+2);
    int kstu  = std::min(ke,kstar+remap+2);
    // In the case of a cylindrical coordinate system, phi axis modified
    if (coord == CYL) {
      jstl = js;
      jstu = je;
    }
    // Loop through cubic region around star
    for (int k=kstl; k<=kstu; ++k) {
      Real zc  = pcoord->x3v(k) - Star->x[2];
      Real zc2 = SQR(zc);
      for (int j=jstl; j<=jstu; ++j) {
        Real yc  = pcoord->x2v(j) - Star->x[1];
        Real yc2 = SQR(yc);
        for (int i=istl; i<=istu; ++i) {
          Real xc  = pcoord->x1v(i) - Star->x[0];
          Real xc2 = SQR(xc);
          Real r2  = xc2 + yc2 + zc2;
          Real r   = std::sqrt(r2);
          // Check to see if cell is within remap zone
          if (r < remap_radius) {
            r  = std::max(1.0,r);
            r2 = std::max(1.0,r2);
            
            Real xy = std::sqrt(xc2 + yc2);
            Real sinphi = xy/r;
            Real cosphi = zc/r;
            
            Real rho = Star->mdot_cgs / (4.0 * PI * r2 * Star->vinf);
            Real pre = (rho / Star->avgm) * kboltz * Star->twnd;
            Real ke  = 0.5 * rho * SQR(Star->vinf);

            // Calculate velocity components based on an outflowing wind of
            // Terminal velocity vinf
            // Cases for cartesian and cylindrical coordinates
            Real u1, u2, u3;
            if (coord == CART) {
              Real costht = xc/xy;
              Real sintht = yc/xy;
              u1 = Star->vinf * sinphi * costht;
              u2 = Star->vinf * sinphi * sintht;
              u3 = Star->vinf * cosphi;
            }
            else if (coord == CYL) {
              u1 = Star->vinf * sinphi;
              u2 = 0.0;
              u3 = Star->vinf * cosphi;
            }

            phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = u1 * rho;
            phydro->u(IM2,k,j,i) = u2 * rho;
            phydro->u(IM3,k,j,i) = u3 * rho;
            phydro->u(IEN,k,j,i) = pre/g1 + ke;

            // Use 1st scalar as marker for each star
            if (NSCALARS > 0) {
              pscalars->s(CLOC,k,j,i) = Star->scal * rho;
              if (dust) {
                pscalars->s(ZLOC,k,j,i) = z_min * rho;
                pscalars->s(ALOC,k,j,i) = a_min * rho;
              }
            }

            Star->ncellsremapped++;
          }
        }
      }
    }
  }

  // if (WR.ncellsremapped > 0 || OB.ncellsremapped > 0) {
  //   printf("Cells remapped:\n");
  //   printf(" + WR = %d\n",WR.ncellsremapped);
  //   printf(" + OB = %d\n",OB.ncellsremapped);
  // }

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief Function called before generating output files
//========================================================================================

void __attribute__((weak)) MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // do nothing
  return;
}
