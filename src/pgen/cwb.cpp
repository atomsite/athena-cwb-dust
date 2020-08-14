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

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"

#include "../hydro/hydro.hpp"
#include "../coordinates/coordinates.hpp"


#include "cwb.hpp"

namespace {
  // Create star objects, see cwb.hpp
  Star WR;  // Object describing properties of Wolf-Rayet star
  Star OB;  // Object describing properties of OB star

  Real a_min;
  Real z_min;
  Real bulk_dens;
  Real stick_eff;
  Real nuc_temp;

  bool orbit;
  bool cool;


  // Constants
  Real kboltz = 1.380649e-16;
}

bool IsInRemap(Real x, Real y, Real z, Star Star) {
  Real x0 = Star.x1;
  Real y0 = Star.x2;
  Real z0 = Star.x3;
  Real rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
  if (rad <= Star.remap_cgs) {return true;}
  return false;
}

// Map winds from both stars, needs some cleanup

void RemapWind(MeshBlock *pmb, const Real time, const Real dt,
               const AthenaArray<Real> &prim,
               const AthenaArray<Real> &bcc,
               AthenaArray<Real> &cons) {
  
  int is = pmb->is, ie = pmb->ie;
  int js = pmb->js, je = pmb->je;
  int ks = pmb->ks, ke = pmb->ke;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x   = pmb->pcoord->x1v(i);
        Real y   = pmb->pcoord->x2v(j);
        Real z   = pmb->pcoord->x3v(k);
        if (IsInRemap(x,y,z,WR) == true) {
          Real rhodot   = WR.mdot_cgs / WR.remap_vol;
          Real mech_lum = 0.5 * rhodot * SQR(WR.vinf);
          cons(IDN,k,j,i) += rhodot * dt;
          cons(IEN,k,j,i) += mech_lum * dt;
        }
        else if (IsInRemap(x,y,z,OB) == true) {
          Real rhodot   = OB.mdot_cgs / OB.remap_vol;
          Real mech_lum = 0.5 * rhodot * SQR(OB.vinf);
          cons(IDN,k,j,i) += rhodot * dt;
          cons(IEN,k,j,i) += mech_lum * dt;
        }
      }
    }
  }
  return;
}

// Not used right now, but a basic cooling function from another code, needs to be revised to fit plasma cooling curves

void CoolingFunction(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Real g = 1.666667;
  Real tau = 0.01;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real temp = (g-1.0)*prim(IEN,k,j,i)/prim(IDN,k,j,i);
        cons(IEN,k,j,i) -= dt*prim(IDN,k,j,i)*(temp - 10.0)/tau/(g-1.0);
      }
    }
  }
  return;
}

int RefineShocks(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  for (int k=pmb->ks; k<=pmb->ke ; k++) {
    for(int j=pmb->js; j<=pmb->je; j++) {
      for(int i=pmb->is; i<=pmb->ie; i++) {
        Real epsr= (std::abs(w(IDN,k,j,i+1)-2.0*w(IDN,k,j,i)+w(IDN,k,j,i-1))
                  +std::abs(w(IDN,k,j+1,i)-2.0*w(IDN,k,j,i)+w(IDN,k,j-1,i)))/w(IDN,k,j,i);
        Real epsp= (std::abs(w(IEN,k,j,i+1)-2.0*w(IEN,k,j,i)+w(IEN,k,j,i-1))
                  +std::abs(w(IEN,k,j+1,i)-2.0*w(IEN,k,j,i)+w(IEN,k,j-1,i)))/w(IEN,k,j,i);
        Real eps = std::max(epsr, epsp);
        maxeps = std::max(maxeps, eps);
      }
    }
  }
  if(maxeps > 0.01) return 1;
  if(maxeps < 0.005) return -1;
  return 0;
}

// int RefineRemaps(MeshBlock *pmb) {
//   return ;
// } 


// 3x members of Mesh class:

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in Mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Initialise stars
  WR.Init(pin,"wr");
  OB.Init(pin,"ob");

  // Initialise dust
  a_min = pin->GetReal("problem","a_min");
  z_min = pin->GetReal("problem","z_min");

  // Get boolean values from problem file
  orbit = pin->GetBoolean("problem","orbit");
  cool  = pin->GetBoolean("problem","cool");

  // Get dust parameters
  bulk_dens = pin->GetReal("problem","bulk_dens");
  stick_eff = pin->GetReal("problem","stick_eff");
  nuc_temp  = pin->GetReal("problem","nuc_temp");

  // Enroll functions based on boolean values
  // Refinement functions
  if (adaptive) {
    EnrollUserRefinementCondition(RefineShocks);
  }
  // Source functions
  EnrollUserExplicitSourceFunction(RemapWind);
  // if (cool) {EnrollUserExplicitSourceFunction(CoolingFunction);}

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
 

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real d_amb = pin->GetReal("problem","d_amb");
  Real p_amb = pin->GetReal("problem","p_amb");

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {

        Real x0 = WR.x1;
        Real y0 = WR.x2;
        Real z0 = WR.x3;
        Real rad;

        if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
          Real x = pcoord->x1v(i);
          Real y = pcoord->x2v(j);
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
          Real x = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j));
          Real z = pcoord->x3v(k);
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        } else { // if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)
          Real x = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::cos(pcoord->x3v(k));
          Real y = pcoord->x1v(i)*std::sin(pcoord->x2v(j))*std::sin(pcoord->x3v(k));
          Real z = pcoord->x1v(i)*std::cos(pcoord->x2v(j));
          rad = std::sqrt(SQR(x - x0) + SQR(y - y0) + SQR(z - z0));
        }

        Real rho_sw = WR.mdot_cgs / (4.0*PI*rad*rad*WR.vinf);
        Real pg_sw  = (rho_sw * kboltz * WR.twnd) / (WR.avgm);

        phydro->u(IDN,k,j,i) = rho_sw * 0.0001;
        phydro->u(IEN,k,j,i) = pg_sw * 0.0001;
      }
    }
  }
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  // do nothing
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
