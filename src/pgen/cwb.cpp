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
}


void RemapWind(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {
  return;
}

void CoolingFunction(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons) {
  return;
}

int RefineShocks(MeshBlock *pmb)
{
  AthenaArray<Real> &w = pmb->phydro->w;
  Real maxeps=0.0;
  int k=pmb->ks;
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
  if(maxeps > 0.01) return 1;
  if(maxeps < 0.005) return -1;
  return 0;
}

int RefineRemaps(MeshBlock *pmb) {
  AthenaArray<Real> &w = pmb->phydro->w;
} 


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
  bulk_dens = 3.0;
  stick_eff = 0.01;
  nuc_temp  = 1e4;

  if(adaptive) {
    EnrollUserRefinementCondition(RefineShocks);
    EnrollUserRefinementCondition(RefineRemaps);
  }

  EnrollUserExplicitSourceFunction(RemapWind);
  EnrollUserExplicitSourceFunction(CoolingFunction);

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
  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  // do nothing
  std::cout << WR.mdot << "\n";
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
