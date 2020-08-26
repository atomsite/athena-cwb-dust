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
#include "../scalars/scalars.hpp"


#include "cwb.hpp"

namespace {
  // Define objects
  // Create star objects, see cwb.hpp
  Star WR;  // Object describing properties of Wolf-Rayet star
  Star OB;  // Object describing properties of OB star
  // Vector containing pointers to both stars
  std::vector<Star*> Stars;

  // Create cooling curve class
  CoolCurve CCurve;
  // Star orbital parameters
  Real dsep;
  Real eta;
  Real rob;
  Real rwr;
  Real phase_offset;
  Real period;
  Real ecc;

  // Dust variables
  Real a_min;
  Real z_min;
  Real bulk_dens;
  Real stick_eff;
  Real nuc_temp;

  // Simulation parameters
  int  remap = 3;     // Remap radius
  bool orbit = false;  // Turn on orbital dynamics
  bool cool  = false;  // Turn on cooling

  // Simulation variables
  int coord; // Coordinate system, 0 for cartesian, 1 for cylindrical

  // Constants, all values are in CGS
  Real massh  = 1.6735575e-21;  // Hydrogen atom mass (g)
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

void RemapWinds(MeshBlock *pmb, const Real time, const Real dt,
                const AthenaArray<Real> &prim,
                const AthenaArray<Real> &bcc,
                AthenaArray<Real> &cons) {

  // First, calculate the remap radius in CGS units, rather than cell widths
  Real remap_radius = Real(remap) * pmb->pcoord->dx1v(0);
  
  // As a meshblock is refined down to the appropriate level it
  // needs to be in an adaptive mesh, this means that the mesh
  // itself is ordered, therefore rather than looping over
  // the entire mesh, only a cubic region around the remap zone
  // is looped through.
  // i/j/kstar: Closest index on axis to star position
  // i/j/kstl:  "Lower" array to access
  // i/j/kstu:  "Upper" array to access

  for (auto Star : Stars) {
    int istar = int((Star->x[0] - pmb->pcoord->x1f(0))/pmb->pcoord->dx1f(0));
    int jstar = int((Star->x[1] - pmb->pcoord->x2f(0))/pmb->pcoord->dx2f(0));
    int kstar = int((Star->x[2] - pmb->pcoord->x3f(0))/pmb->pcoord->dx3f(0));
    int istl = std::max(pmb->is,istar-remap-2);
    int jstl = std::max(pmb->js,jstar-remap-2);
    int kstl = std::max(pmb->ks,kstar-remap-2);
    int istu = std::min(pmb->ie,istar+remap+2);
    int jstu = std::min(pmb->je,jstar+remap+2);
    int kstu = std::min(pmb->ke,kstar+remap+2);
    // In the case of a cylindrical coordinate system, phi axis modified
    if (coord == CYL) {
      jstl = pmb->js;
      jstu = pmb->je;
    }
    // Loop through cubic region around star
    for (int k=kstl; k<=kstu; ++k) {
      Real zc  = pmb->pcoord->x3v(k) - Star->x[2];
      Real zc2 = SQR(zc);
      for (int j=jstl; j<=jstu; ++j) {
        Real yc  = pmb->pcoord->x2v(j) - Star->x[1];
        Real yc2 = SQR(yc);
        for (int i=istl; i<=istu; ++i) {
          Real xc  = pmb->pcoord->x1v(i) - Star->x[0];
          Real xc2 = SQR(xc);
          Real r2  = xc2 + yc2 + zc2;
          Real r   = std::sqrt(r2);
          // Check to see if cell is within remap zone
          if (r < remap_radius) {
            Real xy = xy = std::sqrt(xc2 + yc2);
            Real sinphi = xy/r;
            Real cosphi = zc/r;
            
            Real rho = Star->mdot_cgs / (4.0 * PI * r2 * Star->vinf);
            Real u1, u2, u3;
            Real pre = (rho / Star->avgm) * kboltz * Star->twnd;
            Real ke  = 0.5 * rho * SQR(Star->vinf);
            if (coord == CART) {
              Real costht = xc/xy;
              Real sintht = yc/xy;
              u1 = Star->vinf * sinphi * costht;
              u2 = Star->vinf * sinphi * sintht;
              u3 = Star->vinf * cosphi;
            }
            else if (coord = CYL) {
              u1 = Star->vinf * sinphi;
              u2 = 0.0;
              u3 = Star->vinf * cosphi;
            }
            // Rewrite conserved variables in cell
            cons(IDN,k,j,i) = rho;
            cons(IM1,k,j,i) = rho * u1;
            cons(IM2,k,j,i) = rho * u2;
            cons(IM3,k,j,i) = rho * u3;
            cons(IEN,k,j,i) = pre/g1 + ke;
            // Use 1st scalar as marker for each star
            if (NSCALARS > 0) {
              pmb->pscalars->s(0,k,j,i) = Star->scal * rho;
              pmb->pscalars->r(0,k,j,i) = Star->scal;
            }
          }
        }
      }
    }
  }

  return;
}

void CoolingFunction(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  for (auto Star : Stars) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real rho = prim(IDN,k,j,i);
          Real T   = (prim(IEN,k,j,i) * Star->avgm)/(prim(IDN,k,j,i) * kboltz);
          Real wnd = pmb->pscalars->r(0,k,j,i);
          // Calculate temperature loss due to plasma emission
          Real Lambda = CCurve.FindLambda(T);
          Real eConst = SQR(rho) / SQR(massh);
          Real eLoss  = eConst * Lambda * dt;
          // Remove energy from convserved energy
          cons(IEN,k,j,i) -= eLoss * 1e2;
        }
      }
    }
  }
}

//! \fn void CWBExplicitSourceFunction
//  \brief Function containing smaller functions to be enrolled into meshblock class
//  As only one explicit source function can be enrolled into Athena, and
//  using a single long source function would be clunky, this funciton
//  is used to 

void CWBExplicitSourceFunction(MeshBlock *pmb, const Real time, const Real dt,
                               const AthenaArray<Real> &prim,
                               const AthenaArray<Real> &bcc,
                               AthenaArray<Real> &cons) {
  // Run cooling function
  CoolingFunction(pmb,time,dt,prim,bcc,cons);
  // Map wind onto new regions
  RemapWinds(pmb,time,dt,prim,bcc,cons);
  // Finish up
  return;
}

int CWBRefinementCondition(MeshBlock *pmb) {
  Real dx = pmb->pcoord->dx1f(0);
  Real dy = pmb->pcoord->dx2f(0);
  Real dz = pmb->pcoord->dx3f(0);

  // Refine around stagnation point


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
        phydro->u(IEN,k,j,i) = (pre/g1) + ke;
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  if (orbit) {CalcOrbit(pmy_mesh->time);}
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
