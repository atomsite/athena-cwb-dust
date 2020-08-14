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
//========================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
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
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"   // JMP: needed for scalars

Real gmma1;
Real mass1, mass2;
Real mdot1, mdot2, vinf1, vinf2;
Real xpos1, xpos2, ypos1, ypos2, zpos1, zpos2;
Real xvel1, xvel2, yvel1, yvel2, zvel1, zvel2;
Real eta, rob;
Real dsep;    // current stellar separation
Real dsep0;   // initial stellar separation
Real period;  // orbit period (s)  
Real phaseoff;// phase offset of orbit (from periastron) 
Real ecc;     // orbit eccentricity

int G0_level = 0;

// JMP prototypes
void orbit_calc(Real t);
void RemapWinds(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
void FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);
int RefinementCondition(MeshBlock *pmb);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  mdot1 = pin->GetReal("problem","mdot1");
  mdot2 = pin->GetReal("problem","mdot2");
  vinf1 = pin->GetReal("problem","vinf1");
  vinf2 = pin->GetReal("problem","vinf2");
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
  
  Real Msol = 1.9891e33;
  Real yr = 3.15569e7;

  mdot1 *= Msol/yr;
  mdot2 *= Msol/yr;

  // Which order are these guaranteed in? It might be that the second replaces the first, so that
  // it is not even possible to have two functions...
  EnrollUserExplicitSourceFunction(RemapWinds);
  //EnrollUserExplicitSourceFunction(FixTemperature);

  if (adaptive==true)
      EnrollUserRefinementCondition(RefinementCondition);  

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Problem Generator for the CWB test
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gmma  = peos->GetGamma();
  gmma1 = gmma - 1.0;

  Real Twind = 1.0e4;         // K
  Real avgmass = 1.0e-24;     // g
  Real boltzman = 1.380658e-16;
  Real pi = 2.0*asin(1.0);
  
  dsep = std::sqrt(std::pow(xpos1 - xpos2,2) + std::pow(ypos1 - ypos2,2) + std::pow(zpos1 - zpos2,2));
  eta = mdot2*vinf2/(mdot1*vinf1);                   // wind mtm ratio
  rob = (std::sqrt(eta)/(1.0 + std::sqrt(eta)))*dsep;//distance of stagnation point from star 1 (distance from star 0 is rwr)

  orbit_calc(pmy_mesh->time);
  
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
	    Real u1 = vinf1*sinphi*costhta;
	    Real u2 = vinf1*sinphi*sinthta;
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
	    Real u1 = vinf2*sinphi*costhta;
	    Real u2 = vinf2*sinphi*sinthta;
	    Real u3 = vinf2*cosphi;
	    Real pre = (rho/avgmass)*boltzman*Twind;
	    phydro->u(IDN,k,j,i) = rho;
            phydro->u(IM1,k,j,i) = rho*u1;
            phydro->u(IM2,k,j,i) = rho*u2;
            phydro->u(IM3,k,j,i) = rho*u3;
            phydro->u(IEN,k,j,i) = pre/gmma1 + 0.5*rho*vinf2*vinf2;
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
	  }	  
        }
	
	if (NSCALARS > 0) {
          pscalars->s(0,k,j,i) = 0.0;
        }

      }
    }
  }
  
  return;
}



//========================================================================================
//! \fn void MeshBlock::RemapWinds(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Remap winds
//  NOTE: JMP 11/08/20 - The wind passive scalar for the primary stays within 3% of 1.0.
//                       Therefore use 0.5 to distinguish between primary and secondary wind material.
//========================================================================================
void RemapWinds(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){

  Real g = pmb->peos->GetGamma();
  gmma1 = g - 1.0;

  Real Twind = 1.0e4;         // K
  Real avgmass = 1.0e-24;     // g
  Real boltzman = 1.380658e-16;
  Real pi = 2.0*asin(1.0);
  
  int remapi = 10;
  Real remap_radius = Real(remapi)*pmb->pcoord->dx1v(0); //1.0e11;

  // The following doesn't remap the winds properly...
  //Real remap_radius = 3.75e12;
  //int remapi = int(remap_radius/10.0);


  orbit_calc(time);

  
  Real xpos[2]={xpos1,xpos2};
  Real ypos[2]={ypos1,ypos2};
  Real zpos[2]={zpos1,zpos2};
  Real mdot[2]={mdot1,mdot2};
  Real vinf[2]={vinf1,vinf2};
  Real scalar[2]={1.0,0.0};
  
  for (int nw = 0; nw < 2; ++nw){ // Loop over each wind
    int istar = int((xpos[nw] - pmb->pcoord->x1f(0))/pmb->pcoord->dx1f(0));
    int jstar = int((ypos[nw] - pmb->pcoord->x2f(0))/pmb->pcoord->dx2f(0));
    int kstar = int((zpos[nw] - pmb->pcoord->x3f(0))/pmb->pcoord->dx3f(0));
    int istl = std::max(pmb->is,istar-remapi-2);
    int jstl = std::max(pmb->js,jstar-remapi-2);
    int kstl = std::max(pmb->ks,kstar-remapi-2);
    int istu = std::min(pmb->ie,istar+remapi+2);
    int jstu = std::min(pmb->je,jstar+remapi+2);
    int kstu = std::min(pmb->ke,kstar+remapi+2);
    if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
      jstl = pmb->js; jstu = pmb->je;
    }
    for (int k=kstl; k<=kstu; k++) {
      Real zc = pmb->pcoord->x3v(k) - zpos[nw];
      Real zc2 = zc*zc;
      for (int j=jstl; j<=jstu; j++) {
        Real yc = pmb->pcoord->x2v(j) - ypos[nw];
        Real yc2 = yc*yc;
        for (int i=istl; i<=istu; i++) {
          Real xc = pmb->pcoord->x1v(i) - xpos[nw];
	  Real xc2 = xc*xc;
  	  Real r2 = xc2 + yc2 + zc2;
	  Real r = std::sqrt(r2);
          if (r < remap_radius) {
	    Real xy = xy = std::sqrt(xc2 + yc2);
	    Real sinphi = xy/r;
	    Real cosphi = zc/r;
            Real costhta = xc/xy;
	    Real sinthta = yc/xy;
	    Real rho = mdot[nw]/(4.0*pi*r2*vinf[nw]);
	    Real u1, u2, u3;
	    Real pre = (rho/avgmass)*boltzman*Twind;
            if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
	      u1 = vinf[nw]*sinphi*costhta;
	      u2 = vinf[nw]*sinphi*sinthta;
	      u3 = vinf[nw]*cosphi;
            } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) { // RPZ
	      u1 = vinf1*sinphi;
	      u2 = 0.0;
	      u3 = vinf1*cosphi;
	    }
	    cons(IDN,k,j,i) = rho;
            cons(IM1,k,j,i) = rho*u1;
            cons(IM2,k,j,i) = rho*u2;
            cons(IM3,k,j,i) = rho*u3;
            cons(IEN,k,j,i) = pre/gmma1 + 0.5*rho*vinf[nw]*vinf[nw];
	    if (NSCALARS > 0) {
              pmb->pscalars->s(0,k,j,i) = scalar[nw]*rho;
              pmb->pscalars->r(0,k,j,i) = scalar[nw];
            }

	  }
	}
      }
    }
  }


  
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {

	/*
	// This is worse than just letting it be adiabatic...
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
	*/


	/*
	// This doesn't seem to make the wind temperature any better when using PPM...
	// When using PLM it seems to make it a little worse. Therefore don't use it.
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




//========================================================================================
//! \fn void MeshBlock::FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
//		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
//  \brief Remap winds
//========================================================================================
void FixTemperature(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
		const AthenaArray<Real> &bcc, AthenaArray<Real> &cons){

  Real g = pmb->peos->GetGamma();
  gmma1 = g - 1.0;

  Real Twind = 1.0e4;         // K
  Real avgmass = 1.0e-24;     // g
  Real boltzman = 1.380658e-16;

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
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {	
    stagx = xpos2 - rob;
    stagy = 0.0;
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
void orbit_calc(Real t){

  float xdist,ydist,zdist;
  double time_offset,torbit,phase,phi,E,dE,cosE,sinE;
  double sii,coi,theta,ang;
  double rrel;      // radius vector
  double gamma_ang; // angle between velocity and radius vectors
  double a1,a2;     // semi-major axis of barycentric orbits
  double M;         // effective orbit mass
  double m1,m2,v1,v2;

  Real pi = 2.0*asin(1.0);

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

  Real Msol = 1.9891e33;
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

  // Orbit is in the xy plane, with the semi-major axis along the y-axis

  double sintheta = std::sin(theta);
  double costheta = std::cos(theta);
  double sinang = std::sin(ang);
  double cosang = std::cos(ang);
  
  zpos1 = 0.0;
  zpos2 = 0.0;
  zvel1 = 0.0;
  zvel2 = 0.0;
  xpos1 =  a1*rrel*sintheta;
  ypos1 = -a1*rrel*costheta;
  xpos2 = -a2*rrel*sintheta;
  ypos2 =  a2*rrel*costheta;
  xvel1 =  v1*sinang;
  yvel1 = -v1*cosang;
  xvel2 = -v2*sinang;
  yvel2 =  v2*cosang;

  // Update dsep
  xdist = xpos1 - xpos2;
  ydist = ypos1 - ypos2;
  zdist = zpos1 - zpos2;
  dsep = std::sqrt(xdist*xdist + ydist*ydist + zdist*zdist);

  //std::cout << "xpos1 = " << xpos1 << "; xpos2 = " << xpos2 << "; ypos1 = " << ypos1 << "; ypos2 = " << ypos2 << "; zpos1 = " << zpos1 << "; zpos2 = " << zpos2 << "\n";
  //std::cout << "dsep = " << dsep << "\n";
  //exit(EXIT_SUCCESS);
  
  return;
}