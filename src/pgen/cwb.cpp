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


// CWB global namespace
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

  // Refinement parameters
  int G0_level = 0;

  // Simulation parameters
  int  remap         = 3;      // Remap radius in cell widths
  bool orbit         = false;  // Turn on orbital dynamics
  bool cool          = false;  // Turn on cooling
  bool dust          = false;  // Turn on dust evolution
  bool dust_cool     = false;  // Turn on dust cooling

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

void ThermalRemap(MeshBlock *pmb, const Real dt,
                  const AthenaArray<Real> &prim,
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

    for (int k = kstl; k <= kstu; ++k) {
      Real zc  = pmb->pcoord->x3v(k) - Star->x[2];
      Real zc2 = SQR(zc);
      for (int j = jstl; j <= jstu; ++j) {
        Real yc  = pmb->pcoord->x2v(j) - Star->x[1];
        Real yc2 = SQR(yc);
        for (int i = istl; i <=istu; ++i) {
          Real xc  = pmb->pcoord->x1v(i) - Star->x[0];
          Real xc2 = SQR(xc);
          // Calculate overall distance from star to 
          Real r2  = xc2 + yc2 + zc2;
          Real r   = std::sqrt(r2);

          if (r < Star->remap_cgs) {
            Real MDot = Star->mdot_cgs;
            Real VInf = Star->vinf;
            Real vol  = Star->remap_vol;

            Real rhoDot  = MDot / vol;
            Real mechLum = 0.5 * rhoDot * SQR(VInf);

            cons(IDN,k,j,i) += rhoDot * dt;
            cons(IEN,k,j,i) += mechLum * dt;
          }
        }
      }
    }
  }
  return;
}

void RemapWinds(AthenaArray<Real> &pcoord,
                AthenaArray<Real> &phydro,
                AthenaArray<Real> &pscalars,
                Real is, Real ie,
                Real js, Real je,
                Real ks, Real ke) {

  

  return;
}



// /*!  \brief Restrict the cooling rate at unresolved interfaces between hot 
//  *          diffuse gas and cold dense gas.
//  *
//  *   Replace deltaE with minimum of neighboring deltaEs at the interface.
//  *   Updates dei, which is positive if the gas is cooling.
//  *
//  *   \author Julian Pittard (Original version 13.09.11)
//  *   \version 1.0-stable (Evenstar):
//  *   \date Last modified: 13.09.11 (JMP)
//  */
void RestrictCool(int is,int ie,int js,int je,int ks,int ke,int nd,AthenaArray<Real> &dei,const AthenaArray<Real> &cons){

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
	pre(k,j,i) = g1*ie;
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

void CoolingFunction(MeshBlock *pmb, const Real dt,
                     AthenaArray<Real> &cons)
{
  AthenaArray<Real> dei(pmb->ke+1,pmb->je+1,pmb->ie+1);

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Read conserved variables
        Real rho = cons(IDN,k,j,i);        // Denisty (g cm^-3)
        Real u1  = cons(IM1,k,j,i) / rho;  // U_1 (cm s^-1)
        Real u2  = cons(IM2,k,j,i) / rho;  // U_2 (cm s^-1)
        Real u3  = cons(IM3,k,j,i) / rho;  // U_2 (cm s^-1)
        // Read scalars
        Real wm = pmb->pscalars->r(0,k,j,i);  // Wind mass fraction scalar
        if (dust_cool) {
          Real z = pmb->pscalars->r(ZLOC,k,j,i);  // Dust mass fraction
          Real a = pmb->pscalars->r(ALOC,k,j,i);  // Dust avg. radius (cm)
        }
        // Calculate fluid velocity
        Real v2  = SQR(u1) + SQR(u2) + SQR(u3);  // Square of scalar velocity
        Real v   = std::sqrt(v2);                // Scalar velocity (cm/s)
        // Calculate average mass of particle in cell
        Real avgmass = WR.avgm;
        // Calculate average unshocked wind temperature (K)
        Real TWind   = (wm * WR.twnd) + ((wm - 1.0) * OB.twnd); 
        // Calculate internal energy of cell
        Real ke  = 0.5 * rho * v2;        // Bulk kinetic energy of fluid (erg)
        Real ie  = cons(IEN,k,j,i) - ke;  // Internal energy of fluid (erg)
        Real pre = g1 * ie;               // Running pressure (Ba)
        Real preOld = pre;                // Original pressure (Ba)
        // Calculate Initial temperature
        Real T = (pre * avgmass) / (rho * kboltz);  // Running temperature (K)
        Real TOld = T;  // Original temperature (K)
        // Check for most egregious error, if temperature found to be invalid
        // This prevents the simulation from running on indefinitely
        // Even after it  has failed

        
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
        if (T < 1.5 * TWind) {
          T = TWind; 
        }
        else {
          Real dtInt = 0.0;  // Internal dt, for substepping
          Real lCool = 0.0;  // Cooling length (cm)
          while (dtInt < dt) {
            Real Eint    = pre/g1;
            Real lambda  = CCurve.FindLambda(T);
            Real edotGas = (SQR(rho)/SQR(massh)) * lambda;
            if (dust_cool) {
              // Do some dust cooling here 
            }
            Real totalCoolRate = edotGas;
            Real t_cool = Eint / totalCoolRate;
            // 
            Real tCool  = Eint / totalCoolRate;
            Real dtCool = 0.1 * std::abs(t_cool);
            // Rela 
            if ((dtInt + dtCool) > dt) {
              dtCool = dt - dtInt;
            }
            dtInt += dtCool;
            lCool += v*dtCool;
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

  RestrictCool(pmb->is,pmb->ie,
               pmb->js,pmb->je,
               pmb->ks,pmb->ke,
               pmb->pmy_mesh->ndim,
               dei,cons);
  // // Adjust pressure due to cooling
  AdjustPressureDueToCooling(pmb->is,pmb->ie,
                             pmb->js,pmb->je,
                             pmb->ks,pmb->ke,
                             dei,cons);
  return;
}

//! \fn void EvolveDust
//  \brief Function to evolve dust

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
        Real u1  = cons(IM1,k,j,i);  // Gas velocity x direction (cm s^-1)
        Real u2  = cons(IM2,k,j,i);  // Gas velocity y direction (cm s^-1)
        Real u3  = cons(IM3,k,j,i);  // Gas velocity z direction (cm s^-1)
        // Import scalars
        Real col = pmb->pscalars->s(CLOC,k,j,i) / rho;  // Wind colour
        Real a   = pmb->pscalars->s(ALOC,k,j,i) / rho;  // Grain radius (cm)
        Real z   = pmb->pscalars->s(ZLOC,k,j,i) / rho;  // Dust mass fraction 
        // Calculate temperature
        Real v2   = SQR(u1) + SQR(u2) + SQR(u3);
        Real ke   = 0.5 * rho * v2;
        Real pre  = (cons(IEN,k,j,i) - ke)*g1;
        Real temp = WR.avgm * pre / (rho * kboltz);
        // Replace a with minimum nucleation size and z if needed
        a = std::max(minGrainRadius,a);
        z = std::max(minDustToGasMassRatio,z);
        // Calculate gas number densities
        Real nH = rho * (10.0/14.0) / massh;
        Real nTot = 1.1 * nH;
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
        if (temp < 1.4e4) {
          // Dust growth occurs
          Real wa   = std::sqrt(3.0*kboltz*temp/(A*massh));
          Real dadt = 0.25 * eps_a * rho * wa / grainDensity;
          rhoD_dot  = 4.0 * PI * a2 * grainDensity * nD * dadt;
        }
        if (rhoD_dot != 0.0) {
          Real dRhoD   = rhoD_dot * dt;  // Integrate to find total change
          Real rhoDNew = rhoD + dRhoD;
               rhoDNew = std::max(0.0,rhoDNew);
          Real rhoNew  = rho + (rhoD - rhoDNew);
          // Calculate z
          Real zNew = rhoDNew / rhoNew;
               zNew = std::min(zNew,1.0);
               zNew = std::max(0.0,zNew);

          Real da   = dadt * dt;
          Real aNew = (a + da) / 1.0e-4;

          // Rewrite conserved arrays
          cons(IDN,k,j,i) = rhoNew;
          // Rewrite scalars
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
//  using a single long source function would be clunky, this funciton
//  is used to 

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
