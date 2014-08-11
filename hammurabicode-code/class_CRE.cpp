/*
The "Hammurabi" code simulates Galactic synchrotron emission at all frequencies.
Copyright (C) 2006,2007,2008,2012  Andre Waelkens, Tess Jaffe, Martin Reinecke, Torsten Ensslin

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/
#include "hammurabi.h"
#include "proto_namespace_Vec_Handling.h"
#include "proto_class_CRE.h"

using namespace std;

//-------------------------------------------------------------------------------------------------------
//Routines of class CRE
//-------------------------------------------------------------------------------------------------------

CRE::CRE(paramfile &params)
  {
  Cfield_type=params.find<int>("Cfield_type", 1);

  // Sun position in GC cartesian coordinates.
  SunPosition.x=CGS_U_kpc*params.find<double>("SunPosX", -8.5);
  SunPosition.y=CGS_U_kpc*params.find<double>("SunPosY", 0.);
  SunPosition.z=CGS_U_kpc*params.find<double>("SunPosZ", 0.);

  if(Cfield_type == 1)
    {
    C1_spec_index_p = params.find<double>("C1_p",3.);
    C1_hr = CGS_U_kpc*params.find<double>("C1_hr", 5.);
    C1_hz = CGS_U_kpc*params.find<double>("C1_hz", 1.);
    C1_JE = params.find<double>("C1_JE",0.25);
    }
  else if(Cfield_type ==2)
    {
    C2_spec_index_p1 = params.find<double>("C2_p1", 2.);
    C2_spec_index_p2 = params.find<double>("C2_p2", 3.);
    C2_hr = CGS_U_kpc*params.find<double>("C2_hr", 8.);
    C2_hz = CGS_U_kpc*params.find<double>("C2_hz", 1.);
    }
  else if(Cfield_type ==3)
    {
    C3_spec_index_p1 = params.find<double>("C3_p1", 2.75);
    C3_spec_index_p2 = params.find<double>("C3_p2", 3.);
    C3_hr = CGS_U_kpc*params.find<double>("C3_hr", 8.);
    C3_hz = CGS_U_kpc*params.find<double>("C3_hz", 1.);
    }
  else
    planck_fail("Invalid Cfield_type: "+dataToString(Cfield_type));
  }

// Obtains the C ("spatial term" of the CRE energy density
// distribution, see Ribicky&Lightman).
double CRE::get_C (double R, double THE, double PHI) const
{
  double r_gc, PHI_gc, z_gc;
  vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(R,THE,PHI);
  vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
  Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);

  if (Cfield_type==1)
    return cre_model_WMAP_3th_year (r_gc,z_gc);
  if (Cfield_type==2)
    return cre_model_Sun2008 (r_gc,z_gc,C2_hr,C2_hz);
  return cre_model_Sun2008 (r_gc,z_gc,C3_hr,C3_hz);
}

// this function might change a bunch in future
// it could e.g. depend on R, THE, PHI and obs_freq
double CRE::get_spec_index_p(double obs_freq, double ec_R, double ec_THE, double ec_PHI) const
  {
  if(Cfield_type==1) return C1_spec_index_p;
  if(Cfield_type==2) return (obs_freq>=0.408*CGS_U_GHz) ? C2_spec_index_p2 : C2_spec_index_p1;
  // Cfield_type==3, Sun2008_II, proof of concept case
  double r_gc, PHI_gc, z_gc;
  vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(ec_R,ec_THE,ec_PHI);
  vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
  Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);
  return (ec_R>8.*CGS_U_kpc) ? C3_spec_index_p1 : C3_spec_index_p2;
  }

// ncre model used by Page et al. WMAP 3th year polarization paper.
double CRE::cre_model_WMAP_3th_year (double r, double z) const
  {
  double h_r=C1_hr;
  double h_d=C1_hz;

  double gc_earth_radius=sqrt(dotprod(SunPosition,SunPosition));
  double C_zero=exp(-gc_earth_radius/h_r)*(2./(exp(0./h_d)+exp(-0./h_d)))*(2./(exp(0./h_d)+exp(-0./h_d)));

  // This comes from the Diploma thesis in chapter
  // relativistic electron density
  const double cre_calc_forefactor=4.*CGS_U_pi*(1./CGS_U_C_light)*C1_JE*CGS_U_MEC2/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec);

  const double cre_gamma_10=10.*CGS_U_GeV/CGS_U_MEC2; //We assume, or take, E=10GeV

  double C_earth, C_cre;
  C_earth=cre_calc_forefactor*pow(cre_gamma_10, C1_spec_index_p);       // Value obtained for earth in the thesis. The WMAP paper doesnt give any (not necessary for its calculations).
  C_cre=(C_earth/C_zero)*exp(-r/h_r)*(2./(exp(z/h_d)+exp(-z/h_d)))*(2./(exp(z/h_d)+exp(-z/h_d)));

  return C_cre;
  }

// ncre model used by Sun et al. 2008.
double CRE::cre_model_Sun2008 (double r, double z, double hr, double hz) const
  {
  // Note that Sun et al. point out other options for this 1kpc truncation.
  if(abs(z)>1.*CGS_U_kpc) return 0;

  double C_zero=(6.4e-5)/CGS_U_ccm;
  r=max(3.*CGS_U_kpc,r);
  return C_zero*exp(-((r-SunPosition.Length())/hr)-(abs(z)/hz));
  }
