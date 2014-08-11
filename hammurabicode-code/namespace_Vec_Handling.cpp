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

using namespace std;

// any function of this namespace should be addressed as:
// Vec_Handling::function();
// This is just a way of organizing functions. Nothing more.

namespace Vec_Handling {

/* Returns the perpendicular to LOS component of the vector, given the LOS vector */
double return_perp2LOS (const vec3 &cart_vec, double THE_LOS, double PHI_LOS)
  {
  vec3 unit_vec=return_LOS_unit_vec(THE_LOS, PHI_LOS);
  double sq_cart_vec_LOS=dotprod(unit_vec,cart_vec);
  sq_cart_vec_LOS*=sq_cart_vec_LOS;
  double sq_cart_vec=cart_vec.SquaredLength();

  // there might be numerical imprecision cases were
  // sq_cart_vec-sq_cart_vec_LOS<0, hence the abs.

  return sqrt(abs(sq_cart_vec-sq_cart_vec_LOS));
  }

/* Returns the parallel to LOS component of the vector, given the LOS vector */
double return_par2LOS (const vec3 &cart_vec, double THE_LOS, double PHI_LOS)
  {
  return dotprod(return_LOS_unit_vec(THE_LOS, PHI_LOS),cart_vec);
  }


/* Given LOS direction the input vector is projected on the sphere and then the pol ang is obtained according to R&L standard */
double return_intr_pol_ang(const vec3 &B_vec, double THE, double PHI)
  {
  double intrinsic_pol_ang=vec_inclination_on_sph_plane (THE, PHI, B_vec);
  // Take care that the B_vec is set correctly.
  // Correcting for the R&L standards
  if (intrinsic_pol_ang<0) intrinsic_pol_ang+=2.*CGS_U_pi;
  if (intrinsic_pol_ang<0 || intrinsic_pol_ang> 2.*CGS_U_pi) {cout << " In get_intrinsic_pol_ang, intrinsic pol ang outside definition range " << endl;}

  return intrinsic_pol_ang;
  }


// This routine should get the inclination angle of a
// vec on the plane of the sky.
// It projects the given vector on the ec sph unit vec
// THE and PHI, and uses atan2 to get the inclination.
// For atan2, the angle should be zero if the projection
// on sph_THE=0 and  sph_PHI>0 and pi/2 when sph_PHI=0, sph_THE>0.
// Correcting the angle to a less unconventional coord
// sys is just a matter of summing up a constant.
double vec_inclination_on_sph_plane (double THE_ec, double PHI_ec, const vec3 &input_vec)
  {
  vec3 sph_unit_vec_THE(cos(THE_ec)*cos(PHI_ec),cos(THE_ec)*sin(PHI_ec),-sin(THE_ec));
  vec3 sph_unit_vec_PHI(-sin(PHI_ec),cos(PHI_ec),0.0);

  double y_component=-dotprod(sph_unit_vec_THE, input_vec);// The y_component for atan2, is the projection on sph_unit_THE
  double x_component=-dotprod(sph_unit_vec_PHI, input_vec);// The x_component for atan2, is the projection on sph_unit_PHI

  return atan2(y_component, x_component);
  }

/* Returns the unit vec, given the LOS direction */
vec3 return_LOS_unit_vec (double THE_LOS, double PHI_LOS)
  { return vec3(pointing(THE_LOS,PHI_LOS));  }


vec3 RTHEPHI2cart_vec(double R, double THE, double PHI)
{
  planck_assert(R>=0,"R<0");
  return vec3(pointing(THE,PHI))*R;
}

vec3 ec_cart_vec2gc_cart_vec(const vec3 &ec_cart_vec, const vec3 &SunPosition)
  { return ec_cart_vec+SunPosition; }

void cart_coord2cyl_coord(const vec3 &cart_vec, double & r, double & phi, double & z)
  {
  r=sqrt(cart_vec.x*cart_vec.x+cart_vec.y*cart_vec.y);
  phi=atan2(cart_vec.y,cart_vec.x);
  if(phi<0.) {phi+=2.*CGS_U_pi;}
  z=cart_vec.z;
  }
}
