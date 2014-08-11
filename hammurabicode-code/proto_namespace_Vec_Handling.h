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
#ifndef HAMMURABI_NAMESPC_VEC_HANDLING_H
#define HAMMURABI_NAMESPC_VEC_HANDLING_H

#include "proto_class_B_field2.h"
#include "proto_class_TE_density.h"

// prototype file of namespace Vec_Handling.
// should deal with all vector operations.

namespace Vec_Handling
  {
  double return_perp2LOS (const vec3 &cart_vec, double THE_LOS, double PHI_LOS);
  double return_par2LOS (const vec3 &cart_vec, double THE_LOS, double PHI_LOS);
  double vec_inclination_on_sph_plane (double THE_ec, double PHI_ec, const vec3 &input_vec);
  double return_intr_pol_ang(const vec3 &B_vec, double THE, double PHI);

  vec3 return_LOS_unit_vec (double THE_LOS, double PHI_LOS);
  vec3 RTHEPHI2cart_vec(double R, double THE, double PHI);
  vec3 ec_cart_vec2gc_cart_vec(const vec3 &ec_cart_vec, const vec3 &SunPosition);
  void cart_coord2cyl_coord(const vec3 &cart_vec, double & r, double & phi, double & z);
  }

#endif
