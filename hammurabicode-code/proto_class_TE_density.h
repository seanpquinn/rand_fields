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

#ifndef HAMMURABI_TE_DENSITY_H
#define HAMMURABI_TE_DENSITY_H
#include "arr.h"
#include "vec3.h"
#include "paramfile.h"

#include "arr.h"
#include "vec3.h"

// Protype for class TE_density
class TE_density
{
 private:

  arr3<float> TE_grid;
  double TE_Lx,TE_Ly,TE_Lz;
  int TE_nx,TE_ny,TE_nz;
  bool tegrid;
  double TE_constant;
  vec3 SunPosition;

  void read_grid(const std::string &filename);

public:
  
  TE_density() : tegrid(false){};
  TE_density(paramfile &params);

  double get_density(double R, double THE, double PHI);
  double get_temperature(double R, double THE, double PHI);
  double get_filling_factor(double R, double THE, double PHI);
};

#endif
