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

#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

class CRE
{
 private:
  // Sun position in galactocentric cartesian coordinates
  vec3 SunPosition;

  // C1
  double C1_spec_index_p, C1_hr, C1_hz, C1_JE;

  // C2
  double C2_spec_index_p1, C2_spec_index_p2, C2_hr, C2_hz;

  // C3
  double C3_spec_index_p1, C3_spec_index_p2, C3_hr, C3_hz;

  // number of the C field
  int Cfield_type;

  double cre_model_WMAP_3th_year (double r, double z) const;
  double cre_model_Sun2008 (double r, double z, double hr, double hz) const;

 public:
  CRE(paramfile &params);

  double get_C(double R, double THE, double PHI) const;
  double get_spec_index_p(double obs_freq, double ec_R, double ec_THE, double ec_PHI) const;
};

#endif
