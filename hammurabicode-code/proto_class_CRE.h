/*
The "Hammurabi" code simulates Galactic synchrotron emission at all frequencies.
Copyright (C) 2006,2007,2008  Andre Waelkens, Tess Jaffe, Martin Reinecke, Torsten Ensslin

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

Andre Waelkens can be reached at waelkens@mpa-garching.mpg.de
*/
#ifndef HAMMURABI_CRE_H
#define HAMMURABI_CRE_H

// prototype for class CRE density
class CRE
{
 private:

  void read_CRE_params(paramfile &params);

  // Sun position in galactocentric cartesian
  // coordinates
  vec3 SunPosition;

  // C1
  double C1_spec_index_p;
  double C1_hr;
  double C1_hz;
  double alpha;  // used for C4 as well

  // C2
  double C2_spec_index_p1, C2_spec_index_p2;
  double C2_hr;
  double C2_hz;

  // C3
  double C3_spec_index_p1, C3_spec_index_p2;
  double C3_hr;
  double C3_hz;

  // C4
  double C4_spec_index_p;
  //double alpha; already defined
  double C4_ncre[81][21];
  
  // C5
  double C5_spec_index_p;
  double C5_hr;
  double C5_hz;
  double C5_ha;

  // C6
  double C6_spec_index_p;
  double C6_hr;
  double C6_hz;
  double C6_ha;

  // C7
  double C7_spec_index_p;
  double C7_ncre[81][21];
  double C7_ha;
  double C7_hb;

  // C8
  double C8_spec_index_p;
  double C8_ncre[81][21];
  double C8_zcre;


 public:

  // number of the C field
  int Cfield_type;

  CRE(paramfile &params); // the constructor
  CRE(); // the constructor
  ~CRE(); // the destructor

  double get_C(double R, double THE, double PHI) const;
  double get_spec_index_p(double obs_freq, double ec_R, double ec_THE, double ec_PHI) const; 

  //  double ncre_model_TD_thesis (double r, double z);
  double cre_model_WMAP_3th_year (double r, double z) const;
  double cre_model_Sun2008 (double r, double z) const;
  double cre_model_Sun2008_II (double r, double z) const;
  double cre_model_galprop (double r, double z) const;
  double cre_model_wmap_mod (double r, double z) const;
  double cre_model_JanssonFarrar2011 (double r, double z) const;
  double cre_model_galprop_mod (double r, double z) const;
  double cre_model_galprop_scaled (double r, double z) const;


  void setup_Cfield1(double p, double hr, double hz, double CRE_alpha);
  void setup_Cfield2(double p1, double p2, double hr, double hz);
  void setup_Cfield3(double p1, double p2, double hr, double hz);
  void setup_Cfield4(double p,  double CRE_alpha);
  void setup_Cfield5(double p, double ha,  double hr, double hz, double CRE_alpha);
  void setup_Cfield6(double p, double ha,  double hr, double hz, double CRE_alpha);
  void setup_Cfield7(double p, double ha,  double hb, double CRE_alpha);
  void setup_Cfield8(double p, double z_cre,  double CRE_alpha);
};
#endif
