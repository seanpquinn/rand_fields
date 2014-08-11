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
#ifndef HAMMURABI_INTEGRATOR_H
#define HAMMURABI_INTEGRATOR_H

#include "proto_class_B_field2.h"
#include "proto_class_TE_density.h"
#include "proto_class_CRE.h"
#include "healpix_map.h"

struct struct_observables
{
  double I, Q, U;
  double RM;
  double DM;
  double tau;
  double ff;
  
  struct_observables()
  : I(0), Q(0), U(0), RM(0), DM(0), tau(0), ff(0) {}
  
  void zero()
  { I=Q=U=RM=DM=tau=ff=0; }
  
  struct_observables &operator+= (const struct_observables &other)
  {
    I+=other.I; Q+=other.Q; U+=other.U;
    RM+=other.RM; DM+=other.DM;
    tau+=other.tau;
    ff+=other.ff;
    return *this;
  }
};

struct struct_defl_obs
{
  double Bx;
  double By;
};

// These needed because of the order in which the includes go through;
// in fact, though they are included above, B_field and List use
// Integrator, so this stuff will be looked at before their stuff.  So
// these dummy declarations are needed.
class B_field;
class List;
// proto type file of class Integrator.
// should store all info about the spheres

class Integrator
{
 private:

  // All set in the constructor
  //---------------------------------
  int obs_NSIDE;
  int obs_shell, total_shell;
  int vec_size_R;

  double MAX_RADIUS;
  double delta_R;

  // Sun position in galactocentric cartesian
  // coordinates
  vec3 SunPosition;

  // Simulation box size
  double x_diameter, y_diameter, z_diameter;
  double x_Sun, y_Sun, z_Sun;

  int simulation_boundary_type;

  std::string obs_file_name;
  std::string obs_RM_file_name;
  std::string obs_DM_file_name;
  std::string obs_tau_file_name;
  std::string obs_ff_file_name;
  std::string PF_file_name;

  std::string defl_int_file_name;
  std::string defl_ori_file_name;

  bool pol_map_in_kelvin_flag;

  bool do_sync_emission;
  bool do_rm;
  bool do_dm;
  bool I_PI_PA_flag;
  bool PF_flag;

  // Sun et al. 2008 developments.
  bool do_tau;
  bool do_ff;

  bool do_TE_Bran_coupling;
  double TE_Bran_coupling;

  //---------------------------------
  // really_loud true returns more text
  // output
  bool really_loud;


  //---------------------------------
  // Healpix maps:
  //---------------------------------
  //--galactic emission maps---------
  Healpix_Map<double> obs_I_map;
  Healpix_Map<double> obs_Q_map;
  Healpix_Map<double> obs_U_map;
  Healpix_Map<double> obs_RM_map;
  Healpix_Map<double> obs_DM_map;

  Healpix_Map<double> obs_tau_map;
  Healpix_Map<double> obs_ff_map;

  //--deflection maps----------------
  Healpix_Map<double> obs_defl_Bx_map;
  Healpix_Map<double> obs_defl_By_map;
  //---------------------------------


  int get_subbeam_numb(int target_shell_numb, int current_shell_numb) const;
  int get_shell_NPIX(int shell_numb) const;
  int get_shell_NSIDE(int shell_numb) const;
  double get_delta_R(void) const;
  double get_max_shell_radius(int shell_numb) const;
  double get_min_shell_radius(int shell_numb) const;

  int up_ipix(int up_shell_numb, int current_shell_numb, int current_shell_ipix) const;
  int down_ipix(int below_shell_numb, int current_shell_numb, int current_shell_ipix) const;

  int return_shell_ipix_interval_max(int target_shell_numb, int current_shell_numb, int current_shell_ipix) const;
  int return_shell_ipix_interval_min(int target_shell_numb, int current_shell_numb, int current_shell_ipix) const;

  // radial integration for Galactic emission
  void radial_integration(const double min_shell_radius, const double max_shell_radius, const int current_shell, const int ipix, struct_observables &observables, B_field &mag, TE_density &ne_density, CRE &cre, const double obs_freq);

  // radial integration for deflection map
  void radial_integration(const double min_shell_radius, const double max_shell_radius, const int current_shell, const int ipix, struct_defl_obs &observables, B_field &mag);

  void read_shell_params(paramfile &params);

  bool check_simulation_boundaries (double R, double THE, double PHI) const;

  void get_synchrotron_emissivity (const double cre_C, const double distribution_index_p, const double Bperp, const double omega, double& P_perpendicular_value, double& P_parallel_value) const;

 public:

  Integrator(paramfile & params); // the constuctor
  Integrator() {}; // the empty constuctor, use next function after
  //  void set_shell_params(int obs_shell_numb_set, int total_shell_numb_set, int obs_NSIDE_set, int max_pix_set, int min_pix_set, int vec_size_R_set, double max_radius_set);
  void set_shell_params(int obs_shell_numb_set, int total_shell_numb_set, int obs_NSIDE_set, int vec_size_R_set, double max_radius_set, vec3 sunpos, double x_dia, double y_dia, double z_dia, int simulation_boundary_type_set, std::string obs_file_name_set, std::string obs_RM_file_name_set, std::string obs_DM_file_name_set, std::string obs_tau_file_name_set, std::string PF_file_name_set, std::string defl_int_file_name_set, std::string defl_ori_file_name_set, bool pol_map_in_kelvin_flag_set, bool I_PI_PA_flag, bool PF_flag_set, bool do_sync_emission_set, bool do_rm_set, bool do_dm_set, bool do_tau_set, bool do_ff_set, bool do_TE_Bran_coupling_set, double TE_Bran_coupling_set);

  ~Integrator(); // the destructor

  int get_total_obs_shell_difference() const;


  void integrate2obs_maps(B_field &mag, TE_density &ne_density, CRE &cre, const double obs_freq);
  void integrate2obs_list(List &list, B_field &mag, TE_density &ne_density, CRE &cre);

  void save_maps2file(void);

  // based on radial integration and integrate2obs_maps routines
  void get_deflection_map(B_field &mag);

  void save_defl_map2file(bool do_rigidity, double rigidity);

  friend class List;
  friend class B_field;
};
#endif
