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

#ifndef HAMMURABI_LIST_H
#define HAMMURABI_LIST_H

#include "proto_class_Integrator.h"
#include "arr.h"
#include "paramfile.h"
#include "tess_tools.h"


class List
{
 private:
  struct struct_list_observables
  {
    double obs_freq;
    long ipix;
    // this should be the distance up to the polarized source whose
    // polarization angle obs defines the RM.
    // Useful if your source is not extragalactic.
    double source_dist;
    struct_observables observables;
  };

  double MAX_RADIUS;

  arr<struct_list_observables> data_obs_shell;
  arr<struct_list_observables> data_total_shell;
  int total_to_obs_factor;

  std::string list_in_file_name, list_out_file_name;

  bool pol_map_in_kelvin_flag;
  bool non_eg_pol_sources;
  bool item_obs_freq;
  double in_obs_freq;
	double lon_min;
	double lon_max;
	double lat_min;
	double lat_max;
	bool list_use_lonlat_range;
	int nside;

 public:

  List(paramfile &params, int shell_difference);
  ~List() {}

  List( int length, bool neps, bool iof, double obs_freq, double rmax, int shell_difference, std::string in_file_name, std::string out_file_name, bool kelvin_flag);


  void read_in_list();
  void create_list(paramfile &params, int shell_difference);

  void fill_list(arr<long> &listarr, arr<double> &obs_freq, arr<double> &dist);
  void fill_list(arr<long> &listarr, double obs_freq);

  long return_list_length() const
    { return data_obs_shell.size(); }
  long obs_pix_number (long list_position) const
    { return data_obs_shell[list_position].ipix; }

  double return_list_source_dist(int n) const
    { return data_obs_shell[n].source_dist;}

  // list frequency item
  //------------------------------------------------------------------
  bool list_freq_check(void) const
    {if(item_obs_freq){return true;}else{return false;}}
  double return_list_freq(int n) const
    {
      return data_obs_shell[n].obs_freq;
    }
  void load_full_list_freq(double ext_freq);
  //------------------------------------------------------------------

  void load_observables(long list_position,
                        const struct_observables &observed_values)
    { data_obs_shell[list_position].observables=observed_values; }

  struct_observables get_observables(long list_position)
    { return data_obs_shell[list_position].observables; }

  void add_observables(long list_position,
                       const struct_observables &observed_values)
    { data_obs_shell[list_position].observables+=observed_values; }

  void add_RM_observables(long list_position, int RM_pix, struct_observables observables, double max_shell_radius)
    {
      using namespace std;
      if(list_position>data_total_shell.size())
        {
          cerr << " Too large list_position in void add_RM_observables(long list_position, int RM_pix, double RM) " << endl;
          cerr << " list_position > data_total_shell.size(): " << list_position << " >  " << data_total_shell.size() << endl;
          exit(1);
        }

      if (data_total_shell[list_position].ipix!=RM_pix) {
        cerr<<" Previously determined RM total shell pixel number "<<data_total_shell[list_position].ipix<<" doesn't agree with that just calculated in Integrator: "<<RM_pix<<endl;
        exit(1);
      }

      data_total_shell[list_position].observables+=observables;
      data_total_shell[list_position].source_dist=max_shell_radius;
    }


  void zero() {
    for (long i=0;i<data_obs_shell.size();i++)
      data_obs_shell[i].observables.zero();

    for (long i=0;i<data_total_shell.size();i++)
      data_total_shell[i].observables.zero();
  }

  void save_list2file( std::string filename="" ) ;
  void save_latlon( std::string filename="" ) ;

  void return_list_physical( arr<long> &pixlist, arr<struct_observables> &obslist);

  friend class Integrator;

};

#endif
