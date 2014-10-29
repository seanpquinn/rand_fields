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

// Main for the Hammurabi code

// the header
#include "hammurabi.h"
#include "proto_class_List.h"
#include "CGS_units_file.h"
#include <write_bfield.cpp>

using namespace std;

int main (int argc, const char **argv)
{
  cout << " Starting the Hammurabi code " << endl;
//  module_startup ("hammurabi", argc, argv, 2, "<parameter file>");
  // From the Healpix package. Necessary for parameter
  // file reading. It was suggested by Martin to be left here.

  paramfile params(argv[1]);

  B_field magobj(params);
  TE_density ne_density(params);
  CRE cre(params);

  Integrator intobj(params);

  // Compares random field grid size and integrator cells
  magobj.SanityCheck(intobj); 

  // Does nothing if parameter file has B_field_do_random=F
  magobj.fillRandom();

  double obs_freq=(params.find<double>("obs_freq_GHz",9999.))*CGS_U_GHz;

  if (params.find<bool>("do_defl_map",false))
  {
    intobj.get_deflection_map(magobj);
    double rigidity=(params.find<double>("rigidity", 4.e19));
    double forefactor=((CGS_U_Joule/CGS_U_erg)/(CGS_U_Coulomb/CGS_U_esu));
    cout << " Rigidity is: " << rigidity << " V " << endl;
    cout << "          or: " << rigidity*forefactor << "erg/esu" << endl;
    rigidity*=forefactor*(CGS_U_erg/CGS_U_esu);
    bool do_rigidity=true;
    intobj.save_defl_map2file(do_rigidity, rigidity);
  }
  else if ( (params.find<std::string>("list_in_file_name","")).compare("") == 0 &&  !params.find<bool>("list_use_lonlat_range",false) ) {
    intobj.integrate2obs_maps(magobj, ne_density, cre, obs_freq);
    intobj.save_maps2file();
  } else {
    List listobj(params, intobj.get_total_obs_shell_difference());
    listobj.create_list(params, intobj.get_total_obs_shell_difference());

    // Attention, the following command overwrites list defined
    // frequencies
    //    listobj.load_full_list_freq(obs_freq);
    intobj.integrate2obs_list(listobj, magobj, ne_density, cre);
    listobj.save_list2file();
 
// J.West:	listobj.save_latlon();

  }

  bfield_out(params.find<int>("file_name_num",0), params.find<int>("B_field_type",7)); // [Modified by SPQ 5/8/14]
  Log(" Finishing the Hammurabi code\n");

  return 0;
}
