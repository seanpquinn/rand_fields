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
#include <gsl/gsl_sf_gamma.h>
#include "CGS_units_file.h"
#include "hammurabi.h"
#include "fitshandle.h"
#include "healpix_map_fitsio.h"
#include "proto_class_List.h"
#include "proto_namespace_Vec_Handling.h"

using namespace std;

// this amazing class, will integrate things.

// the constructor
Integrator::Integrator(paramfile & params)
{
  cout << " class Integrator constructor " << endl;
  read_shell_params(params);
  delta_R=MAX_RADIUS/vec_size_R;

  if(do_sync_emission==true && do_rm==false)
    {cout << " Attention! Faraday depolarization is off !!" << endl;}
}

Integrator::~Integrator()
{
        //  cout << " class Integrator destructor " << endl;
}

void Integrator::read_shell_params(paramfile &params)
{
  // shell configuration
  obs_shell=params.find<int>("obs_shell_index_numb", 1);
  total_shell=params.find<int>("total_shell_numb", 1);

  obs_NSIDE=params.find<int>("obs_NSIDE", 16);

  vec_size_R=params.find<int>("vec_size_R", 16);

  MAX_RADIUS=CGS_U_kpc*params.find<double>("max_radius", 32);

  // Keep them in that order!
  //------------------------------------------------------------------
  // Sun position in GC cartesian coordinates.
  SunPosition.x=CGS_U_kpc*params.find<double>("SunPosX", -8.5);
  SunPosition.y=CGS_U_kpc*params.find<double>("SunPosY", 0.);
  SunPosition.z=CGS_U_kpc*params.find<double>("SunPosZ", 0.);

  // Simulation box size and Sun coordinates
  x_diameter=CGS_U_kpc*params.find<double>("x_diameter", 40.);
  y_diameter=CGS_U_kpc*params.find<double>("y_diameter", 40.);
  z_diameter=CGS_U_kpc*params.find<double>("z_diameter", 8.);

  // loading Sun coordinates
  x_Sun=(x_diameter/2.)-SunPosition.Length();
  y_Sun=y_diameter/2.;
  z_Sun=z_diameter/2.;
  //------------------------------------------------------------------

  // case of simulation boundaries. Boundaries should be
  // designed taking into account the sort of fields (so
  // that no information is cut of by the
  // check_simulation_boundaries routine).

  simulation_boundary_type=params.find<int>("simulation_boundary_type", 1);

  // output file names
  obs_file_name=params.find<string>("obs_file_name","IQU.fits");
  obs_RM_file_name=params.find<string>("obs_RM_file_name","RM.fits");
  obs_DM_file_name=params.find<string>("obs_DM_file_name","DM.fits");
  obs_tau_file_name=params.find<string>("obs_tau_file_name","tau.fits");
  obs_ff_file_name=params.find<string>("obs_ff_file_name","ff.fits");
  PF_file_name=params.find<string>("PF_file_name","pol_fraction.fits");

  defl_int_file_name=params.find<string>("defl_int_fname","defl_int.fits");
  defl_ori_file_name=params.find<string>("defl_ori_fname","defl_ori.fits");

  // output units for IQU, default is Kelvin
  pol_map_in_kelvin_flag=params.find<bool>("pol_map_in_kelvin_flag", true);

  // equivalently tho IQU one can use I_PI_PA maps.
  I_PI_PA_flag=params.find<bool>("I_PI_PA_flag", false);

  // it might be convenient to also get a pol. frac. map.
  PF_flag=params.find<bool>("PF_flag", false);

  // superfulous computations can be ignored
  do_sync_emission=params.find<bool>("do_sync_emission", true);
  do_rm=params.find<bool>("do_rm", true);
  do_dm=params.find<bool>("do_dm", false);
  do_tau=params.find<bool>("do_tau", false);
  do_ff=params.find<bool>("do_ff", false);

  if(do_ff==true && do_tau==false)
    {cout << " ATTETION! " << endl;
    cout << " do_ff==true && do_tau==false" << endl;
    cout << " cant compute ff without tau! setting do_tau==true " << endl;
    do_tau=true;}

  if(do_sync_emission==false && do_rm==false && do_dm==false && do_tau==false && params.find<bool>("do_defl_map", false) == false && do_ff==false)
    {cerr << "No IQU,RM,DM,tau,ff flag is true. Stopping Hammurabi." << endl; exit(1);}
  // really_loud==true results in more text output
  really_loud=params.find<bool>("class_Integrator loud", false);

  do_TE_Bran_coupling=params.find<bool>("do_TE_Bran_coupling", false);
  TE_Bran_coupling=params.find<double>("TE_Bran_coupling",0.01);
  if(do_rm==false && do_TE_Bran_coupling==true)
    {cout << " ATTENTION! " << endl;
    cout << " do_rm==false && do_TE_Bran_coupling==true " << endl;
    cout << " cant use TE_Bran_coupling wihtout do_rm==true " << endl;
    cout << " setting do_rm=true " << endl;
    do_rm=true;}

}

void Integrator::set_shell_params(int obs_shell_numb_set, int total_shell_numb_set, int obs_NSIDE_set, int vec_size_R_set, double max_radius_set, vec3 sunpos, double x_dia, double y_dia, double z_dia, int simulation_boundary_type_set, string obs_file_name_set, string obs_RM_file_name_set, string obs_DM_file_name_set, string obs_tau_file_name_set, string PF_file_name_set, string defl_int_file_name_set, string defl_ori_file_name_set, bool pol_map_in_kelvin_flag_set, bool I_PI_PA_flag_set, bool PF_flag_set, bool do_sync_emission_set, bool do_rm_set, bool do_dm_set, bool do_tau_set, bool do_ff_set, bool do_TE_Bran_coupling_set, double TE_Bran_coupling_set)
{
  obs_shell=obs_shell_numb_set;
  total_shell=total_shell_numb_set;
  obs_NSIDE=obs_NSIDE_set;
  vec_size_R=vec_size_R_set;
  MAX_RADIUS=max_radius_set*CGS_U_kpc;
  SunPosition=sunpos;
  x_diameter=x_dia;
  y_diameter=y_dia;
  z_diameter=z_dia;
  x_Sun=(x_diameter/2.)-SunPosition.Length();
  y_Sun=y_diameter/2.;
  z_Sun=z_diameter/2.;

  simulation_boundary_type=simulation_boundary_type_set;
  obs_file_name=obs_file_name_set;
  obs_RM_file_name=obs_RM_file_name_set;
  obs_DM_file_name=obs_DM_file_name_set;
  obs_tau_file_name=obs_tau_file_name_set;
  PF_file_name=PF_file_name_set;

  defl_int_file_name=defl_int_file_name_set;
  defl_ori_file_name=defl_ori_file_name_set;

  pol_map_in_kelvin_flag=pol_map_in_kelvin_flag_set;
  I_PI_PA_flag=I_PI_PA_flag_set;
  PF_flag=PF_flag_set;
  do_sync_emission=do_sync_emission_set;
  do_rm=do_rm_set;
  do_dm=do_dm_set;

  do_tau=do_tau_set;
  do_ff=do_ff_set;
  if(do_ff==true && do_tau==false)
   {cout << " ATTETION! " << endl;
    cout << " do_ff==true && do_tau==false" << endl;
    cout << " cant compute ff without tau! setting do_tau==true " << endl;
    do_tau=true;}

  really_loud=false;
  do_TE_Bran_coupling=do_TE_Bran_coupling_set;
  TE_Bran_coupling=TE_Bran_coupling_set;

  delta_R=MAX_RADIUS/vec_size_R;


}

// Obtains the IQU and RM maps.
// Note, RM map will always have NSIDE
// of total_shell, while IQU maps have obs_NSIDE.
void Integrator::integrate2obs_maps(B_field &mag, TE_density &ne_density, CRE &cre, const double obs_freq)
{
  cout << " integrate2obs_maps() initiated " << endl;

  cout << " Observation shell has NSIDE " << get_shell_NSIDE(obs_shell) << endl;

  //------------------------------------------------------------------
  // sets obs_IQU_maps
  if(do_sync_emission)
    {
      cout << " Allocating memory for obs shell " << endl;
      obs_I_map.SetNside(obs_NSIDE,NEST);
      obs_Q_map.SetNside(obs_NSIDE,NEST);
      obs_U_map.SetNside(obs_NSIDE,NEST);

      cout << " Zeroing obs shell before computation begin " << endl;
      for(int ipix=0;ipix<obs_I_map.Npix();ipix++)
        {
          obs_I_map[ipix]=0.;
          obs_Q_map[ipix]=0.;
          obs_U_map[ipix]=0.;
        }
      cout << " Zeroed!" << endl;
    }


  //------------------------------------------------------------------
  // sets obs_RM_map
  if(do_rm)
    {
      cout << " Observation RM shell has NSIDE " << get_shell_NSIDE(total_shell) << endl;

      cout << " Allocating memory for obs RM shell " << endl;
      // Note, RM map has NSIDE of total_shell. Because it does
      // not make sense to average over RM subbeams, like for
      // the Stokes parameters. To understand this consider how
      // quasi-monochromatic linearly polarized waves, i.e. a
      // bunch of linearly polarized waves with random phases,
      // add up. An average zero RM bunch of 4 subbeams will
      // not produce a linearly dependent on lambda square curve
      // for the measured polarization angle.

      obs_RM_map.SetNside(get_shell_NSIDE(total_shell),NEST);

      cout << " Zeroing obs RM shell before computation begin " << endl;
      for(int ipix=0;ipix<obs_RM_map.Npix();ipix++)
        {
          obs_RM_map[ipix]=0.;
        }
      cout << " Zeroed!" << endl;
    }
  //------------------------------------------------------------------
  // sets obs_DM_map
  if(do_dm)
    {
      cout << " Observation DM shell has NSIDE " << get_shell_NSIDE(total_shell) << endl;

      cout << " Allocating memory for obs DM shell " << endl;
      // like the RM, DM is not computed as an average.
      // DM are obtained by measuring the arrival time of
      // a pulsar pulse at different frequencies. It makes
      // no sense in averaging over several pulsars. DM
      // comes from a point source, there is no such thing as
      // diffuse RM.

      obs_DM_map.SetNside(get_shell_NSIDE(total_shell),NEST);

      cout << " Zeroing obs DM shell before computation begin " << endl;
      for(int ipix=0;ipix<obs_DM_map.Npix();ipix++)
        {
          obs_DM_map[ipix]=0.;
        }
      cout << " Zeroed!" << endl;
    }
  //------------------------------------------------------------------
  // sets obs_tau_map
  if(do_tau)
    {
      cout << " Observation tau shell has NSIDE " << get_shell_NSIDE(total_shell) << endl;

      cout << " Allocating memory for obs tau shell " << endl;
      obs_tau_map.SetNside(get_shell_NSIDE(total_shell),NEST);
      cout << " Zeroing obs tau shell before computation begin " << endl;
      for(int ipix=0;ipix<obs_tau_map.Npix();ipix++)
        {
          obs_tau_map[ipix]=0.;
        }
      cout << " Zeroed!" << endl;
    }
  //------------------------------------------------------------------
  // sets obs_ff_map
  if(do_ff)
    {
      cout << " Observation ff shell has NSIDE " << get_shell_NSIDE(total_shell) << endl;

      cout << " Allocating memory for obs ff shell " << endl;
      obs_ff_map.SetNside(obs_NSIDE,NEST);
      cout << " Zeroing obs ff shell before computation begin " << endl;
      for(int ipix=0;ipix<obs_ff_map.Npix();ipix++)
        {
          obs_ff_map[ipix]=0.;
        }
      cout << " Zeroed!" << endl;
    }

  //------------------------------------------------------------------

  // This for loop goes trough all shells, computes them for the proper
  // shell radiae etc...
  for (int current_shell=1;current_shell<(total_shell+1);current_shell++)
    {
      Log(" shell number: "); cout << current_shell << endl;
      cout << " NSIDE:        " << get_shell_NSIDE(current_shell) << endl;

      double min_shell_radius;
      double max_shell_radius;
      if(current_shell==1)
        {
          min_shell_radius=0.;
        }
      else
        {
          min_shell_radius=get_min_shell_radius(current_shell);
        }
      max_shell_radius=get_max_shell_radius(current_shell);

      cout << " Current shell interval is [" << min_shell_radius/CGS_U_kpc << " , " << max_shell_radius/CGS_U_kpc << "]" << endl;

      // Regardless of whether we allocate memory for the maps, we still declare them
      // The declaration occupies negligible amount of space.
      cout << " Declaring I,Q,U and RM, DM, tau and ff temporary Healpix_Maps " << endl;
      Healpix_Map<double> current_I_map;
      Healpix_Map<double> current_Q_map;
      Healpix_Map<double> current_U_map;
      Healpix_Map<double> current_RM_map;
      Healpix_Map<double> current_DM_map;
      Healpix_Map<double> current_tau_map;
      Healpix_Map<double> current_ff_map;
      cout << " temporary Healpix_Maps declared " << endl;


      //------------------------------------------------------------------
      // Now memory for the maps is allocated. We run trough all the bool
      // flags (do_sync_emission, ....).
      //------------------------------------------------------------------
      // synchrotron
      if(do_sync_emission)
        {
          cout << " Allocating I,Q,U temporary map space " << endl;

          current_I_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          current_Q_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          current_U_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          cout << " allocation done " << endl;
          cout << " Zeroing temporary maps " << endl;
          for(int ipix=0;ipix<current_I_map.Npix(); ipix++)
            {
              current_I_map[ipix]=0.;
              current_Q_map[ipix]=0.;
              current_U_map[ipix]=0.;
            }
          cout << " temporary maps zeroed " << endl;
        }
      //------------------------------------------------------------------
      // rm
      if(do_rm)
        {
          cout << " Allocating RM temporary map space " << endl;
          current_RM_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          cout << " allocation done " << endl;
          cout << " Zeroing temporary maps " << endl;
          for(int ipix=0;ipix<current_RM_map.Npix(); ipix++)
            {
              current_RM_map[ipix]=0.;
            }
          cout << " temporary maps zeroed " << endl;
        }
      //------------------------------------------------------------------
      // dm
      if(do_dm)
        {
          cout << " Allocating DM temporary map space " << endl;
          current_DM_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          cout << " allocation done " << endl;
          cout << " Zeroing temporary maps " << endl;
          for(int ipix=0;ipix<current_DM_map.Npix(); ipix++)
            {
              current_DM_map[ipix]=0.;
            }
          cout << " temporary maps zeroed " << endl;
        }
      //------------------------------------------------------------------
      // tau
      if(do_tau)
        {
          cout << " Allocating tau temporary map space " << endl;
          current_tau_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          cout << " allocation done " << endl;
          cout << " Zeroing temporary maps " << endl;
          for(int ipix=0;ipix<current_tau_map.Npix(); ipix++)
            {
              current_tau_map[ipix]=0.;
            }
          cout << " temporary maps zeroed " << endl;
        }
      //------------------------------------------------------------------
      // ff
      if(do_ff)
        {
          cout << " Allocating ff temporary map space " << endl;
          current_ff_map.SetNside(get_shell_NSIDE(current_shell),NEST);
          cout << " allocation done " << endl;
          cout << " Zeroing temporary maps " << endl;
          for(int ipix=0;ipix<current_ff_map.Npix(); ipix++)
            {
              current_ff_map[ipix]=0.;
            }
          cout << " temporary maps zeroed " << endl;
        }
      // Done! Memory for maps is allocated and they are zeroed.
      //------------------------------------------------------------------


      //------------------------------------------------------------------
      // Setting current_map_Npix.
      int current_map_Npix=0;
      // determining the current_map_Npix. If
      // neither do_sync_emission, nor do_rm, do dm
      // nor do_tau flags are true, nothing will
      // be computed in the next for loop!
      // Note that unlike for the obs maps, all
      // current maps have the same Npix. But since
      // one might not be defined, since no memory
      // was allocated, we run trough all.
      // Note that per construction there can be
      // no do_ff flag without a do_tau flag.
      if(do_sync_emission)
        {
          current_map_Npix=current_I_map.Npix();
        }
      else if(do_rm)
        {
          current_map_Npix=current_RM_map.Npix();
        }
      else if(do_dm)
        {
          current_map_Npix=current_DM_map.Npix();
        }
      else if(do_tau)
        {
          current_map_Npix=current_tau_map.Npix();
        }
      else
        {
          cerr << " Flags (I,Q,U,RM,DM,tau) are all false." << endl;
          cerr << " Nothing to be done. Stopping the code." << endl;
          exit(1);
        }

      // now, using the current_map_Npix, we cycle
      // over all to be integrated pixels.

      struct_observables observables;

      for(int ipix=0;ipix<current_map_Npix; ipix++)
        {

          pointing ptg;
          // if any flag is active ptg should be
          // loaded with something. If all flags are active
          // any of the if statements would have the same effect
          //-------------------------------------------------------
          if(do_sync_emission) {ptg=current_I_map.pix2ang(ipix);}
          else if(do_rm) {ptg=current_RM_map.pix2ang(ipix);}
          else if(do_dm) {ptg=current_DM_map.pix2ang(ipix);}
          else if(do_tau) {ptg=current_tau_map.pix2ang(ipix);}
          //-------------------------------------------------------

          // Get RM and/or tau and/or DM out from inner rings to pass in.
          // Note, this works since we integrate from the inner shell
          // in outward direction.
          if(do_rm) {observables.RM=obs_RM_map.interpolated_value(ptg);}
          if(do_dm) {observables.DM=obs_DM_map.interpolated_value(ptg);}
          if(do_tau) {observables.tau=obs_tau_map.interpolated_value(ptg);}

          radial_integration(min_shell_radius,max_shell_radius, current_shell, ipix, observables, mag, ne_density, cre, obs_freq);

          // Saves observables in respective current maps for each true flag.
          //------------------------------------------------------------------
          if(do_sync_emission)
            {
              current_I_map[ipix]=observables.I;
              current_Q_map[ipix]=observables.Q;
              current_U_map[ipix]=observables.U;
            }
          if(do_rm)
            {
              current_RM_map[ipix]=observables.RM;
            }
          if(do_dm)
            {
              current_DM_map[ipix]=observables.DM;
            }
          if(do_tau)
            {
              current_tau_map[ipix]=observables.tau;
            }
          if(do_ff)
            {
              current_ff_map[ipix]=observables.ff;
            }

          // observables saved in respective current maps
          //------------------------------------------------------------------
        }

      //------------------------------------------------------------------
      // Adds values to obs maps. This is done differently for IQU/ff and RM,tau.

      // Synchrotron and ff emission:
      if((do_sync_emission==true) || (do_ff==true))
        {
          // Add the values to obs_I,Q,U_map or obs_ff_map using whatever interpolation scheme.
          int obs_map_npix;
          if(do_sync_emission){obs_map_npix=obs_I_map.Npix();}
          else if(do_ff){obs_map_npix=obs_ff_map.Npix();}
          if(current_shell<=obs_shell)
            {
              for(int ipix=0;ipix<obs_map_npix;ipix++)
                {
                  pointing ptg;
                  if(do_sync_emission){
                    ptg=obs_I_map.pix2ang(ipix);

                    obs_I_map[ipix]+=current_I_map.interpolated_value(ptg);
                    obs_Q_map[ipix]+=current_Q_map.interpolated_value(ptg);
                    obs_U_map[ipix]+=current_U_map.interpolated_value(ptg);}

                  if(do_ff){
                    ptg=obs_ff_map.pix2ang(ipix);
                    obs_ff_map[ipix]+=current_ff_map.interpolated_value(ptg);}

                  // Alternatively
                  //      if(do_sync_emission){obs_I_map[ipix]+=current_I_map[down_ipix(current_shell,ipix)];}
                  //      if(do_ff){obs_ff_map[ipix]+=current_ff_map[down_ipix(current_shell,ipix)];}
                }
            }
          else if(current_shell>obs_shell)
            {
              for(int ipix=0;ipix<obs_map_npix;ipix++)
                {
                  /*
                    pointing ptg;
                    if(do_sync_emission){
                    ptg=obs_I_map.pix2ang(ipix);
                    obs_I_map[ipix]+=current_I_map.interpolated_value(ptg);
                    obs_Q_map[ipix]+=current_Q_map.interpolated_value(ptg);
                    obs_U_map[ipix]+=current_U_map.interpolated_value(ptg);}
                    if(do_ff){
                    ptg=obs_ff_map.pix2ang(ipix);
                    obs_ff_map[ipix]+=current_ff_map.interpolated_value(ptg);}
                  */
                  // Alternatively
                  int level_factor=1;
                  int ipix_manipulation=ipix;
                  for(int n=0;n<(current_shell-obs_shell);n++)
                    {
                      level_factor*=4;
                      ipix_manipulation*=4;
                    }
                  for(int n=ipix_manipulation;n<(ipix_manipulation+level_factor);n++)
                    {
                      if(do_sync_emission){
                        obs_I_map[ipix]+=current_I_map[n]/level_factor;
                        obs_Q_map[ipix]+=current_Q_map[n]/level_factor;
                        obs_U_map[ipix]+=current_U_map[n]/level_factor;}
                      if(do_ff){
                        obs_ff_map[ipix]+=current_ff_map[n]/level_factor;}
                    }
                }
            }
          else {cerr << " Erro when summing up shells! " << endl; exit(1);}
        }


      // RM and/or tau
      if(do_rm || do_dm || do_tau)
        {
          // RM or tau map is always stored to the highest possible res.
          for(int ipix=0;ipix<get_shell_NPIX(total_shell);ipix++)
            {
              // This interpolation scheme is not the same as in the
              // first hammurabi code, where every shell was interpolated
              // into the subsequent one till this eventualy became the
              // final sphere.
              /*
                pointing ptg;
                if(do_rm)
                {
                ptg=obs_RM_map.pix2ang(ipix);
                obs_RM_map[ipix]+=current_RM_map.interpolated_value(ptg);
                }
                if(do_dm)
                {
                ptg=obs_DM_map.pix2ang(ipix);
                obs_DM_map[ipix]+=current_DM_map.interpolated_value(ptg);
                }
                // and
                if(do_tau)
                {
                ptg=obs_tau_map.pix2ang(ipix);
                obs_tau_map[ipix]+=current_tau_map.interpolated_value(ptg);
                }
              */
              // note that if NSIDE of current_RM_map and obs_RM_map
              // are the same, interpolation corresponds to
              // obs_RM_map[ipix]+=current_RM_map[ipix].
              // Same for DM and tau.

              // Alternatively, with no interpolation
              int current_ipix=ipix;
              for(int n=0;n<(total_shell-current_shell);n++)
                {
                  current_ipix/=4;
                }
              if(do_rm){obs_RM_map[ipix]+=current_RM_map[current_ipix];}
              if(do_dm){obs_DM_map[ipix]+=current_DM_map[current_ipix];}
              if(do_tau){obs_tau_map[ipix]+=current_tau_map[current_ipix];}
            }
        }
      // Done adding values to respective obs shells
      //------------------------------------------------------------------

    }

  cout << " integrate2obs_maps() finished " << endl;
}



// Obtains the IQU, RM and tau lists. Does not use Healpix_Map
// obj, but instead Healpix_Base obj.
// Note, RM, tau lists will always have NSIDE
// of total_shell, while IQU map lists have obs_NSIDE.
// This is because it makes no sense to average
// over RM, tau subbeams. For RM it would yield results deviating from
// the linear lambda square dependence of the pol angle.
// This observation analogy does not work for tau, since I am
// not quite sure we can actualy observe that quantity. Still
// it is used in a very similar way to RM, hence we store it
// in the same way.
// Note, this routine, unlike the one for Healpix maps has no unique
// frequency at which the list pixels are evaluated.
void Integrator::integrate2obs_list(List &list, B_field &mag, TE_density &ne_density, CRE &cre)
{
  cout << " integrate2obs_list() initiated " << endl;

  cout << " Observation shell has NSIDE " << get_shell_NSIDE(obs_shell) << endl;
  cout << " Total shell has NSIDE " << get_shell_NSIDE(total_shell) << endl;

  // This for loop goes trough all shells, computes them for the proper
  // shell radiae etc...
  for (int current_shell=1;current_shell<(total_shell+1);current_shell++)
    {
      Log(" shell number: ");  cout << current_shell << endl;
      cout << " NSIDE:        " << get_shell_NSIDE(current_shell) << endl;

      double min_shell_radius;
      double max_shell_radius;
      if(current_shell==1)
        {
          min_shell_radius=0.;
        }
      else
        {
          min_shell_radius=get_min_shell_radius(current_shell);
        }

      // This loop runs trough the list of pixels
      // stored in the array.
      for(int n=0;n<list.return_list_length(); n++)
        {
          // checking for the obs_freq, present in the list
          // file.
          double obs_freq;
          obs_freq=list.return_list_freq(n);

          // Here we check if max_shell_radius is still smaller
          // than the specified maximum in the list input for
          // each entry n!
          max_shell_radius=get_max_shell_radius(current_shell);
          if(min_shell_radius>list.return_list_source_dist(n))
            {
              // for loop is done!
              continue;
            }
          else if(max_shell_radius>list.return_list_source_dist(n))
            {
              max_shell_radius=list.return_list_source_dist(n);
            }

          // checking if the pixel number is valid
          if(list.obs_pix_number(n)>get_shell_NPIX(obs_shell) || list.obs_pix_number(n)<0)
            {
              cerr << " Erro! ipix is outside allowed boundaries " << endl;
              cerr << " ipix = " << list.obs_pix_number(n) << endl;
              cerr << " does not belong to [0," << get_shell_NPIX(obs_shell) << "], as defined by obs_shell and NSIDE " << endl;
              cerr << " stoping code " << endl;
              exit(1);
            }
          // gets the corresponding pixel interval
          // for this schell. ipix_min==ipix_max for shells
          // <= obs_shell.
          int ipix_min=return_shell_ipix_interval_min(current_shell, obs_shell, list.obs_pix_number(n));
          int ipix_max=return_shell_ipix_interval_max(current_shell, obs_shell, list.obs_pix_number(n));

          struct_observables observables;
          double factor=get_subbeam_numb(current_shell, obs_shell);
          for(int ipix=ipix_min;ipix<ipix_max;ipix++)
            {
              // Get RM and tau out from inner rings to pass in.
              observables=list.get_observables(n);

              radial_integration(min_shell_radius,max_shell_radius,current_shell, ipix, observables, mag, ne_density, cre, obs_freq);
              observables.I/=factor;
              observables.Q/=factor;
              observables.U/=factor;
              observables.RM/=factor;
              observables.DM/=factor;
              // Note, no subbeam number weighting done
              // on the RM observable.
              list.add_observables(n, observables);

              // should be for current shell ipix, not for obs shell ipix.
              // previously: int RM_ipix_min, RM_ipix_max;
              int current_shell_ipix_min, current_shell_ipix_max;

              // Undo the renorm done above for adding IQU to total_shell.
              observables.I*=factor;
              observables.Q*=factor;
              observables.U*=factor;
              observables.RM*=factor;
              observables.DM*=factor;


              if(current_shell<obs_shell)
                {
                  current_shell_ipix_min=return_shell_ipix_interval_min(total_shell, obs_shell, list.obs_pix_number(n));
                  current_shell_ipix_max=return_shell_ipix_interval_max(total_shell, obs_shell, list.obs_pix_number(n));
                }
              else
                {
                  current_shell_ipix_min=return_shell_ipix_interval_min(total_shell, current_shell, ipix);
                  current_shell_ipix_max=return_shell_ipix_interval_max(total_shell, current_shell, ipix);
                }

              for(int m=current_shell_ipix_min;m<current_shell_ipix_max;m++)
                {
                  int list_position_total=n*get_subbeam_numb(total_shell, obs_shell);
                  int current_shell_pix_numb;
                  int extra_increment = (ipix-ipix_min)*(current_shell_ipix_max-current_shell_ipix_min);
                  int increment = m-current_shell_ipix_min+extra_increment;
                  if(really_loud) {cout << " increment " << increment << " extra increment " << extra_increment << endl;}
                  current_shell_pix_numb=list.obs_pix_number(n)*get_subbeam_numb(total_shell, obs_shell)+increment;
                  list.add_RM_observables(list_position_total+increment, current_shell_pix_numb, observables, max_shell_radius);
                }
            }
        }
    }
  cout << " integrate2obs_list() finished " << endl;
}


int Integrator::get_total_obs_shell_difference(void) const
{
  return (total_shell-obs_shell);
}


// returns the number of subbeams.
// 1 if shell_numb is <=obs_shell,
// factors of 4 larger otherwise.
int Integrator::get_subbeam_numb(int target_shell_numb, int current_shell_numb) const
{
  int subbeam_numb=1;
  if(target_shell_numb>current_shell_numb)
    {
      for(int n=current_shell_numb;n<target_shell_numb;n++)
        {
          subbeam_numb*=4;
        }
    }
  return subbeam_numb;
}


// returns the NSIDE of shell_numb
int Integrator::get_shell_NSIDE(int shell_numb) const
{
  if(shell_numb<1 || shell_numb>total_shell)
    {
      cerr << " Error in int Integrator::get_shell_NSIDE(int shell_numb) , shell_numb is outside definition " << endl;
      exit(1);
    }
  int difference =  shell_numb-obs_shell;
  int shell_NSIDE = obs_NSIDE;
  if(difference>0)
    {
      for(int n=0;n<difference;n++)
        {
          shell_NSIDE*=2;
        }
    }
  else if(difference<0)
    {
      for(int n=0;n<abs(difference);n++)
        {
          shell_NSIDE/=2;
        }
    }
  return shell_NSIDE;
}

int Integrator::get_shell_NPIX(int shell_numb) const
{
  int shell_NSIDE = 0;
  shell_NSIDE=get_shell_NSIDE(shell_numb);
  int shell_NPIX=0;
  shell_NPIX=12*shell_NSIDE*shell_NSIDE;
  return shell_NPIX;
}

// Computes the upper radial boundary of the shell
double Integrator::get_max_shell_radius(int shell_numb) const
{
  if(shell_numb<1 || shell_numb>(total_shell))
    {
      cerr << " invalid shell_numb in double Integrator::get_max_shell_radius(int shell_numb) " << endl;
      exit(1);
    }

  double max_shell_radius=MAX_RADIUS;
  for(int n=total_shell;n>shell_numb;n--)
    {
      max_shell_radius/=2.;
    }
  max_shell_radius/=delta_R;
  int step_number;
  step_number=int(max_shell_radius);
  max_shell_radius=step_number*delta_R;

  return max_shell_radius;
}

// Computes the lower radial boundary of the shell
double Integrator::get_min_shell_radius(int shell_numb) const
{
  if(shell_numb<1 || shell_numb>total_shell)
    {
      cerr << " invalid shell_numb in double Integrator::get_min_shell_radius(int shell_numb) " << endl;
      exit(1);
    }
  double min_shell_radius=MAX_RADIUS;
  for(int n=total_shell;n>(shell_numb-1);n--)
    {
      min_shell_radius/=2.;
    }
  // Need to reset the radius to something close to
  // a multiple of delta_R.
  min_shell_radius/=delta_R;
  int step_number;
  step_number=int(min_shell_radius);
  min_shell_radius=step_number*delta_R;

  // The min_shell_radius for the first shell is always zero.
  // First shell has shell_numb=1
  if (shell_numb==1)
    {
      min_shell_radius=0.;
    }
  return min_shell_radius;
}

double Integrator::get_delta_R(void) const
{
  return MAX_RADIUS/vec_size_R;
}

int Integrator::up_ipix(int up_shell_numb, int current_shell_numb, int current_shell_ipix) const
{
  if(up_shell_numb<current_shell_numb)
    {
      cerr << " Error in int Integrator::up_ipix(int up_shell_numb, int current_shell_numb, int current_shell_ipix), up_shell_numb<current_shell_numb " << endl;
      cerr << up_shell_numb << ", " << current_shell_numb << endl;
      exit(1);
    }
  int result=current_shell_ipix;
  for(int n=0;n<(up_shell_numb-current_shell_numb);n++)
    {
      result*=4;
    }
  return result;
}

int Integrator::down_ipix(int below_shell_numb, int current_shell_numb, int current_shell_ipix) const
{
  // might comment this if loop out. not something which should happen anyways.
  if(below_shell_numb>current_shell_numb)
    {
      cerr << " Error in int Integrator::down_ipix(int below_shell_numb, int current_shell_numb, int current_shell_ipix), below_shell_numb>current_shell_numb " << endl;
      cerr << below_shell_numb << ", " << current_shell_numb << endl;
      exit(1);
    }
  int result=current_shell_ipix;
  for(int n=0;n<(current_shell_numb-below_shell_numb);n++)
    {
      result/=4;
    }
  return result;
}

// returns the min ipix for the target shell
// corresponding to the current shell pixel given
int Integrator::return_shell_ipix_interval_min(int target_shell_numb, int current_shell_numb, int current_shell_ipix) const
{
  if(target_shell_numb<0 || current_shell_numb < 0 || target_shell_numb>total_shell || current_shell_numb> total_shell)
    {
      cerr << " Error! Wrong shell numbers in int Integrator::return_shell_ipix_interval_min(...) " << endl;
      exit(1);
    }

  int result;
  if(target_shell_numb>current_shell_numb)
    {
      result=up_ipix(target_shell_numb, current_shell_numb, current_shell_ipix);
    }
  else
    {
      result=down_ipix(target_shell_numb, current_shell_numb, current_shell_ipix);
    }
  return result;
}



// returns the max ipix + 1 (so it's suitable for a
// standard for loop, like for(int n=ipix_min;n<ipix_max;n++)
// for the current shell corresponding to the obs
// shell pixel given
int Integrator::return_shell_ipix_interval_max(int target_shell_numb, int current_shell_numb, int current_shell_ipix) const
{
  int result=0;
  result=return_shell_ipix_interval_min(target_shell_numb, current_shell_numb, current_shell_ipix)+get_subbeam_numb(target_shell_numb, current_shell_numb);
  return result;
}


bool Integrator::check_simulation_boundaries (double R, double THE, double PHI) const
{
  switch(simulation_boundary_type)
    {
    case 1:
      double x,z,y;
      // first sets earth centric coords
      x=R*cos(PHI)*sin(THE);
      y=R*sin(PHI)*sin(THE);
      z=R*cos(THE);
      // then sets array coords
      x=x+x_Sun;
      y=y+y_Sun;
      z=z+z_Sun;
      // then checks if outside sim box
      if (x<0||x>x_diameter) {/*cout << "Outside sim box. Stopping calculation. " << R/CGS_U_kpc << "; " << THE << "; " << PHI <<
                                endl;*/ return true;}
      if (y<0||y>y_diameter) {/*cout << "Outside sim box. Stopping calculation. " << R/CGS_U_kpc << "; " << THE << "; " << PHI <<
                                endl;*/ return true;}
      if (z<0||z>z_diameter) {/*cout << "Outside sim box. Stopping calculation. " << R/CGS_U_kpc << "; " << THE << "; " << PHI <<
                                endl;*/ /*R=loop_end*/ return true;}

      return false;
      break;
    case 2:
      return false;
      break;
    default:
      cout << " Warning no valid case for check_simulation_boundaries chosen !" << endl;
      cout << " Returning false (i.e. simulation boundaries are defined by max radius)" << endl;
      exit(1);
      return false;
      break;
    }
}

// Saves all 4 outputs (I,Q,U and RM).
// Note that obs_freq is necessary to set intensity_to_br_temp_factor,
// it MUST be the same freq at which the maps where first computed!!!
// That's dangerous still gonna fix that!
void Integrator::save_maps2file()
{
  cout << " void Integrator::save_maps2file(void) initiated " << endl;

  cout << " setting physical units to maps " << endl;

  double unit_corr_factor;
  if (pol_map_in_kelvin_flag)
    {
      unit_corr_factor=CGS_U_cm*CGS_U_cm/CGS_U_erg; // Otherwise we display values in code internal untis
    }
  else
    {
      unit_corr_factor=1.0;
    }

  cout << " maps have now physical units " << endl;

  if(do_sync_emission)
    {
      cout << " the IQU obs map " << endl;
      for(int ipix=0;ipix<obs_I_map.Npix();ipix++)
        {
          obs_I_map[ipix]*=(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
          obs_Q_map[ipix]*=(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
          obs_U_map[ipix]*=(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
        }

      if(PF_flag)
        {
          cout << " PF_flag is true, generatin pol. frac. map " << endl;
          Healpix_Map<double> temp_PF_map;
          temp_PF_map.SetNside(obs_I_map.Nside(), NEST);
          for(int ipix=0;ipix<obs_I_map.Npix();ipix++)
            {
              temp_PF_map[ipix]=sqrt(obs_Q_map[ipix]*obs_Q_map[ipix]+obs_U_map[ipix]*obs_U_map[ipix])/obs_I_map[ipix];
            }
          // writing the file
          fitshandle out_PF;
          out_PF.create(string("!")+PF_file_name);
          write_Healpix_map_to_fits(out_PF, temp_PF_map, PLANCK_FLOAT64);
          out_PF.close();
          cout << " pol. fraction map saved to " << PF_file_name << endl;
        }

      // instead of outoputting IQU, one can output
      // I, PI, pol_ang
      if(I_PI_PA_flag)
        {
          cout << " Attention, I_PI_PA_flag is true, " << endl;
          cout << " saving I_PI_PA instead of IQU " << endl;
          Healpix_Map<double> temp_PI_map, temp_PA_map;
          temp_PI_map.SetNside(obs_I_map.Nside(), NEST);
          temp_PA_map.SetNside(obs_I_map.Nside(), NEST);

          for(int ipix=0;ipix<obs_I_map.Npix();ipix++)
            {
              temp_PI_map[ipix]=sqrt(obs_Q_map[ipix]*obs_Q_map[ipix]+obs_U_map[ipix]*obs_U_map[ipix]);
              temp_PA_map[ipix]=(1./2.)*atan2(obs_U_map[ipix], obs_Q_map[ipix]);
              temp_PA_map[ipix]=CGS_U_pi-(temp_PA_map[ipix]+CGS_U_pi/2.);
              //  if(temp_PA_map[ipix]<0.) {temp_PA_map[ipix]+=CGS_U_pi;}


            }
          for(int ipix=0;ipix<obs_I_map.Npix();ipix++)
            {
              obs_Q_map[ipix]=temp_PI_map[ipix];
              obs_U_map[ipix]=temp_PA_map[ipix];
              // uncomment this for fractional polarization map
              // obs_U_map[ipix]=temp_PI_map[ipix]/obs_I_map[ipix];
            }
        }

      // writing the file
      fitshandle out;
      out.create(string("!")+obs_file_name);
      write_Healpix_map_to_fits(out, obs_I_map, obs_Q_map, obs_U_map, PLANCK_FLOAT64);
      cout << " obs map saved to " << obs_file_name << endl;
    }
  if(do_rm)
    {
      cout << " the RM map " << endl;
      double rm_unit_factor=CGS_U_m*CGS_U_m;
      for(int n=0;n<obs_RM_map.Npix();n++)
        {
          obs_RM_map[n]*=rm_unit_factor;
        }
      // writing the file
      fitshandle out_RM;
      out_RM.create(string("!")+obs_RM_file_name);
      write_Healpix_map_to_fits(out_RM,obs_RM_map,PLANCK_FLOAT64);
      cout << " obs RM map saved to " << obs_RM_file_name << endl;
    }
  if(do_dm)
    {
      cout << " the DM map " << endl;
      double dm_unit_factor=CGS_U_ccm/CGS_U_pc;
      for(int n=0;n<obs_DM_map.Npix();n++)
        {
          obs_DM_map[n]*=dm_unit_factor;
        }
      // writing the file
      fitshandle out_DM;
      out_DM.create(string("!")+obs_DM_file_name);
      write_Healpix_map_to_fits(out_DM,obs_DM_map,PLANCK_FLOAT64);
      out_DM.close();
      cout << " obs DM map saved to " << obs_DM_file_name << endl;
    }
  if(do_tau)
    {
      cout << " the tau map " << endl;
      // normalizing
      for(int n=0;n<obs_tau_map.Npix();n++)
        {
          obs_tau_map[n]*=1.;
        }
      // writing the file
      fitshandle out_tau;
      out_tau.create(string("!")+obs_tau_file_name);
      write_Healpix_map_to_fits(out_tau,obs_tau_map,PLANCK_FLOAT64);
      cout << " obs tau map saved to " << obs_tau_file_name << endl;
    }
  if(do_ff)
    {
      cout << " the ff map " << endl;
      // writing the file
      fitshandle out_ff;
      out_ff.create(string("!")+obs_ff_file_name);
      write_Healpix_map_to_fits(out_ff,obs_ff_map,PLANCK_FLOAT64);
      cout << " obs ff map saved to " << obs_ff_file_name << endl;
    }
  cout << " void Integrator::save_maps2file(void) finished " << endl;
}



// performs the integration for one pixel in a shell in [min_shell_radius, max_shell_radius].
// Requires class B_field, class TE_density, class CRE obj.
// Stores I, Q, U and RM and TAU values in the corresponding maps of class Integrator.
//
// Input observables contain only RM/TAU from inner shells out to here.
// This will be erased and replaced with the RM for *just* this shell,
// but that input value will be used to rotate the polarisation
// vectors.
//
void Integrator::radial_integration(const double min_shell_radius, const double max_shell_radius, const int current_shell, const int ipix,  struct_observables & observables, B_field &mag, TE_density &ne_density, CRE &cre, const double obs_freq)
{
  // tau and RM,DM are treated in analogous fashion
  double inner_shells_RM=0.;
  double inner_shells_DM=0.;
  double inner_shells_tau=0.;
  if(do_rm) {inner_shells_RM=observables.RM;}
  else{inner_shells_RM=0.;}
  if(do_dm) {inner_shells_DM=observables.DM;}
  else{inner_shells_DM=0.;}
  if(do_tau) {inner_shells_tau=observables.tau;}
  else{inner_shells_tau=0.;}
  observables.zero();

  pointing ptg_obj;
  Healpix_Base shell_base_obj(get_shell_NSIDE(current_shell), NEST, SET_NSIDE);
  ptg_obj=shell_base_obj.pix2ang(ipix); // Will have the angles theta and phi in the ptg object.

  double THE, PHI;

  THE=ptg_obj.theta;
  PHI=ptg_obj.phi;

  double l,b;

  l=PHI;
  b=(CGS_U_pi/2.)-THE;

  // In our LOS integration we consider the "domain" of the
  // input parameters (ne, ncre, B_field, ...) to range from
  // -delta_r/2. to delta_r/2. of the current distance from
  // the observer. (previous_)delta_rm, (previous_)delta_dm,
  // (previous_)delta_tau refer to that interval. However we
  // evaluate our functions (e.g. synchrotron emissivity) at
  // the distance. Thus the current RM, DM, tau, ... is:
  // RM(distance)=RM(distance - delta_r)+(previous_delta_RM+delta_RM)/2.

  // for computing rm
  // default values have to be zero!
  double delta_RM=0.;
  double previous_delta_RM=0.;

  // for computing dm
  // default values have to be zero!
  double delta_DM=0.;
  double previous_delta_DM=0.;

  // for computing tau
  // default values have to be zero!
  double delta_tau=0.;
  double previous_delta_tau=0.;

  // later delta_tau_average=(previous_delta_tau+delta_tau)/2.
  // This is here explicitly, unlike what happens for RM and DM
  // because we need this value for computing the free-free
  // integral, which is and integral over the optical depth tau.
  double delta_tau_average=0.;

  double P_pol_real=0., P_pol_imaginary=0., P_tot=0.;

  // At every run, this loop uses values of current_step and current_step-1.
  for (double dist=min_shell_radius; dist<(max_shell_radius+(0.1*delta_R)); dist=dist+delta_R)
    {
      // First checks if we are not already far enough.
      if(check_simulation_boundaries(dist, THE,PHI)) {break;}

      double B_per, B_par, intr_pol_ang;

      // Getting galactocentric position from
      // earh-centric dist,THE,PHI
      vec3 B_vec;
      vec3 pos;
      pos=Vec_Handling::RTHEPHI2cart_vec(dist, THE, PHI);

      pos.x=pos.x+SunPosition.x;
      pos.y=pos.y+SunPosition.y;
      pos.z=pos.z+SunPosition.z;

      // this call assumes pos is galactocentric
      B_vec=mag.return_B_cart(pos);
      B_par=Vec_Handling::return_par2LOS(B_vec, THE,PHI);
      B_per=Vec_Handling::return_perp2LOS(B_vec, THE,PHI);
      intr_pol_ang=Vec_Handling::return_intr_pol_ang(B_vec, THE, PHI);

      // for RM and/or DM and/or tau computations
      if(do_rm || do_dm || do_tau)
        {
          double ne;
          ne=ne_density.get_density(dist,THE,PHI);

          if(do_rm)
            {
              const double RM_forefactor=(CGS_U_qe*CGS_U_qe*CGS_U_qe)/(2.0*CGS_U_pi*CGS_U_MEC2*CGS_U_MEC2);

              // Sun et al. 2008, sec. 6.4.2
              // this might become more sophisticated one day.
              //-------------------------------------------------------
              if(do_TE_Bran_coupling){
                double B_par_reg, B_par_ran;
                B_par_reg=Vec_Handling::return_par2LOS(mag.return_Breg_cart(pos), THE,PHI);
                B_par_ran=Vec_Handling::return_par2LOS(mag.return_Bran_cart(pos), THE,PHI);
                delta_RM=ne*(B_par_reg+B_par_ran/sqrt(TE_Bran_coupling))*delta_R*RM_forefactor;
              }
              //-------------------------------------------------------
              else{
              delta_RM=ne*B_par*delta_R*RM_forefactor;
              }
            }
          if(do_dm)
            {
              delta_DM=ne*delta_R;
            }
          if(do_tau)
            {
              double delta_EM=ne*ne*delta_R;
              // filling factor correction from Sun et al. 2008
              delta_EM/=ne_density.get_filling_factor(dist,THE,PHI);

              double Ti=ne_density.get_temperature(dist, THE,PHI);
              // From Rolfs and Wilson 2000
              const double tau_forefactor=8.235e-2;
              // Note! obs_freq has to be in GHz
              // Ti in K and ne in cm^-3.
              // delta_R has to be in pc!
              double delta_EM_in_pc_per_ccm=delta_EM*CGS_U_ccm*CGS_U_ccm/CGS_U_pc;

              delta_tau=delta_EM_in_pc_per_ccm*pow(Ti,-1.35)*tau_forefactor*pow(obs_freq/CGS_U_GHz,-2.1);
              // delta_tau is adimensional
            }


          if (dist==min_shell_radius)
            {
              if(do_rm){previous_delta_RM=delta_RM;}
              if(do_dm){previous_delta_DM=delta_DM;}
              if(do_tau){previous_delta_tau=delta_tau;}

              if (dist==0)
                {
                      observables.RM=0;
                      observables.DM=0;
                      observables.tau=0;
                }
              continue;
            }

          // This is just the RM, DM or tau from the inner edge of *this*
          // shell to the current position.  Will be added to inner shells
          // outside this function.
          if(do_rm)
            {
              observables.RM+=(previous_delta_RM/2+delta_RM/2);
              previous_delta_RM=delta_RM;
            }
          if(do_dm)
            {
              observables.DM+=(previous_delta_DM/2+delta_DM/2);
              previous_delta_DM=delta_DM;
            }
          if(do_tau)
            {
              delta_tau_average=(previous_delta_tau/2+delta_tau/2);
              observables.tau+=delta_tau_average;
              previous_delta_tau=delta_tau;
            }
        }

      // for synchrotron emission computations
      if(do_sync_emission)
        {
          // in case do_rm && do_dm && do tau is false we need
          // to call a continue for dist==min_shell_radius
          if(do_rm==false && do_dm==false && do_tau==false)
            {if (dist==min_shell_radius){continue;}}

          double c_cre;
          c_cre=cre.get_C (dist,THE,PHI);

          double P_perpendicular_value=0.;
          double P_parallel_value=0.;
          double spec_index_p=cre.get_spec_index_p(obs_freq,dist,THE,PHI);

          double omega_obs=2.*CGS_U_pi*obs_freq;

          get_synchrotron_emissivity(c_cre, spec_index_p, B_per, omega_obs,P_perpendicular_value,P_parallel_value);

          double lambda_square;

          lambda_square=((2.*CGS_U_pi*CGS_U_C_light)/omega_obs)*((2.*CGS_U_pi*CGS_U_C_light)/omega_obs);

          double qui;

          // Add inner shells' total RM to the current shell's RM for
          // rotating the vectors.
          // Note that for no do_rm inner_shells_RM+observables.RM=0 .
          qui=(inner_shells_RM+observables.RM)*lambda_square-intr_pol_ang;


          if (P_perpendicular_value<0 || P_parallel_value<0) {cout << " Error, I think, P_perp or P_parallel are < 0 " << P_perpendicular_value << " " <<       P_parallel_value << endl;}

          double pol_emissivity_func_of_R, tot_emissivity_func_of_R;
          pol_emissivity_func_of_R=P_perpendicular_value - P_parallel_value;
          tot_emissivity_func_of_R=P_perpendicular_value + P_parallel_value;

          double real_part, imaginary_part;
          real_part=cos(2.0*qui);
          imaginary_part=sin(2.0*qui);


          P_pol_real=P_pol_real+pol_emissivity_func_of_R*real_part*delta_R;
          P_pol_imaginary=P_pol_imaginary+pol_emissivity_func_of_R*imaginary_part*delta_R;

          P_tot=P_tot+tot_emissivity_func_of_R*delta_R;

          double intensity_to_br_temp_factor;
          if (pol_map_in_kelvin_flag)
            {
              intensity_to_br_temp_factor=CGS_U_C_light*CGS_U_C_light/(2.0*CGS_U_kB*obs_freq*obs_freq);
            }
          else
            {
              intensity_to_br_temp_factor=1.0;
            }

          // for Sun et al. 2008 low frequency attenuation by
          // free-free absorption.
          //-------------------------------------------------------
          double attenuation_factor=1.;
          if(do_tau)
            {
              attenuation_factor=exp(-observables.tau);
            }
          //-------------------------------------------------------

          observables.I+=tot_emissivity_func_of_R*delta_R*intensity_to_br_temp_factor*attenuation_factor;
          observables.Q+=pol_emissivity_func_of_R*real_part*delta_R*intensity_to_br_temp_factor*attenuation_factor;
          observables.U+=pol_emissivity_func_of_R*imaginary_part*delta_R*intensity_to_br_temp_factor*attenuation_factor;
        }
      // Free-free emission
      if(do_ff)
        {
          // Unlike in the synchrotron case, we do not need
          // to call a continue for dist==min_shell_radius.
          // Note that do_tau will always be true if do_ff is.

          double previous_dist=dist-delta_R;
          if(previous_dist<0)
            {
              if(abs(previous_dist)<(delta_R+0.1*delta_R))
                {
                  cerr << " previous_dist<-delta_R!! Stopping code! " << endl; exit(1);
                }
              else
                {
                  previous_dist=0.;
                }
            }

          // Here the convention is that tau grows with distance.
          observables.ff+=0.5*(ne_density.get_temperature(previous_dist, THE, PHI)+ne_density.get_temperature(dist, THE,PHI))*(exp(-(observables.tau)))*(exp(+delta_tau_average)-1.);
        }
    }
}

// overloaded radial_integration .
// used for computing deflection maps.

void Integrator::radial_integration(const double min_shell_radius, const double max_shell_radius, const int current_shell, const int ipix,  struct_defl_obs & observables, B_field &mag)
{
  observables.Bx=observables.By=0;

  pointing ptg_obj;
  Healpix_Base shell_base_obj(get_shell_NSIDE(current_shell), NEST, SET_NSIDE);
  ptg_obj=shell_base_obj.pix2ang(ipix); // Will have the angles theta and phi in the ptg object.

  double THE, PHI;

  THE=ptg_obj.theta;
  PHI=ptg_obj.phi;

  double l,b;

  l=PHI;
  b=(CGS_U_pi/2.)-THE;

  // At every run, this loop uses values of current_step and current_step-1.
  for (double dist=min_shell_radius; dist<(max_shell_radius+(0.1*delta_R)); dist=dist+delta_R)
    {
      // The code is not ready at this stage to use random fields for
      // deflection computations. If there is no attenuation in the random
      // field, so that it is zero at the simulation box borders (currently not
      // the case), then there will be non physical simulation-box-border features.
      if(mag.return_dorand()) {
        cerr << " Error in void Integrator::radial_integration(...) for deflection maps " << endl;
        cerr << " do_rand==true! Unless your random field box is such that the field is attenuated " << endl;
        cerr << " to zero at its borders, this computation will have simulation-box border features! " << endl;
        cerr << " Stopping the code. " << endl;
        exit(1);
      }

      // Note! No simulation boundaries are checked here! Since there is no simulation box, see above.

      double B_per, intr_pol_ang;

      // Getting galactocentric position from
      // earh-centric dist,THE,PHI
      vec3 B_vec;
      vec3 pos;
      pos=Vec_Handling::RTHEPHI2cart_vec(dist, THE, PHI);

      pos.x=pos.x+SunPosition.x;
      pos.y=pos.y+SunPosition.y;
      pos.z=pos.z+SunPosition.z;

      // this call assumes pos is galactocentric
      B_vec=mag.return_B_cart(pos);
      B_per=Vec_Handling::return_perp2LOS(B_vec, THE,PHI);
      intr_pol_ang=Vec_Handling::return_intr_pol_ang(B_vec, THE, PHI);

      if (dist==min_shell_radius){continue;}

      // in the local coordinate system on the plane of
      // the sky, x and y are defined as:
      // x axis points to the right of the GC (i.e. in
      // decreasing l direction, i.e. negative phi direction)
      // , the y axis points towards the galactic north
      // (i.e. negative theta direction).

      observables.Bx+=B_per*cos(intr_pol_ang);
      observables.By+=B_per*sin(intr_pol_ang);
    }
}





void Integrator::get_synchrotron_emissivity (const double cre_C, double distribution_index_p, const double Bperp, const double omega, double& P_perpendicular_value, double& P_parallel_value) const
{
  double common_part=0;
  // The value of Bsinalpha should be abs(Bsinalpha).
  // There is no sense in having P_perpendicular_value
  // and P_parallel_value <0.
  double abs_Bperp=abs(Bperp);

  // notice that we assume m=mass of electron.

  if (Bperp==0 && distribution_index_p<=1) {cout << "Warning Bperp==0, common_part in emissivity.get_P might blow up for distibution_indexes <= 1" << endl;}

  const double common_factor=(sqrt(3.))*(CGS_U_qe*CGS_U_qe*CGS_U_qe)/(4.*CGS_U_pi*CGS_U_MEC2);
  const double common_factor2=(2./3.)*(CGS_U_ME*CGS_U_C_light)/CGS_U_qe;

  common_part=(1./2.)*common_factor*pow((abs_Bperp),((distribution_index_p+1.)/2.))*pow((common_factor2*omega),(-(distribution_index_p-1.)/2.));


  // taking into account the spatial component of the
  // cre distribution:
  common_part*=cre_C;


  if (common_part<0) {cout << "Error, common_part in Integrator::get_synchrotron_emissivity(...) is <0 : " << common_part << endl;}

  const double gamma_function_1=gsl_sf_gamma(distribution_index_p/4.-1./12.);
  const double gamma_function_2=(pow(2.,(distribution_index_p+1.)/2.)/(distribution_index_p+1.))*gsl_sf_gamma(distribution_index_p/4.+19./12.);
  const double gamma_function_3=pow(2.,(distribution_index_p-3.)/2.)*gsl_sf_gamma(distribution_index_p/4.+7./12.);
  const double perpendicular_gamma_part=gamma_function_1*(gamma_function_2+gamma_function_3);
  const double parallel_gamma_part=gamma_function_1*(gamma_function_2-gamma_function_3);

  if (perpendicular_gamma_part<0 || parallel_gamma_part <0) {cout << "Error, perpendicular_gamma_part or parallel_gamma_part <0 in Integrator::get_synchrotron_emissivity(...)" << perpendicular_gamma_part << " " << parallel_gamma_part << endl;}

  P_perpendicular_value=(1./(4.*CGS_U_pi))*common_part*perpendicular_gamma_part;
  P_parallel_value=(1./(4.*CGS_U_pi))*common_part*parallel_gamma_part;

}


// similar to computing QU maps
void Integrator::get_deflection_map(B_field &mag)
{
  cout << " get_deflection_map(...) initiated " << endl;

  cout << " Allocating memory for obs shell " << endl;
  obs_defl_Bx_map.SetNside(obs_NSIDE,NEST);
  obs_defl_By_map.SetNside(obs_NSIDE,NEST);
  cout << " Zeroing obs shell before computation begin " << endl;
  for(int ipix=0;ipix<obs_defl_Bx_map.Npix();ipix++)
    {
      obs_defl_Bx_map[ipix]=0.;
      obs_defl_By_map[ipix]=0.;
    }
  cout << " Zeroed!" << endl;

  // This for loop goes trough all shells, computes them for the proper
  // shell radiae etc...
  for (int current_shell=1;current_shell<(total_shell+1);current_shell++)
    {
      Log(" shell number: "); cout << current_shell << endl;
      cout << " NSIDE:        " << get_shell_NSIDE(current_shell) << endl;

      double min_shell_radius;
      double max_shell_radius;
      if(current_shell==1)
        {
          min_shell_radius=0.;
        }
      else
        {
          min_shell_radius=get_min_shell_radius(current_shell);
        }
      max_shell_radius=get_max_shell_radius(current_shell);

      cout << " Current shell interval is [" << min_shell_radius/CGS_U_kpc << " , " << max_shell_radius/CGS_U_kpc << "]" << endl;

      cout << " Declaring deflection temporary Healpix_Maps " << endl;
      Healpix_Map<double> current_defl_Bx_map;
      Healpix_Map<double> current_defl_By_map;
      cout << " temporary Healpix_Maps declared " << endl;


      //------------------------------------------------------------------
      // Now memory for the maps is allocated.
      //------------------------------------------------------------------

      cout << " Allocating I,Q,U temporary map space " << endl;

      current_defl_Bx_map.SetNside(get_shell_NSIDE(current_shell),NEST);
      current_defl_By_map.SetNside(get_shell_NSIDE(current_shell),NEST);
      cout << " allocation done " << endl;
      cout << " Zeroing temporary maps " << endl;
      for(int ipix=0;ipix<current_defl_Bx_map.Npix(); ipix++)
        {
          current_defl_Bx_map[ipix]=0.;
          current_defl_By_map[ipix]=0.;
        }
      cout << " temporary maps zeroed " << endl;

      // Done! Memory for maps is allocated and they are zeroed.
      //------------------------------------------------------------------

      //------------------------------------------------------------------
      // Setting current_map_Npix.
      int current_map_Npix=0;
      current_map_Npix=current_defl_Bx_map.Npix();

      // now, using the current_map_Npix, we cycle
      // over all to be integrated pixels.

      struct_defl_obs observables;

      for(int ipix=0;ipix<current_map_Npix; ipix++)
        {

          pointing ptg;
          ptg=current_defl_Bx_map.pix2ang(ipix);

          radial_integration(min_shell_radius,max_shell_radius, current_shell, ipix, observables, mag);

          current_defl_Bx_map[ipix]=observables.Bx;
          current_defl_By_map[ipix]=observables.By;
        }

      //------------------------------------------------------------------
      // Adds values to obs maps.

      // Add the values to obs_defl_int/ori_map using whatever interpolation scheme.
      if(current_shell<=obs_shell)
        {
          for(int ipix=0;ipix<obs_defl_Bx_map.Npix();ipix++)
            {
              pointing ptg;
              ptg=obs_defl_Bx_map.pix2ang(ipix);
              obs_defl_Bx_map[ipix]+=current_defl_Bx_map.interpolated_value(ptg);
              obs_defl_By_map[ipix]+=current_defl_By_map.interpolated_value(ptg);

              // Alternatively
              //          obs_defl_Bx_map[ipix]+=current_defl_Bx_map[down_ipix(current_shell,ipix)];
            }
        }
      else if(current_shell>obs_shell)
        {
          for(int ipix=0;ipix<obs_defl_Bx_map.Npix();ipix++)
            {
              /*
                pointing ptg;
                ptg=obs_defl_Bx_map.pix2ang(ipix);
                obs_defl_Bx_map[ipix]+=current_defl_Bx_map.interpolated_value(ptg);
                obs_defl_By_map[ipix]+=current_defl_By_map.interpolated_value(ptg);
              */
              // Alternatively
              int level_factor=1;
              int ipix_manipulation=ipix;
              for(int n=0;n<(current_shell-obs_shell);n++)
                {
                  level_factor*=4;
                  ipix_manipulation*=4;
                }
              for(int n=ipix_manipulation;n<(ipix_manipulation+level_factor);n++)
                {
                  obs_defl_Bx_map[ipix]+=current_defl_Bx_map[n]/level_factor;
                  obs_defl_By_map[ipix]+=current_defl_By_map[n]/level_factor;
                }
            }
        }
      else {cerr << " Erro when summing up shells! " << endl; exit(1);}

      // Done adding values to respective obs shells
      //------------------------------------------------------------------

    }
  cout << " get_deflection_map(...) done " << endl;
}

void Integrator::save_defl_map2file(bool do_rigidity, double rigidity)
{
  cout << " preparing units for output " << endl;
  Healpix_Map<double> temp_Bx_map;
  Healpix_Map<double> temp_By_map;

  temp_Bx_map.SetNside(obs_defl_Bx_map.Nside(), obs_defl_Bx_map.Scheme());
  temp_By_map.SetNside(obs_defl_By_map.Nside(), obs_defl_By_map.Scheme());

  double unit_factor;
  if(do_rigidity){unit_factor=1.;}
  else{unit_factor=CGS_U_Gauss*CGS_U_m;}

  for(int ipix=0;ipix<obs_defl_Bx_map.Npix();ipix++)
    {
      temp_Bx_map[ipix]=obs_defl_Bx_map[ipix]*delta_R/unit_factor;
      temp_By_map[ipix]=obs_defl_By_map[ipix]*delta_R/unit_factor;

      double rad2deg=180./CGS_U_pi;
      obs_defl_Bx_map[ipix]=rad2deg*sqrt(temp_Bx_map[ipix]*temp_Bx_map[ipix]+temp_By_map[ipix]*temp_By_map[ipix]);
      if(do_rigidity){obs_defl_Bx_map[ipix]/=rigidity;}
      obs_defl_By_map[ipix]=atan2(temp_By_map[ipix], temp_Bx_map[ipix]);
    }

  cout << " done " << endl;
  {
    fitshandle out;
    out.create(string("!")+defl_int_file_name);
    write_Healpix_map_to_fits(out, obs_defl_Bx_map, PLANCK_FLOAT64);
    out.close();
    if(do_rigidity) {cout << "do_rigidity is true. Maps computed for rigidity=" << rigidity/(CGS_U_erg/CGS_U_esu) << "erg/esu" << endl;}
    else{cout << " defl int map in units of Gauss*meter " << endl;}
    cout << " deflection intensity map saved to " << defl_int_file_name << endl;
  }
  {
    fitshandle out2;
    out2.create(string("!")+defl_ori_file_name);
    write_Healpix_map_to_fits(out2, obs_defl_By_map, PLANCK_FLOAT64);
    out2.close();
    cout << " deflection orientation map saved to " << defl_ori_file_name << endl;
  }
}
