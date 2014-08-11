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

#include <fstream>
#include "proto_class_List.h"
#include "hammurabi.h"

using namespace std;

List::List(paramfile &params, int shell_difference)
{

  non_eg_pol_sources=params.find<bool>("non_eg_pol_sources_flag", false);
  item_obs_freq=params.find<bool>("item_obs_freq_flag", false);
  in_obs_freq=CGS_U_GHz*params.find<double>("obs_freq_GHz",9999.);
  list_use_lonlat_range = params.find<bool>("list_use_lonlat_range",false);
  lon_min=params.find<double>("lon_min",0.);	
  lon_max=params.find<double>("lon_max",359.99);
  lat_min=params.find<double>("lat_min",-90.);	
  lat_max=params.find<double>("lat_max",90.);	
  nside=params.find<int>("obs_NSIDE",2.);

  MAX_RADIUS=CGS_U_kpc*params.find<double>("max_radius", 32);
	
  list_in_file_name=params.find<std::string>("list_in_file_name", "list_in.txt");
  list_out_file_name=params.find<std::string>("list_out_file_name", "list_out");

  pol_map_in_kelvin_flag=params.find<bool>(" pol_map_in_kelvin_flag",true);

  zero();
}


List::List(int length, bool neps, bool iof, double obs_freq, double rmax, int shell_difference, std::string in_file_name, std::string out_file_name, bool kelvin_flag)
{
  cout << " List constructor " << endl;
  cout << " sono qui " << endl;
  data_obs_shell.alloc(length);
  long data_total_shell_size;

  non_eg_pol_sources=neps;
  item_obs_freq=iof;
  in_obs_freq=obs_freq;
  MAX_RADIUS=CGS_U_kpc*rmax;

  total_to_obs_factor=  ((int) pow(4.,(double) shell_difference));
  data_total_shell_size=data_obs_shell.size()*total_to_obs_factor;
  cout << " la data_total_shell size : " << data_total_shell_size << endl;
  data_total_shell.alloc(data_total_shell_size);
  list_in_file_name=in_file_name;
  list_out_file_name=out_file_name;

  pol_map_in_kelvin_flag=kelvin_flag;

  zero();
}

void List::create_list(paramfile &params, int shell_difference)
{
	//read_in_list();
	
	//if (list_in_file_name.compare("") == 0 &&  list_use_lonlat_range==true) {
	if (list_use_lonlat_range==true) {
				
		int count;

		cout << " Longitude range: " << lon_min << " to " << lon_max << endl;
		cout << " Latitude range: " << lat_min << " to " << lat_max << endl;
		count = count_lonlat_region(nside, lat_min, lat_max, lon_min, lon_max);
		
		cout << " Number of pixels: " << count << endl;	
		data_obs_shell.alloc(count);
		long data_total_shell_size;		
		total_to_obs_factor=  ((int) pow(4.,(double) shell_difference));
		data_total_shell_size=data_obs_shell.size()*total_to_obs_factor;
		cout << " la data_total_shell size : " << data_total_shell_size << endl;
		data_total_shell.alloc(data_total_shell_size);
		
		arr<long> pixlist(count);
		get_lonlat_region(nside, lat_min, lat_max, lon_min, lon_max, pixlist);
		fill_list(pixlist, in_obs_freq);	
		
		
	}
	else {

		data_obs_shell.alloc(params.find<int>("list_length",10));
		long data_total_shell_size;		
		total_to_obs_factor=  ((int) pow(4.,(double) shell_difference));
		data_total_shell_size=data_obs_shell.size()*total_to_obs_factor;
		cout << " la data_total_shell size : " << data_total_shell_size << endl;
		data_total_shell.alloc(data_total_shell_size);
		read_in_list();
	}
	
}

void List::read_in_list()
{
  cout << " Reading in list: " << list_in_file_name << endl;
  ifstream in_file (list_in_file_name.c_str());
  planck_assert(in_file, "error opening file");
  for(long n=0;n<data_obs_shell.size();n++)
    {
      in_file >> data_obs_shell[n].ipix;
      // using list specified obs frequencies
      if(item_obs_freq)
        {
          in_file >> data_obs_shell[n].obs_freq;
          data_obs_shell[n].obs_freq*=CGS_U_GHz;
        }
      else
        {
          data_obs_shell[n].obs_freq=in_obs_freq;
        }

      // checking if polarization sources are
      // not extra-galactic
      if(non_eg_pol_sources)
        {
          in_file >> data_obs_shell[n].source_dist;
          data_obs_shell[n].source_dist*=CGS_U_kpc;
          if(data_obs_shell[n].source_dist>MAX_RADIUS)
            {
              cerr << " Warning, source_dist is very large! " << endl;
              cerr << " Source seems to be outside of the Galaxy. " << endl;
              cerr << " source_dist = " << data_obs_shell[n].source_dist/CGS_U_kpc << endl;
              cerr << " stopping the code " << endl;
              exit(1);
            }
        }
      else
        {
          data_obs_shell[n].source_dist=MAX_RADIUS;
        }

      for (long m=0;m<total_to_obs_factor;m++) {
        data_total_shell[n*total_to_obs_factor+m].ipix=data_obs_shell[n].ipix*total_to_obs_factor+m;
        data_total_shell[n*total_to_obs_factor+m].obs_freq=data_obs_shell[n].obs_freq;
      }
    }

  cout << " List is read! " << endl;
}

void List::fill_list( arr<long> &listarr, double obs_freq){
  arr<double>freq(listarr.size());  freq.fill(obs_freq);
  arr<double>dist(listarr.size());  dist.fill(MAX_RADIUS);
  fill_list(listarr,freq,dist);
}
void List::fill_list( arr<long> &listarr, arr<double> &obs_freq, arr<double> &dist){

  if (listarr.size() != data_obs_shell.size()) {
    cout<<"ERROR:  input list array not the right size ("<<listarr.size()<<") for given List length ("<<data_obs_shell.size()<<")\n";
    std::exit(1);
  }
  cout << " Filling list from input array."<<endl;
  for (long n=0;n<data_obs_shell.size();++n) {
    data_obs_shell[n].ipix=listarr[n];
    data_obs_shell[n].obs_freq=obs_freq[n];
    data_obs_shell[n].source_dist=dist[n];
    for (long m=0;m<total_to_obs_factor;m++) {
      data_total_shell[n*total_to_obs_factor+m].ipix=data_obs_shell[n].ipix*total_to_obs_factor+m;
      data_total_shell[n*total_to_obs_factor+m].obs_freq=data_obs_shell[n].obs_freq;
      data_total_shell[n*total_to_obs_factor+m].source_dist=data_obs_shell[n].source_dist;
    }
  }
  cout << " List is filled." << endl;
}


// Saves all 4 outputs (I,Q,U and RM) from List to file.
void List::save_list2file( std::string filename)
{
  cout << " void List::save_list2file() initiated " << endl;

  bool same = (data_obs_shell.size() == data_total_shell.size());

  string extension = ".txt";
  string infoext = "_info.txt";
  string writefile = (filename.compare("")!=0) ? filename+extension : list_out_file_name+extension;
  string info_writefile = (filename.compare("")!=0) ? filename+infoext : list_out_file_name+infoext;
  ofstream out (writefile.c_str());
  ofstream out_info (info_writefile.c_str());
  planck_assert(out, "error opening file");
  planck_assert(out_info, "error opening file");


  double unit_corr_factor=(pol_map_in_kelvin_flag) ? CGS_U_cm*CGS_U_cm/CGS_U_erg : 1.0; // Otherwise we display values in code internal untis

  cout << " Writing info in " << info_writefile << endl;
  out_info << " ipix  " ;
  out_info << " freq [GHz] " ;
  out_info << " I     " << " Q     " << " U     " << " Int Dist [kpc] " ;
  if(same)
    {
      out_info << "  RM [rad/m^2]";
      out_info << "  DM [pc/cm^3]" << endl;
    }
  else
    {
      out_info << endl;
      out_info << " And for total_data: " << endl;
      out_info << " ipix " ;
      out_info << " freq " ;
      out_info << " I     " << " Q     " << " U     " << " Int Dist " << "   RM " << " DM " << endl;;
    }


  cout << " info written " << endl;

  for(long n=0;n<data_obs_shell.size();n++)
    {
      out << data_obs_shell[n].ipix;
      out << " " << data_obs_shell[n].obs_freq/CGS_U_GHz;
      out << " " << data_obs_shell[n].observables.I*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_obs_shell[n].observables.Q*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_obs_shell[n].observables.U*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_obs_shell[n].source_dist/CGS_U_kpc;

    if (same)
      {
        double rad_per_sq_m = 1./(CGS_U_m*CGS_U_m);
        double pc_per_ccm = CGS_U_pc/CGS_U_ccm;
        out << " " << data_obs_shell[n].observables.RM/rad_per_sq_m;
        out << " " << data_obs_shell[n].observables.DM/pc_per_ccm << endl;
      }
    else out  << endl;
    }

  cout << " List result for obs_shell saved to: " << writefile << endl;

  if (same) return;

  writefile = (filename.compare("")!=0) ? filename+"_total_data"+extension : list_out_file_name+"_total_data"+extension;

  ofstream out2 (writefile.c_str());
  planck_assert(out2, "error opening file");


  for(long n=0;n<data_total_shell.size();++n)
    {
      out2 << data_total_shell[n].ipix;
      out2 << " " << data_total_shell[n].obs_freq/CGS_U_GHz;
      out2 << " " << data_total_shell[n].observables.I*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_total_shell[n].observables.Q*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_total_shell[n].observables.U*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << " " << data_total_shell[n].source_dist/CGS_U_kpc << " " << data_total_shell[n].observables.RM*(CGS_U_m*CGS_U_m) << " " << data_total_shell[n].observables.DM*(CGS_U_ccm/CGS_U_pc) << endl;
    }

  cout << " List result for total_shell saved to: " << writefile << endl;

  cout << " void List::save_list2file() finished " << endl;
}

void List::save_latlon( std::string filename) 
{
	int pixel;
	double lat, lon;
	double unit_corr_factor=(pol_map_in_kelvin_flag) ? CGS_U_cm*CGS_U_cm/CGS_U_erg : 1.0; // Otherwise we display values in code internal untis
	double rad_per_sq_m = 1./(CGS_U_m*CGS_U_m);
	
	cout << " void List::save_latlon() initiated " << endl;
	
	//Stokes I	
	
	string extension_I = "_I.txt";

	string writefile_I = (filename.compare("")!=0) ? filename+extension_I : list_out_file_name+extension_I;
	ofstream out_I (writefile_I.c_str());
	planck_assert(out_I, "error opening file");

	//Stokes Q
	
	string extension_Q = "_Q.txt";
	
	string writefile_Q = (filename.compare("")!=0) ? filename+extension_Q : list_out_file_name+extension_Q;
	ofstream out_Q (writefile_Q.c_str());
	planck_assert(out_Q, "error opening file");	
	
	//Stokes U
	
	string extension_U = "_U.txt";
	
	string writefile_U = (filename.compare("")!=0) ? filename+extension_U : list_out_file_name+extension_U;
	ofstream out_U (writefile_U.c_str());
	planck_assert(out_U, "error opening file");
	
	//Rotation Measure
	
	string extension_RM = "_RM.txt";
	
	string writefile_RM = (filename.compare("")!=0) ? filename+extension_RM : list_out_file_name+extension_RM;
	ofstream out_RM (writefile_RM.c_str());
	planck_assert(out_RM, "error opening file");
		
	for(long n=0;n<data_obs_shell.size();n++)
    {
		
		pixel=data_obs_shell[n].ipix;
		lat = get_lat(nside, pixel);
		lon = get_lon(nside, pixel);
		out_I << lat;
		out_I << "," << lon;
		out_I << "," << data_obs_shell[n].observables.I*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << endl;
		
		out_Q << lat;
		out_Q << "," << lon;
		out_Q << "," << data_obs_shell[n].observables.Q*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << endl;
		
		out_U << lat;
		out_U << "," << lon;
		out_U << "," << data_obs_shell[n].observables.U*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad)) << endl;
		
		out_RM << lat;
		out_RM << "," << lon;
		out_RM << "," << data_obs_shell[n].observables.RM/rad_per_sq_m << endl;
		
		
    }
	
	cout << " List result for obs_shell I saved to: " << writefile_I << endl;	
	cout << " List result for obs_shell Q saved to: " << writefile_Q << endl;
	cout << " List result for obs_shell U saved to: " << writefile_U << endl;
	cout << " List result for obs_shell RM saved to: " << writefile_RM << endl;
	
}



void List::return_list_physical( arr<long> &pixlist, arr<struct_observables> &obslist) {
    if (pixlist.size() != data_obs_shell.size()) { pixlist.alloc(data_obs_shell.size());}
    if (obslist.size() != data_obs_shell.size()) { obslist.alloc(data_obs_shell.size());}

  double unit_corr_factor=(pol_map_in_kelvin_flag) ? CGS_U_cm*CGS_U_cm/CGS_U_erg : 1.0; // Otherwise we display values in code internal untis
  double rad_per_sq_m = 1./(CGS_U_m*CGS_U_m);
  double pc_per_ccm = CGS_U_pc/CGS_U_ccm;


    for (long n=0;n<data_obs_shell.size();n++) {
      pixlist[n]=data_obs_shell[n].ipix;
      //      obslist[n]=data_obs_shell[n].observables;
      obslist[n].I=data_obs_shell[n].observables.I*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
      obslist[n].Q=data_obs_shell[n].observables.Q*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
      obslist[n].U=data_obs_shell[n].observables.U*(1./unit_corr_factor)/(CGS_U_erg/(CGS_U_sec*CGS_U_Hz*CGS_U_cm*CGS_U_cm*CGS_U_sterad));
      obslist[n].RM=data_obs_shell[n].observables.RM/rad_per_sq_m;
      obslist[n].DM=data_obs_shell[n].observables.DM/pc_per_ccm;
    }
}


// loads ext_freq into data_obs_shell[n].obs_freq
void List::load_full_list_freq(double ext_freq)
{
  if(item_obs_freq)
    {
      std::cerr << " overwriting list defined freq! Freq for this run is: "
                << ext_freq <<  std::endl;
    }
  for(int n=0;n<data_obs_shell.size();n++)
    {data_obs_shell[n].obs_freq=ext_freq;}
  for(int n=0;n<data_total_shell.size();n++)
    {data_total_shell[n].obs_freq=ext_freq;}
}
