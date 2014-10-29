/*
The "Hammurabi" code simulates Galactic synchrotron emission at all frequencies.
Copyright (C) 2006,2007,2008,2012 Andre Waelkens, Tess Jaffe, Martin Reinecke, Torsten Ensslin

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

#ifndef HAMMURABI_B_FIELD_H
#define HAMMURABI_B_FIELD_H
#include <ctime>
#include "fftw3.h"
#include "planck_rng.h"
#include "tess_tools.h"
#include "paramfile.h"
#include "proto_class_Integrator.h"

// Tess' new version of proto_class_B_field.h

/* Class B_field meant to be used as follows:

   Integrator integ(params);
   B_field magfield(params);
   ...
   integ.calculate(magfield, n_e_thermal, ncre, map_i, map_q, map_u, map_rm);

   If params indicate only a regular component, it is then calculated
   analytically whenever the Integrator calls it as
   b3=magfield.return_B_cart(coords);

   If params indicate a random component, it is computed and stored in
   a grid at the time the constructor is called with the input
   parameters.  Then, a call to return_B_cart adds the random
   component to the regular component; the former is interpolated from
   the grid and the latter computed analytically at the input
   coordinates.

   To use a model with different parameters in a single code, call the
   constructor with empty brackets
   B_field magfield();
   and then call the setup function explicitely based on whatever
   you want, e.g. looping over pitch angles and calling
   setup_field2(debug,  pitch,  r0,  z0,  b0,  phi0);

   For multiple realisations of the random component, after the first
   run (set up either way as above), clean it and redo the random:
   magfield.cleanRandom();
   magfield.fillRandom();

*/

// Needed so that we can use Integrator's private data below in SanityCheck
class Integrator;

class B_field
{
 private:

  int bfield_quiet;
  long b_cnt;
  bool bfield_debug;
  double gc_r_max,ec_r_max;
  bool grid_interp;
  double grid_lx, grid_ly, grid_lz;
  int grid_nx, grid_ny, grid_nz;

  // Sun position in galactocentric cartesian
  // coordinates
  vec3 SunPosition;

  // Random field:
  bool dorand;
  double *b_of_k_x;
  double *b_of_k_y;
  double *b_of_k_z;
  double bran_alpha;
  double bran_cutoff_kpc;
  int bran_seed;
  double bran_random_rms;
  double bran_random_c0;
  std::string bran_file;
  double bran_rmax;
  double bran_transform;
  double bran_mem_limit;

  // External input file
  std::string bran_inp_file;

  // Random numbers:
  planck_rng bran_rng;

  // For writing total field
  std::string btot_file;
  double tlon;
  double tlat;
  bool Bonly;
	

  // field1:
  double b1_b0, b1_Psi_0,b1_Psi_1,b1_Xsi_0;


  // field2:
  double b2_pitch_deg;
  double b2_r0;
  double b2_z0;
  double b2_b0;
  double b2_phi0;

  // field3:
  // ASS+RING from Sun et al. A&A V.477 2008
  double b3_B0;
  // Note, redundant with SunPosition. However
  // might want to tune it independently anyways.
  double b3_Rsun;
  double b3_R0;
  double b3_z0;
  double b3_Rc;
  double b3_Bc;
  double b3_pitch_deg;

  double b3H_B0;
  double b3H_z0;
  double b3H_z1a;
  double b3H_z1b;
  double b3H_R0;

  // field4:
  // HMR field from Kachelriess et al. APh 2007
  double b4_Rsun;
  double b4_z1;
  double b4_z2;
  double b4_r1;
  double b4_p;
  double b4_epsilon_0;

  // field5:
  // TT field from Kachelriess et al. APh 2007
  double b5_b0;
  double b5_Rsun;
  double b5_r_min;
  double b5_d;
  double b5_z0;
  double b5_p;

  // field6:
  // read from file
  std::string breg_inp_file;
  double *b_file_x;
  double *b_file_y;
  double *b_file_z;

  // field7
  // Random field generation, Sean Quinn Aug 12 2014
  double b7_kmin;
  double b7_Lc;
  unsigned int b7_N_m;
  double b7_nf;
  double b7_dkn_const;
  //long long unsigned int b7_start_seed; Taken out on 10/13/14

  //field8
  // Random field computation, Sean Quinn Aug 15 2014
  double b8_b1;
  double b8_b2;
  double b8_b3;
  double b8_b4;
  double b8_b5;
  double b8_b6;
  double b8_b7;
  double b8_b8;
  double b8_b9;
  double b8_b10;
  double b8_b11;
  double b8_b12;
  double b8_b13;
  double b8_b14;
  double b8_b15;
  double b8_b16;
  double b8_b17;
  double b8_b18;
  double b8_b19;
  double b8_b20;
  double br8_sqrt_beta;
  std::vector<double> b8_bstrix;
  std::vector<double> b8_bstriy;
  std::vector<double> b8_bstriz;

  //field9
  // Random field computation, Sean Quinn Aug 15 2014
  double b9_b1;
  double b9_b2;
  double b9_b3;
  double b9_b4;
  double b9_b5;
  double b9_b6;
  double b9_b7;
  double b9_b8;
  double b9_b9;
  double b9_b10;
  double b9_b11;
  double b9_b12;
  double b9_b13;
  double b9_b14;
  double b9_b15;
  double b9_b16;
  double b9_b17;
  double b9_b18;
  double b9_b19;
  double b9_b20;
  double br9_b1;
  double br9_b2;
  double br9_b3;
  double br9_b4;
  double br9_b5;
  double br9_b6;
  double br9_b7;
  double br9_b8;
  double br9_b9;
  double br9_b10;
  double br9_b11;
  double br9_b12;
  double br9_b13;
  double br9_kappa;
  double br9_sqrt_beta;
  std::vector<double> b9_brandx;
  std::vector<double> b9_brandy;
  std::vector<double> b9_brandz;
  std::vector<double> b9_bstrix;
  std::vector<double> b9_bstriy;
  std::vector<double> b9_bstriz;

  //field10 
  // Random field computation, Sean Quinn Aug 15 2014
  std::vector<double> b10_brandx;
  std::vector<double> b10_brandy;
  std::vector<double> b10_brandz;
  double b10_brms;

  // field11:
  double b11_b1;
  double b11_b2;
  double b11_b3;
  double b11_b4;
  double b11_b5;
  double b11_b6;
  double b11_b7; 
  double b11_b8;
  double b11_b9;
  double b11_b10;
  double b11_b11; 
  double b11_b12;
  double b11_b13;
  double b11_b14;
  double b11_b15;
  double b11_b16;
  double b11_b17;
  double b11_b18;
  double b11_b19; 
  double b11_b20;
  double b11_b21;
  double b11_b22;
  double b11_b23; 
  double b11_b24;
  double b11_b25;
  double b11_b26;
  double b11_b27;
  double b11_b28;
  double b11_b29;
  double b11_b30;

  //field12
  double b12_bxc;
  double b12_byc;
  double b12_bzc;
  std::vector<double> b12_brandx;
  std::vector<double> b12_brandy;
  std::vector<double> b12_brandz;
  double b12_brms;

  // Your field here:
  //vec3 field10(vec3 coords);

  vec3 field1(vec3 coords);
  vec3 field2(vec3 coords);
  vec3 field3(vec3 coords);
  vec3 field4(vec3 coords);
  vec3 field5(vec3 coords);
  vec3 field6(vec3 coords);
  vec3 field7(vec3 coords); //Added by SPQ Aug 12 2014
  vec3 field8(vec3 coords); //Added by SPQ Aug 15 2014
  vec3 field9(vec3 coords); //Added by SPQ Aug 15 2014
  vec3 field10(vec3 coords); //Added by SPQ Aug 15 2014
  vec3 field11(vec3 coords); //Added by SPQ Oct 20 2014
  vec3 field12(vec3 coords); //Added by SPQ Oct 28 2014

  void read_B_params(paramfile &params);

 public:


        int bfield_type;

        B_field(paramfile &params);

        B_field();

	~B_field(){
	  if (dorand) {
	    fftw_free(b_of_k_x);
	    fftw_free(b_of_k_y);
	    fftw_free(b_of_k_z);
	  }
	  if (bfield_type==6) {
	    delete [] b_file_x;
	    delete [] b_file_y;
	    delete [] b_file_z;
	  }	
	}

        bool return_dorand(void) {return dorand;}

        // Functions to read in parameters for each particular field type:

	void setup_field1(double b0,double psi0,double psi1,double xsi0);
	void setup_field2(double pitch_deg, double r0, double z0, double b0, double phi0);
	void setup_field3(double pitch_deg, double B0, double Rsun, double R0, double z0, double Rc, double Bc, double H_B0, double H_z0, double H_z1a, double H_z1b, double H_R0);
	void setup_field4(double Rsun, double z1, double z2, double r1, double p, double epsilon_0);
	void setup_field5(double b0, double Rsun, double r_min, double d, double z0, double p);
	void setup_field6(double lx,double ly,double lz,int nx,int ny,int nz, double tlon, double tlat, bool interp, std::string breg_inp_file);
    void setup_field7(double kmin,double Lc,unsigned int num_modes,double nf,double dkn_const);//,long long unsigned int start_seed); //Added by SPQ Aug 12 2014, removed by SPQ 10/13/14
    void setup_field8(double b1, double b2, double b3, double b4, double b5, 
                            double b6, double b7, double b8, double b9, double b10, 
	                        double b11, double b12, double b13, double b14, double b15, 
                            double b16, double b17, double b18, double b19, double b20, 
                            std::vector<double> bstrix, std::vector<double> bstriy, std::vector<double> bstriz, 
                            double sqrt_beta); //Added by SPQ Aug 15 2014
    void setup_field9(double b1, double b2, double b3, double b4, double b5, 
                            double b6, double b7, double b8, double b9, double b10, 
	                        double b11, double b12, double b13, double b14, double b15, 
                            double b16, double b17, double b18, double b19, double b20,
                            double br1, double br2, double br3, double br4, double br5, double br6,
                            double br7, double br8, double br9, double br10, double br11,
                            double br12, double br13, std::vector<double> brandx, 
                            std::vector<double> brandy, std::vector<double> brandz, double kappa, double sqrt_beta,
                            std::vector<double> bstrix, std::vector<double> bstriy, std::vector<double> bstriz); //Added by SPQ Aug 16 2014
    void setup_field10(std::vector<double> brandx, std::vector<double> brandy, std::vector<double> brandz, double brms_factor); //Added by SPQ Aug 16 2014
    void setup_field11(double b1, double b2, double b3, double b4, double b5, double b6, double b7, double b8, 
                            double b9, double b10,  double b11, double b12, double b13, double b14, double b15, 
                            double b16, double b17, double b18, double b19, double b20, double b21, double b22,  
                            double b23, double b24, double b25, double b26, double b27, double b28, double b29, double b30); //Added by SPQ Oct 20 2014
    void setup_field12(double bx_const, double by_const, double bz_const,std::vector<double> brandx, 
                        std::vector<double> brandy, std::vector<double> brandz, double brms_factor); // Added Oct 28 2014 SPQ
	void setup_random(double alpha,double cutoff_kpc,int seed, double random_rms, double c0, double lx_kpc,double ly_kpc, double lz_kpc, int nx, int ny, int nz, std::string bran_file, double rmax_ran, double tlon, double tlat, bool interp, double mem_lim,bool debug, std::string bran_inp_file);


	void setup_writer(std::string btotfile,int nx,int ny,int nz,double lx,double ly,double lz, double tlon, double tlat, bool interp, bool B_only);
	bool write_B_total();

	// Your field here:
	// void setup_field10();


        // Parameters for random component set by setup_random,
        // setup_random_external above

        // Optionally do *not* alloc fillRandom if doing a second
        // simulation of the same dimensions; waste of time.
        void fillRandom( bool do_alloc = true);
        void fillRandom(bool do_alloc, std::string infile);

        // Free random component.
        void cleanRandom();

        vec3 return_B_cart( vec3 coords);

        vec3 return_Breg_cart( vec3 coords);

        vec3 return_Bran_cart( vec3 coords);



        void SanityCheck(Integrator &intor, double lmin=0, double lmax=0);
        void SanityCheck();


}; // class B_field
#endif
