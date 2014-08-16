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

#include <hammurabi.h>
#include <proto_class_B_field2.h>
#include <vec3.h>
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/algorithm/string.hpp>
#include <gsl/gsl_math.h>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sstream>

using namespace std;

//---------------------------------------------------------------
// Construct and set up from parameter file
//---------------------------------------------------------------
B_field::B_field( paramfile &params){
  read_B_params(params);
  b_cnt=0;
}
//---------------------------------------------------------------


//---------------------------------------------------------------
// Construct empty field to be set up by a subsequent call to setup_field1...
//---------------------------------------------------------------
B_field::B_field(){
  b_cnt=0;
  bfield_debug=false;
  bfield_type=1;
  bfield_quiet=0;
  dorand=false;
  // Sun position in GC cart coordinates default setting:
  SunPosition.x=-CGS_U_kpc*8.5;
  SunPosition.y=CGS_U_kpc*0.;
  SunPosition.z=CGS_U_kpc*0.;
}
//---------------------------------------------------------------

  unsigned int nside_val; // Global value added by SPQ. Required to compute size of rand field arrays. Defaults to Nside=16.
  unsigned int brfa_len; // Global value added by SPQ. Required to compute size of rand field arrays. Defaults to Nside=16.

//---------------------------------------------------------------
//  Read from parameter file and set up given field type.
//---------------------------------------------------------------
void B_field::read_B_params(paramfile &params) {
  // If you skip this function, need to run first the generic_setup
  // and then one of the setups plus maybe the random setup.

  nside_val = params.find<unsigned int>("obs_NSIDE",16); // Added by SPQ Aug 12 2014
  //First step: determine number of Hammurabi steps through grid to build rand field array
  if(nside_val==1) brfa_len=710;
  if(nside_val==2) brfa_len=2934;
  if(nside_val==4) brfa_len=11790;
  if(nside_val==8) brfa_len=47316;
  if(nside_val==16) brfa_len=189276;
  if(nside_val==32) brfa_len=757472;
  if(nside_val==64) brfa_len=3030270;
  if(nside_val==128) brfa_len=12118494;
  if(nside_val==256) brfa_len=48469568;
  //

  bfield_debug=params.find<bool>("B_field_debug",false);
  bfield_type=params.find<int>("B_field_type",3);
  bfield_quiet=params.find<int>("B_field_quiet",0);
  gc_r_max=params.find<double>("B_field_r_max",30.)*CGS_U_kpc;
  ec_r_max=params.find<double>("B_field_rmax_sun",30.)*CGS_U_kpc;

  // Sun position in GC cartesian coordinates.u
  vec3 Sun;
  SunPosition.x=CGS_U_kpc*params.find<double>("SunPosX", -8.5);
  SunPosition.y=CGS_U_kpc*params.find<double>("SunPosY", 0.);
  SunPosition.z=CGS_U_kpc*params.find<double>("SunPosZ", 0.);

  dorand=params.find<bool>("B_field_do_random",true);
    
  // Common parameters for handling grids
  std::string btotfile=params.find<string>("B_field_total_out","");
  int nx,ny,nz;
  double lx,ly,lz,tlon,tlat;

  if(btotfile.compare("") != 0 || bfield_type == 6 || dorand ) {	
    grid_interp=params.find<bool>("B_field_interp",true);
    nx=params.find<int>("B_field_nx",256);
    ny=params.find<int>("B_field_ny",nx);
    nz=params.find<int>("B_field_nz",nx/5);
    lx=params.find<double>("B_field_lx",40);
    ly=params.find<double>("B_field_ly",lx)*CGS_U_kpc;
    lz=params.find<double>("B_field_lz",(nz>1)?lx/5:0)*CGS_U_kpc;
    lx*=CGS_U_kpc;
    if (nz==1 && lz > 0) ErrLog("Don't give both nz==1 and Lz>0\n");
    tlon=params.find<double>("B_field_transform_lon",-999);
    tlat=params.find<double>("B_field_transform_lat",-999);
    
  } 

  if(btotfile.compare("") != 0) { 
    bool B_only=params.find<bool>("Only_write_B_field",false);
    setup_writer(btotfile,nx,ny,nz,lx,ly,lz,tlon,tlat,grid_interp,B_only);
  }


  if (bfield_type == 1) { 

    // All defaults;  nothing to set up
    double b0 =params.find<double>("B_field_b0",6.)*CGS_U_muGauss;
    double psi0=params.find<double>("B_field_psi0_deg",35)*CGS_U_pi/180.;
    double psi1=params.find<double>("B_field_psi1_deg",0.9)*CGS_U_pi/180.;
    double xsi0=params.find<double>("B_field_xsi0_deg",25)*CGS_U_pi/180.;
    setup_field1(b0,psi0,psi1,xsi0);
  }
  else if (bfield_type == 2) {
    double pitch_deg= params.find<double>("b2_field_pitch_deg",-10);
    double r0 =params.find<double>("b2_field_r0",10.55);
    double z0 =params.find<double>("b2_field_z0",1.);
    double b0 =params.find<double>("b2_field_b0",6.);
    b0*=CGS_U_muGauss;
    double phi0=params.find<double>("b2_field_phi0",CGS_U_pi);
    setup_field2(pitch_deg,r0,z0,b0,phi0);
  }
  else if(bfield_type == 3) {
    double pitch_deg=params.find<double>("b3_pitch_deg", -12);
    double B0=CGS_U_muGauss*params.find<double>("b3_B0", 2.);
    double Rsun=CGS_U_kpc*params.find<double>("b3_Rsun", 8.5);
    double R0=CGS_U_kpc*params.find<double>("b3_R0", 10.);
    double z0=CGS_U_kpc*params.find<double>("b3_z0", 1.);
    double Rc=CGS_U_kpc*params.find<double>("b3_Rc", 5.);
    double Bc=CGS_U_muGauss*params.find<double>("b3_Bc", 2.);
    // for the halo field
    double H_B0=CGS_U_muGauss*params.find<double>("b3H_B0", 10.);
    double H_z0=CGS_U_kpc*params.find<double>("b3H_z0", 1.5);
    double H_z1a=CGS_U_kpc*params.find<double>("b3H_z1a", 0.2);
    double H_z1b=CGS_U_kpc*params.find<double>("b3H_z1b", 0.4);
    double H_R0=CGS_U_kpc*params.find<double>("b3H_R0", 4.);

    setup_field3(pitch_deg, B0, Rsun, R0, z0, Rc, Bc, H_B0, H_z0, H_z1a, H_z1b, H_R0);
  }
  else if(bfield_type == 4) {
    double Rsun=CGS_U_kpc*params.find<double>("b4_Rsun", 8.5);
    double z1=CGS_U_kpc*params.find<double>("b4_z1",0.3);
    double z2=CGS_U_kpc*params.find<double>("b4_z2",4.);
    double r1=CGS_U_kpc*params.find<double>("b4_r1",2.);
    double p=(CGS_U_pi/180.)*params.find<double>("b4_p",-10.);
    double epsilon_0=CGS_U_kpc*params.find<double>("epsilon_0",10.55);

    setup_field4(Rsun, z1, z2, r1, p, epsilon_0);
  }
  else if(bfield_type == 5) {
    double b0=CGS_U_muGauss*params.find<double>("b5_b0", 1.4);
    double Rsun=CGS_U_kpc*params.find<double>("b5_Rsun",8.5);
    double r_min=CGS_U_kpc*params.find<double>("b5_r_min",4.);
    double d=CGS_U_kpc*params.find<double>("b5_d",-0.5);
    double z0=CGS_U_kpc*params.find<double>("b5_z0",1.5);
    double p=(CGS_U_pi/180.)*params.find<double>("b5_p",-8.);

    setup_field5(b0, Rsun, r_min, d, z0, p);
  }
	// field 6 - read a coherent field from a file
  else if (bfield_type == 6) {
    std::string inp_file=params.find<string>("B_field_coherent_inp","");
    setup_field6(lx, ly, lz, nx, ny, nz, tlon, tlat, grid_interp, inp_file);
  } 

  //Random field generation function, added by SPQ Aug 12 2014
  else if (bfield_type == 7) {
    long long unsigned int start_seed = params.find<long long unsigned int>("start_seed",800); //Seed to begin the realization
    double lmax = params.find<double>("lmax",100); // Outer scale of turbulence in pc
    unsigned int num_modes = params.find<unsigned int>("N_m",800); // Number of modes to sum for turbulent field.
    //Prepare pre-computed quantities
    double kmin = 2. * CGS_U_pi / lmax;
    double kmax = 1e5 * kmin;
    double Lc = lmax / 5.; // Correlation length as defined by Harari, Mollerach & Roulet, JHEP 03 045 2002
    double dkn_const = pow(10,( log10(kmax) - log10(kmin) ) / num_modes)-1.; //Log spacing constant
    double kn = kmin; // kn begins at kmin
    double nf = 0.; // Norm factor starts at 0
    double delta_kn = 0.; // Change in kn starts at 0
    // Compute power spectrum normalization constant
    for(unsigned int i=0;i < num_modes; i++){
        delta_kn = kn * dkn_const;
        nf = nf + 4. * CGS_U_pi * pow(kn,2) * delta_kn / (1. + pow(kn*Lc,11./3.) );
        kn = kn + delta_kn;
    }
    setup_field7(kmin,Lc,num_modes,dkn_const,nf,start_seed);
  }
  //End field 7 setup
 
  //Field 8: JF12_reg + JF12b_striated
  //This function loads a normalized random field file
  //Requires JF12a and JF12b style optimization parameters
  else if (bfield_type == 8) {
    double sqrt_beta = params.find<double>("sqrt_beta",1.16619);
    double kappa = params.find<double>("kappa",1.);
    string brand_file=params.find<string>("brand_filename","./nfield_out_1.txt");
    double b1       = params.find<double>("B_param_1",1.);
    double b2       = params.find<double>("B_param_2",1.);
    double b3       = params.find<double>("B_param_3",1.);
    double b4       = params.find<double>("B_param_4",1.);
    double b5       = params.find<double>("B_param_5",1.);
    double b6       = params.find<double>("B_param_6",1.);
    double b7       = params.find<double>("B_param_7",1.);
    double b8       = params.find<double>("B_param_8",1.);
    double b9       = params.find<double>("B_param_9",1.);
    double b10       = params.find<double>("B_param_10",1.);
    double b11       = params.find<double>("B_param_11",1.);
    double b12       = params.find<double>("B_param_12",1.);
    double b13       = params.find<double>("B_param_13",1.);
    double b14       = params.find<double>("B_param_14",1.);
    double b15       = params.find<double>("B_param_15",1.);
    double b16       = params.find<double>("B_param_16",1.);
    double b17       = params.find<double>("B_param_17",1.);
    double b18       = params.find<double>("B_param_18",1.);
    double b19       = params.find<double>("B_param_19",1.);
    double b20       = params.find<double>("B_param_20",1.);
    double br1       = params.find<double>("Br_param_1",1.);
    double br2       = params.find<double>("Br_param_2",1.);
    double br3       = params.find<double>("Br_param_3",1.);
    double br4       = params.find<double>("Br_param_4",1.);
    double br5       = params.find<double>("Br_param_5",1.);
    double br6       = params.find<double>("Br_param_6",1.);
    double br7       = params.find<double>("Br_param_7",1.);
    double br8       = params.find<double>("Br_param_8",1.);
    double br9       = params.find<double>("Br_param_9",1.);
    double br10       = params.find<double>("Br_param_10",1.);
    double br11       = params.find<double>("Br_param_11",1.);
    double br12       = params.find<double>("Br_param_12",1.);
    double br13       = params.find<double>("Br_param_13",1.);
    vector<double> brandx;
    vector<double> brandy;
    vector<double> brandz;
    string line;
    ifstream infile (brand_file.c_str());
    int file_length=0;
    if (infile.is_open())
    {
        cout << "Random field loaded successfully" << endl;
        while ( getline(infile,line) ){
            file_length = file_length + 1;
        }    
    }
    else{
        cout << "FAILED to load random field" << endl;
    }
    infile.close();
    ifstream infile2 (brand_file.c_str());
    for (int i=0;i<file_length;i++){
            double bx_temp = 0.;
            double by_temp = 0.;
            double bz_temp = 0.;
            infile2 >> bx_temp;
            infile2 >> by_temp;
            infile2 >> bz_temp;
            brandx.push_back(bx_temp);
            brandy.push_back(by_temp);
            brandz.push_back(bz_temp);
    }
    infile2.close();
    setup_field8(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
                 b15,b16,b17,b18,b19,b20,br1,br2,br3,br4,br5,br6,
                 br7,br8,br9,br10,br11,br12,br13,brandx,brandy,brandz,
                 kappa,sqrt_beta);
}
    //End field8 setup

  //Field 9: JF12_reg + JF12b_striated + JF12b_random
  //This function loads a normalized random field file
  //Requires JF12a and JF12b style optimization parameters
  else if (bfield_type == 9) {
    double sqrt_beta = params.find<double>("sqrt_beta",1.16619);
    double kappa = params.find<double>("kappa",1.);
    string brand_file=params.find<string>("brand_filename","./nfield_out_1.txt");
    double b1       = params.find<double>("B_param_1",1.);
    double b2       = params.find<double>("B_param_2",1.);
    double b3       = params.find<double>("B_param_3",1.);
    double b4       = params.find<double>("B_param_4",1.);
    double b5       = params.find<double>("B_param_5",1.);
    double b6       = params.find<double>("B_param_6",1.);
    double b7       = params.find<double>("B_param_7",1.);
    double b8       = params.find<double>("B_param_8",1.);
    double b9       = params.find<double>("B_param_9",1.);
    double b10       = params.find<double>("B_param_10",1.);
    double b11       = params.find<double>("B_param_11",1.);
    double b12       = params.find<double>("B_param_12",1.);
    double b13       = params.find<double>("B_param_13",1.);
    double b14       = params.find<double>("B_param_14",1.);
    double b15       = params.find<double>("B_param_15",1.);
    double b16       = params.find<double>("B_param_16",1.);
    double b17       = params.find<double>("B_param_17",1.);
    double b18       = params.find<double>("B_param_18",1.);
    double b19       = params.find<double>("B_param_19",1.);
    double b20       = params.find<double>("B_param_20",1.);
    double br1       = params.find<double>("Br_param_1",1.);
    double br2       = params.find<double>("Br_param_2",1.);
    double br3       = params.find<double>("Br_param_3",1.);
    double br4       = params.find<double>("Br_param_4",1.);
    double br5       = params.find<double>("Br_param_5",1.);
    double br6       = params.find<double>("Br_param_6",1.);
    double br7       = params.find<double>("Br_param_7",1.);
    double br8       = params.find<double>("Br_param_8",1.);
    double br9       = params.find<double>("Br_param_9",1.);
    double br10       = params.find<double>("Br_param_10",1.);
    double br11       = params.find<double>("Br_param_11",1.);
    double br12       = params.find<double>("Br_param_12",1.);
    double br13       = params.find<double>("Br_param_13",1.);
    vector<double> brandx;
    vector<double> brandy;
    vector<double> brandz;
    string line;
    ifstream infile (brand_file.c_str());
    int file_length=0;
    if (infile.is_open())
    {
        cout << "Random field loaded successfully" << endl;
        while ( getline(infile,line) ){
            file_length = file_length + 1;
        }    
    }
    else{
        cout << "FAILED to load random field" << endl;
    }
    infile.close();
    ifstream infile2 (brand_file.c_str());
    for (int i=0;i<file_length;i++){
            double bx_temp = 0.;
            double by_temp = 0.;
            double bz_temp = 0.;
            infile2 >> bx_temp;
            infile2 >> by_temp;
            infile2 >> bz_temp;
            brandx.push_back(bx_temp);
            brandy.push_back(by_temp);
            brandz.push_back(bz_temp);
    }
    infile2.close();
    setup_field9(b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13,b14,
                 b15,b16,b17,b18,b19,b20,br1,br2,br3,br4,br5,br6,
                 br7,br8,br9,br10,br11,br12,br13,brandx,brandy,brandz,
                 kappa,sqrt_beta);
}
    //End field9 setup

  //Field10: Pure GJQ random, no scaling by JF12
  //This function loads a normalized random field file
  //The random field file is the only required parameter
  else if (bfield_type == 10) {
    string brand_file=params.find<string>("brand_filename","./nfield_out_1.txt");
    vector<double> brandx;
    vector<double> brandy;
    vector<double> brandz;
    string line;
    ifstream infile (brand_file.c_str());
    int file_length=0;
    if (infile.is_open())
    {
        cout << "Random field loaded successfully" << endl;
        while ( getline(infile,line) ){
            file_length = file_length + 1;
        }    
    }
    else{
        cout << "FAILED to load random field" << endl;
    }
    infile.close();
    ifstream infile2 (brand_file.c_str());
    for (int i=0;i<file_length;i++){
            double bx_temp = 0.;
            double by_temp = 0.;
            double bz_temp = 0.;
            infile2 >> bx_temp;
            infile2 >> by_temp;
            infile2 >> bz_temp;
            brandx.push_back(bx_temp);
            brandy.push_back(by_temp);
            brandz.push_back(bz_temp);
    }
    infile2.close();
    setup_field10(brandx,brandy,brandz);
  } 
  // field 10

  if (dorand) {
    // From Han et al. 2004 ApJ
    double cutoff=CGS_U_kpc*params.find<double>("B_field_cutoff",5.);
    int seed=params.find<int>("B_field_seed",0);
    double rms=CGS_U_muGauss*params.find<double>("B_field_RMS_uG",6.);
    double c0=params.find<double>("B_field_c0_uG",1);

    //Note, -0.37-2 is Han et al. (2004 ApJ) for 3D.  -1 for 2D, -0.37 in 1D.
    double alpha=params.find<double>("B_field_alpha",(nz>1) ? -0.37-2. : -0.37-1.);


    std::string file=params.find<string>("B_field_random_out","");
    double rmax_ran=params.find<double>("B_field_rmax_ran",20.)*CGS_U_kpc;
    double mem_lim=params.find<double>("B_ran_mem_lim",1.);
    std::string inp_file=params.find<string>("B_field_random_inp","");

    setup_random(alpha,cutoff,seed,rms,c0,lx,ly,lz,nx,ny,nz,file, rmax_ran, tlon, tlat, grid_interp, mem_lim,bfield_debug, inp_file);

  } else dorand=false;


} //read_B_params
//---------------------------------------------------------------


//---------------------------------------------------------------
/*
   The purpose of the setup functions is so that the calling program
   does not have to use only the parameter file.  The B_field can be
   set up with the constructor, given the parameter file, or with
   empty brackets, in which case the setup must be called later and
   given all of the parameters explicitely.  This allows you to use
   different parameters in one program, for example.
*/
//---------------------------------------------------------------
void B_field::setup_field1(double b0,double psi0,double psi1,double xsi0) {
  Log("Generating WMAP model magnetic field\n");
  b1_b0=b0;
  b1_Psi_0=psi0;
  b1_Psi_1=psi1;
  b1_Xsi_0=xsi0;
}
// setup_field1 end
//---------------------------------------------------------------


void B_field::setup_field2(double pitch_deg, double r0, double z0, double b0, double phi0){
  Log("Generating Stanev BSS model magnetic field\n");

  b2_pitch_deg=pitch_deg;
  b2_r0=r0;
  b2_z0=z0;
  b2_b0=b0;
  b2_phi0=phi0;
}
// setup_field2 end
//---------------------------------------------------------------


void B_field::setup_field3(double pitch_deg, double B0, double Rsun, double R0, double z0, double Rc, double Bc,
                           double H_B0, double H_z0, double H_z1a, double H_z1b, double H_R0){
  Log("Generating Sun et al. A&A V.477 2008 ASS+RING model magnetic field\n");
  b3_pitch_deg=pitch_deg;
  b3_B0=B0;
  b3_Rsun=Rsun;
  b3_R0=R0;
  b3_z0=z0;
  b3_Rc=Rc;
  b3_Bc=Bc;

  b3H_B0=H_B0;
  b3H_z0=H_z0;
  b3H_z1a=H_z1a;
  b3H_z1b=H_z1b;
  b3H_R0=H_R0;
}
// setup_field3 end
//---------------------------------------------------------------

void B_field::setup_field4(double Rsun, double z1, double z2, double r1, double p, double epsilon_0) {
  Log("Generating Kachelriess et al. APh 2007 HMR model magnetic field\n");

  b4_Rsun=Rsun;
  b4_z1=z1;
  b4_z2=z2;
  b4_r1=r1;
  b4_p=p;
  b4_epsilon_0=epsilon_0;
}

// setup_field4 end
//---------------------------------------------------------------

void B_field::setup_field5(double b0, double Rsun, double r_min, double d, double z0, double p) {
  Log("Generating Kachelriess et al. APh 2007 TT magnetic field model \n");

  b5_b0=b0;
  b5_Rsun=Rsun;
  b5_r_min=r_min;
  b5_d=d;
  b5_z0=z0;
  b5_p=p;
}

// setup_field5 end
//---------------------------------------------------------------



//field 6 used to read in a coherent field from a file;
void B_field::setup_field6(double lx,double ly,double lz,int nx,int ny,int nz, double tlon, double tlat, bool interp, std::string inp_file) {
	Log("Reading regular field from specified file \n");

	grid_lx=lx;
	grid_ly=ly;
	grid_lz=lz;
	grid_nx=nx;
	grid_ny=ny;
	grid_nz=nz;
	grid_interp=interp;

	tlon=tlon*CGS_U_pi/180.;
	tlat=tlat*CGS_U_pi/180.;

	breg_inp_file=inp_file;
	
	Log("Allocating arrays for coherent component\n");
	
	//b_file_x=static_cast<double*>(malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
    //b_file_y=static_cast<double*>(malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
    //b_file_z=static_cast<double*>(malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
	int n = sizeof(double)*grid_nx*grid_ny*grid_nz;
	cout << " n = " << n << "\n";
	b_file_x=new double[n];
	b_file_y=new double[n];
	b_file_z=new double[n];
	
	//Read in the file ----------------------------------------
	if(breg_inp_file.compare("") != 0)
    {
		cout << " Reading in breg input file " << endl;
		std::ifstream reg_inp;
		double tmp;
		reg_inp.open(breg_inp_file.c_str());
		
		for (int i=0;i<grid_nx;i++) 
			for (int j=0;j<grid_ny;j++) 
				for (int l=0;l<grid_nz;l++) {
					// Removed the weird indexing from the random calculation code - 
					// nothing complex here (??)
					int index=Index3d(grid_nx,grid_ny,grid_nz,i,j,l);
					//cout<<"Index="<<index<<" i="<<i<<" j="<<j<<" l="<<l<<"\n";
					//
					if (reg_inp.eof()) {
						Log("B-field ERROR:  input file ended before expected:  i=");
						std::cout<<i<<", j="<<j<<", l="<<l<<"\n";
						std::exit(-1);
					}
					
					// Replace these lines 					
					//reg_inp >> b_file_x[index];
					//reg_inp >> b_file_y[index];
					//reg_inp >> b_file_z[index];
					
					reg_inp.read( reinterpret_cast<char*>(&b_file_x[index]), sizeof( double ) );
					reg_inp.read( reinterpret_cast<char*>(&b_file_y[index]), sizeof( double ) );
					reg_inp.read( reinterpret_cast<char*>(&b_file_z[index]), sizeof( double ) );
					
					
					// Fixing the units:
					b_file_x[index]*=CGS_U_muGauss;
					b_file_y[index]*=CGS_U_muGauss;
					b_file_z[index]*=CGS_U_muGauss;

				}
		
		if (!reg_inp.eof()) reg_inp >> tmp;
		if (!reg_inp.eof()) {
			Log("B-field ERROR:  input file not ended when expected\n");
			std::exit(-1);
		}
		
		reg_inp.close();
    } else {
		Log("Field type 6 requires in input file containing the coherent field. Define param:breg_inp_file\n");
		
	}
	//Done reading the file ----------------------------------------

}
// setup_field6 end
//---------------------------------------------------------------

//Setup field 7
void B_field::setup_field7(double kmin,double Lc,unsigned int num_modes,
                           double nf,double dkn_const,long long unsigned int start_seed) {
    b7_kmin = kmin;
    b7_Lc = Lc;
    b7_N_m = num_modes;
    b7_nf = nf;
    b7_dkn_const = dkn_const;
    b7_start_seed = start_seed;
}
// Setup field 7 end 
//---------------------------------------------------------------

//Setup field 8
void B_field::setup_field8(double b1, double b2, double b3, double b4, double b5, 
                            double b6, double b7, double b8, double b9, double b10, 
	                        double b11, double b12, double b13, double b14, double b15, 
                            double b16, double b17, double b18, double b19, double b20,
                            double br1, double br2, double br3, double br4, double br5, double br6,
                            double br7, double br8, double br9, double br10, double br11,
                            double br12, double br13, vector<double> brandx, 
                            vector<double> brandy, vector<double> brandz, double kappa, double sqrt_beta){
  Log("Generating      field");
  b8_b1=b1;
  b8_b2=b2;
  b8_b3=b3;
  b8_b4=b4;
  b8_b5=b5;
  b8_b6=b6;
  b8_b7=b7;
  b8_b8=b8;
  b8_b9=b9;
  b8_b10=b10;
  b8_b11=b11;
  b8_b12=b12;
  b8_b13=b13;
  b8_b14=b14;
  b8_b15=b15;
  b8_b16=b16;
  b8_b17=b17;
  b8_b18=b18;
  b8_b19=b19;
  b8_b20=b20;
  br8_b1=br1;
  br8_b2=br2;
  br8_b3=br3;
  br8_b4=br4;
  br8_b5=br5;
  br8_b6=br6;
  br8_b7=br7;
  br8_b8=br8;
  br8_b9=br9;
  br8_b10=br10;
  br8_b11=br11;
  br8_b12=br12;
  br8_b13=br13;
  br8_kappa=kappa;
  br8_sqrt_beta=sqrt_beta;
  b8_brandx=brandx;
  b8_brandy=brandy;
  b8_brandz=brandz;
} // setup_field8

//Setup field 9
void B_field::setup_field9(double b1, double b2, double b3, double b4, double b5, 
                            double b6, double b7, double b8, double b9, double b10, 
	                        double b11, double b12, double b13, double b14, double b15, 
                            double b16, double b17, double b18, double b19, double b20,
                            double br1, double br2, double br3, double br4, double br5, double br6,
                            double br7, double br8, double br9, double br10, double br11,
                            double br12, double br13, vector<double> brandx, 
                            vector<double> brandy, vector<double> brandz, double kappa, double sqrt_beta){
  Log("Generating      field");
  b9_b1=b1;
  b9_b2=b2;
  b9_b3=b3;
  b9_b4=b4;
  b9_b5=b5;
  b9_b6=b6;
  b9_b7=b7;
  b9_b8=b8;
  b9_b9=b9;
  b9_b10=b10;
  b9_b11=b11;
  b9_b12=b12;
  b9_b13=b13;
  b9_b14=b14;
  b9_b15=b15;
  b9_b16=b16;
  b9_b17=b17;
  b9_b18=b18;
  b9_b19=b19;
  b9_b20=b20;
  br9_b1=br1;
  br9_b2=br2;
  br9_b3=br3;
  br9_b4=br4;
  br9_b5=br5;
  br9_b6=br6;
  br9_b7=br7;
  br9_b8=br8;
  br9_b9=br9;
  br9_b10=br10;
  br9_b11=br11;
  br9_b12=br12;
  br9_b13=br13;
  br9_kappa=kappa;
  br9_sqrt_beta=sqrt_beta;
  b9_brandx=brandx;
  b9_brandy=brandy;
  b9_brandz=brandz;
} // setup_field9

//Setup field10
void B_field::setup_field10(vector<double> brandx, vector<double> brandy, vector<double> brandz){
  b10_brandx=brandx;
  b10_brandy=brandy;
  b10_brandz=brandz;
} //setup_field10



//Setup B-field writer
void B_field::setup_writer(std::string btotfile,int out_nx, int out_ny, int out_nz, double out_lx,double out_ly,double out_lz, double tlon, double tlat,  bool interp, bool B_only) {
	btot_file=btotfile;
	grid_lx=out_lx;
	grid_ly=out_ly;
	grid_lz=out_lz;
	grid_nx=out_nx;
	grid_ny=out_ny;
	grid_nz=out_nz;
	tlon=tlon*CGS_U_pi/180.;
	tlat=tlat*CGS_U_pi/180.;
	grid_interp=interp;
	Bonly=B_only;
}
// setup_writer end
//---------------------------------------------------------------



//---------------------------------------------------------------
void B_field::setup_random(double alpha,double cutoff,int seed, double rms, double c0, double lx,double ly, double lz, int nx, int ny, int nz, std::string file, double rmax_ran, double transform_lon, double transform_lat, bool interp, double mem_lim,bool debug, std::string inp_file) {
  dorand=true;
  bran_alpha=alpha;
  bran_cutoff_kpc=cutoff;
  if (seed == 0) {
    // set it off the clock, so guaranteed to be different even in
    // cosmomc where called with same paramter file each time and
    // can't remember last time's run
    seed=int(time(NULL));
    Log("DEBUG:  random seed is ");std::cout<<seed<<std::endl<<std::flush;
  }
  bran_seed=seed;
  bran_rng.seed(seed);
  bran_random_rms=rms;
  bran_random_c0=c0;
  grid_lx=lx;
  grid_ly=ly;
  grid_lz=lz;
  grid_nx=nx;
  grid_ny=ny;
  grid_nz=nz;
  bran_file=file;
  bran_rmax=rmax_ran;
  tlon=transform_lon*CGS_U_pi/180.;
  tlat=transform_lat*CGS_U_pi/180.;
  grid_interp=interp;
  bran_mem_limit=mem_lim;
  bfield_debug=debug;
  bran_inp_file=inp_file;
}
//---------------------------------------------------------------

//---------------------------------------------------------------
// WMAP model, no parameters
//---------------------------------------------------------------
vec3 B_field::field1(vec3 coords){

        //      double R=coords.Length();
        double r=sqrt(coords.x*coords.x + coords.y*coords.y);
        double PHI=atan2(coords.y,coords.x);

  //Limits to the field. "r ranges from 3kpc to 20kpc" (Page et al 2006)
  if (r>(20.0*CGS_U_kpc) || r<(3.0*CGS_U_kpc ))
    {
      // To avoid nan problems later
      r=1;
    }

  // Functions given in the paper, Psi(r), and Xsi(z).
  double Psi_r;
  double Psi_0, Psi_1;
  //  Psi_0=0.610865238;        //corresponds to 35deg
  //  Psi_1=0.015707963;        //0.9deg
  Psi_0=b1_Psi_0;
  Psi_1=b1_Psi_1;
  Psi_r=Psi_0+Psi_1*log(r/(8.0*CGS_U_kpc));


  double Xsi_z;
  double Xsi_0;
  //  Xsi_0=0.436332313;        //25deg
  Xsi_0=b1_Xsi_0;
  Xsi_z=Xsi_0*tanh(coords.z/(1.0*CGS_U_kpc));   //I trust tanh() without testing it.

  // B-field in cylindrical coordinates:
  double B_cyl[3]={b1_b0*cos(Psi_r)*cos(Xsi_z), b1_b0*sin(Psi_r)*cos(Xsi_z), b1_b0*sin(Xsi_z) };

  vec3 B_vec3(0,0,0);
  if (r<3.*CGS_U_kpc || r>20.*CGS_U_kpc) { B_vec3*=0;} else {
          double B_cart[3];
          Cyl2Cart(PHI,B_cyl,B_cart);
          B_vec3.x=B_cart[0];
          B_vec3.y=B_cart[1];
          B_vec3.z=B_cart[2];
  }

  return B_vec3;


} // field1
//---------------------------------------------------------------

//---------------------------------------------------------------
// Stanev BSS model (astro-ph/9607086)
//---------------------------------------------------------------
vec3 B_field::field2(vec3 coords) {

  double r=sqrt(coords.x*coords.x + coords.y*coords.y);
  if (r > gc_r_max) { return vec3(0,0,0);}
  double PHI=atan2(coords.y,coords.x);
  double z=coords.z;

  double PHIprime=b2_phi0-PHI;  // PHIprime running clock-wise from neg. x-axis

  double p_ang=b2_pitch_deg*CGS_U_pi/180.;
  double beta=1./tan(p_ang);

  double Rsun=8.5*CGS_U_kpc;
  double r0=b2_r0*CGS_U_kpc;
  double z0=b2_z0*CGS_U_kpc;
  if (abs(z)>0.5*CGS_U_kpc) {z0=4.*CGS_U_kpc;}

  double B_0=b2_b0*Rsun/4./CGS_U_kpc*CGS_U_muGauss;
  if (r>4.*CGS_U_kpc) {B_0=b2_b0*Rsun/r*CGS_U_muGauss;}

  double B_cyl[3]={B_0*cos(PHIprime-beta*log(r/r0))*sin(p_ang)*exp(-abs(z)/z0),
                                   -B_0*cos(PHIprime-beta*log(r/r0))*cos(p_ang)*exp(-abs(z)/z0),
                                   0};
  vec3 B_vec3(0,0,0);
  if (r<1.*CGS_U_kpc) { B_vec3*=0;} else {
          double B_cart[3];
          Cyl2Cart(PHI,B_cyl,B_cart);
          B_vec3.x=B_cart[0];
          B_vec3.y=B_cart[1];
          B_vec3.z=B_cart[2];
  }
  return B_vec3;

}// field2
//---------------------------------------------------------------

//---------------------------------------------------------------
// Sun et al. A&A V.477 2008 ASS+RING model
//---------------------------------------------------------------
vec3 B_field::field3(vec3 coords) {

  double r=sqrt(coords.x*coords.x + coords.y*coords.y);
  double PHI=atan2(coords.y,coords.x);
  double z=coords.z;

  // first we set D2
  // ------------------------------------------------------------
  double D2;
  if(r > 7.5*CGS_U_kpc)
    {
      D2 = 1.;
    }
  else if(r <= 7.5*CGS_U_kpc && r> 6.*CGS_U_kpc)
    {
      D2 = -1.;
    }
  else if(r<=6.*CGS_U_kpc && r>5.*CGS_U_kpc)
    {
      D2 = 1.;
    }
  else if(r<=5.*CGS_U_kpc)
    {
      D2=-1.;
    }
  else {cerr << " Error! r in field3 is:" << r/CGS_U_kpc << "kpc " << endl; exit(1);}
  // ------------------------------------------------------------

  // now we set D1
  // ------------------------------------------------------------
  double D1;
  if(r>b3_Rc)
    {
      D1 = b3_B0*exp(-((r-b3_Rsun)/b3_R0)-(abs(z)/b3_z0));
    }
  else if(r<=b3_Rc)
    {
      D1 = b3_Bc;
    }
  else {cerr << " Error! r in field3 is:" << r/CGS_U_kpc << "kpc " << endl; exit(1);}
  // ------------------------------------------------------------

  double p_ang=b3_pitch_deg*CGS_U_pi/180.;
  double B_cyl[3]={D1*D2*sin(p_ang),
                   -D1*D2*cos(p_ang),
                   0.};

  // Taking into account the halo field
  double halo_field;

  // for better overview
  double b3H_z1_actual;
  if(abs(z)<b3H_z0) {b3H_z1_actual=b3H_z1a;}
  else{b3H_z1_actual=b3H_z1b;}
  double hf_piece1 = (b3H_z1_actual*b3H_z1_actual)/(b3H_z1_actual*b3H_z1_actual+(abs(z)-b3H_z0)*(abs(z)-b3H_z0));
  double hf_piece2 = exp(-(r-b3H_R0)/(b3H_R0));

  halo_field = b3H_B0*hf_piece1*(r/b3H_R0)*hf_piece2;

  B_cyl[1]+=halo_field;

  vec3 B_vec3(0,0,0);

  double B_cart[3];
  Cyl2Cart(PHI,B_cyl,B_cart);
  B_vec3.x=B_cart[0];
  B_vec3.y=B_cart[1];
  B_vec3.z=B_cart[2];

  return B_vec3;

}// field3
//---------------------------------------------------------------


//---------------------------------------------------------------
// Kachelriess et al. APh 2007 HMR model
//---------------------------------------------------------------
vec3 B_field::field4(vec3 coords) {

  double r=sqrt(coords.x*coords.x + coords.y*coords.y);
  double PHI=atan2(coords.y,coords.x);
  double z=coords.z;

  double f_z=(1./(2.*cosh(z/b4_z1)))+(1./(2.*cosh(z/b4_z2)));

  if(r<0.0000000005*CGS_U_kpc) {r=0.5*CGS_U_kpc;}

  double b_r=(3.*b4_Rsun/r)*tanh(r/b4_r1)*tanh(r/b4_r1)*tanh(r/b4_r1)*CGS_U_muGauss;

  double B_r_PHI=b_r*cos(PHI-((1./tan(b4_p))*log(r/b4_epsilon_0)));

  // B-field in cylindrical coordinates:
  double B_cyl[3]={B_r_PHI*sin(b4_p)*f_z, B_r_PHI*cos(b4_p)*f_z , 0.};

  vec3 B_vec3(0,0,0);
  if (r>gc_r_max) { B_vec3*=0;} else {
    double B_cart[3];
    Cyl2Cart(PHI,B_cyl,B_cart);
    B_vec3.x=B_cart[0];
    B_vec3.y=B_cart[1];
    B_vec3.z=B_cart[2];
  }

  return B_vec3;
}// field4
//---------------------------------------------------------------

//---------------------------------------------------------------
// Kachelriess et al. APh 2007 TT model
//---------------------------------------------------------------
vec3 B_field::field5(vec3 coords) {

  double r=sqrt(coords.x*coords.x + coords.y*coords.y);
  double PHI=atan2(coords.y,coords.x);
  double z=coords.z;

  double beta=1./tan(b5_p);
  if(abs(b5_p)<1.e-16){beta=1.;}

  double phase=(beta*log(1.+b5_d/b5_Rsun))-CGS_U_pi/2.;
  //  double epsilon0=(b5_Rsun+b5_d)*exp(-(CGS_U_pi/2.)*tan(b5_p));

  double sign;
  if(z<0){sign=1.;}
  else{sign=-1.;}
  double f_z=sign*exp(-(abs(z)/b5_z0));

  // there is a factor 1/cos(phase) difference between
  // the original TT and Kachelriess.
  //double b_r=b5_b0*(b5_Rsun/(r));
  double b_r=b5_b0*(b5_Rsun/(r*cos(phase)));

  //if(r<b5_r_min) {b_r=b5_b0*(b5_Rsun/(b5_r_min));}
  if(r<b5_r_min) {b_r=b5_b0*(b5_Rsun/(b5_r_min*cos(phase)));}

  double B_r_PHI=b_r*cos(PHI-beta*log(r/b5_Rsun)+phase);
  //  if(r<1.e-26*CGS_U_kpc){B_r_PHI=b_r*cos(PHI+phase);}

  // B-field in cylindrical coordinates:
  double B_cyl[3]={B_r_PHI*sin(b5_p)*f_z, B_r_PHI*cos(b5_p)*f_z , 0.};

  vec3 B_vec3(0,0,0);
  if (r>gc_r_max) { B_vec3*=0;} else {
    double B_cart[3];
    Cyl2Cart(PHI,B_cyl,B_cart);
    B_vec3.x=B_cart[0];
    B_vec3.y=B_cart[1];
    B_vec3.z=B_cart[2];
  }

  return B_vec3;
}// field5
//---------------------------------------------------------------



//---------------------------------------------------------------
// Read a coherent field from a file and return
//---------------------------------------------------------------
vec3 B_field::field6(vec3 coords) {

	//Does the interpolation to the xyz grid
	double r=sqrt(coords.x*coords.x + coords.y*coords.y);
	//cout << " r= " << r << endl;

	//if (r > b6_rmax || (b6_nz > 1 && abs(coords.z) > 0.49*b6_lz_kpc)) { return vec3(0,0,0);}
	
	double r_sun=sqrt((coords.x-SunPosition.x)*(coords.x-SunPosition.x) + coords.y*coords.y);
	if (r_sun > ec_r_max) return vec3(0,0,0);
	
	double dx=grid_lx/grid_nx;
	double dy=grid_ly/grid_ny;
	double dz=grid_lz/grid_nz;
	
	// By default, coords in (x,y,z) where Sun-toward-Galactic Center is
	// x, with zero at GC, Sun to north galactic pole is z, ditto, etc.
	// B-random gridded the same way.  BUT if b6_tlon within
	// +-360, i.e. reasonable instead of default -999, then the grid is
	// meant to start at the Sun and extend only in one direction.  So
	// rotate the input coordinates into the system the grid uses.
	bool transform=false;
	
	if (tlon < 2.*CGS_U_pi && tlon > -2.*CGS_U_pi) {
	  transform=true;
	  vec3 coords2;
	  coords=coords-SunPosition;
	  Cart2LOS(tlat, tlon, coords, coords2);
	  
	  coords=coords2;
	}
	
	int xl,yl,zl;
	double tmp= transform ? coords.x/dx : (coords.x + grid_lx/2.)/dx;
	if (tmp < 0 && abs(tmp) > 1.e-10 ) xl=-1; 
	else xl= int(floor(tmp));  
	
	tmp=(coords.y + grid_ly/2.)/dy;
	if (tmp < 0 && abs(tmp) > 1.e-10) yl=-1; 
	else yl= int(floor(tmp));  
    
	tmp= (coords.z + grid_lz/2.)/dz;
	if (tmp < 0 && abs(tmp) > 1.e-10) zl=-1; 
	else zl= int(floor(tmp));
	
	if (xl < 0 || yl < 0 || zl < 0 ) { return vec3(0,0,0);}
	if (xl > grid_nx-1 || yl>grid_ny-1 || zl>grid_nz-1) { return vec3(0,0,0);}
	
	vec3 B_vec3;
	
	// At the borders, the linear interpolation doesn't work, so there,
	// just return the edge value itself.  Can't make much difference,
	// especially since the field is generally zeroed well before the
	// edge.
	
	// trilinear interpolation (see e.g. Wikipedia):
	//
	// Second check is not for out-of-bounds but for too close to edge for interp.
	if (grid_interp && xl+1<grid_nx && yl+1<grid_ny && zl+1<grid_nz) {
		int idx1,idx2;
		vec3 i1,i2,j1,j2,w1,w2;
		double xd= transform ? (coords.x - xl*grid_lx/grid_nx)/dx  : (coords.x +grid_lx/2 - xl*grid_lx/grid_nx)/dx;
		double yd=(coords.y +grid_ly/2 - yl*grid_ly/grid_ny)/dy;
		double zd=(coords.z +grid_lz/2 - zl*grid_lz/grid_nz)/dz;
		
		idx1=Index3d(grid_nx,grid_ny,grid_nz,xl,yl,zl);
		idx2=Index3d(grid_nx,grid_ny,grid_nz,xl,yl,zl+1);
		i1=vec3(b_file_x[idx1]*(1.-zd) + b_file_x[idx2]*zd, 
				b_file_y[idx1]*(1.-zd) + b_file_y[idx2]*zd, 
				b_file_z[idx1]*(1.-zd) + b_file_z[idx2]*zd );
		idx1=Index3d(grid_nx,grid_ny,grid_nz,xl,yl+1,zl);
		idx2=Index3d(grid_nx,grid_ny,grid_nz,xl,yl+1,zl+1);
		i2=vec3(b_file_x[idx1]*(1.-zd) + b_file_x[idx2]*zd,
				b_file_y[idx1]*(1.-zd) + b_file_y[idx2]*zd,
				b_file_z[idx1]*(1.-zd) + b_file_z[idx2]*zd);
		idx1=Index3d(grid_nx,grid_ny,grid_nz,xl+1,yl,zl);
		idx2=Index3d(grid_nx,grid_ny,grid_nz,xl+1,yl,zl+1);
		j1=vec3(b_file_x[idx1]*(1.-zd) + b_file_x[idx2]*zd,
				b_file_y[idx1]*(1.-zd) + b_file_y[idx2]*zd,
				b_file_z[idx1]*(1.-zd) + b_file_z[idx2]*zd);
		idx1=Index3d(grid_nx,grid_ny,grid_nz,xl+1,yl+1,zl);
		idx2=Index3d(grid_nx,grid_ny,grid_nz,xl+1,yl+1,zl+1);
		j2=vec3(b_file_x[idx1]*(1.-zd) + b_file_x[idx2]*zd,
				b_file_y[idx1]*(1.-zd) + b_file_y[idx2]*zd,
				b_file_z[idx1]*(1.-zd) + b_file_z[idx2]*zd);
		
		w1=i1*(1.-yd) + i2*yd;
		w2=j1*(1.-yd) + j2*yd;
		
		B_vec3= w1*(1.-xd) + w2*xd;
		
		if (B_vec3.Length() > 5000.*CGS_U_muGauss) {
			cout<<"WARNING:  got big B-field for xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl;
			cout<<"B_vec3.Length() = "<< B_vec3.Length()/CGS_U_muGauss << " uG"<<endl;}
		// end of interpolation
	} else {
		
		int idx=Index3d(grid_nx,grid_ny,grid_nz,xl,yl,zl);
		
		B_vec3= vec3(b_file_x[idx],b_file_y[idx],b_file_z[idx]);
		if (B_vec3.Length() > 5000.*CGS_U_muGauss) {
			cout<<"WARNING:  got big B-field for xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl;
			cout<<"B_vec3.Length() = "<< B_vec3.Length()/CGS_U_muGauss << " uG"<< B_vec3.Length() << endl;}
	} // no interpolation
	
	
	//  If transforming to a different longitude, need to rotate.  Note
	//  that this only rotates in x-y plane.  The idea is a field
	//  pointing toward you when at lon=0 (all x) is still pointing
	//  toward you when you move the object to lon=90 (all y).
	//  A sign flip I don't understand but empirically, this mostly works.
	if (transform) {
		vec3 tmp;
		LOS2Cart(tlat, tlon, B_vec3, tmp);
		B_vec3=tmp;
	}
	
	return B_vec3;
		
}// field6
//---------------------------------------------------------------

//field7 Random field generation, Sean Quinn Aug 12 2014
//First step: initialize array to store field and coord values
vector<double> brfx;
vector<double> brfy;
vector<double> brfz;
vector<double> brxc;
vector<double> bryc;
vector<double> brzc;
//Next step: initialize global random seed
long long unsigned int globRandSeed=0;
//Next step: create dummy initialization variable used once in field7 function
unsigned int init_indx = 0;
//Final step: create index for brandfield array
unsigned int brf_indx = 0;

vec3 B_field::field7(vec3 coords){

 //Setup vectors and other globals
 if (init_indx == 0){ 
    brfx.resize(brfa_len); //This line and the 5 following lines resize the vectors to appropriate size.
    brfy.resize(brfa_len);
    brfz.resize(brfa_len);
    brxc.resize(brfa_len);
    bryc.resize(brfa_len);
    brzc.resize(brfa_len);
    globRandSeed = b7_start_seed; // This passes the start seed from the parameter file, used for different realizations
    init_indx = init_indx + 1; //This ensures the conditional is true only at the beginning of the run
 }
 
 double Rmax = 20*CGS_U_kpc; // outer boundary of GMF
 // x,y,z is a Galactocentric cartesian system, with the Sun on the negative x-axis
 double r     = sqrt(coords.x*coords.x + coords.y*coords.y);

 // define boundaries outside of which B is zero. No grid points beyond Rmax are written to the file
 if (r > Rmax)     { return vec3(0,0,0);}

 boost::random::uniform_real_distribution<> unrDistA(0,2. * CGS_U_pi);
 boost::random::uniform_real_distribution<> unrDist(0,2. * CGS_U_pi);

 double zNPrime, kN, deltaKNConst;
 double kMin = b7_kmin;
 int PtsN = b7_N_m;
 deltaKNConst = b7_dkn_const;
 kN = kMin;

 double sum_prefactor;

 complex<double> BRandX(0.,0), BRandY(0.,0.), BRandZ(0.,0.);
 double thetaN, phiN, betaN, alphaN;

 // Two loops required to separate the origin seed from all other points
if (coords.x/CGS_U_kpc==-8.5 and coords.y==0 and coords.z==0){ //If at solar, use solar seed
    int solar_seed = 0; //Solar seeds 
    for( int nRan=0; nRan < PtsN; nRan++ ){
      boost::random::mt19937 gen(solar_seed);
      thetaN = unrDistA(gen);
      solar_seed = solar_seed + 1;
      boost::random::mt19937 hen(solar_seed);
      phiN = unrDist(hen);
      solar_seed = solar_seed + 1;
      boost::random::mt19937 fen(solar_seed);
      betaN = unrDist(fen);
      solar_seed = solar_seed + 1;
      boost::random::mt19937 wren(solar_seed);
      alphaN = unrDist(wren);
      solar_seed = solar_seed + 1;

      complex<double> x_hat(cos(alphaN) * cos(thetaN) * cos(phiN), -sin(alphaN)*sin(phiN));
      complex<double> y_hat(cos(alphaN) * cos(thetaN) * sin(phiN), sin(alphaN)*cos(phiN));
      double z_hat = -cos(alphaN) * sin(thetaN);
      zNPrime = (sin(thetaN)*cos(phiN) * coords.x + sin(thetaN)*sin(phiN)*coords.y + cos(thetaN)*coords.z)/CGS_U_pc;
      sum_prefactor = kN * sqrt(kN * deltaKNConst)/sqrt(1+pow(b7_Lc*kN,11./3.));
      complex<double> exp_term_arg(0,zNPrime * kN + betaN);
      BRandX += sum_prefactor * x_hat * exp(exp_term_arg);
      BRandY += sum_prefactor * y_hat * exp(exp_term_arg);
      BRandZ += sum_prefactor * z_hat * exp(exp_term_arg);

      kN += kN * deltaKNConst;
     }
}
 else{ //For all other coordinates use globRandSeed, which continues to increment
     for( int nRan=0; nRan < PtsN; nRan++ ){

      boost::random::mt19937 gen(globRandSeed);
      thetaN = unrDistA(gen);
      globRandSeed = globRandSeed + 1;
      boost::random::mt19937 hen(globRandSeed);
      phiN = unrDist(hen);
      globRandSeed = globRandSeed + 1;
      boost::random::mt19937 fen(globRandSeed);
      betaN = unrDist(fen);
      globRandSeed = globRandSeed + 1;
      boost::random::mt19937 wren(globRandSeed);
      alphaN = unrDist(wren);
      globRandSeed = globRandSeed + 1;

      complex<double> x_hat(cos(alphaN) * cos(thetaN) * cos(phiN), -sin(alphaN)*sin(phiN));
      complex<double> y_hat(cos(alphaN) * cos(thetaN) * sin(phiN), sin(alphaN)*cos(phiN));
      double z_hat = -cos(alphaN) * sin(thetaN);
      zNPrime = (sin(thetaN)*cos(phiN) * coords.x + sin(thetaN)*sin(phiN)*coords.y + cos(thetaN)*coords.z)/CGS_U_pc;
      sum_prefactor = kN * sqrt(kN * deltaKNConst)/sqrt(1+pow(b7_Lc*kN,11./3.));
      complex<double> exp_term_arg(0,zNPrime * kN + betaN);
      BRandX += sum_prefactor * x_hat * exp(exp_term_arg);
      BRandY += sum_prefactor * y_hat * exp(exp_term_arg);
      BRandZ += sum_prefactor * z_hat * exp(exp_term_arg);

      kN += kN * deltaKNConst;
     }
 }

 double BRandX_final = 2. * sqrt(CGS_U_pi) * real(BRandX) / b7_nf;
 double BRandY_final = 2. * sqrt(CGS_U_pi) * real(BRandY) / b7_nf;
 double BRandZ_final = 2. * sqrt(CGS_U_pi) * real(BRandZ) / b7_nf;

 vec3 B_vec3(0.,0.,0.);

 brfx[brf_indx]=BRandX_final;
 brfy[brf_indx]=BRandY_final;
 brfz[brf_indx]=BRandZ_final;
 brxc[brf_indx]=coords.x/CGS_U_kpc;
 bryc[brf_indx]=coords.y/CGS_U_kpc;
 brzc[brf_indx]=coords.z/CGS_U_kpc;

 brf_indx = brf_indx + 1;

 return B_vec3;
 
}
//field7

//Field 8: JF12a + JF12b*GJQs Added by Sean Quinn Aug 15 2014
unsigned int b8_mmm = 0; //Random field vector index
//---------------------------------------------------------------
vec3 B_field::field8(vec3 coords) {

  // define the input parameters

  // disk parameters
  double b_arm_1 = b8_b1;   // micro Gauss field strength in spiral arm region 1
  double b_arm_2 = b8_b2;
  double b_arm_3 = b8_b3;
  double b_arm_4 = b8_b4;
  double b_arm_5 = b8_b5;
  double b_arm_6 = b8_b6;
  double b_arm_7 = b8_b7;
  double b_ring  = b8_b8;   // field strength in molecular ring
  double h_disk  = b8_b9;   // height of transition between disk and toroidal halo
  double w_disk  = b8_b10;  // transition width between disk and toroidal halo

  //Random Field disk parameters, Sean Quinn 2/12/2014
  double bRandDisk1 = br8_b1;
  double bRandDisk2 = br8_b2;
  double bRandDisk3 = br8_b3;
  double bRandDisk4 = br8_b4;
  double bRandDisk5 = br8_b5;
  double bRandDisk6 = br8_b6;
  double bRandDisk7 = br8_b7;
  double bRandDisk8 = br8_b8;
  double bRandDiskInt = br8_b9;
  double z0RandDisk = br8_b10;
  double B0RandHalo = br8_b11;
  double r0RandHalo = br8_b12;
  double z0RandHalo = br8_b13;
  double beta_rand = br8_sqrt_beta;
  double brand_scale = br8_kappa;

  // toroidal halo parameters
  double Bn = b8_b11;   // field strength in the north
  double Bs = b8_b12;   //            ... in the south
  double rn = b8_b13;   // transition radius in the north
  double rs = b8_b14;   //            ... in the south
  double wh = b8_b15;   // transition width
  double z0 = b8_b16;   // scale height in z

  // X-field parameters
  double B0_X    = b8_b17;
  double Xtheta  = b8_b18;
  double rpc_X   = b8_b19; // called r_X^c in paper
  double r0_X    = b8_b20; // called r_X in paper

  // -------------------------------------------------------------------------------------------------
  // define fixed parameters

  double Rmax   = 20*CGS_U_kpc;  // outer boundary of GMF
  double rho_GC = 1.*CGS_U_kpc;  // interior boundary of GMF


  // fixed disk parameters
  double inc     = 11.5;            // inclination, in degrees
  double rmin    = 5.*CGS_U_kpc;    // outer boundary of the molecular ring region
  double rcent   = 3.*CGS_U_kpc;    // inner boundary of the molecular ring region (field is zero within this radius)
  double f[8]    = {0.130, 0.165, 0.094, 0.122, 0.13, 0.118, 0.084, 0.156}; // fractions of cicumference spanned by each spiral, sums to unity
  double rc_B[8] = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5};             // the radii where the spiral arm boundaries cross the negative x-axis
  // --------------------------------------------------------------------------------------------------
  // x,y,z is a Galactocentric cartesian system, with the Sun on the negative x-axis
  double r     = sqrt(coords.x*coords.x + coords.y*coords.y);
  double rho   = sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z);
  double PHI   = atan2(coords.y,coords.x);
  double z     = coords.z;

  // define boundaries outside of which B is zero
  if (r > Rmax)     { return vec3(0,0,0);}

  // ---------------------------------------------------------------------------------------------------
  // Disk component:  a divergenceless form of Brown et al. (2007). (8 spiral regions, 7 free parameters, the 8th set to conserve flux)

  // set B0 to 1 muG at r=5 kpc (so that the b_i parameters refer to B strength at 5 kpc, after which they go as ~1/r):
  double B0       = rmin/r*CGS_U_muGauss;
  // the logistic equation, to be multiplied to the toroidal halo field, and (1-zprofile) multiplied to the disk:
  double zprofile = 1./(1+exp(-2./w_disk*(abs(z)/CGS_U_kpc-h_disk)) );

  double B_cyl_disk[3]={0,0,0}; // the disk field in cylindrical coordinates
  double B_cyl_Randdisk[3]={0,0,0}; //Random field disk component
  
  if ( (r > rcent) )  // disk field zero elsewhere
  {
    if (r < rmin) {   // circular field in molecular ring
       B_cyl_disk[1] = B0*b_ring*(1-zprofile);
       // Compute the r < 5 kpc disk random field magnitude, Sean Quinn 2/14/2014
    }
    else
    {
       // use flux conservation to calculate the field strength in the 8th spiral arm
       double bv_B [8] = {b_arm_1, b_arm_2, b_arm_3, b_arm_4, b_arm_5, b_arm_6, b_arm_7, 0.};
       double b7 = 0.;
       for (int i=0; i<7; i++){
            b7 -= f[i]*bv_B[i]/f[7];
        } // last spiral strength is set by the others, to conserve flux
       bv_B[7] = b7;

       // iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to.
       double b_disk   = 0;
       double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi));
       if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+CGS_U_pi)); }
       if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }
       // Modified to also select the appropriate random field disk coefficient, Sean Quinn 2/12/2014, original line commented out below.
       // for (int i=7; i>=0; i--){ if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i];} }    // "region 8,7,6,..,2"
       for (int i=7; i>=0; i--){
            if (r_negx < rc_B[i]*CGS_U_kpc){
                b_disk = bv_B[i];
            }
        }

       // the disk field in cylindrical coordinates
       B_cyl_disk[0] =  b_disk*B0*sin(CGS_U_pi/180.*inc)*(1-zprofile);
       B_cyl_disk[1] =  b_disk*B0*cos(CGS_U_pi/180.*inc)*(1-zprofile);
     }
   }
   
 // ---------------------------------------------------------------------------------------------------
 // -------KRF Disk component
	if ( (r < rmin) )
		B_cyl_Randdisk[1] = bRandDiskInt * exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.) * CGS_U_muGauss;
	else{
		double bv_RandB [8] = {bRandDisk1, bRandDisk2, bRandDisk3, bRandDisk4, bRandDisk5, bRandDisk6, bRandDisk7, 0.}; //Random constant array
		double b7rand = 0.;
		for (int i=0; i<7; i++){
			b7rand -= f[i]*bv_RandB[i]/f[7];
		} // last spiral strength is set by the others, to conserve flux
		bv_RandB[7] = b7rand;
		// iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to.
		double b_diskRand = 0; //Random disk constant
		double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi));
		if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+CGS_U_pi)); }
		if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }
		// Modified to also select the appropriate random field disk coefficient, Sean Quinn 2/12/2014, original line commented out below.
		// for (int i=7; i>=0; i--){ if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i];} }    // "region 8,7,6,..,2"
		for (int i=7; i>=0; i--){
			if (r_negx < rc_B[i]*CGS_U_kpc){
				b_diskRand = bv_RandB[i];
			}
		}
		//Random disk field
		B_cyl_Randdisk[0] = b_diskRand * B0 * sin(CGS_U_pi/180.*inc)* exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.);
		B_cyl_Randdisk[1] = b_diskRand * B0 * cos(CGS_U_pi/180.*inc) * exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.);
	 }
 // ---------------------------------------------------------------------------------------------------

 // ---------------------------------------------------------------------------------------------------
 // --- Toroidal halo component
  double b1, rh;
  double B_h = 0;
  double B_hrand = 0;
  if ( z >= 0) {     // NORTH
    b1  = Bn*CGS_U_muGauss;
    rh  = rn;               // transition radius between inner-outer region, (units added later)
  }
  else if ( z < 0 ){  // SOUTH
    b1  = Bs*CGS_U_muGauss;
	rh  = rs;
  }
  B_h = b1*(1. - 1./(1.+exp(-2./wh*(r/CGS_U_kpc-rh))))*exp(-(abs(z))/(z0*CGS_U_kpc)); // vertical exponentialfall-off
  B_hrand = (B0RandHalo*CGS_U_muGauss) * exp(-r/(r0RandHalo*CGS_U_kpc)) * exp(-pow(z/(z0RandHalo*CGS_U_kpc),2)/2.); //spq 2/13/14
  double B_cyl_h[3]={ 0., B_h*zprofile, 0. };
  double B_cyl_hrand[3]={0., B_hrand, 0.}; //spq 2/13/14

 // ---------------------------------------------------------------------------------------------------
 // --- X-field component

 // apply units to input parameters
 B0_X   *= CGS_U_muGauss;
 r0_X   *= CGS_U_kpc;
 rpc_X  *= CGS_U_kpc;
 Xtheta *= CGS_U_pi/180.;

 double    rp_X   =  0.;  // the mid-plane radius for the field line that pass through r
 double    B_X    =  0.;
 double    r_sign =  1.;
 if (z<0){ r_sign = -1.;}

 // dividing line between region with constant elevation angle, and the interior:
 double rc_X = rpc_X + abs(z)/tan(Xtheta);

 if (r<rc_X){ // interior region, with varying elevation angle
   rp_X   = r*rpc_X/rc_X;
   B_X    = B0_X * pow(rpc_X/rc_X ,2.) * exp(-rp_X/r0_X);
   Xtheta = atan( abs(z)/ (r-rp_X) );   // modified elevation angle in interior region
   if (z==0.){Xtheta=CGS_U_pi/2.;}      // to avoid some NaN
 }
 else {       // exterior region with constant elevation angle
   rp_X = r - abs(z)/tan(Xtheta);
   B_X  = B0_X * rp_X/r * exp(-rp_X/r0_X);
 }

 // X-field in cylindrical coordinates
 double B_cyl_X[3]={ B_X*cos(Xtheta)*r_sign, 0.  ,  B_X*sin(Xtheta) };

 // ---------------------------------------------------------------------------------------------------
 // --- add disk + halo components together -------
    
 double B_cyl[3] = {0,0,0};
 B_cyl[0] = B_cyl_disk[0] +  B_cyl_h[0] + B_cyl_X[0];
 B_cyl[1] = B_cyl_disk[1] +  B_cyl_h[1] + B_cyl_X[1];
 B_cyl[2] = B_cyl_disk[2] +  B_cyl_h[2] + B_cyl_X[2];

 double B_cart[3];
 double Branddisk_cart[3];
 double Brandhalo_cart[3];
 vec3 B_vec3(0.,0.,0.);
 Cyl2Cart(PHI,B_cyl,B_cart);
 Cyl2Cart(PHI,B_cyl_Randdisk,Branddisk_cart);
 Cyl2Cart(PHI,B_cyl_hrand,Brandhalo_cart);
 B_vec3.x=B_cart[0];
 B_vec3.y=B_cart[1];
 B_vec3.z=B_cart[2];

 if (rho < rho_GC) { 
 	B_vec3.x=0.;
	B_vec3.y=0.;
	B_vec3.z=0.;
 }

 //Random field magnitude is added in quadrature
 double BRdisk_mag2 = pow(Branddisk_cart[0],2)+pow(Branddisk_cart[1],2)+pow(Branddisk_cart[2],2);
 double BRhalo_mag2 = pow(Brandhalo_cart[0],2)+pow(Brandhalo_cart[1],2)+pow(Brandhalo_cart[2],2);
 double BRandMag = sqrt(BRdisk_mag2+BRhalo_mag2) * brand_scale;

 double BRandX_final = BRandMag * b8_brandx[b8_mmm];
 double BRandY_final = BRandMag * b8_brandy[b8_mmm];
 double BRandZ_final = BRandMag * b8_brandz[b8_mmm];

 double BsrfX = 0.;
 double BsrfY = 0.;
 double BsrfZ = 0.;

 if (rho > rho_GC){
 //Compute prefactor for striated random field (dot products and vector magnitudes)
 double Breg_magsq = pow(B_vec3.x,2)+pow(B_vec3.y,2)+pow(B_vec3.z,2);
 double Bsrf_lmbda = beta_rand * (BRandX_final*B_vec3.x+BRandY_final*B_vec3.y+BRandZ_final*B_vec3.z) / Breg_magsq;
 BsrfX = Bsrf_lmbda * B_vec3.x;
 BsrfY = Bsrf_lmbda * B_vec3.y;
 BsrfZ = Bsrf_lmbda * B_vec3.z;
 }
 
 B_vec3.x = B_vec3.x + BsrfX;
 B_vec3.y = B_vec3.y + BsrfY;
 B_vec3.z = B_vec3.z + BsrfZ;

 b8_mmm = b8_mmm + 1;

 return B_vec3;

} //field8


//Field 9: JF12a + JF12b*GJQs + JF12b*GJQ Added by Sean Quinn Aug 15 2014
unsigned int b9_mmm = 0; //Random field vector index
//---------------------------------------------------------------
vec3 B_field::field9(vec3 coords) {

  // define the input parameters

  // disk parameters
  double b_arm_1 = b9_b1;   // micro Gauss field strength in spiral arm region 1
  double b_arm_2 = b9_b2;
  double b_arm_3 = b9_b3;
  double b_arm_4 = b9_b4;
  double b_arm_5 = b9_b5;
  double b_arm_6 = b9_b6;
  double b_arm_7 = b9_b7;
  double b_ring  = b9_b8;   // field strength in molecular ring
  double h_disk  = b9_b9;   // height of transition between disk and toroidal halo
  double w_disk  = b9_b10;  // transition width between disk and toroidal halo

  //Random Field disk parameters, Sean Quinn 2/12/2014
  double bRandDisk1 = br9_b1;  //variable names are suggestive of JF12b paper
  double bRandDisk2 = br9_b2;
  double bRandDisk3 = br9_b3;
  double bRandDisk4 = br9_b4;
  double bRandDisk5 = br9_b5;
  double bRandDisk6 = br9_b6;
  double bRandDisk7 = br9_b7;
  double bRandDisk8 = br9_b8;
  double bRandDiskInt = br9_b9;
  double z0RandDisk = br9_b10;
  double B0RandHalo = br9_b11;
  double r0RandHalo = br9_b12;
  double z0RandHalo = br9_b13;
  double beta_rand = br9_sqrt_beta;
  double brand_scale = br9_kappa;

  // toroidal halo parameters
  double Bn = b9_b11;   // field strength in the north
  double Bs = b9_b12;   //            ... in the south
  double rn = b9_b13;   // transition radius in the north
  double rs = b9_b14;   //            ... in the south
  double wh = b9_b15;   // transition width
  double z0 = b9_b16;   // scale height in z

  // X-field parameters
  double B0_X    = b9_b17;
  double Xtheta  = b9_b18;
  double rpc_X   = b9_b19; // called r_X^c in paper
  double r0_X    = b9_b20; // called r_X in paper

  // -------------------------------------------------------------------------------------------------
  // define fixed parameters

  double Rmax   = 20*CGS_U_kpc;  // outer boundary of GMF
  double rho_GC = 1.*CGS_U_kpc;  // interior boundary of GMF


  // fixed disk parameters
  double inc     = 11.5;            // inclination, in degrees
  double rmin    = 5.*CGS_U_kpc;    // outer boundary of the molecular ring region
  double rcent   = 3.*CGS_U_kpc;    // inner boundary of the molecular ring region (field is zero within this radius)
  double f[8]    = {0.130, 0.165, 0.094, 0.122, 0.13, 0.118, 0.084, 0.156}; // fractions of cicumference spanned by each spiral, sums to unity
  double rc_B[8] = {5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5};             // the radii where the spiral arm boundaries cross the negative x-axis
  // --------------------------------------------------------------------------------------------------
  // x,y,z is a Galactocentric cartesian system, with the Sun on the negative x-axis
  double r     = sqrt(coords.x*coords.x + coords.y*coords.y);
  double rho   = sqrt(coords.x*coords.x + coords.y*coords.y + coords.z*coords.z);
  double PHI   = atan2(coords.y,coords.x);
  double z     = coords.z;

  // define boundaries outside of which B is zero
  if (r > Rmax)     { return vec3(0,0,0);}

  // ---------------------------------------------------------------------------------------------------
  // Disk component:  a divergenceless form of Brown et al. (2007). (8 spiral regions, 7 free parameters, the 8th set to conserve flux)

  // set B0 to 1 muG at r=5 kpc (so that the b_i parameters refer to B strength at 5 kpc, after which they go as ~1/r):
  double B0       = rmin/r*CGS_U_muGauss;
  // the logistic equation, to be multiplied to the toroidal halo field, and (1-zprofile) multiplied to the disk:
  double zprofile = 1./(1+exp(-2./w_disk*(abs(z)/CGS_U_kpc-h_disk)) );

  double B_cyl_disk[3]={0,0,0}; // the disk field in cylindrical coordinates
  double B_cyl_Randdisk[3]={0,0,0}; //Random field disk component
  
  if ( (r > rcent) )  // disk field zero elsewhere
  {
    if (r < rmin) {   // circular field in molecular ring
       B_cyl_disk[1] = B0*b_ring*(1-zprofile);
       // Compute the r < 5 kpc disk random field magnitude, Sean Quinn 2/14/2014
    }
    else
    {
       // use flux conservation to calculate the field strength in the 8th spiral arm
       double bv_B [8] = {b_arm_1, b_arm_2, b_arm_3, b_arm_4, b_arm_5, b_arm_6, b_arm_7, 0.};
       double b7 = 0.;
       for (int i=0; i<7; i++){
            b7 -= f[i]*bv_B[i]/f[7];
        } // last spiral strength is set by the others, to conserve flux
       bv_B[7] = b7;

       // iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to.
       double b_disk   = 0;
       double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi));
       if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+CGS_U_pi)); }
       if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }
       // Modified to also select the appropriate random field disk coefficient, Sean Quinn 2/12/2014, original line commented out below.
       // for (int i=7; i>=0; i--){ if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i];} }    // "region 8,7,6,..,2"
       for (int i=7; i>=0; i--){
            if (r_negx < rc_B[i]*CGS_U_kpc){
                b_disk = bv_B[i];
            }
        }

       // the disk field in cylindrical coordinates
       B_cyl_disk[0] =  b_disk*B0*sin(CGS_U_pi/180.*inc)*(1-zprofile);
       B_cyl_disk[1] =  b_disk*B0*cos(CGS_U_pi/180.*inc)*(1-zprofile);
     }
   }
   
 // ---------------------------------------------------------------------------------------------------
 // -------KRF Disk component
	if ( (r < rmin) )
		B_cyl_Randdisk[1] = bRandDiskInt * exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.) * CGS_U_muGauss;
	else{
		double bv_RandB [8] = {bRandDisk1, bRandDisk2, bRandDisk3, bRandDisk4, bRandDisk5, bRandDisk6, bRandDisk7, 0.}; //Random constant array
		double b7rand = 0.;
		for (int i=0; i<7; i++){
			b7rand -= f[i]*bv_RandB[i]/f[7];
		} // last spiral strength is set by the others, to conserve flux
		bv_RandB[7] = b7rand;
		// iteratively figure out which spiral arm the current coordinates (r, phi) corresponds to.
		double b_diskRand = 0; //Random disk constant
		double r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI-CGS_U_pi));
		if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+CGS_U_pi)); }
		if (r_negx > rc_B[7]*CGS_U_kpc) {r_negx = r*exp(-1/tan(CGS_U_pi/180.*(90-inc))*(PHI+3*CGS_U_pi)); }
		// Modified to also select the appropriate random field disk coefficient, Sean Quinn 2/12/2014, original line commented out below.
		// for (int i=7; i>=0; i--){ if (r_negx < rc_B[i]*CGS_U_kpc) { b_disk = bv_B[i];} }    // "region 8,7,6,..,2"
		for (int i=7; i>=0; i--){
			if (r_negx < rc_B[i]*CGS_U_kpc){
				b_diskRand = bv_RandB[i];
			}
		}
		//Random disk field
		B_cyl_Randdisk[0] = b_diskRand * B0 * sin(CGS_U_pi/180.*inc)* exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.);
		B_cyl_Randdisk[1] = b_diskRand * B0 * cos(CGS_U_pi/180.*inc) * exp(-pow(z/(z0RandDisk*CGS_U_kpc),2)/2.);
	 }
 // ---------------------------------------------------------------------------------------------------

 // ---------------------------------------------------------------------------------------------------
 // --- Toroidal halo component
  double b1, rh;
  double B_h = 0;
  double B_hrand = 0;
  if ( z >= 0) {     // NORTH
    b1  = Bn*CGS_U_muGauss;
    rh  = rn;               // transition radius between inner-outer region, (units added later)
  }
  else if ( z < 0 ){  // SOUTH
    b1  = Bs*CGS_U_muGauss;
	rh  = rs;
  }
  B_h = b1*(1. - 1./(1.+exp(-2./wh*(r/CGS_U_kpc-rh))))*exp(-(abs(z))/(z0*CGS_U_kpc)); // vertical exponentialfall-off
  B_hrand = (B0RandHalo*CGS_U_muGauss) * exp(-r/(r0RandHalo*CGS_U_kpc)) * exp(-pow(z/(z0RandHalo*CGS_U_kpc),2)/2.); //spq 2/13/14
  double B_cyl_h[3]={ 0., B_h*zprofile, 0. };
  double B_cyl_hrand[3]={0., B_hrand, 0.}; //spq 2/13/14

 // ---------------------------------------------------------------------------------------------------
 // --- X-field component

 // apply units to input parameters
 B0_X   *= CGS_U_muGauss;
 r0_X   *= CGS_U_kpc;
 rpc_X  *= CGS_U_kpc;
 Xtheta *= CGS_U_pi/180.;

 double    rp_X   =  0.;  // the mid-plane radius for the field line that pass through r
 double    B_X    =  0.;
 double    r_sign =  1.;
 if (z<0){ r_sign = -1.;}

 // dividing line between region with constant elevation angle, and the interior:
 double rc_X = rpc_X + abs(z)/tan(Xtheta);

 if (r<rc_X){ // interior region, with varying elevation angle
   rp_X   = r*rpc_X/rc_X;
   B_X    = B0_X * pow(rpc_X/rc_X ,2.) * exp(-rp_X/r0_X);
   Xtheta = atan( abs(z)/ (r-rp_X) );   // modified elevation angle in interior region
   if (z==0.){Xtheta=CGS_U_pi/2.;}      // to avoid some NaN
 }
 else {       // exterior region with constant elevation angle
   rp_X = r - abs(z)/tan(Xtheta);
   B_X  = B0_X * rp_X/r * exp(-rp_X/r0_X);
 }

 // X-field in cylindrical coordinates
 double B_cyl_X[3]={ B_X*cos(Xtheta)*r_sign, 0.  ,  B_X*sin(Xtheta) };

 // ---------------------------------------------------------------------------------------------------
 // --- add disk + halo components together -------
    
 double B_cyl[3] = {0,0,0};
 B_cyl[0] = B_cyl_disk[0] +  B_cyl_h[0] + B_cyl_X[0];
 B_cyl[1] = B_cyl_disk[1] +  B_cyl_h[1] + B_cyl_X[1];
 B_cyl[2] = B_cyl_disk[2] +  B_cyl_h[2] + B_cyl_X[2];

 double B_cart[3];
 double Branddisk_cart[3];
 double Brandhalo_cart[3];
 vec3 B_vec3(0.,0.,0.);
 Cyl2Cart(PHI,B_cyl,B_cart);
 Cyl2Cart(PHI,B_cyl_Randdisk,Branddisk_cart);
 Cyl2Cart(PHI,B_cyl_hrand,Brandhalo_cart);
 B_vec3.x=B_cart[0];
 B_vec3.y=B_cart[1];
 B_vec3.z=B_cart[2];

 if (rho < rho_GC) { 
 	B_vec3.x=0.;
	B_vec3.y=0.;
	B_vec3.z=0.;
 }

 //Random field magnitude is added in quadrature
 double BRdisk_mag2 = pow(Branddisk_cart[0],2)+pow(Branddisk_cart[1],2)+pow(Branddisk_cart[2],2);
 double BRhalo_mag2 = pow(Brandhalo_cart[0],2)+pow(Brandhalo_cart[1],2)+pow(Brandhalo_cart[2],2);
 double BRandMag = sqrt(BRdisk_mag2+BRhalo_mag2) * brand_scale;

 double BRandX_final = BRandMag * b9_brandx[b9_mmm];
 double BRandY_final = BRandMag * b9_brandy[b9_mmm];
 double BRandZ_final = BRandMag * b9_brandz[b9_mmm];

 double BsrfX = 0.;
 double BsrfY = 0.;
 double BsrfZ = 0.;

 if (rho > rho_GC){
 //Compute prefactor for striated random field (dot products and vector magnitudes)
 double Breg_magsq = pow(B_vec3.x,2)+pow(B_vec3.y,2)+pow(B_vec3.z,2);
 double Bsrf_lmbda = beta_rand * (BRandX_final*B_vec3.x+BRandY_final*B_vec3.y+BRandZ_final*B_vec3.z) / Breg_magsq;
 BsrfX = Bsrf_lmbda * B_vec3.x;
 BsrfY = Bsrf_lmbda * B_vec3.y;
 BsrfZ = Bsrf_lmbda * B_vec3.z;
 }
 
 B_vec3.x = B_vec3.x + BsrfX + BRandX_final;
 B_vec3.y = B_vec3.y + BsrfY + BRandY_final;
 B_vec3.z = B_vec3.z + BsrfZ + BRandZ_final;

 b9_mmm = b9_mmm + 1;

 return B_vec3;

} //field9


//Field 10: Pure GJQ (no scaling) Added by Sean Quinn Aug 15 2014
unsigned int b10_mmm = 0; //Random field vector index
//---------------------------------------------------------------
vec3 B_field::field10(vec3 coords) {
  vec3 B_vec3(0.,0.,0.);
  B_vec3.x = b10_brandx[b10_mmm]*CGS_U_muGauss;
  B_vec3.y = b10_brandx[b10_mmm]*CGS_U_muGauss;
  B_vec3.z = b10_brandx[b10_mmm]*CGS_U_muGauss;
  b10_mmm = b10_mmm + 1;
  return B_vec3;
} //field10


void B_field::fillRandom(bool do_alloc, std::string infile){
        bran_inp_file=infile;
        fillRandom(do_alloc);
}

//---------------------------------------------------------------
/*
   To be called in the constructor for the B_field if the parameter
   B_field_do_random is true.  Otherwise, use setup_random and then
   call fillRandom before needed.  (Call cleanRandom() before
   another realisation.)
 */
void B_field::fillRandom( bool do_alloc ){
  if (!dorand) return;

  double k0=grid_nx/(3.*grid_lx) ; 
  double k_cutoff=1./bran_cutoff_kpc;

  if (do_alloc) {
    Log("Allocating arrays for random component\n");
    b_of_k_x=static_cast<double*>(fftw_malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
    b_of_k_y=static_cast<double*>(fftw_malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
    b_of_k_z=static_cast<double*>(fftw_malloc( sizeof(double)*grid_nx*grid_ny*2*(grid_nz/2+1)));
  }

  // Should the condition below be satisfied we load a file
  // with the random field information in real space.
  // After that the routine finishes.
  if(bran_inp_file.compare("") != 0)
    {
      cout << " Reading in bran input file " << endl;
      std::ifstream ran_inp;
      double tmp;
      ran_inp.open(bran_inp_file.c_str());

      for (int i=0;i<grid_nx;i++)
	for (int j=0;j<grid_ny;j++)
	  for (int l=0;l<grid_nz;l++) {
	    // Note the weird indexing (2*(grid_nz/2+1).
	    // This is due to the fact that this array is
	    // also used for complex computations.
	    int index=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);
	    if (ran_inp.eof()) {
	      Log("B-field ERROR:  input file ended before expected:  i="); 
	      std::cout<<i<<", j="<<j<<", l="<<l<<"\n"; 
	      std::exit(-1);
	    }
	    ran_inp.read( reinterpret_cast<char*>(&b_of_k_x[index]), sizeof( double ) );
	    ran_inp.read( reinterpret_cast<char*>(&b_of_k_y[index]), sizeof( double ) );
	    ran_inp.read( reinterpret_cast<char*>(&b_of_k_z[index]), sizeof( double ) );
	    
	    // Fixing the units:
	    b_of_k_x[index]*=CGS_U_muGauss;
	    b_of_k_y[index]*=CGS_U_muGauss;
	    b_of_k_z[index]*=CGS_U_muGauss;
	  }

      if (!ran_inp.eof()) ran_inp >> tmp;
      if (!ran_inp.eof()) {
	Log("B-field ERROR:  input file not ended when expected\n");
	std::exit(-1);
      }
      ran_inp.close();
    } else {

  fftw_complex *dummy_x=reinterpret_cast<fftw_complex*>(b_of_k_x);
  fftw_complex *dummy_y=reinterpret_cast<fftw_complex*>(b_of_k_y);
  fftw_complex *dummy_z=reinterpret_cast<fftw_complex*>(b_of_k_z);

  Log("Making FFTW3 plans\n");
  fftw_plan px = fftw_plan_dft_c2r_3d( grid_nx, grid_ny, grid_nz, dummy_x, b_of_k_x, FFTW_MEASURE);
  fftw_plan py = fftw_plan_dft_c2r_3d( grid_nx, grid_ny, grid_nz, dummy_y, b_of_k_y, FFTW_MEASURE);
  fftw_plan pz = fftw_plan_dft_c2r_3d( grid_nx, grid_ny, grid_nz, dummy_z, b_of_k_z, FFTW_MEASURE);

  // Burn in RNG?
  for (int i=0;i<10000;i++) { bran_rng.rand_gauss();}

  //  Log("brandom:  Filling Fourier space arrays with power-law random component\n",1);

  //  arr<float> flags(n*n*2*(grid_nz/2+1));  flags.fill(0);

  Log("Filling grid in Fourier space\n");

  for (int i=0;i<grid_nx;i++) {
    for (int j=0;j<grid_ny;j++) {
      for (int l=0;l<2*(grid_nz/2+1);l+=2) {
	int idx=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);

	// FFT expects up to n/2 positive while n/2 to n negative
	double kx= (i < grid_nx/2) ? i : i-grid_nx;
        double ky= (j < grid_ny/2) ? j : j-grid_ny;

        // l stepping through z-direction of array, which spans
        // 2*(n/2+1) doubles, representing n/2+1 complex numbers, so the
        // corresponding wavenumber k_z is actually l/2 (l always
        // even)
        double kz= l/2;

	// But then we want to normalize these by the dimensions:
	kx*=(1./grid_lx); 
	ky*=(1./grid_ly); 
	kz*=(1./grid_lz); 

        double k = sqrt( kx*kx + ky*ky + kz*kz );

        vec3 kvec(kx,ky,kz);

        if ( (bran_cutoff_kpc > 1.e-10) && (k < k_cutoff)) {
          b_of_k_x[idx] = 0.;
          b_of_k_x[idx+1] = 0.;
          b_of_k_y[idx] = 0.;
          b_of_k_y[idx+1] = 0.;
          b_of_k_z[idx] = 0.;
          b_of_k_z[idx+1] = 0.;
          // So that two sims with different cutoffs end up with the
          // same random numbers in the same places.
          for (int q=0;q<6;q++) bran_rng.rand_gauss();
        } else {


          // See Phys Rev E 66, 038301 (2002) for where I got this trick:

          // Generate a random complex vector whose length follows
          //  the required power law distribution and whose direction
          //  is random:

          /* P(k)=<B(k)^2>=<B_x(k)^2+B_y(k)^2+B_z(k)^2>propto k^alpha
             therefore real(B_x) propto sqrt(P(k)/6) so that
             B_x(k)^2=real(B_x)^2+imag(B_x)^2 terms, of which there
             will be six total, add up to P(k).  But norm irrelevant,
             set later.  Sqrt means alpha/2.
          */

          // To make debugging easier:
          vec3 rans=vec3( bran_rng.rand_gauss(), bran_rng.rand_gauss(), bran_rng.rand_gauss());
          //vec3 rans=vec3(1.,1.,1.);

	  // Real part:
	  vec3 a(
		 rans.x*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)),
		 rans.y*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)),
		 rans.z*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)) );
	  
	  //  Note: no 1/sqrt(grid_nx*grid_ny*grid_nz) factor.  True,
	  //  you would need it to get correct results for something
	  //  put through both forward and backward transforms, but
	  //  you're not doing that here.  Indeed, here you want to be
	  //  able to give the same c0 at two resolutions and get the
	  //  same amount of large-scale structure.  

          // TO BE FIXED: this is messing up the isotropy for some
          // reason.  Define b as this function of a and k.  The
          // result is that the divergence of the whole brings down a
          // k from the exponential function, which when dotted by
          // this gives zero, i.e. zero divergence.
          vec3 b= a;// - kvec * ( dotprod(a,kvec)/kvec.SquaredLength());
          b_of_k_x[idx]=b.x;
          b_of_k_y[idx]=b.y;
          b_of_k_z[idx]=b.z;


	  // Repeat for imaginary part:
	  rans=vec3( bran_rng.rand_gauss(), bran_rng.rand_gauss(), bran_rng.rand_gauss());
	  //rans=vec3(1.,1.,1.);
	  a=vec3(
		 rans.x*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)), ///sqrt(grid_nx*grid_ny*grid_nz),
		 rans.y*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)), ///sqrt(grid_nx*grid_ny*grid_nz),
		 rans.z*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)) );///sqrt(grid_nx*grid_ny*grid_nz));


          // Repeat for imaginary part:
          rans=vec3( bran_rng.rand_gauss(), bran_rng.rand_gauss(), bran_rng.rand_gauss());
          //rans=vec3(1.,1.,1.);
          a=vec3(
                 rans.x*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)), ///sqrt(bran_nx*bran_ny*bran_nz),
                 rans.y*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)), ///sqrt(bran_nx*bran_ny*bran_nz),
                 rans.z*sqrt(bran_random_c0)*pow(k/k_cutoff,bran_alpha/2)*exp(-0.5*k*k/(k0*k0)) );///sqrt(bran_nx*bran_ny*bran_nz));



          b= a;// - kvec * ( dotprod(a,kvec)/kvec.SquaredLength());
          b_of_k_x[idx+1]=b.x;
          b_of_k_y[idx+1]=b.y;
          b_of_k_z[idx+1]=b.z;



        } // not bigger than cut-off


      } // l
    } // j
  } // i


  if (bfield_debug){
    std::ofstream ran_out_k;
    ran_out_k.open("debug_bran_x_k.bin",std::ios::out|std::ios::binary);

    for (int i=0;i<grid_nx;i++)
      for (int j=0;j<grid_ny;j++)
      for (int l=0;l<2*(grid_nz/2+1);l+=2) {
	int index=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);
	float tmp=float(b_of_k_x[index]);
	ran_out_k.write( reinterpret_cast<char*>(&tmp),sizeof(float)); 
	tmp=float(b_of_k_x[index+1]);
	ran_out_k.write( reinterpret_cast<char*>(&tmp),sizeof(float)); 
      }

    ran_out_k.close();
  }


  Log("Running FFTs\n");
  //----------------
  // Now, do FFTs.
  //----------------
  fftw_execute(px);
  fftw_execute(py);
  fftw_execute(pz);
  fftw_destroy_plan(px);
  fftw_destroy_plan(py);
  fftw_destroy_plan(pz);
  } // else (when not given input file)

  bool error=false; // defunct
  if (bran_file.compare("") != 0 || error){
    if (bran_file.compare("")==0) bran_file="debug_error_bran.bin";
    Log(std::string("Writing random component to ")+bran_file+"\n");
    std::ofstream ran_out;
        // Append N.bin, and increment:
        int i=0,done=0;
        do {
                std::stringstream filename;
                filename<<bran_file<<i<<".bin";
                std::ifstream in(filename.str().data());
                if (in.is_open()) { i++; in.close(); continue;}
                else {
                        ran_out.open(filename.str().data(),std::ios::out|std::ios::binary);
                        for (int i=0;i<grid_nx;i++)
                        for (int j=0;j<grid_ny;j++)
                        for (int l=0;l<grid_nz;l++) {
                        int index=(grid_nz==1) ? Index3d(grid_nx,2*(grid_ny/2+1),1,i,j,l) : Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);
                        //=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);
                        // b_of_k_x is in muGauss, which isn't quite what you think it
                        // is.  Write it out by dividing its units and putting back
                        // the 10^-6.
                        //                      float tmp=float(b_of_k_x[index]*1.e-6/CGS_U_muGauss);
                        double tmp=b_of_k_x[index]*1.e-6/CGS_U_muGauss;
                        ran_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
                        tmp=b_of_k_y[index]*1.e-6/CGS_U_muGauss;
                        ran_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
                        tmp=b_of_k_z[index]*1.e-6/CGS_U_muGauss;
                        ran_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
                        }

                        ran_out.close();
                        done++;
                        if (error) ErrLog("ERROR:  exiting after error above\n");
                } // if it doesn't exist
        } while (done == 0);
  } // if writing


  // Fudge normalization to given RMS:
  double var_x,var_y,var_z;
  double avg_x=0,avg_y=0,avg_z=0;

  var_x=Variance( b_of_k_x, grid_nx*grid_ny*grid_nz,avg_x);
  var_y=Variance( b_of_k_y, grid_nx*grid_ny*grid_nz,avg_y);
  var_z=Variance( b_of_k_z, grid_nx*grid_ny*grid_nz,avg_z);

  // Fudge normalization to given RMS:
  if (bran_random_rms > 0) {
    //    if (bfield_debug) {
      Log("Original RMS of result in x-direction is "); std::cout<<sqrt(var_x)/CGS_U_muGauss<<endl;
      Log("Original RMS of result in y-direction is "); std::cout<<sqrt(var_y)/CGS_U_muGauss<<endl;
      Log("Original RMS of result in z-direction is "); std::cout<<sqrt(var_z)/CGS_U_muGauss<<endl;

          //      bran_rms_orig=(sqrt(var_x)+sqrt(var_y)+sqrt(var_z))/3.;

      Log("Resetting all to "); std::cout<<bran_random_rms/CGS_U_muGauss<<endl<<std::flush;
      //    }


  /*  if (fabs(avg_x)/sqrt(var_x/grid_nx) > 5.) { Log("ERROR:  random component in x direction doesn't appear to have zero mean!\n"); error=true;}
  if (fabs(avg_y)/sqrt(var_y/grid_ny) > 5.) { Log("ERROR:  random component in y direction doesn't appear to have zero mean!\n"); error=true;}
  if (fabs(avg_z)/sqrt(var_z/grid_nz) > 5.) { Log("ERROR:  random component in z direction doesn't appear to have zero mean!\n"); error=true;}
  */

      if (!error) {	
	for (int i=0;i<grid_nx;i++) {
	  for (int j=0;j<grid_ny;j++) {
	    for (int l=0;l<grid_nz;l++) {
	      int index=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),i,j,l);
	      b_of_k_x[index]*=bran_random_rms/sqrt(var_x);
	      b_of_k_y[index]*=bran_random_rms/sqrt(var_y);
	      b_of_k_z[index]*=bran_random_rms/sqrt(var_z);
	      
	    } // l
	  } // j
	} // i

      }
  }
  else {

    Log("RMS of result in x-direction is "); std::cout<<sqrt(var_x)/CGS_U_muGauss<<endl<<std::flush;
    Log("RMS of result in y-direction is "); std::cout<<sqrt(var_y)/CGS_U_muGauss<<endl<<std::flush;
    Log("RMS of result in z-direction is "); std::cout<<sqrt(var_z)/CGS_U_muGauss<<endl<<std::flush;

  }


} // fillRandom
//---------------------------------------------------------------







bool B_field::write_B_total() {

  if (btot_file.compare("") != 0){
    if (btot_file.compare("")==0) btot_file="debug_error_btot.bin";
    Log(std::string("Writing total B to ")+btot_file+"\n");
    
    std::stringstream Bx_name;
    std::stringstream By_name;
    std::stringstream Bz_name;
    
    Bx_name<<btot_file<<"_x.bin";
    By_name<<btot_file<<"_y.bin";
    Bz_name<<btot_file<<"_z.bin";
    
    std::ofstream Bx_out;
    std::ofstream By_out;
    std::ofstream Bz_out;
    
    //step size within the box
    double dx=grid_lx/grid_nx;
    double dy=grid_ly/grid_ny;
    double dz=grid_lz/grid_nz;
    
    vec3 los;
    vec3 point;
    vec3 Bcart;
    
    double tmp;
    
    Bx_out.open(Bx_name.str().data(),std::ios::out|std::ios::binary);
    By_out.open(By_name.str().data(),std::ios::out|std::ios::binary);
    Bz_out.open(Bz_name.str().data(),std::ios::out|std::ios::binary);
    
    
    for (int i=0;i<grid_nx;i++)
      for (int j=0;j<grid_ny;j++)
	for (int l=0;l<grid_nz;l++) {
	  
	  //in x-coordinate start at the front of the box and move to the back
	  //in y- and z- start in the centre
	  int ii=i;
	  int jj=j-grid_ny/2;
	  int ll=l-grid_nz/2;
	  //point within a box centred at a particular los
	  los=vec3(ii*dx,jj*dy,ll*dz);
	  
	  //find the coordinates of that point in the sun-centred frame
	  LOS2Cart(tlat, tlon,los,point);
	  
	  //transform the x-coordinate back to the galactocentric frame
	  //point.x=-1*(point.x+SunPosition.Length());
	  point.x=point.x+SunPosition.x;
	  Bcart = return_B_cart (point);
	  //cout << " point.x=" << point.x << " point.y=" << point.y << " point.z=" << point.z << endl;
	  //cout << " los.x=" << los.x << " los.y=" << los.y << " los.z=" << los.z << endl;
	  //cout << " Bcart.x=" << Bcart.x << " Bcart.y=" << Bcart.y << " Bcart.z=" << Bcart.z << endl;
	  
	  tmp=Bcart.x/1.e-6*CGS_U_muGauss;
	  Bx_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
	  tmp=Bcart.y/1.e-6*CGS_U_muGauss;
	  By_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
	  tmp=Bcart.z/1.e-6*CGS_U_muGauss;
	  Bz_out.write( reinterpret_cast<char*>(&tmp),sizeof(double));
	  
	}
    
    Bx_out.close();
    By_out.close();
    Bz_out.close();
    
    //if (error) ErrLog("ERROR:  exiting after error above\n");
    
  } 
  
  // if file not specified, do nothing
  else {return false;}
  
  //If B_only parameter is set in parameter file returns true
  return Bonly;
  
}




//---------------------------------------------------------------
vec3 B_field::return_B_cart (vec3 coords) {
  if (bfield_quiet != 0 && int(b_cnt) % bfield_quiet == 0) {
    Log("B_field on "); cout<<b_cnt/1e4<<"e4'th call to return_B_cart\n";
  }
  b_cnt++;
  vec3 regular=return_Breg_cart(coords);
  if (!dorand) {
    return regular;
  }// no random
  else {
    vec3 random=return_Bran_cart(coords);
    return random+regular;
  }


}// returnBcartesian
//---------------------------------------------------------------






/* ---------------------------------------------------------------
   Let integer index for each bin refer to the center of the bin.
   So 0th is actually at coordinates dx/2, and ith at (i+1/2)*dx
   (both plus offset of lx/2).

   --------------------------------------------------------------- */
vec3 B_field::return_Bran_cart(vec3 coords) {

  double r=sqrt(coords.x*coords.x + coords.y*coords.y);
  if (r > bran_rmax || (grid_nz > 1 && abs(coords.z) > 0.49*grid_lz)) { return vec3(0,0,0);}

  double r_sun=(coords-SunPosition).Length();   //sqrt((coords.x+SunPosition.Length())*(coords.x+SunPosition.Length()) + coords.y*coords.y);
  if (r_sun > ec_r_max) return vec3(0,0,0);

  double dx=grid_lx/grid_nx;
  double dy=grid_ly/grid_ny;
  double dz=grid_lz/grid_nz;

  // By default, coords in (x,y,z) where Sun-toward-Galactic Center is
  // x, with zero at GC, Sun to north galactic pole is z, ditto, etc.
  // B-random gridded the same way.  BUT if tlon within +-360,
  // i.e. reasonable instead of default -999, then the grid is meant
  // to start at the Sun and extend only in one direction.  So rotate
  // the input coordinates into the system the grid uses, but do the
  // shift to the edge below.
  bool transform=false;

  if (tlon < 2.*CGS_U_pi && tlon > -2.*CGS_U_pi) {
    transform=true;
    vec3 coords2;
    coords=coords-SunPosition;
    Cart2LOS(tlat,tlon,coords,coords2);

    // NOTE: below, also have to rotate the field vectors, because
    // they are still relative to the original galacto-centric system. 

    coords=coords2;
  }

  int xl,yl,zl;
  double tmp= transform ? coords.x/dx : (coords.x + grid_lx/2.)/dx;
  if (tmp < 0 && abs(tmp) > 1.e-10 ) xl=-1;
  else xl= int(floor(tmp));

  tmp=(coords.y + grid_ly/2.)/dy;
  if (tmp < 0 && abs(tmp) > 1.e-10) yl=-1;
  else yl= int(floor(tmp));

  tmp= (coords.z + grid_lz/2.)/dz;
  if (tmp < 0 && abs(tmp) > 1.e-10) zl=-1;
  else zl= int(floor(tmp));

  if (xl < 0 || yl < 0 || zl < 0 ) { return vec3(0,0,0);}
  if (xl > grid_nx-1 || yl>grid_ny-1 || zl>grid_nz-1) { return vec3(0,0,0);}

  vec3 result;

  // At the borders, the linear interpolation doesn't work, so there,
  // just return the edge value itself.  Can't make much difference,
  // especially since the field is generally zeroed well before the
  // edge.

  // trilinear interpolation (see e.g. Wikipedia):
  //
  // Second check is not for out-of-bounds but for too close to edge for interp.
  if (grid_interp && xl+1<grid_nx && yl+1<grid_ny && zl+1<grid_nz) {
    int idx1,idx2;
    vec3 i1,i2,j1,j2,w1,w2;
    double xd= transform ? (coords.x - xl*grid_lx/grid_nx)/dx  : (coords.x +grid_lx/2 - xl*grid_lx/grid_nx)/dx;
    double yd=(coords.y +grid_ly/2 - yl*grid_ly/grid_ny)/dy;
    double zd=(coords.z +grid_lz/2 - zl*grid_lz/grid_nz)/dz;

    idx1=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl,yl,zl);
    idx2=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl,yl,zl+1);
    i1=vec3(b_of_k_x[idx1]*(1.-zd) + b_of_k_x[idx2]*zd,
            b_of_k_y[idx1]*(1.-zd) + b_of_k_y[idx2]*zd,
            b_of_k_z[idx1]*(1.-zd) + b_of_k_z[idx2]*zd );
    idx1=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl,yl+1,zl);
    idx2=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl,yl+1,zl+1);
    i2=vec3(b_of_k_x[idx1]*(1.-zd) + b_of_k_x[idx2]*zd,
            b_of_k_y[idx1]*(1.-zd) + b_of_k_y[idx2]*zd,
            b_of_k_z[idx1]*(1.-zd) + b_of_k_z[idx2]*zd);
    idx1=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl+1,yl,zl);
    idx2=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl+1,yl,zl+1);
    j1=vec3(b_of_k_x[idx1]*(1.-zd) + b_of_k_x[idx2]*zd,
            b_of_k_y[idx1]*(1.-zd) + b_of_k_y[idx2]*zd,
            b_of_k_z[idx1]*(1.-zd) + b_of_k_z[idx2]*zd);
    idx1=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl+1,yl+1,zl);
    idx2=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl+1,yl+1,zl+1);
    j2=vec3(b_of_k_x[idx1]*(1.-zd) + b_of_k_x[idx2]*zd,
            b_of_k_y[idx1]*(1.-zd) + b_of_k_y[idx2]*zd,
            b_of_k_z[idx1]*(1.-zd) + b_of_k_z[idx2]*zd);

    w1=i1*(1.-yd) + i2*yd;
    w2=j1*(1.-yd) + j2*yd;

    result= w1*(1.-xd) + w2*xd;

    if (result.Length() > 5000.*CGS_U_muGauss) {
      cout<<"WARNING:  got big B-field for xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl;
      cout<<"result.Length() = "<< result.Length()/CGS_U_muGauss << " uG"<<endl;}
    // end of interpolation
  } else {

   int idx=Index3d(grid_nx,grid_ny,2*(grid_nz/2+1),xl,yl,zl);

    result= vec3(b_of_k_x[idx],b_of_k_y[idx],b_of_k_z[idx]);
    if (result.Length() > 5000.*CGS_U_muGauss) {
      cout<<"WARNING:  got big B-field for xl="<<xl<<", yl="<<yl<<", zl="<<zl<<endl;
      cout<<"result.Length() = "<< result.Length()/CGS_U_muGauss << " uG"<< result.Length() << endl;}
  } // no interpolation


  //  If transforming to a different longitude, need to rotate.  Note
  //  that this only rotates in x-y plane.  The idea is a field
  //  pointing toward you when at lon=0 (all x) is still pointing
  //  toward you when you move the object to lon=90 (all y).
  //  A sign flip I don't understand but empirically, this mostly works.
  if (transform) {
    vec3 tmp;
    LOS2Cart(tlat,tlon,result,tmp);
    result=tmp;
  }

  return result;

}
//---------------------------------------------------------------



//---------------------------------------------------------------
vec3 B_field::return_Breg_cart( vec3 coords)
{
  vec3 regular;
  switch ( bfield_type) {
  case 1:  regular = field1(coords); break;
  case 2:  regular = field2(coords); break;
  case 3:  regular = field3(coords); break;
  case 4:  regular = field4(coords); break;
  case 5:  regular = field5(coords); break;
  case 6:  regular = field6(coords); break;
  case 7:  regular = field7(coords); break; //Added by SPQ Aug 12 2014
  case 8:  regular = field8(coords); break; //Added by SPQ Aug 15 2014
  case 9:  regular = field9(coords); break; //Added by SPQ Aug 15 2014
  case 10:  regular = field10(coords); break; //Added by SPQ Aug 15 2014
    // Your field here:
    //  case 10: regular=field10(coords); break;
  default: cerr << " No bfield_type specified " << endl; exit(1); break;
  }
  return regular;
}// return_Breg_cart
//---------------------------------------------------------------

// To be called in the case that a single B_field object is used
// with multiple random realisations.  Run this between calls to fillRandom.
void B_field::cleanRandom(){
          if (dorand) {
            Log("B_field freeing up random component\n");
            fftw_free(b_of_k_x);
            fftw_free(b_of_k_y);
            fftw_free(b_of_k_z);
          }
}

void B_field::SanityCheck( ) {
  if (!dorand) return;
  cout<<endl<<endl<<"----------------------"<<endl<<"HAMMURABI Sanity Check:"<<endl<<endl;

  // Check random grid bin sizes:
  double bran_deltax=grid_lx/grid_nx;
  double bran_deltay=grid_ly/grid_ny;
  double bran_deltaz=grid_lz/grid_nz;
  cout<<"Random grid deltaX="<<bran_deltax/CGS_U_kpc<<" kpc"<<endl;
  cout<<"Random grid deltaY="<<bran_deltay/CGS_U_kpc<<" kpc"<<endl;
  cout<<"Random grid deltaZ="<<bran_deltaz/CGS_U_kpc<<" kpc"<<endl;

  //  if (bran_deltax != bran_deltay || bran_deltay != bran_deltaz) cout<<"WARNING:  your random grid does not have the same resolution in all directions."<<endl;
  if (abs(bran_deltax-bran_deltay)/bran_deltax > 0.01  || abs(bran_deltay-bran_deltaz)/bran_deltay > 0.01) cout<<"WARNING:  your random grid does not have the same resolution in all directions."<<endl;


  // If using a transform, check that the random grid covers the
  // longitude range given

  bool transform= (tlon < 2.*CGS_U_pi && tlon > -2.*CGS_U_pi) ? true : false ; 

  if (transform) {
    double deltalon=grid_ly/grid_lx;
    double grid_lmin=(tlon - 0.5*deltalon)*180./CGS_U_pi;
    double grid_lmax=(tlon + 0.5*deltalon)*180./CGS_U_pi;

    cout<<"Transform along longitude "<<tlon*180./CGS_U_pi<<" deg"<<endl;
    cout<<"Grid from longitude "<<grid_lmin<<" to "<<grid_lmax<<" deg"<<endl;

  }

  // Total memory required:
  double bytes = grid_nx*grid_ny*grid_nz*3.*8.;
  if (bytes > bran_mem_limit*1.e9) {
    cout<<"ERROR:  you're asking for "<<bytes/1.e9<<" GB of memory.  Might want to rethink or disable this if you're sure you can"<<endl;
    exit(-1);
  } else  cout<<"Memory for grid:  "<<bytes/1.e6<<" MB"<<endl;


  cout<<endl<<"HAMMURABI Sanity Check done"<<endl<<"----------------------"<<endl<<endl;
}




void B_field::SanityCheck( Integrator &intor, double lmin, double lmax) {
  if (!dorand) return;

  cout<<endl<<endl<<"----------------------"<<endl<<"HAMMURABI Sanity Check:"<<endl<<endl;
  double R_max=intor.MAX_RADIUS;
  int vec_size_R=intor.vec_size_R;
  double delta_R=R_max/vec_size_R;
  int obs_nside=intor.obs_NSIDE;
  int obs_shell=intor.obs_shell;
  int total_shell=intor.total_shell;

  cout<<"DeltaR="<<delta_R/CGS_U_kpc<<" kpc"<<endl;

  // Size perp to LOS of last integration box, which isn't the
  // smallest but is good to check:
  int max_nside=obs_nside; for (int s=obs_shell+1;s<=total_shell;s++) { max_nside*=2;}
  int npix_max=12*max_nside*max_nside;
  double dtheta=sqrt( 4.*CGS_U_pi/npix_max );
  double deltax=dtheta*R_max;
  cout<<"At R_max (max Nside="<<max_nside<<"), DeltaX="<<deltax/CGS_U_kpc<<" kpc"<<endl;

  // Check random grid bin sizes:
  double bran_deltax=grid_lx/grid_nx;
  double bran_deltay=grid_ly/grid_ny;
  double bran_deltaz=grid_lz/grid_nz;
  cout<<"Random grid deltaX="<<bran_deltax/CGS_U_kpc<<" kpc"<<endl;
  cout<<"Random grid deltaY="<<bran_deltay/CGS_U_kpc<<" kpc"<<endl;
  cout<<"Random grid deltaZ="<<bran_deltaz/CGS_U_kpc<<" kpc"<<endl;

  //  if (bran_deltax != bran_deltay || bran_deltay != bran_deltaz) cout<<"WARNING:  your random grid does not have the same resolution in all directions."<<endl;
  if (abs(bran_deltax-bran_deltay)/bran_deltax > 0.01  || abs(bran_deltay-bran_deltaz)/bran_deltay > 0.01) cout<<"WARNING:  your random grid does not have the same resolution in all directions."<<endl;

  // Compare random grid to integration resolution:
  if ( abs( bran_deltax - deltax )/deltax > 2) cout<<"WARNING:  your random grid is much coarser than your integration resolution at R_max"<<endl;

  if ( abs( bran_deltax - deltax )/bran_deltax > 2) cout<<"WARNING:  your random grid is much finer than your integration resolution at R_max"<<endl;


  // If using a transform, check that the random grid covers the
  // longitude range given

  bool transform= (tlon < 2.*CGS_U_pi && tlon > -2.*CGS_U_pi) ? true : false ; 
  if (transform) {
    double deltalon=grid_ly/R_max;
    double grid_lmin=(tlon - 0.5*deltalon)*180./CGS_U_pi;
    double grid_lmax=(tlon + 0.5*deltalon)*180./CGS_U_pi;

    cout<<"Transform along longitude "<<tlon*180./CGS_U_pi<<" deg"<<endl;
    cout<<"Grid from longitude "<<grid_lmin<<" to "<<grid_lmax<<" deg"<<endl;
    if (!(lmin == 0 && lmax==0)) cout<<"Pixel list from "<<lmin<<" to "<<lmax<<" deg"<<endl;

  }

  // Total memory required:
  double bytes = grid_nx*grid_ny*grid_nz*3.*8.;
  if (bytes > bran_mem_limit*1.e9) {
    cout<<"ERROR:  you're asking for "<<bytes/1.e9<<" GB of memory.  Might want to rethink or disable this if you're sure you can"<<endl;
    exit(-1);
  } else  cout<<"Memory for grid:  "<<bytes/1.e6<<" MB"<<endl;



  cout<<endl<<"HAMMURABI Sanity Check done"<<endl<<"----------------------"<<endl<<endl;
    //  cout<<"----------------------"<<endl<<"HAMMURABI Sanity Check done;  exiting"<<endl<<endl;
  //  exit(0);
}











/*---------------------------------------------------------------
  Additional field functions:

  Search for "Your field here" above and in proto_B_field2.cpp for
  where to add your own field definition(s).

  1) The parameters and declarations for the setup and field functions
  must be added to the proto.

  2) Here (or in another file, but then add it appropriately to the
  Makefile), code the setup and field functions for your field(s).

  3) Above, add them to the switch on the field type and to the
  read_B_params function.

  ---------------------------------------------------------------*/

/* Your field here:
//---------------------------
void B_field::setup_field10(){
}
vec3 B_field::field10(vec3 coords) {
 vec3 B_cart;
 ...
 return B_cart;
}

//---------------------------
*/
