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
#include <hammurabi.h>
#include "proto_namespace_Vec_Handling.h"
#include "proto_class_CRE.h"

using namespace std;

// from the NE2001 data.
const double thick_disk_radius=17.*CGS_U_kpc;

//-------------------------------------------------------------------------------------------------------
//Routines of class CRE
//-------------------------------------------------------------------------------------------------------

// class CRE constructor
CRE::CRE(paramfile & params)
{
  cout << " class CRE constructor " << endl;
  read_CRE_params(params);
}

// Empty constructor. Cfield to be set up by later setup_Cfield1...
// call.
CRE::CRE()
{
  cout << " class CRE constructor " << endl;
  Cfield_type = 0; // so routines know of empty constructor

  // Sun position in GC cart coordinates default setting:
  SunPosition.x=-CGS_U_kpc*8.5;
  SunPosition.y=CGS_U_kpc*0.;
  SunPosition.z=CGS_U_kpc*0.;
}

// class CRE destructor
CRE::~CRE(void)
{
  cout << " class CRE destructor " << endl;
}

void CRE::read_CRE_params(paramfile &params)
{
  cout << " Reading parameters for class CRE " << endl;
  
  Cfield_type=params.find<int>("Cfield_type", 1);

  // Sun position in GC cartesian coordinates.
  SunPosition.x=CGS_U_kpc*params.find<double>("SunPosX", -8.5);
  SunPosition.y=CGS_U_kpc*params.find<double>("SunPosY", 0.);
  SunPosition.z=CGS_U_kpc*params.find<double>("SunPosZ", 0.);

  if(Cfield_type == 1)
    {
      double spec_index = params.find<double>("C1_p",3.);
      double hr = CGS_U_kpc*params.find<double>("C1_hr", 5.);
      double hz = CGS_U_kpc*params.find<double>("C1_hz", 1.);
      double CRE_alpha = params.find<double>("alpha", 1.);
       setup_Cfield1(spec_index, hr, hz, CRE_alpha);
    }
  else if(Cfield_type ==2)
    {
      double spec_index_p1 = params.find<double>("C2_p1", 2.);
      double spec_index_p2 = params.find<double>("C2_p2", 3.);
      double hr = CGS_U_kpc*params.find<double>("C2_hr", 8.);
      double hz = CGS_U_kpc*params.find<double>("C2_hz", 1.);
      setup_Cfield2(spec_index_p1, spec_index_p2, hr, hz);
    }
  else if(Cfield_type ==3)
    {
      double spec_index_p1 = params.find<double>("C3_p1", 2.75);
      double spec_index_p2 = params.find<double>("C3_p2", 3.);
      double hr = CGS_U_kpc*params.find<double>("C3_hr", 8.);
      double hz = CGS_U_kpc*params.find<double>("C3_hz", 1.);
      setup_Cfield3(spec_index_p1, spec_index_p2, hr, hz);
    }
  else if(Cfield_type ==4)
    {
      double spec_index = params.find<double>("C4_p",3.);
      double CRE_alpha = params.find<double>("alpha", 1.);
      setup_Cfield4(spec_index, CRE_alpha);
    }
  else if(Cfield_type == 5)
    {
      double spec_index = params.find<double>("C5_p",3.);
      double ha = CGS_U_kpc*params.find<double>("C5_ha", 0.);
      double hr = CGS_U_kpc*params.find<double>("C5_hr", 5.);
      double hz = CGS_U_kpc*params.find<double>("C5_hz", 1.);
      double CRE_alpha = params.find<double>("alpha", 1.);
      setup_Cfield5(spec_index, ha, hr, hz, CRE_alpha);
    }
  else if(Cfield_type == 6)
    {
      double spec_index = params.find<double>("C6_p",3.);
      double ha = CGS_U_kpc*params.find<double>("C6_ha", 5.);
      double hr = CGS_U_kpc*params.find<double>("C6_hr", 3.);
      double hz = CGS_U_kpc*params.find<double>("C6_hz", 1.);
      double CRE_alpha = params.find<double>("alpha", 1.);
      setup_Cfield6(spec_index, ha, hr, hz, CRE_alpha);
    }
  else if(Cfield_type == 7)
    {
      double spec_index = params.find<double>("C7_p",3.);
      double ha = 1./CGS_U_kpc*params.find<double>("C7_ha", 0.);
      double hb = 1./CGS_U_kpc*params.find<double>("C7_hb", 0.);
      double CRE_alpha = params.find<double>("alpha", 1.);
      setup_Cfield7(spec_index, ha, hb, CRE_alpha);
    }
  else if(Cfield_type == 8)
    {
      double spec_index = params.find<double>("C8_p",3.);
      double CRE_alpha = params.find<double>("alpha", 1.);
      double z_cre = CGS_U_kpc*params.find<double>("C8_zcre", 1.);
      setup_Cfield8(spec_index, z_cre, CRE_alpha);
    }
  else
    {
      cerr << "Invalid Cfield_type: " << Cfield_type << endl;
      exit(1);
    }

  cout << " Done " << endl;
}

// Obtains the C ("spatial term" of the CRE energy density
// distribution, see Ribicky&Lightman).
double CRE::get_C (double R, double THE, double PHI) const
{
  double r_gc, PHI_gc, z_gc;
  vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(R,THE,PHI);
  vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
  Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);

  if(Cfield_type==1){return cre_model_WMAP_3th_year (r_gc,z_gc);}
  else if(Cfield_type==2){return cre_model_Sun2008 (r_gc,z_gc);}
  else if(Cfield_type==3){return cre_model_Sun2008_II (r_gc,z_gc);}
  else if(Cfield_type==4){return cre_model_galprop (r_gc,z_gc);}
  else if(Cfield_type==5){return cre_model_wmap_mod (r_gc,z_gc);}
  else if(Cfield_type==6){return cre_model_JanssonFarrar2011 (r_gc,z_gc);}
  else if(Cfield_type==7){return cre_model_galprop_mod (r_gc,z_gc);}
  else if(Cfield_type==8){return cre_model_galprop_scaled (r_gc,z_gc);}
  else if(Cfield_type==0)
    {
      cerr << " Cfield_type==0! CRE class running on empty constructor!" << endl;
      cerr << " Call to setup_Cfield1... needed first! " << endl;
      exit(1);
    }
  else 
    {
      cerr << "Invalid Cfield_type: " << Cfield_type << endl;
    }
}

// this function might change a bunch in future
// it could e.g. depend on R, THE, PHI and obs_freq
double CRE::get_spec_index_p(double obs_freq, double ec_R, double ec_THE, double ec_PHI) const
{
  if(Cfield_type==1)  {return C1_spec_index_p;}
  else if(Cfield_type==4)  {return C4_spec_index_p;}
  else if(Cfield_type==5)  {return C5_spec_index_p;}
  else if(Cfield_type==6)  {return C6_spec_index_p;}
  else if(Cfield_type==7)  {return C7_spec_index_p;}
  else if(Cfield_type==8)  {return C8_spec_index_p;}
  else if(Cfield_type==2) 
    {
      if(obs_freq>=0.408*CGS_U_GHz)
	{return C2_spec_index_p2;}
      else if(obs_freq<0.408*CGS_U_GHz)
	{return C2_spec_index_p1;}
      else
	{
	  cerr << " Error in double CRE::get_spec_index_p(double obs_freq) const " << endl;
	  cerr << " something wrong with obs_freq:" << obs_freq << endl;
	  exit(1);
	}
    }
  else if(Cfield_type==3) // Sun2008_II, proof of concept case
    {
      double r_gc, PHI_gc, z_gc;
      vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(ec_R,ec_THE,ec_PHI);
      vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
      Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);
      if(ec_R>8.*CGS_U_kpc)
	{return C3_spec_index_p1;}
      else
	{return C3_spec_index_p2;}
    }
  else
    {
      cerr << " Error in double CRE::get_spec_index_p(double obs_freq) const "
	   << endl;
      cerr << " Invalid Cfield_type: " << Cfield_type << endl;
      exit(1);
    }
}

// ncre model used by Page et al. WMAP 3th year polarization paper.
double CRE::cre_model_WMAP_3th_year (double r, double z) const
{
  //double h_r=5.0*CGS_U_kpc;
  //double h_d=1.0*CGS_U_kpc;
  

  double h_r = C1_hr;  // RJ MOD!!!
  double h_d = C1_hz;  
  double CRE_alpha = alpha; 
 
  double gc_earth_radius=sqrt(dotprod(SunPosition,SunPosition));
  double C_zero=exp(-gc_earth_radius/h_r)*(2./(exp(0./h_d)+exp(-0./h_d)))*(2./(exp(0./h_d)+exp(-0./h_d)));

  // This comes from the Diploma thesis in chapter
  // relativistic electron density
  const double cre_calc_forefactor=4.*CGS_U_pi*(1./CGS_U_C_light)*0.25*CGS_U_MEC2/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec); 

  const double cre_gamma_10=10.*CGS_U_GeV/CGS_U_MEC2; //We assume, or take, E=10GeV 

  double C_earth, C_cre;
  C_earth=cre_calc_forefactor*pow(cre_gamma_10, C1_spec_index_p);	// Value obtained for earth in the thesis. The WMAP paper doesnt give any (not necessary for its calculations).
  //C_cre=(C_earth/C_zero)*exp(-r/h_r)*(2./(exp(z/h_d)+exp(-z/h_d)))*(2./(exp(z/h_d)+exp(-z/h_d)));
  C_cre=CRE_alpha*(C_earth/C_zero)*exp(-r/h_r)*(2./(exp(z/h_d)+exp(-z/h_d)))*(2./(exp(z/h_d)+exp(-z/h_d))); // RJ mod
  
  //cout << C_cre << endl;

  return C_cre;
}

// ncre model used by Sun et al. 2008.
double CRE::cre_model_Sun2008 (double r, double z) const
{
  double C_zero=(6.4e-5)/CGS_U_ccm;
  double C_cre = C_zero*exp(-((r-SunPosition.Length())/(C2_hr))-(abs(z)/(C2_hz)));

  if(r<3.*CGS_U_kpc) 
    {
      C_cre=C_zero*exp(-((3.*CGS_U_kpc-SunPosition.Length())/(C2_hr))-(abs(z)/(C2_hz)));
    }
  // Note that Sun et al. point out other 
  // options for this 1kpc truncation.
  if(abs(z)>1.*CGS_U_kpc) {C_cre=0.;}

  //cout << C_cre << endl;

  return C_cre;
}

// ncre model used by Sun et al. 2008.
double CRE::cre_model_Sun2008_II (double r, double z) const
{
  double C_zero=(6.4e-5)/CGS_U_ccm;
  double C_cre = C_zero*exp(-((r-SunPosition.Length())/(C3_hr))-(abs(z)/(C3_hz)));

  if(r<3.*CGS_U_kpc) 
    {
      C_cre=C_zero*exp(-((3.*CGS_U_kpc-SunPosition.Length())/(C3_hr))-(abs(z)/(C3_hz)));
    }
  // Note that Sun et al. point out other 
  // options for this 1kpc truncation.
  if(abs(z)>1.*CGS_U_kpc) {C_cre=0.;}

  return C_cre;
}

double CRE::cre_model_galprop (double r, double z) const
{
  double CRE_alpha = alpha; 
  double C_cre;

  // find r,z indices [hardcoded size]
  int zi = int( (abs(z)/CGS_U_kpc)*10 + 40. + 0.05 ); // +40 because of symmetry in z
  int rj = int( (r/CGS_U_kpc) + 0.5 ); 

  C_cre = CRE_alpha*C4_ncre[zi][rj]/CGS_U_ccm;

  return C_cre;
}

// modified wmap ncre model, for parameters  ha=5.62,  hr=5.34, hz=1.99 it is a pretty good spatial fit to Galprop model
double CRE::cre_model_wmap_mod (double r, double z) const
{

  double CRE_alpha = alpha; 
  double gc_earth_radius=sqrt(dotprod(SunPosition,SunPosition));
  double C_zero=exp(-abs(gc_earth_radius-C5_ha)/C5_hr)*(2./(exp(0./C5_hz)+exp(-0./C5_hz)))*(2./(exp(0./C5_hz)+exp(-0./C5_hz)));

  // This comes from the Diploma thesis in chapter
  // relativistic electron density
  const double cre_calc_forefactor=4.*CGS_U_pi*(1./CGS_U_C_light)*0.25*CGS_U_MEC2/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec); 

  const double cre_gamma_10=10.*CGS_U_GeV/CGS_U_MEC2; //We assume, or take, E=10GeV 

  double C_earth, C_cre;
  C_earth=cre_calc_forefactor*pow(cre_gamma_10, C5_spec_index_p);	// Value obtained for earth in the thesis. The WMAP paper doesnt give any (not necessary for its calculations).
  C_cre=alpha*(C_earth/C_zero)*exp(-abs(r-C5_ha)/C5_hr)*(2./(exp(z/C5_hz)+exp(-z/C5_hz)))*(2./(exp(z/C5_hz)+exp(-z/C5_hz))); 


  return C_cre;
}


// model with exponential radial fall off after r=ha, and constant for r<ha.
double CRE::cre_model_JanssonFarrar2011 (double r, double z) const
{

  double CRE_alpha = alpha; 
  double gc_earth_radius=sqrt(dotprod(SunPosition,SunPosition));
  double C_zero = 0.;
  if (gc_earth_radius>C6_ha)
  {
    C_zero=exp(-(gc_earth_radius-C6_ha)/C6_hr)*(2./(exp(0./C6_hz)+exp(-0./C6_hz)))*(2./(exp(0./C6_hz)+exp(-0./C6_hz)));
  }
  else 
  {
    C_zero=(2./(exp(0./C6_hz)+exp(-0./C6_hz)))*(2./(exp(0./C6_hz)+exp(-0./C6_hz)));
  }

  // This comes from the Diploma thesis in chapter
  // relativistic electron density
  const double cre_calc_forefactor=4.*CGS_U_pi*(1./CGS_U_C_light)*0.25*CGS_U_MEC2/(CGS_U_GeV*100.*CGS_U_cm*100.*CGS_U_cm*CGS_U_sec); 

  const double cre_gamma_10=10.*CGS_U_GeV/CGS_U_MEC2; //We assume, or take, E=10GeV 

  double C_earth, C_cre;
  C_earth=cre_calc_forefactor*pow(cre_gamma_10, C6_spec_index_p);	// Value obtained for earth in the thesis. The WMAP paper doesnt give any (not necessary for its calculations).

  if (r>C6_ha)
  {
  C_cre=alpha*(C_earth/C_zero)*exp(-(r-C6_ha)/C6_hr)*(2./(exp(z/C6_hz)+exp(-z/C6_hz)))*(2./(exp(z/C6_hz)+exp(-z/C6_hz))); 
  }
  else
  {
  C_cre=alpha*(C_earth/C_zero)*(2./(exp(z/C6_hz)+exp(-z/C6_hz)))*(2./(exp(z/C6_hz)+exp(-z/C6_hz))); 
  }

  return C_cre;
}

double CRE::cre_model_galprop_mod (double r, double z) const
{
  double CRE_alpha = alpha; 
  double C_cre;
  double h_a = C7_ha;  
  double h_b = C7_hb;  
  double ncre_scaling = exp(h_a*(abs(SunPosition.x)-r)+h_b*abs(z)); // scaling that leaves n_cre(Earth) unchanged.

  // find r,z indices [hardcoded size]
  int zi = int( (abs(z)/CGS_U_kpc)*10 + 40. + 0.05 ); // +40 because of symmetry in z
  int rj = int( (r/CGS_U_kpc) + 0.5 ); 

  C_cre = CRE_alpha*ncre_scaling*C4_ncre[zi][rj]/CGS_U_ccm;
  return C_cre;
}

// modified GALPROP model, by a scaling factor
double CRE::cre_model_galprop_scaled (double r, double z) const
{
  double CRE_alpha = alpha; 
  double C_cre; 

  double ncre_scaling = exp(-abs(z)/C8_zcre); 

  // find r,z indices [hardcoded size]
  int zi = int( (abs(z)/CGS_U_kpc)*10 + 40. + 0.05 ); // +40 because of symmetry in z
  int rj = int( (r/CGS_U_kpc) + 0.5 ); 

  C_cre = CRE_alpha*ncre_scaling*C8_ncre[zi][rj]/CGS_U_ccm;
  return C_cre;
}




// Routines to set up the Cfield values:

void CRE::setup_Cfield1(double p, double hr, double hz, double CRE_alpha)
{
  Log("Setting up Page et. al 2007 C field model");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=1;}
  C1_spec_index_p=p;
  C1_hr=hr;
  C1_hz=hz;
  alpha=CRE_alpha;
}

void CRE::setup_Cfield2(double p1, double p2, double hr, double hz)
{
  Log("Setting up Sun et. al 2008 C field model");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=2;}
  C2_spec_index_p1=p1;
  C2_spec_index_p2=p2;
  C2_hr=hr;
  C2_hz=hz;
}

void CRE::setup_Cfield3(double p1, double p2, double hr, double hz)
{
  Log("Setting up Sun et. al 2008 II C field model");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=3;}
  C3_spec_index_p1=p1;
  C3_spec_index_p2=p2;
  C3_hr=hr;
  C3_hz=hz;
}
void CRE::setup_Cfield4(double p, double CRE_alpha)
{
  Log("Setting up n_cre model from GALPROP input\n");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=4;}
  C4_spec_index_p=p;
  alpha=CRE_alpha;
  
  // create look-up table for GALPROP data
  string galprop_data_filename = "ncre_galprop.txt";
  ifstream galprop_data(galprop_data_filename.c_str() ); 
  if (!galprop_data) {
    cerr << "   Could not open galprop data file: " << galprop_data_filename << endl;
    exit(EXIT_FAILURE);
  }
  int zi = 81; // list indices, in r (z comes in 0.1 kpc steps from -4 to 4 kpc, r in 1 kpc steps from 0 to 20 kpc)
  int rj = 21;
  for (int i=0; i<zi; i++){
    for (int j=0; j<rj; j++){
      galprop_data >> C4_ncre[i][j];
    }
  }
}

void CRE::setup_Cfield5(double p, double ha, double hr, double hz, double CRE_alpha)
{
  Log("Setting up modified WMAP n_cre model");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=5;}
  C5_spec_index_p=p;
  C5_ha=ha;
  C5_hr=hr;
  C5_hz=hz;
  alpha=CRE_alpha;
}

void CRE::setup_Cfield6(double p, double ha, double hr, double hz, double CRE_alpha)
{
  Log("Setting up Jansson-Farrar (2011) ncre model ");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=6;}
  C6_spec_index_p=p;
  C6_ha=ha;
  C6_hr=hr;
  C6_hz=hz;
  alpha=CRE_alpha;
}


void CRE::setup_Cfield7(double p, double ha, double hb, double CRE_alpha)
{
  Log("Setting up (modified) n_cre model from GALPROP input\n");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=7;}
  C7_spec_index_p=p;
  C7_ha=ha;
  C7_hb=hb;
  alpha=CRE_alpha;
  
  // create look-up table for GALPROP data
  string galprop_data_filename = "ncre_galprop.txt";
  ifstream galprop_data(galprop_data_filename.c_str() ); 
  if (!galprop_data) {
    cerr << "   Could not open galprop data file: " << galprop_data_filename << endl;
    exit(EXIT_FAILURE);
  }
  int zi = 81; // list indices, in r (z comes in 0.1 kpc steps from -4 to 4 kpc, r in 1 kpc steps from 0 to 20 kpc)
  int rj = 21;
  for (int i=0; i<zi; i++){
    for (int j=0; j<rj; j++){
      galprop_data >> C4_ncre[i][j];
    }
  }
}

void CRE::setup_Cfield8(double p, double z_cre, double CRE_alpha)
{
  Log("Setting up re-scaled GALPROP n_cre\n");
  // for empty constuctor Cfield_type has to be set.
  if(Cfield_type==0) {Cfield_type=8;}
  C8_spec_index_p=p;
  alpha=CRE_alpha;
  C8_zcre=z_cre;
  // create look-up table for GALPROP data
  string galprop_data_filename = "ncre_galprop.txt";
  ifstream galprop_data(galprop_data_filename.c_str() ); 
  if (!galprop_data) {
    cerr << "   Could not open galprop data file: " << galprop_data_filename << endl;
    exit(EXIT_FAILURE);
  }
  int zi = 81; // list indices, in r (z comes in 0.1 kpc steps from -4 to 4 kpc, r in 1 kpc steps from 0 to 20 kpc)
  int rj = 21;
  for (int i=0; i<zi; i++){
    for (int j=0; j<rj; j++){
      galprop_data >> C8_ncre[i][j];
    }
  }
}



