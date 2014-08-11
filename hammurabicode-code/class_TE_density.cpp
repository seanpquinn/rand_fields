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

#include "proto_class_TE_density.h"
#include <cmath>
#include "hammurabi.h"
#include "proto_namespace_Vec_Handling.h"

/*  NEW WAY:  include cfortran.h in header.  Then following instructions at
    http://wwwasd.web.cern.ch/wwwasd/cernlib/cfortran.html
                 PROTOCCALLSFSUB2(SUB_NAME,sub_name,STRING,PINT)
#define SUB_NAME(A,B) CCALLSFSUB2(SUB_NAME,sub_name,STRING,PINT, A,B)

                                ^     -                                       -
   number of arguments _____|    |   STRING   BYTE    PBYTE       BYTEV(..)|
                              /  |   STRINGV  DOUBLE  PDOUBLE   DOUBLEV(..)|
                             /   |  PSTRING   FLOAT   PFLOAT     FLOATV(..)|
    types of arguments ____ /    | PNSTRING   INT     PINT         INTV(..)|
                            \    | PPSTRING   LOGICAL PLOGICAL LOGICALV(..)|
                             \   |  PSTRINGV  LONG    PLONG       LONGV(..)|
                              \  |   ZTRINGV  SHORT   PSHORT     SHORTV(..)|
                                 |  PZTRINGV  ROUTINE PVOID      SIMPLE    |

*/
#include "cfortran.h"

PROTOCCALLSFSUB10(DMDSM,dmdsm,\
                  FLOAT,FLOAT,INT,PFLOAT,PFLOAT,PSTRING,PFLOAT,PFLOAT,PFLOAT,PFLOAT)
#define DMDSM(L,B,NDIR,DM,R,LIMIT,SM,SUNGCDISTANCE,SMTHETA,SMISO) CCALLSFSUB10(DMDSM,dmdsm, \
                  FLOAT,FLOAT,INT,PFLOAT,PFLOAT,PSTRING,PFLOAT,PFLOAT,PFLOAT,PFLOAT, \
                  L,B,NDIR,DM,R,LIMIT,SM,SUNGCDISTANCE,SMTHETA,SMISO)

using namespace std;

//should get ne calling directly the fortran routine.
double TE_density::get_density(double R, double THE, double PHI)
 {

  if (TE_constant > 1.e-20) return TE_constant;

  if (!tegrid) {

	// fortran code variables
	int ndir=-1;
	char limit[] = {0,0};
	float sungcdistance=SunPosition.Length()/CGS_U_kpc;

	float l=PHI;
	float b=(CGS_U_pi/2.)-THE;
	float r=R/CGS_U_kpc; //r is the radius in kpc to be forwarded to the fortran routine
	//-1 factor is to call the proper dmdsm routine. Otherwise it will try to return a dispersion measure DM.	
	//	dmdsm_(&l,&b,&ndir,&dm,&r,&limit,&sm,&sungcdistance,&smtheta,&smiso); 
        float dm, sm, smtheta, smiso;
	DMDSM(l,b,ndir,dm,r,limit,sm,sungcdistance,smtheta,smiso);

        return sm/CGS_U_ccm;


  } // if no input grid
  else {


    // Coordinates given (r,theta,phi) in Sun-centric and hammurabi
    // units.  Convert to galacto-centric (x,y,z) in kpc and then interpolate
    // from grid.
    double lat=CGS_U_pi/2. - THE;

    double x=(R*cos(PHI)*cos(lat)-SunPosition.Length())/CGS_U_kpc;
    double y=R*sin(PHI)*cos(lat)/CGS_U_kpc;
    double z=R*sin(lat)/CGS_U_kpc;

    // Return zero if out of bounds of the grid.  (Why doesn't fabs() work?)
    if (x > 0.5*TE_Lx || y > 0.5*TE_Ly || z > 0.5*TE_Lz)  return 0.;
    if (x < -0.5*TE_Lx || y < -0.5*TE_Ly || z < -0.5*TE_Lz)  return 0.;

    double dx=TE_Lx/TE_nx;
    double dy=TE_Ly/TE_ny;
    double dz=TE_Lz/TE_nz;

    int xl,yl,zl;
    double tmp=(x + 0.5*TE_Lx)/dx;
    if (tmp < 0 && fabs(tmp) > 1.e-10) xl=-1;
    else xl=int(round(tmp));

    tmp=(y + 0.5*TE_Ly)/dy;
    if (tmp < 0 && fabs(tmp) > 1.e-10) yl=-1;
    else yl=int(round(tmp));

    tmp=(z + 0.5*TE_Lz)/dz;
    if (tmp < 0 && fabs(tmp) > 1.e-10) zl=-1;
    else zl=int(round(tmp));

    if (xl < 0 || yl < 0 || zl < 0 ) { return 0;}
    if (xl > TE_nx-1 || yl> TE_ny-1 || zl> TE_nz-1) { return 0;}


    // trilinear interpolation (see e.g. Wikipedia):

    // Second check is not for out-of-bounds but for too close to edge for interp
    if ( xl+1<TE_nx && yl+1<TE_ny && zl+1<TE_nz) {

      double i1,i2,j1,j2,w1,w2;
      double xd= (x+0.5*TE_Lx - xl*TE_Lx/TE_nx)/dx;
      double yd= (y+0.5*TE_Ly - yl*TE_Ly/TE_ny)/dy;
      double zd= (z+0.5*TE_Lz - zl*TE_Lz/TE_nz)/dz;

      i1=TE_grid(xl,yl,zl)*(1.-zd) + TE_grid(xl,yl,zl+1)*zd;
      i2=TE_grid(xl,yl+1,zl)*(1-zd) + TE_grid(xl,yl+1,zl+1)*zd;
      j1=TE_grid(xl+1,yl,zl)*(1-zd) + TE_grid(xl+1,yl,zl+1)*zd;
      j2=TE_grid(xl+1,yl+1,zl)*(1-zd) + TE_grid(xl+1,yl,zl+1)*zd;

      w1=i1*(1-yd)+i2*yd;
      w2=j1*(1-yd)+j2*yd;

      return (w1*(1-xd)+w2*xd)/CGS_U_ccm;

    } // if not edge
    else
      return TE_grid(xl,yl,zl)/CGS_U_ccm;

  } // if given grid file


} // get_density


TE_density::TE_density(paramfile &params) {

  TE_constant=params.find<double>("TE_constant_pccm",0.)/CGS_U_ccm;

  // Sun position in GC cartesian coordinates.
  SunPosition.x=CGS_U_kpc*params.find<double>("SunPosX", -8.5);
  SunPosition.y=CGS_U_kpc*params.find<double>("SunPosY", 0.);
  SunPosition.z=CGS_U_kpc*params.find<double>("SunPosZ", 0.);

  string TE_grid_filename=params.find<string>("TE_grid_filename","");
  if (TE_grid_filename.compare("") != 0) {
    tegrid=true;
    TE_Lx=params.find<double>("TE_lx_kpc",40);
    TE_Ly=params.find<double>("TE_ly_kpc",40);
    TE_Lz=params.find<double>("TE_lz_kpc",8);
    TE_nx=params.find<int>("TE_nx",200);
    TE_ny=params.find<int>("TE_ny",200);
    TE_nz=params.find<int>("TE_nz",40);

    cout<<"Allocating memory for TE grid\n";
    TE_grid.alloc(TE_nx,TE_ny,TE_nz);
    TE_grid.fill(0.);

    cout<<"Reading TE grid\n";
    read_grid(TE_grid_filename);

  } // if grid filename given
  else { tegrid=false;}
}



void TE_density::read_grid(const std::string &filename) {

  std::ifstream input(filename.c_str(), std::ios::in|std::ios::binary);
  planck_assert(input,"Could not open input TE file "+filename);
  float tmp;

  for (int i=0; i<TE_nx; i++)
    for (int j=0; j<TE_ny; j++)
      for (int k=0; k<TE_nz; k++) {
	if (input.eof()) { std::cout<<"TE file ended sooner than expected at i="<<i<<",j="<<j<<",k="<<k<<"!\n"; exit(1);}
	input.read(reinterpret_cast<char *>(&tmp),sizeof(float));
	TE_grid(i,j,k)=tmp;
      }
  if (!input.eof()) {
    input.read(reinterpret_cast<char *>(&tmp),sizeof(float));
    if (!input.eof()) { cout<<"TE file larger than expected at i="<<TE_nx<<",j="<<TE_ny<<",k="<<TE_nz<<"!\n"; exit(1);}
  }

  cout<<"TE file successfully read\n";

}

// See Sun et al. 2008 and references therein
double TE_density::get_temperature(double R, double THE, double PHI)
{
  double r_gc, PHI_gc, z_gc;
  vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(R,THE,PHI);
  vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
  Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);

  double Ti;
  // for the formulae in Sun et al. 2008 we require
  // r_gc, and z_gc to be in units of kpc.
  z_gc=abs(z_gc/CGS_U_kpc);
  r_gc/=CGS_U_kpc;
  Ti=5780.+287.*r_gc-562.*abs(z_gc)+1770.*abs(z_gc*z_gc);
  return Ti;
}

// See Sun et al. 2008 and references therein
double TE_density::get_filling_factor(double R, double THE, double PHI)
{
double r_gc, PHI_gc, z_gc;
  vec3 cart_vec = Vec_Handling::RTHEPHI2cart_vec(R,THE,PHI);
  vec3 gc_cart_vec = Vec_Handling::ec_cart_vec2gc_cart_vec(cart_vec, SunPosition);
  Vec_Handling::cart_coord2cyl_coord(gc_cart_vec,r_gc,PHI_gc, z_gc);

  double filling_factor=0.07*exp(2.*abs(z_gc/CGS_U_kpc));
  if(abs(z_gc)>0.75*CGS_U_kpc){filling_factor=0.32;}
  return filling_factor;
}
