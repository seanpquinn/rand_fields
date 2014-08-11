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

#include "tess_tools.h"

void Log ( std::string message, bool flushit) {
  char stime[81];
  time_t system_time = time(NULL);
  char time_text[81];
  strftime(time_text, 80, "%Y-%m-%dT%R:%S", localtime(&system_time));
  sprintf(stime,"[%s]  ", time_text);
  message = stime+message;
  std::cout << message;
  if (flushit) std::cout<<std::flush;
}

void ErrLog( std::string message , int errval) {
  Log(message);
  std::exit(errval);
}


void Cyl2Cart( double phi, double invec[3], double outvec[3]){
        outvec[0]=0.;outvec[1]=0.;outvec[2]=0.;
        arr2<double> cyl_unit_vec_array(3,3);
        cyl_unit_vec_array[0][0]=cos(phi);
        cyl_unit_vec_array[1][0]=sin(phi);
        cyl_unit_vec_array[2][0]=0.;
        cyl_unit_vec_array[0][1]=-sin(phi);
        cyl_unit_vec_array[1][1]=cos(phi);
        cyl_unit_vec_array[2][1]=0.;
        cyl_unit_vec_array[0][2]=0.;
        cyl_unit_vec_array[1][2]=0.;
        cyl_unit_vec_array[2][2]=1.;

        for (int n=0;n<3;n++){
                for (int m=0;m<3;m++) {
                        outvec[n]=outvec[n]+cyl_unit_vec_array[n][m]*invec[m];
                }
        }
}// Cyl2Cart

long Index3d(long /* n1 */, long n2, long n3, long i, long j, long l)
{ return (i*n2*n3 + j*n3 + l);}


/*
template <typename T> T Mean ( T inarray[], int size) {
  T avg=0;
  for (int m=0; m<size; ++m) avg += inarray[m];
  avg/=size;
  return avg;
}

template <typename T> T Variance ( T inarray[], int size)
{
  T avg=Mean(inarray,size);
  T var = 0.;
  for (int m=0; m<size; ++m) var += pow(inarray[m]-avg,2);
  var/=size;
  return var;
}
*/


double Mean ( double inarray[], int size) {
  double avg=0;
  for (int m=0; m<size; ++m) avg += inarray[m];
  avg/=size;
  return avg;
}

double Variance ( double inarray[], int size)
{
  double avg=Mean(inarray,size);
  double var = 0.;
  for (int m=0; m<size; ++m) var += (inarray[m]-avg)*(inarray[m]-avg);
  var/=size;
  return var;
}

// If you want the average back, set it to zero first, else its
// current value will be used instead of calculated!
double Variance ( double inarray[], int size, double &avg)
{
  if (fabs(avg) < 1e-20) avg=Mean(inarray,size);
  double var = 0.;
  for (int m=0; m<size; ++m) var +=(inarray[m]-avg)*(inarray[m]-avg);
  var/=size;
  return var;
}

//Tess' original Cart2LOS - MODIFIED!
void Cart2LOS(double lat, double lon, vec3 infield, vec3 &outfield){

        /* Tess' original transform
        outfield.x=     infield.x*cos(lat)*cos(lon)             + infield.y*cos(lat)*sin(lon) + infield.z*sin(lat);
        outfield.y=     -1.*infield.x*cos(lat)*sin(lon) + infield.y*cos(lat)*cos(lon);
        outfield.z=     infield.x*sin(lat)*cos(lon)             + infield.y*sin(lat)*sin(lon) - infield.z*cos(lat);
         */
        
        ///* My modified transform
        outfield.x=     infield.x*cos(lat)*cos(lon)             + infield.y*cos(lat)*sin(lon) + infield.z*sin(lat);
        outfield.y=     -1.*infield.x*sin(lon)                  + infield.y*cos(lon);
        outfield.z=     -1.*infield.x*sin(lat)*cos(lon) - infield.y*sin(lat)*sin(lon) + infield.z*cos(lat);
        // */
}//Cart2LOS

//Transpose of Tess' modified Cart2LOS
void LOS2Cart(double lat, double lon, vec3 infield, vec3 &outfield){
        
        /* Tess' original transform
        outfield.x= infield.x*cos(lat)*cos(lon)         - infield.y*cos(lat)*sin(lon) + infield.z*sin(lat)*cos(lon); 
        outfield.y= infield.x*cos(lat)*sin(lon)         + infield.y*cos(lat)*cos(lon) + infield.z*sin(lat)*sin(lon);
        outfield.z= infield.x*sin(lat)                                                                                    - infield.z*cos(lat);
        */
        
        ///* My modified transform
        outfield.x=     infield.x*cos(lat)*cos(lon)             - infield.y*sin(lon)    - infield.z*sin(lat)*cos(lon); 
        outfield.y=     infield.x*cos(lat)*sin(lon)             + infield.y*cos(lon)    - infield.z*sin(lat)*sin(lon);
        outfield.z=     infield.x*sin(lat)                                                                              + infield.z*cos(lat);
        //*/
}



bool FileExists( const std::string filename) {
  std::fstream fin;
  fin.open(filename.data(), std::ifstream::in);
  if( fin.is_open() ) {
    fin.close();
    return true;
  }
  else return false;
}



// My polint that is a copy of the numerical recipes version but with
// the inputs normal C arrays instead of their f'ed up "unit-offset"
// arrays.  They recommend using the original code and calling it with
// array-1 instead of array, but I think that's dangerous.
void polint(double inx[], double iny[], int n, double x, double *y, double *dy)
{
  // Point xa and ya one element before input arrays so that you can
  // then use xa[1] to xa[n] etc.  First value set to dummy.
        double *xa=inx-1;
        double *ya=iny-1;

        int i,m,ns=1;
        double den,dif,dift,ho,hp,w;
        //      double *c,*d;
        double cc[n],dd[n];
        double *c=cc-1;
        double *d=dd-1;

        dif=fabs(x-xa[1]);
        //      c=vector(1,n);
        //      d=vector(1,n);
        for (i=1;i<=n;i++) {
                if ( (dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        *y=ya[ns--];
        for (m=1;m<n;m++) {
                for (i=1;i<=n-m;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+m]-x;
                        w=c[i+1]-d[i];
                        if ( (den=ho-hp) == 0.0)
                                ErrLog("Error in routine polint");
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
        }

} // polint


//////////////////////////////////////////////////////////////////////
//  Min function
//////////////////////////////////////////////////////////////////////
double Min ( double *inarray, int n, int &p )
{

  double min = inarray[0];
  for (long i=0; i<n;++i) {
    if (min > inarray[i] )  { min = inarray[i]; p = i; }
  }
  return min;

}

int count_lonlat_region( int nside, double minLat, double maxLat, double minLon, double maxLon)
{
	int npix;
	double THE, PHI;
	double maxTheta, maxPhi, minTheta, minPhi;
	int count=0;

	double ratio = 180.0 / CGS_U_pi;
	
	minPhi=minLon/ratio;
	maxPhi=maxLon/ratio;
	minTheta=(CGS_U_pi/2.)-(maxLat/ratio);
	maxTheta=(CGS_U_pi/2.)-(minLat/ratio);
	pointing ptg_obj;
	Healpix_Base map(nside, NEST, SET_NSIDE);
	npix=map.Npix();

	
	for (int i=0; i<npix; i++) {
		
		ptg_obj=map.pix2ang(i); // Will have the angles theta and phi in the ptg object.		
		
		THE=ptg_obj.theta;
		PHI=ptg_obj.phi;
		
		if (THE>minTheta && THE<maxTheta && PHI>minPhi && PHI<maxPhi) {

			count++;
		}
		
	}
	return count;
	
}


int get_lonlat_region( int nside, double minLat, double maxLat, double minLon, double maxLon, arr<long> &pixlist)
{
	int npix;
	double THE, PHI;
	double maxTheta, maxPhi, minTheta, minPhi;
	int count=0;
	double ratio = 180.0 / CGS_U_pi;

	minPhi=minLon/ratio;
	maxPhi=maxLon/ratio;
	minTheta=(CGS_U_pi/2.)-(maxLat/ratio);
	maxTheta=(CGS_U_pi/2.)-(minLat/ratio);
	pointing ptg_obj;
	Healpix_Base map(nside, NEST, SET_NSIDE);
	npix=map.Npix();
	
	
	for (int i=0; i<npix; i++) {
		
		ptg_obj=map.pix2ang(i); // Will have the angles theta and phi in the ptg object.		
		
		THE=ptg_obj.theta;
		PHI=ptg_obj.phi;
		
		if (THE>minTheta && THE<maxTheta && PHI>minPhi && PHI<maxPhi) {
			
			pixlist[count] = i;
			count++;
		}
		
	}
	return count;
	
}

double get_lat( int nside, int pixel)
{

	double THE;
	double lat;

	double ratio = 180.0 / CGS_U_pi;

	pointing ptg_obj;
	Healpix_Base map(nside, NEST, SET_NSIDE);
	ptg_obj=map.pix2ang(pixel); // Will have the angles theta and phi in the ptg object.		
		
	THE=ptg_obj.theta;

	lat=((CGS_U_pi/2.)-THE)*ratio;
	
	return lat;
	
}

double get_lon( int nside, int pixel)
{

	double PHI;
	double lon;
	
	double ratio = 180.0 / CGS_U_pi;
	
	pointing ptg_obj;
	Healpix_Base map(nside, NEST, SET_NSIDE);
	ptg_obj=map.pix2ang(pixel); // Will have the angles theta and phi in the ptg object.		
	
	PHI=ptg_obj.phi;

	lon=PHI*ratio;
	
	return lon;
	
}

