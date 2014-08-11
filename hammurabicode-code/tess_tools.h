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

#ifndef TESS_TOOLS
#define TESS_TOOLS

#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include "arr.h"
#include "vec3.h"
#include "CGS_units_file.h"
#include "pointing.h"
#include "healpix_base.h"



void Log ( std::string message, bool flushit = true);
void ErrLog( std::string message , int errval = -1);
void Cyl2Cart( double phi, double invec[3], double outvec[3]);
long Index3d(long n1, long n2, long n3, long i, long j, long l);
//template <typename T> T Mean ( T inarray[], int size);
//template <typename T> T Variance ( T inarray[], int size);

double Mean ( double inarray[], int size);
double Variance ( double inarray[], int size);
double Variance ( double inarray[], int size, double &avg );
void LOS2Cart(double lat, double lon, vec3 infield, vec3 &outfield);
void CoordsCart2LOS(double lat, double lon, vec3 infield, vec3 &outfield);
void Cart2LOS(double lat, double lon, vec3 infield, vec3 &outfield);
bool FileExists( const std::string filename);
void polint(double inx[], double iny[], int n, double x, double *y, double *dy);
double Min ( double *inarray, int n, int &p );
int get_lonlat_region( int nside, double lon_min, double lon_max, double lat_min, double lat_max, arr<long> &pixlist);
int count_lonlat_region( int nside, double lon_min, double lon_max, double lat_min, double lat_max);
double get_lat(int nside, int pixel);
double get_lon(int nside, int pixel);

#endif
