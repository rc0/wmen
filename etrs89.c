/*
* Copyright (c) 2012, Richard P. Curnow
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
*     * Redistributions of source code must retain the above copyright
*       notice, this list of conditions and the following disclaimer.
*     * Redistributions in binary form must reproduce the above copyright
*       notice, this list of conditions and the following disclaimer in the
*       documentation and/or other materials provided with the distribution.
*     * Neither the name of the <organization> nor the
*       names of its contributors may be used to endorse or promote products
*       derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
* ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* */

/* Produce contour plots of E,N deltas between ITRS and ETRS89
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "contour.h"
#include "osxx02.h"

static int get_itrs(double lat, double lon, double *Ei, double *Ni)
{
  struct llh itrs;
  struct en en_itrs;
  struct en en_osgb;

  itrs.lat = lat;
  itrs.lon = lon;
  itrs.h = 0;
  wgs84_to_grid(&itrs, &en_itrs);
  if (apply_ostn02(&en_itrs, &en_osgb) == 0) return 0;
  *Ei = en_osgb.E;
  *Ni = en_osgb.N;
  return 1;
}

static int get_etrs89(double lat, double lon, int year, double *Ee89, double *Ne89)
{
  struct llh itrs;
  struct llh etrs89;
  struct en en_etrs;
  struct en en_osgb;

  itrs.lat = lat;
  itrs.lon = lon;
  itrs.h = 0;
  itrs2005_to_etrs89(&itrs, year, &etrs89);
#if 0
  printf("%.10f %.10f -> %.10f %.10f\n",
      itrs.lat, itrs.lon,
      etrs89.lat, etrs89.lon);
  fflush(stdout);
#endif
  wgs84_to_grid(&etrs89, &en_etrs);
  if (apply_ostn02(&en_etrs, &en_osgb) == 0) return 0;
  *Ee89 = en_osgb.E;
  *Ne89 = en_osgb.N;
  return 1;
}

/* Convert the coordinates to web mercator */
static void remap_clines(const struct cline *lines, double lonmin, double lonstep, double latmin, double latstep)/*{{{*/
{
  const struct cline *l;
  for (l=lines; l; l=l->next) {
    struct cpoint *p;
    for (p=l->points; p; p=p->next) {
      struct llh llh;
      struct mxy mxy;
      llh.lat = latmin + latstep * p->y;
      llh.lon = lonmin + lonstep * p->x;
      llh.h = 0;
      wgs84_to_mxy(&llh, &mxy);
      p->x = mxy.X;
      p->y = mxy.Y;
    }
  }
}
/*}}}*/

#define C0 "black", "#000000"
#define C7 "red", "#ff0000"
#define THICK1 "thin", 1.0
#define THICK2 "thick", 1.4


#define N_EAST 12

struct level east_levels[N_EAST] = 
{
  {C7, -0.23, THICK2},
  {C7, -0.24, THICK2},
  {C7, -0.25, THICK2},
  {C7, -0.26, THICK2},
  {C7, -0.27, THICK2},
  {C7, -0.28, THICK2},
  {C7, -0.29, THICK2},
  {C7, -0.30, THICK2},
  {C7, -0.31, THICK2},
  {C7, -0.32, THICK2},
  {C7, -0.33, THICK2},
  {C7, -0.34, THICK2}
};

#define N_NORTH 6
struct level north_levels[N_NORTH] = 
{
  {C7, -0.40, THICK2},
  {C7, -0.41, THICK2},
  {C7, -0.42, THICK2},
  {C7, -0.43, THICK2},
  {C7, -0.44, THICK2},
  {C7, -0.45, THICK2},
};

int main (int argc, char **argv)
{

  double lonmin, lonmax;
  double latmin, latmax;
  double step;
  double lat, lon;
  struct cdata *outer_e, *outer_n;
  struct crow *inner_e, *inner_n;
  struct cline *lines;
  int i;

  load_osxx02();

  step = 0.025;

  latmin = 54.0;
  latmax = 56.0;
  lonmin = -1.0;
  lonmax = +1.0;
#if 1
  latmin = 49.0;
  latmax = 61.0;
  lonmin = -8.0;
  lonmax = +2.0;
#endif

  outer_e = new_cdata();
  outer_n = new_cdata();
  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    inner_e = next_crow(outer_e);
    inner_n = next_crow(outer_n);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      double Ei, Ni, Ee89, Ne89;
      double de, dn;
      int si, se;
      si = get_itrs(lat, lon, &Ei, &Ni);
      se = get_etrs89(lat, lon, 2012, &Ee89, &Ne89);
      if (si && se) {
        de = Ee89 - Ei;
        dn = Ne89 - Ni;
        printf("%9.2f,%9.2f  %9.2f,%9.2f  %6.2f,%6.2f\n",
            Ei, Ni,
            Ee89, Ne89,
            de, dn);
        fflush(stdout);
        add_cnode(inner_e, de);
        add_cnode(inner_n, dn);
      } else {
        add_empty_cnode(inner_e);
        add_empty_cnode(inner_n);
      }
    }
  }

  start_svg("itrs_etrs_east.svg");
  for (i=0; i<N_EAST; i++) {
    lines = generate_isolines(outer_e, east_levels[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    emit_svg(east_levels[i].level, east_levels[i].svg_colour, east_levels[i].width_scale, lines);
#if 0
    emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
#endif
    free_clines(lines);
  }
  finish_svg();

  start_svg("itrs_etrs_north.svg");
  for (i=0; i<N_NORTH; i++) {
    lines = generate_isolines(outer_n, north_levels[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    emit_svg(north_levels[i].level, north_levels[i].svg_colour, north_levels[i].width_scale, lines);
#if 0
    emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
#endif
    free_clines(lines);
  }
  finish_svg();

  free_cdata(outer_e);
  free_cdata(outer_n);

  return 0;
}
