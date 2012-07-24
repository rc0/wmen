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

/* Generate conversions between ITRS2005 and ETRS89 -based Eastings, Northings,
 * as a function of WM.  These corrections could then be rolled into the
 * formulae out of wm_to_grid, if desired.
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "contour.h"
#include "osxx02.h"

static int get_itrs(double lat, double lon, double *Ei, double *Ni)/*{{{*/
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
/*}}}*/
static int get_etrs89(double lat, double lon, double year, double *Ee89, double *Ne89)/*{{{*/
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
/*}}}*/
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
#define N_NORTH 6

struct level east_levels[N_EAST] = /*{{{*/
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
/*}}}*/
struct level north_levels[N_NORTH] = /*{{{*/
{
  {C7, -0.40, THICK2},
  {C7, -0.41, THICK2},
  {C7, -0.42, THICK2},
  {C7, -0.43, THICK2},
  {C7, -0.44, THICK2},
  {C7, -0.45, THICK2},
};
/*}}}*/
static void contours_2012(void)/*{{{*/
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
#if 0
        printf("%9.2f,%9.2f  %9.2f,%9.2f  %6.2f,%6.2f\n",
            Ei, Ni,
            Ee89, Ne89,
            de, dn);
#endif
        fflush(stdout);
        add_cnode(inner_e, de);
        add_cnode(inner_n, dn);
      } else {
        add_empty_cnode(inner_e);
        add_empty_cnode(inner_n);
      }
    }
  }

  start_svg("itrs_etrs_2012_east.svg");
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

  start_svg("itrs_etrs_2012_north.svg");
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
}
/*}}}*/

#define OX 2
#define OY 2
#define OT 2
#define OXY (OX*OY)
#define OXYT (OX*OY*OT)

struct model/*{{{*/
{
  double X0, Y0, T0;
  double SX, SY, ST;
  /* TODO : not even clear that the same order polynomials are appropriate for
   * both of these? */
  double e[OX][OY][OT];
  double n[OX][OY][OT];
};

/*}}}*/
static void contours_err(const struct model *m)/*{{{*/
{
  double lonmin, lonmax;
  double latmin, latmax;
  double step;
  double lat, lon;
  struct llh wgs;
  struct mxy mxy;
  struct cdata *outer_e, *outer_n;
  struct crow *inner_e, *inner_n;
  struct cline *lines;
  int i;
  double t = 2012.0;

  load_osxx02();

  step = 0.025;

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
      double Ep, Np;
      double de, dn;
      double dx, dy, dt;
      int si, se;
      si = get_itrs(lat, lon, &Ei, &Ni);
      se = get_etrs89(lat, lon, 2012, &Ee89, &Ne89);
      if (si && se) {
        wgs.lat = lat;
        wgs.lon = lon;
        wgs.h = 0;
        wgs84_to_mxy(&wgs, &mxy);
        dx = m->SX * (mxy.X - m->X0);
        dy = m->SY * (mxy.Y - m->Y0);
        dt = m->ST * (t - m->T0);
        Ep =  0.04989 - 0.01455*dt - 0.00190*dt*dy;
        Np = -0.06564 - 0.01602*dt - 0.00112*dt*dx;

        de = 1000.0 * ((Ee89 - Ei) - Ep);
        dn = 1000.0 * ((Ne89 - Ni) - Np);
        fflush(stdout);
        add_cnode(inner_e, de);
        add_cnode(inner_n, dn);
      } else {
        add_empty_cnode(inner_e);
        add_empty_cnode(inner_n);
      }
    }
  }

  start_svg("itrs_etrs_error_2012_east.svg");
  for (i=0; i<N_LEVELS_2; i++) {
    lines = generate_isolines(outer_e, levels_2side[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    emit_svg(levels_2side[i].level, levels_2side[i].svg_colour, levels_2side[i].width_scale, lines);
#if 0
    emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
#endif
    free_clines(lines);
  }
  finish_svg();

  start_svg("itrs_etrs_error_2012_north.svg");
  for (i=0; i<N_LEVELS_2; i++) {
    lines = generate_isolines(outer_n, levels_2side[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    emit_svg(levels_2side[i].level, levels_2side[i].svg_colour, levels_2side[i].width_scale, lines);
#if 0
    emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
#endif
    free_clines(lines);
  }
  finish_svg();

  free_cdata(outer_e);
  free_cdata(outer_n);
}
/*}}}*/
static void fit_wm(struct model *model)/*{{{*/
{
  double latmin = 49.8;
  double latmax = 61.0;
  double lonmin, lonmax;
  double tmin, tmax, tstep;
  double step = 0.1;

  double lat, lon;
  struct llh wgs;
  struct mxy mxy;
  struct en en;
  double dx, dy, dt;
  long double xx[2*OX-1], yy[2*OY-1], tt[2*OT-1];
  long double LE[OXYT], LN[OXYT], RE[OXYT], RN[OXYT];
  Matrix M;
  int i, j, k, m, n, u, v;
  double t;

  tmin = 1989.0;
  tmax = 2039.00001;
  tstep = 0.05;

  for (i=0; i<OXYT; i++) {
    for (j=0; j<OXYT; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<OXYT; i++) {
    RE[i] = 0.0;
    RN[i] = 0.0;
  }

#if 1
  model->X0 = 0.494440093;
  model->Y0 = 0.312663855; /* 0.3187 */;
  model->T0 = 1989.0;
  model->SX = 61.0;
  model->SY = 36.0;
  model->ST = 1.0;
#endif

  for (t = tmin; t < tmax; t += tstep) {
    for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
      uk_range(lat, &lonmin, &lonmax);
      for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {

        double Ei, Ni, Ee89, Ne89;
        double de, dn;
        int si, se;

        wgs.lat = lat;
        wgs.lon = lon;
        wgs.h = 0;
        wgs84_to_mxy(&wgs, &mxy);
        si = get_itrs(lat, lon, &Ei, &Ni);
        se = get_etrs89(lat, lon, t, &Ee89, &Ne89);

        if (si && se) {
          de = Ee89 - Ei;
          dn = Ne89 - Ni;
          dx = model->SX * (mxy.X - model->X0);
          dy = model->SY * (mxy.Y - model->Y0);
          dt = model->ST * (t - model->T0);
#if 0
          printf("%f %f %f %f\n", lat, lon, dx, dy);
          fflush(stdout);
#endif
          xx[0] = 1.0;
          for (i=1; i<2*OX-1; i++) {
            xx[i] = dx * xx[i-1];
          }
          yy[0] = 1.0;
          for (i=1; i<2*OY-1; i++) {
            yy[i] = dy * yy[i-1];
          }
          tt[0] = 1.0;
          for (i=1; i<2*OT-1; i++) {
            tt[i] = dt * tt[i-1];
          }

          for (i=0; i<OY; i++) {
            for (j=0; j<OY; j++) {
              for (m=0; m<OX; m++) {
                for (n=0; n<OX; n++) {
                  for (u=0; u<OT; u++) {
                    for (v=0; v<OT; v++) {
                      M[OXY*u+OX*i+m][OXY*v+OX*j+n] += xx[m+n]*yy[i+j]*tt[u+v];
                    }
                  }
                }
              }
            }
          }
          for (i=0; i<OX; i++) {
            for (j=0; j<OY; j++) {
              for (k=0; k<OT; k++) {
                RE[k*OXY+i+OX*j] += de * xx[i] * yy[j] * tt[k];
                RN[k*OXY+i+OX*j] += dn * xx[i] * yy[j] * tt[k];
              }
            }
          }
        }
      }
    }
  }

  solve2(M, OXYT, LE, LN, RE, RN);

  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      for (k=0; k<OT; k++) {
        model->e[i][j][k] = LE[k*OXY+i+OX*j];
        model->n[i][j][k] = LN[k*OXY+i+OX*j];
      }
    }
  }
}
/*}}}*/
static void print_model(const struct model *model) {/*{{{*/
  int i, j, k;
  printf("X0=%.10f\n", model->X0);
  printf("Y0=%.10f\n", model->Y0);
  printf("T0=%.10f\n", model->T0);
  printf("SX=%.3f\n", model->SX);
  printf("SY=%.3f\n", model->SY);
  printf("ST=%.3f\n", model->ST);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      for (k=0; k<OT; k++) {
        printf("e[%d][%d][%d] = %15.5f   n[%d][%d][%d] = %15.5f\n",
            i, j, k,
            model->e[i][j][k],
            i, j, k,
            model->n[i][j][k]
            );
      }
    }
  }
}
/*}}}*/

int main (int argc, char **argv)/*{{{*/
{
  struct model model;
  struct model model2;
  load_osxx02();
#if 0
  fit_wm(&model);
  print_model(&model);
#endif

#if 1
  contours_2012();
#endif

  model2.X0 = 0.494440093;
  model2.Y0 = 0.312663855; /* 0.3187 */;
  model2.T0 = 1989.0;
  model2.SX = 61.0;
  model2.SY = 36.0;
  model2.ST = 1.0;
  memset(model2.e, 0, sizeof(model2.e));
  memset(model2.n, 0, sizeof(model2.n));
  contours_err(&model2);

  return 0;
}
/*}}}*/
