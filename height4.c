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

/* Find an approximation to map lat/long to the difference between WGS84
 * elevation and height above mean sea level in UK (ODN) */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "contour.h"
#include "osxx02.h"

#if 0
/* Accurate to <0.5m across whole country */
#define OP 3
#define OQ 2
#else
/* Adds about a 1cm error to the Helmert transform formula in the OS document,
 * but we cull 3 terms to get a 6 term formula so resulting additional error is
 * less than about 10cm. */
#define OP 2
#define OQ 2
#endif
#define OPQ (OP*OQ)

struct model
{
  double lat0;
  double lon0;
  double slat;
  double slon;
  double c[OP][OQ];
};

static void fit(struct model *model)
{
  /* p = lat, q = lon */
#if 1
  double lat0 = 54.5;
  double lon0 = -2.0;
#elif 0
  double lat0 = 52.0;
  double lon0 = -2.0;
#else
  double lat0 = 49.0;
  double lon0 = -6.0;
#endif
  double latmin = 50.0;
  double latmax = 61.0;
  double lonmin = -8.0;
  double lonmax = 2.0;
  double slon = 1.0/4.0;
  double slat = 1.0/4.5;
  double step = 0.01;

  double lat, lon;
  double p, q;
  long double pp[2*OP-1], qq[2*OQ-1];

  Matrix M;
  long double L[OPQ], R[OPQ];
  int i, j, m, n;

  for (i=0; i<OPQ; i++) {
    for (j=0; j<OPQ; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<OPQ; i++) {
    R[i] = 0.0;
  }

  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    uk_range(lat, &lonmin, &lonmax);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      struct llh wgs;
      struct en en_wgs;
      double corr;
      wgs.lat = lat;
      wgs.lon = lon;
      wgs.h = 0;
      wgs84_to_grid(&wgs, &en_wgs);
      if (apply_osgm02(&en_wgs, 1, &corr) == 0) continue;
      corr = -corr;
      p = slat * (lat - lat0);
      q = slon * (lon - lon0);
      pp[0] = 1.0;
      qq[0] = 1.0;
      for (i=1; i<2*OP-1; i++) {
        pp[i] = p * pp[i-1];
      }
      for (i=1; i<2*OQ-1; i++) {
        qq[i] = q * qq[i-1];
      }
      for (i=0; i<OQ; i++) {
        for (j=0; j<OQ; j++) {
          for (m=0; m<OP; m++) {
            for (n=0; n<OP; n++) {
              M[OP*i+m][OP*j+n] += pp[m+n]*qq[i+j];
            }
          }
        }
      }
      for (i=0; i<OP; i++) {
        for (j=0; j<OQ; j++) {
          R[i+OP*j] += corr * pp[i] * qq[j];
        }
      }
    }
  }

  solve(M, OPQ, L, R);

  model->lat0 = lat0;
  model->lon0 = lon0;
  model->slat = slat;
  model->slon = slon;
  for (i=0; i<OP; i++) {
    for (j=0; j<OQ; j++) {
      model->c[i][j] = L[i+OP*j];
    }
  }
  return;
}

static double predict(double lat, double lon, const struct model *model)
{
  double dlat, dlon;
  double pp[OP], qq[OQ];
  int i, j;
  double result;
  pp[0] = 1.0;
  qq[0] = 1.0;
  dlat = model->slat * (lat - model->lat0);
  dlon = model->slon * (lon - model->lon0);
  for (i=1; i<OP; i++) {
    pp[i] = dlat * pp[i-1];
  }
  for (i=1; i<OQ; i++) {
    qq[i] = dlon * qq[i-1];
  }
  result = 0.0;
  for (i=0; i<OP; i++) {
    for (j=0; j<OQ; j++) {
      result += model->c[i][j] * pp[i] * qq[j];
    }
  }
  return result;
}

static void puncture_model(struct model *model)
{
  int i, j;
  for (i=0; i<OP; i++) {
    for (j=0; j<OQ; j++) {
      if (fabs(model->c[i][j]) < 1.00) {
        model->c[i][j] = 0.0;
      }
    }
  }
  return;
}

static void round_model(struct model *model)
{
  int i, j;
  for (i=0; i<OP; i++) {
    for (j=0; j<OQ; j++) {
      model->c[i][j] = 0.1 * round(10.0 * model->c[i][j]);
    }
  }
}

static void print_model(const struct model *model)
{
  int i, j;
  printf("lat0=%f\n", model->lat0);
  printf("lon0=%f\n", model->lon0);
  printf("slat=%f\n", model->slat);
  printf("slon=%f\n", model->slon);
  for (i=0; i<OP; i++) {
    for (j=0; j<OQ; j++) {
      printf("c[%d][%d] = %f\n", i, j, model->c[i][j]);
    }
  }
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

static double predict2(double lat, double lon)
{
  double M = lat;
  double L = lon;

  return -24.64 - ((L+4.93)*(M-72.8))/12.0;
}

int main(int argc, char **argv)
{
  double latmin;
  double latmax;
  double lonmin;
  double lonmax;
  double lat, lon;
  double step;
  struct model model;
  struct cdata *cdata;
  struct crow *row;
  struct cline *lines;
  int i;
  int is_svg;

  is_svg = (argc > 1);

  load_osxx02();

  fit(&model);
  print_model(&model);
  puncture_model(&model);
  round_model(&model);
  print_model(&model);

  latmin = 49.0;
  latmax = 61.0;
  lonmin = -8.0;
  lonmax = +2.0;
  cdata = new_cdata();
  step = 0.025;
  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    row = next_crow(cdata);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      double actual;
      double predicted;
      double predicted2;
      double error;
      struct llh wgs;
      struct en en_wgs;
      predicted = predict(lat, lon, &model);
      predicted2 = predict2(lat, lon);
      wgs.lat = lat;
      wgs.lon = lon;
      wgs.h = 0;
      wgs84_to_grid(&wgs, &en_wgs);
      if (apply_osgm02(&en_wgs, 0, &actual)) {
        actual = -actual;
        error = actual - predicted;
        add_cnode(row, error);
        printf("%8.3f %8.3f %8.2f %8.2f %8.2f %8.2f\n", lat, lon, actual, predicted, error, predicted2);
      } else {
        add_empty_cnode(row);
      }
    }
  }

  if (is_svg) {
    start_svg("height4_new.svg");
  } else {
    start_tikz("height4.tex");
  }
  for (i=0; i<N_LEVELS_2; i++) {
    double level;
    level = levels_2side[i].level;
    if ((fabs(level) > 0.0)  && (fabs(level) < 0.5)) continue;
    lines = generate_isolines(cdata, level);
    remap_clines(lines, lonmin, step, latmin, step);
    if (is_svg) {
      emit_svg(level, levels_2side[i].svg_colour, levels_2side[i].width_scale, lines);
    } else {
      emit_tikz(level, levels_2side[i].colour, levels_2side[i].thickness, lines);
    }
    free_clines(lines);
  }
  if (is_svg) {
    finish_svg();
  } else {
    finish_tikz();
  }

  generate_contour_key_tikz(stdout);

  return 0;
}

