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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "contour.h"

static FILE *out;
/* pngs are loaded into inkscape at 1 pixel per point and each tile is 256
 * pixels in the png */

static inline double radof(double deg)/*{{{*/
{
  return deg * (M_PI/180.0);
}
/*}}}*/
static void mxy_to_en(double x, double y, double *E, double *N)/*{{{*/
{
  x = 72.0 * (x - 0.494444);
  y = 36.0 * (y - 0.318700);

  double y2 = y*y;
  double y3 = y*y2;
  double y4 = y*y3;
  double y5 = y*y4;
  double y6 = y*y5;

  double x2 = x*x;
  double x3 = x*x2;
  double x4 = x*x3;
  double x5 = x*x4;

  double E0 = +400087.5 -3.9*y;
  double E1 = +323810.0 +45906.2*y +1579.8*y2 -241.1*y3 -37.6*y4 -1.9*y5;
  double E3 = -133.5 +60.0*y +18.9*y2 +1.7*y3;
  double E5 = -0.5;

  *E = E0 +E1*x +E3*x3 +E5*x5;

  double N0 = +511732.7 -646151.6*y -45594.0*y2 -1016.5*y3 +122.3*y4 +14.8*y5 +0.6*y6;
  double N1 = -1.9;
  double N2 = +11502.6 +799.1*y -180.0*y2 -37.7*y3 -2.5*y4;
  double N4 = +7.5 +4.7*y +0.6*y2;

  *N = N0 +N1*x +N2*x2 +N4*x4;
}
/*}}}*/
void en_to_wm(double E, double N, double *X, double *Y)/*{{{*/
{
  double e = (E - 400000.0) / 400000.0;
  double n = (N - 650000.0) / 650000.0;
  double e2 = e*e;
  double e3 = e2*e;
  double e4 = e3*e;
  double e5 = e4*e;
  double n2 = n*n;
  double n3 = n2*n;
  double n4 = n3*n;
  double n5 = n4*n;
  double n6 = n5*n;

  double X0 = +12361001.62 -19.66*n -2.24*n2;
  double X1 = +442441.04 +66066.30*n +12171.89*n2 +2057.38*n3 +364.45*n4 +62.11*n5;
  double X2 = -2.52 -23.19*n -23.76*n2 +3.72*n3;
  double X3 = -1540.23 -778.77*n -273.35*n2 -63.49*n3 -25.26*n4 -12.21*n5;
  double X4 = +7.76*n3 +6.93*n4;
  double X5 = +10.01 +12.09*n4 +10.61*n5;
  *X = (X0 +X1*e +X2*e2 +X3*e3 +X4*e4 +X5*e5) / 25000000.0;
  double Y0 = +7816596.43 -720504.09*n -53575.89*n2 -6596.44*n3 -830.91*n4 -117.01*n5 -15.94*n6;
  double Y1 = -13.88 -5.83*n +4.10*n2 +4.32*n3;
  double Y2 = +20368.69 +7485.47*n +1909.20*n2 +454.99*n3 +98.62*n4 +23.58*n5;
  double Y3 = -9.34 -12.10*n +26.34*n2 +23.54*n3 -5.75*n4;
  double Y4 = -119.10 -84.19*n -41.68*n2 -7.53*n3 +20.95*n4 +15.10*n5;
  double Y5 = +14.47*n4 +30.69*n5 +14.37*n6;
  *Y = (Y0 +Y1*e +Y2*e2 +Y3*e3 +Y4*e4 +Y5*e5) / 25000000.0;
}
/*}}}*/

#define X0 30.0
#define Y0 22.0
#define SCALE 3.5

void XY_to_xy(double X, double Y, double *x, double *y)/*{{{*/
{
  *x = SCALE * (X*REF6 - X0);
  *y = SCALE * (Y0 - Y*REF6);
}
/*}}}*/
#define MAXP 32768
static void inner_reduce(double *e, double *n, unsigned char *k, int i0, int i1, double thresh)/*{{{*/
{
  /* implement Douglas-Peucker algorithm to simplify line shapes */
  double de, dn, d;
  double u, v;

  if ((i1 - i0) < 2) return;

  de = e[i1] - e[i0];
  dn = n[i1] - n[i0];
  d = sqrt(de*de + dn*dn);
  if (d > 0.0) {
    int i;
    int maxi = -1;
    double maxd;
    u = -dn/d;
    v =  de/d;
    for (i=i0+1; i<i1; i++) {
      de = e[i] - e[i0];
      dn = n[i] - n[i0];
      d = fabs(de*u + dn*v);
      if ((maxi < 0) || (d > maxd)) {
        maxi = i;
        maxd = d;
      }
    }
    if (maxd > thresh) {
      k[maxi] = 1;
      inner_reduce(e, n, k, i0, maxi, thresh);
      inner_reduce(e, n, k, maxi, i1, thresh);
    } else {
      /* don't keep any */
    }

  } else {
    /* degenerate path : find the point with the greatest distance from the start */
    int i;
    int maxi = -1;
    double maxd;
    for (i=i0+1; i<i1; i++) {
      de = e[i] - e[i0];
      dn = n[i] - n[i0];
      d = sqrt(de*de + dn*dn);
      if ((maxi < 0) || (d > maxd)) {
        maxi = i;
        maxd = d;
      }
    }
    k[maxi] = 1;
    inner_reduce(e, n, k, i0, maxi, thresh);
    inner_reduce(e, n, k, maxi, i1, thresh);
  }
}
/*}}}*/
static int reduce(int np, double *e, double *n, unsigned char *k, double thresh)/*{{{*/
{
  int is_cycle;
  int count, i;

  double de, dn, d;

  de = e[np-1] - e[0];
  dn = n[np-1] - n[0];
  d = sqrt(de*de + dn*dn);
  is_cycle = (d < thresh) ? 1 : 0;

  memset(k, 0, np);
  k[0] = k[np-1] = 1;
  inner_reduce(e, n, k, 0, np-1, thresh);

  count = 0;
  for (i=0; i<np; i++) {
    if (k[i]) ++count;
  }

  return is_cycle;
}
/*}}}*/


struct cline *add_cline(struct cline *last)/*{{{*/
{
  struct cline *result;
  result = malloc(sizeof(struct cline));
  result->next = last;
  result->points = NULL;
  return result;
}
/*}}}*/
struct cpoint *new_point(double x, double y)/*{{{*/
{
  struct cpoint *result;
  result = malloc(sizeof(struct cpoint));
  result->x = x;
  result->y = y;
  result->next = NULL;
  return result;
}
/*}}}*/
static void append_point(struct cline *cl, struct cpoint *p)/*{{{*/
{
  p->next = cl->points;
  cl->points = p;
}
/*}}}*/
struct cline *load_map(const char *filename)/*{{{*/
{
  FILE *in;
  char line[256];
  double lat, lon;
  double x, yy, y;
  struct cline *result, *current;

  result = current = NULL;

  in = fopen(filename, "r");
  if (!in) {
    fprintf(stderr, "Cannot open map information at %s\n", filename);
    exit(2);
  }

  while (fgets(line, 256, in)) {
    if (sscanf(line, "%lf%lf", &lon, &lat) == 2) {
      struct cpoint *point;
      x = radof(lon);
      yy = radof(lat);
      y = log(tan(yy) + 1.0/cos(yy));
      x = 0.5 * (1.0 + x/M_PI);
      y = 0.5 * (1.0 - y/M_PI);
      point = new_point(x, y);
      if (!current) {
        current = result = add_cline(result);
      }
      append_point(current, point);
    } else {
      current = NULL;
    }
  }
  return result;
}
/*}}}*/
#define THICKNESS "thin"
void draw_map_tikz(FILE *out, const struct cline *map)/*{{{*/
{
  double E, N;
  const struct cline *l;
  const struct cpoint *p;
  double ee[MAXP], nn[MAXP];
  unsigned char keep[MAXP];
  int q, i;
  int is_cycle;
  int limit;

  for (l=map; l; l=l->next) {
    for (q=0, p=l->points; p; p=p->next, q++) {
      mxy_to_en(p->x, p->y, &E, &N);

      ee[q] = E / 100000.0;
      nn[q] = N / 100000.0;
    }

    is_cycle = reduce(q, ee, nn, keep, 0.01);
    fprintf(out, "\\path [draw,%s,brown!70!black]", THICKNESS);
    fprintf(out, " (%.2fcm,%.2fcm)\n", ee[0], nn[0]);
    limit = is_cycle ? q-1 : q;
    for (i=1; i<limit; i++) {
      if (keep[i]) {
        fprintf(out, " -- (%.2fcm,%.2fcm)\n", ee[i], nn[i]);
      }
    }
    if (is_cycle) {
      fprintf(out, " -- cycle\n");
    }
    fprintf(out, ";\n");
  }
}
/*}}}*/
void draw_map_mxy_tikz(FILE *out, const struct cline *map)/*{{{*/
{
  const struct cline *l;
  const struct cpoint *p;
  double xx[MAXP], yy[MAXP];
  unsigned char keep[MAXP];
  int q, i;
  int is_cycle;
  int limit;
  
  fprintf(stderr, "Got here!\n");

  for (l=map; l; l=l->next) {
    for (q=0, p=l->points; p; p=p->next, q++) {
      XY_to_xy(p->x, p->y, xx+q, yy+q);
    }

    is_cycle = reduce(q, xx, yy, keep, 0.01);
    fprintf(out, "\\path [draw,%s,brown!70!black]", THICKNESS);
    fprintf(out, " (%.2fcm,%.2fcm)\n", xx[0], yy[0]);
    limit = is_cycle ? q-1 : q;
    for (i=1; i<limit; i++) {
      if (keep[i]) {
        fprintf(out, " -- (%.2fcm,%.2fcm)\n", xx[i], yy[i]);
      }
    }
    if (is_cycle) {
      fprintf(out, " -- cycle\n");
    }
    fprintf(out, ";\n");
  }
}
/*}}}*/
void start_tikz(char *name)/*{{{*/
{
  struct cline *map;
  out = fopen(name, "w");
  fprintf(out, "\\begin{tikzpicture}\n");
  fprintf(out, "\\begin{scope}\n");

  /* Draw UK map */
  map = load_map("uk_coast/uk_coast.txt");
  draw_map_tikz(out, map);
  free_clines(map);
}
/*}}}*/

#define FONT "\\tiny"
void emit_tikz(double level, const char *colour, const char *thickness, struct cline *lines)/*{{{*/
{
  struct cline *line;
  double ee[MAXP], nn[MAXP];
  unsigned char keep[MAXP];
  int q, i;

  for (line=lines; line; line=line->next) {
    struct cpoint *p;
    struct cpoint *fp, *np;
    double minE, minN, maxE, maxN;
    double mhd;
    double E, N;
    int is_cycle;
    int limit;
    q = 0;
    fp = line->points;
    for (p=line->points; p; p=np) {
      np = p->next;
      mxy_to_en(p->x, p->y, &E, &N);
      ee[q] = E/100000.0;
      nn[q] = N/100000.0;
      q++;
    }

    is_cycle = reduce(q, ee, nn, keep, 0.01);
    fprintf(out, "\\path [draw,%s,%s] ", thickness, colour);
    fprintf(out, "(%.2fcm,%.2fcm)\n", ee[0], nn[0]);
    minE = maxE = ee[0];
    minN = maxN = nn[0];
    limit = is_cycle ? q-1 : q;
    for (i=1; i<limit; i++) {
      if (keep[i]) {
        fprintf(out, " -- (%.2fcm,%.2fcm)\n", ee[i], nn[i]);
        if (ee[i] < minE) minE = ee[i];
        if (ee[i] > maxE) maxE = ee[i];
        if (nn[i] < minN) minN = nn[i];
        if (nn[i] > maxN) maxN = nn[i];
      }
    }
    if (is_cycle) {
      fprintf(out, " -- cycle\n");
    }
    fprintf(out, ";\n");


    mhd = (maxE - minE) + (maxN - minN);
    if (mhd > 0.7) {
      double E, N;
      char buffer[32];
      mxy_to_en(fp->x, fp->y, &E, &N);
      E /= 100000.0;
      N /= 100000.0;
      /* Draw caption on the line */
      if (level < 1.0) {
        sprintf(buffer, "%.2f", level);
      } else if (level < 10.0) {
        sprintf(buffer, "%.1f", level);
      } else {
        sprintf(buffer, "%.0f", level);
      }
      fprintf(out, "\\path (%.2fcm,%.2fcm) node "
          "[fill=white,fill opacity=0.7,"
          " text opacity=1,font=%s,text=%s,"
          "inner sep=1pt,draw=%s,thin,draw opacity=1] {%s};\n",
            E, N, FONT, colour, colour, buffer);
    }
  }

  return;
}
/*}}}*/
void finish_tikz(void)/*{{{*/
{
  fprintf(out, "\\end{scope}\n");
  fprintf(out, "\\end{tikzpicture}\n");
  fclose(out);
}
/*}}}*/
static void key_heading(FILE *out)/*{{{*/
{
  fprintf(out, "\\begin{tabular}{r l}\n");
  fprintf(out, "\\toprule\n");
  fprintf(out, "\\textbf{Level} & \\textbf{Sample} \\\\\n");
  fprintf(out, "\\midrule\n");
}
/*}}}*/
static void key_footing(FILE *out)/*{{{*/
{
  fprintf(out, "\\bottomrule\n");
  fprintf(out, "\\end{tabular}\n");
}
/*}}}*/
void generate_contour_key_tikz(FILE *out)/*{{{*/
{
  int i;
  int percol;
  percol = (N_LEVELS_1 + 3) >> 2;

  key_heading(out);
  for (i=0; i<N_LEVELS_1; i++) {
    char buffer[32];
    double level = levels_1side[i].level;
    if (level < 1.0) {
      sprintf(buffer, "%.2f", level);
    } else if (level < 10.0) {
      sprintf(buffer, "%.1f", level);
    } else {
      sprintf(buffer, "%.0f", level);
    }

    fprintf(out, "%.2f & {\\tikz \\draw [%s,%s] (0cm,1ex) ",
        level,
        levels_1side[i].colour,
        levels_1side[i].thickness);
    fprintf(out, "node[fill=white,fill opacity=0.7,text opacity=1,font=%s,text=%s,inner sep=1pt,draw=%s,thin,draw opacity=1] {%s}",
        FONT, levels_1side[i].colour, levels_1side[i].colour, buffer);
    fprintf(out, "-- ++(1cm,0cm);}\\\\\n");
    if (((i+1) % percol) == 0) {
      key_footing(out);
      fprintf(out, "\\quad\n");
      key_heading(out);
    }
  }
  key_footing(out);
}
/*}}}*/
