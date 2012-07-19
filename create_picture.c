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

/* Generate picture of OS grid, with WM grid drawn distorted on top. */

#include <stdio.h>

#include "contour.h"

static void write_prologue(void)
{
  struct cline *map;
  printf("\\begin{tikzpicture}\n");
  printf("\\begin{scope}\n");

  map = load_map("uk_coast/uk_coast.txt");
  draw_map_tikz(stdout, map);
  free_clines(map);
}

static void write_epilogue(void)
{
  printf("\\end{scope}\n");
  printf("\\end{tikzpicture}\n");
}

#define RIGHT 7
#define TOP 13
#define FSIZE_GRID "tiny"
#define COL "green!80!black"

#define XORIG 0.494440093
#define YORIG 0.312663855
#define XSCALE (1.0/61.0)
#define YSCALE (1.0/36.0)

static void write_verticals(void)
{
  int i;
  for (i=0; i<=RIGHT; i++) {
    printf("\\path [draw,thin,black] (%dcm,0cm) -- (%dcm,%dcm);\n", i, i, TOP);
    printf("\\path (%dcm,-0.3cm) node [font=\\%s] {%d};\n", i, FSIZE_GRID, i*100000);
  }
}

static void write_horizontals(void)
{
  int i;
  for (i=0; i<=TOP; i++) {
    printf("\\path [draw,thin,black] (0cm,%dcm) -- (%dcm,%dcm);\n", i, RIGHT, i);
    printf("\\path (-0.6cm,%dcm) node [font=\\%s] {%d};\n", i, FSIZE_GRID, i*100000);
  }
}

static char *s0 = "VQLFAWRMGBXSNHCYTOJDZUPKE";
static char *s1 = "SNHTOJ";

static void write_squares(void)
{
  int i, j;
  for (i=0; i<RIGHT; i++) {
    for (j=0; j<TOP; j++) {
      printf("\\path (%.1fcm,%.1fcm) node {%c%c};\n",
          (double)i+0.5,
          (double)j+0.5,
          s1[3*(i/5) + (j/5)],
          s0[5*(i%5) + (j%5)]
          );
    }
  }
}

static void wm_to_en(double X, double Y, double *E, double *N) {
    double x, y;
    x = 72.0 * (X - 0.494444);
    y = 36.0 * (Y - 0.318700);

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


static void write_wm_horiz(void)
{
  int i;
  double j;
  double X, Y, E, N;
  int first;
  double step=0.1;
  for (i=18; i<=22; i++) {
    first = 1;
    printf("\\path [draw,thick," COL "]\n");
    for (j=30; j<=33.001; j+=step) {
      if (first) {
        first = 0;
      } else {
        printf(" --");
      }
      X = j / 64.0;
      Y = (double) i / 64.0;
      wm_to_en(X, Y, &E, &N);
      printf(" (%.2fcm,%.2fcm)\n", E/100000.0, N/100000.0);
    }
    printf(";\n");
    if (i == 22) continue;
    X = 30.0 / 64.0;
    Y = (0.5 + (double) i) / 64.0;
    wm_to_en(X, Y, &E, &N);
    printf("\\path (%.1fcm,%.1fcm) node[" COL "] {%d};\n",
        -0.3 + E/100000.0, N/100000.0, i);
    X = 33.0 / 64.0;
    Y = (0.5 + (double) i) / 64.0;
    wm_to_en(X, Y, &E, &N);
    printf("\\path (%.1fcm,%.1fcm) node[" COL "] {%d};\n",
        0.3 + E/100000.0, N/100000.0, i);
  }
}

static void write_wm_vert(void)
{
  int i;
  double j;
  double X, Y, E, N;
  int first;
  double step=0.1;
  for (i=30; i<=33; i++) {
    first = 1;
    printf("\\path [draw,thick," COL "]\n");
    for (j=18; j<=22.001; j+=step) {
      if (first) {
        first = 0;
      } else {
        printf(" --");
      }
      X = (double) i / 64.0;
      Y = (double) j / 64.0;
      wm_to_en(X, Y, &E, &N);
      printf(" (%.2fcm,%.2fcm)\n", E/100000.0, N/100000.0);
    }
    printf(";\n");
    if (i == 33) continue;
    X = (0.5 + (double)i) / 64.0;
    Y = 18.0 / 64.0;
    wm_to_en(X, Y, &E, &N);
    printf("\\path (%.1fcm,%.1fcm) node[" COL "] {%d};\n",
        E/100000.0, 0.2 + N/100000.0, i);
    X = (0.5 + (double) i) / 64.0;
    Y = 22.0 / 64.0;
    wm_to_en(X, Y, &E, &N);
    printf("\\path (%.1fcm,%.1fcm) node[" COL "] {%d};\n",
        E/100000.0, -0.2 + N/100000.0, i);
  }
}

static void write_ab_orig(void)
{
  double X, Y, E, N;
  X = XORIG;
  Y = YORIG;
  wm_to_en(X, Y, &E, &N);
  printf("\\draw [very thick,blue] (%.3fcm,%.3fcm) circle(0.1cm);\n",
      E/100000.0, N/100000.0);
}


int main (int argc, char **argv)
{
  write_prologue();

  write_verticals();
  write_horizontals();
  write_squares();

  write_wm_horiz();
  write_wm_vert();

  write_ab_orig();

  write_epilogue();
  return 0;
}

