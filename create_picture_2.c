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

/* Generate picture which is linear over WM space, with the OS grid and map
 * drawn distorted on top. */

#include <stdio.h>

#include "contour.h"

static void write_prologue(void)/*{{{*/
{
  struct cline *map;
  printf("\\begin{tikzpicture}\n");
  printf("\\begin{scope}\n");

  map = load_map("uk_coast/uk_coast.txt");
  draw_map_mxy_tikz(stdout, map);
  free_clines(map);
}
/*}}}*/
static void write_epilogue(void)/*{{{*/
{
  printf("\\end{scope}\n");
  printf("\\end{tikzpicture}\n");
}
/*}}}*/

#define RIGHT 7
#define TOP 13
#define FSIZE_GRID "tiny"
#define COL "green!80!black"

static void write_verticals(void)/*{{{*/
{
  int i;
  double X, Ya, Yb, x0, x1, y0, y1;
  Ya = 22/REF6;
  Yb = 18/REF6;
  for (i=30; i<=33; i++) {
    X = i/REF6;
    XY_to_xy(X, Ya, &x0, &y0);
    XY_to_xy(X, Yb, &x1, &y1);
    printf("\\path [draw,thick," COL "] (%.2fcm,%.2fcm) -- (%.2fcm,%.2fcm);\n",
        x0, y0, x1, y1);
  }

  Ya = 22.1/REF6;
  Yb = 17.9/REF6;
  for (X = 0.47; X<0.51001; X+=0.01) {
    XY_to_xy(X, Ya, &x0, &y0);
    XY_to_xy(X, Yb, &x1, &y1);
    printf("\\path [draw,thick,red] (%.2fcm,%.2fcm) -- (%.2fcm,%.2fcm);\n",
        x0, y0, x1, y1);
    printf("\\path (%.2fcm,%.2fcm) node [red,font=\\%s] {%.2f};\n",
        x0, y0-0.2, FSIZE_GRID, X);
    printf("\\path (%.2fcm,%.2fcm) node [red,font=\\%s] {%.2f};\n",
        x1, y1+0.2, FSIZE_GRID, X);
  }

}
/*}}}*/
static void write_horizontals(void)/*{{{*/
{
  int i;
  double Y, Xa, Xb, x0, x1, y0, y1;
  Xa = 30/REF6;
  Xb = 33/REF6;
  for (i=18; i<=22; i++) {
    Y = i/REF6;
    XY_to_xy(Xa, Y, &x0, &y0);
    XY_to_xy(Xb, Y, &x1, &y1);
    printf("\\path [draw,thick," COL "] (%.2fcm,%.2fcm) -- (%.2fcm,%.2fcm);\n",
        x0, y0, x1, y1);
  }
  Xa = 29.9/REF6;
  Xb = 33.1/REF6;
  for (Y = 0.29; Y<0.34001; Y+=0.01) {
    XY_to_xy(Xa, Y, &x0, &y0);
    XY_to_xy(Xb, Y, &x1, &y1);
    printf("\\path [draw,thick,red] (%.2fcm,%.2fcm) -- (%.2fcm,%.2fcm);\n",
        x0, y0, x1, y1);
    printf("\\path (%.2fcm,%.2fcm) node [red,font=\\%s] {%.2f};\n",
        x0-0.3, y0, FSIZE_GRID, Y);
    printf("\\path (%.2fcm,%.2fcm) node [red,font=\\%s] {%.2f};\n",
        x1+0.3, y1, FSIZE_GRID, Y);
  }
}
/*}}}*/

#define XORIG 0.494440093
#define YORIG 0.312663855
#define XSCALE (1.0/61.0)
#define YSCALE (1.0/36.0)

static void draw_limits(void)
{
  double X0, Y0, X1, Y1, x0, x1, y0, y1;
  X0 = XORIG - XSCALE;
  X1 = XORIG + XSCALE;
  Y0 = YORIG - YSCALE;
  Y1 = YORIG + YSCALE;
  XY_to_xy(X0, Y0, &x0, &y0);
  XY_to_xy(X1, Y1, &x1, &y1);
  printf("\\path [draw,thick,blue] (%.2fcm,%.2fcm)\n"
      "-- (%.2fcm,%.2fcm)\n"
      "-- (%.2fcm,%.2fcm)\n"
      "-- (%.2fcm,%.2fcm)\n"
      "-- cycle;\n",
      x0, y0, x0, y1, x1, y1, x1, y0);
}

static char *s0 = "VQLFAWRMGBXSNHCYTOJDZUPKE";
static char *s1 = "SNHTOJ";

static void write_squares(void)/*{{{*/
{
  int i, j;
  double E, N, X, Y, x, y;
  for (i=0; i<RIGHT; i++) {
    for (j=0; j<TOP; j++) {
      E = i*100000 + 50000;
      N = j*100000 + 50000;
      en_to_wm(E, N, &X, &Y);
      XY_to_xy(X, Y, &x, &y);
      printf("\\path (%.2fcm,%.2fcm) node {%c%c};\n",
          x, y,
          s1[3*(i/5) + (j/5)],
          s0[5*(i%5) + (j%5)]
          );
    }
  }
}
/*}}}*/

static void write_en_horiz(void)/*{{{*/
{
  int i;
  double j;
  double X, Y, E, N;
  double x, y;
  int first;
  double step=0.1;
  for (i=0; i<=13; i++) {
    first = 1;
    printf("\\path [draw,thin,black]\n");
    for (j=0; j<=7.001; j+=step) {
      if (first) {
        first = 0;
      } else {
        printf(" --");
      }
      E = j*100000.0;
      N = i*100000.0;
      en_to_wm(E, N, &X, &Y);
      XY_to_xy(X, Y, &x, &y);
      printf(" (%.2fcm,%.2fcm)\n", x, y);
    }
    printf(";\n");
    E = 750000;
    N = i*100000.0;
    en_to_wm(E, N, &X, &Y);
    XY_to_xy(X, Y, &x, &y);
    printf("\\path (%.2fcm,%.2fcm) node [font=\\%s] {%d};\n", x, y, FSIZE_GRID, i*100000);
    E = -50000;
    N = i*100000.0;
    en_to_wm(E, N, &X, &Y);
    XY_to_xy(X, Y, &x, &y);
    printf("\\path (%.2fcm,%.2fcm) node [font=\\%s] {%d};\n", x, y, FSIZE_GRID, i*100000);
  }
}
/*}}}*/
static void write_en_vert(void)/*{{{*/
{
  int i;
  double j;
  double X, Y, E, N;
  double x, y;
  int first;
  double step=0.1;
  for (i=0; i<=7; i++) {
    first = 1;
    printf("\\path [draw,thin,black]\n");
    for (j=0; j<=13.001; j+=step) {
      if (first) {
        first = 0;
      } else {
        printf(" --");
      }
      E = i*100000.0;
      N = j*100000.0;
      en_to_wm(E, N, &X, &Y);
      XY_to_xy(X, Y, &x, &y);
      printf(" (%.2fcm,%.2fcm)\n", x, y);
    }
    printf(";\n");
    E = i*100000;
    N = 13.2*100000.0;
    en_to_wm(E, N, &X, &Y);
    XY_to_xy(X, Y, &x, &y);
    printf("\\path (%.2fcm,%.2fcm) node [font=\\%s] {%d};\n", x, y, FSIZE_GRID, i*100000);
    E = i*100000;
    N = -0.2*100000.0;
    en_to_wm(E, N, &X, &Y);
    XY_to_xy(X, Y, &x, &y);
    printf("\\path (%.2fcm,%.2fcm) node [font=\\%s] {%d};\n", x, y, FSIZE_GRID, i*100000);
  }
}
/*}}}*/
static void write_ab_orig(void)/*{{{*/
{
  double X, Y, x, y;
  X = XORIG;
  Y = YORIG;
  XY_to_xy(X, Y, &x, &y);
  printf("\\draw [very thick,blue] (%.3fcm,%.3fcm) circle(0.1cm);\n",
      x, y);
}
/*}}}*/

int main (int argc, char **argv)/*{{{*/
{
  write_prologue();

  write_verticals();
  write_horizontals();
  draw_limits();
  write_squares();

  write_en_horiz();
  write_en_vert();

  write_ab_orig();

  write_epilogue();
  return 0;
}
/*}}}*/
