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

/* Find some approximations to go from Web-Mercator X,Y to OS national grid
 * eastings, northings */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "contour.h"
#include "osxx02.h"

#if defined(ONE_METRE_MODEL)
/* This is accurate to <1m over all the OSTN02 defined area, but the
 * coefficients can't be dropped in any kind of graceful way to get
 * significantly simpler but less accurate models. */
#define OX 8
#define OY 8

#define QX 6
#define QY 6
#else
/* The 'regular' mode which gracefully degrades as terms are dropped */
#define OX 6
#define OY 6

#define QX 4
#define QY 4
#endif

#define OXY (OX*OY)
#define QXY (QX*QY)



struct model/*{{{*/
{

  double X0, Y0;
  double SX, SY;
  /* TODO : not even clear that the same order polynomials are appropriate for
   * both of these? */
  double e[OX][OY];
  double n[OX][OY];
};
/*}}}*/
static void fit(struct model *model)/*{{{*/
{
  double lat0 = 54.5;
  double lon0 = -2.0;
  double SX, SY;
  double latmin = 49.8;
  double latmax = 61.0;
  double lonmin, lonmax;
  double step = 0.025;

  double lat, lon;
  struct llh osgb, wgs;
  struct mxy mxy;
  double xmin, ymin, xmax, ymax;
  struct en en;
  double dx, dy;
  long double xx[2*OX-1], yy[2*OY-1];
  long double LE[OXY], LN[OXY], RE[OXY], RN[OXY];
  Matrix M;
  int i, j, m, n;

  for (i=0; i<OXY; i++) {
    for (j=0; j<OXY; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<OXY; i++) {
    RE[i] = 0.0;
    RN[i] = 0.0;
  }

  osgb.lat = lat0;
  osgb.lon = lon0;
  osgb.h = 0.0;
  osgb36_to_wgs84(&osgb, &wgs);
  wgs84_to_mxy(&wgs, &mxy);
  model->X0 = mxy.X;
  model->Y0 = mxy.Y;

  osgb.lat = 50;
  osgb.lon = -8;
  osgb.h = 0.0;
  osgb36_to_wgs84(&osgb, &wgs);
  wgs84_to_mxy(&wgs, &mxy);
  xmin = mxy.X;
  ymin = mxy.Y;
  osgb.lat = 61;
  osgb.lon = +2;
  osgb.h = 0.0;
  osgb36_to_wgs84(&osgb, &wgs);
  wgs84_to_mxy(&wgs, &mxy);
  xmax = mxy.X;
  ymax = mxy.Y;

  model->SX = SX = 2.0 / fabs(xmax - xmin);
  model->SY = SY = 2.0 / fabs(ymax - ymin);

#if 1
  printf("Original X0 : %.12f\n", model->X0);
  printf("Original Y0 : %.12f\n", model->Y0);
  printf("Original SX : %.12f\n", model->SX);
  printf("Original SY : %.12f\n", model->SY);
  model->X0 = 0.494440093;
  model->Y0 = 0.312663855; /* 0.3187 */;
  model->SX = 61.0;
  model->SY = 36.0;
#endif

  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    uk_range(lat, &lonmin, &lonmax);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      osgb.lat = lat;
      osgb.lon = lon;
      osgb.h = 0;
      osgb36_to_grid(&osgb, &en);
      osgb36_to_wgs84(&osgb, &wgs);
      wgs84_to_mxy(&wgs, &mxy);
      dx = model->SX * (mxy.X - model->X0);
      dy = model->SY * (mxy.Y - model->Y0);
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
      for (i=0; i<OY; i++) {
        for (j=0; j<OY; j++) {
          for (m=0; m<OX; m++) {
            for (n=0; n<OX; n++) {
              M[OX*i+m][OX*j+n] += xx[m+n]*yy[i+j];
            }
          }
        }
      }
      for (i=0; i<OX; i++) {
        for (j=0; j<OY; j++) {
          RE[i+OX*j] += en.E * xx[i] * yy[j];
          RN[i+OX*j] += en.N * xx[i] * yy[j];
        }
      }
    }
  }

  solve2(M, OXY, LE, LN, RE, RN);

  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->e[i][j] = LE[i+OX*j];
      model->n[i][j] = LN[i+OX*j];
    }
  }
}
/*}}}*/
static void predict(double x, double y,/*{{{*/
    const struct model *model,
    double *E, double *N)
{
  double dx, dy;
  double xx[2*OX-1], yy[2*OY-1];
  int i, j;
  double ee, nn;
  dx = model->SX * (x - model->X0);
  dy = model->SY * (y - model->Y0);
  xx[0] = 1.0;
  for (i=1; i<OX; i++) {
    xx[i] = dx * xx[i-1];
  }
  yy[0] = 1.0;
  for (i=1; i<OY; i++) {
    yy[i] = dy * yy[i-1];
  }
  ee = nn = 0.0;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      ee += model->e[i][j] * xx[i] * yy[j];
      nn += model->n[i][j] * xx[i] * yy[j];
    }
  }
  *E = ee;
  *N = nn;
}
/*}}}*/
static void predict_ll(double lat, double lon,/*{{{*/
    const struct model *model,
    double *E, double *N,
    double *e_err, double *n_err)
{
  struct llh osgb, wgs;
  struct mxy mxy;
  struct en en;

  osgb.lat = lat;
  osgb.lon = lon;
  osgb.h = 0;
  osgb36_to_wgs84(&osgb, &wgs);
  osgb36_to_grid(&osgb, &en);
  wgs84_to_mxy(&wgs, &mxy);
  predict(mxy.X, mxy.Y, model, E, N);
  *e_err = en.E - *E;
  *n_err = en.N - *N;
}
/*}}}*/
static int predict_ll_ostn(double lat, double lon,/*{{{*/
    const struct model *model,
    double *E, double *N,
    double *e_err, double *n_err)
{
  struct llh wgs;
  struct mxy mxy;
  struct en en_wgs;
  struct en en_osgb;

  wgs.lat = lat;
  wgs.lon = lon;
  wgs.h = 0;
  wgs84_to_grid(&wgs, &en_wgs);
  if (apply_ostn02(&en_wgs, &en_osgb) == 0) return 0;
  wgs84_to_mxy(&wgs, &mxy);
  predict(mxy.X, mxy.Y, model, E, N);
  *e_err = en_osgb.E - *E;
  *n_err = en_osgb.N - *N;
  return 1;
}
/*}}}*/

static void trim_model(struct model *model)/*{{{*/
{
  /* Create a low-order fit to the residual between the initial model and the
   * OSTN02 model.  Solve this and use it to trim the coefficients */

  double latmin = 49.8;
  double latmax = 61.0;
  double lonmin, lonmax;
  double step = 0.025;
  double lat, lon;
  struct llh wgs;
  struct mxy mxy;
  double dx, dy;
  double min_dx, min_dy;
  double max_dx, max_dy;

  long double xx[2*QX-1], yy[2*QY-1];
  long double LE[QXY], LN[QXY], RE[QXY], RN[QXY];
  Matrix M;
  double E, N, e_err, n_err;
  int i, j, m, n;

  min_dx = 1.0;
  max_dx = -1.0;
  min_dy = 1.0;
  max_dy = -1.0;

  for (i=0; i<QXY; i++) {
    for (j=0; j<QXY; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<QXY; i++) {
    RE[i] = 0.0;
    RN[i] = 0.0;
  }

  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    uk_range(lat, &lonmin, &lonmax);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      wgs.lat = lat;
      wgs.lon = lon;
      wgs.h = 0;
      if(predict_ll_ostn(lat, lon, model, &E, &N, &e_err, &n_err) == 0) continue;
      wgs84_to_mxy(&wgs, &mxy);
      dx = model->SX * (mxy.X - model->X0);
      dy = model->SY * (mxy.Y - model->Y0);
      if (dx < min_dx) min_dx = dx;
      if (dx > max_dx) max_dx = dx;
      if (dy < min_dy) min_dy = dy;
      if (dy > max_dy) max_dy = dy;
      xx[0] = 1.0;
      for (i=1; i<2*QX-1; i++) {
        xx[i] = dx * xx[i-1];
      }
      yy[0] = 1.0;
      for (i=1; i<2*QY-1; i++) {
        yy[i] = dy * yy[i-1];
      }

      for (i=0; i<QY; i++) {
        for (j=0; j<QY; j++) {
          for (m=0; m<QX; m++) {
            for (n=0; n<QX; n++) {
              M[QX*i+m][QX*j+n] += xx[m+n]*yy[i+j];
            }
          }
        }
      }

      for (i=0; i<QX; i++) {
        for (j=0; j<QY; j++) {
          RE[i+QX*j] += e_err * xx[i] * yy[j];
          RN[i+QX*j] += n_err * xx[i] * yy[j];
        }
      }

    }
  }

  solve2(M, QXY, LE, LN, RE, RN);

  printf("Trim coefficients:\n");
  for (i=0; i<QX; i++) {
    for (j=0; j<QY; j++) {
      model->e[i][j] += LE[i+QX*j];
      model->n[i][j] += LN[i+QX*j];
      printf("delta_e[%d][%d] = %15.5Lf   delta_n[%d][%d] = %15.5Lf\n",
          i, j,
          LE[i+QX*j],
          i, j,
          LN[i+QX*j]);
    }
  }
  printf("\\begin{tabular}{l l r r}\n");
  printf("\\toprule\n");
  printf("\\textbf{i} & \\textbf{j} & $\\Delta{}e_{ij}$ & $\\Delta{}n_{ij}$ \\\\\n");
  printf("\\midrule\n");
  for (i=0; i<QX; i++) {
    for (j=0; j<QY; j++) {
      printf("%d & %d & %15.5Lf & %15.5Lf \\\\\n",
          i, j,
          LE[i+QX*j],
          LN[i+QX*j]
          );
    }
  }
  printf("\\bottomrule\n");
  printf("\\end{tabular}\n");
  printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
  printf("X and Y ranges encountered over OSTN area:\n");
  printf("  X : %.12f -> %.12f\n", min_dx, max_dx);
  printf("  Y : %.12f -> %.12f\n", min_dy, max_dy);
  printf("=-=-=-=-=-=-=-=-=-=-=-=-=-=-\n");
}
/*}}}*/
static void print_model(const struct model *model) {/*{{{*/
  int i, j;
  int nze=0, nzn=0;
  printf("X0=%.10f\n", model->X0);
  printf("Y0=%.10f\n", model->Y0);
  printf("SX=%.3f\n", model->SX);
  printf("SY=%.3f\n", model->SY);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      printf("e[%d][%d] = %15.5f   n[%d][%d] = %15.5f\n",
          i, j,
          model->e[i][j],
          i, j,
          model->n[i][j]
          );
      if (fabs(model->e[i][j]) > 0.000001) nze++;
      if (fabs(model->n[i][j]) > 0.000001) nzn++;
    }
  }
  printf("%2d eastings terms, %2d northings terms\n", nze, nzn);
}
/*}}}*/
static void print_model_tex(const struct model *model) {/*{{{*/
  int i, j;
  printf("\\begin{tabular}{l l r r}\n");
  printf("\\toprule\n");
  printf("\\textbf{i} & \\textbf{j} & $e_{ij}$ & $n_{ij}$ \\\\\n");
  printf("\\midrule\n");
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      printf("%d & %d & %15.5f & %15.5f \\\\\n",
          i, j,
          model->e[i][j],
          model->n[i][j]
          );
    }
  }
  printf("\\bottomrule\n");
  printf("\\end{tabular}\n");
}
/*}}}*/
static void print_model_java_inner(char foo, const double m[OX][OY]) {/*{{{*/
  int i, j;
  unsigned char any[16];
  for (i=0; i<OX; i++) {
    any[i] = 0;
    for (j=0; j<OY; j++) {
      if (fabs(m[i][j]) > 1.0e-10) {
        if (any[i] == 0) {
          printf("    double %c%d =", foo, i);
        }
        any[i] = 1;
        printf(" %+.2f", m[i][j]);
        if (j > 0) printf("*beta");
        if (j > 1) printf("%d", j);
      }
    }
    if (any[i]) {
      printf(";\n");
    }
  }
  printf("    double %c = ", foo);
  for (i=0; i<OX; i++) {
    if (any[i]) {
      if (i > 0) printf(" +");
      printf("%c%d", foo, i);
      if (i > 0) printf("*alpha");
      if (i > 1) printf("%d", i);
    }
  }
  printf(";\n");
}
/*}}}*/
static void print_model_java(const struct model *model) {/*{{{*/
  printf("double alpha, beta;\n");
  printf("alpha = %.3f * (X - %.10f);\n", model->SX, model->X0);
  printf("beta = %.3f * (Y - %.10f);\n", model->SY, model->Y0);
  print_model_java_inner('E', model->e);
  print_model_java_inner('N', model->n);
}
/*}}}*/
static void print_model_estrin_inner(const double m[OX][OY], const char *res)/*{{{*/
{
  SMatrix s;
  int i;
  int j;
  memset(s, 0, sizeof(s));
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      double term = m[i][j];
      if (fabs(term) > 0.0) {
        char buffer[32];
        sprintf(buffer, "%.2f", term);
        s[i][j] = strdup(buffer);
      }
    }
  }
  estrin2(s, OX, OY, "alpha", "beta", res);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      if (s[i][j]) free(s[i][j]);
    }
  }
}
/*}}}*/
static void print_model_estrin(const struct model *model)/*{{{*/
{
  printf("double alpha, beta;\n");
  printf("alpha = %.3f * (X - %.10f);\n", model->SX, model->X0);
  printf("beta = %.3f * (Y - %.10f);\n", model->SY, model->Y0);
  print_model_estrin_inner(model->e, "E");
  print_model_estrin_inner(model->n, "N");
}
/*}}}*/
static void print_model_latex_inner(char foo, const double m[OX][OY]) {/*{{{*/
  int i, j;
  unsigned char any[16];
  for (i=0; i<OX; i++) {
    any[i] = 0;
    for (j=0; j<OY; j++) {
      if (fabs(m[i][j]) > 1.0e-10) {
        if (any[i] == 0) {
          printf("%c_%d &=", foo, i);
        }
        if (any[i]) {
          printf(" %+.2f", m[i][j]);
        } else {
          printf(" %.2f", m[i][j]);
        }
        any[i] = 1;
        if (j > 0) printf("\\beta");
        if (j > 1) printf("^%d", j);
      }
    }
    if (any[i]) {
      printf("\\nonumber \\\\\n");
    }
  }
  printf("%c &= ", foo);
  for (i=0; i<OX; i++) {
    if (any[i]) {
      if (i > 0) printf(" +");
      printf("%c_%d", foo, i);
      if (i > 0) printf("\\alpha");
      if (i > 1) printf("^%d", i);
    }
  }
  printf("\n");
}
/*}}}*/
static void print_model_latex(const struct model *model) {/*{{{*/
  printf("\\begin{align}\n");
  print_model_latex_inner('E', model->e);
  printf("\\\\[1ex]\n");
  print_model_latex_inner('N', model->n);
  printf("\\end{align}\n");
}
/*}}}*/
static void drop_coef(int bw, int bi, int bj, double c[OX][OY], double accum[OX][OY])/*{{{*/
{
  double coef;
  coef = c[bi][bj];
  switch (bw ? bj : bi) {
    case 0:
    case 1:
      break; /* Just discard the coefficient */
    case 2:
      if (bw) {
        c[bi][0] += 0.5*coef;
      } else {
        c[0][bj] += 0.5*coef;
      }
      break;
    case 3:
      if (bw) {
        c[bi][1] += 0.75*coef;
      } else {
        c[1][bj] += 0.75*coef;
      }
      break;
    case 4:
      if (bw) {
        c[bi][2] += coef;
        c[bi][0] -= 0.125*coef;
      } else {
        c[2][bj] += coef;
        c[0][bj] -= 0.125*coef;
      }
      break;
    case 5:
      if (bw) {
        c[bi][3] += (20.0/16.0)*coef;
        c[bi][1] -= (5.0/16.0)*coef;
      } else {
        c[3][bj] += (20.0/16.0)*coef;
        c[1][bj] -= (5.0/16.0)*coef;
      }
      break;
    case 6:
      if (bw) {
        c[bi][4] += (48.0/32.0)*coef;
        c[bi][2] -= (18.0/32.0)*coef;
        c[bi][0] += (1.0/32.0)*coef;
      } else {
        c[4][bj] += (48.0/32.0)*coef;
        c[2][bj] -= (18.0/32.0)*coef;
        c[0][bj] += (1.0/32.0)*coef;
      }
      break;
    case 7:
      if (bw) {
        c[bi][5] += (112.0/64.0)*coef;
        c[bi][3] -= (56.0/64.0)*coef;
        c[bi][1] += (7.0/64.0)*coef;
      } else {
        c[5][bj] += (112.0/64.0)*coef;
        c[3][bj] -= (56.0/64.0)*coef;
        c[1][bj] += (7.0/64.0)*coef;
      }
      break;
    case 8:
      if (bw) {
#if (OY > 6)
        c[bi][6] += (256.0/128.0)*coef;
        c[bi][4] -= (160.0/128.0)*coef;
        c[bi][2] += (32.0/128.0)*coef;
        c[bi][0] -= (1.0/128.0)*coef;
#endif
      } else {
#if (OX > 6)
        c[6][bj] += (256.0/128.0)*coef;
        c[4][bj] -= (160.0/128.0)*coef;
        c[2][bj] += (32.0/128.0)*coef;
        c[0][bj] -= (1.0/128.0)*coef;
#endif
      }
      break;
    case 9:
      if (bw) {
#if (OY > 7)
        c[bi][7] += (576.0/256.0)*coef;
        c[bi][5] -= (432.0/256.0)*coef;
        c[bi][3] += (120.0/256.0)*coef;
        c[bi][1] -= (9.0/256.0)*coef;
#endif
      } else {
#if (OX > 7)
        c[7][bj] += (576.0/256.0)*coef;
        c[5][bj] -= (432.0/256.0)*coef;
        c[3][bj] += (120.0/256.0)*coef;
        c[1][bj] -= (9.0/256.0)*coef;
#endif
      }
      break;
    default:
      fprintf(stderr, "Can't handle a polynomial with order this high (bi=%d bj=%d bw=%d), giving up\n", bi, bj, bw);
      exit(2);
  }
  c[bi][bj] = 0.0;
  if (accum) {
    accum[bi][bj] += coef;
  }
}
/*}}}*/

struct checkpoint/*{{{*/
{
  struct checkpoint *next;
  double X, Y;
  double E, N;
};
/*}}}*/
static struct checkpoint *make_checkpoints(void)/*{{{*/
{
  double lonmin, lonmax;
  double latmin, latmax;
  double lat, lon;
  double step;
  struct checkpoint *result;
  struct llh wgs;
  struct mxy mxy;
  struct en en_wgs;
  struct en en_osgb;

  result = NULL;
  latmin = 49.0;
  latmax = 61.0;
  lonmin = -8.0;
  lonmax = +2.0;
  step = 0.01;

  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      struct checkpoint *new_point;
      wgs.lat = lat;
      wgs.lon = lon;
      wgs.h = 0;
      wgs84_to_grid(&wgs, &en_wgs);
      wgs84_to_mxy(&wgs, &mxy);
      if (apply_ostn02(&en_wgs, &en_osgb) == 0) continue;
      new_point = malloc(sizeof(struct checkpoint));
      new_point->X = mxy.X;
      new_point->Y = mxy.Y;
      new_point->E = en_osgb.E;
      new_point->N = en_osgb.N;
      new_point->next = result;
      result = new_point;
    }
  }
  return result;
}
/*}}}*/
static void free_checkpoints(struct checkpoint *list)/*{{{*/
{
  struct checkpoint *next;
  for (; list; list=next) {
    next = list->next;
    free(list);
  }
}
/*}}}*/
struct victim/*{{{*/
{
  int bi, bj, bw;
  double score;
};
/*}}}*/
static void find_victim(double coef[OX][OY], double accum[OX][OY], struct victim *v)/*{{{*/
{
  int which;
  double score;
  int i, j;

  v->bi = v->bj = v->bw = -1;
  v->score = 0.0;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      double score0;
      if (fabs(coef[i][j]) == 0.0) continue;
      score0 = fabs(coef[i][j] + accum[i][j]) - fabs(accum[i][j]);
      if ((i < 2) && (j < 2)) {
        score = score0;
        which = 0;
      } else if (i > j) {
        score = score0 / (double)(1<<(i-1));
        which = 0;
      } else if (i < j) {
        score = score0 / (double)(1<<(j-1));
        which = 1;
      } else { /* It's arbitrary */
        score = score0 / (double)(1<<(j-1));
        which = 1;
      }

      if ((v->bw < 0) || (score < v->score)) {
        v->bi = i;
        v->bj = j;
        v->bw = which;
        v->score = score;
      }
    }
  }
}
/*}}}*/

#define N_LEVELS 8
static double levels[N_LEVELS] = {
  1.0, 1.5, 2.0, 2.5, 5.0, 10.0, 35.0, 100.0
};

static void scoring(struct model *m, int index, int its_e, struct victim *v, struct checkpoint *expected, const char *tag)/*{{{*/
{
  int count[N_LEVELS+1];
  struct checkpoint *pt;
  int band;
  memset(count, 0, sizeof(count));
  for (pt=expected; pt; pt=pt->next) {
    double E, N;
    double de, dn, err;

    predict(pt->X, pt->Y, m, &E, &N);
    de = E - pt->E;
    dn = N - pt->N;
    err = sqrt(de*de + dn*dn);
    for (band=0; band<N_LEVELS; band++) {
      if (err <= levels[band]) {
        break;
      }
    }
    ++count[band];
  }
  printf("%3d: %4s : %c[%d][%d] in %c (%15.6f) :",
      index,
      tag ? tag : "",
      its_e ? 'E' : 'N',
      v->bi, v->bj,
      v->bw ? 'j' : 'i',
      v->score);
  for (band=0; band<=N_LEVELS; band++) {
    printf(" %8d", count[band]);
  }
  printf("\n");
}
/*}}}*/
static void scoring_header(void)/*{{{*/
{
  int band;
  int i;
  char buffer[32];
  for (i=0; i<45; i++) putchar(' ');
  for (band=0; band<N_LEVELS; band++) {
    sprintf(buffer, "<=%.1f", levels[band]);
    printf("%9s", buffer);
  }
  sprintf(buffer, ">%.1f", levels[N_LEVELS-1]);
  printf("%9s", buffer);
  printf("\n");

}
/*}}}*/

struct save
{
  int index;
  char *tag;
  double min_contour;
  struct model m;
};

static void stepwise_reduction(struct model *m, int n_saves, struct save *saves)/*{{{*/
{
  /* Gradually replace terms involving high powers of x and y with lower power
   * terms, using Chebyshev polynomial coefficients. */
  double ae[OX][OY];
  double an[OX][OY];
  int i, j;
  struct model temp;
  struct checkpoint *expected;
  int index;
  int sp;
  char *tag;
  expected = make_checkpoints();


  memcpy(&temp, m, sizeof(struct model));
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      ae[i][j] = 0.0;
      an[i][j] = 0.0;
    }
  }

  scoring_header();

  sp = 0;
  for (index=0; ; index++) {
    struct victim ve, vn;
    find_victim(temp.e, ae, &ve);
    find_victim(temp.n, an, &vn);

    if ((ve.bw < 0) && (vn.bw < 0)) {
      printf("No eligible candidates\n");
      break;
    } else {
      int its_e;
      if (vn.bw < 0) {
        its_e = 1;
      } else if (ve.bw < 0) {
        its_e = 0;
      } else if (ve.score < vn.score) {
        its_e = 1;
      } else {
        its_e = 0;
      }
      if (its_e) {
        drop_coef(ve.bw, ve.bi, ve.bj, temp.e, ae);
      } else {
        drop_coef(vn.bw, vn.bi, vn.bj, temp.n, an);
      }

      if ((sp < n_saves) && (saves[sp].index == index)) {
        memcpy(&saves[sp].m, &temp, sizeof(struct model));
        tag = saves[sp].tag;
        sp++;
      } else {
        tag = NULL;
      }

      scoring(&temp, index, its_e, its_e ? &ve : &vn, expected, tag);
    }

  }

  free_checkpoints(expected);

}
/*}}}*/
static void puncture_model_4(const char *name, double threshold, double c[OX][OY])/*{{{*/
{
  /* Gradually replace terms involving high powers of x and y with lower power
   * terms, using Chebyshev polynomial coefficients. */
  double total;
  double accum[OX][OY];
  int i, j;

  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      accum[i][j] = 0.0;
    }
  }

  total = 0.0;
  do {
    int which;
    double best_score;
    double score;
    int bi, bj;
    int bw;

    bi = bj = bw = -1;
    best_score = 0.0;
    for (i=0; i<OX; i++) {
      for (j=0; j<OY; j++) {
        double score0;
        if (fabs(c[i][j]) == 0.0) continue;
        score0 = fabs(c[i][j] + accum[i][j]) - fabs(accum[i][j]);
        if ((i < 2) && (j < 2)) {
          score = score0;
          which = 0;
        } else if (i > j) {
          score = score0 / (double)(1<<(i-1));
          which = 0;
        } else if (i < j) {
          score = score0 / (double)(1<<(j-1));
          which = 1;
        } else { /* It's arbitrary */
          score = score0 / (double)(1<<(j-1));
          which = 1;
        }

        if ((bw < 0) || (score < best_score)) {
          bi = i;
          bj = j;
          bw = which;
          best_score = score;
        }
      }
    }

    if (bw < 0) {
      printf("No eligible candidates\n");
      break;
    }
    if ((best_score + total > threshold)) {
      printf("Remaining scores too large, best is %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f\n",
          bw ? "Y" : "X", name,
          bi, bj, c[bi][bj], best_score);
      break;
    }

    total += best_score;
    printf("Drop %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f, total=%15.8f\n",
        bw ? "Y" : "X",
        name,
        bi, bj, c[bi][bj], best_score, total);

    drop_coef(bw, bi, bj, c, accum);
  } while(1);

}
/*}}}*/
#if 0
static void puncture_model_3(const char *name, double threshold, double c[OX][OY])/*{{{*/
{
  /* Gradually replace terms involving high powers of x and y with lower power
   * terms, using Chebyshev polynomial coefficients. */
  double total;

  total = 0.0;
  do {
    int which;
    double best_score;
    double score;
    double coef;
    int i, j;
    int bi, bj;
    int bw;

    bw = -1;
    best_score = 0.0;
    for (i=0; i<OX; i++) {
      for (j=0; j<OY; j++) {
        if (fabs(c[i][j]) == 0.0) continue;
        if ((i < 2) && (j < 2)) {
          score = fabs(c[i][j]);
          which = 0;
        } else if (i > j) {
          score = fabs(c[i][j]) / (double)(1<<(i-1));
          which = 0;
        } else if (i < j) {
          score = fabs(c[i][j]) / (double)(1<<(j-1));
          which = 1;
        } else { /* It's arbitrary */
          score = fabs(c[i][j]) / (double)(1<<(j-1));
          which = 1;
        }

        if ((bw < 0) || (score < best_score)) {
          bi = i;
          bj = j;
          bw = which;
          best_score = score;
        }
      }
    }

    if (bw < 0) {
      printf("No eligible candidates\n");
      break;
    }
    if ((best_score + total > threshold)) {
      printf("Remaining scores too large, best is %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f\n",
          bw ? "Y" : "X", name,
          bi, bj, c[bi][bj], best_score);
      break;
    }

    total += best_score;
    printf("Drop %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f, total=%15.8f\n",
        bw ? "Y" : "X",
        name,
        bi, bj, c[bi][bj], best_score, total);

    coef = c[bi][bj];
    switch (bw ? bj : bi) {
      case 0:
      case 1:
        break; /* Just discard the coefficient */
      case 2:
        if (bw) {
          c[bi][0] += 0.5*coef;
        } else {
          c[0][bj] += 0.5*coef;
        }
        break;
      case 3:
        if (bw) {
          c[bi][1] += 0.75*coef;
        } else {
          c[1][bj] += 0.75*coef;
        }
        break;
      case 4:
        if (bw) {
          c[bi][2] += coef;
          c[bi][0] -= 0.125*coef;
        } else {
          c[2][bj] += coef;
          c[0][bj] -= 0.125*coef;
        }
        break;
      case 5:
        if (bw) {
          c[bi][3] += (20.0/16.0)*coef;
          c[bi][1] -= (5.0/16.0)*coef;
        } else {
          c[3][bj] += (20.0/16.0)*coef;
          c[1][bj] -= (5.0/16.0)*coef;
        }
        break;
      case 6:
        if (bw) {
          c[bi][4] += (48.0/32.0)*coef;
          c[bi][2] -= (18.0/32.0)*coef;
          c[bi][0] += (1.0/32.0)*coef;
        } else {
          c[4][bj] += (48.0/32.0)*coef;
          c[2][bj] -= (18.0/32.0)*coef;
          c[0][bj] += (1.0/32.0)*coef;
        }
        break;
      case 7:
        if (bw) {
          c[bi][5] += (112.0/64.0)*coef;
          c[bi][3] -= (56.0/64.0)*coef;
          c[bi][1] += (7.0/64.0)*coef;
        } else {
          c[5][bj] += (112.0/64.0)*coef;
          c[3][bj] -= (56.0/64.0)*coef;
          c[1][bj] += (7.0/64.0)*coef;
        }
        break;
      case 8:
        if (bw) {
          c[bi][6] += (256.0/128.0)*coef;
          c[bi][4] -= (160.0/128.0)*coef;
          c[bi][2] += (32.0/128.0)*coef;
          c[bi][0] -= (1.0/128.0)*coef;
        } else {
          c[6][bj] += (256.0/128.0)*coef;
          c[4][bj] -= (160.0/128.0)*coef;
          c[2][bj] += (32.0/128.0)*coef;
          c[0][bj] -= (1.0/128.0)*coef;
        }
        break;
      case 9:
        if (bw) {
          c[bi][7] += (576.0/256.0)*coef;
          c[bi][5] -= (432.0/256.0)*coef;
          c[bi][3] += (120.0/256.0)*coef;
          c[bi][1] -= (9.0/256.0)*coef;
        } else {
          c[7][bj] += (576.0/256.0)*coef;
          c[5][bj] -= (432.0/256.0)*coef;
          c[3][bj] += (120.0/256.0)*coef;
          c[1][bj] -= (9.0/256.0)*coef;
        }
        break;
      default:
        fprintf(stderr, "Can't handle a polynomial with order this high, giving up\n");
        exit(2);
    }
    c[bi][bj] = 0.0;
  } while(1);

}
/*}}}*/
#endif
static void puncture_model_2(const char *name, double threshold, double c[OX][OY])/*{{{*/
{
  /* Gradually replace terms involving high powers of x and y with lower power
   * terms, using Chebyshev polynomial coefficients. */
  double total;
  int i, j, ii, jj;
  int flags[OX][OY];
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      flags[i][j] = 1;
    }
  }

  total = 0.0;
  do {
    int best;
    int which;
    double score;
    double coef;

    best = -1;
    for (ii=0; ii<OX; ii++) {
      for (jj=0; jj<OY; jj++) {
        if (flags[ii][jj] == 0) continue;
        if ((ii < 2) && (jj < 2)) continue;
        if (ii + jj > best) {
          best = ii + jj;
          i = ii;
          j = jj;
        }
      }
    }

    if (best < 0) break;

    /* Can't replace constant or linear terms */
    if (i > j) {
      score = fabs(c[i][j]) / (double)(1<<(i-1));
      which = 0;
    } else if (i < j) {
      score = fabs(c[i][j]) / (double)(1<<(j-1));
      which = 1;
    } else { /* It's arbitrary */
      score = fabs(c[i][j]) / (double)(1<<(j-1));
      which = 1;
    }

    if (score + total > threshold) {
      printf("Can't drop %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f\n",
          which ? "Y" : "X",
          name,
          i, j, c[i][j],
          score);
      flags[i][j] = 0;
      continue;
    } else {
      total += score;
      printf("Drop %s coefficient at %s[%d][%d]=%15.8f : score=%15.8f, total=%15.8f\n",
          which ? "Y" : "X",
          name,
          i, j, c[i][j], score, total);

      flags[i][j] = 0;
    }

    coef = c[i][j];
    switch (which ? j : i) {
      case 2:
        if (which) {
          c[i][0] += 0.5*coef;
        } else {
          c[0][j] += 0.5*coef;
        }
        break;
      case 3:
        if (which) {
          c[i][1] += 0.75*coef;
        } else {
          c[1][j] += 0.75*coef;
        }
        break;
      case 4:
        if (which) {
          c[i][2] += coef;
          c[i][0] -= 0.125*coef;
        } else {
          c[2][j] += coef;
          c[0][j] -= 0.125*coef;
        }
        break;
      case 5:
        if (which) {
          c[i][3] += (20.0/16.0)*coef;
          c[i][1] -= (5.0/16.0)*coef;
        } else {
          c[3][j] += (20.0/16.0)*coef;
          c[1][j] -= (5.0/16.0)*coef;
        }
        break;
      default:
        fprintf(stderr, "Can't handle a polynomial with order this high, giving up\n");
        exit(2);
    }
    c[i][j] = 0.0;
  } while(1);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      if (fabs(c[i][j]) < 0.1*threshold) c[i][j] = 0.0;
    }
  }

}
/*}}}*/
#if 0
static void puncture_model(struct model *model)/*{{{*/
{
  int i, j;
  double thresh = 0.25;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      if (fabs(model->e[i][j]) < thresh) model->e[i][j] = 0.0;
      if (fabs(model->n[i][j]) < thresh) model->n[i][j] = 0.0;
    }
  }
}
/*}}}*/
#endif
static void round_model_2(struct model *model)/*{{{*/
{
  int i, j;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->e[i][j] = 0.01 * round(100.0 * model->e[i][j]);
      model->n[i][j] = 0.01 * round(100.0 * model->n[i][j]);
    }
  }
}
/*}}}*/
#if 0
static void round_model_1(struct model *model)/*{{{*/
{
  int i, j;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->e[i][j] = 0.1 * round(10.0 * model->e[i][j]);
      model->n[i][j] = 0.1 * round(10.0 * model->n[i][j]);
    }
  }
}
/*}}}*/
static void round_model_0(struct model *model)/*{{{*/
{
  int i, j;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->e[i][j] = round(model->e[i][j]);
      model->n[i][j] = round(model->n[i][j]);
    }
  }
}
/*}}}*/
#endif

#if 0
static void emit_clines(const struct cline *lines, double lonmin, double lonstep, double latmin, double latstep)/*{{{*/
{
  const struct cline *l;
  for (l=lines; l; l=l->next) {
    struct cpoint *p;
    for (p=l->points; p; p=p->next) {
      printf("%f, %f\n",
          lonmin + lonstep * p->x,
          latmin + latstep * p->y);
    }
    printf("\n");
  }
}
/*}}}*/
#endif
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

struct accuracy/*{{{*/
{
  double e;
  double n;
  double minc;
  char *base;
};
/*}}}*/
#define N_ACCURACIES 5

static struct accuracy accuracies[N_ACCURACIES] =/*{{{*/
{
  {   1.0,   1.9,   0.35, "wm_to_grid_1m%s.%s" },
  {   2.0,   2.0,   0.35, "wm_to_grid_2m%s.%s" },
  {   4.0,   4.0,   0.5,  "wm_to_grid_4m%s.%s" },
  {  20.0,  20.0,   0.7,  "wm_to_grid_20m%s.%s" },
  {  90.0, 100.0,   5.0,  "wm_to_grid_100m%s.%s" }
};
/*}}}*/


#define N_SAVES 7
static struct save saves[N_SAVES] = {
  {  36,   "1m", 0.25, },
  {  47,   "2m", 0.35, },
  {  51,   "2m5", 0.5, },
  {  57,   "5m",  0.7, },
  {  61,  "10m",  1.0, },
  {  67,  "35m",  2.5, },
  {  69, "100m", 10.0, },
};

int main (int argc, char **argv)/*{{{*/
{
  struct model model;
  struct model rounded;
  double step;
  double lat, lon;
  struct cdata *cdata_notrim;
  struct cdata *cdata_trim;
  struct cdata *cdata_total;
  struct cdata *cdata_east;
  struct cdata *cdata_north;
  struct crow *crow_notrim;
  struct crow *crow_trim;
  struct crow *crow_total;
  struct crow *crow_east;
  struct crow *crow_north;
  struct cline *lines;
  double lonmin, lonmax;
  double latmin, latmax;
  double E, N;
  int i, acc;
  int is_svg;

  is_svg = (argc > 1);

  load_osxx02();

  fit(&model);
  step = 0.025;
  latmin = 49.0;
  latmax = 61.0;
  lonmin = -8.0;
  lonmax = +2.0;

  printf("===============================\n");
  printf("INITIAL MODEL FIT:\n");
  print_model(&model);
  print_model_tex(&model);

  rounded = model;
  puncture_model_4("E", 1.0, rounded.e);
  puncture_model_4("N", 1.0, rounded.n);
  round_model_2(&rounded);
  print_model_latex(&rounded);

  /* Plot out the merit of the least squares fit - compared against its
   * original objective function */
  if (!is_svg) {
    cdata_notrim = new_cdata();
    for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
      crow_notrim = next_crow(cdata_notrim);
      for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
        double e_err, n_err;
        double t_err;
        predict_ll(lat, lon, &rounded, &E, &N, &e_err, &n_err);
        t_err = sqrt(e_err*e_err + n_err*n_err);
        add_cnode(crow_notrim, t_err);
      }
    }

    start_tikz("wm_to_grid_notrim_7param.tex");
    for (i=0; i<N_LEVELS_1; i++) {
      if (levels_1side[i].level < 0) continue;
      lines = generate_isolines(cdata_notrim, levels_1side[i].level);
      remap_clines(lines, lonmin, step, latmin, step);
      emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
      free_clines(lines);
    }
    finish_tikz();

    free_cdata(cdata_notrim);
  }

  /* Plot out the merit of the least squares fit - compared against OSTN02 */
  cdata_notrim = new_cdata();
  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    crow_notrim = next_crow(cdata_notrim);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      double e_err, n_err;
      double t_err;
      if (predict_ll_ostn(lat, lon, &rounded, &E, &N, &e_err, &n_err) == 0) {
        add_empty_cnode(crow_notrim);
      } else {
        t_err = sqrt(e_err*e_err + n_err*n_err);
        add_cnode(crow_notrim, t_err);
      }
    }
  }

  if (is_svg) {
    start_svg("wm_to_grid_notrim_ostn.svg");
  } else {
    start_tikz("wm_to_grid_notrim_ostn.tex");
  }
  for (i=0; i<N_LEVELS_1; i++) {
    if (levels_1side[i].level < 0.5) continue;
    lines = generate_isolines(cdata_notrim, levels_1side[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    if (is_svg) {
      emit_svg(levels_1side[i].level, levels_1side[i].svg_colour, levels_1side[i].width_scale, lines);
    } else {
      emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
    }
    free_clines(lines);
  }
  if (is_svg) {
    finish_svg();
  } else {
    finish_tikz();
  }

  free_cdata(cdata_notrim);

  printf("===============================\n");
  trim_model(&model);
  printf("AFTER TRIMMING:\n");
  print_model(&model);
  printf("===============================\n");

  /* Plot out the merit of the least squares fit - compared against OSTN02 after trim */
  cdata_trim = new_cdata();
  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    crow_trim = next_crow(cdata_trim);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      double e_err, n_err;
      double t_err;
      if (predict_ll_ostn(lat, lon, &model, &E, &N, &e_err, &n_err) == 0) {
        add_empty_cnode(crow_trim);
      } else {
        t_err = sqrt(e_err*e_err + n_err*n_err);
        add_cnode(crow_trim, t_err);
      }
    }
  }

  if (is_svg) {
    start_svg("wm_to_grid_trim_ostn.svg");
  } else {
    start_tikz("wm_to_grid_trim_ostn.tex");
  }
  for (i=0; i<N_LEVELS_1; i++) {
    if (levels_1side[i].level < 0.25) continue;
    lines = generate_isolines(cdata_trim, levels_1side[i].level);
    remap_clines(lines, lonmin, step, latmin, step);
    if (is_svg) {
      emit_svg(levels_1side[i].level, levels_1side[i].svg_colour, levels_1side[i].width_scale, lines);
    } else {
      emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
    }
    free_clines(lines);
  }
  if (is_svg) {
    finish_svg();
  } else {
    finish_tikz();
  }
  free_cdata(cdata_trim);

  printf("===============================\n");
  stepwise_reduction(&model, N_SAVES, saves);
  printf("===============================\n");
  for (i=0; i<N_SAVES; i++) {
    int count_e, count_n;
    int j, k;
    struct cdata *cdata;
    struct crow *crow;
    char filename_svg[64], filename_tex[64];

    round_model_2(&saves[i].m);
    count_e = count_n = 0;
    for (j=0; j<OX; j++) {
      for (k=0; k<OY; k++) {
        if (fabs(saves[i].m.e[j][k]) > 0.0) {
          ++count_e;
        }
        if (fabs(saves[i].m.n[j][k]) > 0.0) {
          ++count_n;
        }
      }
    }

    printf("REDUCED MODEL FOR ACCURACY LEVEL : %s\n", saves[i].tag);
    printf("%d non-zero east, %d non-zero north\n", count_e, count_n);
    print_model(&saves[i].m);
    print_model_latex(&saves[i].m);
    print_model_java(&saves[i].m);
    print_model_estrin(&saves[i].m);

    cdata = new_cdata();
    for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
      crow = next_crow(cdata);
      for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
        double e_err, n_err;
        double t_err;
        if (predict_ll_ostn(lat, lon, &saves[i].m, &E, &N, &e_err, &n_err) == 0) {
          add_empty_cnode(crow);
        } else {
          t_err = sqrt(e_err*e_err + n_err*n_err);
          add_cnode(crow, t_err);
        }
      }
    }

    sprintf(filename_svg, "wm_to_grid_%s.svg", saves[i].tag);
    sprintf(filename_tex, "wm_to_grid_%s.tex", saves[i].tag);
    start_svg(filename_svg);
    start_tikz(filename_tex);
    for (j=0; j<N_LEVELS_1; j++) {
      if (levels_1side[j].level < saves[i].min_contour) continue;
      lines = generate_isolines(cdata, levels_1side[j].level);
      remap_clines(lines, lonmin, step, latmin, step);
      emit_svg(levels_1side[j].level, levels_1side[j].svg_colour, levels_1side[j].width_scale, lines);
      emit_tikz(levels_1side[j].level, levels_1side[j].colour, levels_1side[j].thickness, lines);
      free_clines(lines);
    }
    finish_svg();
    finish_tikz();
    free_cdata(cdata);
    printf("===============================\n");
  }

  exit(0);

  for (acc=0; acc<N_ACCURACIES; acc++) {
    struct model copy;
    char filename[64];
    memcpy(&copy, &model, sizeof(struct model));
    printf("===============================\n");
    printf("PUNCTURE : x_thresh=%f y_thresh=%f\n",
        accuracies[acc].e,
        accuracies[acc].n);
    puncture_model_4("E", accuracies[acc].e, copy.e);
    puncture_model_4("N", accuracies[acc].n, copy.n);
    printf("=================\n");
    printf("AFTER PUNCTURING:\n");
    print_model(&copy);

    round_model_2(&copy);
    printf("=================\n");
    printf("AFTER ROUNDING:\n");
    print_model(&copy);
    print_model_latex(&copy);
    print_model_java(&copy);
    print_model_estrin(&copy);

    cdata_total = new_cdata();
    cdata_east = new_cdata();
    cdata_north = new_cdata();
    for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
      crow_total = next_crow(cdata_total);
      crow_east = next_crow(cdata_east);
      crow_north = next_crow(cdata_north);
      for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
        double e_err, n_err;
        double t_err;
        if (predict_ll_ostn(lat, lon, &copy, &E, &N, &e_err, &n_err) == 0) {
          add_empty_cnode(crow_total);
          add_empty_cnode(crow_east);
          add_empty_cnode(crow_north);
        } else {
          t_err = sqrt(e_err*e_err + n_err*n_err);
          add_cnode(crow_total, t_err);
          add_cnode(crow_east, e_err);
          add_cnode(crow_north, n_err);

        }
      }
    }

    sprintf(filename, accuracies[acc].base, "", is_svg ? "svg" : "tex");
    if (is_svg) {
      start_svg(filename);
    } else {
      start_tikz(filename);
    }
    for (i=0; i<N_LEVELS_1; i++) {
      if (levels_1side[i].level < accuracies[acc].minc) continue;
      lines = generate_isolines(cdata_total, levels_1side[i].level);
      remap_clines(lines, lonmin, step, latmin, step);
      if (is_svg) {
        emit_svg(levels_1side[i].level, levels_1side[i].svg_colour, levels_1side[i].width_scale, lines);
      } else {
        emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
      }
      free_clines(lines);
    }
    if (is_svg) {
      finish_svg();
    } else {
      finish_tikz();
    }
    free_cdata(cdata_total);

    sprintf(filename, accuracies[acc].base, "_east", is_svg ? "svg" : "tex");
    if (is_svg) {
      start_svg(filename);
    } else {
      start_tikz(filename);
    }
    for (i=0; i<N_LEVELS_2; i++) {
      double level = levels_2side[i].level;
      if ((fabs(level) > 0) && (fabs(level) < accuracies[acc].minc)) continue;
      lines = generate_isolines(cdata_east, level);
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
    free_cdata(cdata_east);

    sprintf(filename, accuracies[acc].base, "_north", is_svg ? "svg" : "tex");
    if (is_svg) {
      start_svg(filename);
    } else {
      start_tikz(filename);
    }
    for (i=0; i<N_LEVELS_2; i++) {
      double level = levels_2side[i].level;
      if ((fabs(level) > 0) && (fabs(level) < accuracies[acc].minc)) continue;
      lines = generate_isolines(cdata_north, level);
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
    free_cdata(cdata_north);
  }

  return 0;
}
/*}}}*/
