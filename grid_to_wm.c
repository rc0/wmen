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

/* Find some approximations to go from OS national grid E,N to Web-Mercator X,Y
 *
 * Flow is like this:
 * - loops are over WGS84 lat/long space : easy to fix the range of interest
 *
 * */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"
#include "contour.h"
#include "osxx02.h"

#define OX 7
#define OY 7
#define OXY (OX*OY)

#define QX 4
#define QY 4
#define QXY (QX*QY)

#define FIT_STEP 0.01


struct model/*{{{*/
{

  /* Origin and scalings in OSGB grid space */
  double E0, N0;
  double SE, SN;
  double SX, SY;

  /* TODO : not even clear that the same order polynomials are appropriate for
   * both of these? */
  double x[OX][OY];
  double y[OX][OY];
};
/*}}}*/
static void fit(struct model *model)/*{{{*/
{
  double lat0 = 54.5;
  double lon0 = -2.0;
  double SE, SN;
  double latmin = 49.8;
  double latmax = 61.0;
  double lonmin, lonmax;
  double step = FIT_STEP;

  double lat, lon;
  struct llh osgb, wgs;
  struct mxy mxy;
  double emin, emax, nmin, nmax;
  struct en en;
  double de, dn;
  long double ee[2*OX-1], nn[2*OY-1];
  long double LX[OXY], LY[OXY], RX[OXY], RY[OXY];
  Matrix M;
  int i, j, m, n;

  for (i=0; i<OXY; i++) {
    for (j=0; j<OXY; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<OXY; i++) {
    RX[i] = 0.0;
    RY[i] = 0.0;
  }

  osgb.lat = lat0;
  osgb.lon = lon0;
  osgb.h = 0.0;
  osgb36_to_grid(&osgb, &en);
  model->E0 = en.E;
  model->N0 = en.N;

  osgb.lat = 50;
  osgb.lon = -8;
  osgb.h = 0.0;
  osgb36_to_grid(&osgb, &en);
  emin = en.E;
  nmin = en.N;

  osgb.lat = 61;
  osgb.lon = +2;
  osgb.h = 0.0;
  osgb36_to_grid(&osgb, &en);
  emax = en.E;
  nmax = en.N;

  model->SE = SE = 2.0 / fabs(emax - emin);
  model->SN = SN = 2.0 / fabs(nmax - nmin);
  /* 2*pi*6378140*cos(51) */
  model->SX = 25 * 1.0e6;
  model->SY = 25 * 1.0e6;

#if 1
  printf("Original E0 : %.12f\n", model->E0);
  printf("Original N0 : %.12f\n", model->N0);
  printf("Original SE : %.12f (1.0/%.12f)\n", model->SE, 1.0/model->SE);
  printf("Original SN : %.12f (1.0/%.12f)\n", model->SN, 1.0/model->SN);
  printf("Original SX : %.12f\n", model->SX);
  printf("Original SY : %.12f\n", model->SY);
  model->E0 = 400000.0;
  model->N0 = 650000.0;
  model->SE = 1.0 / 400000.0;
  model->SN = 1.0 / 650000.0;
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
      de = model->SE * (en.E - model->E0);
      dn = model->SN * (en.N - model->N0);
      ee[0] = 1.0;
      for (i=1; i<2*OX-1; i++) {
        ee[i] = de * ee[i-1];
      }
      nn[0] = 1.0;
      for (i=1; i<2*OY-1; i++) {
        nn[i] = dn * nn[i-1];
      }
      for (i=0; i<OY; i++) {
        for (j=0; j<OY; j++) {
          for (m=0; m<OX; m++) {
            for (n=0; n<OX; n++) {
              M[OX*i+m][OX*j+n] += ee[m+n]*nn[i+j];
            }
          }
        }
      }
      for (i=0; i<OX; i++) {
        for (j=0; j<OY; j++) {
          RX[i+OX*j] += (model->SX * mxy.X) * ee[i] * nn[j];
          RY[i+OX*j] += (model->SY * mxy.Y) * ee[i] * nn[j];
        }
      }
    }
  }

  solve2(M, OXY, LX, LY, RX, RY);

  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->x[i][j] = LX[i+OX*j];
      model->y[i][j] = LY[i+OX*j];
    }
  }

}
/*}}}*/

static void predict(double e, double n,/*{{{*/
    const struct model *model,
    double *X, double *Y)
{
  double de, dn;
  double ee[2*OX-1], nn[2*OY-1];
  int i, j;
  double xx, yy;
  de = model->SE * (e - model->E0);
  dn = model->SN * (n - model->N0);
  ee[0] = 1.0;
  for (i=1; i<OX; i++) {
    ee[i] = de * ee[i-1];
  }
  nn[0] = 1.0;
  for (i=1; i<OY; i++) {
    nn[i] = dn * nn[i-1];
  }
  xx = yy = 0.0;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      xx += model->x[i][j] * ee[i] * nn[j];
      yy += model->y[i][j] * ee[i] * nn[j];
    }
  }
  *X = xx / model->SX;
  *Y = yy / model->SY;
}
/*}}}*/
static int predict_ll_ostn(double lat, double lon,/*{{{*/
    const struct model *model,
    double *X, double *Y,
    double *e_err, double *n_err)
{
  struct llh wgs;
  struct mxy mxy;
  struct en en_wgs;
  struct en en_osgb0, en_osgb1;

  wgs.lat = lat;
  wgs.lon = lon;
  wgs.h = 0;
  wgs84_to_grid(&wgs, &en_wgs);
  if (apply_ostn02(&en_wgs, &en_osgb0) == 0) return 0;
  predict(en_osgb0.E, en_osgb0.N, model, &mxy.X, &mxy.Y);
  *X = mxy.X;
  *Y = mxy.Y;
  mxy_to_wgs84(&mxy, &wgs);
  wgs84_to_grid(&wgs, &en_wgs);
  if (apply_ostn02(&en_wgs, &en_osgb1) == 0) return 0;

  *e_err = en_osgb1.E - en_osgb0.E;
  *n_err = en_osgb1.N - en_osgb0.N;
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
  struct mxy mxy0, mxy1;
  double de, dn;
  struct en en_wgs;
  struct en en_osgb;

  long double ee[2*QX-1], nn[2*QY-1];
  long double LX[QXY], LY[QXY], RX[QXY], RY[QXY];
  Matrix M;
  double x_err, y_err;
  int i, j, m, n;

  for (i=0; i<QXY; i++) {
    for (j=0; j<QXY; j++) {
      M[i][j] = 0.0;
    }
  }
  for (i=0; i<QXY; i++) {
    RX[i] = 0.0;
    RY[i] = 0.0;
  }

  for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
    uk_range(lat, &lonmin, &lonmax);
    for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
      wgs.lat = lat;
      wgs.lon = lon;
      wgs.h = 0;
      wgs84_to_grid(&wgs, &en_wgs);
      if (apply_ostn02(&en_wgs, &en_osgb) == 0) continue;
      wgs84_to_mxy(&wgs, &mxy0);
      predict(en_osgb.E, en_osgb.N, model, &mxy1.X, &mxy1.Y);
      x_err = mxy0.X - mxy1.X;
      y_err = mxy0.Y - mxy1.Y;
      de = model->SE * (en_osgb.E - model->E0);
      dn = model->SN * (en_osgb.N - model->N0);
      ee[0] = 1.0;
      for (i=1; i<2*QX-1; i++) {
        ee[i] = de * ee[i-1];
      }
      nn[0] = 1.0;
      for (i=1; i<2*QY-1; i++) {
        nn[i] = dn * nn[i-1];
      }

      for (i=0; i<QY; i++) {
        for (j=0; j<QY; j++) {
          for (m=0; m<QX; m++) {
            for (n=0; n<QX; n++) {
              M[QX*i+m][QX*j+n] += ee[m+n]*nn[i+j];
            }
          }
        }
      }

      for (i=0; i<QX; i++) {
        for (j=0; j<QY; j++) {
          RX[i+QX*j] += model->SX * x_err * ee[i] * nn[j];
          RY[i+QX*j] += model->SY * y_err * ee[i] * nn[j];
        }
      }

    }
  }

  solve2(M, QXY, LX, LY, RX, RY);

  printf("Trim coefficients:\n");
  for (i=0; i<QX; i++) {
    for (j=0; j<QY; j++) {
      model->x[i][j] += LX[i+QX*j];
      model->y[i][j] += LY[i+QX*j];
      printf("delta_e[%d][%d] = %15.5Lf   delta_n[%d][%d] = %15.5Lf\n",
          i, j,
          LX[i+QX*j],
          i, j,
          LY[i+QX*j]);
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
          LX[i+QX*j],
          LY[i+QX*j]
          );
    }
  }
  printf("\\bottomrule\n");
  printf("\\end{tabular}\n");
}
/*}}}*/
static void print_model(const struct model *model) {/*{{{*/
  int i, j;
  int nzx=0, nzy=0;
  printf("E0=%f\n", model->E0);
  printf("N0=%f\n", model->N0);
  printf("SE=%f (1/%f)\n", model->SE, 1.0/model->SE);
  printf("SN=%f (1/%f)\n", model->SN, 1.0/model->SN);
  printf("SX : %.12f\n", model->SX);
  printf("SY : %.12f\n", model->SY);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      printf("x[%d][%d] = %15.5f   y[%d][%d] = %15.5f\n",
          i, j,
          model->x[i][j],
          i, j,
          model->y[i][j]
          );
      if (fabs(model->x[i][j]) > 0.000001) nzx++;
      if (fabs(model->y[i][j]) > 0.000001) nzy++;
    }
  }
  printf("%2d X terms, %2d Y terms\n", nzx, nzy);
}
/*}}}*/
#if 0
static void print_model_tabular(const struct model *model) {/*{{{*/
  int i, j;
  printf("\\begin{tabular}{l l r r}\n");
  printf("\\toprule\n");
  printf("\\textbf{i} & \\textbf{j} & $x_{ij}$ & $y_{ij}$ \\\\\n");
  printf("\\midrule\n");
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      printf("%d & %d & %15.5f & %15.5f \\\\\n",
          i, j,
          model->x[i][j],
          model->y[i][j]
          );
    }
  }
  printf("\\bottomrule\n");
  printf("\\end{tabular}\n");
}
/*}}}*/
#endif
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
        if (j > 0) printf("*n");
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
      if (i > 0) printf("*e");
      if (i > 1) printf("%d", i);
    }
  }
  printf(";\n");
}
/*}}}*/
static void print_model_java(const struct model *model) {/*{{{*/
  print_model_java_inner('X', model->x);
  print_model_java_inner('Y', model->y);
}
/*}}}*/
static void print_model_estrin_inner(const double m[OX][OY], const char *res)
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
  estrin2(s, OX, OY, "e", "n", res);
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      if (s[i][j]) free(s[i][j]);
    }
  }

}

static void print_model_estrin(const struct model *model)
{
  print_model_estrin_inner(model->x, "X");
  print_model_estrin_inner(model->y, "Y");
}


/* Convert the coordinates to web mercator */
static void print_model_latex_inner(char foo, const double m[OX][OY], double output_scale) {/*{{{*/
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
        if (j > 0) printf("\\nu");
        if (j > 1) printf("^%d", j);
      }
    }
    if (any[i]) {
      printf("\\nonumber \\\\\n");
    }
  }
  printf("%.0f %c &= ", (output_scale+0.5), foo);
  for (i=0; i<OX; i++) {
    if (any[i]) {
      if (i > 0) printf(" +");
      printf("%c_%d", foo, i);
      if (i > 0) printf("\\eta");
      if (i > 1) printf("^%d", i);
    }
  }
  printf("\n");
}
/*}}}*/
static void print_model_latex(const struct model *model) {/*{{{*/
  /* TODO : show the translation and scaling on the domain side. */
  printf("\\begin{align}\n");
  print_model_latex_inner('X', model->x, model->SX);
  printf("\\\\[1ex]\n");
  print_model_latex_inner('Y', model->y, model->SY);
  printf("\\end{align}\n");
}
/*}}}*/
#if 0
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
    double coef;
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
      case 10:
        if (bw) {
#if (OY > 8)
          c[bi][8] += (1280.0/512.0)*coef;
          c[bi][6] -= (1120.0/512.0)*coef;
          c[bi][4] += (400.0/512.0)*coef;
          c[bi][2] -= (50.0/512.0)*coef;
          c[bi][0] += (1.0/512.0)*coef;
#endif
        } else {
#if (OX > 8)
          c[8][bj] += (1280.0/512.0)*coef;
          c[6][bj] -= (1120.0/512.0)*coef;
          c[4][bj] += (400.0/512.0)*coef;
          c[2][bj] -= (50.0/512.0)*coef;
          c[0][bj] += (1.0/512.0)*coef;
#endif
        }
        break;
      default:
        fprintf(stderr, "Can't handle a polynomial with order this high, giving up\n");
        exit(2);
    }
    c[bi][bj] = 0.0;
    accum[bi][bj] += coef;
  } while(1);

}
/*}}}*/
#if 0
static void round_model_1(struct model *model)/*{{{*/
{
  int i, j;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->x[i][j] = 0.1 * round(10.0 * model->x[i][j]);
      model->y[i][j] = 0.1 * round(10.0 * model->y[i][j]);
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
      model->x[i][j] = 0.01 * round(100.0 * model->x[i][j]);
      model->y[i][j] = 0.01 * round(100.0 * model->y[i][j]);
    }
  }
}
/*}}}*/
static void round_model_3(struct model *model)/*{{{*/
{
  int i, j;
  for (i=0; i<OX; i++) {
    for (j=0; j<OY; j++) {
      model->x[i][j] = 0.001 * round(1000.0 * model->x[i][j]);
      model->y[i][j] = 0.001 * round(1000.0 * model->y[i][j]);
    }
  }
}
/*}}}*/
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

#if 0
static void test_coord(void)/*{{{*/
{
  struct llh wgs0, wgs1;
  struct mxy xy;

  wgs0.lat = 50.0;
  wgs0.lon = -8.0;
  wgs0.h = 0.0;
  wgs84_to_mxy(&wgs0, &xy);
  mxy_to_wgs84(&xy, &wgs1);
  printf("(%.12f,%.12f) -> (%.12f,%.12f) -> (%.12f,%.12f)\n",
      wgs0.lat, wgs0.lon,
      xy.X, xy.Y,
      wgs1.lat, wgs1.lon);
  wgs0.lat = 61.0;
  wgs0.lon = 2.0;
  wgs0.h = 0.0;
  wgs84_to_mxy(&wgs0, &xy);
  mxy_to_wgs84(&xy, &wgs1);
  printf("(%.12f,%.12f) -> (%.12f,%.12f) -> (%.12f,%.12f)\n",
      wgs0.lat, wgs0.lon,
      xy.X, xy.Y,
      wgs1.lat, wgs1.lon);
  wgs0.lat = 54.5012345678;
  wgs0.lon = -0.0123456789;
  wgs0.h = 0.0;
  wgs84_to_mxy(&wgs0, &xy);
  mxy_to_wgs84(&xy, &wgs1);
  printf("(%.12f,%.12f) -> (%.12f,%.12f) -> (%.12f,%.12f)\n",
      wgs0.lat, wgs0.lon,
      xy.X, xy.Y,
      wgs1.lat, wgs1.lon);
}
/*}}}*/
#endif

struct accuracy
{
  double x;
  double y;
  double minc;
  char *base;
};

#define N_ACCURACIES 4

static struct accuracy accuracies[N_ACCURACIES] =
{
  {   2.5,   2.0,   0.35, "grid_to_wm_1m%s.%s" },
  {  30.0,  28.0,   1.0,  "grid_to_wm_10m%s.%s" },
  { 100.0, 150.0,   3.5,  "grid_to_wm_35m%s.%s" },
  { 200.0, 150.0,  10.0,  "grid_to_wm_100m%s.%s" }
};

int main (int argc, char **argv)/*{{{*/
{
  struct model model;
  double step;
  double lat, lon;
  struct cdata *cdata_notrim, *cdata_total, *cdata_east, *cdata_north;
  struct crow *crow_notrim, *crow_total, *crow_east, *crow_north;
  struct cline *lines;
  double lonmin, lonmax;
  double latmin, latmax;
  int i, acc;
  double X, Y;
  int is_svg;

#if 0
  test_coord();
  exit (0);
#endif

  is_svg = (argc > 1);

  load_osxx02();

  fit(&model);

  printf("===============================\n");
  printf("INITIAL MODEL FIT:\n");
  print_model(&model);

  printf("===============================\n");
  step = 0.025;
  latmin = 49.0;
  latmax = 61.0;
  lonmin = -8.0;
  lonmax = +2.0;

  if (!is_svg) {
    cdata_notrim = new_cdata();
    for (lat=latmin; lat<=latmax+step/2.0; lat+=step) {
      crow_notrim = next_crow(cdata_notrim);
      for (lon=lonmin; lon<=lonmax+step/2.0; lon+=step) {
        double e_err, n_err;
        double t_err;
        if (predict_ll_ostn(lat, lon, &model, &X, &Y, &e_err, &n_err) == 0) {
          add_empty_cnode(crow_notrim);
        } else {
          t_err = sqrt(e_err*e_err + n_err*n_err);
          add_cnode(crow_notrim, t_err);
        }
      }
    }

    start_tikz("grid_to_wm_notrim.tex");
    for (i=0; i<N_LEVELS_1; i++) {
      if (levels_1side[i].level < 0.5) continue;
      lines = generate_isolines(cdata_notrim, levels_1side[i].level);
      remap_clines(lines, lonmin, step, latmin, step);
      emit_tikz(levels_1side[i].level, levels_1side[i].colour, levels_1side[i].thickness, lines);
      free_clines(lines);
    }
    finish_tikz();

    free_cdata(cdata_notrim);
  }

  printf("===============================\n");
  trim_model(&model);
  printf("AFTER TRIMMING:\n");
  print_model(&model);

  for (acc=0; acc<N_ACCURACIES; acc++) {
    struct model copy;
    char filename[64];
    memcpy(&copy, &model, sizeof(struct model));
    printf("===============================\n");
    printf("PUNCTURE : x_thresh=%f y_thresh=%f\n",
        accuracies[acc].x,
        accuracies[acc].y);
    puncture_model_4("X", accuracies[acc].x, copy.x);
    puncture_model_4("Y", accuracies[acc].y, copy.y);

    printf("=================\n");
    printf("AFTER PUNCTURING:\n");
    print_model(&copy);

    round_model_3(&copy);
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
        if (predict_ll_ostn(lat, lon, &copy, &X, &Y, &e_err, &n_err) == 0) {
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

    if (is_svg) {
      sprintf(filename, accuracies[acc].base, "", "svg");
      start_svg(filename);
    } else {
      sprintf(filename, accuracies[acc].base, "", "tex");
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

    if (is_svg) {
      sprintf(filename, accuracies[acc].base, "_east", "svg");
      start_svg(filename);
    } else {
      sprintf(filename, accuracies[acc].base, "_east", "tex");
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

    if (is_svg) {
      sprintf(filename, accuracies[acc].base, "_north", "svg");
    } else {
      sprintf(filename, accuracies[acc].base, "_north", "tex");
    }
    start_tikz(filename);
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
