/*
* Copyright (c) 2013, Richard P. Curnow
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

/* Find approximation to distance between two points.  Actually, the scaling
 * per unit web-mercator space, as a function of the Y-coordinate

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "tool.h"

static double a = 6378137.0;
static double b = 6356752.3141;

static double latofy(double y) {
  double t = (1.0 - 2.0*y) * M_PI;
  double M = (2.0 * atan(exp(t))) - 0.5*M_PI;
  return M;
}

static double latdegofy(double y) {
  double M = latofy(y);
  return M * (180.0/M_PI);
}

static double fy(double y)
{
  /* given y in [0,1] as the Y-coord of the web mercator position,
   * compute the 'radius', equivalent to Re * cos(lat) */

  double M = latofy(y);
  double e2 = 1.0 - (b*b)/(a*a);
  double SM = sin(M);
  double nu = a / sqrt(1.0 - e2*SM*SM);
  return nu * cos(M);
}

static void do_fit(int n, double *coef)
{
  Matrix m;
  long double c[32], r[32];
  double nodes[32];
  int i, j;

  for (i=1; i<=n; i++) {
    /* In range [-1,1] */
    double node = cos(M_PI*(double)(i+i-1)/(double)(n+n));

    /* Only cover the range 71N-71S */
    node = (0.5-0.215)*node + 0.5;

#if 0
    /* Scale/shift to range [0,1] */
    node = 0.5*node + 0.5;
#endif
    printf("%2d : %10.5f\n", i, node);
    nodes[i-1] = node;
  }

  for (i=0; i<n; i++) {
    double nn = nodes[i];
    double prod = 1.0;
    double z = nn - 0.5;
    for (j=0; j<n; j++) {
      m[i][j] = prod;
      prod *= z;
    }
    r[i] = a/fy(nn);
  }

  solve(m, n, c, r);
  for (i=0; i<n; i++) {
    coef[i] = c[i];
  }
  return;

}

int main (int argc, char **argv)
{
#if 1
  double coef[32];
  int i;
  int n = 7;

  do_fit(n, coef);
  for (i=0; i<n; i++) {
    printf("%2d : %15.8f\n", i, coef[i]);
  }
#endif

#if 1
  double y;
  double z, z2, z4, z6, fz;
  double fyy;
  double err;
#if 0
  double k2 = 19.74765212;
  double k4 = 61.09368648;
  double k6 = 114.75590338;
#endif
#if 0
  double k2 = 19.7;
  double k4 = 61.0;
  double k6 = 115;
#endif
  double k2 = coef[2];
  double k4 = coef[4];
  double k6 = coef[6];
  y = 0.0;
  while (y < 1.0) {
    z = (y - 0.5);
    z2 = z*z;
    z4 = z2*z2;
    z6 = z4*z2;
    fz = a / (1.0 + k2*z2 + k4*z4 + k6*z6);
    fyy = fy(y);
    err = (fyy - fz)/fyy;

    printf("%15.8f %10.2f %15.8f %15.8f %15.8f %8.2f\n",
        y, latdegofy(y), fyy, fz, err, err*1000);
    y += 0.01;
  }
#endif

  return 0;
}


