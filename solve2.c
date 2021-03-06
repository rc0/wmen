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
#include <math.h>

#include "tool.h"

static void reduce(Matrix m, int n, long double *r0, long double *r1)
{
  int i, j, k;
  int ii;
  long double t, scale;
  for (i=0; i<n; i++) {
    ii = -1;
    for (j=i+1; j<n; j++) {
      if (fabs(m[j][i] > fabs(m[i][i]))) {
        ii = j;
      }
    }
    if (ii > i) {
      /* swap */
      for (k=i; k<n; k++) {
        t = m[ii][k];
        m[ii][k] = m[i][k];
        m[i][k] = t;
      }
      t = r0[ii];
      r0[ii] = r0[i];
      r0[i] = t;
      t = r1[ii];
      r1[ii] = r1[i];
      r1[i] = t;
    }
    /* Now subtract multiples of row i from the rows below */
    for (j=i+1; j<n; j++) {
      if (m[i][i] == 0.0) {
        printf("scale is 0 at %d,%d\n", i, i);
        exit(1);
      }
      scale = m[j][i] / m[i][i];
      for (k=i; k<n; k++) {
        m[j][k] -= scale * m[i][k];
      }
      r0[j] -= scale * r0[i];
      r1[j] -= scale * r1[i];
    }
  }
}

static void backfill(Matrix m, int n,
    long double *l0, long double *l1,
    long double *r0, long double *r1)
{
  int i, j;
  for (i=n-1; i>=0; i--) {
    long double t0, t1;
    t0 = r0[i];
    t1 = r1[i];
    for (j=i+1; j<n; j++) {
      t0 -= m[i][j] * l0[j];
      t1 -= m[i][j] * l1[j];
    }
    l0[i] = t0 / m[i][i];
    l1[i] = t1 / m[i][i];
  }
}

void solve2(Matrix m, int n, long double *l0, long double *l1,
    long double *r0, long double *r1)
{
  /* Row reduction with partial pivot */
  reduce(m, n, r0, r1);
  backfill(m, n, l0, l1, r0, r1);
}


