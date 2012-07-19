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

/* Read in the OSTN02, OSGM02 models from file and provide a lookup function */
#include <stdio.h>
#include <stdlib.h>

#include "tool.h"
#include "osxx02.h"

struct point/*{{{*/
{
  /* one point in the ostn02/osgm02 model */
  int flag;
  double e_shift;
  double n_shift;
  double g_height;
};
/*}}}*/

#define N_EAST 701
#define N_NORTH 1251

static struct point ostn[N_NORTH][N_EAST];

static inline struct point *lookup_point(int ee, int nn)/*{{{*/
{
  return &ostn[nn][ee];
}
/*}}}*/
void load_osxx02(void)/*{{{*/
{
  int e, n;
  char line[256];
  int R, E, N, F;
  double ES, NS, GH;
  FILE *in;

  in = fopen("OSTN02_OSGM02_GB.txt", "r");
  if (!in) {
    fprintf(stderr, "Could not open OSTN02_OSGM02_GB.txt\n");
    exit(2);
  }

  e = n = 0;

  while (fgets(line, sizeof(line), in)) {
    if (sscanf(line, "%d,%d,%d,%lf,%lf,%lf,%d",
          &R, &E, &N,
          &ES, &NS, &GH,
          &F) != 7) {
      fprintf(stderr, "Failed to read record at e=%d, n=%d\n", e, n);
      exit(2);
    }
    if (E != 1000*e) {
      fprintf(stderr, "Eastings mismatch at e=%d, n=%d\n", e, n);
      exit(2);
    }
    if (N != 1000*n) {
      fprintf(stderr, "Northings mismatch at e=%d, n=%d\n", e, n);
      exit(2);
    }
    ostn[n][e].e_shift = ES;
    ostn[n][e].n_shift = NS;
    ostn[n][e].g_height = GH;
    ostn[n][e].flag = F;
    e++;
    if (e == N_EAST) {
      e = 0;
      n++;
    }
  }

  fclose(in);

}
/*}}}*/
int apply_ostn02(const struct en *en_wgs, struct en *en_osgb)/*{{{*/
{
  int e0, n0;
  int e1, n1;
  double de, dn;
  double se, sn;
  e0 = (int)(0.001 * en_wgs->E);
  n0 = (int)(0.001 * en_wgs->N);
  e1 = e0 + 1;
  n1 = n0 + 1;
  de = (0.001 * en_wgs->E) - (double)e0;
  dn = (0.001 * en_wgs->N) - (double)n0;
  if (e0 < 0) return 0;
  if (n0 < 0) return 0;
  if (e1 >= N_EAST) return 0;
  if (n1 >= N_NORTH) return 0;
  if (lookup_point(e0, n0)->flag == 0) return 0;
  if (lookup_point(e1, n0)->flag == 0) return 0;
  if (lookup_point(e0, n1)->flag == 0) return 0;
  if (lookup_point(e1, n1)->flag == 0) return 0;
  se = (1.0-de) * (1.0-dn) * lookup_point(e0, n0)->e_shift
     + (    de) * (1.0-dn) * lookup_point(e1, n0)->e_shift
     + (1.0-de) * (    dn) * lookup_point(e0, n1)->e_shift
     + (    de) * (    dn) * lookup_point(e1, n1)->e_shift;
  sn = (1.0-de) * (1.0-dn) * lookup_point(e0, n0)->n_shift
     + (    de) * (1.0-dn) * lookup_point(e1, n0)->n_shift
     + (1.0-de) * (    dn) * lookup_point(e0, n1)->n_shift
     + (    de) * (    dn) * lookup_point(e1, n1)->n_shift;
  en_osgb->E = en_wgs->E + se;
  en_osgb->N = en_wgs->N + sn;
  return 1;
}
/*}}}*/

int apply_osgm02(const struct en *en_wgs, int only_newlyn, double *corr)
{
  int e0, n0;
  int e1, n1;
  double de, dn;
  e0 = (int)(0.001 * en_wgs->E);
  n0 = (int)(0.001 * en_wgs->N);
  e1 = e0 + 1;
  n1 = n0 + 1;
  de = (0.001 * en_wgs->E) - (double)e0;
  dn = (0.001 * en_wgs->N) - (double)n0;
  if (e0 < 0) return 0;
  if (n0 < 0) return 0;
  if (e1 >= N_EAST) return 0;
  if (n1 >= N_NORTH) return 0;

  if (only_newlyn) {
    /* Only want ODN points */
    if (lookup_point(e0, n0)->flag != 1) return 0;
    if (lookup_point(e1, n0)->flag != 1) return 0;
    if (lookup_point(e0, n1)->flag != 1) return 0;
    if (lookup_point(e1, n1)->flag != 1) return 0;
  } else {
    if (lookup_point(e0, n0)->flag == 0) return 0;
    if (lookup_point(e1, n0)->flag == 0) return 0;
    if (lookup_point(e0, n1)->flag == 0) return 0;
    if (lookup_point(e1, n1)->flag == 0) return 0;
  }
  *corr = (1.0-de) * (1.0-dn) * lookup_point(e0, n0)->g_height
        + (    de) * (1.0-dn) * lookup_point(e1, n0)->g_height
        + (1.0-de) * (    dn) * lookup_point(e0, n1)->g_height
        + (    de) * (    dn) * lookup_point(e1, n1)->g_height;
  return 1;
}
